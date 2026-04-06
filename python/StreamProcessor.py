#!/usr/bin/env python3
# Main script for real-time DAS data processing and streaming

import socket, zmq
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import sys
import atexit
from DASPacket_Mod import *
from DAS_RealTimeMod import (
    RingBuffer,
    real_time_picking_async,
    start_picking_thread,
    merge_stream_picks,
    time_format,
    init_pga_state,
    update_real_time_pga,
)
import time
from datetime import datetime, timedelta, timezone
from obspy.core.inventory import inventory
import asyncio
import PyEW


def _ensure_stomp_client_importable(project_root):
    """Add the read-only stomp_client dependency to sys.path if needed."""
    stomp_path = Path(project_root) / 'external' / 'stomp_client'
    if not stomp_path.exists():
        raise ImportError(f'stomp_client path not found: {stomp_path}')
    stomp_path_str = str(stomp_path)
    if stomp_path_str not in sys.path:
        sys.path.insert(0, stomp_path_str)


def _validate_pga2finder_args(args):
    """Validate PGA-to-FinDer runtime configuration."""
    if not args.PGA2FinDer:
        return

    if args.workInterval <= 0.0:
        raise ValueError('--PGA2FinDer requires --workInterval > 0')
    if args.xmlmeta is None:
        raise ValueError('--PGA2FinDer requires --xmlmeta')
    if args.pickingChannel is None:
        raise ValueError('--PGA2FinDer requires --pickingChannel')
    if args.finderConfig is None:
        raise ValueError('--PGA2FinDer requires --finderConfig')
    if args.nChPGASmooth < 0 or args.nChPGASmooth % 2 != 0:
        raise ValueError('--nChPGASmooth must be a non-negative even integer')

    picking_channel_df = pd.read_csv(args.pickingChannel, nrows=1)
    required_cols = {'Channel', 'PGA/PSR-Ratio'}
    missing_cols = required_cols.difference(picking_channel_df.columns)
    if missing_cols:
        raise ValueError(
            'pickingChannel file is missing required columns: %s'
            % ', '.join(sorted(missing_cols))
        )

    _ensure_stomp_client_importable(Path(__file__).resolve().parents[1])


def init_finder_sender(config_path, verbose_level=1):
    """Initialize and connect a FinDer sender using the read-only stomp client dependency."""
    _ensure_stomp_client_importable(Path(__file__).resolve().parents[1])
    from finder_sender import FinderMessageSender
    from stomp_client import StompConfig, StompConnection, LoggingListener

    cfg = StompConfig(config_path)
    conn = StompConnection(cfg)
    listener = LoggingListener(conn, verbose_level=verbose_level)
    conn.set_listener(listener)
    sender = FinderMessageSender(conn, verbose_level=verbose_level)
    sender.connect()
    return sender


def _coerce_timestamp_to_epoch(value):
    """Convert various timestamp-like objects to Unix epoch seconds."""
    ts = pd.Timestamp(value)
    if ts.tzinfo is None:
        ts = ts.tz_localize(timezone.utc)
    else:
        ts = ts.tz_convert(timezone.utc)
    return ts.timestamp()


def build_finder_pga_payload(pga_rows, pga_state):
    """Build a FinDer-compatible payload and metadata from PGA rows."""
    if pga_rows is None or len(pga_rows) == 0:
        raise ValueError('No PGA rows available for FinDer payload')

    message_timestamp = datetime.now(timezone.utc).timestamp()
    stations = {}
    missing_geo = []
    for _, row in pga_rows.iterrows():
        static_station = pga_state.finder_station_template[row['sncl']]
        if static_station['lat'] is None or static_station['lon'] is None:
            missing_geo.append(row['channel_info'])
            continue
        station_timestamp = _coerce_timestamp_to_epoch(row['pga_time'])
        stations[row['sncl']] = {
            **static_station,
            'timestamp': station_timestamp,
            'PGA': f"{float(row['pga']):.6e}",
        }

    if missing_geo:
        raise ValueError(
            'Missing latitude/longitude for channels: %s'
            % ', '.join(missing_geo)
        )

    payload = {'stations': stations}
    metadata = dict(pga_state.finder_metadata)
    return payload, metadata, message_timestamp


def build_finder_debug_message(payload, metadata, message_timestamp, config_path):
    """Build the exact outbound FinDer FFD2 text without requiring a live broker connection."""
    _ensure_stomp_client_importable(Path(__file__).resolve().parents[1])
    from finder_sender import FinderMessageSender
    from stomp_client import StompConfig, StompConnection

    cfg = StompConfig(config_path)
    conn = StompConnection(cfg)
    sender = FinderMessageSender(conn, verbose_level=0)
    envelope_meta = dict(sender.metadata)
    envelope_meta.update(metadata)
    return sender.build_ffd2_format(
        envelope_meta,
        payload['stations'],
        timestamp=message_timestamp,
    )


def write_pga_debug_outputs(pga_rows, finder_message, output_dir, timestamp):
    """Write PGA CSV and exact outbound FinDer text payload for debugging."""
    if output_dir is None:
        return

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = pd.Timestamp(datetime.fromtimestamp(timestamp, tz=timezone.utc)).strftime('%Y-%m-%dT%H%M%S.%fZ')

    csv_path = outdir / f'{stamp}_pga.csv'
    txt_path = outdir / f'{stamp}_finder.txt'

    rows_to_write = pga_rows.copy()
    if 'peak_strain_rate_time' in rows_to_write.columns:
        rows_to_write['peak_strain_rate_time'] = rows_to_write['peak_strain_rate_time'].apply(
            lambda x: pd.Timestamp(x).isoformat() if pd.notna(x) else None
        )
    if 'pga_time' in rows_to_write.columns:
        rows_to_write['pga_time'] = rows_to_write['pga_time'].apply(
            lambda x: pd.Timestamp(x).isoformat() if pd.notna(x) else None
        )
    rows_to_write.to_csv(csv_path, index=False)
    txt_path.write_text(finder_message)


def _disconnect_finder_sender(sender):
    """Disconnect the FinDer sender if it exists."""
    if sender is None:
        return
    try:
        sender.disconnect()
    except Exception as exc:
        print(f'Error disconnecting FinDer sender: {exc}', flush=True)


def _schedule_finder_reconnect_cooldown(cooldown=None):
    """Delay the next FinDer reconnect attempt to avoid blocking every work interval."""
    global finderReconnectNotBefore
    if cooldown is None:
        cooldown = FINDER_RECONNECT_COOLDOWN
    finderReconnectNotBefore = time.monotonic() + cooldown


def _finder_reconnect_ready():
    """Return True when the FinDer sender reconnect cooldown has expired."""
    return time.monotonic() >= finderReconnectNotBefore


def doWork(strmRdr, args, waveRing=None, pickRing=None, loop=None, minimumPhaseNetTime=30.0):
    global finderSender, finderReconnectNotBefore

    try:
        workInterval = args.workInterval
        strainRate = args.strainRate
        ringbuff_size = args.ringbuffer
        filelength = args.filelength
        filepath = args.filepath
        pickOutInt = args.pickOutputInterval
        if args.xmlmeta is None and filelength > 0.0:
            raise ValueError('User must provide XML meta information')
        if args.xmlmeta is not None:
            inventor = inventory.read_inventory(args.xmlmeta)
        else:
            inventor = None
        if workInterval > ringbuff_size:
            raise ValueError('workInterval (%s) must be smaller than or equal to ringbuffer size (%s)' % (workInterval, ringbuff_size))
        taskPicking = None  # Picking task to check if previous task is done
        TT_picksBuf = None
        pgaState = None
        pickingChannel_df = None

        ii = 0  # packet counter
        packet = strmRdr.getNextPacket()
        fs = strmRdr.getFs(packet)
        deltaStrainRate = timedelta(seconds=float(0.5 / fs))
        nch = strmRdr.getNumChannel(packet)
        pickingChannel = np.arange(nch)
        if args.pickingChannel is not None:
            pickingChannel_df = pd.read_csv(args.pickingChannel)
            pickingChannel = pickingChannel_df['Channel'].values.astype(int)
            print(f'Picking will be performed on {pickingChannel.shape[0]} channels: {pickingChannel}', flush=True)

        # Picking output folder
        pickOutput = args.pickOutput
        startPicking = time.time()

        ringbuff = RingBuffer(int(ringbuff_size * fs), pickingChannel)
        ringbuff.setObspyTraceHeader(inventor)
        streamCh = None
        chCodes = None
        if ringbuff.channels_info is not None:
            streamCh = np.array(ringbuff.chIds)
            chCodes = ringbuff.statNames

        if args.PGA2FinDer:
            pgaState = init_pga_state(
                pickingChannel_df,
                ringbuff.good_ch,
                ringbuff.channels_info,
                ringbuff.chIds,
                ringbuff.rawChIds,
                ringbuff.statNames,
                ringbuff.chIds,
                ringbuff.stats,
                args.nChPGASmooth,
            )

        # Factors to convert phase to strain
        conv_factor = strmRdr.getConversionFactor(packet)
        if strainRate:
            conv_factor *= fs  # Division by dt

        OldtimeSample = strmRdr.getPayloadRad(packet) * conv_factor
        if strainRate:
            packet = strmRdr.getNextPacket()
            currtimeSample = strmRdr.getPayloadRad(packet) * conv_factor
            ringbuff.append(currtimeSample - OldtimeSample, timestamps=strmRdr.getPacketTimestamp(packet) - deltaStrainRate)
            OldtimeSample = currtimeSample
        else:
            ringbuff.append(OldtimeSample, timestamps=strmRdr.getPacketTimestamp(packet))
        ii += strmRdr.getNumTimeSamples(packet)
        while True:
            packet = strmRdr.getNextPacket()
            if packet == b'':
                break

            # Computing time step difference between last sample and latest sample received to check for data gaps
            if strmRdr.getNumTimeSamples(packet) > 1:
                curr_timestamp = strmRdr.getPacketTimestamp(packet)[0]
                last_timestamp = ringbuff.getTimeStamps()[-1]
            else:
                curr_timestamp = strmRdr.getPacketTimestamp(packet)
                last_timestamp = ringbuff.getTimeStamps()[-1]
            timediff = (curr_timestamp - last_timestamp).total_seconds()
            if strainRate:
                timediff -= float(0.5 / fs)  # to account for origin shift of the strain rate

            if ~np.isclose(1.0 / fs, timediff, atol=1e-3):
                print('Skipped samples', flush=True)
                print(f'Packet current and last timestamps: {last_timestamp}, {curr_timestamp}; time gap: {timediff}', flush=True)
                if filelength > 0.0 and filepath is not None:
                    print('Writing file at %s' % ringbuff.getTimeStamps()[-1].strftime(time_format), flush=True)
                    ringbuff.writeObsPyTraces(fs, filepath)
                time.sleep(0.1)
                return

            # Filling the buffer
            currtimeSample = strmRdr.getPayloadRad(packet) * conv_factor
            if strainRate:
                ringbuff.append(currtimeSample - OldtimeSample, timestamps=strmRdr.getPacketTimestamp(packet) - deltaStrainRate)
                OldtimeSample = currtimeSample
            else:
                ringbuff.append(currtimeSample, timestamps=strmRdr.getPacketTimestamp(packet))

            ii += strmRdr.getNumTimeSamples(packet)

            t_each_print = 1.0  # seconds
            if ii % int(t_each_print * fs) == 0:
                print(
                    f'Packet first and last timestamp and shape of ringbuffer: {ringbuff.getTimeStamps()[0]}, '
                    f'{ringbuff.getTimeStamps()[-1]}, {ringbuff.getData().shape}',
                    flush=True,
                )

            # Do some processing every workInterval
            if workInterval > 0.0 and ii % int(workInterval * fs) == 0 and ii > 0:
                if args.PGA2FinDer and pgaState is not None:
                    pga_rows = update_real_time_pga(ringbuff, pgaState)
                    payload, metadata, message_timestamp = build_finder_pga_payload(pga_rows, pgaState)

                    finder_message = None
                    if args.PGAOutput is not None:
                        try:
                            finder_message = build_finder_debug_message(
                                payload,
                                metadata,
                                message_timestamp,
                                args.finderConfig,
                            )
                        except Exception as exc:
                            print(f'Failed to build local FinDer debug message: {exc}', flush=True)

                    if finderSender is None and _finder_reconnect_ready():
                        try:
                            finderSender = init_finder_sender(args.finderConfig, verbose_level=max(0, args.debug))
                            finderReconnectNotBefore = 0.0
                        except Exception as exc:
                            print(f'Failed to initialize FinDer sender: {exc}', flush=True)
                            finderSender = None
                            _schedule_finder_reconnect_cooldown()

                    if finderSender is not None:
                        try:
                            finderSender.send_finder(
                                payload=payload,
                                content_type='text/plain',
                                metadata=metadata,
                                timestamp=message_timestamp,
                            )
                        except Exception as exc:
                            print(f'Failed to send FinDer message: {exc}', flush=True)
                            _disconnect_finder_sender(finderSender)
                            finderSender = None
                            _schedule_finder_reconnect_cooldown()

                    if args.PGAOutput is not None and finder_message is not None:
                        write_pga_debug_outputs(pga_rows, finder_message, args.PGAOutput, message_timestamp)

                # Picking has been requested if loop is not None
                if loop is not None:
                    # Performing some checks before submitting picking task
                    if taskPicking is not None:
                        # Previous picking done?
                        if not taskPicking.done():
                            continue
                        TT_picksNew = taskPicking.result()
                        # Picking streaming task
                        TT_picksBuf = merge_stream_picks(
                            TT_picksBuf,
                            TT_picksNew,
                            delta_t_thres=2.0,
                            maxBuf=pickOutInt,
                            pickRing=pickRing,
                            streamCh=streamCh,
                            chCodes=chCodes,
                        )
                        # Checking if writing of picking data base was requested
                        if time.time() - startPicking >= pickOutInt and pickOutput is not None:
                            if TT_picksBuf is not None:
                                filename = pickOutput + '/%s.csv' % datetime.fromtimestamp(startPicking).astimezone(timezone.utc).strftime(time_format)
                                print('Writing picking file %s' % filename, flush=True)
                                TT_picksBuf = TT_picksBuf.drop_duplicates(subset=['station_id', 'phase_time'])
                                TT_picksBuf.to_csv(filename, index=False)
                            startPicking = time.time()
                    if ringbuff.getData().shape[1] < minimumPhaseNetTime * fs:
                        continue
                    taskPicking = real_time_picking_async(ringbuff.getData().copy(), 1.0 / fs, ringbuff.getTimeStamps())
                    taskPicking = asyncio.run_coroutine_threadsafe(taskPicking, loop)

            if filelength > 0.0 and ii % int(ringbuff_size * fs) == 0 and ii > 0 and filepath is not None:
                print('Writing file at %s' % ringbuff.getTimeStamps()[-1].strftime(time_format), flush=True)
                ringbuff.writeObsPyTraces(fs, filepath)

            if waveRing is not None and ii % int(ringbuff_size * fs) == 0:
                ringbuff.send2ew(fs, waveRing)

    except Exception as e:
        print(e)
        pass


waveRing = None
pickRing = None
loop = None
finderSender = None
finderReconnectNotBefore = 0.0
FINDER_RECONNECT_COOLDOWN = 10.0


def main():
    """Main function to perform picking and earthworm picks and data streaming."""
    global waveRing, pickRing, loop, finderSender, finderReconnectNotBefore

    _validate_pga2finder_args(args)

    if not (1 <= args.port <= 65535):
        error_exit('Invalid input port number')

    strmType = args.streamType

    if strmType == 'OptaSense':
        try:
            socket.inet_aton(args.host)
        except socket.error:
            error_exit('Invalid input host IP address')
    elif strmType == 'ASN':
        try:
            context = zmq.Context()
        except Exception as e:
            error_exit(f'Failed to initialize ZMQ context: {e}')
    else:
        raise ValueError(f'ERROR! streamType provided ({strmType}) not supported!')

    if args.wavering[0] and waveRing is None:
        ringNumber, modID, inst_id, hb_freq, db_ew = args.wavering
        waveRing = PyEW.EWModule(ringNumber, modID, inst_id, float(hb_freq), bool(db_ew))
        waveRing.add_ring(ringNumber)

    if args.pickring[0] and pickRing is None:
        ringNumber, modID, inst_id, hb_freq, db_ew = args.pickring
        pickRing = PyEW.EWModule(ringNumber, modID, inst_id, float(hb_freq), bool(db_ew))
        pickRing.add_ring(ringNumber)

    if args.picking and loop is None:
        loop = start_picking_thread(args.device)

    if args.PGA2FinDer and finderSender is None and _finder_reconnect_ready():
        try:
            finderSender = init_finder_sender(args.finderConfig, verbose_level=max(0, args.debug))
            finderReconnectNotBefore = 0.0
        except Exception as exc:
            print(f'Failed to initialize FinDer sender at startup: {exc}', flush=True)
            finderSender = None
            _schedule_finder_reconnect_cooldown()

    atexit.register(lambda: _disconnect_finder_sender(finderSender))

    print(f'Starting data streaming for {strmType}...', flush=True)

    while True:
        try:
            if strmType == 'OptaSense':
                with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as inp:
                    print(f'Connecting to OptaSense data stream {args.host}:{args.port}...', flush=True)
                    inp.connect((args.host, args.port))
                    print('Connected to OptaSense', flush=True)
                    strmRdr = OptaSenseStreamReader(inp)
                    doWork(strmRdr, args, loop=loop, waveRing=waveRing, pickRing=pickRing)

            elif strmType == 'ASN':
                with context.socket(zmq.SUB) as inp:
                    print(f'Connecting to ASN data stream on {args.host}:{args.port}...', flush=True)
                    inp.connect(f'tcp://{args.host}:{args.port}')
                    inp.setsockopt(zmq.SUBSCRIBE, b'')
                    print('Connected to ASN', flush=True)
                    strmRdr = ASN_StreamReader(inp)
                    doWork(strmRdr, args, loop=loop, waveRing=waveRing, pickRing=pickRing)

        except zmq.ZMQError as e:
            print(f'ASN stream error: {e}', flush=True)
            time.sleep(5)

        except Exception as e:
            print(f'Error connecting to {strmType} stream: {e}', flush=True)
            time.sleep(5)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Main program to process real-time DAS data streams.\n')
    parser.add_argument('--host', metavar='HOST', required=True, help='a hostname of the input data stream')
    parser.add_argument('--port', metavar='PORT', type=int, required=True, help='a port of the input data stream')
    parser.add_argument('--streamType', '-strTp', metavar='streamType', type=str, required=True, help='Stream data type. Currently supported: OptaSense, ASN')
    parser.add_argument('--strainRate', '-strnRt', metavar='strainRate', type=int, default=1, help='Save strain rate or strain. Default 1')
    parser.add_argument('--workInterval', '-wrkint', metavar='workInterval', type=float, default=1.0, help='Work interval for data processing. Default 1.0 [s]')
    parser.add_argument('--picking', '-pick', metavar='picking', type=int, default=0, help='Flag to run picking using PhaseNet-DAS (Zhu et al., 2023)')
    parser.add_argument('--pickingChannel', '-pch', metavar='pickingChannel', type=str, default=None, help='CSV file containing the reduced DAS channel subset in its `Channel` column plus any per-channel metadata used during processing.')
    parser.add_argument('--pickOutput', '-pOut', metavar='pickOutput', type=str, default=None, help='Folder where picking database are written. Currently hourly output.')
    parser.add_argument('--pickOutputInterval', '-pOutInt', metavar='pickOutInt', type=float, default=3600.0, help='Interval to write picks into .csv file. Default 3600.0 [s]')
    parser.add_argument('--pickring', '-pring', metavar='pickring', type=int, nargs=5, default=[None, 8, 141, 30, 0], help='PyEarthworm module parameters for pick ring (ring number,module ID, INST_ID,HeartBeats frequency [s], Debug Flag); see parameters for PyEW.EWModule')
    parser.add_argument('--device', '-dev', metavar='device', type=str, default='cuda', help='Device on which to run PhaseNet-DAS; picking must be 1 to work. To use different GPU card use cuda:1 to select GPU card ID 1')
    parser.add_argument('--ringbuffer', '-Rbfsz', metavar='ringbuffer', type=float, default=60.0, help='Ring buffer size in seconds to be stored for processing. Default 60.0 [s]')
    parser.add_argument('--filelength', '-flng', metavar='filelength', type=float, default=0.0, help='Interval to be writing data to disk in seconds. Default 0.0 [s], meaning no file writing')
    parser.add_argument('--filepath', '-flpt', metavar='filepath', type=str, default=None, help='Path for writing the mseed files')
    parser.add_argument('--xmlmeta', '-xml', metavar='xmlmeta', type=str, default=None, help='Path to the XML metadata channel info')
    parser.add_argument('--wavering', '-wring', metavar='wavering', type=int, nargs=5, default=[None, 8, 141, 30, 0], help='PyEarthworm module parameters for wave ring (ring number,module ID, INST_ID,HeartBeats frequency [s], Debug Flag); see parameters for PyEW.EWModule')
    parser.add_argument('--PGA2FinDer', '-pga2fd', metavar='PGA2FinDer', type=int, default=0, help='Flag to enable PGA-to-FinDer processing. Default 0')
    parser.add_argument('--finderConfig', '-fdCfg', metavar='finderConfig', type=str, default=None, help='Path to the STOMP/FinDer config file')
    parser.add_argument('--nChPGASmooth', '-pgaNchSm', metavar='nChPGASmooth', type=int, default=0, help='Number of neighbor channels used for PGA smoothing, split evenly left/right around the center channel. Must be even. Default 0')
    parser.add_argument('--PGAOutput', '-pgaOut', metavar='PGAOutput', type=str, default=None, help='Optional output directory for PGA debug CSV and FinDer message text files')
    parser.add_argument('--debug', '-dbg', metavar='debug', type=int, default=0, help='Debug flag for asyncio module')
    args = parser.parse_args()
    main()
