# Module containing real-time processing classes and functions
import numpy as np
from obspy.core.trace import Trace 
from obspy.core.trace import Stats
from obspy.core.utcdatetime import UTCDateTime
import os

# Necessary for multi-threading 
import asyncio

# Necessary to run PhaseNet-DAS in real time
import threading
import dateutil.parser
from datetime import datetime, timezone
import atexit
from dataclasses import dataclass
import pandas as pd
# DAS utilities related to picking process
import DAS_ML


############################################################################################################
# DAS RingBuffer
############################################################################################################

time_format = "%Y-%m-%dT%H%M%SZ"

class RingBuffer:
    """ class that implements a not-yet-full buffer """
    def __init__(self, buff_size, good_ch):
        """
        Input:
        buff_size [int]: maximum number of time samples in the ring buffer  
        """
        if buff_size < 0:
            raise ValueError("buff_size must be a positive integer")
        self.max = buff_size
        self.data = []
        self.good_ch = good_ch
        self.channels_info = None
        self.timeStamps = [] # Rolling buffer of timestamp axis
        
    def append(self, x, timestamps=None):
        """append an element at the end of the buffer"""
        if x.ndim == 2:
            self.data.extend([np.expand_dims(row, axis=1) for row in x[:, self.good_ch]])
        else:
            self.data.append(np.expand_dims(x[self.good_ch], axis=1))
        if timestamps is not None:
            self.timeStamps += list(timestamps) if x.ndim == 2 else [timestamps]
        if len(self.data) >= self.max:
            self.cur = 0
            self.max = len(self.data) # To allow for non-divisible data packets
            # Permanently change self's class from non-full to full
            self.__class__ = self.__Full

    def getData(self):
        """ Return array of elements from the oldest to the newest. """
        return np.concatenate(self.data, axis=1)
    
    def getTimeStamps(self):
        """ Return array of related time stamps """
        return np.array(self.timeStamps)
    
    def setObspyTraceHeader(self, inventory=None):
        """Method to set channel header from inventory.

        The XML station labels carry a reduced-subset channel index in the form
        `<label>/<reduced_subset_index>`. Only the reduced-subset index is used
        here: it is a positional index into `self.good_ch`, not a raw DAS
        channel id.
        """
        if inventory is None:
            self.stats = None
            self.channels_info = None
            self.chIds = None
            self.rawChIds = None
            self.statNames = None
            self.latitudes = None
            self.longitudes = None
            return

        station_dict = inventory.get_contents()
        xml_channels_info = station_dict['channels']
        network_codes = [stat.split(".")[0] for stat in xml_channels_info]
        station_codes = [stat.split(".")[1] for stat in xml_channels_info]
        channel_codes = [stat.split(".")[3] for stat in xml_channels_info]
        xml_subset_idx = [int(stat.split(" ")[-1][:-1].split("/")[1]) for stat in station_dict['stations']]
        good_ch = np.asarray(self.good_ch, dtype=int)

        self.stats = []
        self.statNames = [] # Necessary for streaming traveltime picks
        self.latitudes = []
        self.longitudes = []
        kept_channels_info = []
        kept_subset_idx = []
        kept_raw_ch_ids = []
        skipped_subset_idx = []

        for idx, subset_idx in enumerate(xml_subset_idx):
            subset_idx = int(subset_idx)
            if subset_idx < 0 or subset_idx >= len(good_ch):
                skipped_subset_idx.append(subset_idx)
                continue

            stat = Stats()
            stat.network = network_codes[idx]
            stat.station = station_codes[idx]
            stat.channel = channel_codes[idx]
            channel_meta = inventory.select(
                network=network_codes[idx],
                station=station_codes[idx],
                channel=channel_codes[idx],
            )
            latitude = None
            longitude = None
            if len(channel_meta.networks) > 0 and len(channel_meta.networks[0].stations) > 0:
                station_meta = channel_meta.networks[0].stations[0]
                if len(station_meta.channels) > 0:
                    chan = station_meta.channels[0]
                    latitude = getattr(chan, 'latitude', None)
                    longitude = getattr(chan, 'longitude', None)
                if latitude is None:
                    latitude = getattr(station_meta, 'latitude', None)
                if longitude is None:
                    longitude = getattr(station_meta, 'longitude', None)
            stat.latitude = latitude
            stat.longitude = longitude

            self.stats.append(stat)
            self.latitudes.append(latitude)
            self.longitudes.append(longitude)
            self.statNames.append("%s.%s.%s.--" % (network_codes[idx], station_codes[idx], channel_codes[idx]))
            kept_channels_info.append(xml_channels_info[idx])
            kept_subset_idx.append(subset_idx)
            kept_raw_ch_ids.append(int(good_ch[subset_idx]))

        self.channels_info = np.array(kept_channels_info)
        self.chIds = np.array(kept_subset_idx, dtype=int)
        self.rawChIds = np.array(kept_raw_ch_ids, dtype=int)
        self.statNames = np.array(self.statNames)
        self.latitudes = np.array(self.latitudes, dtype=object)
        self.longitudes = np.array(self.longitudes, dtype=object)

        if skipped_subset_idx:
            print(
                f'Skipping {len(skipped_subset_idx)} XML-selected channels with out-of-range reduced-subset indices: {skipped_subset_idx}',
                flush=True,
            )
        return
    
    def writeObsPyTraces(self, fs, datapath, scaling=1e6):
        """Method to write Obspy traces"""
        if self.channels_info is None:
            raise ValueError("Call method setObspyTraceHeader to set trace info before using writeObsPyTraces")
        # Getting channel data to write and timestamps
        traceData = (self.getData()[self.chIds,:]*scaling).astype(np.int32)
        timeStamps = self.getTimeStamps()
        # Loop for writing traces
        for idx in range(len(self.chIds)):
            self.stats[idx].sampling_rate = fs
            self.stats[idx].npts = traceData.shape[1]
            self.stats[idx].starttime = UTCDateTime(timeStamps[0])
            self.stats[idx].delta = 1.0/fs
            tr = Trace(traceData[idx,:], header=self.stats[idx])
            file_path = "%s/%s_%s.mseed"%(datapath,self.channels_info[idx], timeStamps[0].strftime(time_format))
            tr.write(file_path, format="MSEED", reclen=512, encoding='STEIM2')
            os.chmod(file_path, 0o644)
        return
    
    def send2ew(self, fs, waveMod, ringID=0, scaling=1e6):
        """Function to send buffered data to Earthworm wavering"""
        if scaling > 1.0:
            traceData = (self.getData()[self.chIds,:]*scaling).astype(np.int32)
            dataFormat = 'i4'
        else:
            # Not currently tested
            traceData = self.getData()[self.chIds,:].astype(np.float32)
            dataFormat = 'f4'
        npts = traceData.shape[1]
        timeStamps = self.getTimeStamps()
        # Ensure the timestamp is in UTC
        timestamp_init = timeStamps[0].astimezone(timezone.utc).timestamp()
        print(timeStamps, timestamp_init)
        for idx in range(len(self.chIds)):
            wave = {
                'station': self.stats[idx].station, 
                'network': self.stats[idx].network, 
                'channel': self.stats[idx].channel, 
                'location': '--', 
                'nsamp': npts, 
                'samprate': fs, 
                'startt': timestamp_init,
                'endt': timestamp_init+(npts-1)/fs,
                'datatype': dataFormat,
                'data': traceData[idx,:]
            }
            waveMod.put_wave(ringID, wave)
        return
    
    class __Full:
            """ class that implements a full buffer """
            def append(self, x, timestamps=None):
                """Append an element overwriting the oldest one."""
                if x.ndim == 2:
                    ntimes = x.shape[0]
                    for idx in range(ntimes):
                        self.data[self.cur] = np.expand_dims(x[idx, self.good_ch], axis=1)
                        if timestamps is not None:
                            self.timeStamps[self.cur] = timestamps[idx]
                        self.cur = (self.cur + 1) % self.max
                else:
                    self.data[self.cur] = np.expand_dims(x[self.good_ch], axis=1)
                    if timestamps is not None:
                        self.timeStamps[self.cur] = timestamps
                    # Update pointer
                    self.cur = (self.cur + 1) % self.max
            
            def getData(self):
                """ Return array of elements in correct order """
                return np.concatenate(self.data[self.cur:]+self.data[:self.cur], axis=1)
            
            def getTimeStamps(self):
                """ Return array of timestamps in correct order """
                return np.array(self.timeStamps[self.cur:]+self.timeStamps[:self.cur])
            
            def writeObsPyTraces(self, fs, datapath, scaling=1e6):
                """Method to write Obspy traces"""
                if self.channels_info is None:
                    raise ValueError("Call method setObspyTraceHeader to set trace info before using writeObsPyTraces")
                # Getting channel data to write and timestamps
                traceData = (self.getData()[self.chIds,:]*scaling).astype(np.int32)
                timeStamps = self.getTimeStamps()
                # Loop for writing traces
                for idx in range(len(self.chIds)):
                    self.stats[idx].sampling_rate = fs
                    self.stats[idx].npts = traceData.shape[1]
                    self.stats[idx].starttime = UTCDateTime(timeStamps[0])
                    self.stats[idx].delta = 1.0/fs
                    tr = Trace(traceData[idx,:], header=self.stats[idx])
                    file_path = "%s/%s_%s.mseed"%(datapath,self.channels_info[idx], timeStamps[0].strftime(time_format))
                    tr.write(file_path, format="MSEED", reclen=512, encoding='STEIM2')
                    os.chmod(file_path, 0o644)
                return
            
            def send2ew(self, fs, waveMod, ringID=0, scaling=1e6):
                """
                    Sends buffered seismic data to Earthworm wavering.

                    Parameters:
                        - fs (float): Sampling frequency of the data.
                        - waveMod (object): Earthworm wave module instance to send data to.
                        - ringID (int, optional): Identifier for the Earthworm ring. Defaults to 0.
                        - scaling (float, optional): Scaling factor for the data. Defaults to 1e6.

                    Returns:
                        - None
                """
                if scaling > 1.0:
                    traceData = (self.getData()[self.chIds,:]*scaling).astype(np.int32)
                    dataFormat = 'i4'
                else:
                    # Not currently tested
                    traceData = self.getData()[self.chIds,:].astype(np.float32)
                    dataFormat = 'f4'
                npts = traceData.shape[1]
                timeStamps = self.getTimeStamps()
                timestamp_init = timeStamps[0].timestamp()
                for idx in range(len(self.chIds)):
                    wave = {
                        'station': self.stats[idx].station, 
                        'network': self.stats[idx].network, 
                        'channel': self.stats[idx].channel, 
                        'location': '--', 
                        'nsamp': npts, 
                        'samprate': fs, 
                        'startt': timestamp_init,
                        'endt': timestamp_init+(npts-1)/fs,
                        'datatype': dataFormat,
                        'data': traceData[idx,:]
                    }
                    waveMod.put_wave(ringID, wave)
                return


############################################################################################################
# PhaseNet-DAS realtime functions and utilities
############################################################################################################

# Function to run the event loop in a separate thread
def start_event_loop(loop):
    asyncio.set_event_loop(loop)
    loop.run_forever()

# Function to stop correctly thread at the exit
def close_picking_thread(thread, loop):
    loop.call_soon_threadsafe(loop.stop)
    thread.join()
    return

def start_picking_thread(device="cuda"):
    """Function to start picking thread and loading ML model on proper device"""
    # Create an event loop for the separate thread
    loop = asyncio.new_event_loop()

    # Start the event loop in a new thread
    thread = threading.Thread(target=start_event_loop, args=(loop,))
    thread.start()

    # Registering the closure of the loop once the main program has ended
    atexit.register(lambda: close_picking_thread(thread, loop))

    # Loading ML model on device
    print("Loading ML model...", end=" ", flush=True)
    
    # Make sure `DAS_ML.preload_model_async` is a coroutine
    task1 = asyncio.run_coroutine_threadsafe(DAS_ML.preload_model_async(device=device), loop)
    
    # Wait for task1 to complete
    try:
        task1.result()  # This will block until task1 completes
        print("DONE", flush=True)
    except Exception as e:
        print(f"Failed to load ML model: {e}", flush=True)
        # Ensure the thread and loop are closed if there's an error
        close_picking_thread(thread, loop)
        raise

    return loop

async def real_time_picking_async(DASdata, dt, timeStamps, minbuf=2.0, maxbuf=68.0):
    """Function performing real-time"""
    time_format_picking = "%Y-%m-%dT%H:%M:%S.%f+00:00" 
    fs = 1.0/dt
    first_timestamp = dateutil.parser.parse(timeStamps[0].strftime(time_format_picking))
    TT_picks = DAS_ML.phasenet_das(DASdata, first_timestamp.strftime(time_format_picking), 0, dt)
    if len(TT_picks) == 0:
        TT_picks = None
    else:
        TT_picks["station_id"] = TT_picks["station_id"].apply(lambda x: int(x))
        TT_picks['phase_time'] = pd.to_datetime(TT_picks['phase_time'])
        TT_picks["phase_time_seconds"] = TT_picks["phase_time"].apply(lambda x: (x - first_timestamp).total_seconds())
        TT_picks = TT_picks[(TT_picks['phase_time_seconds'] >= minbuf) & (TT_picks['phase_time_seconds'] <= maxbuf)]
        if len(TT_picks) == 0:
            return None
        TT_picks["begin_time"] = first_timestamp
        TT_picks["peak Strain rate [nm/m/s]"] = DAS_ML.extract_peak_amp(DASdata, TT_picks["phase_time_seconds"].to_numpy(), 
                                  fs , 0.0, 0.0, 25.0, chIDs=TT_picks["station_id"].to_numpy())
    return TT_picks

Pickcounter = 0
def merge_stream_picks(TT_picksBuf, TT_picksNew, delta_t_thres=2.0, maxBuf=3600.0, pickRing=None, streamCh=None, chCodes=None, minPhaseScore=0.6, qualityThresHold=[1.0, 0.97, 0.9, 0.0], quality_values=[0, 1, 2, 3]):
    """Function to buffer traveltime picks and stream them using a pyEarthworm pick ring"""
    global Pickcounter  # Declare Pickcounter as global inside the function

    # Get current time in UTC
    curTimeUTC = datetime.now(timezone.utc)

    if TT_picksBuf is None and TT_picksNew is None:
        return None

    # Ensure 'phase_time' is a datetime object in both buffers
    if TT_picksBuf is not None:
        TT_picksBuf = TT_picksBuf.copy()
        TT_picksBuf["phase_time"] = pd.to_datetime(TT_picksBuf["phase_time"])

        # Remove old picks beyond maxBuf
        TT_picksBuf["phase_time_seconds"] = TT_picksBuf["phase_time"].apply(lambda x: (curTimeUTC - x).total_seconds())
        TT_picksBuf = TT_picksBuf[TT_picksBuf["phase_time_seconds"] < maxBuf]
        
        if len(TT_picksBuf) == 0:
            return None

    if TT_picksNew is not None:
        TT_picksNew = TT_picksNew.copy()
        TT_picksNew["phase_time"] = pd.to_datetime(TT_picksNew["phase_time"])

    # Combine the buffers
    if TT_picksBuf is not None and TT_picksNew is not None:
        TT_picks = pd.concat([TT_picksBuf, TT_picksNew])
    elif TT_picksBuf is None and TT_picksNew is not None:
        TT_picks = TT_picksNew
    else:
        return TT_picksBuf  # Nothing new, return the buffer

    # Sort by station_id, phase_type, and phase_time
    TT_picks = TT_picks.sort_values(by=['station_id', 'phase_type', 'phase_time'])

    # Compute time difference between consecutive picks within the same station and phase
    TT_picks['time_diff'] = TT_picks.groupby(['station_id', 'phase_type'])['phase_time'].diff().dt.total_seconds().abs()

    # Keep only new picks that are sufficiently spaced apart
    TT_picks = TT_picks[(TT_picks['time_diff'] > delta_t_thres) | (TT_picks['time_diff'].isna())]

    # Select new picks for streaming
    TT_picksNew = TT_picks[(TT_picks['time_diff'] > delta_t_thres) | (TT_picks['time_diff'].isna())]

    if TT_picksBuf is not None:
        # Find picks in TT_picksNew that are NOT in TT_picksBuf
        TT_picksNew = TT_picksNew.merge(TT_picksBuf[['station_id', 'phase_type', 'phase_time']], 
                                        on=['station_id', 'phase_type', 'phase_time'], 
                                        how='left', 
                                        indicator=True)
        # Keep only the new picks
        TT_picksNew = TT_picksNew[TT_picksNew['_merge'] == 'left_only'].drop(columns=['_merge'])
    
    if len(TT_picksNew) > 0:
        # Ensure 'phase_score' is numeric
        TT_picksNew.loc[:, 'phase_score'] = pd.to_numeric(TT_picksNew['phase_score'], errors='coerce')

        # Filter picks based on phase score threshold
        TT_picksNew = TT_picksNew[TT_picksNew["phase_score"] > minPhaseScore]

        # Assign quality values based on thresholds
        conditions = [
            (TT_picksNew['phase_score'] >= qualityThresHold[0]),
            (TT_picksNew['phase_score'] >= qualityThresHold[1]) & (TT_picksNew['phase_score'] < qualityThresHold[0]),
            (TT_picksNew['phase_score'] >= qualityThresHold[2]) & (TT_picksNew['phase_score'] < qualityThresHold[1]),
            (TT_picksNew['phase_score'] < qualityThresHold[2])
        ]
        TT_picksNew = TT_picksNew.copy()  # Explicitly create a copy
        TT_picksNew.loc[:, 'Q'] = np.select(conditions, quality_values)

        # Stream picks through pickRing if available
        if pickRing is not None and streamCh is not None:
            # Filter only picks from the selected stream channels
            TT_picksNew = TT_picksNew[TT_picksNew["station_id"].isin(streamCh)]
            if len(TT_picksNew) > 0:
                # Add picks to pickRing
                for _, pick in TT_picksNew.iterrows():
                    chidx = np.where(pick["station_id"] == streamCh)[0]
                    picktime = pick["phase_time"].strftime("%Y-%m-%dT%H%M%S.%f+00:00")[:-6].replace("-", "").replace("T", "").replace(":", "")
                    pickString = "8 99 4 %s %s ?%s %s 0 0 0" % (Pickcounter, chCodes[chidx][0], pick['Q'], picktime)
                    pickRing.put_msg(1, 8, pickString)
                    print(pickString)
                    Pickcounter += 1
        else:
            print("Cannot stream picks through pick ring without a running pickRing and a provided streamCh")

    # Drop the 'time_diff' column as it's no longer needed
    TT_picks = TT_picks.drop(columns=['time_diff'])
    
    return TT_picks


############################################################################################################
# PGA2FinDer realtime functions and utilities
############################################################################################################


@dataclass
class PGAState:
    conversion_df: pd.DataFrame
    buffer_channels: np.ndarray
    channels_info: np.ndarray
    xml_subset_idx: np.ndarray
    raw_channel_ids: np.ndarray
    stat_names: np.ndarray
    export_row_idx: np.ndarray
    n_ch_smooth: int
    conversion_factors: np.ndarray
    channel_to_factor: dict
    finder_station_template: dict
    finder_metadata: dict


def init_pga_state(conversion_df, buffer_channels, channels_info, xml_subset_idx, raw_channel_ids, stat_names, export_row_idx, stats, n_ch_smooth):
    """Initialize PGA conversion metadata for the current ring-buffer channels and XML export subset.

    `xml_subset_idx` and `export_row_idx` are positional indices into the reduced
    channel subset defined by `buffer_channels`. They are not raw DAS channel ids,
    and the leading XML label component is not used for indexing.
    """
    if channels_info is None or xml_subset_idx is None or raw_channel_ids is None or stat_names is None or export_row_idx is None or stats is None:
        raise ValueError('PGA processing requires XML-derived export metadata')
    if conversion_df is None:
        raise ValueError('PGA processing requires a conversion dataframe')

    buffer_channels = np.asarray(buffer_channels, dtype=int)
    export_row_idx = np.asarray(export_row_idx, dtype=int)
    if export_row_idx.size == 0:
        raise ValueError('No XML-selected channels overlap the input channel list')

    conversion_df = conversion_df.copy()
    conversion_df['Channel'] = conversion_df['Channel'].astype(int)
    conversion_df['PGA/PSR-Ratio'] = pd.to_numeric(
        conversion_df['PGA/PSR-Ratio'], errors='coerce'
    )

    if conversion_df['PGA/PSR-Ratio'].isna().any():
        bad_rows = conversion_df[conversion_df['PGA/PSR-Ratio'].isna()]['Channel'].tolist()
        raise ValueError(
            'Invalid PGA/PSR-Ratio values for channels: %s'
            % ', '.join(map(str, bad_rows))
        )

    channel_to_factor = dict(zip(conversion_df['Channel'], conversion_df['PGA/PSR-Ratio']))
    missing_buffer_channels = [int(ch) for ch in buffer_channels if int(ch) not in channel_to_factor]
    if missing_buffer_channels:
        raise ValueError(
            'Missing PGA/PSR-Ratio for input channels: %s'
            % ', '.join(map(str, missing_buffer_channels[:20]))
        )

    conversion_factors = np.array([channel_to_factor[int(ch)] for ch in buffer_channels], dtype=float)
    finder_station_template = {}
    for idx, sncl in enumerate(stat_names):
        stat = stats[idx]
        latitude = getattr(stat, 'latitude', None)
        longitude = getattr(stat, 'longitude', None)
        finder_station_template[sncl] = {
            'lat': f"{float(latitude):.3f}" if latitude is not None else None,
            'lon': f"{float(longitude):.3f}" if longitude is not None else None,
        }

    return PGAState(
        conversion_df=conversion_df,
        buffer_channels=buffer_channels,
        channels_info=np.array(channels_info),
        xml_subset_idx=np.array(xml_subset_idx, dtype=int),
        raw_channel_ids=np.array(raw_channel_ids, dtype=int),
        stat_names=np.array(stat_names),
        export_row_idx=export_row_idx,
        n_ch_smooth=int(n_ch_smooth),
        conversion_factors=conversion_factors,
        channel_to_factor=channel_to_factor,
        finder_station_template=finder_station_template,
        finder_metadata={'columns': ['lat', 'lon', 'sncl', 'timestamp', 'HSZ', 'HS1', 'HS2']},
    )


def compute_peak_strain_rate_window(data, time_stamps):
    """Compute peak absolute strain rate and the corresponding time per channel."""
    if data.ndim != 2:
        raise ValueError('data must be a 2D array with shape [nch, nt]')
    if len(time_stamps) != data.shape[1]:
        raise ValueError('time_stamps length must match data.shape[1]')

    peak_indices = np.argmax(np.abs(data), axis=1)
    peak_strain_rate = np.abs(data[np.arange(data.shape[0]), peak_indices])
    peak_times = np.asarray(time_stamps)[peak_indices]
    return peak_strain_rate, peak_indices, peak_times


def convert_peak_strain_rate_to_pga(peak_strain_rate, peak_indices, time_stamps, conversion_factors):
    """Convert peak strain rate values to pseudo-PGA values using channel factors."""
    peak_strain_rate = np.asarray(peak_strain_rate, dtype=float)
    peak_indices = np.asarray(peak_indices, dtype=int)
    conversion_factors = np.asarray(conversion_factors, dtype=float)
    if peak_strain_rate.shape != conversion_factors.shape:
        raise ValueError('peak_strain_rate and conversion_factors must have the same shape')

    pga_values = peak_strain_rate * conversion_factors * 100.0
    peak_times = np.asarray(time_stamps)[peak_indices]
    return pga_values, peak_times


def smooth_pga_values_median_subset(pga_values, pga_times, export_row_idx, n_ch_smooth):
    """Apply spatial median smoothing only for the requested export channels."""
    pga_values = np.asarray(pga_values, dtype=float)
    pga_times = np.asarray(pga_times)
    export_row_idx = np.asarray(export_row_idx, dtype=int)
    if pga_values.ndim != 1:
        raise ValueError('pga_values must be a 1D array')
    if len(pga_values) != len(pga_times):
        raise ValueError('pga_values and pga_times must have the same length')
    if n_ch_smooth < 0 or n_ch_smooth % 2 != 0:
        raise ValueError('n_ch_smooth must be a non-negative even integer')
    if export_row_idx.ndim != 1:
        raise ValueError('export_row_idx must be a 1D array')

    half_window = n_ch_smooth // 2
    smoothed_values = np.empty(export_row_idx.shape[0], dtype=float)
    smoothed_times = np.empty(export_row_idx.shape[0], dtype=object)

    for out_idx, ich in enumerate(export_row_idx):
        left = max(0, ich - half_window)
        right = min(len(pga_values), ich + half_window + 1)
        window_values = pga_values[left:right]
        window_times = pga_times[left:right]
        order = np.argsort(window_values, kind='stable')
        median_idx = order[len(order) // 2]
        smoothed_values[out_idx] = window_values[median_idx]
        smoothed_times[out_idx] = window_times[median_idx]

    return smoothed_values, smoothed_times


def update_real_time_pga(ringbuff, pga_state):
    """Compute export-ready PGA rows from the current ring buffer window."""
    if ringbuff.channels_info is None:
        raise ValueError('PGA processing requires ring buffer metadata from XML')

    data = ringbuff.getData()
    time_stamps = ringbuff.getTimeStamps()
    peak_strain_rate, peak_indices, peak_times = compute_peak_strain_rate_window(data, time_stamps)
    pga_values_raw, pga_times_raw = convert_peak_strain_rate_to_pga(
        peak_strain_rate,
        peak_indices,
        time_stamps,
        pga_state.conversion_factors,
    )
    pga_values, pga_times = smooth_pga_values_median_subset(
        pga_values_raw,
        pga_times_raw,
        pga_state.export_row_idx,
        pga_state.n_ch_smooth,
    )

    rows = []
    for idx, channel_info in enumerate(pga_state.channels_info):
        export_idx = int(pga_state.export_row_idx[idx])
        network_code, station_code, _, channel_code = channel_info.split('.')
        stats = ringbuff.stats[idx]
        rows.append({
            'xml_subset_index': int(pga_state.xml_subset_idx[idx]),
            'channel_index': int(pga_state.xml_subset_idx[idx]),
            'csv_row_index': export_idx,
            'raw_channel_id': int(pga_state.raw_channel_ids[idx]),
            'channel_info': channel_info,
            'sncl': pga_state.stat_names[idx],
            'network': network_code,
            'station': station_code,
            'channel': channel_code,
            'location': '--',
            'latitude': getattr(stats, 'latitude', None),
            'longitude': getattr(stats, 'longitude', None),
            'peak_strain_rate': float(peak_strain_rate[export_idx]),
            'peak_strain_rate_time': peak_times[export_idx],
            'pga': float(pga_values[idx]),
            'pga_time': pga_times[idx],
            'HSZ': 0.0,
            'HS1': float(pga_values[idx]),
            'HS2': 0.0,
        })

    return pd.DataFrame(rows)
