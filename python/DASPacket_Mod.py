import sys
from datetime import datetime, timedelta, timezone
import numpy as np
import json, struct

def error_exit(message, code=1):
    sys.stderr.write("Error: {0}\n".format(str(message)))
    sys.exit(code)

class StreamReader:
    """Template class for any Stream data reader"""

    def getNextPacket(self):
        raise NotImplementedError("getNextPacket must be overridden")
        
    def getPacketTimestamp(self, header):
        """Returns a timestamp from the header."""
        raise NotImplementedError("getPacketTimestamp must be overridden")

    def getNumTimeSamples(self, packet):
        raise NotImplementedError("getNumTimeSamples must be overridden")

    def getSampleCount(self, packet):
        raise NotImplementedError("getSampleCount must be overridden")

    def getNumChannel(self, packet):
        raise NotImplementedError("getNumChannel must be overridden")

    def getGaugeLengthProc(self, packet):
        raise NotImplementedError("getGaugeLengthProc must be overridden")

    def getHeader(self, packet):
        raise NotImplementedError("getHeader must be overridden")

    def getPayload(self, packet):
        raise NotImplementedError("getPayload must be overridden")

    def getPayloadRad(self, packet):
        raise NotImplementedError("getPayloadRad must be overridden")

    def getDecfactorORsamplesPerPacket(self, packet):
        raise NotImplementedError("getDecfactorORsamplesPerPacket must be overridden")

    def getFs(self, packet):
        raise NotImplementedError("getFs must be overridden")

    def getConversionFactor(self, packet):
        raise NotImplementedError("getConversionFactor must be overridden")
    

class ASN_StreamReader(StreamReader):
    """Streamer Reader for ASN data"""
    def __init__(self, socket):
        self.socket = socket
        self.enable_data_delay_diag = False
        self.last_packet_header_time = None
        self.last_packet_receive_time = None
        self.last_packet_delay_seconds = None
        self.header = None
        self.rois = []
        self.times = []
        self.data = []
        self.SampleCount = 0
        self.unpackFormat = ""
        self.readHeader(self.socket.recv().decode('utf-8'))

    def readHeader(self, text):
        header = json.loads(text)
        unpackFormat = f"{header['nChannels'] * header['nPackagesPerMessage']}{'f' if header['dataType']=='float' else 'h'}"
        rois = [] # List with roi slices. 
        channel = 0
        for roi in header['roiTable']:
            channels = (roi['roiEnd'] -  roi['roiStart']) //  roi['roiDec'] +1
            rois.append(slice(channel, channel+channels, 1))
            channel += channels
        self.header = header
        self.rois = rois
        self.unpackFormat = unpackFormat
        return 

    def getHeader(self, packet):
        return self.header
       
    def getROIs(self, packet):
        return self.rois

    def timeString(self, timestamp):
        return datetime.datetime.fromtimestamp(timestamp).strftime('%Y-%m-%d %H:%M:%S.%fZ')

    def getNextPacket(self):
        """Function to obtain next packet; for ASN unpacking everything directly"""
        packet = self.socket.recv()
        if len(packet) != self.header['bytesPerPackage']*self.header['nPackagesPerMessage'] + 8:
            print('Parameters Changed. Updated header, ROIs and format; restarting data stream')
            return b''
        # Getting timestamp
        timestamp = struct.unpack('Q' , packet[:8])[0]*1e-9 # Nanoseconds to seconds
        time_utc = datetime.fromtimestamp(timestamp, tz=timezone.utc)
        dt =  self.header["dt"]
        # Reading the data raw
        self.data = np.array(struct.unpack(self.unpackFormat, packet[8:])).reshape((self.header['nPackagesPerMessage'], self.header['nChannels']))
        # Updating internal sample count
        self.SampleCount += self.data.shape[0]
        self.times = [time_utc + timedelta(seconds=(it*dt)) for it in np.arange(self.data.shape[0])]
        if self.enable_data_delay_diag:
            self.last_packet_header_time = time_utc
            self.last_packet_receive_time = datetime.now(timezone.utc)
            self.last_packet_delay_seconds = (
                self.last_packet_receive_time - self.last_packet_header_time
            ).total_seconds()
        return packet
        
    def getPacketTimestamp(self, header):
        """Returns the timestamps of the data points"""
        return self.times

    def getNumTimeSamples(self, packet):
        return self.data.shape[0]

    def getSampleCount(self, packet):
        return self.SampleCount

    def getNumChannel(self, packet):
        return self.header['nChannels']

    def getGaugeLengthProc(self, packet):
        return self.header['gaugeLength']

    def getPayload(self, packet):
        raise NotImplementedError("getPayload must be overridden")

    def getPayloadRad(self, packet):
        return self.data

    def getDecfactorORsamplesPerPacket(self, packet):
        return self.header["nPackagesPerMessage"]

    def getFs(self, packet):
        return 1.0/self.header["dt"]

    def getConversionFactor(self, packet):
        return 1.0/(self.header["sensitivities"][0]["factor"])



class OptaSenseStreamReader(StreamReader):
    """Streamer Reader for OptaSense data"""
    def __init__(self, socket):
        self.buffer = b''
        self.socket = socket
        self.enable_data_delay_diag = False
        self.last_packet_header_time = None
        self.last_packet_receive_time = None
        self.last_packet_delay_seconds = None

    def getPacketHeaderSize(self, packet):
        return int(self.from_bytes(packet[8:9]))
    
    def getPacketPayloadSize(self, packet):
        return int(self.from_bytes(packet[12:15]))
    
    def getPacketFooterSize(self, packet):
        return int(self.from_bytes(packet[9:10]))

    def getPacketSize(self, header):
        SW0 = header[0:4]
        SW1 = header[4:8]
        if SW0 != b'\xff\xff\xff\x7f' or SW1 != b'\x00\x00\x00\x80':
            raise RuntimeError("No Sync Word found in the data. Wrong format?")
        return (self.getPacketHeaderSize(header) +
                self.getPacketPayloadSize(header) +
                self.getPacketFooterSize(header))

    def getNextPacket(self):
        BUFF_SIZE = 2048
        LEN_SIZE = 16
        
        chunks = []
        bytes_recd = len(self.buffer)
        while bytes_recd < LEN_SIZE:
            chunk = self.socket.recv(min(LEN_SIZE - bytes_recd, BUFF_SIZE))
            if chunk == b'':
                raise RuntimeError("socket connection broken")
            chunks.append(chunk)
            bytes_recd += len(chunk)
        self.buffer += b''.join(chunks)
        
        chunks = []
        PACKET_LEN = self.getPacketSize(self.buffer)
        bytes_recd = len(self.buffer)
        while bytes_recd < PACKET_LEN:
            chunk = self.socket.recv(min(PACKET_LEN - bytes_recd, BUFF_SIZE))
            if chunk == b'':
                raise RuntimeError("socket connection broken")
            chunks.append(chunk)
            bytes_recd += len(chunk)
        self.buffer += b''.join(chunks)
        
        packet = self.buffer[:PACKET_LEN]
        self.buffer = self.buffer[PACKET_LEN:]
        if self.enable_data_delay_diag:
            self.last_packet_header_time = self.getPacketTimestamp(packet)
            self.last_packet_receive_time = datetime.now(timezone.utc)
            self.last_packet_delay_seconds = (
                self.last_packet_receive_time - self.last_packet_header_time
            ).total_seconds()
        return packet
    
    def from_bytes(self, data, big_endian=False):
        if isinstance(data, str):
            data = bytearray(data)
        if big_endian:
            data = reversed(data)
        num = 0
        for offset, byte in enumerate(data):
            num += byte << (offset * 8)
        return num
    
    def getPacketTimestamp(self, header):
        """Function to obtain timestamp sample"""
        '''Returns a timestamp from the header in the following format: yyyy-MM-ddThh:mm:ss.uuuuuu.'''
        seconds = int(self.from_bytes(header[96:100])) #Seconds since 0:00 Jan 01 2000
        milliseconds = int(self.from_bytes(header[100:104]))/1000000 #Nanoseconds in the header
        base_datetime = datetime(2000, 1, 1, tzinfo=timezone.utc)
        delta = timedelta(0, seconds, 0, milliseconds)
        timestamp = (base_datetime + delta)
        return timestamp
    
    def getNumTimeSamples(self, packet):
        return int(self.from_bytes(packet[80:81]))
    
    def getSampleCount(self, packet):
        return np.int64(self.from_bytes(packet[88:96]))
    
    def getNumChannel(self, packet):
        return int(self.from_bytes(packet[76:79]))
    
    def getGaugeLengthProc(self, packet, nFiber=1.4682, c=299792458.0):
        u16ADCClock = 1.0 / float(self.from_bytes(packet[52:54]))
        return int(self.from_bytes(packet[60:62])) * c * 1e-6 * u16ADCClock / nFiber * 0.5
    
    def getTimePerPacket(self, packet):
        u16ADCClock = int(self.from_bytes(packet[52:54]))
        u32PingPeriodCsu = int(self.from_bytes(packet[44:47]))
        u16Decimation = self.getDecfactorORsamplesPerPacket(packet)
        u8NumTimeSamples = self.getNumTimeSamples(packet)
        return int(1 / (u16ADCClock / u32PingPeriodCsu / u16Decimation)) * u8NumTimeSamples
    
    def getHeader(self, packet):
        return packet[:self.getPacketHeaderSize(packet)]
    
    def getPayload(self, packet):
        header_size = self.getPacketHeaderSize(packet)
        payload_size = self.getPacketPayloadSize(packet)
        return packet[header_size:header_size + payload_size]
    
    def getFooter(self, packet):
        return packet[-self.getPacketFooterSize(packet):]
    
    def getPayloadRad(self, packet):
        dataWidth = self.getDataWidth(packet) // 8
        dtype = '<i2' if dataWidth == 2 else '<i4'
        return np.frombuffer(self.getPayload(packet), dtype=dtype)
    
    def getPhaseLSB(self, packet):
        return int(self.from_bytes(packet[62:63]))
    
    def getDataWidth(self, packet):
        return int(self.from_bytes(packet[75:76]))
    
    def getDecfactorORsamplesPerPacket(self, packet):
        return int(self.from_bytes(packet[82:84]))
    
    def getOCP(self, packet):
        return int(self.from_bytes(packet[81:82]))
    
    def getFs(self, packet):
        u16ADCClock = int(self.from_bytes(packet[52:54]))
        u32PingPeriodCsu = int(self.from_bytes(packet[44:47]))
        u16Decimation = self.getDecfactorORsamplesPerPacket(packet)
        return u16ADCClock / u32PingPeriodCsu / u16Decimation * 1000000
    
    def getConversionFactor(self, packet):
        GaugeL = self.getGaugeLengthProc(packet)
        nFiber = 1.4682
        lamdLaser = 1550.0
        eta = 0.78
        factor = 4.0 * np.pi * eta * nFiber * GaugeL / lamdLaser
        radconv = 1.0
        phaseLSB = self.getPhaseLSB(packet)
        scaleFactor = 2 * np.pi / (2 ** phaseLSB)
        return 1.0 / factor / radconv * scaleFactor