#!/usr/bin/env python
import logging
from math import pow
from itertools import izip

from netCDF4 import Dataset
import numpy as N

sigmet_log = logging

class FileFormatError(Exception):
    pass

def sigmet16_to_float(short_data):
    return convLUT[short_data]

def makeLUT():
    shortData = N.arange(2**16, dtype=N.int32)
    mask = (shortData & 0xF000) > 0
    man = shortData[mask] & 0x7FF
    exp = (shortData[mask] >> 12) & 0x00F
    mask2 = (shortData[mask] & 0x0800) > 0
    man[mask2] = man[mask2] | 0xFFFFF000
    man[~mask2] = man[~mask2] | 0x00000800
    convLUT[mask] = man.astype(N.float32) * (1 << exp) / 3.3554432E7
    convLUT[~mask] = (shortData[~mask] << 20) / 1.759218603E13

convLUT = N.zeros((2**16,), dtype=N.float32)
makeLUT()

def sigmet_to_netcdf(filename, data):
    nc = Dataset(filename, 'w')
    sigmet_log.debug('Opened file %s for writing.', filename)
    
    #Set Global Attributes
    nc.SiteID = data['SiteID']
    sigmet_log.debug('Set global attribute SiteID to: %s', data['SiteID'])
    
    #Create Dimensions
    nc.createDimension('Range', len(data['BinRanges']))
    nc.createDimension('Pulse', len(data['Azimuth']))
    sigmet_log.debug('Created nc dimensions.')
    
    #Create variable to go along with range
    var = nc.createVariable('Range', 'f', ('Range',))
    var.Units = 'm'
    var[:] = data['BinRanges']
    sigmet_log.debug('Set range values.')
    
    #Single-valued variables
    sing_vars = [('Wavelength', 'm'), ('PulseWidth', 's'),
                 ('NoisePowerH', 'dBm'), ('NoisePowerV', 'dBm')]
    for varname, units in sing_vars:
        var = nc.createVariable(varname, 'f', ())
        var.Units = units
        var[:] = data.get(varname)
    sigmet_log.debug('Set single valued variables.')

    #Variables that are a function of pulse
    pulse_vars = [('Azimuth', 'deg'), ('Elevation', 'deg'),
                  ('Time', 's'), ('TimeMSec', 'msec'), ('PRT', 's')]
    for varname, units in pulse_vars:
        var = nc.createVariable(varname, 'f', ('Pulse',))
        var.Units = units
        var[:] = data.get(varname)

    #Time series data - horizontal first
    dims = ('Pulse', 'Range')
    iq = N.array(data['IQDataH'])
    var = nc.createVariable('I_Horizontal', 'f', dims)
    var.Units = 'mW'
    var[:] = iq.real
    var = nc.createVariable('Q_Horizontal', 'f', dims)
    var.Units = 'mW'
    var[:] = iq.imag
    
    #Now vertical
    iq = N.array(data['IQDataV'])
    var = nc.createVariable('I_Vertical', 'f', dims)
    var.Units = 'mW'
    var[:] = iq.real
    var = nc.createVariable('Q_Vertical', 'f', dims)
    var.Units = 'mW'
    var[:] = iq.imag

    nc.close()

def read_iq_data(datafile):
    data = dict()
    #File should start with this string
    header = datafile.readline()
    if header.strip()!= 'rvp8PulseInfo start':
        sigmet_log.warning('rvp8PulseInfo not found. Continuing....')
    else:
        #Get all the information from the rvp8PulseInfo structure
        #While loop used so we don't mix reads with iterators
        sigmet_log.debug('rvp8PulseInfo start found.')
        while True:
            line = datafile.readline()
            if line.strip() == 'rvp8PulseInfo end':
                sigmet_log.debug('rvp8PulseInfo end found.')
                break
            else:
                attrib, val = line.split('=')
                data[attrib] = val.strip()
                sigmet_log.debug('Info Attribute: %s->%s', attrib, val)
    
    data['pulse_data'] = list()
    
    #We should now have a sequence of pulses
    while True:
        header = datafile.readline()
        if header.strip()!= 'rvp8PulseHdr start':
            if not header:
                #We hit the end of the file
                sigmet_log.debug('End of file detected.')
                break
            else:
                sigmet_log.error('rvp8PulseHdr not found!')
                raise FileFormatError('rvp8PulseHdr not found!')
        else:
            #Get all the information from the rvp8PulseHdr structure
            #While loop used so we don't mix reads with iterators
            sigmet_log.debug('rvp8PulseHdr start found')
            pulse = dict()
            while True:
                line = datafile.readline()
                if line.strip() == 'rvp8PulseHdr end':
                    sigmet_log.debug('rvp8PulseHdr end found.')
                    break
                else:
                    attrib, val = line.split('=')
                    pulse[attrib] = val.strip()
                    sigmet_log.debug('Hdr Attribute: %s->%s', attrib, val)
            #Now read in the binary data - size is 2 * iNumVecs * iVIQPerBin
            #2 for I+Q
            num_bins = int(pulse['iNumVecs'])
            num_pol = int(pulse['iVIQPerBin'])

            #Int16 isn't correct, but will at least give 16-bit chunks
            count = 2 * num_bins * num_pol
            buf = datafile.read(2 * count)
            pulse['iqdata'] = N.frombuffer(buf, dtype=N.int16, count=count)
            data['pulse_data'].append(pulse)
            sigmet_log.debug('Finished decoding HighSNR data.')

    return data

def convert_sigmet(data):
    ndat = dict()
    ndat['SiteID'] = data['sSiteName'][:-2] #Takes TS off of name
    ndat['Wavelength'] = float(data['fWavelengthCM']) / 100.0 # cm->m
    ndat['PulseWidth'] = float(data['fPWidthUSec']) * 1e-6 # usec->sec

    npwr_H, npwr_V = data['fNoiseDBm'].split()
    ndat['NoisePowerH'] = float(npwr_H)
    ndat['NoisePowerV'] = float(npwr_V)
    
    #The range mask is a sequence of 16-bit integers which work as bitmask,
    #showing which bins at resolution range_res were actually generated.
    range_res = float(data['fRangeMaskRes'])
    range_mask = [int(i) for i in data['iRangeMask'].split()]
    ranges = []
    cur_rng = range_res
    for mask in range_mask:
        for i in xrange(16):
            if mask & (1<<i):
                ranges.append(cur_rng)
            cur_rng += range_res
    
    sys_clock = float(data['fSyClkMHz']) * 1.0e6 #Get system clock in Hz
    max_power = pow(10.0, float(data['fSaturationDBM']) / 10.0) / 2 # dBm -> mW
    for pulse in data['pulse_data']:
        #Angles are given in a binary angle format that maps 0-360 to 16-bit int
        ndat.setdefault('Azimuth', []).append(int(pulse['iAz']) * 360./65535)
        ndat.setdefault('Elevation', []).append(int(pulse['iEl']) * 360./65535)

        #Time is in epoch time, MSec is the number of milliseconds beyond that
        ndat.setdefault('Time', []).append(int(pulse['iTimeUTC']))
        ndat.setdefault('TimeMSec', []).append(int(pulse['iMSecUTC']))
        
        #Convert to PRT in seconds from the number of system clock ticks
        ndat.setdefault('PRT', []).append(int(pulse['iPrevPRT']) / sys_clock)
        
        #Covert packed format to floating point values and then to complex
        iqdata = sigmet16_to_float(pulse['iqdata']).astype(N.complex64)
        iq_comb = iqdata[::2] * max_power + 1.0j * iqdata[1::2] * max_power
        
        #Split into horizontal and vertical channels -- skip first sample, as it
        #is actually the sample of the burst pulse itself
        num_bins = int(pulse['iNumVecs'])
        ndat.setdefault('IQDataH', []).append(iq_comb[1:num_bins])
        ndat.setdefault('IQDataV', []).append(iq_comb[num_bins + 1:])
    # Make sure the number of bin ranges matches the amount of data
    ndat['BinRanges'] = ranges[-num_bins + 1:]
    
    return ndat

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print "Usage:\nsigmetiq <filename>"
        sys.exit(-1)
    
    sigmet_log.basicConfig(level=logging.WARN)
    
    filename = sys.argv[1]
    sigmet_log.info('Reading data from file: %s', filename)
    data = read_iq_data(file(filename, 'rb'))
    sigmet_log.info('Done reading.')

    sigmet_log.info('Parsing...')
    parsed_data = convert_sigmet(data)
    sigmet_log.info('Parsing done.')

    sigmet_log.info('Writing to netcdf file: %s', filename + '.nc')
    sigmet_to_netcdf(filename + '.nc', parsed_data)
    sigmet_log.info('Done writing.')
