#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 11:49:55 2023

@author: apdcam
"""

import os
import sys  # For error messages
import warnings
import datetime  # For handling the timing
import calendar  # For handling the timing
import urllib.request  # For getting data from Web API (JSON)
import json  # For handling downloaded JSON files
import logging
import numpy as np
import decimal

class ArchiveSignal(object):

    """**Bases:** :class:`object`

    This class handles the data readout from the Wendelstein 7-X ArchiveDB data
    acquisition system via the Web API interface.

    """

    def __init__(self, signal_id=None, exp_id="20221201.001"):
        """ Initializes the basic signal information
        """
        if signal_id is not None:
            self.stream = stream_gen(signal_id, exp_id=exp_id)  # The address of the datastream
        else:
            self.stream=None
        self.time_query = ""  # The time for which to look for the data
        ''' Resampling. For high time resolution data, the data
        for one shot can be otherwise unreasonably large to handle:'''
        self.resample_query = ""
        self.url = ""

    @staticmethod
    def url_request(url_str):
        """
        This function handles an url request for the W7X archive system.

        :param url_str: the requested url
        :type url_str: str
        :return: the requested data
        """
        request = urllib.request.Request(url_str)
        request.add_header('Content-Type', 'application/json')

        try:
            url = urllib.request.urlopen(request)
        except Exception as e:
#            logger.error("During the URL request, the following problem occured:")
#            logger.error(e)
#            logger.error("The requested URL:")
#            logger.error(request.get_full_url())
            return None

        try:
            data = json.loads(url.read().decode(url.info().get_param('charset') or 'utf-8'))
        except Exception as e:
#            logger.error("There were some problems fetching the JSON object!")
#            logger.error(e)
#            logger.error(sys.exc_info()[0])
#            logger.error("The requested URL:")
#            logger.error(request.get_full_url())
            url.close()
            return None

        url.close()
        if data == None:
            warnings.warn('Data missing at '+url_str)

        return data

    def archive_pull(self):
        """
        This function creates the url and pulls the data from the W7X archive system.

        :stream: the signal string
        :time_query: the time query string
        :reduction: the string for the downsampling
        :return: the requested data
        """
        base = 'http://archive-webapi.ipp-hgw.mpg.de/'
        self.url = base + self.stream + self.time_query + self.resample_query
        data = self.url_request(self.url)
        return data

    def last_version(self):
        """This function gives back the latest data version for the requested data
    
        :return: an integer for the latest data version
    
        """
        base = 'http://archive-webapi.ipp-hgw.mpg.de/'
        self.url = base + self.stream + self.time_query + self.resample_query
        data = self.url_request(self.url.replace("signal", "versions"))
        return data["versionInfo"][0]["number"]

    def time_query_gen(self, start, stop):
        """This function gives back the time_query string for the requested data

        :param start: contains the start time of the signal request
        :type start: datetime.datetime object
        :param stop: contains the end time of the signal request
        :type stop: a datetime.datetime object
        :return: the time_query part of the url

        """
        try:
            start_ns = self.convert_nanosecs(start.year, start.month, start.day,
                                         start.hour, start.minute, start.second, start.microsecond)
        except AttributeError:
            start_ns = start
        self.from_time = start_ns
        try:
            stop_ns = self.convert_nanosecs(stop.year, stop.month, stop.day, stop.hour,
                                            stop.minute, stop.second, stop.microsecond)
        except AttributeError:
            stop_ns = stop
        self.time_query = "_signal.json?from=" + str(start_ns) + '&upto=' + str(stop_ns)
        return self.time_query
    
    def time_query_gen_ns(self, start_ns, stop_ns):
        """This function gives back the time_query string for the requested data if the start
        and stop time are given in ns

        :param start: contains the start time of the signal request
        :type start: int
        :param stop: contains the end time of the signal request
        :type stop: int
        :return: the time_query part of the url

        """
        self.from_time = start_ns
        self.time_query = "_signal.json?from=" + str(start_ns) + '&upto=' + str(stop_ns)
        return self.time_query

    def shotid_time_query_gen(self, shotid):
        """ This one creates the whole the time query for the data using
        the shotid. It sends a request to the archive to get the shotid
        timerange.

        :shotid: The ID of the shot given as a string.
        """
        url = "http://archive-webapi.ipp-hgw.mpg.de/programs.json?from=" + shotid
        response = self.url_request(url)
        if (response is None):
            raise ConnectionError("Cannot access W7X webapi.")
        from_time = response['programs'][0]['trigger']['1'][0]
        self.from_time = from_time
        upto_time = response['programs'][0]['trigger']['4'][0]
        self.time_query = "_signal.json?from=" + str(from_time-1000000000) + '&upto=' + str(upto_time+1000000000)
        return self.time_query

    def convert_nanosecs(self, year, month, day, hour, minute, second=0, microsecond=0):
        """Convert the given inputs to nanoseconds based on the UTC time format.

        :param year: the year of the date
        :type year: int
        :param month: the month of the date
        :type month: int
        :param day: the day of the date
        :type day: int
        :param hour: the hour of the date
        :type hour: int
        :param minute: the minute of the date
        :type minute: int
        :param second: second of the date (optional)
        :type second: int
        :param microsecond: microsecond of the date (optional)
        :type microsecond: int

        :return: the time in nanoseconds (UTC format)

        """
        d = datetime.datetime(year, month, day, hour, minute, second, microsecond)
        time_ns = decimal.Decimal(round(calendar.timegm(d.timetuple()) * 1000) + d.microsecond)*1000000
        return time_ns

    def downsamp_gen(self, sample_num):
        """ Creates the string for obtaining a server-side reducted data.
        Currently the server can do a min/max reduction:  the original data
        array will be devided into nSamples/2 sections, and the reduced array
        consists of the minimum and maximum values of each section.
        Seemingly however, this does not work perfectly at the moment due to
        problems on the W7-X network side.

        :sample_num: number of samples, integer

        :return: reduction string part of the url
        """
        self.resample_query = '&reduction=minmax&nSamples='+str(sample_num)
        return self.resample_query


def stream_gen(signal_id, exp_id="20221201.001"):
    """ Creates the string for obtaining a number of specific signals of OP1

    :signal_id: signal id string
    :return: reduction string part of the url, if there is no signal_id
    string in stream_dict, then just return signal_id
    """
    #for campaign OP1, the location of a number of signals was different
    if int(exp_id[:8]) < 20221201:
        return stream_gen_OP1(signal_id)
    elif int(exp_id[:8]) < 20231201:
        return stream_gen_OP2(signal_id)
    else:
        return stream_gen_curr(signal_id)

def stream_gen_OP1(signal_id):
    stream_dict = dict()
    # This is Pixelfy at the AEQ21 port:
    stream_dict['Pixelfly1'] = \
        'ArchiveDB/raw/W7X/CoDaStationDesc.14823/DataModuleDesc.14833_DATASTREAM/0/image/unscaled'

    # This is Pixelfly at the AEQ51 port:
    stream_dict['Pixelfly2'] = \
        'ArchiveDB/raw/W7X/CoDaStationDesc.14858/DataModuleDesc.14866_DATASTREAM/0/image/unscaled'

    # This is EDICAM at the AEQ11 port:
    stream_dict['EDICAM-AEQ11'] = \
        'ArchiveDB/raw/W7X/QSV_VideoCameras_Edicam_AEQ11/ROIP1Data_DATASTREAM/0/image/unscaled'

    # Total ECRH power:
    stream_dict['ECRH'] = \
        'ArchiveDB/codac/W7X/CBG_ECRH/TotalPower_DATASTREAM/V1/0/Ptot_ECRH/scaled/'
    # NBI source 7 beam current:
    stream_dict['NBI-7'] = \
        'ArchiveDB/codac/W7X/CoDaStationDesc.31/DataModuleDesc.24119_DATASTREAM/5/HGV_3%20Monitor%20I/scaled/'
    # NBI source 8 beam current:
    stream_dict['NBI-8'] = \
        'ArchiveDB/codac/W7X/CoDaStationDesc.31/DataModuleDesc.24119_DATASTREAM/7/HGV_4%20Monitor%20I/scaled/'

    # AUG2 source:
    stream_dict['AUG-2'] = \
        'Test/raw/W7X/QSK_CXRS/AUG2_DATASTREAM/'

    # Control coils
    stream_dict['ACM11'] = \
        'ArchiveDB/raw/W7X/CoDaStationDesc.12/DataModuleDesc.16_DATASTREAM/0/ACM11%20_%20ACG%2038%20Module%201%20upper/unscaled/'
    stream_dict['ACM19'] = \
        'ArchiveDB/raw/W7X/CoDaStationDesc.12/DataModuleDesc.16_DATASTREAM/1/ACM19%20_%20ACG%2039%20Module%201%20lower/unscaled/'
    stream_dict['ACM21'] = \
        'ArchiveDB/raw/W7X/CoDaStationDesc.12/DataModuleDesc.16_DATASTREAM/2/ACM21%20_%20ACG%2028%20Module%202%20upper/unscaled/'
    stream_dict['ACM29'] = \
        'ArchiveDB/raw/W7X/CoDaStationDesc.12/DataModuleDesc.16_DATASTREAM/3/ACM29%20_%20ACG%2029%20Module%202%20lower/unscaled/'
    stream_dict['ACM31'] = \
        'ArchiveDB/raw/W7X/CoDaStationDesc.12/DataModuleDesc.16_DATASTREAM/4/ACM31%20_%20ACG%2018%20Module%203%20upper/unscaled/'
    stream_dict['ACM39'] = \
        'ArchiveDB/raw/W7X/CoDaStationDesc.12/DataModuleDesc.16_DATASTREAM/5/ACM39%20_%20ACG%2019%20Module%203%20lower/unscaled/'
    stream_dict['ACM41'] = \
        'ArchiveDB/raw/W7X/CoDaStationDesc.12/DataModuleDesc.16_DATASTREAM/6/ACM41%20_%20ACG%2048%20Module%204%20upper/unscaled/'
    stream_dict['ACM49'] = \
        'ArchiveDB/raw/W7X/CoDaStationDesc.12/DataModuleDesc.16_DATASTREAM/7/ACM49%20_%20ACG%2049%20Module%204%20lower/unscaled/'
    stream_dict['ACM51'] = \
        'ArchiveDB/raw/W7X/CoDaStationDesc.12/DataModuleDesc.16_DATASTREAM/8/ACM51%20_%20ACG%2058%20Module%205%20upper/unscaled/'
    stream_dict['ACM59'] = \
        'ArchiveDB/raw/W7X/CoDaStationDesc.12/DataModuleDesc.16_DATASTREAM/8/ACM51%20_%20ACG%2058%20Module%205%20upper/unscaled/'
        
    if signal_id in stream_dict:
        return stream_dict[signal_id]
    elif signal_id[0:4] == 'Zeff':
        # signal format expected Zeff-v20-data
        split_out = signal_id.split('-')
        version = split_out[1]
        version = version[1:]
        stream_type = split_out[2]
        if stream_type == 'data':
            location_string = \
                'Test/raw/Minerva2/Minerva.Zeff.AEZ40-AET40.USB-Spectrometer/Zeff_mean_DATASTREAM/V'+version+"/0/Zeff_mean/"
        elif stream_type == 'param':
            location_string = \
                'Test/raw/Minerva2/Minerva.Zeff.AEZ40-AET40.USB-Spectrometer/Zeff_mean_PARLOG/V'+version+"/0/Zeff_mean/"
        elif stream_type == 'err':
            location_string = \
                'Test/raw/Minerva2/Minerva.Zeff.AEZ40-AET40.USB-Spectrometer/Zeff_std_DATASTREAM/V'+version+"/0/Zeff_std/"      
        return location_string
    elif signal_id[0:2] == 'TS':
        # Thomson scattering, expected string format e.g. TS-v20-ne_map-vol1
        # That would correspond to the version 20, volume 1 data for the
        # electron density
        split_out = signal_id.split('-')
        version = split_out[1]
        version = version[1:]
        quantity = split_out[2]
        stream_type = split_out[3]
        volume = split_out[4]
        volume = volume[3:]
        thomson_stream_type_dict = dict()
        thomson_stream_type_dict['data'] = 'DATASTREAM'
        thomson_stream_type_dict['param'] = 'PARLOG'
        thomson_dict = dict()
        thomson_dict['Te_map'] = '0/Te_map/unscaled/'
        thomson_dict['ne_map'] = '1/ne_map/unscaled/'
        thomson_dict['pe_map'] = '2/pe_map/unscaled/'
        thomson_dict['Te_mean'] = '3/Te_mean/unscaled/'
        thomson_dict['ne_mean'] = '4/ne_mean/unscaled/'
        thomson_dict['pe_mean'] = '5/pe_mean/unscaled/'
        thomson_dict['Te_low97.5'] = '6/Te_low97.5/unscaled/'
        thomson_dict['Te_high97.5'] = '7/Te_high97.5/unscaled/'
        thomson_dict['ne_low97.5'] = '8/ne_low97.5/unscaled/'
        thomson_dict['ne_high97.5'] = '9/ne_high97.5/unscaled/'
        thomson_dict['pe_low97.5'] = '10/pe_low97.5/unscaled/'
        thomson_dict['pe_high97.5'] = '11/pe_high97.5/unscaled/'
        thomson_dict['Te_s'] = '12/Te_s/unscaled/'
        thomson_dict['ne_s'] = '13/ne_s/unscaled/'
        thomson_dict['pe_s'] = '14/pe_s/unscaled/'
        thomson_dict['signal_map_1'] = '15/signal_map_1/unscaled/'
        thomson_dict['signal_map_2'] = '16/signal_map_2/unscaled/'
        thomson_dict['signal_map_3'] = '17/signal_map_3/unscaled/'
        thomson_dict['signal_map_4'] = '18/signal_map_4/unscaled/'
        thomson_dict['signal_map_5'] = '19/signal_map_5/unscaled/'
        thomson_dict['signal_mean_1'] = '20/signal_mean_1/unscaled/'
        thomson_dict['signal_mean_2'] = '21/signal_mean_2/unscaled/'
        thomson_dict['signal_mean_3'] = '22/signal_mean_3/unscaled/'
        thomson_dict['signal_mean_4'] = '23/signal_mean_4/unscaled/'
        thomson_dict['signal_mean_5'] = '24/signal_mean_5/unscaled/'
        thomson_dict['x_coord'] = 'position/x_m/'
        thomson_dict['y_coord'] = 'position/y_m/'
        thomson_dict['z_coord'] = 'position/z_m/'
        if stream_type == 'data':
            location_string = 'Test/raw/W7X/QTB_Profile/volume_' + volume + '_DATASTREAM/V' + version + \
                '/' + thomson_dict[quantity]
        elif stream_type == 'param':
            location_string = 'Test/raw/W7X/QTB_Profile/volume_' + volume + '_PARLOG/V' + version \
                + '/parms/' + thomson_dict[quantity]
        return location_string
    
    elif signal_id[0:4] == 'ABES':
        # Thomson scattering, expected string format e.g. TS-v20-ne_map-vol1
        # That would correspond to the version 20, volume 1 data for the
        # electron density
        split_out = signal_id.split('-')
        version = split_out[1]
        version = version[1:]
        if version[-1] == 'S' or version[-1] == 's':
            version = version[:-1]
            sandbox = True
        else:
            sandbox = False
        quantity = split_out[2]
        stream_type = split_out[3]
        volume = split_out[4]
        volume = volume[3:]
        abes_stream_type_dict = dict()
        abes_stream_type_dict['data'] = 'DATASTREAM'
        abes_stream_type_dict['param'] = 'PARLOG'
        abes_dict = dict()
        abes_dict['density'] = '0/density/'
        abes_dict['density_error_low'] = '1/density_error_low/'
        abes_dict['density_error_high'] = '2/density_error_high/'
        abes_dict['meas_light'] = '3/meas_light/'
        abes_dict['meas_light_error'] = '4/meas_light_error/'
        abes_dict['x_coord'] = '5/x_coord/'
        abes_dict['y_coord'] = '6/y_coord/'
        if stream_type == 'data':
            location_string = 'Test/raw/W7X/QSI/vol_' + volume + '_DATASTREAM/V' + version + \
                '/' + abes_dict[quantity]
        elif stream_type == 'param':
            location_string = 'Test/raw/W7X/QSI/vol_' + volume + '_PARLOG/V' + version \
                + '/parms/' + abes_dict[quantity]
        
        if sandbox is True:
            location_string = location_string.replace("Test", "Sandbox")
        return location_string

    else:
        # print('unkown signal_id, using it as url string directly')
        return signal_id

def stream_gen_OP2(signal_id):
    """ Creates the string for obtaining a number of specific signals of OP2

    :signal_id: signal id string
    :return: reduction string part of the url, if there is no signal_id
    string in stream_dict, then just return signal_id
    """
    stream_dict = dict()
    # Total ECRH power:
    stream_dict['ECRH'] = \
        'Test/codac/W7X/CBG_ECRH/TotalPower_DATASTREAM/V1/0/Ptot_ECRH/scaled/'
    # NBI source 7 beam current:
    stream_dict['NBI-7'] = \
        'ArchiveDB/codac/W7X/CoDaStationDesc.31/DataModuleDesc.24119_DATASTREAM/5/HGV_3%20Monitor%20I/scaled/'
    # NBI source 8 beam current:
    stream_dict['NBI-8'] = \
        'ArchiveDB/codac/W7X/CoDaStationDesc.31/DataModuleDesc.24119_DATASTREAM/7/HGV_4%20Monitor%20I/scaled/'
    # QSI-CXRS data
    stream_dict["QSI-CXRS"] = \
        "Test/raw/W7X/QSI/cxrs_DATASTREAM/0/Images/"
    
    
    if signal_id in stream_dict:
        return stream_dict[signal_id]
    elif signal_id[0:4] == 'Zeff':
        # signal format expected Zeff-v20-data
        split_out = signal_id.split('-')
        version = split_out[1]
        version = version[1:]
        stream_type = split_out[2]
        if stream_type == 'data':
            location_string = \
                'Test/raw/Minerva2/Minerva.Zeff.AEZ40-AET40.USB-Spectrometer/Zeff_mean_DATASTREAM/V'+version+"/0/Zeff_mean/"
        elif stream_type == 'param':
            location_string = \
                'Test/raw/Minerva2/Minerva.Zeff.AEZ40-AET40.USB-Spectrometer/Zeff_mean_PARLOG/V'+version+"/0/Zeff_mean/"
        elif stream_type == 'err':
            location_string = \
                'Test/raw/Minerva2/Minerva.Zeff.AEZ40-AET40.USB-Spectrometer/Zeff_std_DATASTREAM/V'+version+"/0/Zeff_std/"      
        return location_string
    elif signal_id[0:2] == 'TS':
        # Thomson scattering, expected string format e.g. TS-v20-ne
        # That would correspond to the version 20, data for the
        # electron density profile
        split_out = signal_id.split('-')
        version = split_out[1][1:]
        quantity = split_out[2]
        thomson_dict={
            "rho": "0/rho/",
            "ne": "1/ne/",
            "te": "2/te/",
            "ne_std": "3/ne_std/",
            "te_std": "4/te_std/",
            "R": "5/R/",
            "Phi": "6/Phi/",
            "Z": "7/Z/"}
        location_string = f"Test/raw/Minerva2/Minerva.ThomsonScattering.Profiles/Preliminary_Fast_DATASTREAM/V{version}" + \
                '/' + thomson_dict[quantity]
        return location_string
    elif signal_id[0:3] == 'ECE':
        # ECE scattering, expected string format e.g. ECE-te-ch01 or ECE-rho
        # That would correspond to the version 1, data for the
        # electron temperature profile
        split_out = signal_id.split('-')
        quantity = split_out[1]
        if quantity == "te":
            location_string = "ArchiveDB/raw/Minerva/Minerva.ECE.DownsampledRadiationTemperatureTimetraces/"\
                +f"relative_calibrated_signal_DATASTREAM/V1/{int(split_out[2][2:])-1}/QME-{split_out[2]}/"
        elif quantity == "rho":
            location_string = "ArchiveDB/raw/Minerva/Minerva.ECE.ColdResonance/"\
                +"cold_resonance_DATASTREAM/V1/0/rho/"
        return location_string
    else:
        return signal_id

def stream_gen_curr(signal_id):
    """ Creates the string for obtaining a number of specific signals of OP2

    :signal_id: signal id string
    :return: reduction string part of the url, if there is no signal_id
    string in stream_dict, then just return signal_id
    """
    stream_dict = dict()
    # Total ECRH power:
    stream_dict['ECRH'] = \
        'Test/codac/W7X/CBG_ECRH/TotalPower_DATASTREAM/V1/0/Ptot_ECRH/scaled/'
    # NBI source 7 beam current:
    stream_dict['NBI-7'] = \
        'ArchiveDB/codac/W7X/CoDaStationDesc.31/DataModuleDesc.24119_DATASTREAM/5/HGV_3%20Monitor%20I/scaled/'
    # NBI source 8 beam current:
    stream_dict['NBI-8'] = \
        'ArchiveDB/codac/W7X/CoDaStationDesc.31/DataModuleDesc.24119_DATASTREAM/7/HGV_4%20Monitor%20I/scaled/'
    # QSI-CXRS data
    stream_dict["QSI-CXRS"] = \
        "Test/raw/W7X/QSS_DivertorSpectroscopy/PI_CCD_06_1-QSS60OC095_DATASTREAM/0/Images/"
    # WDia data
    stream_dict["Wdia"] = \
            "ArchiveDB/raw/W7X/ControlStation.2161/Wdia_DATASTREAM/0/Diamagnetic%20Energy/"
    # WDia data
    stream_dict["Wdia"] = \
            "ArchiveDB/raw/W7X/ControlStation.2161/Wdia_DATASTREAM/0/Diamagnetic%20Energy/"

    if signal_id in stream_dict:
        return stream_dict[signal_id]
    elif signal_id[0:4] == 'Zeff':
        # signal format expected Zeff-v20-data
        split_out = signal_id.split('-')
        version = split_out[1]
        version = version[1:]
        stream_type = split_out[2]
        if stream_type == 'data':
            location_string = \
                'Test/raw/Minerva2/Minerva.Zeff.AEZ40-AET40.USB-Spectrometer/Zeff_mean_DATASTREAM/V'+version+"/0/Zeff_mean/"
        elif stream_type == 'param':
            location_string = \
                'Test/raw/Minerva2/Minerva.Zeff.AEZ40-AET40.USB-Spectrometer/Zeff_mean_PARLOG/V'+version+"/0/Zeff_mean/"
        elif stream_type == 'err':
            location_string = \
                'Test/raw/Minerva2/Minerva.Zeff.AEZ40-AET40.USB-Spectrometer/Zeff_std_DATASTREAM/V'+version+"/0/Zeff_std/"      
        return location_string
    elif signal_id[0:2] == 'TS':
        # Thomson scattering, expected string format e.g. TS-v20-ne
        # That would correspond to the version 20, data for the
        # electron density profile
        split_out = signal_id.split('-')
        version = split_out[1][1:]
        quantity = split_out[2]
        thomson_dict={
            "rho": "0/rho/",
            "ne": "1/ne/",
            "te": "2/te/",
            "ne_std": "3/ne_std/",
            "te_std": "4/te_std/",
            "R": "5/R/",
            "Phi": "6/Phi/",
            "Z": "7/Z/"}
        location_string = f"Test/raw/Minerva2/Minerva.ThomsonScattering.Profiles/Preliminary_Fast_DATASTREAM/V{version}" + \
                '/' + thomson_dict[quantity]
        return location_string
    elif signal_id[0:3] == 'ECE':
        # ECE scattering, expected string format e.g. ECE-te-ch01 or ECE-rho
        # That would correspond to the version 1, data for the
        # electron temperature profile
        split_out = signal_id.split('-')
        quantity = split_out[1]
        if quantity == "te":
            location_string = "ArchiveDB/raw/Minerva/Minerva.ECE.DownsampledRadiationTemperatureTimetraces/"\
                +f"relative_calibrated_signal_DATASTREAM/V1/{int(split_out[2][2:])-1}/QME-{split_out[2]}/"
        elif quantity == "rho":
            location_string = "ArchiveDB/raw/Minerva/Minerva.ECE.ColdResonance/"\
                +"cold_resonance_DATASTREAM/V1/0/rho/"
        return location_string
    else:
        return signal_id