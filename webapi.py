"""

.. module:: FLAP/***
   :platform: Windows/Ubuntu
   :synopsis: Signal reading from Wendelstein 7-X ArchiveDB

.. moduleauthor:: G. Cseh and M. Vecsei


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

import flap

if (flap.VERBOSE):
    print("Importing w7x_webapi")

logger = logging.getLogger('flapLogger')

# ----------------------------------------------------------------------------------------------------------------------
# -------------------------------------------Downloading via webapi-----------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# -------------------------------------------FLAP INDEPENDENT CODE------------------------------------------------------

class GetSignal(object):

    """**Bases:** :class:`object`

    This class handles the data readout from the Wendelstein 7-X ArchiveDB data
    acquisition system via the Web API interface.

    """

    def __init__(self, signal_id):
        """ Initializes the basic signal information
        """
        self.stream = self.stream_gen(signal_id)  # The address of the datastream
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
        return self.url_request(self.url)

    def time_query_gen(self, start, stop):
        """This function gives back the time_query string for the requested data

        :param start: contains the start time of the signal request
        :type start: datetime.datetime object
        :param stop: contains the end time of the signal request
        :type stop: a datetime.datetime object
        :return: the time_query part of the url

        """
        start_ns = self.convert_nanosecs(start.year, start.month, start.day,
                                         start.hour, start.minute, start.second, start.microsecond)
        self.from_time = start_ns
        stop_ns = self.convert_nanosecs(stop.year, stop.month, stop.day, stop.hour,
                                        stop.minute, stop.second, stop.microsecond)
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

    def stream_gen(self, signal_id):
        """ Creates the string for obtaining a number of specific signals

        :signal_id: signal id string
        :return: reduction string part of the url, if there is no signal_id
        string in stream_dict, then just return signal_id
        """
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


class ScalarData:
    """**

    This class handles time evolution of scalar data, e.g. ECRH power and so on.

    """

    __module__ = "customfiles.InternalData"

    def __init__(self, data, from_time=None, data_source_obj=None):
        ""
        self.time = data["dimensions"]
        self.val = data['values']
        self.label = data['label']
        self.unit = data['unit']
        if from_time is not None:
            self.from_time = from_time
        elif hasattr(data_source_obj, 'from_time'):
            self.from_time = data.from_time
        else:
            self.from_time = self.time[0]


    def time_norm(self):
        """
        Normalizes the time axis of the data. OTherwise it would be given in nanoseconds
        """
        self.time[:] = [float(i - self.from_time) * 1.0e-9 for i in self.time]


class SignalList():
    """

    The class is only used to handle multiple signals.

    """

    def __init__(self, signal_names):
            ""
            self.list = [GetSignal(signal) for signal in signal_names]

    def __getattr__(self, name):
        """
        When we try to access some function that doesn't exist in the class
        definition, treat it as a class method and call it on every signal element.
        """
        def method(*args):
            for signal in self.list:
                method_to_call = getattr(signal, name)
                method_to_call(*args)
        return method


class ScalarDataList():
    """

    This class handles multiple ScalarData objects

    """
    def __init__(self, data_list, from_time=None, data_source_obj=None):
        self.list = [ScalarData(data, from_time=from_time, data_source_obj=data_source_obj) for data in data_list]

    def __getattr__(self, name):
        """
        When we try to access some function that doesn't exist in the class
        definition, treat it as a class method and call it on every ScalarData
        element.
        """
        def method(*args):
            for data in self.list:
                method_to_call = getattr(data, name)
                method_to_call(*args)
        return method


class VectorData():
    """"
    Obtains the more complicated data structures from the archive
    """

    def __init__(self, signal_id, signal_types=None, param_types=None, volumes=None):
        """
        By using the SignalList class, the necessary structure for reading the
        data and the parameters are defined
        """
        if not (signal_types is None):
            self.signal_types = signal_types
            signal_names = [[signal_id+'-'+signal_type+'-data-vol'+str(vol_index) for vol_index in volumes]
                            for signal_type in self.signal_types]
            self.signals = [SignalList(signal) for signal in signal_names]
        if not (param_types is None):
            self.param_types = param_types
            param_names = [[signal_id+'-'+param_type+'-vol'+str(vol_index) for vol_index in volumes]
                           for param_type in self.param_types]
            self.params = [SignalList(param) for param in param_names]

    def get_data(self):
        """
        Obtaining the data and storing it in ScalarDataList
        """
        if hasattr(self, "signals"):
            archive_pull_out = [[signal_id.archive_pull() for signal_id in signal.list] for signal in self.signals]
            self.data = [ScalarDataList(data, from_time=self.signals[0].list[0].from_time)
                         for data in archive_pull_out]
        if hasattr(self, "params"):
            archive_pull_out = [[param_id.archive_pull() for param_id in param.list] for param in self.params]
            self.params = [ScalarDataList(data, from_time=self.signals[0].list[0].from_time)
                           for data in archive_pull_out]

    def squeeze(self):
        '''
        This function changes all the acquired data to a single structure to
        allow working with it
        '''
        self.data_struct = dict()
        index = 0
        if hasattr(self, "param_types"):
            for param in self.param_types:
                self.data_struct[param] = np.asarray([data.val[0] for data in self.params[index].list])
                index = index + 1
        self.data_struct["time"] = np.asarray(self.data[0].list[0].time)
        index = 0
        for quantity in self.signal_types:
            self.data_struct[quantity] = np.vstack([np.asarray(data.val) for data in self.data[index].list])
            index = index+1


class ThomsonData(VectorData):
        """"
        Currently the Thomson scattering data is uploaded in a very complciated manner:
            The output of every data of every mesurement point in space is in a
            different datastream. E.g. The time evolution of the electron temperature
            n point in space is given by one separate datastream. Moreover, the location
            of on measurement volume in space is stored in a separata parlog datastream.
            For more info look up the description in https://w7x-logbook.ipp-hgw.mpg.de/components?id=QTB#
            This class is made to simplify the processing of Thomson Scattering data
        """
        
        def __init__(self, signal_id, signal_types=None, volumes=None):
            # signal_id is expected to be in the form 'TS-v20'
            if signal_types is None:
                signal_types = ['Te_map', 'ne_map', 'pe_map', 'Te_mean', 'ne_mean', 'pe_mean',
                                'Te_low97.5', 'Te_high97.5', 'ne_low97.5', 'ne_high97.5',
                                'pe_low97.5', 'pe_high97.5', 'Te_s', 'ne_s', 'pe_s', 'signal_map_1',
                                'signal_map_2', 'signal_map_3', 'signal_map_4', 'signal_map_5',
                                'signal_mean_1', 'signal_mean_2', 'signal_mean_3', 'signal_mean_4',
                                'signal_mean_5']

            param_types = ['x_coord-param', 'y_coord-param', 'z_coord-param']

            if volumes is None:
                volumes = np.linspace(1, 16, num=16, dtype=int)

            VectorData.__init__(self, signal_id, signal_types=signal_types, param_types=param_types,
                                volumes=volumes)


class ABESData(VectorData):
        """"
        Obtains the ABESData from the archive.
        For more info look up the description in https://w7x-logbook.ipp-hgw.mpg.de/components?id=QSI#
        This class is made to simplify the processing of data
        """
        
        def __init__(self, signal_id, signal_types=None, volumes=None):
            # signal_id is expected to be in the form 'TS-v20'
            if signal_types is None:
                signal_types = ['density', 'density_error_low', 'density_error_high', 'meas_light', 'meas_light_error',
                                'x_coord', 'y_coord']

            param_types = None # ['x_coord-data', 'y_coord-data']

            if volumes is None:
                volumes = np.linspace(1, 30, num=30, dtype=int)

            VectorData.__init__(self, signal_id, signal_types=signal_types, param_types=param_types,
                                volumes=volumes)

        def get_data(self):
            """
            Obtaining the data and storing it in ScalarDataList
            """
            if hasattr(self, "signals"):
                archive_pull_out = [[signal_id.archive_pull() for signal_id in signal.list] for signal in self.signals]
                # Additional line to remove nones as not every volume is uploaded for every abes shot
                for element in range(0, len(archive_pull_out)):
                    archive_pull_out[element] = [archive_pull_out[element][i] for i, e in
                                                 enumerate(archive_pull_out[element]) if e is not None]
                    self.signals[element].list = [self.signals[element].list[i] for i, e in
                                                  enumerate(archive_pull_out[element]) if e is not None]
                self.data = [ScalarDataList(data, from_time=self.signals[0].list[0].from_time)
                             for data in archive_pull_out]

            if hasattr(self, "params"):
                archive_pull_out = [[param_id.archive_pull() for param_id in param.list] for param in self.params]              
                # Additional line to remove nones as not every volume is uploaded for every abes shot
                for element in range(0, len(archive_pull_out)):
                    archive_pull_out[element] = [archive_pull_out[element][i] for i, e in
                                                 enumerate(archive_pull_out[element]) if e is not None]
                    self.params[element].list = [self.params[element].list[i] for i, e in
                                                 enumerate(archive_pull_out[element]) if e is not None]
                self.params = [ScalarDataList(data, from_time=self.signals[0].list[0].from_time)
                               for data in archive_pull_out]

# ----------------------------------------FLAP RELATED CODE-------------------------------------------------------------


def conv_scalar_to_dataobj(data, url, data_name, exp_id, options):
    """Converts the ScalarData obtained by the code to DataObject class of flap
    """
    if not data['dimensionCount'] == 1:
        raise TypeError('conv_scalar_to_dataobj is incompatible with higher than one dimensional data, ' +
                        'current data dimension: '+str(data['dimensionCount']))

    info = "Data downloaded from: "+url+os.linesep

    # Constructing the time flap Coordinate
    name = 'Time'
    # getting the unit
    unit = 'ns'
    if 'Scale Time' in options.keys():
        if options['Scale Time'] is True:
            unit = 'Second'
    # Checking equidistance of time
    equidistant = False
    if 'Check Time Equidistant' in options.keys():
        if options['Check Time Equidistant'] is True:
            timesteps = np.array(data["dimensions"][1:])-np.array(data["dimensions"][0:-1])
            equidistance = np.linalg.norm(timesteps-timesteps[0])/np.linalg.norm(timesteps)
            if equidistance < 1e-6:
                info = info+"Time variable is taken as equidistant to an accuracy of "+str(equidistance)+os.linesep
                equidistant = True
                start = data["dimensions"][0]
                step = np.mean(timesteps)
            else:
                    info = info+"Time variable is not equidistant, deviation: "+str(equidistance)+os.linesep
    shape = [data['dimensionSize']]

    if equidistant is True:
        time_coord = flap.Coordinate(name=name, unit=unit, mode=flap.CoordinateMode(equidistant=True), shape=shape,
                                     start=start, step=[step], dimension_list=[0])
    else:
        time_coord = flap.Coordinate(name=name, unit=unit, mode=flap.CoordinateMode(equidistant=False), shape=shape,
                                     values=np.array(data["dimensions"]), dimension_list=[0])

    # Temporal sample coord
    time_sample = flap.Coordinate(name='Sample', unit=1, mode=flap.CoordinateMode(equidistant=True), shape=shape,
                                  start=1, step=[1], dimension_list=[0])

    coords = [time_coord, time_sample]

    # Constructing the DataObject
    data_unit = flap.Unit(name=data['label'], unit=data['unit'])
    data_title = data_name
    data_shape = list(np.shape(data['values']))
    d = flap.DataObject(data_array=np.array(data['values']), error=None, data_unit=data_unit, coordinates=coords,
                        exp_id=exp_id, data_title=data_title, info=info, data_shape=data_shape)
    return d

def conv_aug2_to_dataobj(data, url, data_name, exp_id, options):
    """Converts the ScalarData obtained by the code to DataObject class of flap
    """

    info = "Data downloaded from: "+url+os.linesep

    # Constructing the time flap Coordinate
    name = 'Time'
    # getting the unit
    unit = 'ns'
    if 'Scale Time' in options.keys():
        if options['Scale Time'] is True:
            unit = 'Second'
    # Checking equidistance of time
    equidistant = False
    if 'Check Time Equidistant' in options.keys():
        if options['Check Time Equidistant'] is True:
            timesteps = np.array(data["dimensions"][1:])-np.array(data["dimensions"][0:-1])
            equidistance = np.linalg.norm(timesteps-timesteps[0])/np.linalg.norm(timesteps)
            if equidistance < 1e-6:
                info = info+"Time variable is taken as equidistant to an accuracy of "+str(equidistance)+os.linesep
                equidistant = True
                start = data["dimensions"][0]
                step = np.mean(timesteps)
            else:
                    info = info+"Time variable is not equidistant, deviation: "+str(equidistance)+os.linesep
    shape = [len(data['dimensions'])]

    if equidistant is True:
        time_coord = flap.Coordinate(name=name, unit=unit, mode=flap.CoordinateMode(equidistant=True), shape=shape,
                                     start=start, step=[step], dimension_list=[0])
    else:
        time_coord = flap.Coordinate(name=name, unit=unit, mode=flap.CoordinateMode(equidistant=False), shape=shape,
                                     values=np.asarray(data["dimensions"]), dimension_list=[0])

    # Temporal sample coord
    time_sample = flap.Coordinate(name='Sample', unit='1', mode=flap.CoordinateMode(equidistant=True), shape=shape,
                                  start=1, step=[1], dimension_list=[0])
    channel_coord = flap.Coordinate(name='Channel', unit='', mode=flap.CoordinateMode(equidistant=True),
                                    shape=[np.shape(data['values'])[1]], start=1, step=[1], dimension_list=[1])
    wavelength_coord = flap.Coordinate(name='Wavelength', unit='', mode=flap.CoordinateMode(equidistant=True),
                                       shape=[np.shape(data['values'])[2]], start=0, step=np.asarray([1]), dimension_list=[2])

    coords = [time_coord, time_sample, channel_coord, wavelength_coord]

    # Constructing the DataObject
    data_unit = flap.Unit(name=data['label'], unit='None')
    data_title = data_name
    data_shape = list(np.shape(data['values']))
    d = flap.DataObject(data_array=np.array(data['values']), error=None, data_unit=data_unit, coordinates=coords,
                        exp_id=exp_id, data_title=data_title, info=info, data_shape=data_shape)
    return d

def conv_vector_to_dataobj(data, url, data_name, exp_id, options):
    """Converts the VectorData obtained by the code to DataObject class of flap
    """

    info = "Data downloaded from: "+url+os.linesep

    # Constructing the time flap Coordinate
    name = 'Time'
    # getting the unit
    unit = 'ns'
    if 'Scale Time' in options.keys():
        if options['Scale Time'] is True:
            unit = 'Second'
    # Checking equidistance of time
    equidistant = False
    if 'Check Time Equidistant' in options.keys():
        if options['Check Time Equidistant'] is True:
            timesteps = np.array(data.data[0].list[0].time[1:])-np.array(data.data[0].list[0].time[0:-1])
            equidistance = np.linalg.norm(timesteps-timesteps[0])/np.linalg.norm(timesteps)
            if equidistance < 1e-6:
                info = info+"Time variable is taken as equidistant to an accuracy of "+str(equidistance)+os.linesep
                equidistant = True
                start = data["dimensions"][0]
                step = np.mean(timesteps)
            else:
                    info = info+"Time variable is not equidistant, deviation: "+str(equidistance)+os.linesep

    shape = list(np.shape(data.data_struct[data.signal_types[0]]))

    if equidistant is True:
        time_coord = flap.Coordinate(name=name, unit=unit, mode=flap.CoordinateMode(equidistant=True), shape=shape,
                                     start=start, step=[step], dimension_list=[0])
    else:
        time_coord = flap.Coordinate(name=name, unit=unit, mode=flap.CoordinateMode(equidistant=False), shape=shape,
                                     values=np.asarray(data.data[0].list[0].time), dimension_list=[0])

    # Temporal sample coord
    time_sample = flap.Coordinate(name='Sample', unit='1', mode=flap.CoordinateMode(equidistant=True), shape=shape,
                                  start=1, step=[1], dimension_list=[0])

    coords = [time_coord, time_sample]

    index = 0
    for params in data.param_types:
        param_coord = flap.Coordinate(name=data.params[index].list[0].label, unit=data.params[index].list[0].unit,
                                      mode=flap.CoordinateMode(equidistant=False), shape=shape,
                                      values=np.asarray(data.data_struct[params]), dimension_list=[1])
        index = index+1
        coords.append(param_coord)

    # Constructing the DataObject
    data_unit = flap.Unit(name=data.data[0].list[0].label, unit=data.data[0].list[0].unit)
    data_title = data_name
    data_array = np.transpose(data.data_struct[data.signal_types[0]])
    data_shape = list(np.shape(data_array))
    d = flap.DataObject(data_array=np.asarray(data_array), error=None,
                        data_unit=data_unit, coordinates=coords, exp_id=exp_id, data_title=data_title,
                        info=info, data_shape=data_shape)
    return d


def get_data(exp_id=None, data_name=None, no_data=False, options=None, coordinates=None, data_source=None):
    """ Data read function for the W7-X Alkali BES diagnostic
    data_name: should be a string, currently either:
        ECRH
        NBI-7
        NBI-8
        TS-v20 - thomson data
        or a stream string: ArchiveDB/codac/W7X/CBG_ECRH/TotalPower_DATASTREAM/V1/0/Ptot_ECRH/scaled/
    exp_id: Experiment ID, YYYYMMDD.xxx, e.g. '20180912.012'
        if not given, then a start and a stop time should be given in the options
        coordinates: List of flap.Coordinate() or a single flap.Coordinate
                 Defines read ranges. The following coordinates are interpreted:
                     'Time [s]': The read times
                     Only a single equidistant range is interpreted.
    options:
        Time Start and Time Stop times can be given (see exp_id) e.g.
         options={'Time Start':datetime.datetime(2018, 8, 22, 10, 46, 53),
                  'Time Stop': datetime.datetime(2018, 8, 22, 10, 47, 8)}
        Downsample - the sampling of the data
        Scale Time - if set and True, then the data is normalized to seconds and 0s point is set to
                     the start of the shot
        Check Time Equidistant - if set and True, then the data is checked if it is equidistant in time.
                              it is assumed to be equidistant, if the difference between the timesteps relatively small
    """

    if no_data:
        raise NotImplementedError('no data input is not implemented yet.')

    get_vector = False
    if data_name[0:3] == "TS-":
        data_setup = GetSignal(data_name+"-data-vol1")
        get_vector = True
    elif data_name[0:5] == "ABES-":
        data_setup = GetSignal(data_name+"-data-vol15")
        get_vector = True
    else:
        data_setup = GetSignal(data_name)
    # Generating the data relevant time intervals
    if exp_id is None:
        if not ('Time Start' in options.keys() and 'Time Stop' in options.keys()):
            raise ValueError("either exp_id should be given, or options should have a start and stop key")
        else:
            data_setup.time_query_gen(options['Time Start'], options['Time Stop'])
    else:
        data_setup.shotid_time_query_gen(exp_id)

    start_shift = decimal.Decimal(0)
    if coordinates is not None:
        coord = coordinates[0]
        if (coord.unit.name == 'Time [s]') and (coord.mode.equidistant):
                read_range = [float(coord.c_range[0]), float(coord.c_range[1])]
        else:
            raise ValueError("Only timeranges may be defined for reading")
        orig_times = data_setup.time_query.split('=')
        orig_start = decimal.Decimal((orig_times[1].split('&'))[0])
        start_shift = decimal.Decimal((read_range[0]+1.0)*1e9)
        new_start = orig_start + start_shift
        new_stop = orig_start + decimal.Decimal((read_range[1]+1.0)*1e9)
        data_setup.time_query = "_signal.json?from=" + str(new_start) + '&upto=' + str(new_stop)

    # Downsampling the data if requested
    if 'Downsample' in options.keys():
        warnings.warn('Downsampling requests do not seem to work properly on the W7-X webapi')
        data_setup.downsamp_gen(options['Downsample'])

    # Converting data to flap.DataObject format
    if get_vector is True:
        if data_name[0:3] == "TS-":
            # data is expected in the form 'TS-v20-ne_map'
            split_out = data_name.split('-')
            version_data = split_out[0]+"-"+split_out[1]
            signal_type = split_out[2]
            vector_data = ThomsonData(version_data, signal_types=[signal_type])
        elif data_name[0:5] == "ABES-":
            # data is expected in the form 'TS-v20-ne_map'
            split_out = data_name.split('-')
            version_data = split_out[0]+"-"+split_out[1]
            signal_type = split_out[2]
            vector_data = ABESData(version_data, signal_types=[signal_type])
        else:
            raise NotImplementedError('Only Thomson and Alkali beam data can be currently read')

        times = data_setup.time_query.split('=')
        start = decimal.Decimal((times[1].split('&'))[0])
        end = decimal.Decimal((times[2].split('&'))[0])
        [x.time_query_gen_ns(start, end) for x in vector_data.signals]
        [x.time_query_gen_ns(start, end) for x in vector_data.params]
        vector_data.get_data()
        vector_data.squeeze()

        # Normalizing the time
        if 'Scale Time' in options.keys():
            if options['Scale Time'] is True:
                [x.time_norm() for x in vector_data.data]
                if coordinates is not None:
                    coord = coordinates[0]
                    if (coord.unit.name == 'Time [s]') and (coord.mode.equidistant):
                        vector_data.data[0].list[0].time = [element + read_range[0] for element in
                                                            vector_data.data[0].list[0].time]
                        for x in vector_data.params:                    
                            x.list[0].time = [element + read_range[0] for element in x.list[0].time]

        d = conv_vector_to_dataobj(vector_data, data_setup.url, data_name, exp_id, options)
    else:
        # Downloading data
        data = data_setup.archive_pull()

        # Normalizing the time
        if 'Scale Time' in options.keys():
            if options['Scale Time'] is True:
                mintime = min(data["dimensions"])+1000000000-start_shift
                data["dimensions"][:] = [float(i - mintime) * 1.0e-9 for i in data["dimensions"]]
        if data_name[0:4] == "AUG-":
            d = conv_aug2_to_dataobj(data, data_setup.url, data_name, exp_id, options)
        else:
            d = conv_scalar_to_dataobj(data, data_setup.url, data_name, exp_id, options)
    return d


def add_coordinate(data_object, new_coordinates, options=None):
    raise NotImplementedError("Coordinate conversions not implemented yet.")


# ----------------------------------------------------------------------------------------------------------------------
# ---------------------------------------Uploading via webapi-----------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

class WriteABESSignal(GetSignal):
    def __init__(self, datalist, shotid, minindex=0, version=1.0):
        self.dens = datalist[0]
        self.light_meas = datalist[1]
#        self.light_recon = datalist[2]
        self.shotid = shotid
        self.minindex = minindex
        self.version = version
        self.json={}
#        self.set_params(data,shotid)
#        self.set_data(data)

    def shotid_time_gen(self):
        """ This one creates the whole the time query for the data using
        the shotid. It sends a request to the archive to get the shotid
        timerange.

        :shotid: The ID of the shot given as a string.
        """
        url = "http://archive-webapi.ipp-hgw.mpg.de/programs.json?from=" + self.shotid
        response = self.url_request(url)
        from_time = decimal.Decimal(response['programs'][0]['trigger']['1'][0])
        upto_time = decimal.Decimal(response['programs'][0]['trigger']['4'][0])
        return [from_time, upto_time]

    def set_data(self):
        #setting the time
        [from_time, upto_time] = self.shotid_time_gen()
        time_vec=np.unique(self.dens.coordinate('Time')[0])
        self.dimensions = [from_time + decimal.Decimal(int(time_vec[index]*1e9)) for index in range(0,len(time_vec))]

        temp = {}
        spat_coord = np.unique(np.unique(self.dens.coordinate('Beam axis')[0]))
        index=0
        for location in spat_coord:
            index = index + 1
            dens_loc = self.dens.slice_data(slicing={'Beam axis':location})
            density = dens_loc.data*1e19
            density_low_error = (density-dens_loc.error*1e19).clip(0)
            density_high_error = density+dens_loc.error*1e19
            light_meas_loc = self.light_meas.slice_data(slicing={'Beam axis':location})
            self.data = {
                    "datatype" : "float",
                    "values": [density.tolist(),
                               density_low_error.tolist(),
                               density_high_error.tolist(),
                               light_meas_loc.data.tolist(),
                               light_meas_loc.error.tolist(),
                               dens_loc.coordinate('Device x')[0].tolist(),
                               dens_loc.coordinate('Device y')[0].tolist(),
                               dens_loc.coordinate('Device z')[0].tolist()],
                    "dimensions": [int(element) for element in self.dimensions]}
            volume_name = "vol_"+str(index+1+self.minindex)
            temp[volume_name] = json.dumps(self.data).encode('utf-8')
        self.json["data"] = temp
        
    def set_params(self):
        self.params = {
                "label": "parms",
                "values": [{
                        "chanDescs": {
                                "[0]": {"name": "density",
                                        "active": 1,
                                        "description": "Reconstructed density profile",
                                        "physicalQuantity": {"from": 0.0, "upto": 0.0, "type":"m^-3"}},
                                "[1]": {"name": "density_error_low",
                                        "active": 1,
                                        "description": "Lower limit of reconstructed density profile error",
                                        "physicalQuantity": {"from": 0.0, "upto": 0.0, "type":"m^-3"}},
                                "[2]": {"name": "density_error_high",
                                        "active": 1,
                                        "description": "Upper limit of reconstructed density profile error",
                                        "physicalQuantity": {"from": 0.0, "upto": 0.0, "type":"m^-3"}},
                                "[3]": {"name": "meas_light",
                                        "active": 1,
                                        "description": "measured light profile",
                                        "physicalQuantity": {"from": 0.0, "upto": 0.0, "type":"V"}},
                                "[4]": {"name": "meas_light_error",
                                        "active": 1,
                                        "description": "measured light profile error",
                                        "physicalQuantity": {"from": 0.0, "upto": 0.0, "type":"V"}},
                                "[5]": {"name": "x_coord",
                                        "active": 1,
                                        "description": "x position",
                                        "physicalQuantity": {"from": 0.0, "upto": 0.0, "type":"m"}},
                                "[6]": {"name": "y_coord",
                                        "active": 1,
                                        "description": "y position",
                                        "physicalQuantity": {"from": 0.0, "upto": 0.0, "type":"m"}},
                                "[7]": {"name": "z_coord",
                                        "active": 1,
                                        "description": "z position",
                                        "physicalQuantity": {"from": 0.0, "upto": 0.0, "type":"m"}}
                                }
                            }],
                "dimensions": [int(element) for element in self.dimensions]
                }
        self.json["params"] = json.dumps(self.params).encode('utf-16')
    
    def create_version(self, sandbox=False):
            version = { 
                    "versionInfo" : [ 
                            {
                                    "reason": "updated reconstruction and error approximation method",
                                    "producer": "mive",
                                    "code_release": "v"+str(self.version),
                                    "analysis_environment": "flap"
                                    }
                            ]
                    }
            json_version = json.dumps(version).encode('utf-16')
            if sandbox is True:
                url='http://archive-webapi.ipp-hgw.mpg.de/Sandbox/raw/W7X/QSI/'
            else:
                url='http://archive-webapi.ipp-hgw.mpg.de/Test/raw/W7X/QSI/'
            
            for key in self.json["data"].keys():
                try:
                    url_data = url + key  + "_DATASTREAM"+"/_versions.json"
                    url_parms = url + key  + "_PARLOG"+"/_versions.json"
                    
                    
                    req = urllib.request.Request(url_data)
                    req.add_header('Content-Type', 'application/json; charset=utf-16')
    #                req.get_method = lambda: "POST"
                    urllib.request.urlopen(req, json_version)
                    req = urllib.request.Request(url_parms)
                    req.add_header('Content-Type', 'application/json; charset=utf-16')
                    urllib.request.urlopen(req, json_version)
                except urllib.error.HTTPError as e:
                    print('error for creating version for '+key)
                    print(e)
                    print(e.headers)

    def upload_json(self, sandbox=None):
            if hasattr(self, 'sandbox') and sandbox == None:
                sandbox = self.sandbox
            if sandbox is True:
                url = 'http://archive-webapi.ipp-hgw.mpg.de/Sandbox/raw/W7X/QSI/'
            else:
                url = 'http://archive-webapi.ipp-hgw.mpg.de/Test/raw/W7X/QSI/'
    
            for key in self.json["data"].keys():
                print('uploading '+key)
                try:
                    url_parms = url + key+"_PARLOG/V"+str(self.version)
                    req = urllib.request.Request(url_parms)
                    req.add_header('Content-Type', 'application/json; charset=utf-8')
                    urllib.request.urlopen(req, self.json["params"])
                    url_data = url + key+"_DATASTREAM/V"+str(self.version)
                    req = urllib.request.Request(url_data)
                    req.add_header('Content-Type', 'application/json; charset=utf-8')
                    urllib.request.urlopen(req, self.json["data"][key])
                except urllib.error.HTTPError as e:
                    print(e)
                    print(e.headers)


def register(data_source=None):
    flap.register_data_source('W7X_WEBAPI', get_data_func=get_data, add_coord_func=add_coordinate)

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------


class Points3D:
    '''
    Used for storing a number of points and to perform coordinate transformations
    If a point is added in either xyz of radius - theta - z coordinates, an automatic transformation to
    the other coordinate system is performed. Not yet compatible with FLAP coordinates
    '''

    __module__ = "customfiles.InternalData"

    def __init__(self):
        self.xyz = np.empty([0, 3])  # point coordinates in stellarator xyz
        self.rtz = np.empty([0, 3])  # point coordinates in stellarator radius-theta-z coordinates
        self.reff = np.empty([0])  # a vector storing the effective radius of every point
        self.shotID = None
        self.ref_eq = None  # the w7x vmec equilibrium reference ID of the shot

    def append_xyz(self, xyz_list):
        '''
        Add points to the object based on their xyz coordinates. xyz_list has to be
        either a numpy array of length 3 defining the point location, or a two-dimensional
        numpy array.
        '''
        if xyz_list.ndim == 1:
            xyz_list = np.asarray([xyz_list])
        self.xyz = np.append(self.xyz, xyz_list, axis=0)
        self.rtz = self.xyz_to_rtz(self.xyz)

    def append_rtz(self, rtz_list):
        '''
        Add points to the object based on their rtz coordinates. rtz_list has to be
        either a numpy array of length 3 defining the point location, or a two-dimensional
        numpy array.
        '''
        if rtz_list.ndim == 1:
            rtz_list = np.asarray([rtz_list])
        self.rtz = np.append(self.rtz, rtz_list, axis=0)
        self.xyz = self.rtz_to_xyz(self.rtz)

    @staticmethod
    def rtz_to_xyz(rtz_list):
        '''
        Convert rtz coordinates to xyz and returns the value. rtz_list has to be
        either a numpy array of length 3 defining the point location, or a two-dimensional
        numpy array. Practically list of points.
        '''
        if rtz_list.ndim == 1:
            rtz_list = np.asarray([rtz_list])
        x_list = rtz_list[:, 0]*np.cos(rtz_list[:, 1])
        y_list = rtz_list[:, 0]*np.sin(rtz_list[:, 1])
        return np.squeeze(np.transpose([x_list, y_list, rtz_list[:, 2]]))

    @staticmethod
    def xyz_to_rtz(xyz_list):
        '''
        Convert xyz coordinates to rtz and returns the value. xyz_list has to be
        either a numpy array of length 3 defining the point location, or a two-dimensional
        numpy array. Practically list of points.
        '''
        if xyz_list.ndim == 1:
            xyz_list = np.asarray([xyz_list])
        r_list = np.linalg.norm(xyz_list[:, 0:2], axis=1)
        t_list = np.arctan2(xyz_list[:, 1], xyz_list[:, 0])
        return [np.squeeze(np.transpose([r_list, t_list, xyz_list[:, 2]]))]

    def xyz_to_reff(self, shotID=None):
        '''
        Obtain the effective radius based on self.xyz of the object and vmec results
        '''
        if shotID:
            self.shotID = shotID
            self.ref_eq = get_refeq(shotID)
            if self.ref_eq:
                url = 'http://svvmec1.ipp-hgw.mpg.de:8080/vmecrest/v1/'+self.ref_eq+'/reff.json?'
            else:
                info = 'Points3D ' + shotID + 'error: No data at ' + url + "\n" + \
                        "The information at http://svvmec1.ipp-hgw.mpg.de:8080/vmecrest/v1/w7x_ref may help"
                logger.error(info)
                print(info)
        self.reff = np.empty([0])
        for point in self.xyz:
            url_str = url + 'x='+str(point[0])+'&y='+str(point[1])+'&z='+str(point[2])
            data = GetSignal.url_request(url_str)
            self.reff = np.append(self.reff, data['reff'][0])


def get_lcfs(shotID, phi):
    """
    Used for getting the location of the last closed flux surface at a given phi location.
    Phi should be given in degrees
    Return the r and z coordinates
    """
    ref_eq = get_refeq(shotID)
    if ref_eq:
        url = 'http://svvmec1.ipp-hgw.mpg.de:8080/vmecrest/v1/'+ref_eq+'/lcfs.json?phi='+str(phi)
    else:
        info = 'No refrence ID for shot' + str(shotID) + "\n" + \
                "The information at http://svvmec1.ipp-hgw.mpg.de:8080/vmecrest/v1/w7x_ref may help"
        logger.error(info)
        print(info)
        return -1
    data = GetSignal.url_request(url)
    r = data['lcfs'][0]['x1']
    z = data['lcfs'][0]['x3']
    return r, z


def get_refeq(shotID):
    """
    Used for getting the reference number for a w7x configurations from a shotID
    Input of shotID is a string, output is a sting in the format e.g. w7x_ref_120
    """
    source_stream = 'ArchiveDB/raw/W7XAnalysis/Equilibrium/RefEq_PARLOG/V1/parms/equilibriumID/'
    refeq = GetSignal(source_stream)
    refeq.shotid_time_query_gen(shotID)
    data = refeq.archive_pull()
    ref_eq = data['values'][0]
    return ref_eq
