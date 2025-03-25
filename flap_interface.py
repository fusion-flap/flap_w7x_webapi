"""

.. module:: FLAP/***
   :platform: Windows/Ubuntu
   :synopsis: Signal reading from Wendelstein 7-X ArchiveDB

.. moduleauthor:: G. Cseh and M. Vecsei


"""
import os
import io
import pickle
import sys  # For error messages
import warnings
import datetime  # For handling the timing
import calendar  # For handling the timing
import urllib.request  # For getting data from Web API (JSON)
import json  # For handling downloaded JSON files
import logging
import numpy as np
import decimal
import copy
import configparser

import flap

from . import archive_signal
from . import spec_data

if (flap.VERBOSE):
    print("Importing w7x_webapi")

logger = logging.getLogger('flapLogger')

def standard_dataobject(dataobject):
    '''Converts the input dataobject to a standard flap dataobject so that it can be saved and used by other users of
    flap without having to use flap_w7x_webapi
    INPUT: data - the dataobject to be converted
    OUTPUT: a standardized flap dataobject
    '''
    new_data = flap.DataObject(data_array=dataobject.data, data_unit=dataobject.data_unit,
                                        coordinates=dataobject.coordinates, exp_id=dataobject.exp_id,
                                        data_title=dataobject.data_title, data_shape=dataobject.shape,
                                        error=dataobject.error, info=dataobject.info)
    new_data.shape = dataobject.shape
    return copy.deepcopy(new_data)

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
    shape = [len(data["dimensions"])]

    if equidistant is True:
        time_coord = flap.Coordinate(name=name, unit=unit, mode=flap.CoordinateMode(equidistant=True), shape=shape,
                                     start=start, step=[step], dimension_list=[0])
    else:
        time_coord = flap.Coordinate(name=name, unit=unit, mode=flap.CoordinateMode(equidistant=False), shape=shape,
                                     values=np.array(data["dimensions"]), dimension_list=[0])

    # Temporal sample coord
    time_sample = flap.Coordinate(name='Sample', unit='Sample number', mode=flap.CoordinateMode(equidistant=True), shape=shape,
                                  start=1, step=[1], dimension_list=[0])

    coords = [time_coord, time_sample]

    # Constructing the DataObject
    data_unit = flap.Unit(name=data['label'], unit=data['unit'])
    data_title = data_name
    data_shape = list(np.shape(data['values']))
    d = flap.DataObject(data_array=np.array(data['values']), error=None, data_unit=data_unit, coordinates=coords,
                        exp_id=exp_id, data_title=data_title, info=info, data_shape=data_shape)
    return d

def conv_matrix_to_dataobj(data, url, data_name, exp_id, options):
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
    shape = [len(data["dimensions"])]

    if equidistant is True:
        time_coord = flap.Coordinate(name=name, unit=unit, mode=flap.CoordinateMode(equidistant=True), shape=shape,
                                     start=start, step=[step], dimension_list=[0])
    else:
        time_coord = flap.Coordinate(name=name, unit=unit, mode=flap.CoordinateMode(equidistant=False), shape=shape,
                                     values=np.array(data["dimensions"]), dimension_list=[0])

    # Temporal sample coord
    time_sample = flap.Coordinate(name='Sample', unit='Sample number', mode=flap.CoordinateMode(equidistant=True), shape=shape,
                                  start=1, step=[1], dimension_list=[0])

    coords = [time_coord, time_sample]

    # Constructing the DataObject
    data_unit = flap.Unit(name=data['label'], unit=data['unit'])
    data_title = data_name
    data_shape = list(np.shape(data['values']))
    for buffer_coordinate in range(len(data_shape)-1):
        buffer_coord = flap.Coordinate(name=f"Coord {buffer_coordinate+1}", unit="n.a.",
                                       mode=flap.CoordinateMode(equidistant=True), 
                                       shape=data_shape[buffer_coordinate+1], start=0, step=[1],
                                       dimension_list=[buffer_coordinate+1])
        coords += [buffer_coord]
    #raise ValueError("stop")
    d = flap.DataObject(data_array=np.array(data['values']), error=None, data_unit=data_unit, coordinates=coords,
                        exp_id=exp_id, data_title=data_title, info=info, data_shape=data_shape)
    return d

def conv_thomson_to_dataobj(data, url, data_name, exp_id, options):
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
            timesteps = np.array(data[0]["dimensions"][1:])-np.array(data[0]["dimensions"][0:-1])
            equidistance = np.linalg.norm(timesteps-timesteps[0])/np.linalg.norm(timesteps)
            if equidistance < 1e-6:
                info = info+"Time variable is taken as equidistant to an accuracy of "+str(equidistance)+os.linesep
                equidistant = True
                start = data[0]["dimensions"][0]
                step = np.mean(timesteps)
            else:
                    info = info+"Time variable is not equidistant, deviation: "+str(equidistance)+os.linesep
    shape = [len(data[0]['dimensions'])]

    if equidistant is True:
        time_coord = flap.Coordinate(name=name, unit=unit, mode=flap.CoordinateMode(equidistant=True), shape=shape,
                                     start=start, step=[step], dimension_list=[0])
    else:
        time_coord = flap.Coordinate(name=name, unit=unit, mode=flap.CoordinateMode(equidistant=False), shape=shape,
                                     values=np.asarray(data[0]["dimensions"]), dimension_list=[0])

    # Temporal sample coord
    time_sample = flap.Coordinate(name='Sample', unit='1', mode=flap.CoordinateMode(equidistant=True), shape=shape,
                                  start=1, step=[1], dimension_list=[0])
    reff_coord = flap.Coordinate(name='Normalized effective radius', unit='1', mode=flap.CoordinateMode(equidistant=False), shape=shape,
                                  values=data[2]["values"])
    r_coord = flap.Coordinate(name='Device R', unit='m', mode=flap.CoordinateMode(equidistant=False), shape=shape,
                                  values=data[3]["values"])
    phi_coord = flap.Coordinate(name='Device phi', unit='rad', mode=flap.CoordinateMode(equidistant=False), shape=shape,
                                  values=data[4]["values"])
    z_coord = flap.Coordinate(name='Device z', unit='m', mode=flap.CoordinateMode(equidistant=False), shape=shape,
                                  values=data[5]["values"])

    coords = [time_coord, time_sample, reff_coord, r_coord, phi_coord, z_coord]

    # Constructing the DataObject
    if data_name.split("-")[-1] =="ne":
        data_unit = flap.Unit(name="Electron density", unit='1e19m-3')
    else:
        data_unit = flap.Unit(name="Electron temperature", unit='keV')
    data_title = data_name
    data_shape = list(np.shape(data[0]['values']))
    d = flap.DataObject(data_array=np.array(data[0]['values']), error=np.array(data[1]['values']), data_unit=data_unit, coordinates=coords,
                        exp_id=exp_id, data_title=data_title, info=info, data_shape=data_shape)
    return d

def conv_ece_to_dataobj(data, url, data_name, exp_id, options):
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
            timesteps = np.array(data[1]["dimensions"][1:])-np.array(data[1]["dimensions"][0:-1])
            equidistance = np.linalg.norm(timesteps-timesteps[0])/np.linalg.norm(timesteps)
            if equidistance < 1e-6:
                info = info+"Time variable is taken as equidistant to an accuracy of "+str(equidistance)+os.linesep
                equidistant = True
                start = data[1]["dimensions"][0]
                step = np.mean(timesteps)
            else:
                    info = info+"Time variable is not equidistant, deviation: "+str(equidistance)+os.linesep
    shape = [len(data[1]['dimensions'])]

    if equidistant is True:
        time_coord = flap.Coordinate(name=name, unit=unit, mode=flap.CoordinateMode(equidistant=True), shape=shape,
                                     start=start, step=[step], dimension_list=[0])
    else:
        time_coord = flap.Coordinate(name=name, unit=unit, mode=flap.CoordinateMode(equidistant=False), shape=shape,
                                     values=np.asarray(data[1]["dimensions"]), dimension_list=[0])

    # Temporal sample coord
    time_sample = flap.Coordinate(name='Sample', unit='1', mode=flap.CoordinateMode(equidistant=True), shape=shape,
                                  start=1, step=[1], dimension_list=[0])
    order = np.argsort(data[0]["values"][0]) # order the data according to the location
    reff_coord = flap.Coordinate(name='Normalized effective radius', unit='1', mode=flap.CoordinateMode(equidistant=False), shape=[len(data[0]["values"][0])],
                                  values=np.asarray(data[0]["values"][0])[order], dimension_list=[1])
    
    temp_data = np.zeros([len(data[1]['dimensions']), len(data[0]["values"][0])])
    for channel in np.arange(len(data[0]["values"][0])):
        index = np.where(order==channel)[0][0]
        temp_data[:,index] = data[channel+1]["values"]

    coords = [time_coord, time_sample, reff_coord]

    # Constructing the DataObject
    data_unit = flap.Unit(name="Electron temperature", unit='keV')
    data_title = data_name
    data_shape = temp_data.shape
    d = flap.DataObject(data_array=temp_data, data_unit=data_unit, coordinates=coords,
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

    shape = list(np.shape(np.transpose(data.data_struct[data.signal_types[0]])))

    if equidistant is True:
        time_coord = flap.Coordinate(name=name, unit=unit, mode=flap.CoordinateMode(equidistant=True), shape=[shape[0]],
                                     start=start, step=[step], dimension_list=[0])
    else:
        time_coord = flap.Coordinate(name=name, unit=unit, mode=flap.CoordinateMode(equidistant=False), shape=[shape[0]],
                                     values=np.asarray(data.data[0].list[0].time), dimension_list=[0])

    # Temporal sample coord
    time_sample = flap.Coordinate(name='Sample', unit='1', mode=flap.CoordinateMode(equidistant=True), shape=[shape[0]],
                                  start=1, step=[1], dimension_list=[0])

    coords = [time_coord, time_sample]

    index = 0
    for params in data.param_types:
        param_coord = flap.Coordinate(name=data.params[index].list[0].label, unit=data.params[index].list[0].unit,
                                      mode=flap.CoordinateMode(equidistant=False), shape=[shape[1]],
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

def get_data_v1(exp_id=None, data_name=None, no_data=False, options={}, coordinates=None, data_source=None):
    """ Data read function for the W7-X Alkali BES diagnostic
    data_name: should be a string, currently either:
        A full webapi name 
            e.g. ArchiveDB/codac/W7X/CBG_ECRH/TotalPower_DATASTREAM/V1/0/Ptot_ECRH/scaled/
        A webapi name with some unix-style regular expression. 
            e.g. ArchiveDB/raw/Minerva/Minerva.ECE.RadiationTemperatureTimetraces/relative_calibrated_signal_DATASTREAM/V1/0/QME-ch[01-12]//
            This example will read 12 channels from the ECE diagnostic into one data object.
        A list of strings
            The listed signals will be read into one data object, if the time scales are the same.
        A name or list of names which are listed in <Virtual name file>
            If the virtual signal is a complex signal it will generate a complex data object.
            If the virtual signal is a list of signls (...,...,...) these will be put as multiple signals in the data object.       
        Regular expression processing will be applied both at the data_name level and after translation using the
        <Virtual name file>
    exp_id: Experiment ID, YYYYMMDD.xxx, e.g. '20180912.012'
        if not given, then a start and a stop time should be given in the options
    coordinates: List of flap.Coordinate() or a single flap.Coordinate
        Defines read ranges. The following coordinates are interpreted:
            'Time [s]': The read times
                        Only a single equidistant range is interpreted.
    options:
        Time Start and Time Stop : datetime.datetime
            Times can be given as an alternative to exp_id e.g.
            options={'Time Start':datetime.datetime(2018, 8, 22, 10, 46, 53),
                     'Time Stop': datetime.datetime(2018, 8, 22, 10, 47, 8)}
        Virtual name file : str 
            A file name to translate virtual names to webapi entries. 
            For format see webapi_virtual_names()

        Downsample : float, default is None
            The sampling of the data to this frequency.
        Scale Time : bool, default is False
            if set and True, then the data is normalized to seconds and 0s point is set to
            the start of the shot
        Check Time Equidistant : bool , default is False
            If set and True, then the data is checked if it is equidistant in time.
            it is assumed to be equidistant, if the difference between the timesteps relatively small
        Cache Data : bool, default is None
            Store a local copy of the data so as next time it is read from there.
        Cache Directory : str
            Cache directory name. If None will store data under the source file directory.
        Renew Cache : bool
            Read from webapi even if data is present in cache.
        V1 : bool, default is False
            If True the new version of the program will be used: get_data_v1
        Verbose : bool
            Print diagnostic messages.
    """

    options_default = {'Cache Data': True,
                       'Cache Directory': None,
                       'Renew Cache': False,
                       'Scale Time':False,
                       'Time Start': None,
                       'Time Stop': None,
                       'Downsample': None,
                       'Check Time Equidistant': False,
                       'V1' : False,
                       'Virtual name file' : None,
                       'Verbose' : False
                       }
  
    _options = flap.config.merge_options(options_default,options,data_source='W7X_WEBAPI')
 
    # This is a date until which a 1 second shift in the time axis is applied
    start_shift_date = '20251231'
    print("Start shift date needs to be discussed!!!")
    
    if ((exp_id is not None) and (type(exp_id) is not str) or (len(exp_id) != 12)):
        raise TypeError("exp_id should be string in format YYYYMMDD.xxx")
        
    if ((exp_id is None) and (_options['Time Start'] is None) and (_options['Time Stop'] is None)):
        raise ValueError("Either exp_id or start and stop times should be set.")
    
    data_description = flap.interpret_signals(data_name,exp_id=exp_id,virtual_signal_file=_options['Virtual name file'])
    
    if coordinates is not None:
        try:
            coord = coordinates[0]
            if (coord.unit.name == 'Time') and (coord.mode.equidistant):
                read_range = [float(coord.c_range[0]), float(coord.c_range[1])]
            else:
                raise ValueError("Only timeranges may be defined for reading")
        except IndexError:
            read_range = None
    else:
        read_range = None
            
    # Assembling a list of signals and list of webapi request
    # This will collect the data of all virtual signals
    signal_list = []
    # This will collect the names of all signals.
    signal_name_list = []
    # This will collect the common timescale for all signals
    common_time = None    
    # This will contain the data type (complex or float) of the data
    common_dtype = None
    # This will contain the data unit returned by webapi
    common_unit = None
    

    if (data_description.data_type == 'real'):
        complex_data = False
    elif (data_description.data_type == 'complex'):
        complex_data = True
    else:
        raise TypeError("Invalid signal type. This is an internal error.") 
                 
    # This cycle goes through the singals listed in data_description   
    for i_signal,signal in enumerate(data_description.signal_list):
        # Assembling a list of webapi nodes needed for this data
        if (complex_data):
            webapi_request_list = signal
            signal_name_list.append('complex({:s},{:s})'.format(signal[0],signal[1]))
        else:
            webapi_request_list = [signal]
            signal_name_list.append(signal)

        # This variable collects data for the constituents of one data_name
        webapi_data_list = []
        # Common parameters for this data _name
        data_ndim = None
        data_shape = None
        data_unit = None
        
        # Checing from small possible bugs in signal names
        for i in range(len(webapi_request_list)):
            while (webapi_request_list[i][0] == '/'):
                webapi_request_list[i] = webapi_request_list[i][1:]
            if (webapi_request_list[i][-1] != '/'):
                webapi_request_list[i] += '/'
                
        for webapi_name in webapi_request_list: 

            # Trying to read cached data
            if (_options['Cache Data'] and (_options['Cache Directory'] != None)):
                filename = str(exp_id)+'_'+webapi_name
                for c in ['\\',':','/']:
                   filename = filename.replace(c,'_')
                directory=os.path.join(_options['Cache Directory'],str(exp_id))
                filename = os.path.join(directory,filename+'.pickle')
                if (not _options['Renew Cache']):
                    while True:
                        if (_options['Verbose']):
                            print("Reading cache file: '{:s}'.".format(filename),flush=True)
                        try:
                            f = io.open(filename,'rb')
                        except:
                            data_cached = False
                            if (_options['Verbose']):
                                print("Cannot open cache file.",flush=True)
                            break
                        try:
                            webapi_pickle = pickle.load(f)
                            this_data = webapi_pickle[0]
                            this_time = webapi_pickle[1]
                            exp_id_from_cache = webapi_pickle[2]
                            shot_ref_time_from_cache = webapi_pickle[3]
                            this_unit = webapi_pickle[4]
                        except Exception:
                            print("Invalid cache file: {:s}".format(filename),flush=True)
                            data_cached = False
                            f.close()
                            break
                        f.close()
                        data_cached = True
                        if (exp_id_from_cache is not None):
                            ref_time = shot_ref_time_from_cache
                            shot_ref_time = shot_ref_time_from_cache
                        else:
                            ref_time = 0
                            shot_ref_time = None
                        # Checking whether the cached data covers the requested time
                        if (read_range is not None):
                            if (exp_id is None):
                                # In this case the read_range is an absolute time
                                this_tres = this_time[1] - this_time[0]
                                if ((read_range[0] < this_time[0] - this_tres) or (read_range[1] > this_time[-1] + this_tres)):
                                    if (_options['Verbose']):
                                        print("Requested time range is outside of time range in cache file.",flush=True)
                                    data_cached = False
                                else:
                                    ind = np.nonzero(np.logical_and(this_time >= read_range[0],
                                                                    this_time <= read_range[1]
                                                                    )[0]
                                                    )
                                    if (len(ind) == 0):
                                        if (_options['Verbose']):
                                            print("No data in cache file in requested time range.",flush=True)
                                        data_cached = False
                                    else:
                                        this_time = this_time[ind]
                                        this_data = this_data[ind]
                            else:
                                if (exp_id_from_cache is None):
                                    # The cached data is not linked to an exp_id. We don't know here the 
                                    # time of the experiment, we don't use the cached data
                                    if (_options['Verbose']):
                                        print("Data in cache file is not linked to experiment.",flush=True)
                                    data_cached = False
                                else:
                                    if (read_range[0] < 0):
                                        read_range[0] = 0
                                    # In this case the read_range is a relative time in the experiment
                                    this_tres = this_time[1] - this_time[0]
                                    if ((read_range[0] < (this_time[0] - shot_ref_time_from_cache - this_tres) /1e9) 
                                        or (read_range[1] > (this_time[-1] - shot_ref_time_from_cache + this_tres) / 1e9)):
                                            data_cached = False
                                            if (_options['Verbose']):
                                                print("Requested time range ({:7.3f}-{:7.3f}) is outside of time range in cache file ({:7.3f}-{:7.3f}).".format(read_range[0],
                                                                                                                                                read_range[1],
                                                                                                                                                (this_time[0] - shot_ref_time_from_cache - this_tres) /1e9,
                                                                                                                                                (this_time[-1] - shot_ref_time_from_cache + this_tres) / 1e9),
                                                                                                                                            flush=True
                                                      )
                                    else:
                                        ind = np.nonzero(np.logical_and((this_time - shot_ref_time_from_cache) /1e9 >= read_range[0],
                                                                        (this_time - shot_ref_time_from_cache) /1e9 <= read_range[1])
                                                         )[0]
                                        if (len(ind) == 0):
                                            data_cached = False
                                        else:
                                            this_time = this_time[ind]
                                            this_data = this_data[ind]

                        break    
                else:
                    data_cached = False
            else:
                data_cached = False
                
            if (not data_cached):
                data_setup = archive_signal.ArchiveSignal(webapi_name, exp_id=exp_id)
                # Generating the data relevant time intervals
                if exp_id is None:
                    if not ((_options['Time Start'] is not None) and (_options['Time Stop'] is not None)):
                        raise ValueError("either exp_id should be given, or options should have a 'Time Start' and 'Time Stop' key.")
                    else:
                        data_setup.time_query_gen(_options['Time Start'], _options['Time Stop'])
                        ref_time = 0 
                        shot_ref_time = None
                else:
                    try:
                        data_setup.shotid_time_query_gen(exp_id)
                        orig_times = data_setup.time_query.split('=')
                        orig_start = int((orig_times[1].split('&'))[0])
                        orig_stop = int(orig_times[2])
                        shot_ref_time = int(data_setup.time_query.split('=')[1].split('&')[0])
                        if int(exp_id[:8]) <= int(start_shift_date):
                            start_shift = 1000000000
                        else:
                            start_shift = 0
                        new_start = orig_start + start_shift
                        new_stop = orig_stop + start_shift
                        shot_ref_time = new_start
                        ref_time= shot_ref_time
                        data_setup.time_query = "_signal.json?from=" + str(int(new_start)) + '&upto=' + str(int(new_stop))   
                    except ConnectionError:
                        raise ConnectionError("Cannot access webapi to read {:s} from {:s}.".format(webapi_name,exp_id))

                start_shift =int(0)
                if coordinates is not None:
                    coord = coordinates[0]
                    if (coord.unit.name == 'Time') and (coord.mode.equidistant):
                            read_range = [float(coord.c_range[0]), float(coord.c_range[1])]
                    else:
                        raise ValueError("Only timeranges may be defined for reading")
                    orig_times = data_setup.time_query.split('=')
                    orig_start = int((orig_times[1].split('&'))[0])
                    # Increasing the time interval a bit to ensure that the required data is in the range
                    new_start = orig_start + read_range[0] * 1000000000 - 100000
                    new_stop = orig_start + read_range[1] * 1000000000 + 100000                     
                    data_setup.time_query = "_signal.json?from=" + str(int(new_start)) + '&upto=' + str(int(new_stop))
                                                         
                # Downsampling the data if requested
                if _options['Downsample'] is not None:
                    warnings.warn('Downsampling requests do not seem to work properly on the W7-X webapi')
                    data_setup.downsamp_gen(_options['Downsample'])

                if (_options['Verbose']):
                    print("Reading '{:s}'".format(webapi_name),flush=True)
                
                # Reading data
                webapi_data = data_setup.archive_pull()
                if (webapi_data is None):
                    raise ValueError("Cannot read data from webapi:{:s}".format(webapi_name))
                if (_options['Verbose']):
                    print("Read {:d} samples.".format(len(webapi_data['values'])),flush=True)
                if (webapi_data['dimensionCount'] != 1):
                    raise NotImplementedError("Only one-dimensional data read is supported now. This data has {:d} dimensions.".format(webapi_data['dimensionCount']))
                this_data = np.array(webapi_data['values'])
                this_time = np.array(webapi_data['dimensions'])
                this_unit = webapi_data['unit']

            if (len(webapi_data_list) == 0):
                # If this is the first signal in this signal storing common parameters of the virtual signal
                if (exp_id is not None):
                    ref_time = shot_ref_time
                else:
                    ref_time = 0
                webapi_data_list.append(this_data)
                data_ndim = this_data.ndim
                data_shape = this_data.shape
                data_time = this_time
                data_unit = this_unit
            else:
                # If this is not the first signal in the virtual signal checking if common parameters agree with other signals
                if (data_ndim != this_data.ndim):
                    raise TypeError("Multiple webapi data cannot be combined into one DataObject: different number of dimensions.")
                for i_d in range(data_ndim):
                    if (data_shape[i_d] != this_data.shape[i_d]):                       
                        raise TypeError("Multiple webapi data cannot be combined into one DataObject: different data shape.")
                if (data_unit != this_unit):
                    raise TypeError("Multiple webapi data cannot be combined into one DataObject: different data unit.")
                if (len(np.nonzero(data_time != this_time)[0]) != 0):
                    raise TypeError("Multiple webapi data cannot be combined into one DataObject: different timescale.")                        
                webapi_data_list.append(this_data)
                
            if (not data_cached and _options['Cache Data'] and (_options['Cache Directory'] is not None)):
                webapi_data_pickle = [webapi_data_list[-1],data_time,exp_id,shot_ref_time,data_unit]
                while True:
                    if not (os.path.exists(directory)):
                        try:
                            os.mkdir(directory)
                        except:
                            raise SystemError("The shot folder cannot be created. Cache directory might not be present.")
                    try:          
                        f = io.open(filename,"wb")
                    except:
                        print("Warning: Cannot open cache file: "+filename)
                        break
                    try:
                        pickle.dump(webapi_data_pickle,f)
                    except Exception as e:
                        print("Warning: Cannot write cache file: "+filename)
                        break
                    try:
                        f.close()
                        del webapi_data_pickle
                    except Exception as e:
                        print("Warning: Cannot close cache file: "+filename)
                    break
        # End of collecting all data for a data_name constituents
        
        if (len(signal_list) == 0):
            # This is the first data_name processed
            common_time = data_time
            common_dtype = webapi_data_list[0].dtype
            common_unit = data_unit
        else:
            if (len(common_time) != len (data_time)):
                raise ValueError("Different data length for signals. Try renewing cache.")
            if (len(np.nonzero(np.absolute(common_time - data_time) \
                               > ((common_time[1] - common_time[0]) / len(common_time + 1)) * 0.5
                               )[0]
                    ) != 0
                ):
                raise ValueError("Signals have different timescale. Cannot construct one flap.DataObject.")
            if (common_dtype != webapi_data_list[0].dtype):
                raise ValueError("Signals have different data type (complex, float). Cannot construct one flap.DataObject.")
            if (common_unit != data_unit):
                raise ValueError("Signals have different data unit. Cannot construct one flap.DataObject.")
        if (complex_data):
            signal_list.append(webapi_data_list[0] + 1j * webapi_data_list[1])
        else:
            signal_list.append(webapi_data_list[0])
        if (data_description.invert is not None):
            if (data_description.invert[i_signal]):
                signal_list[-1] *= -1
               


    if (_options['Check Time Equidistant']):
        mean_tstep = (data_time[-1] - data_time[0]) / (len(data_time) - 1)
        if (len(np.nonzero(np.absolute(np.diff(data_time) - mean_tstep) > mean_tstep / 100)[0]) != 0):
            # Time scale is non-equidistant
            coord_time = flap.Coordinate(name='Time',unit='Second',mode=flap.coordinate.CoordinateMode(equidistant=False),
                                         shape=[data_time.shape],values=(data_time-ref_time) / 1e9,dimension_list=[0])
        else:
            # Time scale is equidistant
            coord_time = flap.Coordinate(name='Time',unit='Second',mode=flap.coordinate.CoordinateMode(equidistant=True),
                                        shape=[data_time.shape],start=(data_time[0] - ref_time) / 1E9,step=mean_tstep / 1E9,dimension_list=[0])
    coord_sample = flap.Coordinate(name='Sample',unit='a.u.',mode=flap.coordinate.CoordinateMode(equidistant=True),
                                   shape=[data_time.shape],start=0,step=1,dimension_list=[0])

    if (len(signal_list) == 1):
        coord_signal = flap.Coordinate(name='Signal name',mode=flap.coordinate.CoordinateMode(equidistant=False),shape=[],
                                       values=signal_name_list,dimension_list=[])
        d = np.array(signal_list[0])
        coordinates=[coord_time,coord_sample,coord_signal]
    else:
        d = np.ndarray((len(signal_list[0]),len(signal_list)),dtype=signal_list[0].dtype)
        for i,s in enumerate(signal_list):
            d[:,i] = s
        coord_signal = flap.Coordinate(name='Signal name',mode=flap.coordinate.CoordinateMode(equidistant=False),shape=len(signal_name_list),
                                       values=signal_name_list,dimension_list=[1])
        coord_sn = flap.Coordinate(name='Signal number',mode=flap.coordinate.CoordinateMode(equidistant=True),shape=len(signal_name_list),
                                   start=1,step=1,dimension_list=[1])
        coordinates=[coord_time,coord_sample,coord_signal,coord_sn]
        
    return flap.DataObject(data_array=d,
                           data_unit=flap.Unit(name='Signal',unit=data_unit),
                           coordinates=coordinates,
                           exp_id=exp_id,data_title='W7-X data'
                           )

def get_data(exp_id=None, data_name=None, no_data=False, options={}, coordinates=None, data_source=None):
    """ Data read function for the W7-X webapi database
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
        Time Start and Time Stop : datetime.datetime
            Times can be given as an alternative to exp_id e.g.
            options={'Time Start':datetime.datetime(2018, 8, 22, 10, 46, 53),
                     'Time Stop': datetime.datetime(2018, 8, 22, 10, 47, 8)}
        Downsample : float, default is None
            The sampling of the data to this frequency.
        Scale Time : bool, default is False
            if set and True, then the data is normalized to seconds and 0s point is set to
            the start of the shot
        Check Time Equidistant : bool , default is False
            If set and True, then the data is checked if it is equidistant in time.
            it is assumed to be equidistant, if the difference between the timesteps relatively small
        Cache Data : bool, default is None
            Store a local copy of the data so as next time it is read from there.
        Cache Directory : str
            Cache directory name. If None will store data under the source file directory.
        V1 : bool, default is False
            If True the new version of the program will be used: get_data_v1
    """

    options_default = {'Cache Data': True,
                       'Cache Directory': None,
                       'Scale Time':False,
                       'Time Start': None,
                       'Time Stop': None,
                       'Downsample': None,
                       'Check Time Equidistant': False,
                       'V1' : False
                       }
    try:
        options['V1']
    except KeyError:
        options['V1'] = False
    if (options['V1']):
        return get_data_v1(exp_id=exp_id, data_name=data_name, no_data=no_data, options=options, coordinates=coordinates, data_source=data_source)
            

    _options = flap.config.merge_options(options_default,options,data_source='W7X_WEBAPI')
 
    if (_options['V1']):
        get_data_v1(exp_id=exp_id, data_name=data_name, no_data=no_data, options=options, coordinates=coordinates, data_source=data_source)
        

    note = "" # any notes that are added to the dataobject info
    
    #checking if the data is already in the cache
    curr_path = os.path.dirname(os.path.abspath(__file__))
    location = os.path.sep.join(curr_path.split(os.path.sep))
    source = data_name.lower().split("/")
    if (_options['Cache Directory'] is None):
        filename = os.path.join(os.path.sep.join([location,"cached"]),
                                ('_'.join(source)+"-"+str(exp_id)+".hdf5").lower()
                                )
    else:
        filename = os.path.join(_options['Cache Directory'],
                                ('_'.join(source)+"-"+str(exp_id)+".hdf5").lower()
                                )
    if os.path.exists(filename):
        return flap.load(filename)

    if no_data:
        raise NotImplementedError('no data input is not implemented yet.')

    get_vector = False
    if data_name[0:3] == "TS-":
        #old data was uploaded in a differen manner
        if int(exp_id[:8]) <= 20221001:
            data_setup = archive_signal.ArchiveSignal(data_name+"-data-vol1", exp_id=exp_id)
            get_vector = True
        else:
            # checking if a version number was given, if not, the version number temporarily set to 1 and later the latest version is downloaded
            if len(data_name.split("-")) == 2:
                print("-".join(["TS-v1", data_name.split("-")[1]]))
                get_version = spec_data.ThomsonData("-".join(["TS-v1", data_name.split("-")[1]]), exp_id=exp_id)
                # Generating the data relevant time intervals
                if exp_id is None:
                    if not ((_options['Time Start'] is not None) and (_options['Time Stop'] is not None)):
                        raise ValueError("either exp_id should be given, or options should have a 'Time Start' and 'Time Stop' key.")
                    else:
                        get_version.time_query_gen(_options['Time Start'], _options['Time Stop'])
                else:
                    get_version.shotid_time_query_gen(exp_id)
                version = get_version.list[0].last_version() 
                data_setup = spec_data.ThomsonData("-".join([f"TS-v{version}", data_name.split("-")[1]]), exp_id=exp_id)
            else:
                data_setup = spec_data.ThomsonData(data_name, exp_id=exp_id)
    if data_name == "ECE-te":
        data_setup = spec_data.ECEData(data_name)
    elif data_name[0:5] == "ABES-":
        data_setup = archive_signal.ArchiveSignal(data_name+"-data-vol15", exp_id=exp_id)
        get_vector = True
    else:
        data_setup = archive_signal.ArchiveSignal(data_name, exp_id=exp_id)
    # Generating the data relevant time intervals
    if exp_id is None:
        if not ((_options['Time Start'] is not None) and (_options['Time Stop'] is not None)):
            raise ValueError("either exp_id should be given, or options should have a 'Time Start' and 'Time Stop' key.")
        else:
            data_setup.time_query_gen(_options['Time Start'], _options['Time Stop'])
    else:
        data_setup.shotid_time_query_gen(exp_id)

    start_shift = decimal.Decimal(0)
    if coordinates is not None:
        coord = coordinates[0]
        if (coord.unit.name == 'Time') and (coord.mode.equidistant):
                read_range = [float(coord.c_range[0]), float(coord.c_range[1])]
        else:
            raise ValueError("Only timeranges may be defined for reading")
        orig_times = data_setup.time_query.split('=')
        orig_start = decimal.Decimal((orig_times[1].split('&'))[0])
        if int(exp_id[:8]) <= 20221001:
            start_shift = decimal.Decimal((read_range[0]+1.0)*1e9)
        else:
            start_shift = decimal.Decimal(read_range[0]*1e9)
        new_start = orig_start + start_shift
        if int(exp_id[:8]) <= 20221001:
            new_stop = orig_start + decimal.Decimal((read_range[1]+1.0)*1e9)
        else:
            new_stop = orig_start + decimal.Decimal((read_range[1])*1e9)
        data_setup.time_query = "_signal.json?from=" + str(new_start) + '&upto=' + str(new_stop)

    # Downsampling the data if requested
    if _options['Downsample'] is not None:
        warnings.warn('Downsampling requests do not seem to work properly on the W7-X webapi')
        data_setup.downsamp_gen(_options['Downsample'])

    # Converting data to flap.DataObject format
    if get_vector is True:
        if data_name[0:3] == "TS-":
            # data is expected in the form 'TS-v20-ne_map'
            split_out = data_name.split('-')
            version_data = split_out[0]+"-"+split_out[1]
            signal_type = split_out[2]
            vector_data = spec_data.ThomsonDataOld(version_data, signal_types=[signal_type])
        elif data_name[0:5] == "ABES-":
            split_out = data_name.split('-')
            version_data = split_out[0]+"-"+split_out[1]
            signal_type = split_out[2]
            vector_data = spec_data.ABESData(version_data, signal_types=[signal_type], exp_id=exp_id)
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
        if _options['Scale Time'] is True:
            [x.time_norm() for x in vector_data.data]
            if coordinates is not None:
                coord = coordinates[0]
                if (coord.unit.name == 'Time [s]') and (coord.mode.equidistant):
                    vector_data.data[0].list[0].time = [element + read_range[0] for element in
                                                        vector_data.data[0].list[0].time]
                    for x in vector_data.params:                    
                        x.list[0].time = [element + read_range[0] for element in x.list[0].time]

        d = conv_vector_to_dataobj(vector_data, data_setup.url, data_name, exp_id, _options)
    else:
        # Downloading data
        data = data_setup.archive_pull()
        
        if str(type(data)) != "<class 'list'>":
            if data is None:
             raise IOError(f"Error reading data from webapi at {data_setup.url}")
        else:
            for index in range(len(data)):
                if (data[index] is None):
                    raise IOError(f"Error reading data from webapi at {data_setup.list[index].url}")

        # Normalizing the time
        if _options['Scale Time'] is True:
            if data_name == 'QSI-CXRS':
                note += f"\n{data_name} timing was incorrect during OP2.1, correcting this automatically"
                warnings.warn(note)
                mintime = min(data["dimensions"])-start_shift
                data["dimensions"][:] = [float(i - mintime) * 1.0e-9 for i in data["dimensions"]]
            elif data_name[:3] == "TS-":
                for d in data:
                    mintime = min(d["dimensions"])+1000000000-start_shift
                    d["dimensions"][:] = [float(i - mintime) * 1.0e-9 for i in d["dimensions"]]
            elif data_name == "ECE-te":
                for d in data:
                    mintime = min(d["dimensions"])+1000000000-start_shift
                    d["dimensions"][:] = [float(i - mintime) * 1.0e-9 for i in d["dimensions"]]
            else:
                mintime = min(data["dimensions"])+1000000000-start_shift
                data["dimensions"][:] = [float(i - mintime) * 1.0e-9 for i in data["dimensions"]]
        if data_name[0:4] == "AUG-":
            d = conv_aug2_to_dataobj(data, data_setup.url, data_name, exp_id, _options)
        elif data_name[0:3] == "TS-":
            d = conv_thomson_to_dataobj(data, data_setup.list[0].url, data_name, exp_id, _options)
        elif data_name == "ECE-te":
            d = conv_ece_to_dataobj(data, data_setup.list[0].url, data_name, exp_id, _options)
        else:
            if data['dimensionCount'] == 1:
                d = conv_scalar_to_dataobj(data, data_setup.url, data_name, exp_id, _options)
            else:
                d = conv_matrix_to_dataobj(data, data_setup.url, data_name, exp_id, _options)

    d.info = f"Source: {data_setup.stream}\n"
    d.info += note

    if _options["Cache Data"] is True:
        d = standard_dataobject(d)
        flap.save(d, filename)
    return d


def add_coordinate(data_object, new_coordinates, options=None):
    raise NotImplementedError("Coordinate conversions not implemented yet.")




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
        self.rphiz = np.empty([0, 3])  # point coordinates in stellarator radius-phi-z coordinates
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
        self.rphiz = self.xyz_to_rphiz(self.xyz)

    def append_rphiz(self, rphiz_list):
        '''
        Add points to the object based on their rphiz coordinates. rphiz_list has to be
        either a numpy array of length 3 defining the point location, or a two-dimensional
        numpy array.
        '''
        if rphiz_list.ndim == 1:
            rphiz_list = np.asarray([rphiz_list])
        self.rphiz = np.append(self.rphiz, rphiz_list, axis=0)
        self.xyz = self.rphiz_to_xyz(self.rphiz)

    @staticmethod
    def rphiz_to_xyz(rphiz_list):
        '''
        Convert rphiz coordinates to xyz and returns the value. rphiz_list has to be
        either a numpy array of length 3 defining the point location, or a two-dimensional
        numpy array. Practically list of points.
        '''
        if rphiz_list.ndim == 1:
            rphiz_list = np.asarray([rphiz_list])
        x_list = rphiz_list[:, 0]*np.cos(rphiz_list[:, 1])
        y_list = rphiz_list[:, 0]*np.sin(rphiz_list[:, 1])
        return np.squeeze(np.transpose([x_list, y_list, rphiz_list[:, 2]]))

    @staticmethod
    def xyz_to_rphiz(xyz_list):
        '''
        Convert xyz coordinates to rphiz and returns the value. xyz_list has to be
        either a numpy array of length 3 defining the point location, or a two-dimensional
        numpy array. Practically list of points.
        '''
        if xyz_list.ndim == 1:
            xyz_list = np.asarray([xyz_list])
        r_list = np.linalg.norm(xyz_list[:, 0:2], axis=1)
        t_list = np.arctan2(xyz_list[:, 1], xyz_list[:, 0])
        return [np.squeeze(np.transpose([r_list, t_list, xyz_list[:, 2]]))]

    def xyz_to_reff(self, shotID=None, normalize=True):
        '''
        Obtain the effective radius based on self.xyz of the object and vmec results
        '''
        if shotID:
            self.shotID = shotID
            self.ref_eq = get_refeq(shotID)
        if hasattr(self, "ref_eq") is False:
            self.ref_eq = get_refeq(shotID)
        if self.ref_eq:
            url = 'http://svvmec1.ipp-hgw.mpg.de:8080/vmecrest/v1/'+self.ref_eq+'/reff.json?'
        else:
            info = 'Points3D ' + shotID + 'error: No data at ' + url + "\n" + \
                    "The information at http://svvmec1.ipp-hgw.mpg.de:8080/vmecrest/v1/w7x_ref may help"
            logger.error(info)
            print(info)
        self.reff = np.empty([0])
        info = "Obtaining effective radius of points"
        logger.info(info)
        pointnum = 0
        for point in self.xyz:
            info = f"{int(pointnum/len(self.xyz)*1000)*0.1}%"
            logger.info(info)
            url_str = url + 'x='+str(point[0])+'&y='+str(point[1])+'&z='+str(point[2])
            data = archive_signal.ArchiveSignal.url_request(url_str)
            self.reff = np.append(self.reff, data['reff'][0])
            pointnum += 1
        if normalize is True:
            self.normalize_reff()
    
    def reff_at_lcfs(self, shotID=None):
        if shotID is not None:
            self.shotID = shotID
        r,z = get_lcfs(self.shotID, 72)
        r = r[0]
        z = z[0]
        if hasattr(self, "ref_eq") is False:
            self.ref_eq = get_refeq(self.shotID)
        url = 'http://svvmec1.ipp-hgw.mpg.de:8080/vmecrest/v1/'+self.ref_eq+'/reff.json?'
        x = r*np.cos(72/180*np.pi)
        y = r*np.sin(72/180*np.pi)
        x=int(x*1000)*0.001
        y=int(y*1000)*0.001
        z=int(z*1000)*0.001

        url_str = url + 'x='+str(x)+'&y='+str(y)+'&z='+str(z)
        data = archive_signal.ArchiveSignal.url_request(url_str)
        reff = data['reff'][0]
        return reff

    def normalize_reff(self):
        self.reff = self.reff.astype("float")/self.reff_at_lcfs()


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
    data = archive_signal.ArchiveSignal.url_request(url)
    r = data['lcfs'][0]['x1']
    z = data['lcfs'][0]['x3']
    return r, z




def get_refeq(shotID):
    """
    Used for getting the reference number for a w7x configurations from a shotID
    Input of shotID is a string, output is a sting in the format e.g. w7x_ref_120
    """
    source_stream = 'ArchiveDB/raw/W7XAnalysis/Equilibrium/RefEq_PARLOG/V1/parms/equilibriumID/'
    refeq = archive_signal.ArchiveSignal(source_stream)
    refeq.shotid_time_query_gen(shotID)
    data = refeq.archive_pull()
    ref_eq = data['values'][0]
    return ref_eq
