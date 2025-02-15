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

def webapi_virtual_names(data_name, exp_id, channel_config_file):

    """
    Translates virtual data names to webapi entries.
    Returns a list of virtual names and translated full names.
    The list is read from the configuration file. The file is a standard
    configuration file with one section [Virtual names]. Each entry looks
    like one of the following templates.
        Simple translation:
            <name> = <webapi name>
        Translation in an expID range, including ends:
            <name>(startID-endID) = <webapi name>
          One of startID or endID may be omitted.
        <webapi name> can be either a single name or complex(name1, name2) which
        will construct a complex signal from the two entries. 
        It can also be a list of signals (aa,bb,...)
    Wildcards are allowed in data names. Data names can be a single string or
        list of strings.
    Return value:
        virt_names, webapi_str, webapi_names
        virt names is a list of the interpreted names.
        virt_str is the description from the virtual name translation file
        webapi_names is a list of the same length as virt_names. Elements can be the following:
        Elements can be the following:
            string: a single webapi name
            list: [<type>, <name1>, <name2>, ...]
                type can be
                  'complex': two webapi entries are expected (real, imag), complex
                  signal will be created from them.

    """
    if ((type(exp_id) is not str) or (len(exp_id) != 12)):
        raise TypeError("exp_id should be string in format YYYYMMDD.xxx")
    exp_id_num = exp_id
        
    config = configparser.ConfigParser()
    config.optionxform = str
    read_ok = config.read(channel_config_file)
    if (read_ok == []):
        raise ValueError("Invalid webapi virtual name file "+channel_config_file) 
    try:
        entries = config.items('Virtual names')
    except Exception as e:
        raise ValueError("Invalid webapi virtual name file "+channel_config_file)
    entry_descr = [e[0] for e in entries]
    entry_values = [e[1] for e in entries]
    webapi_names = []
    entry_names = []
    # Finding etries with matching exp_id
    for e,ev in zip(entry_descr, entry_values):
        start_brace = e.find('(')
        stop_brace = e.find(')')
        if (start_brace * stop_brace < 0):
            raise ValueError("Invalid entry '{:s}' in virtual name file: {:s}".format(e, channel_config_file))
        if (start_brace > 0):
            # ExpID range given
            exp_id_range = e[start_brace+1:stop_brace].split('-')
            if (len(exp_id_range) != 2):
                raise ValueError("Invalid exp_id range in entry '{:s}' in virtual name file: {:s}".format(e, channel_config_file))
            if (exp_id_range[0].strip() == ''):
                exp_id_start = None
            else:
                try:
                    exp_id_start = int(exp_id_range[0])
                except ValueError:
                    raise ValueError("Invalid exp_id start in entry '{:s}' in virtual name file: {:s}".format(e, channel_config_file))
            if (exp_id_range[1].strip() == ''):
                exp_id_stop = None
            else:
                try:
                    exp_id_stop = int(exp_id_range[1])
                except ValueError:
                    raise ValueError("Invalid exp_id stop in entry '{:s}' in virtual name file: {:s}".format(e, channel_config_file))
            if ((exp_id_start is not None) and (exp_id_num < exp_id_start) or \
                (exp_id_stop is not None) and (exp_id_num > exp_id_stop)) :
                continue
            entry_names.append(e[:start_brace])
        else:
            entry_names.append(e)
        webapi_names.append(ev)
        
    if (type(data_name ) is not list):
        _data_name = [data_name]
    else:
        _data_name = data_name
    select_list = []
    select_webapi_list = []
    for i,dn in enumerate(_data_name):    
        try:
            sl, si = flap.select_signals(entry_names, dn)
            select_list += sl
            select_webapi_list += [webapi_names[i] for i in si]
        except ValueError as e:
            select_list.append(None)
            select_webapi_list.append(dn)
            
    webapi_descr = []
    for descr in select_webapi_list:
        start_brace = descr.find('(')
        stop_brace = descr.find(')')
        if (start_brace * stop_brace < 0):
            raise ValueError("Invalid value '{:s}' in virtual name file: {:s}".format(descr, channel_config_file))
        if (start_brace >= 0):
            webapi_list = descr[start_brace+1:stop_brace].split(',')
            webapi_type = descr[:start_brace]
            webapi_descr.append([webapi_type] + webapi_list)
        else:
            webapi_descr.append(descr)    
    return select_list, select_webapi_list, webapi_descr        

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
        V1 : bool, default is False
            If True the new version of the program will be used: get_data_v1
        Verbose : bool
            Print diagnostic messages.
    """

    options_default = {'Cache Data': True,
                       'Cache Directory': None,
                       'Scale Time':False,
                       'Time Start': None,
                       'Time Stop': None,
                       'Downsample': None,
                       'Check Time Equidistant': False,
                       'V1' : False,
                       'Virtual name file' : 'Virtual_names_webapi.cfg',
                       'Verbose' : False
                       }
  
    _options = flap.config.merge_options(options_default,options,data_source='W7X_WEBAPI')
 
    if ((type(exp_id) is not str) or (len(exp_id) != 12)):
        raise TypeError("exp_id should be string in format YYYYMMDD.xxx")
    
    if (_options['Virtual name file'] is not None):
        try:
            virt_names, virt_webapi_txt, virt_webapi = webapi_virtual_names(data_name, exp_id, _options['Virtual name file'])
        except Exception as e:
            raise e
    
    # Assembling a list of signals and list of webapi request
    signal_list = []
    common_time = None    
    for name, webapi_descr in zip(virt_names,virt_webapi):
        # Assembling a list of webapi nodes needed for this data
        webapi_request_list = []
        if (name is None):
            # This was not recognized as virtual name
            webapi_request_list = [webapi_descr]
            signal_list.append(webapi_descr)
            readtype = 0
        else:
            # This was recongnized as a virtual signal
            signal_list.append(name)
            if (type(webapi_descr) is not list):
                # This was recognized as a single webapii node
                webapi_request_list = [webapi_descr]
                readtype = 0
            else:
                # This is a composite virtual signal
                webapi_request_list = webapi_descr[1:]
                readtype = 1 
        
        # Processing reqgular expressions
        webapi_request_list = flap.select_signals(None,webapi_request_list)[0]
        
        for webapi_name in webapi_request_list: 
 
            # Trying to read cached data
            if ((_options['Cache Data']) and (_options['Cache Directory'] != None)):
                filename = str(exp_id)+'_'+webapi_name
                for c in ['\\',':']:
                   filename = filename.replace(c,'_')
                #filename = os.path.join(_options['Cache directory'],filename+'.pickle')
                directory=os.path.join(_options['Cache directory'],str(exp_id))
                if not (os.path.exists(directory)):
                    try:
                        os.mkdir(directory)
                    except:
                        raise SystemError("The shot folder cannot be created. Cache directory might not be present.")
                                      
                filename = os.path.join(directory,filename+'.pickle')
                try:
                    f = io.open(filename,'rb')
                    webapi_pickle = pickle.load(f)
                    f.close()
                    try:
                        print("Get data from pickle file.")
                        del webapi_pickle
                        data_cached = True
                    except:
                        data_cached = False
                except:
                    data_cached = False
            else:
                data_cached = False
                
            if (not data_cached):
                webapi_data = [1,2,3]
                
                # data_setup = archive_signal.ArchiveSignal(webapi_name, exp_id=exp_id)
                # # Generating the data relevant time intervals
                # if exp_id is None:
                #     if not ((_options['Time Start'] is not None) and (_options['Time Stop'] is not None)):
                #         raise ValueError("either exp_id should be given, or options should have a 'Time Start' and 'Time Stop' key.")
                #     else:
                #         data_setup.time_query_gen(_options['Time Start'], _options['Time Stop'])
                # else:
                #     data_setup.shotid_time_query_gen(exp_id)
     
                # start_shift = decimal.Decimal(0)
                # if coordinates is not None:
                #     coord = coordinates[0]
                #     if (coord.unit.name == 'Time') and (coord.mode.equidistant):
                #             read_range = [float(coord.c_range[0]), float(coord.c_range[1])]
                #     else:
                #         raise ValueError("Only timeranges may be defined for reading")
                #     orig_times = data_setup.time_query.split('=')
                #     orig_start = decimal.Decimal((orig_times[1].split('&'))[0])
                #     if int(exp_id[:8]) <= 20221001:
                #         start_shift = decimal.Decimal((read_range[0]+1.0)*1e9)
                #     else:
                #         start_shift = decimal.Decimal(read_range[0]*1e9)
                #     new_start = orig_start + start_shift
                #     if int(exp_id[:8]) <= 20221001:
                #         new_stop = orig_start + decimal.Decimal((read_range[1]+1.0)*1e9)
                #     else:
                #         new_stop = orig_start + decimal.Decimal((read_range[1])*1e9)
                #     data_setup.time_query = "_signal.json?from=" + str(new_start) + '&upto=' + str(new_stop)
     
                # # Downsampling the data if requested
                # if _options['Downsample'] is not None:
                #     warnings.warn('Downsampling requests do not seem to work properly on the W7-X webapi')
                #     data_setup.downsamp_gen(_options['Downsample'])

                if (_options['Verbose']):
                    print("Reading {:s}".format(webapi_name),flush=True)

                webapi_data = [1,2,3]

                if (not data_cached and (_options['Cache Data']) and (_options['Cache Directory'] is not None)):
                    webapi_data_pickle = [webapi_data]
                    while True:
                        try:          
                            f = io.open(filename,"wb")
                        except:
                            print("Warning: Cannot open cache file: "+filename)
                            break
                        try:
                            pickle.dump(webapi_data_pickle,f)
                            del webapi_data_pickle
                        except Exception as e:
                            print("Warning: Cannot write cache file: "+filename)
                            break
                        try:
                            f.close()
                        except Exception as e:
                            print("Warning: Cannot write cache file: "+filename)
                        break

        print("webapi_descr[0]:'{:s}'".format(webapi_descr[0]))


    return flap.DataObject(data_array=np.array(webapi_data))

    

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
        #old data was uploadid in a differen manner
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
