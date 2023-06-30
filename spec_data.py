#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 14:30:16 2023

@author: apdcam
"""
import os
import sys  # For error messages
import warnings
import decimal
import json
import urllib
import numpy as np
from . import archive_signal

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

    def __init__(self, signal_names, exp_id="20221201.001"):
            ""
            self.list = [archive_signal.ArchiveSignal(signal, exp_id=exp_id) for signal in signal_names]

    def __getattr__(self, name):
        """
        When we try to access some function that doesn't exist in the class
        definition, treat it as a class method and call it on every signal element.
        """
        def method(*args):
            out = []
            for signal in self.list:
                method_to_call = getattr(signal, name)
                out += [method_to_call(*args)]
            return out
        return method


class ScalarDataList():
    """

    This class handles multiple ScalarData objects

    """
    def __init__(self, data_list, from_time=None, data_source_obj=None):
        self.list = []
        for data in data_list:
            if data is not None:
                self.list.append(ScalarData(data, from_time=from_time, data_source_obj=data_source_obj))

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

    def __init__(self, signal_id, signal_types=None, param_types=None, volumes=None,
                 exp_id="20221201.001"):
        """
        By using the SignalList class, the necessary structure for reading the
        data and the parameters are defined
        """
        self.signals = []
        self.params = []
        self.signal_types = []
        self.param_types = []
        if not (signal_types is None):
            self.signal_types = signal_types
            signal_names = [[signal_id+'-'+signal_type+'-data-vol'+str(vol_index) for vol_index in volumes]
                            for signal_type in self.signal_types]
            self.signals = [SignalList(signal, exp_id=exp_id) for signal in signal_names]
        if not (param_types is None):
            self.param_types = param_types
            param_names = [[signal_id+'-'+param_type+'-vol'+str(vol_index) for vol_index in volumes]
                           for param_type in self.param_types]
            self.params = [SignalList(param, exp_id=exp_id) for param in param_names]

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



class ThomsonDataOld(VectorData):
        """"
        Data before OP2.1: The Thomson scattering data is uploaded in a very complciated manner:
            The output of every data of every mesurement point in space is in a
            different datastream. E.g. The time evolution of the electron temperature
            at a point in space is given by one unique datastream. Moreover, the location
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


class ThomsonData(SignalList):
        """"
        
        """
        
        def __init__(self, signal_id, exp_id="20221201.001"):
            # signal_id is expected to be in the form 'TS-v1-ne or TS-v1-ne, where v1 is version 1'
            uncertainty = f"{signal_id}_std"
            additional_data = ["rho","R","Phi","Z"]
            additional_signals = ["-".join(signal_id.split("-")[:2]+[data]) for data in additional_data]
            signal_names = [signal_id, uncertainty]+additional_signals
            self.list  = [archive_signal.ArchiveSignal(signal, exp_id=exp_id) for signal in signal_names]

class ECEData(SignalList):
        """"
        
        """
        
        def __init__(self, signal_id, exp_id="20221201.001"):
            te = [f"te-ch{i:02d}" for i in np.arange(32)+1]
            get_data = ["rho"]+te
            signal_names = ["-".join(["ECE"]+[data]) for data in get_data]
            self.list  = [archive_signal.ArchiveSignal(signal, exp_id=exp_id) for signal in signal_names]

class ABESData(VectorData):
        """"
        Obtains the ABESData from the archive.
        For more info look up the description in https://w7x-logbook.ipp-hgw.mpg.de/components?id=QSI#
        This class is made to simplify the processing of data
        """
        
        def __init__(self, signal_id, signal_types=None, volumes=None, exp_id="20221201.001"):
            # signal_id is expected to be in the form 'TS-v20'
            if signal_types is None:
                signal_types = ['density', 'density_error_low', 'density_error_high', 'meas_light', 'meas_light_error',
                                'x_coord', 'y_coord']

            param_types = None # ['x_coord-data', 'y_coord-data']

            if volumes is None:
                volumes = np.linspace(1, 30, num=30, dtype=int)

            VectorData.__init__(self, signal_id, signal_types=signal_types, param_types=param_types,
                                volumes=volumes, exp_id=exp_id)
            self.data = []
            self.params = []

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

# ----------------------------------------------------------------------------------------------------------------------
# ---------------------------------------Uploading via webapi-----------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

class WriteABESSignal(archive_signal.ArchiveSignal):
    def __init__(self, datalist, shotid, minindex=0, version=1.0):
        self.dens = datalist[0]
        self.light_meas = datalist[1]
#        self.light_recon = datalist[2]
        self.shotid = shotid
        self.minindex = minindex
        self.version = int(version)
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
                    #parameters
                    url_loc = url + key+"_PARLOG/V"+str(self.version)
                    req = urllib.request.Request(url_loc)
                    req.add_header('Content-Type', 'application/json; charset=utf-8')
                    urllib.request.urlopen(req, self.json["params"])
                    #data
                    url_loc = url + key+"_DATASTREAM/V"+str(self.version)
                    req = urllib.request.Request(url_loc)
                    req.add_header('Content-Type', 'application/json; charset=utf-8')
                    urllib.request.urlopen(req, self.json["data"][key])
                except urllib.error.HTTPError as e:
                    print(e)
                    print(url_loc)
                    print(e.headers)
