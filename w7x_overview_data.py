# -*- coding: utf-8 -*-
"""
Created on Tue Oct  7 12:31:57 2025

Determine plasma global parameters from W7-X HDF5 file
@author: Zoletnik
"""

import os
import h5py
import numpy as np

import flap


def get_overview_data(exp_id,options=None,mean=False):
    """
    Read shot overview data from the HDF file of the global plot.

    Parameters
    ----------
    exp_id : str
        The experiment ID, YYYYMMDD.XXXX.
    options : dict, optional
        Dictionary of options. The default is None.
        Defaults can also be entered in the flap_default.cfg file under section [W7-X_Overview]
        Possible options:
            'HDF path': The path where the files are. Under this path there must be directories
                       with names YYMMDD. Under the shot directories the HDF files of that day
                       can be either directly or under a w7x subdir.
            'Parameters': list of parameters to read. The names must be identical to the 
                           HDF file entry names. Default is
                           ['ECRH','ICRH','NBI','Prad','Wdia','nedl','ECE_core']
                           (nedl is translated to interferometer)
            
    mean : bool, optional
        If Ture returns the mean value in the time interval when Wdia is above 0.5. The default is False.

    Raises
    ------
    TypeError
        Bad exp_id type, not string.
    ValueError
        Invalid exp_id format.

    Returns
    -------
    retval : dict
        A dictionary with the data.
        keys:
            'Shot timerange': list with two times in second.
            other keys refer to the parameters requested
        If mean is True then the values for the keys give the mean value in the shot timerange.
        If mean is False then the values are a dict with 't' as time vector and 'val' as value.

    """
    
    default_options = {'HDF path': 'W7-X_Overview',
                       'Parameters': ['ECRH','ICRH','NBI','Prad','Wdia','nedl','ECE_core']
                       }
    if (type(exp_id) is not str):
        raise TypeError("exp_id should be a string.")
    if ((len(exp_id) != 12) or (exp_id[8] != '.')):
        raise ValueError("Invalid exp_id. Should be YYYYMMDD.xxx")
                      
    _options = flap.config.merge_options(default_options,options,section='W7-X_Overview')
    
    day_dir = os.path.join(_options['HDF path'],exp_id[2:8],'w7x')
    filename = exp_id[2:8] + exp_id[9:]+'.h5'
    full_filename = os.path.join(day_dir,filename)
    try:
        datafile = h5py.File(full_filename, 'r')
    except Exception:
        day_dir = os.path.join(_options['HDF path'],exp_id[2:8])
        full_filename = os.path.join(day_dir,filename)
        try:
            datafile = h5py.File(full_filename, 'r')
        except Exception as e:
            raise ValueError("Error opening file: {:s}. ({:s})".format(full_filename,str(e)))

    retval = {}
    t = np.array(datafile["/traces/Wdia/dims"])
    wdia = np.array(datafile["/traces/Wdia/vals"])
    ind = np.nonzero(wdia > 0.5)[0]
    if (len(ind) != 0):
        shot_timerange = [t[ind[0]],t[ind[-1]]]
    else:
        shot_timerange = [0,0]
    retval['Shot timerange'] = shot_timerange
    for para in _options['Parameters']:
        if (para == 'nedl'):
            _para = 'interferometer'
        else:
            _para= para
        try:
            t = np.array(datafile["/traces/"+_para+"/dims"])
            y = np.array(datafile["/traces/"+_para+"/vals"])
        except KeyError:
            raise ValueError("Parameter '{:s}' not found in HDF file.".format(para))
        if (mean):
            ind = np.nonzero(np.logical_and(t > shot_timerange[0],
                                            t < shot_timerange[1]
                                            )
                             )[0]
            retval[para] = np.mean(y[ind])
        else:
            retval[para] = {'t':t,'val':y}
    return retval
        
    

    
    
    