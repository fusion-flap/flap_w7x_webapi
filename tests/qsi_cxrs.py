#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 15:29:55 2023

@author: apdcam
"""

import os
import unittest
import datetime
import numpy as np
import flap
import flap_w7x_webapi as webapi

flap.register_data_source('W7X_WEBAPI',
                          get_data_func=webapi.get_data,
                          add_coord_func=webapi.add_coordinate)

# get data by shotID
# qsi_cxrs = flap.get_data('W7X_WEBAPI', name='Test/raw/W7X/QSI/spec_interface_test4_DATASTREAM/0/Images/',
#                          exp_id='20230222.005',
#                          options={'Scale Time': True,
#                                   'Cache Data': True},
#                          object_name='QSI_CXRS_data',
#                          coordinates={'Time': [2, 3]})

#get data by time interval
qsi_cxrs = flap.get_data('W7X_WEBAPI', name='Test/raw/W7X/QSI/spec_interface_test4_DATASTREAM/0/Images/',
                          options={'Scale Time': True,
                                  'Cache Data': True,
                                  'Time Start':datetime.datetime(2023, 2, 21, 13, 42, 00),
                                  'Time Stop': datetime.datetime(2023, 2, 21, 13, 43, 00)},
                          object_name='QSI_CXRS_data')

#slice by ROI:
ROI0 = qsi_cxrs.slice_data(slicing={"Coord 1" : 0})
ROI0.data.shape
#slice by time: (if the data is not acquired by exp_id then time will be in ns)
mintime = ROI0.coordinate("Time")[0][0,0]
ROI0_0 = ROI0.slice_data(slicing={"Time" : mintime})
ROI0_0.plot(axes="Coord 2")