#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:22:43 2023

@author: apdcam
"""

import flap_w7x_webapi as webapi
import flap
from matplotlib import pyplot as plt

flap.register_data_source('W7X_WEBAPI',
                          get_data_func=webapi.get_data,
                          add_coord_func=webapi.add_coordinate)

# laser = "ArchiveDB/raw/Minerva/Minerva.Thomson.Lasers/laser_DATASTREAM/V3/0/laser/"
# th_laser = flap.get_data('W7X_WEBAPI', name=laser,
#                           exp_id="20230323.055",
#                           options={'Scale Time': True,
#                                    'Cache Data': True},
#                           object_name='laser ID',
#                           coordinates={'Time': [0, 10]})
# th_laser.plot()

shotID = "20230323.063"
Te_s="Test/raw/Minerva2/Minerva.ThomsonScattering.Profiles/Preliminary_Fast_DATASTREAM/V1/2/te/"
th_laser = flap.get_data('W7X_WEBAPI', name=Te_s,
                          exp_id=shotID,
                          options={'Scale Time': True,
                                   'Cache Data': True},
                          object_name='temp data',
                          coordinates={'Time': [0, 10]})
plt.subplot(1,2,1)
plt.imshow(th_laser.data)
plt.subplot(1,2,2)
plt.plot(th_laser.data[50,:])

dTe_s=r"Test/raw/Minerva2/Minerva.ThomsonScattering.Profiles/Preliminary_Fast_DATASTREAM/V1/4/te_std/"