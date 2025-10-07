# -*- coding: utf-8 -*-
"""
Created on Tue Oct  7 15:30:43 2025

@author: Zoletnik
"""
import matplotlib.pyplot as plt

import flap_w7x_webapi

retval = flap_w7x_webapi.get_overview_data('20250311.079',mean=True)
print(retval)
retval = flap_w7x_webapi.get_overview_data('20250311.079',mean=False)
plt.close('all')
plt.plot(retval['nedl']['t'],retval['nedl']['val'])
plt.plot([retval['Shot timerange'][0]]*2,plt.ylim(),color='red')
plt.plot([retval['Shot timerange'][1]]*2,plt.ylim(),color='red')