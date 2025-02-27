# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 12:21:28 2025

@author: zoletnik
"""

def w7x_summary(exp_ID,datapath='https://w7x-logbook.ipp-hgw.mpg.de/api//log/history/XP/'):

    try:
        import requests, json
    except ModuleNotFoundError:
        raise RuntimeError("Not on W7-X machine?")

#base = 'https://w7x-logbook.ipp-hgw.mpg.de/api'
#path = '/log/history/XP/XP_20180718.22'


    path = datapath+'XP_'+exp_ID[:9]+"{:d}".format(int(exp_ID[9:]))
#    print(path)
    res = requests.get(path)

    if res.status_code == 200:
      r = res.json()
      for c in r['hits']['hits'][-1]['_source']['tags']:
          if (c['description'] == 'Configuration main coils'):
              print(c['Main field']) 

w7x_summary('20181018.003')