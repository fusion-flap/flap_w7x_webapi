#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:22:43 2023

@author: apdcam
"""

import flap_w7x_webapi as webapi
import flap
from matplotlib import pyplot as plt


class ReadThomsonProfiles():
    def __init__(self, shotID):
        flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
        self.shotID = shotID
        self.searchpath=r"Test/raw/Minerva2/Minerva.ThomsonScattering.Profiles/Preliminary_Fast_DATASTREAM/V1/"
        # self.datanamearr=["rho","ne","te","ne_std","te_std","R","Phi","Z","Volume Number"]
        self.datanamearr=["rho","ne","te","ne_std","te_std","R","Phi","Z"]
        self.d=self.ReadData()
        self.ReorganizeData()
        
    def ReadData(self):
        data_arr={}
        for i,dataname in zip(range(len(self.datanamearr)),self.datanamearr):
            name=self.searchpath+str(i)+'/'+dataname+'/'
            data_arr[dataname]= flap.get_data('W7X_WEBAPI', name=name,exp_id=self.shotID)
        return data_arr
    
    def ReorganizeData(self):
        time=self.d["te"].coordinates[0]
        coordinate_names=["rho","R","Phi","Z"]
        coordinates=[]
        coordinates.append(time)
        for coordinate_name in coordinate_names:
            coordinates.append(flap.Coordinate(name=self.d[coordinate_name].data_unit.name,
                                unit=self.d[coordinate_name].data_unit.unit,
                                shape=self.d[coordinate_name].shape,
                                values=self.d[coordinate_name].data,
                                mode=flap.coordinate.CoordinateMode(equidistant=False)))
        
            
        self.Te=flap.DataObject(exp_id=self.shotID,
                        data_array=self.d["te"].data,
                        error=self.d["te_std"].data,
                        data_unit=self.d["te"].data_unit,
                        coordinates=coordinates)
        
        self.ne=flap.DataObject(exp_id=self.shotID,
                        data_array=self.d["ne"].data,
                        error=self.d["ne_std"].data,
                        data_unit=self.d["ne"].data_unit,
                        coordinates=coordinates)
                        # self.dataarr["te"].coordinate("Time")[0][:,0]       
        
if __name__=="__main__":
    shotID = "20230323.063"
    D=ReadThomsonProfiles(shotID)

def test():
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
    Te_s=r"Test/raw/Minerva2/Minerva.ThomsonScattering.Profiles/Preliminary_Fast_DATASTREAM/V1/2/te/"
    dTe_s=r"Test/raw/Minerva2/Minerva.ThomsonScattering.Profiles/Preliminary_Fast_DATASTREAM/V1/4/te_std/"
    ne_stupid_s=r"Test/raw/Minerva2/Minerva.ThomsonScattering.Profiles/Preliminary_Fast_DATASTREAM/V1/1/ne/"
    R_s=r"Test/raw/Minerva2/Minerva.ThomsonScattering.Profiles/Preliminary_Fast_DATASTREAM/V1/5/R/"
    Rho_s=r"Test/raw/Minerva2/Minerva.ThomsonScattering.Profiles/Preliminary_Fast_DATASTREAM/V1/0/rho/"
    V_s=r"Test/raw/Minerva2/Minerva.ThomsonScattering.Profiles/Preliminary_Fast_DATASTREAM/V1/8/Volume Number/"
    Te = flap.get_data('W7X_WEBAPI', name=Te_s, exp_id=shotID,
                            options={'Scale Time': True,
                                     'Cache Data': True},
                            object_name='temp data',
                            coordinates={'Time': [0, 10]})
    dTe = flap.get_data('W7X_WEBAPI', name=dTe_s,
                              exp_id=shotID,
                              options={'Scale Time': True,
                                       'Cache Data': True},
                              object_name='temp data',
                              coordinates={'Time': [0, 10]})
    R = flap.get_data('W7X_WEBAPI', name=R_s,
                              exp_id=shotID,
                              options={'Scale Time': True,
                                       'Cache Data': True},
                              object_name='temp data',
                              coordinates={'Time': [0, 10]})
    Rho = flap.get_data('W7X_WEBAPI', name=Rho_s,
                              exp_id=shotID,
                              options={'Scale Time': True,
                                       'Cache Data': True},
                              object_name='temp data',
                              coordinates={'Time': [0, 10]})
    
    plt.subplot(1,2,1)
    plt.imshow(Te.data)
    plt.subplot(1,2,2)
    plt.plot(Rho.data[50,:],Te.data[50,:])
    
