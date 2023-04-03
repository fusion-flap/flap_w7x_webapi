import os
import unittest
import datetime
import numpy as np
import flap
import flap_w7x_webapi as webapi

# -------------------------------------INDEPENDENT FUNCTIONS FROM FLAP--------------------------------------------------


def test_basic():
    # Getting ECRH data by giving the exact time
    ecrh = webapi.ArchiveSignal('ECRH')
    start = datetime.datetime(2023, 3, 30, 10, 10, 00)
    end = datetime.datetime(2023, 3, 30, 10, 30, 00)
    ecrh.time_query_gen(start, end)
    ecrh_data = ecrh.archive_pull()
    return int(np.mean(ecrh_data["values"]))


def test_shotid():
    # Getting NBI data with shotid
    nbi7 = webapi.ArchiveSignal('NBI-7')
    nbi7.shotid_time_query_gen('20180912.012')
    nbi7.downsamp_gen(512)
    nbi7_data = webapi.ScalarData(nbi7.archive_pull())
    return [int(np.mean(nbi7_data.val)*1000000), int(min(nbi7_data.time)/1e9)]

def test_op1_shotID():
    # Getting ECRH data with shotid
    shotID = '20180912.012'
    ecrh = webapi.ArchiveSignal('ECRH', exp_id=shotID)
    ecrh.shotid_time_query_gen(shotID)
    ecrh.downsamp_gen(512)
    ecrh_data = webapi.ScalarData(ecrh.archive_pull())
    return [int(np.mean(ecrh_data.val)*1000000), int(min(ecrh_data.time)/1e9)]

def test_get_last_dataversion():
    shotID = '20230216.067'
    ecrh = webapi.ArchiveSignal('ECRH', exp_id=shotID)
    ecrh.shotid_time_query_gen(shotID)
    return ecrh.last_version()

def test_streamlink():
    # Getting ECRH data via the streaming link, can be used for arbitrary data
    ecrh2 = webapi.ArchiveSignal('ArchiveDB/codac/W7X/CBG_ECRH/TotalPower_DATASTREAM/V1/0/Ptot_ECRH/scaled/')
    ecrh2.shotid_time_query_gen('20180912.012')
    ecrh2.downsamp_gen(512)
    ecrh2_data = webapi.ScalarData(ecrh2.archive_pull())
    return True


class StandaloneTest(unittest.TestCase):
    def test_basic(self):
        self.assertEqual(test_basic(), 871)

    def test_shotid(self):
        self.assertEqual(test_shotid(), [-26, 1536744303])

    def test_streamlink(self):
        self.assertTrue(test_streamlink())
    
    def test_op1_shotid(self):
        self.assertEqual(test_op1_shotID(), [1245711860, 1536744303])
        
    def test_get_last_dataversion(self):
        self.assertEqual(test_get_last_dataversion(), 5)

# -------------------------------------INTERFACE WITH FLAP----------------------------------------------------------


def test_register():
    flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
    return 'W7X_WEBAPI' in flap.list_data_sources()


def test_reading():
    import time
    starttime = time.time()
    flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
    d = flap.get_data('W7X_WEBAPI', name='ECRH',
                      exp_id='20181016.032',
                      options={'Scale Time': True,
                               'Check Time Equidistant': True,
                               'Cache Data': False},
                      object_name='ECRH_data')
    print(f"\nDownload time for ECRH data for a single shot: {time.time()-starttime}")
    return [d.coordinates[0].mode.equidistant, d.coordinates[0].start]


def test_downsampling():
    flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
    d = flap.get_data('W7X_WEBAPI', name='NBI-7',
                      exp_id='20181016.037',
                      options={'Scale Time': True,
                               'Downsample': 1024,
                               'Cache Data': False},
                      object_name='ECRH_data')
    return d.coordinates[0].shape


def test_timerange():
    flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
    d = flap.get_data('W7X_WEBAPI', name='ArchiveDB/codac/W7X/CBG_ECRH/TotalPower_DATASTREAM/V1/0/Ptot_ECRH/scaled/',
                      exp_id='20181016.037',
                      options={'Scale Time': True,
                               'Cache Data': False},
                      object_name='ECRH_data',
                      coordinates={'Time': [2, 3]})
    return [d.coordinates[0].values[0], d.coordinates[0].values[-1]]


def test_imagedata():
    flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
    d = flap.get_data('W7X_WEBAPI', name='QSI-CXRS',
                      exp_id='20230330.028',
                      options={'Scale Time': True,
                               'Cache Data': False},
                      object_name='QSI CXRS data',
                      coordinates={'Time': [2, 3]})
    return int(np.mean(d.data)*100)


def test_vectordata():
    flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
    d = flap.get_data('W7X_WEBAPI', name='ABES-v2-density',
                      exp_id='20181016.037',
                      options={'Scale Time': True,
                               'Cache Data': False},
                      object_name='ECRH_data',
                      coordinates={'Time': [0,3]})
    return np.mean(d.data)

def test_cache():
    cached_file = os.sep.join(webapi.__file__.split(os.sep)[:-1] +\
                  ['cached','archivedb_codac_w7x_cbg_ecrh_totalpower_datastream'+
                   '_v1_0_ptot_ecrh_scaled_-20181016.037.hdf5'])
    existed = os.path.exists(cached_file)
    if existed is True:
        print("Can't properly test data caching, file already exists")
    flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
    d = flap.get_data('W7X_WEBAPI', name='ArchiveDB/codac/W7X/CBG_ECRH/TotalPower_DATASTREAM/V1/0/Ptot_ECRH/scaled/',
                      exp_id='20181016.037',
                      options={'Scale Time': True,
                               'Cache Data': True},
                      object_name='ECRH_data',
                      coordinates={'Time': [2, 3]})
    exists = os.path.exists(cached_file)
    if not existed:
        os.remove(cached_file)
    return exists

class FLAPTest(unittest.TestCase):
    def test_register(self):
        self.assertTrue(test_register())

    def test_reading(self):
        self.assertEqual(test_reading(), [True, -1.0])

    # def test_downsampling(self):
    #     self.assertEqual(test_downsampling(), [841])

    def test_timerange(self):
        self.assertEqual(test_timerange(),  [2.0, 2.99996])

    def test_imagedata(self):
        self.assertEqual(test_imagedata(),  60308)

    def test_vectordata(self):
        self.assertEqual(test_vectordata(),  2.614503530033091e+19)
        
    def test_cache(self):
        self.assertTrue(test_cache())

#--------------------------------------------POINTS3D check-------------------------------------------------------------
def test_reff():
    point = webapi.Points3D()
    point.append_xyz(np.array([0,6.0,0]))
    point.xyz_to_reff('20181018.003')
    return point.reff


class Points3DTest(unittest.TestCase):
    def test_vectordata(self):
        self.assertEqual(test_reff(),  0.45680578471016137)


#---------------------------------------TESTING SPECIFIC SIGNALS--------------------------------------

def test_thomson():
    shotID = "20230330.059"
    d = flap.get_data('W7X_WEBAPI', name="TS-v1-ne", # if v1 is the latest data version, this is equal to setting "TS-ne" as name
                      exp_id=shotID,
                      options={'Scale Time': True,
                               'Check Time Equidistant': True,
                               'Cache Data': False},
                      object_name='Thomson')
    return int(np.mean(d.slice_data(slicing={"Time":4}).data)*100)
    
class SignalTest(unittest.TestCase):
    def test_thomson(self):
        self.assertEqual(test_thomson(),  30)


if __name__ == '__main__':
    unittest.main()
