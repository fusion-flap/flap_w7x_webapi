import unittest
import datetime
import numpy as np
import flap
import flap.modules.w7x_webapi as webapi

# -------------------------------------INDEPENDENT FUNCTIONS FROM FLAP--------------------------------------------------


def test_basic():
    # Getting ECRH data by giving the exact time
    ecrh = webapi.GetSignal('ECRH')
    start = datetime.datetime(2018, 8, 22, 10, 46, 53)
    end = datetime.datetime(2018, 8, 22, 10, 47, 8)
    ecrh.time_query_gen(start, end)
    ecrh_data = ecrh.archive_pull()
    return int(np.mean(ecrh_data["values"]))


def test_shotid():
    # Getting NBI data with shotid
    nbi7 = webapi.GetSignal('NBI-7')
    nbi7.shotid_time_query_gen('20180912.012')
    nbi7.downsamp_gen(512)
    nbi7_data = webapi.ScalarData(nbi7.archive_pull())
    return [int(np.mean(nbi7_data.val)*1000000), int(min(nbi7_data.time))]


def test_streamlink():
    # Getting ECRH data via the streaming link, can be used for arbitrary data
    ecrh2 = webapi.GetSignal('ArchiveDB/codac/W7X/CBG_ECRH/TotalPower_DATASTREAM/V1/0/Ptot_ECRH/scaled/')
    ecrh2.shotid_time_query_gen('20180912.012')
    ecrh2.downsamp_gen(512)
    ecrh2_data = webapi.ScalarData(ecrh2.archive_pull())
    return int(np.mean(ecrh2_data.val))


class StandaloneTest(unittest.TestCase):
    def test_basic(self):
        self.assertEqual(test_basic(), 1136)

    def test_shotid(self):
        self.assertEqual(test_shotid(), [-26, 1536744303613702400])

    def test_streamlink(self):
        self.assertEqual(test_streamlink(), 1248)

# -------------------------------------COMPATIBILITY WITH FLAP----------------------------------------------------------


def test_register():
    flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
    return 'W7X_WEBAPI' in flap.list_data_sources()


def test_reading():
    flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
    d = flap.get_data('W7X_WEBAPI', name='ECRH',
                      exp_id='20180912.012',
                      options={'Scale Time': True,
                               'Check Time Equidistant': True},
                      object_name='ECRH_data')
    return [d.coordinates[0].mode.equidistant, d.coordinates[0].start]


def test_downsampling():
    flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
    d = flap.get_data('W7X_WEBAPI', name='NBI-7',
                      exp_id='20181016.037',
                      options={'Scale Time': True,
                               'Downsample': 1024},
                      object_name='ECRH_data')
    return d.coordinates[0].shape


def test_timerange():
    flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
    d = flap.get_data('W7X_WEBAPI', name='ArchiveDB/codac/W7X/CBG_ECRH/TotalPower_DATASTREAM/V1/0/Ptot_ECRH/scaled/',
                      exp_id='20181016.037',
                      options={'Scale Time': True},
                      object_name='ECRH_data',
                      coordinates={'Time [s]': [2, 3]})
    return [d.coordinates[0].values[0], d.coordinates[0].values[-1]]


def test_imagedata():
    flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
    d = flap.get_data('W7X_WEBAPI', name='AUG-2',
                      exp_id='20181016.037',
                      options={'Scale Time': True},
                      object_name='AUG2_data',
                      coordinates={'Time [s]': [2, 3]})
    return np.mean(d.data)


def test_vectordata():
    flap.register_data_source('W7X_WEBAPI',
                              get_data_func=webapi.get_data,
                              add_coord_func=webapi.add_coordinate)
#    d = flap.get_data('W7X_WEBAPI', name='TS-v8-ne_map',

    d = flap.get_data('W7X_WEBAPI', name='ABES-v2-density',
                      exp_id='20181016.037',
                      options={'Scale Time': True},
                      object_name='ECRH_data',
                      coordinates={'Time [s]': [0, 3]})
    return np.mean(d.data)


class FLAPTest(unittest.TestCase):
    def test_register(self):
        self.assertTrue(test_register())

    def test_reading(self):
        self.assertEqual(test_reading(), [True, [-1.0]])

    def test_downsampling(self):
        self.assertEqual(test_downsampling(), [841])

    def test_timerange(self):
        self.assertEqual(test_timerange(),  [2.0, 2.99996])

    def test_imagedata(self):
        self.assertEqual(test_imagedata(),  2356.9879264322917)

    def test_vectordata(self):
        self.assertEqual(test_vectordata(),  2.614503530033091e+19)


if __name__ == '__main__':
    unittest.main()
