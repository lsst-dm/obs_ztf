from lsst.obs.monocam.ingest import MonocamCalibsParseTask
config.parse.retarget(MonocamCalibsParseTask)

config.register.columns = {'filter': 'text',
                           'ccd': 'int',
                           'calibDate': 'text',
                           'validStart': 'text',
                           'validEnd': 'text',
                           }

config.parse.translators = {'ccd': 'translate_ccd',
                            'filter': 'translate_filter',
                            'calibDate': 'translate_calibDate',
                            }

config.register.unique = ['filter', 'ccd', 'calibDate']
config.register.tables = ['bias', 'dark', 'flat', 'fringe']
config.register.visit = ['calibDate', 'filter']
