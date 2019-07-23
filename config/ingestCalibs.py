from lsst.obs.ztf.ingest import ZtfCalibsParseTask
config.parse.retarget(ZtfCalibsParseTask)

config.register.columns = {'filter': 'text',
                           'ccd': 'int',
                           'calibDate': 'text',
                           'validStart': 'text',
                           'validEnd': 'text',
                           }

config.parse.translation = {
    #'field': 'FRAMENUM',
    'expTime': 'EXPTIME',
    'ccd': 'CCD_ID',
    #'filename': 'ORIGNAME',
    'taiObs': 'DATE-OBS',
}

config.parse.translators = {
    #'ccd': 'translate_ccd',
    'filter': 'translate_filter',
    'calibDate': 'translate_calibDate',
}

config.register.unique = ['filter', 'ccd', 'calibDate']
config.register.tables = ['bias', 'dark', 'flat', 'fringe']
config.register.visit = ['calibDate', 'filter']
