from lsst.obs.ztf.ingest import ZtfParseTask
config.parse.retarget(ZtfParseTask)
config.parse.translation = {
    # translate_visit
    'field': 'FIELDID',
    # translate_filter
    # translate_dateObs
    ### 'date': 'DATE-OBS',
    'expTime': 'AEXPTIME',
    'ccd': 'CCD_ID',
    # translate_imageType
    'filename': 'ORIGNAME',
}
config.parse.translators = {
    'visit': 'translate_visit',
    'filter': 'translate_filter',       # defined in base class
    'dateObs': 'translate_dateObs',
    'imageType': 'translate_imageType',
}
config.parse.defaults = {
}

config.register.columns = {
    'visit': 'int',
    'field': 'int',
    'filter': 'text',
    'dateObs': 'text',
    'expTime': 'double',
    'ccd': 'int',
    'imageType': 'text',
    'filename': 'text',
}
config.register.visit = list(config.register.columns.keys())
