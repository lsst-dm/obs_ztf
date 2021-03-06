from lsst.obs.ztf.ingest import ZtfParseTask
config.parse.retarget(ZtfParseTask)
config.parse.hdu = 0

config.parse.translation = {
    # translate_visit
    'field': 'FRAMENUM',
    # translate_filter
    # translate_dayObs
    ### 'date': 'DATE-OBS',
    'expTime': 'EXPTIME',
    'ccd': 'CCD_ID',
    # translate_imageType
    'filename': 'ORIGNAME',
    'taiObs': 'DATE-OBS',
}
config.parse.translators = {
    'visit': 'translate_visit',
    'filter': 'translate_filter',       # defined in base class
    'dateObs': 'translate_dayObs',
    'imageType': 'translate_imageType',
}
config.parse.defaults = {
}

config.register.columns = {
    'visit': 'int',
    'field': 'int',
    'filter': 'text',
    'dateObs': 'text',
    'taiObs': 'text',
    'expTime': 'double',
    'ccd': 'int',
    'imageType': 'text',
    'filename': 'text',
}
config.register.visit = list(config.register.columns.keys())
