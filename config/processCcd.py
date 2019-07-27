import os.path

from lsst.utils import getPackageDir

ObsConfigDir = os.path.join(getPackageDir("obs_ztf"), "config")

for sub in ("isr", "charImage", "calibrate"):
    path = os.path.join(ObsConfigDir, sub + ".py")
    if os.path.exists(path):
        getattr(config, sub).load(path)


for loader in [config.calibrate.astromRefObjLoader, config.calibrate.photoRefObjLoader]:
    loader.ref_dataset_name = 'gaia_DR2'
    loader.filterMap = dict(ZTF_g = 'phot_g_mean_mag',
                            ZTF_r = 'phot_g_mean_mag',
                            ZTF_i = 'phot_g_mean_mag',
    )


config.charImage.repair.cosmicray.nCrPixelMax = 1000000
config.charImage.repair.interp.modelPsf.defaultFwhm = 0.1

config.isr.fwhm=1.00
config.charImage.installSimplePsf.fwhm=1.01
config.charImage.repair.interp.modelPsf.defaultFwhm=1.02

# PSF determination
config.charImage.measurePsf.reserve.fraction = 0.2
config.charImage.measurePsf.starSelector["objectSize"].widthMin = 0.9
config.charImage.measurePsf.starSelector["objectSize"].fluxMin = 4000
try:
    import lsst.meas.extensions.psfex.psfexPsfDeterminer
    config.charImage.measurePsf.psfDeterminer["psfex"].spatialOrder = 2
    config.charImage.measurePsf.psfDeterminer["psfex"].psfexBasis = 'PIXEL_AUTO'
    config.charImage.measurePsf.psfDeterminer["psfex"].samplingSize = 0.5
    config.charImage.measurePsf.psfDeterminer["psfex"].kernelSize = 81
    config.charImage.measurePsf.psfDeterminer.name = "psfex"
except ImportError as e:
    print("WARNING: Unable to use psfex: %s" % e)
    config.charImage.measurePsf.psfDeterminer.name = "pca"

config.charImage.detection.thresholdValue=70

config.calibrate.detection.thresholdValue=80
config.calibrate.astrometry.referenceSelector.doSignalToNoise=False # all errors are NaN
config.calibrate.astrometry.referenceSelector.signalToNoise.minimum=10
config.calibrate.astrometry.referenceSelector.signalToNoise.fluxField='phot_g_mean_mag_flux'
config.calibrate.astrometry.referenceSelector.signalToNoise.errField='phot_g_mean_mag_fluxErr'

config.calibrate.astrometry.referenceSelector.doMagLimit=True
config.calibrate.astrometry.referenceSelector.magLimit.fluxField='phot_g_mean_mag_flux'
config.calibrate.astrometry.referenceSelector.magLimit.maximum = 16.5 # AB

config.calibrate.astrometry.wcsFitter.order = 5
