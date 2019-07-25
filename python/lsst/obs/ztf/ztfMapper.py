#
# LSST Data Management System
# Copyright 2016 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import os
import lsst.utils
import lsst.geom as geom
import lsst.log
import lsst.daf.persistence as dafPersistence
import lsst.afw.image as afwImage
import lsst.afw.image.utils as afwImageUtils
import lsst.afw.geom as afwGeom
import lsst.afw.cameraGeom as cameraGeom
import lsst.obs.base.yamlCamera as yamlCamera
from astro_metadata_translator import fix_header
from lsst.obs.base import CameraMapper, MakeRawVisitInfoViaObsInfo

from .translators.ztf import ZtfTranslator

__all__ = ["ZtfMapper", "ZtfMakeRawVisitInfo"]


def getWcsFromDetector(detector, boresight, rotation=0*geom.degrees, flipX=False):
    """Given a detector and (boresight, rotation), return that detector's WCS

    Parameters
    ----------
    camera : `lsst.afw.cameraGeom.Camera`
        The camera containing the detector.
    detector : `lsst.afw.cameraGeom.Detector`
        A detector in a camera.
    boresight : `lsst.afw.geom.SpherePoint`
       The boresight of the observation.
    rotation : `lsst.afw.geom.Angle`, optional
        The rotation angle of the camera.
        The rotation is "rotskypos", the angle of sky relative to camera
        coordinates (from North over East).
    flipX : `bool`, optional
        Flip the X axis?

    Returns
    -------
    wcs : `lsst::afw::geom::SkyWcs`
        The calculated WCS.
    """
    trans = detector.getTransform(detector.makeCameraSys(cameraGeom.PIXELS),
                                  detector.makeCameraSys(cameraGeom.FIELD_ANGLE))

    wcs = afwGeom.makeSkyWcs(trans, rotation, flipX, boresight)

    return wcs


class ZtfMakeRawVisitInfo(MakeRawVisitInfoViaObsInfo):
    """Make a VisitInfo from the FITS header of a raw image."""


class ZtfMapper(CameraMapper):
    packageName = 'obs_ztf'

    MakeRawVisitInfoClass = ZtfMakeRawVisitInfo
    translatorClass = ZtfTranslator

    def __init__(self, **kwargs):
        policyFile = dafPersistence.Policy.defaultPolicyFile(self.packageName, "ztfMapper.yaml", "policy")
        policy = dafPersistence.Policy(policyFile)
        #
        # Look for the calibrations root "root/CALIB" if not supplied
        #
        if kwargs.get('root', None) and not kwargs.get('calibRoot', None):
            calibSearch = [os.path.join(kwargs['root'], 'CALIB')]
            if "repositoryCfg" in kwargs:
                calibSearch += [os.path.join(cfg.root, 'CALIB') for cfg in kwargs["repositoryCfg"].parents if
                                hasattr(cfg, "root")]
                calibSearch += [cfg.root for cfg in kwargs["repositoryCfg"].parents if hasattr(cfg, "root")]
            for calibRoot in calibSearch:
                if os.path.exists(os.path.join(calibRoot, "calibRegistry.sqlite3")):
                    kwargs['calibRoot'] = calibRoot
                    break
            if not kwargs.get('calibRoot', None):
                lsst.log.Log.getLogger("ZtfCamMapper").warn("Unable to find valid calib root directory")

        super(ZtfMapper, self).__init__(policy, os.path.dirname(policyFile), **kwargs)

        afwImageUtils.resetFilters()
        afwImageUtils.defineFilter('NONE', 0.0)
        afwImageUtils.defineFilter('ZTF_g', 0.0, alias=['g'])
        afwImageUtils.defineFilter('ZTF_r', 0.0, alias=['r'])
        afwImageUtils.defineFilter('ZTF_i', 0.0, alias=['i'])

        ZtfMapper._nbit_tract = 16
        ZtfMapper._nbit_patch = 5
        ZtfMapper._nbit_filter = 6

        ZtfMapper._nbit_id = 64 - (ZtfMapper._nbit_tract + 2*ZtfMapper._nbit_patch + ZtfMapper._nbit_filter)

    @classmethod
    def _makeCamera(cls, policy=None, repositoryDir=None, cameraYamlFile=None):
        """Make a camera  describing the camera geometry.

        policy : ignored
        repositoryDir : ignored
        cameraYamlFile : `str`
           The full path to a yaml file to be passed to `yamlCamera.makeCamera`

        Returns
        -------
        camera : `lsst.afw.cameraGeom.Camera`
            Camera geometry.
        """

        if not cameraYamlFile:
            cameraYamlFile = os.path.join(lsst.utils.getPackageDir(cls.packageName), "policy", "ztf.yaml")

        return yamlCamera.makeCamera(cameraYamlFile)

    def _extractDetectorName(self, dataId):
        return "%02d" % (dataId["ccd"])

    def _computeCcdExposureId(self, dataId):
        """Compute the 64-bit (long) identifier for a CCD exposure.

        @param dataId (dict) Data identifier with visit
        """
        visit = dataId['visit']
        return int(visit)

    def bypass_raw(self, datasetType, pythonType, location, dataId):
        """Magic method that is called automatically if it exists.
        """
        from lsst.ip.isr import AssembleCcdTask

        config = AssembleCcdTask.ConfigClass()
        config.doTrim = False

        assembleTask = AssembleCcdTask(config=config)
        logger = lsst.log.Log.getLogger("ZtfMapper")

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        assert len(location.getLocations()) == 1  # KTL says that this is always true, but check
        relativeFileName = location.getLocations()[0]
        if location.storage.exists(relativeFileName):
            fileName = location.storage.locationWithRoot(relativeFileName)
            data = afwImage.ImageF(fileName, hdu=1)
            data = pythonType(afwImage.MaskedImageF(data))
        else:
            raise IOError(f"Unable to find {relativeFileName}")

        # Get the Detector
        detector = self.camera[self._extractDetectorName(dataId)]
        #
        # Read all the data
        #
        ampDict = {}
        for i, amp in enumerate(detector, 1):
            ampExp = pythonType(amp.getRawBBox())
            ampExp.setDetector(detector)
            ampDict[amp.getName()] = ampExp

            hdu = amp.get('hdu')
            if i == 1:
                assert i == hdu         # don't read twice
                data = data.image
            else:
                data = afwImage.ImageF(fileName, hdu=hdu)

            ampExp[amp.getRawDataBBox()].image = data

            bias = afwImage.ImageF(fileName, hdu=hdu + 4)
            ampExp[amp.getRawHorizontalOverscanBBox()].image = bias

        exposure = assembleTask.assembleCcd(ampDict)

        md = afwImage.readMetadata(fileName, hdu=0)
        fix_header(md, translator_class=self.translatorClass)
        exposure.setMetadata(md)

        visitInfo = ZtfMakeRawVisitInfo(logger)(md)
        exposure.getInfo().setVisitInfo(visitInfo)

        boresight = visitInfo.getBoresightRaDec()
        rotangle = visitInfo.getBoresightRotAngle()

        if boresight.isFinite():
            exposure.setWcs(getWcsFromDetector(exposure.getDetector(), boresight,
                                               90*geom.degrees - rotangle))
        else:
            logger.warn(f"Unable to set WCS for {dataId} from header as RA/Dec/Angle are unavailable")

        return exposure

    def _computeCcdExposureId(self, dataId):
        """Compute the 64-bit (long) identifier for a CCD exposure.

        Parameters
        ----------
        dataId : `dict`
            Data identifier including visit

        Returns
        -------
        id : `int`
            Integer identifier for a CCD exposure.
        """
        try:
            visit = dataId["visit"]
            detector = dataId["ccd"]
        except Exception:
            raise KeyError(f"Require a visit and ccd to calculate detector exposure ID. Got: {dataId}")

        return 100*visit + detector

    def bypass_ccdExposureId(self, datasetType, pythonType, location, dataId):
        return self._computeCcdExposureId(dataId)

    def bypass_ccdExposureId_bits(self, datasetType, pythonType, location, dataId):
        """How many bits are required for the maximum exposure ID"""
        return 31  # max visit ~ 21400000 (a 30s cadence every 10-hour night for 20 years)
    
    def _computeCoaddExposureId(self, dataId, singleFilter):
        """Compute the 64-bit (long) identifier for a coadd.

        Parameters
        ----------
        dataId : `dict`
            Data identifier with tract and patch.
        singleFilter : `bool`
            True means the desired ID is for a single-filter coadd, in which
            case ``dataId`` must contain filter.
        """

        tract = int(dataId['tract'])
        if tract < 0 or tract >= 2**ZtfMapper._nbit_tract:
            raise RuntimeError('tract not in range [0,%d)' % (2**ZtfMapper._nbit_tract))
        patchX, patchY = [int(patch) for patch in dataId['patch'].split(',')]
        for p in (patchX, patchY):
            if p < 0 or p >= 2**ZtfMapper._nbit_patch:
                raise RuntimeError('patch component not in range [0, %d)' % 2**ZtfMapper._nbit_patch)
        oid = (((tract << ZtfMapper._nbit_patch) + patchX) << ZtfMapper._nbit_patch) + patchY
        if singleFilter:
            return (oid << ZtfMapper._nbit_filter) + afwImage.Filter(dataId['filter']).getId()
        return oid

    def bypass_deepCoaddId_bits(self, *args, **kwargs):
        """The number of bits used up for patch ID bits."""
        return 64 - ZtfMapper._nbit_id

    def bypass_deepCoaddId(self, datasetType, pythonType, location, dataId):
        return self._computeCoaddExposureId(dataId, True)

    def bypass_dcrCoaddId_bits(self, datasetType, pythonType, location, dataId):
        return self.bypass_deepCoaddId_bits(datasetType, pythonType, location, dataId)

    def bypass_dcrCoaddId(self, datasetType, pythonType, location, dataId):
        return self.bypass_deepCoaddId(datasetType, pythonType, location, dataId)

    def bypass_deepMergedCoaddId_bits(self, *args, **kwargs):
        """The number of bits used up for patch ID bits."""
        return 64 - ZtfMapper._nbit_id

    def bypass_deepMergedCoaddId(self, datasetType, pythonType, location, dataId):
        return self._computeCoaddExposureId(dataId, False)

    def bypass_dcrMergedCoaddId_bits(self, *args, **kwargs):
        """The number of bits used up for patch ID bits."""
        return self.bypass_deepMergedCoaddId_bits(*args, **kwargs)

    def bypass_dcrMergedCoaddId(self, datasetType, pythonType, location, dataId):
        return self.bypass_deepMergedCoaddId(datasetType, pythonType, location, dataId)
