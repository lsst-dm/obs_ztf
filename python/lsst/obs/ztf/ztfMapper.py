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
import lsst.afw.image.utils as afwImageUtils
import lsst.afw.geom as afwGeom
import lsst.afw.cameraGeom as cameraGeom
import lsst.obs.base.yamlCamera as yamlCamera
from astro_metadata_translator import fix_header
from lsst.obs.base import CameraMapper, MakeRawVisitInfoViaObsInfo, bboxFromIraf

from .translators.ztf import ZtfTranslator

__all__ = ["ZtfMapper", "ZtfMakeRawVisitInfo", "assemble_raw"]


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


def assemble_raw(dataId, componentInfo, cls):
    """Called by the butler to construct the composite type "raw".

    Note that we still need to define "_raw" and copy various fields over.

    Parameters
    ----------
    dataId : `lsst.daf.persistence.dataId.DataId`
        The data ID.
    componentInfo : `dict`
        dict containing the components, as defined by the composite definition
        in the mapper policy.
    cls : 'object'
        Unused.

    Returns
    -------
    exposure : `lsst.afw.image.Exposure`
        The assembled exposure.
    """
    from lsst.ip.isr import AssembleCcdTask

    config = AssembleCcdTask.ConfigClass()
    config.doTrim = False

    assembleTask = AssembleCcdTask(config=config)

    ampExps = componentInfo['raw_amp'].obj
    if len(ampExps) == 0:
        raise RuntimeError("Unable to read raw_amps for %s" % dataId)

    ccd = ampExps[0].getDetector()      # the same (full, CCD-level) Detector is attached to all ampExps
    #
    # Check that the geometry in the metadata matches cameraGeom
    #
    logger = lsst.log.Log.getLogger("ZtfMapper")
    warned = False
    for i, (amp, ampExp) in enumerate(zip(ccd, ampExps)):
        ampMd = ampExp.getMetadata().toDict()

        if amp.getRawBBox() != ampExp.getBBox():  # Oh dear. cameraGeom is wrong -- probably overscan
            if amp.getRawDataBBox().getDimensions() != amp.getBBox().getDimensions():
                raise RuntimeError("Active area is the wrong size: %s v. %s" %
                                   (amp.getRawDataBBox().getDimensions(), amp.getBBox().getDimensions()))
            if not warned:
                logger.warn("amp.getRawBBox() != data.getBBox(); patching. (%s v. %s)",
                            amp.getRawBBox(), ampExp.getBBox())
                warned = True

            w, h = ampExp.getBBox().getDimensions()
            ow, oh = amp.getRawBBox().getDimensions()  # "old" (cameraGeom) dimensions
            #
            # We could trust the BIASSEC keyword, or we can just assume that
            # they've changed the number of overscan pixels (serial and/or
            # parallel).  As Jim Chiang points out, the latter is safer
            #
            bbox = amp.getRawHorizontalOverscanBBox()
            hOverscanBBox = geom.BoxI(bbox.getBegin(),
                                      geom.ExtentI(w - bbox.getBeginX(), bbox.getHeight()))
            bbox = amp.getRawVerticalOverscanBBox()
            vOverscanBBox = geom.BoxI(bbox.getBegin(),
                                      geom.ExtentI(bbox.getWidth(), h - bbox.getBeginY()))

            amp.setRawBBox(ampExp.getBBox())
            amp.setRawHorizontalOverscanBBox(hOverscanBBox)
            amp.setRawVerticalOverscanBBox(vOverscanBBox)
            #
            # This gets all the geometry right for the amplifier, but the size
            # of the untrimmed image will be wrong and we'll put the amp
            # sections in the wrong places, i.e.
            #   amp.getRawXYOffset()
            # will be wrong.  So we need to recalculate the offsets.
            #
            xRawExtent, yRawExtent = amp.getRawBBox().getDimensions()

            x0, y0 = amp.getRawXYOffset()
            ix, iy = x0//ow, y0/oh
            x0, y0 = ix*xRawExtent, iy*yRawExtent
            amp.setRawXYOffset(geom.ExtentI(ix*xRawExtent, iy*yRawExtent))
        #
        # Check the "IRAF" keywords, but don't abort if they're wrong
        #
        if False:
            #
            # Only warn about the first amp, use debug for the others
            #
            detsec = bboxFromIraf(ampMd["DETSEC"]) if "DETSEC" in ampMd else None
            datasec = bboxFromIraf(ampMd["DATASEC"]) if "DATASEC" in ampMd else None
            biassec = bboxFromIraf(ampMd["BIASSEC"]) if "BIASSEC" in ampMd else None

            logCmd = logger.warn if i == 0 else logger.debug
            if detsec and amp.getBBox() != detsec:
                logCmd("DETSEC doesn't match for %s (%s != %s)", dataId, amp.getBBox(), detsec)
            if datasec and amp.getRawDataBBox() != datasec:
                logCmd("DATASEC doesn't match for %s (%s != %s)", dataId, amp.getRawDataBBox(), detsec)
            if biassec and amp.getRawHorizontalOverscanBBox() != biassec:
                logCmd("BIASSEC doesn't match for %s (%s != %s)",
                       dataId, amp.getRawHorizontalOverscanBBox(), detsec)

    ampDict = {}
    for amp, ampExp in zip(ccd, ampExps):
        ampDict[amp.getName()] = ampExp

    exposure = assembleTask.assembleCcd(ampDict)

    md = componentInfo['raw_hdu'].obj
    fix_header(md)  # No mapper so cannot specify the translator class
    exposure.setMetadata(md)
    visitInfo = ZtfMakeRawVisitInfo(logger)(md)
    exposure.getInfo().setVisitInfo(visitInfo)

    boresight = visitInfo.getBoresightRaDec()
    rotangle = visitInfo.getBoresightRotAngle()

    if boresight.isFinite():
        exposure.setWcs(getWcsFromDetector(exposure.getDetector(), boresight,
                                           90*geom.degrees - rotangle))
    else:
        # Should only warn for science observations but VisitInfo does not know
        logger.warn("Unable to set WCS for %s from header as RA/Dec/Angle are unavailable" %
                    (dataId,))

    return exposure


class ZtfMakeRawVisitInfo(MakeRawVisitInfoViaObsInfo):
    """Make a VisitInfo from the FITS header of a raw image."""


class ZtfMapper(CameraMapper):
    packageName = 'obs_ztf'

    MakeRawVisitInfoClass = ZtfMakeRawVisitInfo
    translatorClass = ZtfTranslator

    def __init__(self, **kwargs):
        policyFile = dafPersistence.Policy.defaultPolicyFile(self.packageName, "ztfMapper.yaml", "policy")
        policy = dafPersistence.Policy(policyFile)

        super(ZtfMapper, self).__init__(policy, os.path.dirname(policyFile), **kwargs)
        #
        # The composite objects don't seem to set these
        #
        for d in (self.mappings, self.exposures):
            d['raw'] = d['_raw']
        
        afwImageUtils.resetFilters()
        afwImageUtils.defineFilter('ZTF_g', 0.0, alias=['g'])
        afwImageUtils.defineFilter('ZTF_r', 0.0, alias=['r'])
        afwImageUtils.defineFilter('ZTF_i', 0.0, alias=['i'])

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

    def query_raw_amp(self, format, dataId):
        """Return a list of tuples of values of the fields specified in
        format, in order.

        Parameters
        ----------
        format : `list`
            The desired set of keys.
        dataId : `dict`
            A possible-incomplete ``dataId``.

        Returns
        -------
        fields : `list` of `tuple`
            Values of the fields specified in ``format``.

        Raises
        ------
        ValueError
            The channel number requested in ``dataId`` is out of range.
        """
        nChannel = 4                    # number of possible channels, 1..nChannel

        if "channel" in dataId:         # they specified a channel
            dataId = dataId.copy()
            channel = dataId.pop('channel')  # Do not include in query below
            if channel > nChannel or channel < 1:
                raise ValueError(f"Requested channel is out of range 0 < {channel} <= {nChannel}")
            channels = [channel]
        else:
            channels = range(1, nChannel + 1)  # we want all possible channels

        if "channel" in format:           # they asked for a channel, but we mustn't query for it
            format = list(format)
            channelIndex = format.index('channel')  # where channel values should go
            format.pop(channelIndex)
        else:
            channelIndex = None

        dids = []                       # returned list of dataIds
        for value in self.query_raw(format, dataId):
            if channelIndex is None:
                dids.append(value)
            else:
                for c in channels:
                    did = list(value)
                    did.insert(channelIndex, c)
                    dids.append(tuple(did))

        return dids
    #
    # The composite type "raw" doesn't provide e.g. query_raw, so we defined
    # type _raw in the .paf file with the same template, and forward requests
    # as necessary
    #

    def query_raw(self, *args, **kwargs):
        """Magic method that is called automatically if it exists.

        This code redirects the call to the right place, necessary because of
        leading underscore on ``_raw``.
        """
        return self.query__raw(*args, **kwargs)

    def map_raw_md(self, *args, **kwargs):
        """Magic method that is called automatically if it exists.

        This code redirects the call to the right place, necessary because of
        leading underscore on ``_raw``.
        """
        return self.map__raw_md(*args, **kwargs)

    def map_raw_filename(self, *args, **kwargs):
        """Magic method that is called automatically if it exists.

        This code redirects the call to the right place, necessary because of
        leading underscore on ``_raw``.
        """
        return self.map__raw_filename(*args, **kwargs)

    def bypass_raw_filename(self, *args, **kwargs):
        """Magic method that is called automatically if it exists.

        This code redirects the call to the right place, necessary because of
        leading underscore on ``_raw``.
        """
        return self.bypass__raw_filename(*args, **kwargs)

    def map_raw_visitInfo(self, *args, **kwargs):
        """Magic method that is called automatically if it exists.

        This code redirects the call to the right place, necessary because of
        leading underscore on ``_raw``.
        """
        return self.map__raw_visitInfo(*args, **kwargs)

    def bypass_raw_visitInfo(self, datasetType, pythonType, location, dataId):
        fileName = location.getLocationsWithRoot()[0]
        mat = re.search(r"\[(\d+)\]$", fileName)
        if mat:
            hdu = int(mat.group(1))
            md = readMetadata(fileName, hdu=hdu)
        else:
            md = readMetadata(fileName)  # or hdu = INT_MIN; -(1 << 31)

        makeVisitInfo = self.MakeRawVisitInfoClass(log=self.log)
        fix_header(md, translator_class=self.translatorClass)
        return makeVisitInfo(md)

    def std_raw_amp(self, item, dataId):
        return self._standardizeExposure(self.exposures['raw_amp'], item, dataId,
                                         trimmed=False, setVisitInfo=False)

    def std_raw(self, item, dataId, filter=True):
        """Standardize a raw dataset by converting it to an
        `~lsst.afw.image.Exposure` instead of an `~lsst.afw.image.Image`."""

        return self._standardizeExposure(self.exposures['raw'], item, dataId, trimmed=False,
                                         setVisitInfo=False,  # it's already set, and the metadata's stripped
                                         filter=filter)
