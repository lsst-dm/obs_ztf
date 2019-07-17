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
import lsst.daf.persistence as dafPersistence
import lsst.afw.image.utils as afwImageUtils
import lsst.obs.base.yamlCamera as yamlCamera
from lsst.obs.base import CameraMapper, MakeRawVisitInfoViaObsInfo

from .translators.ztf import ZtfTranslator

__all__ = ["ZtfMapper", "ZtfMakeRawVisitInfo"]


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

        afwImageUtils.defineFilter('g', 0.0, alias=[])
        afwImageUtils.defineFilter('R', 0.0, alias=['r'])
        afwImageUtils.defineFilter('i', 0.0, alias=[])

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
