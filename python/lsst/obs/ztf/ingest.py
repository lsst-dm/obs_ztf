from __future__ import print_function
import os
import re
import lsst.daf.base
from lsst.pipe.tasks.ingest import ParseTask
from lsst.pipe.tasks.ingestCalibs import CalibsParseTask

EXTENSIONS = ["fits", "gz", "fz"]  # Filename extensions to strip off

class ZtfParseTask(ParseTask):
    """Parser suitable for lab data"""

    def _getInfo(self, filename):
        # Grab the basename
        phuInfo, infoList = ParseTask.getInfo(self, filename)
        basename = os.path.basename(filename)
        while any(basename.endswith("." + ext) for ext in EXTENSIONS):
            basename = basename[:basename.rfind('.')]
        phuInfo['basename'] = basename
        return phuInfo, infoList

    def translate_visit(self, md):
        """Generate a unique visit from the timestamp

        Parameters
        ----------
        md : `lsst.daf.base.PropertyList or PropertySet`
            image metadata

        Returns
        -------
        visit_num : `int`
            Visit number, as translated
        """
        mjd = md.get("OBSMJD")
        mmjd = mjd - 55197              # relative to 2010-01-01, just to make the visits a tiny bit smaller
        return int(1e5*mmjd)            # 86400s per day, so we need this resolution

    def translate_imageType(self, md):
        """Return the type of the image (e.g. dark)"""
        if False:
            imageType = None            # maybe this really comes from MODE_NUM??
        else:
            fileName = md.get("ORIGNAME").strip()
            mat = re.search(r"_(.)\.fits$", fileName)
            imageType = mat.group(1) if mat else None

            # from ztf_pipelines_deliverables.pdf as provided by Eric Bellm, Nov. 2017
            imageType = dict(o = 'onsky',
                             b = 'bias',
                             d = 'dark',
                             f = 'flat',
                             c = 'focus',
                             g = 'guider').get(imageType)

        return imageType

    def translate_dateObs(self, md):
        """Generate the day of the observation (UT) in the form "2017-11-06"

        Parameters
        ----------
        md : `lsst.daf.base.PropertyList or PropertySet`
            image metadata

        Returns
        -------
            The day of the observation
        """
        mjd = md.get("DATE-OBS")

        return mjd.split("T")[0]

##############################################################################################################

class ZtfCalibsParseTask(CalibsParseTask):
    """Parser for calibs"""

    def _translateFromCalibId(self, field, md):
        """Get a value from the CALIB_ID written by constructCalibs"""
        data = md.get("CALIB_ID")
        match = re.search(".*%s=(\S+)" % field, data)
        return match.groups()[0]

    def _translate_ccd(self, md):
        return self._translateFromCalibId("ccd", md)

    def _translate_filter(self, md):
        return self._translateFromCalibId("filter", md)

    def _translate_calibDate(self, md):
        return self._translateFromCalibId("calibDate", md)
