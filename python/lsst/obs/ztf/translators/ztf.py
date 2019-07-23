# This file is currently part of obs_lsst but is written to allow it
# to be migrated to the astro_metadata_translator package at a later date.
#
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the LICENSE file in this directory for details of code ownership.
#
# Use of this source code is governed by a 3-clause BSD-style
# license that can be found in the LICENSE file.

"""Metadata translation code for LSST auxTel collimation camera headers"""

__all__ = ("ZtfTranslator", )

import logging

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation

from astro_metadata_translator.translators.helpers import is_non_science, tracking_from_degree_headers

from astro_metadata_translator import cache_translation, FitsTranslator


log = logging.getLogger(__name__)

# Date instrument is taking data at telescope
# Prior to this date many parameters are automatically nulled out
# since the headers have not historically been reliable
TSTART = Time("2019-06-01T00:00", format="isot", scale="utc")

# Define the sensor and group name for AuxTel globally so that it can be used
# in multiple places. There is no raft but for consistency with other LSST
# cameras we define one.
_DETECTOR_GROUP_NAME = "RXX"
_DETECTOR_NAME = "S00"


def is_non_science_or_lab(self):
    """Pseudo method to determine whether this is a lab or non-science
    header.

    Raises
    ------
    KeyError
        If this is a science observation and on the mountain.
    """
    if is_non_science(self):
        return
    if not self._is_on_mountain():
        return
    raise KeyError("Required key is missing and this is a mountain science observation")


class ZtfTranslator(FitsTranslator):
    """Metadata translator for the ZTF Mosaic Camera data.
    """

    name = "ZTF/Mosaic"
    """Name of this translation class"""

    _const_map = {
        "instrument": "ZTF/Mosaic",     # Supports the Mosaic Camera on ZTF
        "telescope": "ZTF",
        "detector_group": _DETECTOR_GROUP_NAME,
        "detector_num": 0,
        "detector_name": _DETECTOR_NAME,  # Single sensor
        "detector_serial": "ASI1600MM_Pro",
        "boresight_rotation_coord": "unknown",
        "science_program": "unknown",
        "relative_humidity": None,
        "pressure": None,
        "temperature": None,
        "altaz_begin": None,
    }

    _trivial_map = {
        "observation_id": ("OBSID", dict(default=None, checker=is_non_science)),
        "boresight_airmass": ("AIRMASS", dict(checker=is_non_science_or_lab)),
        "physical_filter": "FILTER",
        "object": ("OBJECT", dict(checker=is_non_science_or_lab, default="UNKNOWN")),
        "boresight_rotation_angle": ("TELRAD", dict(checker=is_non_science_or_lab,
                                                    default=float("nan"), unit=u.deg)),
    }

    DETECTOR_GROUP_NAME = _DETECTOR_GROUP_NAME
    """Fixed name of detector group."""

    DETECTOR_NAME = _DETECTOR_NAME
    """Fixed name of single sensor."""

    @classmethod
    def can_translate(cls, header, filename=None):
        """Indicate whether this translation class can translate the
        supplied header.

        Parameters
        ----------
        header : `dict`-like
            Header to convert to standardized form.
        filename : `str`, optional
            Name of file being translated.

        Returns
        -------
        can : `bool`
            `True` if the header is recognized by this class. `False`
            otherwise.
        """
        return True if ("INSTRUME" in header and header["INSTRUME"] == "ZTF/MOSAIC") else False

        # Calibration files strip important headers at the moment so guess
        if "DETNAME" in header and header["DETNAME"] == "RXX_S00":
            return True
        return False

    def _is_on_mountain(self):
        date = self.to_datetime_begin()
        if date is None or date > TSTART:
            return True
        return False

    @cache_translation
    def to_detector_exposure_id(self):
        # Docstring will be inherited. Property defined in properties.py
        return self.to_exposure_id()

    @cache_translation
    def to_location(self):
        # Docstring will be inherited. Property defined in properties.py
        return EarthLocation.from_geodetic(
            self._header["OBSLON"], self._header["OBSLAT"], self._header["OBSALT"])

    @cache_translation
    def to_dark_time(self):
        return self.to_exposure_time()

    @staticmethod
    def compute_exposure_id(dayobs, seqnum):
        """Helper method to calculate the AuxTel exposure_id.

        Parameters
        ----------
        dayobs : `str`
            Day of observation in either YYYYMMDD or YYYY-MM-DD format.
        seqnum : `int` or `str`
            Sequence number.

        Returns
        -------
        exposure_id : `int`
            Exposure ID in form YYYYMMDDnnnnn form.
        """
        dayobs = dayobs.replace("-", "")

        if len(dayobs) != 8:
            raise ValueError(f"Malformed dayobs: {dayobs}")

        # Expect no more than 99,999 exposures in a day
        maxdigits = 5
        if seqnum >= 10**maxdigits:
            raise ValueError(f"Sequence number ({seqnum}) exceeds limit")

        # Form the number as a string zero padding the sequence number
        idstr = f"{dayobs}{seqnum:0{maxdigits}d}"
        return int(idstr)

    @cache_translation
    def to_exposure_id(self):
        """Return a unique exposure ID number

        Returns
        -------
        exposure_id : `int`
            Unique exposure number.
        """
        if "CALIB_ID" in self._header:
            self._used_these_cards("CALIB_ID")
            return None

        return self._header["FRAMENUM"]

    # For now "visits" are defined to be identical to exposures.
    to_visit_id = to_exposure_id

    @cache_translation
    def to_exposure_time(self):
        # Docstring will be inherited. Property defined in properties.py
        # Some data is missing a value for EXPTIME.
        # Have to be careful we do not have circular logic when trying to
        # guess
        if self.is_key_ok("EXPTIME"):
            return self.quantity_from_card("EXPTIME", u.s)

        # A missing or undefined EXPTIME is problematic. Set to -1
        # to indicate that none was found.
        log.warning("Insufficient information to derive exposure time. Setting to -1.0s")
        return -1.0 * u.s

    @cache_translation
    def to_tracking_radec(self):
        # Docstring will be inherited. Property defined in properties.py
        radecsys = ("RADESYS",)
        radecpairs = (("TELRA", "TELDEC"),)
        return tracking_from_degree_headers(self, radecsys, radecpairs, unit=(u.hourangle, u.deg))

    @cache_translation
    def to_observation_type(self):
        """Determine the observation type.

        In the absence of an ``IMGTYPE`` header, assumes lab data is always a
        dark if exposure time is non-zero, else bias.

        Returns
        -------
        obstype : `str`
            Observation type.
        """

        # AuxTel observation type is documented to appear in OBSTYPE
        # but for historical reasons prefers IMGTYPE.  Some data puts
        # it in GROUPID (which is meant to be for something else).
        # Test the keys in order until we find one that contains a
        # defined value.
        obstype_keys = ["OBSTYPE", "IMGTYPE"]

        # For now, hope that GROUPID does not contain an obs type value
        # when on the mountain.
        if not self._is_on_mountain():
            obstype_keys.append("GROUPID")

        for k in obstype_keys:
            if self.is_key_ok(k):
                obstype = self._header[k]
                self._used_these_cards(k)
                return obstype.lower()

        # In the absence of any observation type information, return
        # unknown unless we think it might be a bias.
        exptime = self.to_exposure_time()
        if exptime == 0.0:
            obstype = "bias"
        else:
            obstype = "unknown"
        log.warning("Unable to determine observation type. Guessing '%s'", obstype)
        return obstype
