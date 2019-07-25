import os.path

from lsst.utils import getPackageDir

ObsConfigDir = os.path.join(getPackageDir("obs_ztf"), "config")

for sub in ("isr", "charImage", "calibrate"):
    path = os.path.join(ObsConfigDir, sub + ".py")
    if os.path.exists(path):
        getattr(config, sub).load(path)
