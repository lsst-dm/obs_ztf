config.doConvertIntToFloat=False
config.doBias=False
config.doDark=False
config.doLinearize=False
config.doDefect=False

#
# we get
#   3631e6*15.09/R
# photons/s per m^2 for a 0-th mag star. The clear aperture is c. 3.6m^2, and at R = 5 and a QE of 0.9
# this is c. 3.5e10 photons/s
# 
fluxPerSec = 3.5e10
config.fluxMag0T1={
    'ZTF_g': fluxPerSec,
    'ZTF_r': fluxPerSec,
    'ZTF_i': fluxPerSec,
}
