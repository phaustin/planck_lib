# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:light
#     notebook_metadata_filter: all,-language_info,-toc,-latex_envs
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.1
#   latex_metadata:
#     chead: center head
#     lhead: left head
# ---

import numpy as np
from scipy.integrate import quad
import pdb

clight = 2.99792458e+08  #m/s -- speed of light in vacumn
h = 6.62606876e-34  #J s  -- Planck's constant
kb = 1.3806503e-23  # J/K  -- Boltzman's constant
c1 = 2. * h * clight**2.
c2 = h * clight / kb
c3 = 2 * h / clight**2.
c4 = h / kb

sigma = 2. * np.pi**5. * kb**4. / (15 * h**3. * clight**2.)


def planckDeriv(wavel, Temp):
    """
       input: wavel in m, Temp in K
       output: dBlambda/dlambda  W/m^2/m/sr/m
    """
    expterm = np.exp(c2 / (wavel * Temp))
    deriv = c1 / np.pi * wavel**(-6.) * (expterm -
                                         1)**(-2.) * c2 / Temp**2. * expterm
    return deriv


def planckwavelen(wavel, Temp):
    """
       input: wavelength (m), Temp (K)
       output: planck function W/m^2/m/sr
    """
    Blambda = c1 / (wavel**5. * (np.exp(c2 / (wavel * Temp)) - 1))
    return Blambda


def planckfreq(freq, Temp):
    """
      input: freq (Hz), Temp (K)
      output: planck function in W/m^2/Hz/sr
    """
    Bfreq = c3 * freq**3. / (np.exp(c4 * freq / Temp) - 1)
    #pdb.set_trace()
    return Bfreq


def planckwavenum(waven, Temp):
    """
      input: wavenumber (m^{-1}), Temp (K)
      output: planck function in W/m^2/m^{-1}/sr
    """
    Bwaven = c1 * waven**3. / (np.exp(c2 * waven / Temp) - 1)
    return Bwaven


def planckInvert(wavel, Blambda):
    """input wavelength in m and Blambda in W/m^2/m/sr, output
    output brightness temperature in K
    """
    Tbright = c2 / (wavel * np.log(c1 / (wavel**5. * Blambda) + 1.))
    return Tbright


def planckInt(planck_fn, Temp, lower, upper):
    """Integrate planckwavelen given temperatue Temp (K) from lower (m) to upper (m) wavelengths

       output: integrated radiance in W/m^2/sr
       see http://docs.scipy.org/doc/scipy-0.14.0/reference/integrate.html#module-scipy.integrate
    """
    args = (Temp)
    integ = quad(planck_fn, lower, upper, args)
    return integ[0]


def goodInvert(T0, bbr, wavel):
    B0 = planckwavelen(wavel, T0)
    theDeriv = planckDeriv(wavel, T0)
    delB = bbr - B0
    delT = delB / theDeriv
    theT = T0 + delT
    return theT


def rootfind(T0, bbrVec, wavel):
    bbrVec = np.asarray(bbrVec)
    guess = planckwavelen(T0, wavel)
    out = []
    for bbr in bbrVec:
        while np.fabs(bbr - guess) > 1.e-8:
            delB = bbr - guess
            deriv = planckDeriv(wavel, T0)
            delT = delB / deriv
            T0 = T0 + delT
            guess = planckwavelen(wavel, T0)
        out.append(T0)
    return out


def test_planck_wavelen():
    """
       test planck function for several wavelengths
       and Temps
    """
    #
    # need Temp in K and wavelen in m
    #
    the_temps = [200., 250., 350.]
    the_wavelens = np.array([8., 10., 12.]) * 1.e-6
    out = []
    for a_temp in the_temps:
        for a_wavelen in the_wavelens:
            #
            # convert to W/m^2/micron/sr
            #
            the_bbr = planckwavelen(a_wavelen, a_temp) * 1.e-6
            out.append(the_bbr)
    answer = [
        0.4521, 0.8954, 1.1955, 2.7324, 3.7835, 3.9883, 21.4495, 19.8525,
        16.0931
    ]
    np.testing.assert_array_almost_equal(out, answer, decimal=4)
    return None


def test_planck_wavenum():
    """
       test planck function for several wavelengths
       and Temps
    """
    #
    # need Temp in K and wavelen in m
    #
    the_temps = [200., 250., 350.]
    the_wavelens = np.array([8., 10., 12.]) * 1.e-6
    the_wavenums = 1 / the_wavelens
    out = []
    for a_temp in the_temps:
        for a_wavenum in the_wavenums:
            #
            # convert to W/m^2/micron/sr
            #
            the_bbr = planckwavenum(a_wavenum, a_temp)
            out.append(the_bbr)
    answer = [
        2.8932e-05, 8.9535e-05, 1.7215e-04, 1.7487e-04, 3.7835e-04, 5.7431e-04,
        1.3728e-03, 1.9852e-03, 2.3174e-03
    ]
    np.testing.assert_array_almost_equal(out, answer, decimal=4)
    return None


def test_planck_freq():
    """
       test planck function for several wavelengths
       and Temps
    """
    #
    # need Temp in K and wavelen in m
    #
    the_temps = [200., 250., 350.]
    the_wavelens = np.array([8., 10., 12.]) * 1.e-6
    the_wavenums = 1 / the_wavelens
    the_freqs = the_wavenums * clight
    out = []
    for a_temp in the_temps:
        for a_freq in the_freqs:
            the_bbr = planckfreq(a_freq, a_temp)
            out.append(the_bbr)
    answer = [
        9.6508e-14, 2.9866e-13, 5.7424e-13, 5.8331e-13, 1.2620e-12, 1.9157e-12,
        4.5791e-12, 6.6221e-12, 7.7300e-12
    ]
    np.testing.assert_array_almost_equal(out, answer, decimal=4)
    return None


def test_planck_inverse():
    """
       test planck inverse for several round trips
       and Temps
    """
    #
    # need Temp in K and wavelen in m
    #
    the_temps = [200., 250., 350.]
    the_wavelens = np.array([8., 10., 12.]) * 1.e-6
    out = []
    for a_temp in the_temps:
        for a_wavelen in the_wavelens:
            #
            # convert to W/m^2/micron/sr
            #
            the_bbr = planckwavelen(a_wavelen, a_temp)
            out.append((a_wavelen, the_bbr))

    brights = []
    for wavelen, bbr in out:
        brights.append(planckInvert(wavelen, bbr))
    answer = [200.0, 200.0, 200.0, 250.0, 250.0, 250.0, 350.0, 350.0, 350.0]
    np.testing.assert_array_almost_equal(brights, answer, decimal=10)
    return None


def test_planck_wavelen_integral():
    """
       integrage and compare with stefan-boltzman
    """
    Temp = 300.
    stefan = sigma / np.pi * Temp**4.
    totrad = planckInt(planckwavelen, Temp, 1.e-7, 8000.e-6)
    np.testing.assert_almost_equal(totrad, stefan, decimal=5)
    return None


def test_planck_wavenum_integral():
    """
       integrage and compare with stefan-boltzman
    """
    Temp = 300.
    stefan = sigma / np.pi * Temp**4.
    left = 1. / 8000.e-6
    right = 1 / 1.e-7
    totrad = planckInt(planckwavenum, Temp, left, right)
    np.testing.assert_almost_equal(totrad, stefan, decimal=5)
    return None


def test_planck_freq_integral():
    """
       integrage and compare with stefan-boltzman
    """
    Temp = 300.
    stefan = sigma / np.pi * Temp**4.
    left = (1. / 8000.e-6) * clight
    right = (1 / 1.e-7) * clight
    print(left, right)
    totrad = planckInt(planckfreq, Temp, left, right)
    np.testing.assert_almost_equal(totrad, stefan, decimal=5)
    return None

# this trick will run  the following script if
# the file planck.py is run as a program, but won't
# if  planck.py is imported from another  module


if __name__ == '__main__':

    test_planck_wavelen()
    test_planck_wavenum()
    test_planck_freq()
    test_planck_inverse()
    test_planck_wavelen_integral()
    test_planck_wavenum_integral()
    test_planck_freq_integral()
