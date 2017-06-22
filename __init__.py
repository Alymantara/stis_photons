# -*- coding: utf-8 -*-

__version__ = "0.1.0.dev0"
__author__ = "Juan V. Hernandez-Santisteban (pantro@gmail.com)"
__copyright__ = "Copyright 2016-2017 Juan v. Hernandez-Santisteban"
__contributors__ = [
    # Alphabetical by first name.
    "Just me @alymantara",
]


try:
    __STIS_PHOTONS_SETUP__
except NameError:
    __STIS_PHOTONS_SETUP__ = False

if not __STIS_PHOTONS_SETUP__:
    __all__ = ["localise", "hist2d", "quantile"]

    from .stis_photons import localise
