#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import os
import subprocess

__version__ = "0.0.0"

# Only import the module if not run from the setup script.
try:
    __FSPS_SETUP__
except NameError:
    __FSPS_SETUP__ = False


