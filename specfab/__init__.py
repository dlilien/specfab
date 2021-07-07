#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Nichlas Rathmann and David Lilien
#
# Distributed under terms of the MIT license.

"""

"""

from ._specfabpy import specfabpy as sf


def a2_ij(nlm):
    return sf.get_a2_ij(nlm)
