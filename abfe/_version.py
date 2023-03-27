#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# We will use semantic version (major, minor, patch)
__main_version_tuple__ = (0, 0, 0)
__pre_version_tuple__ = ('alpha',1) # or None or an empty tuple, the first char is alpha, beta, rc, etc...

if __pre_version_tuple__:
    __version_tuple__ = tuple(list(__main_version_tuple__) + list(__pre_version_tuple__))
    __version__ = '.'.join([str(i) for i in __main_version_tuple__]) + f'{__pre_version_tuple__[0][0]}' + '.'.join([str(i) for i in __pre_version_tuple__[1:]])
else:
    __version_tuple__ = __main_version_tuple__
    __version__ = '.'.join([str(i) for i in __main_version_tuple__])