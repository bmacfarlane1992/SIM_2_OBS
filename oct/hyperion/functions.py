from __future__ import print_function, division

import string
import random
import os
import shutil
import sys
import tempfile

import h5py
import numpy as np

from astropy import log as logger

MAX_FLOAT = np.log(np.finfo('d').max)

def link_or_copy(group, name, link, copy, absolute_paths=False):
    '''
    Link or copy a dataset or group

    Parameters
    ----------
    group : h5py.Group
        The group to create the link, dataset, or group in
    name : str
        The name of the link, dataset, or group in the new file
    link : h5py.ExternalLink
        A link to the group or dataset to include
    copy : bool
        Whether to copy or link to the dataset
    absolute_paths : bool
        If copy=False, then if absolute_paths is True, absolute filenames
        are used in the link, otherwise the path relative to the file is used.
    '''
    if copy:
        f = h5py.File(link.filename, 'r')
        f.copy(link.path, group, name=name)
        f.close()
    else:
        if absolute_paths:
            group[name] = h5py.ExternalLink(os.path.abspath(link.filename), link.path)
        else:
            group[name] = h5py.ExternalLink(os.path.relpath(link.filename, os.path.dirname(group.file.filename)), link.path)
        try:
            group[name]
        except KeyError:  # indicates linking failed (h5py < 2.1.0)
            logger.warn("Linking failed, copying instead (indicates an outdated version of h5py)")
            del group[name]
            f = h5py.File(link.filename, 'r')
            f.copy(link.path, group, name=name)
            f.close()


class FreezableClass(object):

    _frozen = False
    _final = False
    _attributes = []

    def __init__(self):
        super(FreezableClass, self).__init__()

    def _freeze(self):
        object.__setattr__(self, '_frozen', True)

    def _finalize(self):
        object.__setattr__(self, '_final', True)

    def isfrozen(self):
        return self._frozen

    def isfinal(self):
        return self._final

    def __setattr__(self, key, value):
        if self._final:
            raise Exception("Attribute %s can no longer be changed" % key)
        if self._frozen and not key in self._attributes:
            raise AttributeError("Attribute %s does not exist" % key)
        self._attributes.append(key)
        object.__setattr__(self, key, value)

def is_numpy_array(variable):
    return issubclass(variable.__class__, (np.ndarray,
                                           np.core.records.recarray,
                                           np.ma.core.MaskedArray))
