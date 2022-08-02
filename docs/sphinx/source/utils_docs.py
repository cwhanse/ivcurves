"""
Wrapper of utils.py that is used in the ivcurves scripts
Allows useful functions in utils.py to be used by the Sphinx docs scripts.
"""
import os
import sys

IVCURVES = '../../../ivcurves'

sys.path.insert(0, os.path.abspath(IVCURVES))

import utils # utils.py used by ivcurves scripts


TEST_SETS_DIR = utils.TEST_SETS_DIR


def mp_nstr_precision_func(num_mpf):
    return utils.mp_nstr_precision_func(num_mpf)


def read_iv_curve_parameter_sets(filename):
    return utils.read_iv_curve_parameter_sets(filename)


def make_iv_curve_name(test_set_name, index):
    return utils.make_iv_curve_name(test_set_name, index)


def get_filenames_in_directory(directory_path):
    return utils.get_filenames_in_directory(directory_path)

