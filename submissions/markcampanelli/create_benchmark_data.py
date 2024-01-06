"""Run unit test from pvfit package to generate fit data."""

import pvfit.tests.modeling.dc.single_diode.equation.inference_test
import pytest

test_file = pvfit.tests.modeling.dc.single_diode.equation.inference_test.__file__

pytest.main(["-vvsx", "--cache-clear", "-k", "not test_fit_benchmark", test_file])
