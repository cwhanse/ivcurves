"""Run unit test from pvfit package to generate fit data in the current directory."""

import pvfit.modeling.dc.single_diode.equation.inference_test
import pytest

if __name__ == "__main__":
    test_file = pvfit.modeling.dc.single_diode.equation.inference_test.__file__
    pytest.main(["-vvsx", "--cache-clear", f"{test_file}::test_fit"])
