Test Cases
==========

IV curve fitting algorithms are scored using three benchmark cases.

Case 1 and Case 2
-----------------

Case 1 and Case 2 represent, respectively, a module composed of 72 series-connected single junction crystalline silicon (cSi) cells, and a module composed of 140 series-connected thin-film cells. Each case contains 32 sets of IV curves computed exactly for all combinations of high and low values of the five parameters in the single diode equation. Tables of the parameters are given below.

For both Case 1 and Case 2, the scoring engine expects to recieve the five parameter values for each IV curve. The test score is the sum over all IV curves of a metric that measures the separation between the exactly-known curve and the curve computed from the inputted five parameters. The metric is the sum of the length of the line segments lying between each pair of curves the rays projected from the origin through 100 points on the known curve, where these points are equally-spaced between 0 V and V\ :sub:`OC`\ .

Case 3 comprises 4 sets of 50 IV curves each. Case 3a represents a good quality, 60-cell cSe module at high power condition, and Case 3b represents a poor quality cSi module in low light conditions. Case 3c represents a good quality, 140-cell thin film module at high power conditions, and Case 3d represents a poor quality thin-film module at low light conditions. For each of these four cases, 50 IV curves are formed by computing the exact curve from the assumed parameters, then adding simulated measurement noise. Voltage noise is simulated by uncorrelated random samples from a uniform distribution on [-0.05%, 0.05%]. Current noise is generated from a zero-mean normal distribution with 0.1% standard deviation, with an autocorrelation that decays exponentially. Details are provided in the file build_case3.py.

For Case 3, the scoring engine expects to recieve one set of five parameter values for each part of Case 3 (four sets in total). The test score is the sum (over the parts of Case 3) of the sum (over the five parameters) of the absolute difference in parameter value, as a percent of the known value.


.. datatemplate:nodata::

   {% for test_set_name, data in config.html_context.test_cases.test_case_data.items() %}

   {{ test_set_name }}
   =============================

   {{ make_list_table([
         'Index',
         'Photocurrent (A)',
         'Saturation Current (A)',
         'Series Resistance (Ω)',
         'Shunt Resistance (Ω)',
         'Diode Factor (-)'
      ],
      data,
      title='')
   }}

   {% endfor %}

