Test Cases
==========

Here we describe the test cases.


.. datatemplate:nodata::

   {% for test_set_name, data in config.html_context.test_cases.test_case_data.items() %}

   Test Set: {{ test_set_name }}
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

