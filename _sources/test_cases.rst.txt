Test Cases
==========

.. datatemplate:nodata::

   {% for test_set_name, test_case in config.html_context.test_cases.test_case_data.items() %}

   Test Set: {{ test_set_name }}
   =============================

   {% for Index, data in test_case.items() %}

   {{ data.title }}
   ----------------

   {{ make_list_table([
         'Photocurrent (A)',
         'Saturation Current (A)',
         'Series Resistance (Ω)',
         'Shunt Resistance (Ω)',
         'Diode Factor (-)'
      ],
      [data.parameters],
      title='')
   }}

   .. image:: {{ data.image_path }}
      :width: 600

   {% endfor %}
   {% endfor %}

