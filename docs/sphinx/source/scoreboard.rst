.. _scoreboard:

Scoreboard
==========

Submissions are given a score for some or all test sets, and the sum of these scores is the submission's overall score.
If a submission is not scored on a test set, that test set's score will be blank (--).
Test sets **case1** and **case2** are scored by the distance between the known IV curve and the submission's fitted IV curve (see :func:`ivcurves.compare_curves.score_curve`).
Test sets **case3a** through **case3d** are scored by the difference between the known and fitted single diode equation parameters (see :func:`ivcurves.compare_curves.score_parameters`).

.. datatemplate:nodata::

   .. list-table::
      :header-rows: 1

   {% for row in config.html_context.scoreboard.table_rows %}
      * - {{ row[0] }}
        {% for val in row[1:] %}
        - {{ val }}
        {% endfor %}
   {% endfor %}

