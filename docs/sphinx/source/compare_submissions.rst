:orphan:

.. _compare_submissions:

Compare Submissions
===================

.. datatemplate:nodata::

   {{ make_list_table_from_mappings([
         ('Submission', 'submission'),
         ('Username', 'username'),
         ('Case 1', 'case1'),
         ('Case 2', 'case2'),
         ('Case 3', 'case3'),
         ('Case 4', 'case4'),
         ('Case 5', 'case5'),
         ('Case 6', 'case6'),
         ('Case 7', 'case7'),
         ('Case 8', 'case8')
      ],
      config.html_context.compare_submissions.entries,
      title='')
   }}

