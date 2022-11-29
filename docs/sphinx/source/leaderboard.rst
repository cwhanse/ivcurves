.. _leaderboard:

Leaderboard
===========

Submissions are ranked by their overall score, which is the sum of the scores for each test set.
To see how the submissions rank for the each test set, see :ref:`compare_submissions`.

.. datatemplate:nodata::

   {{ make_list_table_from_mappings([
         ('Rank', 'rank'),
         ('Submission', 'submission'),
         ('Username', 'username'),
         ('Overall Score', 'overall_score'),
         ('Submission Date', 'submission_date')
      ],
      config.html_context.leaderboard.entries,
      title='')
   }}

