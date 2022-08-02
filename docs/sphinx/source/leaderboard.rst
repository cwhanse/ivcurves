.. _leaderboard:

Leaderboard
===========

.. datatemplate:nodata::

   {{ make_list_table_from_mappings([
         ('Rank', 'rank'),
         ('Pull Request', 'pr_number'),
         ('Username', 'username'),
         ('Overall Score', 'overall_score'),
         ('Submission Date', 'submission_date')
      ],
      config.html_context.leaderboard.leaderboard_entries,
      title='')
   }}

