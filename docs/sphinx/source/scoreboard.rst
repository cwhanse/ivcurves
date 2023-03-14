.. _scoreboard:

Scoreboard
==========

.. datatemplate:nodata::

   .. list-table::
      :header-rows: 1

   {% for row in config.html_context.scoreboard.table_rows %}
      * - {{ row[0] }}
        {% for val in row[1:] %}
        - {{ val }}
        {% endfor %}
   {% endfor %}

