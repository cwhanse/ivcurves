.. _participating:

How to Participate
==================


Submission Workflow
-------------------

These are the steps to make a submission as a competitor.
You may as many submissions as you like, and each will be a separate entry on the leaderboard.

#. Create a fork of the ivcurves repository.
#. Create a new folder ``submissions/<your_GitHub_username>``.
   This folder will store all of the files you create.
#. Write or copy your submission's code inside ``submissions/<your_GitHub_username>``.
#. Create a file ``pr_config.json`` with three entries:

   .. code-block:: json

      {
          "RUN_SCORER": true,
          "REQUIREMENTS": "<path_to_requirements.txt>",
          "SUBMISSION_MAIN": "<path_to_submission_entrypoint.py>"
      }

   ``REQUIREMENTS`` and ``SUBMISSION_MAIN`` must be relative paths to ``submissions/<your_GitHub_username>``.
   For example, if your Pip requirements file is located at ``submissions/<your_GitHub_username>/requirements.txt``, ``REQUIREMENTS`` should be ``./requirements.txt``.
#. Push your changes to your fork, and then GitHub will automatically run the scorer with your code.
   When it is finished, you will see either a green check mark or red x icon next to your commit indicating whether the scorer succeeded or not.
   Click on the icon  to see the results of the scorer.
   The scorer will also provide a CSV file containing your code's score.

   .. note::
      GitHub will execute your ``SUBMISSION_MAIN`` from the same folder that contains it.
      Therefore, your code may use relative paths when reading and writing files.

#. When you are ready, and there is a green check mark next to your latest commit, create a pull request into ``cwhanse/ivcurves/main``.
   GitHub will again automatically run and score your code.
   An ivcurves maintainer will be notified of your pull request.
#. After an ivcurves maintainer reviews and approves your pull request, it will be merged.
   GitHub will post your GitHub username and score to the leaderboard.

Submission Workflow Diagram
---------------------------

.. mermaid::

   sequenceDiagram
       participant Competitor
       participant GitHub Actions
       participant ivcurves Maintainer

       activate Competitor
       Note over Competitor: Forks ivcurves repository.
       Note over Competitor: Creates folder submissions/<your_GitHub_username> to put their files into.
       Note over Competitor: Creates pr_config.json.
       Note over Competitor: Writes code for submission.

       Competitor->>GitHub Actions: Pushes code to fork on GitHub
       deactivate Competitor
       activate GitHub Actions
       Note over GitHub Actions: Reads pr_config.json and scores submission.
       GitHub Actions->>Competitor: Reports results with CSV of scores
       deactivate GitHub Actions

       Competitor->>GitHub Actions: Creates pull request
       activate GitHub Actions
       Note over GitHub Actions: Reads pr_config.json and scores submission code.
       GitHub Actions->>Competitor: Reports results with CSV of scores
       GitHub Actions->>ivcurves Maintainer: Notifies of pull request
       deactivate GitHub Actions
       activate ivcurves Maintainer
       Note over ivcurves Maintainer: Reviews pull request.
       ivcurves Maintainer->>GitHub Actions: Merges competitor's pull request
       activate GitHub Actions
       deactivate ivcurves Maintainer
       GitHub Actions->>Competitor: Notifies of pull request merge
       Note over GitHub Actions: Records submission's score in database.
       Note over GitHub Actions: Updates leaderboard and submission documentation.
       deactivate GitHub Actions


Submission Requirements
-----------------------

In the ``test_sets`` folder, there are multiple JSON files whose contents has the following structure:

.. code-block::

  {
    "Manufacturer": "",
    "Sandia ID": "",
    "Material": "",
    "IV Curves": [
      {
        "Index": 1,
        "Voltages": [
          "0.0",
          ...
        ],
        "Currents": [
          "0.9996667777132812",
          ...
        ],
        "v_oc": "39.7481074783976643",
        "i_sc": "0.9996667777132812",
        "v_mp": "33.9368943991155301",
        "i_mp": "0.8461238606639279",
        "p_mp": "28.7148161079236649",
        "Temperature": "298.15",
        "Irradiance": null,
        "Sweep direction": null,
        "Datetime": null
      },
      ...
    ]
  }

Under the ``"IV Curves"`` key is a list of IV curve data sets each with an ``"Index"`` value.
The ``"Index"`` value is the test case number of the test set.

For each JSON file ``<test_set_name>.json`` in ``test_sets``, your code must write a CSV file ``<test_set_name>.csv`` in ``submissions/<your_GitHub_username>``.
Each CSV file must have these columns:

.. datatemplate:nodata::

  {{ make_list_table([
        'Index',
        'photocurrent',
        'saturation_current',
        'resistance_series',
        'resistance_shunt',
        'n',
        'cells_in_series'
     ],
     [['#','#','#','#','#','#','#']],
     title='<test_set_name>.csv')
  }}

Each row the CSV file will contain your code's fitted parameters for each test case in its corresponding test set.


Here is some Python code that may be useful for getting a set all of the JSON filenames in ``test_sets``:

.. code-block:: python

   import json
   import pathlib


   def get_test_set_filenames():
       path_to_test_sets = pathlib.Path('../../test_sets')
       return {f'{path_to_test_sets}/{entry.stem}.json' for entry in path_to_test_sets.iterdir()
                   if entry.is_file()}


   def json_file_to_dict(filepath):
       with open(filepath, 'r') as file:
           return json.load(file)


Documenting Your Submission
---------------------------

These steps will cover how to add your submission to the ivcurves Submissions documentation.
The ivcurves documentation uses numpy-sytle docstrings.

#. In the ``docs/sphinx/source/submissions`` folder, make a new folder ``<your_GitHub_username>``.
   All documentation files you create will go in this folder.
#. For each of your submission's ``.py`` files in the top level of the ``submissions/<your_GitHub_username>`` folder, create a file ``<your_py_filename.py>.rst`` containing the following:

   .. |autosummary| replace:: autosummary

   .. code-block:: rst
      :substitutions:


      .. currentmodule:: submissions.<your_GitHub_username>.<your_py_filename>

      .. |autosummary|::
         :toctree: generated/

         ..
            write the name of each function in <your_py_filename>.py

         <function_name1>
         <function_name2>

   ..
      note to documentation writer: the rst in the code-block above
      is still interpreted by Sphinx. To prevent autosummary from executing,
      it must be substituted in (using sphinx_substitution_extensions)

#. The following steps are for registering your submission's ``.py`` files that are in subfolders under ``submissions/<your_GitHub_username>``.

   #. Create a folder ``<your_subfolder_name>``. This will contain all the documentation files you create in this set of steps.
   #. Inside that folder, for each ``.py`` file under ``submissions/<your_GitHub_username>/<your_subfolder_name>`` create a file ``<your_py_filename.py>.rst``.
   #. Create a file ``index.rst`` containing the following:

      .. code-block:: rst

         <your_subfolder_name>
         =====================

         .. toctree::
            :maxdepth: 2

            ..
               write the name of each .rst file you created here
               the .rst extension should be ommitted

            <your_py_filename1>
            <your_py_filename2>

#. Back in ``docs/sphinx/source/submissions/<your_GitHub_username>``, create a file ``index.rst`` containing the following:

   .. code-block:: rst

      <your_GitHub_username>
      ======================

      .. toctree::
         :maxdepth: 2

         ..
            write the name of each .rst file you created for your .py files in the top level of ``submissions/<your_GitHub_username>``
            the .rst extension should be ommitted

         <your_py_filename1>
         <your_py_filename2>

         ..
            for each subfolder in ``submissions/<your_GitHub_username>``, write the following lines

         <your_subfolder_name1>/index
         <your_subfolder_name2>/index

#. If your submission had a folder structure like this

   .. code-block:: bash

      submissions/<your_GitHub_username>
          |- pr_config.json
          |- requirements.txt
          |- <your_py_filename1>.py
          |- <your_subfolder_name1>/
               |- <your_py_filename1>.py

   then after following the previous steps you should have this folder structure:

   .. code-block:: bash

      docs/sphinx/source/submissions/
        |- index.rst
        |- <your_GitHub_username>/
             |- index.rst
             |- <your_py_filename1>.rst
             |- <your_subfolder_name1>/
                  |- index.rst
                  |- <your_py_filename1>.rst

#. Finally, inside ``submissions/index.rst`` like in the highlighted line:

   .. code-block:: rst
      :emphasize-lines: 8

      Submissions
      ===========

      .. toctree::
         :maxdepth: 2

         <other_GitHub_username1>/index
         <your_GitHub_username>/index

To help describe or contextualize your code, you may create links to external sites using this Sphinx rst directive in your docstrings:

   .. code-block:: rst

      Link to an `external site`_.

      .. _external site: <url>

      ..
         Example:

      Link to `Sphinx documentation`_.

      .. _Sphinx documentation: https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html#hyperlinks

