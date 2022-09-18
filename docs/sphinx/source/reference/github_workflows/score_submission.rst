Score Submission
================

.. mermaid::

   flowchart
   Push(A commit is pushed)
   WorkflowCall(Another workflow calls this workflow)

   Push --> StartWorkflow(Start this workflow)
   WorkflowCall --> StartWorkflow

   StartWorkflow --> StartReadPRConfigJob(Start the read-pr-config job)

   StartReadPRConfigJob --> CheckoutReadPRConfig(Checkout Competitor's pull request branch)
   CheckoutReadPRConfig --> InstallPython310ReadPRConfig(Install Python 3.10)

   InstallPython310ReadPRConfig --> ReadPRConfig(Read and validate the Competitor's pr_config.json)

   ReadPRConfig -->|RUN_SCORER is true| StartScoreSubmissionJob(Start the score-submission job)
   ReadPRConfig -->|RUN_SCORER is false| TerminateWorkflow(End the workflow)

   StartScoreSubmissionJob --> CheckoutScoreSubmission(Checkout Competitor's pull request branch)
   CheckoutScoreSubmission --> InstallPython310ScoreSubmission(Install Python 3.10)
   InstallPython310ScoreSubmission --> InstallIVCurvesDependencies(Install ivcurves Python dependencies)
   InstallIVCurvesDependencies --> InstallCompetitorDependencies(Install Competitor's Python dependencies)

   InstallCompetitorDependencies --> RunCompetitorSubmission(Run the Competitor's submission)
   RunCompetitorSubmission --> ScoreCompetitorOutput(Run ivcurve's compare_curves.py to score competitor's CSV output)

   ScoreCompetitorOutput --> ValidateScores(Validate the Competitor's scores)
   ValidateScores --> UploadOverallScores(Upload the scores to the workflow's artifacts)

