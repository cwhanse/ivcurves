Record Scores
=============

.. mermaid::

   flowchart
   PullRequestMerged(The Competitor's pull request is merged into ivcurves' main branch) --> StartWorkflow(Start this workflow)

   StartWorkflow --> CallScoreSubmission(Call the score-submission workflow)
   CallScoreSubmission --> ReadPRConfig(Read score-submission's run-scorer output variable)

   ReadPRConfig -->|run-scorer is true| StartRecordScoresJob(Start the record-scores job)
   ReadPRConfig -->|run-scorer is false| TerminateWorkflow(End the workflow)

   StartRecordScoresJob --> CheckoutRecordScores(Checkout ivcurves' main branch)
   CheckoutRecordScores --> InstallPython310RecordScores(Install Python 3.10)
   InstallPython310RecordScores --> DownloadOverallScores(Download overall scores from the score-submission workflow)
   DownloadOverallScores --> ValidateAndRecordScores(Validate the overall scores CSV and record the scores to the database)

   ValidateAndRecordScores --> CommitModifiedDatabase(Commit and push updated database to GitHub)

