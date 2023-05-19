Update Scores
=============

.. mermaid::

   flowchart
   MaintainerCallsWorkflow(An Maintainer requests to update all submission scores) --> StartWorkflow(Start this workflow)

   StartWorkflow --> StartGetPRDataJob(Start the get-pr-data job)
   StartGetPRDataJob --> CheckoutGetPRData(Checkout ivcurves' main branch)
   CheckoutGetPRData --> InstallPython310GetPRData(Install Python 3.10)
   InstallPython310GetPRData --> ReadScoresDatabasePullRequestData(Read scores_database.json and create a 2-D array with these columns: pr_number, username, and submission_datetime)

   subgraph collect-pr-submissions
   CheckoutCollectPRSubmissions(Checkout ivcurves' main branch with its entire commit history) --> FindPRMergeCommit(Use GitHub CLI to find the merge commit of the pull request)
   FindPRMergeCommit --> CheckoutPRMergeCommit(Checkout the merge commit of the pull request)
   CheckoutPRMergeCommit --> RenameSubmissionFolder(Rename the pull request author's submission to pr_number--username)
   RenameSubmissionFolder --> UploadRenamedSubmissionFolder(Upload the renamed submission folder to GitHub artifacts)
   end

   ReadScoresDatabasePullRequestData -->|For every pull requeset number in scores_database.json| collect-pr-submissions
   collect-pr-submissions --> StartScoreAllSubmissionsJob(Start the score-all-submissions job)

   StartScoreAllSubmissionsJob --> CheckoutScoreAllSubmissions(Checkout ivcurves' main branch)
   CheckoutScoreAllSubmissions --> InstallPython310ScoreAllSubmissions(Install Python 3.10)
   InstallPython310ScoreAllSubmissions --> InstallIVCurvesDependencies(Install ivcurves Python dependencies)
   InstallIVCurvesDependencies --> DeleteAllSubmissions(Delete all submissions to make room for the ones uploaded by the collect-pr-submissions job)
   DeleteAllSubmissions --> DownloadGitHubArtifactsSubmissions(Download all submissions uploaded by the collect-pr-submissions job)
   DownloadGitHubArtifactsSubmissions --> BeginBashScriptToRunAllSubmissions(Begin a Bash script to run and score all submissions)

   subgraph Bash
   BashCreateVirtualEnv(Bash creates a Python virtual environment for the submission) --> BashInstallSubmissionDependencies(Bash installs the submission's Python dependencies)
   BashInstallSubmissionDependencies --> BashRunSubmission(Bash runs the submission)
   BashRunSubmission --> BashScoreSubmission(Bash scores the submission)
   BashScoreSubmission --> BashValidateRecordScores(Bash validates and records the scores, marking the submission broken if validation fails)
   BashValidateRecordScores --> BashRemoveVirtualEnv(Bash removes the virtual environment for the submission)
   end

   BeginBashScriptToRunAllSubmissions -->|For every submission| Bash
   Bash --> CommitModifiedDatabase(Commit and push the udpated database to GitHub)

