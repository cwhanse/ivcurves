GitHub Workflows
================

.. mermaid::

   flowchart
   ContributorPushes(Maintainer or Competitor pushes code) --> RunScoreSubmissionPush(Run score-submission)
   RunScoreSubmissionPush --> ContributorPullRequest(A pull request is made)
   ContributorPullRequest --> PullRequestReviewAndMerge(The pull request is merged)
   PullRequestReviewAndMerge --> RunRecordScores(Run record-scores for the merged pull request)
   PullRequestReviewAndMerge --> RunScoreSubmissionMergePush(Run the submission of the maintainer who merged the pull request)

.. toctree::
   :maxdepth: 2

   score_submission
   record_scores
   update_scores
   pytest
   utils/index

