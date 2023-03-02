GitHub Workflows
================

.. mermaid::

   flowchart
   ContributorPushes(Maintainer or Competitor pushes code) --> RunScoreSubmissionPush(Run score-submission)
   RunScoreSubmissionPush --> ContributorPullRequest(A pull request is made)
   ContributorPullRequest --> PullRequestReviewAndMerge(The pull request is merged)
   PullRequestReviewAndMerge --> RunRecordScores(Run record-scores, calling score-submission and build-sphinx-docs)
   PullRequestReviewAndMerge --> RunScoreSubmissionMergePush(Run score-submission merger's submission)

.. toctree::
   :maxdepth: 2

   score_submission
   record_scores
   update_scores
   build_sphinx_docs
   pytest
   utils/index

