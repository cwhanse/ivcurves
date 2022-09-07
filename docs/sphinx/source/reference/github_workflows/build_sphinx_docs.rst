Build Sphinx Docs
=================

.. mermaid::

   flowchart
   CalledByAnotherWorkflow(Another workflow calls this workflow) --> StartWorkflow(Start this workflow)
   CalledByMaintainer(A Maintainer manually runs this workflow) --> StartWorkflow

   StartWorkflow --> CheckoutMain(Checkout ivcurves' main branch into the folder 'main')
   CheckoutMain --> InstallPython310(Install Python 3.10)
   InstallPython310 --> InstallIVCurvesDependencies(Install ivcurves Python dependencies)
   InstallIVCurvesDependencies --> GenerateTestCasePlotImages(Generate test case plot images)

   GenerateTestCasePlotImages --> BuildSphinx(Build the Sphinx HTML)

   BuildSphinx --> CheckoutGHPages(Checkout ivcurves' gh-pages branch into the folder 'gh-pages')

   CheckoutGHPages --> CopyBuildFiles(Copy build files in 'main' into 'gh-pages')

   CopyBuildFiles --> CommitNewBuild(Commit and push new build files to ivcurves' gh-pages branch)

   CommitNewBuild --> PublishGHPages(GitHub runs another workflow to publish the changes to the GitHub pages website)


