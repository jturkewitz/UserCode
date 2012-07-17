Instructions to run HSCP analysis
Seth I. Cooper
June 27 2012


Preliminaries:
(1) Obtain "edm" skim files for data and Monte Carlo from the standard HSCP analysis.
These are produced from the central HSCP skim [1] running on the data from the muon and MET
triggers, and contain the dE/dx and isolation collections.  The 2011 files are listed on the
twiki page [2].  I have copies of the files at Minnesota on the HDFS cluster [3].

[1]
http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/SUSYBSMAnalysis/Skimming/python/EXOHSCP_cff.py?revision=1.16&view=markup
[2] https://twiki.cern.ch/twiki/bin/viewauth/CMS/Hscp2011Analysis#EDM_files
[3] /hdfs/cms/user/cooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Nov11/

(2) Using CMSSW_4_2_8, check out the following tag list:
V00-04-07    AnalysisDataFormats/SUSYBSMObjects
V08-03-09    PhysicsTools/Utilities
V03-03-07    RecoLuminosity/LumiDB
V03-00-09    RecoMuon/TrackingTools
V00-06-00    RecoTracker/DeDx 
V00-04-19    SUSYBSMAnalysis/HSCP
V00-04-06    SUSYBSMAnalysis/HSCP/src/BetaCalculatorECAL.cc
V00-07-10    SUSYBSMAnalysis/Skimming

You will also need my package: UserCode/SCooper/HSCP2011/ToyMC.  The HEAD should be in good shape.
Compile.  To run the limit setting code, you will need a recent version of root, 5.32 or above.
At Minnesota, you will need to download and compile the source by hand (I had to disable HDFS
support or it wouldn't compile, and you have to enable RooFit)
in a directory of your choice.

(3) You will need to generate a sort of ntuple (in the format of RooDataSet) from the skim
files.  First, cd into the HSCP2011/ToyMC/bin folder.  Then, edit the SignalDefinitions.py*
file with the location of the signal skim files (at the bottom).  Then you can edit the
RunHSCParticlePlots.py file with the locations of the data skim files, and finally, execute the
RunHSCParticlePlots.py (python RunHSCParticlePlots.py), which will submit jobs to condor which
will run the makeHSCParticlePlots.cpp binary on all the signal and data input files.  This binary
will create output root files which contain the RooDataSets needed to do the analysis, along with
some plots.  The preselection is also applied at this step.  Each signal model has two input
files, one for the tight and one for the loose RPC trigger regimes in 2011.  The data files
(NOT signal sample files) from the output of the condor jobs should be combined into a single data
file using the binary mergeHSCParticlePlotsRooDataSets.cpp.  Edit this file with the appropriate
input file names, recompile, and then run it from the command line (mergeHSCParticlePlotsRooDataSets).
Inputs: EDM skim files from data and signal (2 for each signal model)
Outputs: root files containing the RooDataSets and plots

* This is a file that contains information for each model/mass, such as the standard analysis
counting experiment Pt/Ias cuts and the location of the signal files.  Groups of signals (such as
gluinos, staus, etc.) are defined as well as an overall modelList that is used in the Launch.py
script (see below).  New signal definitions can be added here.


Analysis:
(1) The analysis is controlled via the Launch.py* script and is divided into steps.
Each step can be launched by running Launch.py with the step number (e.g., python Launch.py 2).
Running Launch.py without a step number will give a help message.  The first step of the analysis
is to make the background prediction from the data sidebands.  The binary in charge here is called
makeIasPredictions.cpp.  It creates the A/B/C regions for each dE/dx number-of-measurements (NoM)
and eta slice and the predicted Ias distribution after the mass cut (bin-by-bin ABCD method).
The mass cut is unique to each (different mass) signal model, defined as the signal's nominal
mass - 2*sigma, where sigma is the expected mass resolution.  Therefore one background prediction
must be created for each unique mass.  This step also creates the histogram for the real data
in the search region.  The Launch.py script divides the prediction-making into a few slices per
job and submits the jobs to condor.  The individual prediction files for each mass cut are then
merged together in the next step.  The first invocation of Launch.py will create the directory
structure, with the top level directory name specified as the FarmDirectory variable.
Inputs: root files for data (not signal) containing RooDataSets
Outputs: one root file for each set of slices, containing: background prediction for limits,
background prediction for discovery, background prediction with exponential slope variation
of +/- 1 sigma; data observed in search region; other diagnostic plots

* This is the driver script for the analysis.  There are many things that must be properly
configured here, such as: FarmDirectory name, input root file for data (from step (3) above),
input root file name beginning for signal (for files create from step (3) above); integrated
luminosity; the templates of the python cfgs for each binary; the template XML files for the
HistFactory (see below); and the root macros for the limit setting machinery.  If running the
analysis from scratch, you will probably need to change the FarmDirectory and input root file
names for data and signal, and possibly the integrated luminosity.  The rest should be fine as is.
There are auxiliary python modules (HSCPMakeIasPredictionsLaunch, HSCPMakeScaledPredictionsLaunch,
HSCPDoLimitsLaunch, HSCPDoSignificanceLaunch) which have the lower-level code to divide up
the jobs, write the cfg and sh files for each job, and submit to the batch cluster for each of
these steps.

(2) Merge the prediction files for each mass point: python Launch.py 3.  You should see an error
if there are one or more files missing for any of the models/masses.
Inputs: root files from step (1)
Outputs: 1 root file for each model/mass, containing all the background predictions
(and data distributions) for all NoM/eta slices.

(3) Make the predictions for the signal.  This is done by invoking step 4 (python Launch.py 4)
which will submit jobs (one for each signal model that run the binary makeScaledPredictions.cpp.
This binary looks at the signal, does lumi reweighting for the two different trigger regimes in the
data, pileup reweighting, and application of Ias shape and NoM distribution data/MC systematics.
The signal predictions are all scaled to 1 pb for later usage with the limit setting tools.
All NoM/eta slices (for signal and background) are then combined into a single "unrolled"
histogram, each slice taking up its own set of bins.  This is necessary for the limit setting tools
as they are currently being used.  The final part of this step is the running of the HistFactory,
a part of RooStats [1][2] which takes XML files and the root files containing the histograms
and produces RooStats models which can be fed to the limit-setting macros.  This requires a
recent version of root, so you will need to edit HSCPMakeScaledPredictionsLaunch.py with the
location of your installation of root (if running at Minnesota).
Inputs: Signal skim files, output prediction root file from step (2)
Outputs: Root files containing signal and background predictions and data distributions;
root files containing RooStats models

[1] https://twiki.cern.ch/twiki/bin/view/RooStats/WebHome
[2] https://twiki.cern.ch/twiki/bin/view/CMS/RooStatsTprime

(4) Submit the limit-setting jobs.  You will need to edit HSCPDoLimitsLaunch.py with
the location of your installation of root (if running at Minnesota).  Each job runs the
StandardHypoTestInvDemo.C macro, which is a RooStats tutorial macro, which scans trial cross
section values and computes the expected/observed CLs (expected is from toy generation).
These are then combined together in the next step for difference cross section points or
to sum up the toys for an individual cross section point.  Each job can take hours to run.
Invoked with step 5 (python Launch.py 5).
Inputs: Root files with RooStats model from step (3)
Outputs: Root files with the HypoTestInverter information

(5) Combine the limit-setting jobs.  This is done via Launch.py with a standalone binary called
combineLimitResults.cpp.  I recommend updating the file setupRoot.sh with your installation
of root.  Then you can execute makeCombineLimitResults (source makeCombineLimitResults) which
will compile the standalone binary for you.  Once you have compiled the binary, you can invoke
step 6 (python Launch.py 6) to combine the limit results.  You can also make sure in the output
on the screen that the number of files being merged together is correct, and that a sensible
limit is computed.
Inputs: Root files with the HypoTestInverter information from step (4)
Outputs: Root files with all the HypoTestInverter information combined

(6) Submit the significance jobs.  This will compute the observed signal significance and the
expected significance for several trial cross sections.  The latter is later used to compute
the expected minimum cross section at which we can claim a 5-sigma discovery.  The significance
calculation is done using a different binary (StandardHypoTestDemo.C).  Invoke step 7 (python
Launch.py 7).  Since the asymptotic approximation is used here, the jobs run very quickly.
Inputs: Root files with the HypoTestInverter information from step (4)
Outputs: Root files with HypoTest information, one for each trial cross section

(7) Plot results and tables.  This is done using a standalone script called PlotResults.py,
where you must specify the farm directory you were using with Launch.py (here called BaseDir).
This script will create p-value plots for each signal model, along with text and latex tables
with some information and results, like that in the analysis note, and a plot of the mass limits,
taken from interpolating the limit results from each mass and finding the intersection with
the theoretical cross section.  The output is saved in a folder called results just underneath
the farm directory.  Note that it is not necessary to have both the significance results and
the limit results to run this script.
Inputs: Combined limit root files from step (5) and significance files from step (6)
Outputs: PNG, EPS formats of p-value of mass limit plots; text files containing tables in
text and latex format

