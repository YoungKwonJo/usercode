#!/bin/bash
source /afs/cern.ch/cms/sw/cmsset_default.sh
project CMSSW
cmsrel CMSSW CMSSW_3_1_1
cd CMSSW_3_1_1/src/
eval `scramv1 runtime -sh`
cvs co -r V00-00-04 UserCode/youngjo/copyFromCastorToWeb
cd UserCode/youngjo/copyFromCastorToWeb

./checkgif.sh


# example
# > bsub -q 1nd copyjob.sh
