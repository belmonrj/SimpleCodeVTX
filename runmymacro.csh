#!/bin/csh

set number = $1

source ~/.cshrc

cd /phenix/scratch/belmonrj/SimpleCodeVTX/

#root -b -l -q RunMyMacro.C\(\"Run_SimpleTreeSVXCNT.C\",\"simple_output.root\",1000,\"Run15pp200MBPro104\"\)
root -b -l -q RunMyMacro.C\(\"Run_SimpleTreeSVXCNT.C\",\"simple_output.root\",$number,\"Run15pp200MBPro104\"\)

