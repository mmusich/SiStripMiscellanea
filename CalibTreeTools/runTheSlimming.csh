#!/bin/tcsh

set COUNT=0
foreach inputfile (`eos ls /eos/cms/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR17_Aag/`)
    echo "=================================================" 
    if ($inputfile =~ *"root"*) then
	set myLFN="root://eoscms.cern.ch//eos/cms/store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR17_Aag/$inputfile"
	echo $myLFN
	root -q -b .L treeCloner.C++\(\"$myLFN\"\)
	@ COUNT +=1
    endif
    if($COUNT > 100) then
	break
    endif
end

echo "Total checked files: $COUNT"
