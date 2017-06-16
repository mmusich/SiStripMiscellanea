#!/usr/bin/env python

import os,sys
import getopt
import commands
import time
import ROOT
import urllib
import string
from optparse import OptionParser

PCLDATASET = '/StreamExpressCosmics/Run2017A-SiStripPCLHistos-Express-v1/ALCARECO' 
runs = commands.getstatusoutput("dasgoclient -query='run dataset="+PCLDATASET+"'")[1].split(",")

#######################################################
def batchScriptCERN(runindex):
#######################################################
    '''prepare the LSF version of the batch script, to run on LSF'''
    script = """
#!/bin/bash 
CMSSW_DIR=$CMSSW_BASE/src
OUT_DIR=$CMSSW_DIR/harvest 
LOG_DIR=$CMSSW_DIR/out
LXBATCH_DIR=`pwd`  
cd $CMSSW_DIR
eval `scram runtime -sh`
cd $LXBATCH_DIR 
cp -pr $CMSSW_DIR/cfg/step5_ALCAHARVEST_{runindex}_cfg.py .
cmsRun step5_ALCAHARVEST_{runindex}_cfg.py >& jobLog.out
ls -lh . 
for payloadOutput in $(ls *db ); do cp $payloadOutput $OUT_DIR/Calibration_run{runindex}.db ; done 
for logOutput in $(ls *out ); do cp $logOutput $LOG_DIR/log_run{runindex}.out; done 
""".format(runindex=runindex)
   
    return script


##############################################
def main():
##############################################

    desc="""This is a description of %prog."""
    parser = OptionParser(description=desc,version='%prog version 0.1')
    parser.add_option('-s','--submit',    help='job submitted',    dest='submit',     action='store_true',  default=False)
    (opts, args) = parser.parse_args()

    if(not os.path.exists("cfg")):
        os.system("mkdir cfg")
        os.system("mkdir lsf")
        os.system("mkdir harvest")
        os.system("mkdir out")

    #count=0
    for run in runs:
        #count=count+1
        #if(count>10): 
        #    continue
        run = run.strip("[").strip("]")
        #print run
        files = commands.getstatusoutput("dasgoclient -query='file dataset="+PCLDATASET+" run="+run+"'")[1].split("\n")
        #print run, files
        listOfFiles='['
        for ffile in files:
            listOfFiles=listOfFiles+"\""+str(ffile)+"\","
        listOfFiles+="]"

        cwd = os.getcwd()
        os.system("cp step5_ALCAHARVEST.py ./cfg/step5_ALCAHARVEST_"+run+"_cfg.py")
        os.system("sed -i 's|XXX_FILES_XXX|"+listOfFiles+"|g' "+cwd+"/cfg/step5_ALCAHARVEST_"+run+"_cfg.py")

        lsfdir = os.path.join(cwd,"lsf")
        scriptFileName = os.path.join(lsfdir,"batchHarvester_"+run+".lsf")
        scriptFile = open(scriptFileName,'w')
        scriptFile.write(batchScriptCERN(run)) 
        scriptFile.close()
        os.system('chmod +x %s' % scriptFileName)

        if opts.submit:
            print "submit harvesting job for run",run                      
            os.system("bsub -o tmp.tmp -q cmscaf1nd "+scriptFileName) 

if __name__ == "__main__":        
    main()
