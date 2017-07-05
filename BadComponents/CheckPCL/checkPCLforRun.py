#!/usr/bin/python

'''Script to check CMS SiStrip BadComponents PCL 
'''

__author__ = 'Suvankar Roy Chowdhury'
__copyright__ = 'Copyright 2017, CERN CMS'
__license__ = 'Unknown'
__maintainer__ = 'Suvankar Roy Chowdhury'
__email__ = 'suvankar.roy.chowdhury@cern.ch'
__version__ = 1


import ConfigParser
import glob
import os
import ROOT
import string
import subprocess
import sys
import time
import json
import codecs
import datetime
import getpass
from datetime import datetime
from ROOT import *

def getGridCertificate():
    vinfo = subprocess.check_output(['voms-proxy-info'])
    #Without the extra splitting at the end, the file couldn't be parsed by subprocess calls
    return vinfo.split('\n')[5].split(":")[1].split(" ")[1]

def colorText(val):
    return '\033[1;31m'+str(val)+'\033[1;m' if not val else '\033[1;32m'+str(val)+'\033[1;m'

def checkLogfilesforRun(runNumber, sfLog):
    flistAll = subprocess.Popen(("ls", "/eos/cms/store/express/tier0_harvest/"), stdout=subprocess.PIPE)
    flistRun=[]
    try:
        flistRun = subprocess.check_output(('grep', runNumber), stdin=flistAll.stdout)
    except subprocess.CalledProcessError, e:
        flistRun=[]
    hasdb=False
    hasTag=False
    insertTime="\tNA\t\t"
    if flistRun:
        for f in flistRun.split('\n'):
            if not hasdb:
                if ".db" and "SiStripBadStrip_pcl" in f: 
                    hasdb=True
                    ##check if new cond created
                    checkDBoutput = subprocess.check_output(["conddb","--db", "/eos/cms/store/express/tier0_harvest/"+f,"listTags"], stderr=FNULL)
                    #print checkDBoutput
                    #print len(checkDBoutput.split("\n"))
                    if len(checkDBoutput.split("\n")) > 3:
                        if "SiStripBadStrip_pcl" in checkDBoutput.split("\n")[2]:
                            hasTag=True
                            insertTime=checkDBoutput.split("\n")[2].split(" ")[38]+","+checkDBoutput.split("\n")[2].split(" ")[39]
    sfLog['hasdb'] = hasdb
    sfLog['hasPayload'] = hasTag
    sfLog['insertTime'] = insertTime


def getTimeElapsedFromRunEnd(run):
    out = subprocess.check_output(["das_client", "--limit", "0", "--query", "run={} | grep run.end_time".format(run)])
    if out == "[]\n":
        return 0
    else:
        result = out[:-1].replace('"', '')
        if result != 'null':
            #print result,'\t',datetime.utcnow()
            endTime = datetime.strptime(result,"%Y-%m-%d %H:%M:%S")
            delta = datetime.utcnow() - endTime
            return round(delta.days*24 + delta.seconds/3600.,1)

def setRunDirectory(runNumber):
    dirDict = { 294644:['Data2017', 'Run2017'],\
                290123:['Data2017', 'Commissioning2017'],\
                284500:['Data2016', 'PARun2016'],\
                271024:['Data2016', 'Run2016'],\
                264200:['Data2016', 'Commissioning2016'],\
                246907:['Data2015', 'Run2015'],\
                232881:['Data2015', 'Commissioning2015'],\
                211658:['Data2013', 'Run2013'],\
                209634:['Data2013', 'HIRun2013'],\
                190450:['Data2012', 'Run2012']}
    runKey=0
    for key in sorted(dirDict):
        if runNumber > key:
            runKey=key
    return dirDict[runKey]     

def checkDQMhisto(run):
    runDirval=setRunDirectory(run)
    DataLocalDir=runDirval[0]
    DataOfflineDir=runDirval[1]
    nnn=str(run)[:4]
    uname = getpass.getuser()
    tpath = "/tmp/" + uname + "/"  
    File_Name = ''  
    runClasses=["StreamExpress", "StreamExpressCosmics", "StreamHIExpress"]
    for Run_type in runClasses: 
        url = 'https://cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/'+DataOfflineDir+'/'+Run_type+'/000'+str(nnn)+'xx/'
        subprocess.check_output(["curl", "-k", "--cert", proxyName, "--key", proxyName, "-X", "GET", url,\
        "-o", tpath + "index.html"], stderr=FNULL)
        f=codecs.open(tpath + "index.html", 'r')
        index = f.readlines()
        for s in range(len(index)): 
            if str(run) in index[s]:
                if str("__DQMIO.root") in index[s]:
                    File_Name=str(str(index[s]).split("xx/")[1].split("'>DQM")[0])
                    break
        if File_Name:
           break
    #print run, 'Downloading DQM file:'+url+File_Name
    if File_Name:
        subprocess.check_output(["curl", "-k", "--cert", proxyName, "--key", proxyName, "-X", "GET", url+File_Name,\
        "-o", tpath + File_Name], stderr=FNULL)
        fin = TFile(tpath + File_Name)
        if fin:
            hclus=fin.Get("DQMData/Run "+ str(run) + "/AlCaReco/Run summary/SiStrip/MechanicalView/TIB/TotalNumberOfCluster__TIB")
            if(hclus):
                nclus = hclus.GetEntries()
                fin.Close()
                return nclus
            else:
                return -3
        else:
            return -2
    else:
       return -1

def getFinalStatusofRun(sFlag):
    if sFlag['hasdb']:
        if sFlag['hasPayload']:
            sFlag['status'] = '\033[1;32mOK\033[1;m'
        else:
            if int(statusFlags['nclusTib']) < 10000:
                sFlag['status'] = '\033[1;32mOK\033[1;m'
            else: 
                sFlag['status'] = '\033[1;31mALERT\033[1;m'  if  statusFlags['eTime'] > 20 else '\033[1;35mWAIT\033[1;m'
    else:
        sFlag['status'] = '\033[1;31mALERT\033[1;m'  if  statusFlags['eTime'] > 12 else '\033[1;35mWAIT\033[1;m'

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: python checkPCLforRun.py <Run No.>....."
        sys.exit(1)
    global FNULL
    FNULL=open(os.devnull, 'w')
    global proxyName
    proxyName = getGridCertificate()
    if not proxyName:
        print "Generate grid certificate!!"
        sys.exit(1)
    print "Using user certificate ", proxyName   
    print "Run\tSiStripBadPCLdb\t\tPayload\t\tInsert Time\t\t\tNClsTIB\tTimeFromRunEnd(Hrs)\tStatus"  
    ##loop over runs
    for run in range(1,len(sys.argv)):
        statusFlags = {'hasdb' : False, 'hasPayload' : False, 'insertTime' : 0, 'nclusTib' : 'NA', 'eTime' : 'NA', 'status' : False}
        checkLogfilesforRun(sys.argv[run],statusFlags)
        statusFlags['eTime'] = getTimeElapsedFromRunEnd(sys.argv[run])
        if  statusFlags['hasdb'] and not statusFlags['hasPayload']:
            statusFlags['nclusTib'] = checkDQMhisto(int(sys.argv[run]))
        getFinalStatusofRun(statusFlags)
        print sys.argv[run], "\t", colorText(statusFlags['hasdb']), "\t\t\t",colorText(statusFlags['hasPayload']),"\t\t",\
        statusFlags['insertTime'],"\t", statusFlags['nclusTib'],"\t",statusFlags['eTime'],\
        '\t\t\t',statusFlags['status']
