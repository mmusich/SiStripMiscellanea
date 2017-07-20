#/env/python 

#!/usr/bin/env python

import os,sys
import getopt
import commands
import time
import ROOT
import urllib
import string
from optparse import OptionParser

pattern = '/StreamExpress/Run2017B-PromptCalibProdSiStripGainsAAG-Express-v2/ALCAPROMPT ' 
runs = [298653,298678,298674,298996,298997,298998,299000,299042,299061,299062,299065,299096,299149,299183,299185,299326,299064,299178,299180,299067,299184,299316,299317,299318,299324,299325]

filesToHarvest=[]

fout=open("InputFiles_cff.py",'w+b')

for run in runs:
    ds = commands.getstatusoutput("dasgoclient -query='file dataset="+pattern+" run="+str(run)+"'")[1].split("\n")
    filesToHarvest+=ds

for m_file in filesToHarvest:
    fout.write("'"+m_file+"',\n")

print len(filesToHarvest)


