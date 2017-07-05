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

#pattern = '/ZeroBias*/Run2016*SiStripCalMinBias-18Apr2017*-v*/ALCARECO'
pattern = '/StreamExpress/Run2016*SiStripCalMinBias-Express*-v*/ALCARECO' 
ds = commands.getstatusoutput("dasgoclient -query='dataset dataset="+pattern+"'")[1].split("\n")
ds.sort()
for d in ds:
    events = commands.getstatusoutput("das_client.py --limit=0 --query='summary dataset="+d+" | grep summary.nevents'")[1]
    print d,events


