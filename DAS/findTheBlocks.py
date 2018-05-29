import os,sys
import getopt
import commands
import time
import ROOT
import urllib
import string
from optparse import OptionParser

runs=[314755,314756]

blocksToKeep=[]

for run in runs:
    for i in xrange(1,5):
        ds="/ZeroBias%d/Commissioning2018-SiStripCalMinBias-PromptReco-v1/ALCARECO" % i
        blocks = commands.getstatusoutput("dasgoclient -query='block dataset="+ds+" run="+str(run)+"'")[1].split("\n")
        for block in blocks:
            blocksToKeep.append(block)
        #print run,ds,blocks

output=[]
for x in blocksToKeep:
    if x not in output:
        output.append(x)
    else:
        pass
        #print x,"already stored" 

total_size=0
for block in output:
    block_size = commands.getstatusoutput("dasgoclient -query='summary block="+block+"| grep summary.file_size'")[1]
    total_size+=float(block_size)
    print block+","

print "total size:",total_size/(1000*1000*1000*1000),"TB"
