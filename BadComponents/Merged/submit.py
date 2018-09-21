import os
import sys
import calendar
import optparse
import importlib
import sqlalchemy
import subprocess
import CondCore.Utilities.conddblib as conddb

##############################################
def execme(command, dryrun=False):
##############################################
    '''Wrapper for executing commands.
    '''
    if dryrun:
        print command
    else:
        print " * Executing: %s ..." % (command)
        os.system(command)
        print " * Executed!"

##############################################
def main():
##############################################

    defaultFile='DQM_V0001_R000291659__Cosmics__Commissioning2017-PromptReco-v1__DQMIO.root'
    defaultGT='90X_dataRun2_Prompt_v3'
    defaultRun = 291659
    
    parser = optparse.OptionParser(usage = 'Usage: %prog [options] <file> [<file> ...]\n')
    
    parser.add_option('-f', '--inputDQMFile',
                      dest = 'inputDQMFile',
                      default = defaultFile,
                      help = 'DQM File for FED Errors',
                      )

    parser.add_option('-G', '--inputGT',
                      dest = 'inputGT',
                      default = defaultGT,
                      help = 'Global Tag to get conditions',
                      )

    parser.add_option('-r', '--inputRun',
                      dest = 'inputRun',
                      default = defaultRun,
                      help = 'run to be used',
                      )

    (options, arguments) = parser.parse_args()

    print "Input configuration"
    print "globalTag: ",options.inputGT
    print "runNumber: ",options.inputRun
    print "dqmFile:   ",options.inputDQMFile

    con = conddb.connect(url = conddb.make_url())
    session = con.session()
    RunInfo = session.get_dbtype(conddb.RunInfo)
    
    bestRun = session.query(RunInfo.run_number,RunInfo.start_time, RunInfo.end_time).filter(RunInfo.run_number >= options.inputRun).first()
    if bestRun is None:
        raise Exception("Run %s can't be matched with an existing run in the database." %options.runNumber)
    
    start= bestRun[1]
    stop = bestRun[2]
    
    bestRunStartTime = calendar.timegm( bestRun[1].utctimetuple() ) << 32
    bestRunStopTime  = calendar.timegm( bestRun[2].utctimetuple() ) << 32
    
    print "run start time:",start,"(",bestRunStartTime,")"
    print "run stop time: ",stop,"(",bestRunStopTime,")"
    
    command='cmsRun mergeBadChannel_Template_cfg.py globalTag='+options.inputGT+' runNumber='+str(options.inputRun)+' dqmFile='+options.inputDQMFile+' runStartTime='+str(bestRunStopTime)
    
    execme(command)

if __name__ == "__main__":        
    main()
