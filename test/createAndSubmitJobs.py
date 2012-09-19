#!/usr/bin/env python

import os, sys, optparse, subprocess, string, re


def make_filenamelist(inputDir):
    if not os.path.isdir(inputDir):
        print ('%s is not a directory'%(inputDir))
        sys.exit(1)

    filenamelist = []
    for filename in os.listdir(inputDir):
        if not os.path.isfile(os.path.join(inputDir,filename)):
            continue
        filenamelist.append(filename)

    return filenamelist


def process_input_dir(inputDir, match, filelist):
    inputDir = inputDir.rstrip('/')+'/'
    filenamelist = make_filenamelist(inputDir)

    path = inputDir;
    name = ''
    jobdict = {}

    for filename in filenamelist:
        if( not re.search('.root$', filename) ):
            continue
        if ( match!=None and not re.search(match, filename) ):
            continue
        m1 = re.search('_\d+_\d+_\w+.root', filename)
        if name=='':
            name = re.split('_\d+_\d+_\w+.root', filename)[0]
        jobstring = filename[m1.start():].lstrip('_').replace('.root','').split('_')
        job = int(jobstring[0])
        if job not in jobdict.keys():
            jobdict[job] = []
            jobdict[job].append([int(jobstring[1])])
            jobdict[job].append([jobstring[2]])
        else:
            jobdict[job][0].append(int(jobstring[1]))
            jobdict[job][1].append(jobstring[2])

    jobs = jobdict.keys()
    if( len(jobs)==0 ):
        print 'No matching .root files found'
        sys.exit()

    jobs.sort()
    for job in jobs:
        maxsub = max(jobdict[job][0])
        filename = (path+name+'_%i_%i_%s.root')%(job, maxsub, jobdict[job][1][jobdict[job][0].index(maxsub)])
        filelist.append(filename)

    return


jdl_template = """universe = vanilla
Notification = never
Executable = DUMMY_OUTPUTDIR/CMSSW.sh
Output = DUMMY_OUTPUTDIR/logs/condor_$(Cluster)_$(Process).stdout
Error = DUMMY_OUTPUTDIR/logs/condor_$(Cluster)_$(Process).stderr
Log = DUMMY_OUTPUTDIR/logs/condor_$(Cluster)_$(Process).log
Arguments = $(Process) DUMMY_OUTPUTDIR
Queue DUMMY_NJOBS
"""


bash_template = """#!/bin/bash

CONDOR_PROCESS=$1
OUTPUTDIR=$2

START_TIME=`/bin/date`
echo "Started at $START_TIME"
echo ""

export SCRAM_ARCH=slc5_amd64_gcc462
source /cms/base/cmssoft/cmsset_default.sh
cd $OUTPUTDIR
eval `scramv1 runtime -sh`

LOG="${OUTPUTDIR}/logs/CMSSW_Job_${CONDOR_PROCESS}.log"

echo "Running CMSSW job $CONDOR_PROCESS"
cmsRun -j ${OUTPUTDIR}/output/CMSSW_Job_${CONDOR_PROCESS}_fjr.xml -p ${OUTPUTDIR}/cfgs/CMSSW_Job_${CONDOR_PROCESS}_cfg.py > $LOG 2>&1
exitcode=$?

echo ""
END_TIME=`/bin/date`
echo "Finished at $END_TIME"
exit $exitcode
"""


usage = """Usage: %prog [-m MATCH] -i INPUTDIR -c CMSSW_CFG -o OUTPUTDIR -n NJOBS [--create_only]\n
Example: ./createAndSubmitJobs.py -i /cms/ferencek/store/ferencek/WW_TuneZ2star_8TeV_pythia6_tauola/Summer12-PU_S7_START52_V9-v1_PATTuple/c356f3ff60b3ed8401082f13b541664f -c rutgersjetanalyzer_template_cfg.py -o WW_test -n 50"""


def main():

    parser = optparse.OptionParser(usage=usage)

    parser.add_option('-m', '--match', metavar='MATCH', action='store', help='Only files containing the MATCH string in their names will be considered (This parameter is optional)')
    parser.add_option('-i', '--inputdir', metavar='INPUTDIR', dest='inputdir', action='store', help='Input directory containing .root files')
    parser.add_option('-c', '--cmssw_cfg', metavar="CMSSW_CFG", dest='cmssw_cfg', action='store', help='CMSSW configuration file template')
    parser.add_option('-o', '--outputdir', metavar='OUTPUTDIR', dest='outputdir', action='store', help='Main output directory for Condor jobs')
    parser.add_option('-n', '--njobs', metavar='NJOBS', action='store', dest='njobs', help='Number of jobs')
    parser.add_option('--create_only', action="store_true", dest="create_only", default=False, help="Create the necessary configuration files but skip the job submission (This parameter is optional)")

    (options, args) = parser.parse_args(args=None)

    if ( options.inputdir==None or options.cmssw_cfg==None or options.outputdir==None or options.njobs==None ):
        print ('\nOptions -i, -c, -o, and -n  are required\n')
        parser.print_help()
        sys.exit()

    filelist = []

    process_input_dir(options.inputdir, options.match, filelist)

    # read CMSSW configuration file template
    cmssw_cfg = open(options.cmssw_cfg,'r')
    cmssw_cfg_content = cmssw_cfg.read()

    if( not re.search('DUMMY_INPUTFILES',cmssw_cfg_content) or not re.search('DUMMY_OUTPUTFILE',cmssw_cfg_content) ):
        print "CMSSW configuration file needs to contain DUMMY_INPUTFILES and DUMMY_OUTPUTFILE variables"
        print "Aborting"
        sys.exit()

    ################################################
    pwd = os.getcwd()
    # redefine main output directory as an absolute path (if not defined in such form already)
    outputmain = options.outputdir.rstrip('/')
    if not re.search('^/', outputmain):
	outputmain = pwd + '/' + outputmain

    os.mkdir(outputmain)
    os.mkdir(os.path.join(outputmain,'cfgs'))
    os.mkdir(os.path.join(outputmain,'logs'))
    os.mkdir(os.path.join(outputmain,'output'))
    #################################################
    numfiles = len(filelist)
    ijobmax=int(options.njobs)
    if ijobmax > numfiles:
	ijobmax=numfiles
	print 'Number of jobs requested exceeds the total number of input files.\nThe number of jobs set to ' + str(ijobmax) + ' to match the number of input files'
    filesperjob = int(numfiles/ijobmax)
    if numfiles%ijobmax!=0:
	filesperjob = filesperjob + 1
	ijobmax = int(numfiles/filesperjob)
	if numfiles%filesperjob!=0:
	    ijobmax = ijobmax + 1
        print 'Could not create the exact number of jobs requested.\nFor optimal job splitting, the number of jobs set to '+ str(ijobmax)
    #################################################
    # create jdl file
    jdl_file = open(os.path.join(outputmain,'CMSSW.jdl'),'w')
    jdl_content = re.sub('DUMMY_OUTPUTDIR',outputmain,jdl_template)
    jdl_content = re.sub('DUMMY_NJOBS',str(ijobmax),jdl_content)
    jdl_file.write(jdl_content)
    jdl_file.close()
    # create Bash script
    bash_script = open(os.path.join(outputmain,'CMSSW.sh'),'w')
    bash_script.write(bash_template)
    bash_script.close()
    os.system('chmod +x '+ os.path.join(outputmain,'CMSSW.sh'))
    ################################################
    ifile = 0
    for ijob in range(ijobmax):
	# prepare the list of input files
	inputfiles = '\'file:' + filelist[ifile] + '\''
	for i in range(filesperjob-1):
	    if ifile>(numfiles-2):
	        break
	    ifile = ifile + 1
	    inputfiles += (',\n\'file:' + filelist[ifile] + '\'')
	ifile = ifile + 1

	# create CMSSW configuration file
	job_cfg = open(os.path.join(outputmain,'cfgs','CMSSW_Job_' + str(ijob) + '_cfg.py'),'w')
	job_cfg_content = re.sub('DUMMY_INPUTFILES',inputfiles,cmssw_cfg_content)
	job_cfg_content = re.sub('DUMMY_OUTPUTFILE','\'' + os.path.join(outputmain,'output','CMSSW_Job_' + str(ijob) + '.root\''),job_cfg_content)
	job_cfg.write(job_cfg_content)
	job_cfg.close()

    if not options.create_only:
      # submit Condor jobs
      os.system('condor_submit ' + os.path.join(outputmain,'CMSSW.jdl'))


if __name__ == '__main__':
    main()
