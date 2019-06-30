#! /usr/bin/env python

import os, multiprocessing, math, sys
import ROOT as rt

def submitJobs(jobname, jobflavour, nEpochs, title, data_path):
    path = os.getcwd()
    os.makedirs("/afs/cern.ch/work/m/msommerh/public/ANN_models/"+title)
    workdir = "tmp"+jobname
    os.makedirs(workdir)
    os.chdir(workdir)
       
    #write executable file for submission
    with open('job_{}.sh'.format(1), 'w') as fout:
        fout.write("#!/bin/sh\n")
        fout.write("echo\n")
        fout.write("echo\n")
        fout.write("echo 'START---------------'\n")
        fout.write("echo 'WORKDIR ' ${PWD}\n")

        fout.write("cd "+str(path)+"\n")
        fout.write("export SCRAM_ARCH=slc6_amd64_gcc630\n" )
        fout.write("if [ -r CMSSW_10_2_5_patch1/src ] ; then\n")
        fout.write("    echo 'release CMSSW_10_2_5_patch1 already exists'\n")
        fout.write("else\n")
        fout.write("    scram p CMSSW CMSSW_10_2_5_patch1\n")
        fout.write("fi\n")
        fout.write("cd CMSSW_10_2_5_patch1/src\n")
        fout.write("eval `scram runtime -sh`\n")
        fout.write("cd -\n" )
        fout.write("echo 'cmssw release = ' $CMSSW_BASE\n")
        fout.write("python ANN_training.py {} {} {}\n".format(title, data_path, nEpochs))
        fout.write("echo 'STOP---------------'\n")
        fout.write("echo\n")
        fout.write("echo\n")
    
    #submit job
    os.system("chmod 755 job_{}.sh".format(1) )
    os.system("mv job_*.sh "+jobname+"_"+str(1)+".sh")
    makeSubmitFileCondor(jobname+"_"+str(1)+".sh", jobname, jobflavour)
    os.system("condor_submit submit.sub")
    print "job {} nr {} submitted".format(jobname, 1)
    os.chdir("../..")

def makeSubmitFileCondor(exe, jobname, jobflavour): 
    print "make options file for condor job submission"
    submitfile = open("submit.sub", "w")
    submitfile.write("executable  = "+exe+"\n")
    submitfile.write("arguments             = $(ClusterID) $(ProcId)\n")
    submitfile.write("output                = "+jobname+".$(ClusterId).$(ProcId).out\n")
    submitfile.write("error                 = "+jobname+".$(ClusterId).$(ProcId).err\n")
    submitfile.write("log                   = "+jobname+".$(ClusterId).log\n")
    submitfile.write('+JobFlavour           = "'+jobflavour+'"\n')
    submitfile.write("queue")
    submitfile.close()
        

if __name__ == "__main__":
 
  jobname = "ANN_training_updated_ptcut_ZPrime-matching_quitelong"

  ##choose priority
  #jobflavour = 'espresso' #max 30min
  #jobflavour = 'microcentury' #max 1h
  #jobflavour = 'longlunch' #max 2h
  #jobflavour = 'workday' #max 8h
  #jobflavour = 'tomorrow' #max 1d
  jobflavour = 'testmatch' #max 3d

  nEpochs = 200
  title = "updated_ptcut_ZPrime-matching_quitelong"
  data_path = "/afs/cern.ch/work/m/msommerh/public/ANN_data/updated_ptcut_ZPrime-matching"
  submitJobs(jobname, jobflavour, nEpochs, title, data_path)
 
  print
  print "your jobs:"
  os.system("condor_q")
  userName=os.environ['USER']
  
  print
  print 'Done submitting jobs!'
  print
  
