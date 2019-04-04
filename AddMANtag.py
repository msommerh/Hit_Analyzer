import ROOT as rt
import numpy as np
from array import array
import sys

import keras as kr
from ANN_support import Predict_from_best_model

if len(sys.argv)>3:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        model_path = sys.argv[3]
else:
        print "using default settings:"
        print "input_file = flatTuple_MANtest.root"
        print "output_file = AddMANtag2_test.root"
        print "model_path = ANN/submitted_models_fullrun1_PU_pt_both.h5"
        input_file = "flatTuple_MANtest.root"
        output_file = "AddMANtag2_test.root"
        model_path = "ANN/submitted_models_fullrun1_PU_pt_both.h5"


model = kr.models.load_model(model_path)
np.seterr(divide='ignore', invalid='ignore')


if output_file != input_file:
	f1 = rt.TFile.Open(input_file, "READ")
else:
	f1 = rt.TFile.Open(input_file, "UPDATE")
tree = f1.Get("demo/tree")

maxnJ =150 

_MANtag = np.zeros(maxnJ, dtype=np.double)

jet_branch = tree.GetBranch("nJets")
MANtag_branch = tree.Branch("MANtag",_MANtag, "MANtag[nJets]/D")

N = tree.GetEntries()

for i in xrange(N):
	if i%100==0: print "scanning event nr.{}".format(i)
	tree.GetEntry(i)
	for j in xrange(tree.nJets):
                NC = np.array([tree.nClusters_L1004[j],tree.nClusters_L1006[j],tree.nClusters_L1008[j],tree.nClusters_L1010[j],tree.nClusters_L1016[j],tree.nClusters_L2004[j],tree.nClusters_L2006[j],tree.nClusters_L2008[j],tree.nClusters_L2010[j],tree.nClusters_L2016[j],tree.nClusters_L3004[j],tree.nClusters_L3006[j],tree.nClusters_L3008[j],tree.nClusters_L3010[j],tree.nClusters_L3016[j],tree.nClusters_L4004[j],tree.nClusters_L4006[j],tree.nClusters_L4008[j],tree.nClusters_L4010[j],tree.nClusters_L4016[j]], dtype=np.int32)
                _MANtag[j] = Predict_from_best_model(NC, tree.jet_pt[j], model)

	MANtag_branch.Fill()

if output_file != input_file:
	f2 = rt.TFile(output_file,"RECREATE")
	tree2 = tree.CloneTree()
	tree2.Write()
	f2.Close()
	print "MANtag variable added and saved in:",output_file
else:
	f1.cd("demo")
	tree.Write()
	rt.gDirectory.Delete("tree;1")
	print "MANtag variable added directly to the file:", output_file

f1.Close()
