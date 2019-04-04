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
	print "input_file = flatTuple.root"
        print "output_file = AddMANtag_test3.root"
        print "model_path = ANN/submitted_models_fullrun1_PU_pt_both.h5" 
	input_file = "flatTuple.root"
        output_file = "AddMANtag_test3.root"
        model_path = "ANN/submitted_models_fullrun1_PU_pt_both.h5" 

	
model = kr.models.load_model(model_path)
np.seterr(divide='ignore', invalid='ignore')

#reading

f1 = rt.TFile.Open(input_file)
tree1 = f1.Get("demo/tree")
N1 = tree1.GetEntries()
N = N1

nJets = []
jet_pt = []
jet_eta = []
jet_phi = []
jet_mass = []
jet_bTag = []
nPV = []
PV_x = []
PV_y = []
PV_z = []
nClusters = []
MANtag = []

for i in xrange(N):
	if i%100==0: print "scanning event nr.{}".format(i)
	tree1.GetEntry(i)
	_nJets = tree1.nJets
	_jet_pt = []
	_jet_eta = []
	_jet_phi = []
	_jet_mass = []
	_jet_bTag = []
	_nPV = tree1.nPV
	_PV_x = []
	_PV_y = []
	_PV_z = []
	_nClusters = np.ndarray((0,20), dtype=np.int32)
	_MANtag = []

	for j in xrange(tree1.nJets):
		_jet_pt 	.append(tree1.jet_pt  [j]) 	
                _jet_eta 	.append(tree1.jet_eta [j])
                _jet_phi 	.append(tree1.jet_phi [j])
                _jet_mass 	.append(tree1.jet_mass[j])
                _jet_bTag 	.append(tree1.jet_bTag[j])
		NC = np.array([tree1.nClusters_L1004[j],tree1.nClusters_L1006[j],tree1.nClusters_L1008[j],tree1.nClusters_L1010[j],tree1.nClusters_L1016[j],tree1.nClusters_L2004[j],tree1.nClusters_L2006[j],tree1.nClusters_L2008[j],tree1.nClusters_L2010[j],tree1.nClusters_L2016[j],tree1.nClusters_L3004[j],tree1.nClusters_L3006[j],tree1.nClusters_L3008[j],tree1.nClusters_L3010[j],tree1.nClusters_L3016[j],tree1.nClusters_L4004[j],tree1.nClusters_L4006[j],tree1.nClusters_L4008[j],tree1.nClusters_L4010[j],tree1.nClusters_L4016[j]], dtype=np.int32)
		_nClusters = np.vstack((_nClusters,NC))
		_MANtag.append(Predict_from_best_model(NC, tree1.jet_pt[j], model))

	for j in xrange(tree1.nPV):
		_PV_x.append(tree1.PV_x[j]) 
                _PV_y.append(tree1.PV_y[j])
                _PV_z.append(tree1.PV_z[j])

	nJets   	.append(_nJets   	)   
	jet_pt  	.append(_jet_pt  	)
	jet_eta 	.append(_jet_eta 	) 
	jet_phi 	.append(_jet_phi 	) 
	jet_mass	.append(_jet_mass	)  
	jet_bTag	.append(_jet_bTag	)
	nPV		.append(_nPV		)
	PV_x		.append(_PV_x		) 
	PV_y		.append(_PV_y		)
	PV_z		.append(_PV_z		)
	nClusters	.append(_nClusters	)
	MANtag		.append(_MANtag		)

f1.Close()


#writing

f2 = rt.TFile(output_file,"RECREATE")
tree2 = rt.TTree("tree2", "tree2")

maxnJ = max(nJets)
maxnPV = max(nPV)

_nJets 		 = np.zeros(1, dtype=np.int32)
_jet_pt 	 = np.zeros(maxnJ, dtype=np.double) 
_jet_eta 	 = np.zeros(maxnJ, dtype=np.double)
_jet_phi 	 = np.zeros(maxnJ, dtype=np.double)
_jet_mass 	 = np.zeros(maxnJ, dtype=np.double)
_jet_bTag 	 = np.zeros(maxnJ, dtype=np.double)
_nClusters_L1004 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L1006 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L1008 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L1010 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L1016 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L2004 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L2006 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L2008 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L2010 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L2016 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L3004 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L3006 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L3008 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L3010 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L3016 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L4004 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L4006 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L4008 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L4010 = np.zeros(maxnJ, dtype=np.int32)
_nClusters_L4016 = np.zeros(maxnJ, dtype=np.int32)
_MANtag		 = np.zeros(maxnJ, dtype=np.double)
_nPV 		 = np.zeros(1, dtype=np.int32)
_PV_x	         = np.zeros(maxnPV, dtype=np.double)
_PV_y	         = np.zeros(maxnPV, dtype=np.double)
_PV_z	         = np.zeros(maxnPV, dtype=np.double)

tree2.Branch("nJets",_nJets	       , "nJets/I") 
tree2.Branch("jet_pt",_jet_pt 	       , "jet_pt[nJets]/D")   
tree2.Branch("jet_eta",_jet_eta 	       , "jet_eta [nJets]/D")   
tree2.Branch("jet_phi",_jet_phi 	       , "jet_phi [nJets]/D")   
tree2.Branch("jet_mass",_jet_mass	       , "jet_mass[nJets]/D")   
tree2.Branch("jet_bTag",_jet_bTag	       , "jet_bTag[nJets]/D")   
tree2.Branch("nClusters_L1004",_nClusters_L1004 , "nClusters_L1004[nJets]/I")
tree2.Branch("nClusters_L1006",_nClusters_L1006 , "nClusters_L1006[nJets]/I")
tree2.Branch("nClusters_L1008",_nClusters_L1008 , "nClusters_L1008[nJets]/I")
tree2.Branch("nClusters_L1010",_nClusters_L1010 , "nClusters_L1010[nJets]/I")
tree2.Branch("nClusters_L1016",_nClusters_L1016 , "nClusters_L1016[nJets]/I")
tree2.Branch("nClusters_L2004",_nClusters_L2004 , "nClusters_L2004[nJets]/I")
tree2.Branch("nClusters_L2006",_nClusters_L2006 , "nClusters_L2006[nJets]/I")
tree2.Branch("nClusters_L2008",_nClusters_L2008 , "nClusters_L2008[nJets]/I")
tree2.Branch("nClusters_L2010",_nClusters_L2010 , "nClusters_L2010[nJets]/I")
tree2.Branch("nClusters_L2016",_nClusters_L2016 , "nClusters_L2016[nJets]/I")
tree2.Branch("nClusters_L3004",_nClusters_L3004 , "nClusters_L3004[nJets]/I")
tree2.Branch("nClusters_L3006",_nClusters_L3006 , "nClusters_L3006[nJets]/I")
tree2.Branch("nClusters_L3008",_nClusters_L3008 , "nClusters_L3008[nJets]/I")
tree2.Branch("nClusters_L3010",_nClusters_L3010 , "nClusters_L3010[nJets]/I")
tree2.Branch("nClusters_L3016",_nClusters_L3016 , "nClusters_L3016[nJets]/I")
tree2.Branch("nClusters_L4004",_nClusters_L4004 , "nClusters_L4004[nJets]/I")
tree2.Branch("nClusters_L4006",_nClusters_L4006 , "nClusters_L4006[nJets]/I")
tree2.Branch("nClusters_L4008",_nClusters_L4008 , "nClusters_L4008[nJets]/I")
tree2.Branch("nClusters_L4010",_nClusters_L4010 , "nClusters_L4010[nJets]/I")
tree2.Branch("nClusters_L4016",_nClusters_L4016 , "nClusters_L4016[nJets]/I")
tree2.Branch("MANtag"	      ,_MANtag        	, "MANtag[nJets]/D")
tree2.Branch("nPV" ,_nPV,  "nPV/I")
tree2.Branch("PV_x",_PV_x, "PV_x[nPV]/D")
tree2.Branch("PV_y",_PV_y, "PV_y[nPV]/D")
tree2.Branch("PV_z",_PV_z, "PV_z[nPV]/D")

for i in xrange(N):
	_nJets[0] = nJets[i]
	_nPV[0]  = nPV[i]
	
	for j in range(_nJets[0]):
		_jet_pt[j] = jet_pt  [i][j] 
                _jet_eta[j] = jet_eta [i][j] 
                _jet_phi[j] = jet_phi [i][j] 
                _jet_mass[j] = jet_mass[i][j] 
                _jet_bTag[j] = jet_bTag[i][j]
		_MANtag[j]   = MANtag[i][j] 

	        _nClusters_L1004[j] = nClusters[i][j][0]  
                _nClusters_L1006[j] = nClusters[i][j][1]
                _nClusters_L1008[j] = nClusters[i][j][2]
                _nClusters_L1010[j] = nClusters[i][j][3]
                _nClusters_L1016[j] = nClusters[i][j][4]
                _nClusters_L2004[j] = nClusters[i][j][5]
                _nClusters_L2006[j] = nClusters[i][j][6]
                _nClusters_L2008[j] = nClusters[i][j][7]
                _nClusters_L2010[j] = nClusters[i][j][8]
                _nClusters_L2016[j] = nClusters[i][j][9]
                _nClusters_L3004[j] = nClusters[i][j][10]
                _nClusters_L3006[j] = nClusters[i][j][11]
                _nClusters_L3008[j] = nClusters[i][j][12]
                _nClusters_L3010[j] = nClusters[i][j][13]
                _nClusters_L3016[j] = nClusters[i][j][14]
                _nClusters_L4004[j] = nClusters[i][j][15]
                _nClusters_L4006[j] = nClusters[i][j][16]
                _nClusters_L4008[j] = nClusters[i][j][17]
                _nClusters_L4010[j] = nClusters[i][j][18]
                _nClusters_L4016[j] = nClusters[i][j][19]
	
	for j in range(_nPV[0]):
		_PV_x[j] = PV_x[i][j] 
        	_PV_y[j] = PV_y[i][j]
        	_PV_z[j] = PV_z[i][j]

	tree2.Fill()

f2.Write()
f2.Close()

