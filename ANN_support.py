import numpy as np
import keras as kr
import ROOT as rt

np.seterr(divide='ignore', invalid='ignore')

def discriminants(x_data):
        L1_d = x_data[:,0]
        L4_d = x_data[:,3]
        L1_r = x_data[:,12]
        L2_r = x_data[:,13]
        L3_r = x_data[:,14]
        L4_r = x_data[:,15]
        L2_L1 = L2_r/L1_r
        L3_L2 = L3_r/L2_r
        L4_L3 = L4_r/L3_r
        L4_L1 = L4_r/L1_r
        L2_L1[np.isnan(L2_L1)]=1
        L3_L2[np.isnan(L3_L2)]=1
        L4_L3[np.isnan(L4_L3)]=1
        L4_L1[np.isnan(L4_L1)]=1
        L2_L1[L2_L1>100]=10
        L3_L2[L3_L2>100]=10
        L4_L3[L4_L3>100]=10
        L4_L1[L4_L1>100]=10
        return np.vstack((L4_d-L1_d,L4_L1,L2_L1,L3_L2,L4_L3)).transpose()

def ANN_functional_shape(x_data):
        """takes the 5-cones hits data and puts it into the right input shape for the functional ANN model"""
        x_data_Li_Lj= discriminants(x_data)
        x_data_Li = np.reshape(x_data[:,:20].flatten(),(-1,5,4,1))
        return [x_data_Li, x_data_Li_Lj]

def Predict_from_best_model(nClusters, pt, model):
	"""nClusters is an numpy array of 20 entries corresponding to the number of clusters in each of the 4 layers in each of the 5 cone sizes"""
	nClusters = np.vstack((np.ndarray((0,20)),nClusters))
	pT = np.array([pt/200.])
	X = ANN_functional_shape(nClusters)
        X.append(pT)
	return model.predict(X).item()

def Create_datasets(title, signal_files, bg_files):
	signal_X = np.ndarray((0,23), dtype=np.float64)
	for file_path in signal_files:
		signal_X = np.vstack((signal_X, Extract_Data_From_File(file_path, y=1)))
	bg_X = np.ndarray((0,23), dtype=np.float64)
	for file_path in bg_files:
		bg_X = np.vstack((bg_X, Extract_Data_From_File(file_path, y=0)))

	np.random.shuffle(signal_X)
	np.random.shuffle(bg_X)
	n_s = signal_X.shape[0]
	n_bg = bg_X.shape[0]
	print "There are {} signal jets and {} bg jets.".format(n_s, n_bg)
        if n_bg >1.5*n_s:
                bg_X = bg_X[:int(1.5*n_s),:]
		print "reduced the background to {} jets.".format(int(1.5*n_s))
	elif n_s > 1.5* n_bg:
                signal_X = signal_X[:int(1.5*n_bg),:]
		print "reduced the signal to {} jets.".format(int(1.5*n_bg))

	full_X = np.vstack((signal_X, bg_X))
	np.random.shuffle(full_X)
	m = full_X.shape[0]
	r = int(0.2*m)

	np.save("ANN_data/{}/test_x.npy".format(title), full_X[:r,:20])
        np.save("ANN_data/{}/test_y.npy".format(title), full_X[:r,22].astype(int))
        np.save("ANN_data/{}/test_CSV.npy".format(title), full_X[:r,21])
        np.save("ANN_data/{}/test_pT.npy".format(title), full_X[:r,20])
        np.save("ANN_data/{}/train_x.npy".format(title), full_X[r:,:20])
        np.save("ANN_data/{}/train_y.npy".format(title), full_X[r:,22].astype(int))
        np.save("ANN_data/{}/train_pT.npy".format(title), full_X[r:,20])
	print "Saved data as ANN_data/{}/*_*.npy.".format(title)

def Extract_Data_From_File(input_file, y=1):
	X = np.ndarray((0,23), dtype=np.float64)
	print "opening",input_file
	f1 = rt.TFile.Open(input_file, "READ")
	try:
		tree = f1.Get("demo/tree")
	except ReferenceError:
		print "skipping file"
		return X
	N = tree.GetEntries()
	for i in xrange(N):
                #if i%1000==0: print "scanning event nr.{}".format(i)
                tree.GetEntry(i)
                for j in xrange(tree.nJets):
			if y==1 and tree.jet_MC_bTag[j] != 1: continue
			if y==0 and tree.jet_MC_bTag[j] != 0: continue
                        _X = np.array([tree.nClusters_L1004[j],tree.nClusters_L1006[j],tree.nClusters_L1008[j],tree.nClusters_L1010[j],tree.nClusters_L1016[j],tree.nClusters_L2004[j],tree.nClusters_L2006[j],tree.nClusters_L2008[j],tree.nClusters_L2010[j],tree.nClusters_L2016[j],tree.nClusters_L3004[j],tree.nClusters_L3006[j],tree.nClusters_L3008[j],tree.nClusters_L3010[j],tree.nClusters_L3016[j],tree.nClusters_L4004[j],tree.nClusters_L4006[j],tree.nClusters_L4008[j],tree.nClusters_L4010[j],tree.nClusters_L4016[j], tree.jet_pt[j], tree.jet_bTag[j], y], dtype=np.float64)
			X = np.vstack((X, _X))
	return X

def reweight(data_path, train_pt, train_y, bin_size, scale_correction=1):
	maximum = 5000

	import json
	with open("{}/weights.json".format(data_path), "r") as fp:
		weight_function = json.load(fp)

        weights = np.zeros(len(train_pt))
        for n,entry in enumerate(train_y):
                if entry == 0:
                        weights[n] = 1
                else:
                        if train_pt[n] > maximum:
                                weights[n] = 0
                        else:
                                weights[n] = scale_correction*weight_function[str((int(train_pt[n])/bin_size)*bin_size)]
        return weights

def find_weight_function(data_path, Numerator_pT, Denominator_pT, binsize):

	maximum = 5000
        bins = range(0,maximum+1,binsize)
        nbins = len(bins)-1
        import array
        bins_ = array.array('d',bins)
        Numerator_hist = rt.TH1D("Numerator_hist","Numerator_hist",nbins,bins_)
        Denominator_hist = rt.TH1D("Denominator_hist","Denominator_hist",nbins,bins_)

        for pt in Numerator_pT:
                Numerator_hist.Fill(pt)
        for pt in Denominator_pT:
                Denominator_hist.Fill(pt)

        norm = 1

        Numerator_hist.Scale(norm/Numerator_hist.Integral())
        Denominator_hist.Scale(norm/Denominator_hist.Integral())
        Ratio_hist = Numerator_hist.Clone()
        Ratio_hist.SetName("Ratio_hist")
        Ratio_hist.SetTitle("Ratio_hist")
        Ratio_hist.Divide(Denominator_hist)

        ratio_dict = {}
        for bin_nr in range(nbins):
                ratio_dict[bins[bin_nr]] = Ratio_hist.GetBinContent(bin_nr+1)

        import json
	with open("{}/weights.json".format(data_path), "w") as fp:
		json.dump(ratio_dict,fp)

        #test
        Denominator_test_hist = rt.TH1D("Denominator_test_hist","Denominator_test_hist",nbins,bins_)
        for pt in Denominator_pT:
                if pt>maximum: continue
                Denominator_test_hist.Fill(pt,ratio_dict[(int(pt)/binsize)*binsize])

        Denominator_test_hist.Scale(1.1/Denominator_test_hist.Integral())

        canvas = rt.TCanvas('canvas','canvas',600,600)
        legend = rt.TLegend(0.78,0.9,0.54,0.79)
        Numerator_hist.SetLineColor(2)
        Denominator_hist.SetLineColor(3)
        Denominator_test_hist.SetLineColor(4)
        legend.AddEntry(Numerator_hist, 'num')
        legend.AddEntry(Denominator_hist, 'denom')
        legend.AddEntry(Denominator_test_hist, 'denom_test')

        Denominator_test_hist.Draw()
        Denominator_hist.Draw("SAME")
        Numerator_hist.Draw("SAME")
        #Denominator_test_hist.Draw("SAME")
        legend.Draw()
	f1 = rt.TFile("Plots/test_hist.root", "RECREATE")
        canvas.Write()
	f1.Close()

if __name__ == "__main__":

	
	signal_path = "/eos/user/m/msommerh/HitAnalyzer4_collected/{}/M{}/flatTuple_{}.root"
	bg_path = "/eos/user/m/msommerh/HitAnalyzer4_collected/QCD/{}/flatTuple_{}.root"
	
	#M0_list = ['1000', '1200', '1400', '1600', '1800', '2000', '2500', '3000', '3500', '4000', '4500', '5000', '5500', '6000']
	M0_list = ['2000', '2500', '3000', '3500', '4000', '4500', '5000', '5500', '6000']
	bin_list = [('170','300'), ('300','470'), ('470','600'), ('600','800'), ('800','1000'), ('1000','1400'), ('1400','1800'), ('1800','2400'), ('2400','3200'), ('3200', 'Inf')]

	signal_files = []
	for M0 in M0_list:
		for i in range(1,41):
			signal_files.append(signal_path.format('2017',M0,i))
			signal_files.append(signal_path.format('2018',M0,i))
	bg_files = []
	for bin_ in bin_list:
		for i in range(1,61):
			bg_files.append(bg_path.format(bin_[0]+"to"+bin_[1],i))
	Create_datasets("lessFilter2", signal_files, bg_files)	
		
	train_pT = np.load("ANN_data/lessFilter2/train_pT.npy")
	train_y = np.load("ANN_data/lessFilter2/train_y.npy")

	bg_pT = train_pT[train_y==0]
	signal_pT = train_pT[train_y==1]
	print "bg_pT =",bg_pT[:20]
	print "signal_pT =",signal_pT[:20]	
	print "train_pT =",train_pT[:40]
	print "train_y =",train_y[:40]

	find_weight_function("ANN_data/lessFilter2", bg_pT, signal_pT, 10)
	
