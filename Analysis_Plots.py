import ROOT as rt
import numpy as np
import matplotlib.pyplot as plt

#import keras as kr

from ANN_support import ANN_functional_shape

def ANN_bin_selection(pT,bins):
        """numpy array version of the bin_selection() function: takes an array of all pT values and the
pT-bins and returns an array of same length as pT with the corresponing bins index at each entry. pT val
ues outside the bins are labeled with -100"""
        bin_numbers = np.zeros(len(pT))
        for n in range(len(bins)-1):
                bin_numbers += (n+100)*(pT>bins[n])*(pT<bins[n+1])
        bin_numbers -=100
        return bin_numbers.astype(int)

def MANtag_pT_CSV_from_file(input_file, MC=True, BG=False):
	MANtag = []
	pT = []
	CSV = []
	print "opening file:", input_file
	f1 = rt.TFile.Open(input_file, "READ")
	try:
		tree = f1.Get("demo/tree")
	except ReferenceError:
		print "skipping file"
		return MANtag, pT, CSV
	N = tree.GetEntries()
	print N, "events"
	for i in xrange(N):
		tree.GetEntry(i)
		for j in xrange(tree.nJets):
			if MC:
				if BG and tree.jet_MC_bTag[j] != 0: continue
				if not BG and tree.jet_MC_bTag[j] != 1: continue
			MANtag.append(tree.MANtag[j])
			pT.append(tree.jet_pt[j])
			CSV.append(tree.jet_bTag[j]) ##CSVv2
			#CSV.append(tree.jet_deepCSV_probb[j]) ##deepCSV
	return MANtag, pT, CSV

def MANtag_pT_CSV_from_single_file(input_file, MC=True, BG=False):
	MANtag, pT, CSV = MANtag_pT_CSV_from_file(input_file, MC=MC, BG=BG)
	return np.array(MANtag), np.array(pT), np.array(CSV)

def MANtag_pT_CSV_from_multiple_files(file_list, MC=True, BG=False):
	MANtag, pT, CSV = [], [], []
	for f in file_list:
		_MANtag, _pT, _CSV = MANtag_pT_CSV_from_file(f, MC=MC, BG=BG)
		MANtag += _MANtag
		pT += _pT
		CSV += _CSV
	return np.array(MANtag), np.array(pT), np.array(CSV)

def ANN_Make_Binned_ROC_histograms(title, MANtag, pT, CSV, bins):
        """makes binned ROC histograms for an ANN. takes as input a keras model, the necessary ANN data, pT, CSV and the desired pT-binning"""
	N = len(MANtag)
	assert len(pT) == N and len(CSV) == N

        nbins = 60
        ANN_hist_list = []
        CSV_hist_list = []
        for bin_ in range(len(bins)-1):
                ANN_hist_list.append(rt.TH1D("ANN_"+str(bins[bin_])+"_"+str(bins[bin_+1]),"ANN_"+str(bins[bin_])+"_"+str(bins[bin_+1]),nbins,0,1))
                CSV_hist_list.append(rt.TH1D("CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1]),"CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1]),nbins,0,1))

        bin_numbers = ANN_bin_selection(pT,bins)

        for n in xrange(N):
                if bin_numbers[n] == -100: continue
                ANN_hist_list[int(bin_numbers[n])].Fill(MANtag[n])
                CSV_hist_list[int(bin_numbers[n])].Fill(CSV[n])

        tfile = rt.TFile("ROC_histograms/{}_ROC_histograms.root".format(title),"recreate")
        for hist in ANN_hist_list:
                hist.Write()
        for hist in CSV_hist_list:
                hist.Write()
        print "saved histograms in ROC_histograms/{}_ROC_histograms.root".format(title)

def Get_ROC_Efficiencies(histogram,ran,nCuts,print_cut=False):
        """Helper function used in Make_ROC_Curves(). Given a discriminant histogram, it finds the cut corresponding most closely to a 10% mistag rate"""
        Cuts = np.linspace(ran[0],ran[1],nCuts+1)
        bin_ran = (histogram.GetXaxis().FindBin(ran[0]),histogram.GetXaxis().FindBin(ran[1]))
        Efficiencies = np.zeros(nCuts+1)
        FullIntegral = histogram.Integral(bin_ran[0],bin_ran[1])
        for n,cut in enumerate(Cuts):
                bin_cut = histogram.GetXaxis().FindBin(cut)
                Efficiencies[n] = histogram.Integral(bin_cut,bin_ran[1])/FullIntegral
        diff = 1
        closest = 0
        if print_cut:
                for n,eff in enumerate(Efficiencies):
                        if abs(eff - 0.1) < diff:
                                closest = n
                                diff = abs(eff - 0.1)
                print "Mistag rate:",Efficiencies[closest], "corresponding to a cut at", Cuts[closest]
        return Efficiencies, Cuts[closest]

def find_cuts(file_path,bins):
        """given a pre-made root file containing binned ROC histograms (made by efficient_Make_Binned_ROC_histograms()) and the corresponding list of bins, it returns the 10% mistag rate cut on L4-L1, L4/L1 and CSV for each bin"""
        File = rt.TFile.Open(file_path)
        nbins = 60
	ANN_cuts, CSV_cuts = [], []

        for bin_ in range(len(bins)-1):
                print "\n"
                print "cuts for bin {}_{}: \n".format(bins[bin_],bins[bin_+1])
                print "ANN:"
                _ANN_eff, _ANN_cut =  Get_ROC_Efficiencies(File.Get("ANN_"+str(bins[bin_])+"_"+str(bins[bin_+1])),(0,1),nbins,print_cut=True)
                print "CSV:"
                _CSV_eff, _CSV_cut = Get_ROC_Efficiencies(File.Get("CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1])),(0,1),nbins,print_cut=True)
                print "\n"
		ANN_cuts.append(_ANN_cut)
		CSV_cuts.append(_CSV_cut)
	print "MANtag_cuts =",ANN_cuts
	print "CSV_cuts =",CSV_cuts
	return ANN_cuts, CSV_cuts

def Stacked_tagged_jets(title, Discriminant_Name, AllJets, Discriminant_and_not_CSV, CSV_and_not_Discriminant, Discriminant_and_CSV, y_max):

        AllJets_hist = AllJets.Clone()
        Discriminant_and_not_CSV_hist = Discriminant_and_not_CSV.Clone()
        CSV_and_not_Discriminant_hist = CSV_and_not_Discriminant.Clone()
        Discriminant_and_CSV_hist = Discriminant_and_CSV.Clone()

        stack = rt.THStack("stack", "stack")
        AllJets_hist.SetLineColor(1)
        Discriminant_and_not_CSV_hist.SetFillColor(2)
        CSV_and_not_Discriminant_hist.SetFillColor(3)
        Discriminant_and_CSV_hist.SetFillColor(4)
        stack.Add(Discriminant_and_CSV_hist)
        stack.Add(CSV_and_not_Discriminant_hist)
        stack.Add(Discriminant_and_not_CSV_hist)

        canvas = rt.TCanvas("canvas","canvas",600,600)
        rt.gStyle.SetOptTitle(0)
        rt.gStyle.SetOptStat(0)
        legend = rt.TLegend(0.9,0.9,0.65,0.75)
        legend.AddEntry(AllJets_hist,"all jets")
        legend.AddEntry(Discriminant_and_not_CSV_hist,Discriminant_Name+" and not CSV")
        legend.AddEntry(CSV_and_not_Discriminant_hist,"CSV and not "+Discriminant_Name)
        legend.AddEntry(Discriminant_and_CSV_hist, Discriminant_Name+" and CSV")
        AllJets_hist.GetXaxis().SetTitle("jet p_{T} (GeV)")
        AllJets_hist.GetYaxis().SetTitle('# jets')
        AllJets_hist.GetYaxis().SetTitleOffset(1.5)
        if y_max >0: AllJets_hist.SetMaximum(y_max)
        AllJets_hist.Draw()
        stack.Draw("SAME")
        rt.gPad.RedrawAxis()
        legend.Draw()
	f2 = rt.TFile("Plots/stacked_"+title+".root", "RECREATE")
	canvas.Write()
	f2.Close()
        #canvas.SaveAs("Plots/stacked_"+title+".png")

def ANN_binned_tagged_jets_hist(datalist, discriminant_cuts, CSV_cuts, bins, nbins):
        """creates a histogram for each dataset given as list of tuples (MANtag, pT, CSV, title, range) of all the jets that were b-tagged by passing a given a list of cut values corresponding to the given pT-bins for CSV and ANN versus jet-pT. The histograms are saved to a root file for further use."""
        title = "tagged_jets_vs_pT"
        discriminant = "ANN"
        AllJetsHistlist = []
        CSVHistlist = []
        DiscriminantHistlist = []
        for n,data in enumerate(datalist):
                datatitle = data[3]
                print "working on",datatitle
                ran = data[4]
                CSV = data[2]
                pT = data[1]
                MANtag = data[0]
                AllJetsHistlist.append(rt.TH1D(datatitle+"_AllJets",datatitle+"_"+title,nbins,ran[0],ran[1]))
                AllJetsHistlist[n].SetLineColor(4)
                CSVHistlist.append(rt.TH1D(datatitle+"_CSV",datatitle+"_"+title,nbins,ran[0],ran[1]))
                CSVHistlist[n].SetLineColor(3)
                DiscriminantHistlist.append(rt.TH1D(datatitle+"_Discriminant",datatitle+"_"+title,nbins,ran[0],ran[1]))
                DiscriminantHistlist[n].SetLineColor(2)

                bin_numbers = ANN_bin_selection(pT,bins)

                for i,pT_value in enumerate(pT):
                        if bin_numbers[i] == -100: continue
                        AllJetsHistlist[n].Fill(pT_value)
                        if MANtag[i] >= discriminant_cuts[bin_numbers[i]]: DiscriminantHistlist[n].Fill(pT_value)
                        if CSV[i] >= CSV_cuts[bin_numbers[i]]: CSVHistlist[n].Fill(pT_value)

        legendlist = []
        Tfilelist = []
        for n,data in enumerate(datalist):
                datatitle = data[3]
                rt.gStyle.SetOptStat(0)
                legendlist.append(rt.TLegend(0.9,0.9,0.65,0.75))
                legendlist[n].AddEntry(AllJetsHistlist[n], "All jets")
                legendlist[n].AddEntry(CSVHistlist[n], "CSV")
                legendlist[n].AddEntry(DiscriminantHistlist[n], discriminant)
                AllJetsHistlist[n].GetXaxis().SetTitle("jet_pT")
                AllJetsHistlist[n].GetYaxis().SetTitle('# jets')
                AllJetsHistlist[n].GetYaxis().SetTitleOffset(1.5)
                Tfilelist.append(rt.TFile("tagged_jets_histograms/"+title+"_"+datatitle+discriminant+".root","recreate"))
                print "saved histogram as tagged_jets_histograms/"+title+"_"+datatitle+discriminant+".root"
                AllJetsHistlist[n].Write()
                CSVHistlist[n].Write()
                DiscriminantHistlist[n].Write()


def ANN_exclusive_tagged_jets_hist(signal_title, MANtag, pT, CSV, discriminant_cuts, CSV_cuts, bins, ran, nbins, AllJets, y_max=0):
        """creates three histograms when given a signal dataset: -all jets tagged by CSV but not by ANN, -all jets tagged by ANN but not by CSV, -all jets tagged by both. the amount of jets is plotted versus the feature given as string to 'mode' (currently only jet-pT working) (see FeatureDict)"""
        title = signal_title+"_tagged_jets_vs_pT_exclusive"
        discriminant = "ANN"
        CSV_and_not_Discriminant = rt.TH1D(signal_title+"_CSV_not_Discriminant",title,nbins,ran[0],ran[1])
        CSV_and_not_Discriminant.SetLineColor(3)
        Discriminant_and_not_CSV = rt.TH1D(signal_title+"_Discriminant_and_not_CSV",title,nbins,ran[0],ran[1])
        Discriminant_and_not_CSV.SetLineColor(2)
        Discriminant_and_CSV = rt.TH1D(signal_title+"_Discriminant_and_CSV",title,nbins,ran[0],ran[1])
        Discriminant_and_CSV.SetLineColor(4)

        bin_numbers = ANN_bin_selection(pT,bins)

        for i,pT_value in enumerate(pT):
                if bin_numbers[i] == -100: continue
                CSV_tag, Disc_tag = False, False
                if CSV[i] >= CSV_cuts[bin_numbers[i]]:
                        CSV_tag = True
                if MANtag[i] >= discriminant_cuts[bin_numbers[i]]:
                	Disc_tag = True

                if Disc_tag and not CSV_tag: Discriminant_and_not_CSV.Fill(pT_value)
                if CSV_tag and not Disc_tag: CSV_and_not_Discriminant.Fill(pT_value)
                if CSV_tag and Disc_tag: Discriminant_and_CSV.Fill(pT_value)

        assert AllJets != None, "need to give a histogram containing full amount of jets into AllJets parameter"
        Stacked_tagged_jets(title, discriminant, AllJets, Discriminant_and_not_CSV, CSV_and_not_Discriminant, Discriminant_and_CSV, y_max)

        Tfile= rt.TFile("exclusive_tagged_histograms/"+title+".root","recreate")
        Discriminant_and_not_CSV.Write()
        CSV_and_not_Discriminant.Write()
        Discriminant_and_CSV.Write()
        print "saved as exclusive_tagged_histograms/"+title+".root"

def Relative_Gains(Discriminant_not_CSV, CSV_not_Discriminant, Discriminant_and_CSV):
        ran = (0,2500)
        res = 60
        gain_list = []
        errors = []
        CSV = CSV_not_Discriminant.Clone()
        CSV.Add(Discriminant_and_CSV)
        thresholds = np.linspace(ran[0],ran[1],60)
        for threshold in thresholds:
                N_Disc_excl = Discriminant_not_CSV.Integral(Discriminant_not_CSV.GetXaxis().FindBin(threshold),Discriminant_not_CSV.GetXaxis().FindBin(ran[1]))
                N_CSV = CSV.Integral(CSV.GetXaxis().FindBin(threshold),CSV.GetXaxis().FindBin(ran[1]))
                if N_CSV != 0:
                        gain_list.append(N_Disc_excl/N_CSV)
                        errors.append((np.sqrt(N_Disc_excl)/N_CSV) * np.sqrt(1+ N_Disc_excl/N_CSV))
                else:
                        gain_list.append(np.nan)
                        errors.append(np.nan)
        return (np.array(gain_list), np.array(errors), thresholds)

def Relative_Gain_Plots(input_file, signal_title, pT_range):

        ANN_withPT_file = rt.TFile.Open(input_file)

        withPT_Discriminant_not_CSV = ANN_withPT_file.Get(signal_title+"_Discriminant_and_not_CSV")
        withPT_CSV_not_Discriminant = ANN_withPT_file.Get(signal_title+"_CSV_not_Discriminant")
        withPT_Discrimant_and_CSV = ANN_withPT_file.Get(signal_title+"_Discriminant_and_CSV")

        ANN_withPT_gains, ANN_withPT_errors, thresholds = Relative_Gains(withPT_Discriminant_not_CSV, withPT_CSV_not_Discriminant, withPT_Discrimant_and_CSV)

        plt.figure()
        plt.errorbar(thresholds, ANN_withPT_gains, yerr=ANN_withPT_errors, fmt='magenta', label=r"MANtag")
        plt.xlim(pT_range[0],pT_range[1])
        plt.ylim(0,3.5)
        plt.xlabel(r"jet $p_T$ threshold (GeV)")
        plt.ylabel(r"gain")
        plt.legend(loc=2)
        plt.savefig("Plots/Relative_Gain_Plot.png")
        print "saved figure as Plots/Relative_Gain_Plot.png"
        #plt.show()

def Efficiency_vs_pT(title,histlist, hist_all_jets,y_max,legend_shift=False,BG=False, LargeLegend=False):
        """plots for each histogram of tagged jets given in a list of tuples (histogram, title, colorindex(optional)) the efficiency for each bin, where the x-axis corresponds to the feature given as string (see FeatureDict)."""
        canvas = rt.TCanvas('canvas','canvas',600,600)
        if legend_shift:
                if LargeLegend:
                        #legend = rt.TLegend(0.1,0.1,0.4,0.3)
                        legend = rt.TLegend(0.6,0.9,0.7,0.9)
                else:
                        #legend = rt.TLegend(0.1,0.1,0.35,0.25)
                        legend = rt.TLegend(0.65,0.75,0.9,0.9)
        else:
                if LargeLegend:
                        legend = rt.TLegend(0.1,0.9,0.4,0.7)
                else:
                        legend = rt.TLegend(0.1,0.9,0.35,0.75)
        graphlist = []
        for n,hist in enumerate(histlist):
                graphlist.append(rt.TGraphAsymmErrors())
                #if n==0: graphlist[n].SetTitle(title+"_vs_jet-pT")
                graphlist[n].Divide(hist[0],hist_all_jets,"cl=0.683 b(1,1) mode")
                legend.AddEntry(graphlist[n], histlist[n][1],"LEP")
                if len(hist) > 2:
                        graphlist[n].SetLineColor(hist[2])
                else:
                        if n < 3:
                                graphlist[n].SetLineColor(n+2)
                        else:
                                graphlist[n].SetLineColor(n+3)
                if n<1:
                        graphlist[n].GetXaxis().SetTitle("jet p_{T} (GeV)")
                        if BG:
                                graphlist[n].GetYaxis().SetTitle('mistag rate')
                        else:
                                graphlist[n].GetYaxis().SetTitle('efficiency')
                        graphlist[n].GetYaxis().SetTitleOffset(1.5)
                        graphlist[n].SetMinimum(0.)
                        graphlist[n].SetMaximum(y_max)
                        graphlist[n].Draw()
                else:
                        graphlist[n].Draw("SAME")
        legend.Draw()
        #canvas.SaveAs("Plots/"+title+"_vs_pT.png")
	f2 = rt.TFile("Plots/"+title+"_vs_pT.root","RECREATE")
	canvas.Write()
	f2.Close()

def AnalysisStep0_LoadData(signal_files, bg_files):
	signal_MANtag_data, signal_pT_data, signal_CSV_data = MANtag_pT_CSV_from_multiple_files(signal_files, MC=True, BG=False)
        bg_MANtag_data, bg_pT_data, bg_CSV_data = MANtag_pT_CSV_from_multiple_files(bg_files, MC=True, BG=True)
	return (signal_MANtag_data, signal_pT_data, signal_CSV_data), (bg_MANtag_data, bg_pT_data, bg_CSV_data)

def AnalysisStep1_DeriveCuts(cut_bins, bg_data):
        bg_MANtag_data, bg_pT_data, bg_CSV_data = bg_data

        ANN_Make_Binned_ROC_histograms("BG_large", bg_MANtag_data, bg_pT_data, bg_CSV_data, cut_bins)
        MANtag_cuts, CSV_cuts = find_cuts('ROC_histograms/BG_large_ROC_histograms.root',cut_bins)
	return MANtag_cuts, CSV_cuts

def AnalysisStep2_GeneratePlots(cut_bins, signal_data, bg_data, MANtag_cuts, CSV_cuts):
	signal_MANtag_data, signal_pT_data, signal_CSV_data = signal_data
	bg_MANtag_data, bg_pT_data, bg_CSV_data = bg_data
	ANN_binned_tagged_jets_hist([(signal_MANtag_data, signal_pT_data, signal_CSV_data, "signal_large", (200,3000))], MANtag_cuts, CSV_cuts, cut_bins, 60)
	ANN_binned_tagged_jets_hist([(bg_MANtag_data, bg_pT_data, bg_CSV_data, "BG_large", (200,3000))], MANtag_cuts, CSV_cuts, cut_bins, 60)
	
	signal_tagged_jets_file =   rt.TFile.Open("tagged_jets_histograms/tagged_jets_vs_pT_signal_largeANN.root")
        signal_AllJets_hist = signal_tagged_jets_file.Get('signal_large_AllJets')
	signal_ANN_hist = signal_tagged_jets_file.Get('signal_large_Discriminant')
	signal_CSV_hist = signal_tagged_jets_file.Get('signal_large_CSV')
	BG_tagged_jets_file =   rt.TFile.Open("tagged_jets_histograms/tagged_jets_vs_pT_BG_largeANN.root")
	BG_AllJets_hist = BG_tagged_jets_file.Get('BG_large_AllJets')
	BG_ANN_hist = BG_tagged_jets_file.Get('BG_large_Discriminant')
        BG_CSV_hist = BG_tagged_jets_file.Get('BG_large_CSV')
	
	ANN_exclusive_tagged_jets_hist("signal_large", signal_MANtag_data, signal_pT_data, signal_CSV_data, MANtag_cuts, CSV_cuts, cut_bins, (200,3000), 60, signal_AllJets_hist)

	Relative_Gain_Plots("exclusive_tagged_histograms/signal_large_tagged_jets_vs_pT_exclusive.root", "signal_large", (200,3000))	

	Efficiency_vs_pT("signal_large",[(signal_ANN_hist, "MANtag", 3), (signal_CSV_hist, "CSV", 4)], signal_AllJets_hist,1,BG=False)
	Efficiency_vs_pT("BG_large",[(BG_ANN_hist, "MANtag", 3), (BG_CSV_hist, "CSV", 4)], BG_AllJets_hist,0.3,BG=True)




if __name__ == "__main__":

	#model_path = "ANN/test/model_large.h5"
	#from AddMANtag_list import AddMANtag_to_dataset(model_path)
	
	signal_path = "/afs/cern.ch/work/m/msommerh/public/HitAnalyzer2/{}/M{}/flatTuple_{}.root"
        bg_path = "/afs/cern.ch/work/m/msommerh/public/HitAnalyzer2/QCD/{}/flatTuple_{}.root"

        M0_list = ['1000', '1200', '1400', '1600', '1800', '2000', '2500', '3000', '3500', '4000', '4500', '5000', '5500', '6000']
        bin_list = [('170','300'), ('300','470'), ('470','600'), ('600','800'), ('800','1000'), ('1000','1400'), ('1400','1800'), ('1800','2400'), ('2400','3200'), ('3200', 'Inf')]

        signal_files = []
        for M0 in M0_list[5:]:
                for i in range(11,16):
			#if i == 14 and M0 == '1000': continue
                        signal_files.append(signal_path.format('2017',M0,i))
	for M0 in M0_list[5:]:
                for i in range(11,16):
                        signal_files.append(signal_path.format('2018',M0,i))
        bg_files = []
        for bin_ in bin_list:
                for i in range(31,42):
			#if i == 38 and bin_[0] == '470': continue
			#if i == 35 and bin_[0] == '3200': continue
			#if i == 38 and bin_[0] == '3200': continue
                        bg_files.append(bg_path.format(bin_[0]+"to"+bin_[1],i))

	cut_bins = [0, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250, 3500, 3750, 4000, 5000]
	
	signal_data, bg_data = AnalysisStep0_LoadData(signal_files, bg_files)

	MANtag_cuts, CSV_cuts = AnalysisStep1_DeriveCuts(cut_bins, bg_data)

	#MANtag_cuts = [0.7166666666666667, 0.7166666666666667, 0.7166666666666667, 0.7166666666666667, 0.7166666666666667, 0.7166666666666667, 0.7166666666666667, 0.7166666666666667, 0.7166666666666667, 0.7, 0.6833333333333333, 0.6166666666666667, 0.6166666666666667, 0.55, 0.5, 0.3833333333333333]
	#CSV_cuts = [0.5333333333333333, 0.6, 0.6333333333333333, 0.65, 0.65, 0.6333333333333333, 0.6833333333333333, 0.7333333333333333, 0.7666666666666666, 0.75, 0.75, 0.7833333333333333, 0.7833333333333333, 0.7666666666666666, 0.7333333333333333, 0.75]

	AnalysisStep2_GeneratePlots(cut_bins, signal_data, bg_data, MANtag_cuts, CSV_cuts)









	
	'''	
	#signal_MANtag_data, signal_pT_data, signal_CSV_data = MANtag_pT_CSV_from_single_file("flatTuples/flatTuple_signal.root", MC=True, BG=False)	
	#bg_MANtag_data, bg_pT_data, bg_CSV_data = MANtag_pT_CSV_from_single_file("flatTuples/flatTuple_QCD0toinf.root", MC=True, BG=True)

	signal_MANtag_data, signal_pT_data, signal_CSV_data = MANtag_pT_CSV_from_multiple_files(signal_files, MC=True, BG=False)
	bg_MANtag_data, bg_pT_data, bg_CSV_data = MANtag_pT_CSV_from_multiple_files(bg_files, MC=True, BG=True)
	
	#ANN_Make_Binned_ROC_histograms("BG_large", bg_MANtag_data, bg_pT_data, bg_CSV_data, cut_bins)
	#find_cuts('ROC_histograms/BG_large_ROC_histograms.root',cut_bins)
	

	MANtag_cuts = [0.9166666666666666, 0.9, 0.8666666666666667, 0.85, 0.85, 0.8166666666666667, 0.7833333333333333, 0.7833333333333333, 0.75, 0.6833333333333333, 0.7166666666666667, 0.6833333333333333, 0.65, 0.6666666666666666, 0.5666666666666667, 0.6333333333333333]
	CSV_cuts = [0.5333333333333333, 0.6, 0.6333333333333333, 0.65, 0.65, 0.6333333333333333, 0.6833333333333333, 0.7333333333333333, 0.7666666666666666, 0.75, 0.75, 0.7833333333333333, 0.7833333333333333, 0.7666666666666666, 0.7333333333333333, 0.75]
	
	ANN_binned_tagged_jets_hist([(signal_MANtag_data, signal_pT_data, signal_CSV_data, "signal_large", (200,5000))], MANtag_cuts, CSV_cuts, cut_bins, 60)
	ANN_binned_tagged_jets_hist([(bg_MANtag_data, bg_pT_data, bg_CSV_data, "BG_large", (200,5000))], MANtag_cuts, CSV_cuts, cut_bins, 60)
	
	signal_tagged_jets_file =   rt.TFile.Open("tagged_jets_histograms/tagged_jets_vs_pT_signal_largeANN.root")
        signal_AllJets_hist = signal_tagged_jets_file.Get('signal_large_AllJets')
	signal_ANN_hist = signal_tagged_jets_file.Get('signal_large_Discriminant')
	signal_CSV_hist = signal_tagged_jets_file.Get('signal_large_CSV')
	BG_tagged_jets_file =   rt.TFile.Open("tagged_jets_histograms/tagged_jets_vs_pT_BG_largeANN.root")
	BG_AllJets_hist = BG_tagged_jets_file.Get('BG_large_AllJets')
	BG_ANN_hist = BG_tagged_jets_file.Get('BG_large_Discriminant')
        BG_CSV_hist = BG_tagged_jets_file.Get('BG_large_CSV')
	
	ANN_exclusive_tagged_jets_hist("signal_large", signal_MANtag_data, signal_pT_data, signal_CSV_data, MANtag_cuts, CSV_cuts, cut_bins, (200,5000), 60, signal_AllJets_hist)

	Relative_Gain_Plots("exclusive_tagged_histograms/signal_large_tagged_jets_vs_pT_exclusive.root", "signal_large", (200,5000))	

	Efficiency_vs_pT("signal_large",[(signal_ANN_hist, "MANtag", 3), (signal_CSV_hist, "CSV", 4)], signal_AllJets_hist,1,BG=False)
	Efficiency_vs_pT("BG_large",[(BG_ANN_hist, "MANtag", 3), (BG_CSV_hist, "CSV", 4)], BG_AllJets_hist,0.3,BG=True)
	'''
	
