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

def MANtag_pT_CSV_from_file(input_file, MC=True, BG=False, ZPrime_matching=False):
	MANtag = []
	pT = []
	CSV = []
	deep_CSV = []
	PV = []
	#print "opening file:", input_file
	f1 = rt.TFile.Open(input_file, "READ")
	try:
		tree = f1.Get("demo/tree")
	except ReferenceError:
		print "skipping file"
		return MANtag, pT, CSV, deep_CSV, PV
	N = tree.GetEntries()
	#print N, "events"
	for i in xrange(N):
		tree.GetEntry(i)
		for j in xrange(tree.nJets):
			if MC:
				if BG and tree.jet_MC_bTag[j] != 0: continue
				if not BG and tree.jet_MC_bTag[j] != 1: continue
				if not BG and ZPrime_matching and abs(tree.jet_BMotherID[j]) != 32: continue
			MANtag.append(tree.MANtag[j])
			pT.append(tree.jet_pt[j])
			CSV.append(tree.jet_bTag[j]) ##CSVv2
			deep_CSV.append(tree.jet_deepCSV_probb[j]) ##deepCSV
			PV.append(tree.nPV)
	return MANtag, pT, CSV, deep_CSV, PV

def MANtag_pT_CSV_from_single_file(input_file, MC=True, BG=False, ZPrime_matching=False):
	MANtag, pT, CSV, deep_CSV, PV= MANtag_pT_CSV_from_file(input_file, MC=MC, BG=BG, ZPrime_matching=ZPrime_matching)
	return np.array(MANtag), np.array(pT), np.array(CSV), np.array(deep_CSV), np.array(PV)

def MANtag_pT_CSV_from_multiple_files(file_list, MC=True, BG=False, ZPrime_matching=False):
	MANtag, pT, CSV, deep_CSV, PV = [], [], [], [], []
	for f in file_list:
		_MANtag, _pT, _CSV, _deep_CSV, _PV = MANtag_pT_CSV_from_file(f, MC=MC, BG=BG, ZPrime_matching=ZPrime_matching)
		MANtag += _MANtag
		pT += _pT
		CSV += _CSV
		deep_CSV += _deep_CSV
		PV += _PV
	return np.array(MANtag), np.array(pT), np.array(CSV), np.array(deep_CSV), np.array(PV)

def ANN_Make_Binned_ROC_histograms(title, MANtag, pT, CSV, deep_CSV, bins):
        """makes binned ROC histograms for an ANN. takes as input a keras model, the necessary ANN data, pT, CSV and the desired pT-binning"""
	rt.gROOT.SetBatch(True)
	N = len(MANtag)
	assert len(pT) == N and len(CSV) == N and len(deep_CSV) == N

        nbins = 60
        ANN_hist_list = []
        CSV_hist_list = []
	deep_CSV_hist_list = []
        for bin_ in range(len(bins)-1):
                ANN_hist_list.append(rt.TH1D("ANN_"+str(bins[bin_])+"_"+str(bins[bin_+1]),"ANN_"+str(bins[bin_])+"_"+str(bins[bin_+1]),nbins,0,1))
                CSV_hist_list.append(rt.TH1D("CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1]),"CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1]),nbins,0,1))
                deep_CSV_hist_list.append(rt.TH1D("deep_CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1]),"deep_CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1]),nbins,0,1))

        bin_numbers = ANN_bin_selection(pT,bins)

        for n in xrange(N):
                if bin_numbers[n] == -100: continue
                ANN_hist_list[int(bin_numbers[n])].Fill(MANtag[n])
                CSV_hist_list[int(bin_numbers[n])].Fill(CSV[n])
                deep_CSV_hist_list[int(bin_numbers[n])].Fill(deep_CSV[n])

        tfile = rt.TFile("ROC_histograms/{}_ROC_histograms.root".format(title),"recreate")
        for hist in ANN_hist_list:
                hist.Write()
        for hist in CSV_hist_list:
                hist.Write()
        for hist in deep_CSV_hist_list:
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
        for n,eff in enumerate(Efficiencies):
                if abs(eff - 0.1) < diff:
                        closest = n
                        diff = abs(eff - 0.1)
        if print_cut:
                print "Mistag rate:",Efficiencies[closest], "corresponding to a cut at", Cuts[closest]
        return Efficiencies, Cuts[closest]

def find_cuts(file_path,bins):
        """given a pre-made root file containing binned ROC histograms (made by efficient_Make_Binned_ROC_histograms()) and the corresponding list of bins, it returns the 10% mistag rate cut on L4-L1, L4/L1 and CSV for each bin"""
        File = rt.TFile.Open(file_path)
        nbins = 60
	ANN_cuts, CSV_cuts, deep_CSV_cuts = [], [], []

        for bin_ in range(len(bins)-1):
                _ANN_eff, _ANN_cut =  Get_ROC_Efficiencies(File.Get("ANN_"+str(bins[bin_])+"_"+str(bins[bin_+1])),(0,1),nbins,print_cut=False)
                _CSV_eff, _CSV_cut = Get_ROC_Efficiencies(File.Get("CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1])),(0,1),nbins,print_cut=False)
                _deep_CSV_eff, _deep_CSV_cut = Get_ROC_Efficiencies(File.Get("deep_CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1])),(0,1),nbins,print_cut=False)
		ANN_cuts.append(_ANN_cut)
		CSV_cuts.append(_CSV_cut)
		deep_CSV_cuts.append(_deep_CSV_cut)
	print "MANtag_cuts =",ANN_cuts
	print "CSV_cuts =",CSV_cuts
	print "deep_CSV_cuts =",deep_CSV_cuts
	return ANN_cuts, CSV_cuts, deep_CSV_cuts

def Stacked_tagged_jets(title, Discriminant_Name, AllJets, Discriminant_and_not_CSV, CSV_and_not_Discriminant, Discriminant_and_CSV, y_max):
	rt.gROOT.SetBatch(True)
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

def ANN_binned_tagged_jets_hist(datalist, discriminant_cuts, CSV_cuts, deep_CSV_cuts, bins, nbins):
        """creates a histogram for each dataset given as list of tuples (MANtag, pT, CSV, title, range) of all the jets that were b-tagged by passing a given a list of cut values corresponding to the given pT-bins for CSV and ANN versus jet-pT. The histograms are saved to a root file for further use."""
	rt.gROOT.SetBatch(True)
        title = "tagged_jets_vs_pT"
        discriminant = "ANN"
        AllJetsHistlist = []
        CSVHistlist = []
        deep_CSVHistlist = []
        DiscriminantHistlist = []
        for n,data in enumerate(datalist):
                datatitle = data[4]
                print "working on",datatitle
                ran = data[5]
                CSV = data[2]
		deep_CSV = data[3]
                pT = data[1]
                MANtag = data[0]
                AllJetsHistlist.append(rt.TH1D(datatitle+"_AllJets",datatitle+"_"+title,nbins,ran[0],ran[1]))
                AllJetsHistlist[n].SetLineColor(4)
                CSVHistlist.append(rt.TH1D(datatitle+"_CSV",datatitle+"_"+title,nbins,ran[0],ran[1]))
                CSVHistlist[n].SetLineColor(3)
                deep_CSVHistlist.append(rt.TH1D(datatitle+"_deep_CSV",datatitle+"_"+title,nbins,ran[0],ran[1]))
                deep_CSVHistlist[n].SetLineColor(6)
                DiscriminantHistlist.append(rt.TH1D(datatitle+"_Discriminant",datatitle+"_"+title,nbins,ran[0],ran[1]))
                DiscriminantHistlist[n].SetLineColor(2)

                bin_numbers = ANN_bin_selection(pT,bins)

                for i,pT_value in enumerate(pT):
                        if bin_numbers[i] == -100: continue
                        AllJetsHistlist[n].Fill(pT_value)
                        if MANtag[i] >= discriminant_cuts[bin_numbers[i]]: DiscriminantHistlist[n].Fill(pT_value)
                        if CSV[i] >= CSV_cuts[bin_numbers[i]]: CSVHistlist[n].Fill(pT_value)
                        if deep_CSV[i] >= deep_CSV_cuts[bin_numbers[i]]: deep_CSVHistlist[n].Fill(pT_value)

        legendlist = []
        Tfilelist = []
        for n,data in enumerate(datalist):
                datatitle = data[4]
                rt.gStyle.SetOptStat(0)
                legendlist.append(rt.TLegend(0.9,0.9,0.65,0.75))
                legendlist[n].AddEntry(AllJetsHistlist[n], "All jets")
                legendlist[n].AddEntry(CSVHistlist[n], "CSV")
                legendlist[n].AddEntry(deep_CSVHistlist[n], "deep_CSV")
                legendlist[n].AddEntry(DiscriminantHistlist[n], discriminant)
                AllJetsHistlist[n].GetXaxis().SetTitle("jet_pT")
                AllJetsHistlist[n].GetYaxis().SetTitle('# jets')
                AllJetsHistlist[n].GetYaxis().SetTitleOffset(1.5)
                Tfilelist.append(rt.TFile("tagged_jets_histograms/"+title+"_"+datatitle+discriminant+".root","recreate"))
                print "saved histogram as tagged_jets_histograms/"+title+"_"+datatitle+discriminant+".root"
                AllJetsHistlist[n].Write()
                CSVHistlist[n].Write()
                deep_CSVHistlist[n].Write()
                DiscriminantHistlist[n].Write()


def ANN_exclusive_tagged_jets_hist(signal_title, MANtag, pT, CSV, discriminant_cuts, CSV_cuts, bins, ran, nbins, AllJets, y_max=0):
        """creates three histograms when given a signal dataset: -all jets tagged by CSV but not by ANN, -all jets tagged by ANN but not by CSV, -all jets tagged by both. the amount of jets is plotted versus the feature given as string to 'mode' (currently only jet-pT working) (see FeatureDict)"""
	rt.gROOT.SetBatch(True)
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
        plt.ylim(0,2.0)
        plt.xlabel(r"jet $p_T$ threshold (GeV)")
        plt.ylabel(r"gain")
        plt.legend(loc=2)
        plt.savefig("Plots/Relative_Gain_Plot.png")
        print "saved figure as Plots/Relative_Gain_Plot.png"
        #plt.show()

def Relative_Gain_Plots2(title, input_file_CSV2, input_file_deepCSV, signal_title_CSV2, signal_title_deepCSV, pT_range):

        CSV2_ANN_withPT_file = rt.TFile.Open(input_file_CSV2)

        withPT_Discriminant_not_CSV2 = CSV2_ANN_withPT_file.Get(signal_title_CSV2+"_Discriminant_and_not_CSV")
        withPT_CSV2_not_Discriminant = CSV2_ANN_withPT_file.Get(signal_title_CSV2+"_CSV_not_Discriminant")
        withPT_Discrimant_and_CSV2 = CSV2_ANN_withPT_file.Get(signal_title_CSV2+"_Discriminant_and_CSV")

        CSV2_ANN_withPT_gains, CSV2_ANN_withPT_errors, CSV2_thresholds = Relative_Gains(withPT_Discriminant_not_CSV2, withPT_CSV2_not_Discriminant, withPT_Discrimant_and_CSV2)

        deepCSV_ANN_withPT_file = rt.TFile.Open(input_file_deepCSV)

        withPT_Discriminant_not_deepCSV = deepCSV_ANN_withPT_file.Get(signal_title_deepCSV+"_Discriminant_and_not_CSV")
        withPT_deepCSV_not_Discriminant = deepCSV_ANN_withPT_file.Get(signal_title_deepCSV+"_CSV_not_Discriminant")
        withPT_Discrimant_and_deepCSV = deepCSV_ANN_withPT_file.Get(signal_title_deepCSV+"_Discriminant_and_CSV")

        deepCSV_ANN_withPT_gains, deepCSV_ANN_withPT_errors, deepCSV_thresholds = Relative_Gains(withPT_Discriminant_not_deepCSV, withPT_deepCSV_not_Discriminant, withPT_Discrimant_and_deepCSV)

        plt.figure()
        plt.errorbar(CSV2_thresholds, CSV2_ANN_withPT_gains, yerr=CSV2_ANN_withPT_errors, fmt='blue', label=r"MANtag w.r.t. CSVv2")
	plt.errorbar(deepCSV_thresholds, deepCSV_ANN_withPT_gains, yerr=deepCSV_ANN_withPT_errors, fmt='magenta', label=r"MANtag w.r.t. deepCSV")
        plt.xlim(pT_range[0],pT_range[1])
        plt.ylim(0,2.0)
        plt.xlabel(r"jet $p_T$ threshold (GeV)")
        plt.ylabel(r"gain")
        plt.legend(loc=2)
	plt.grid(True)
        plt.savefig("Plots/Relative_Gain_Plot2_{}.png".format(title))
        print "saved figure as Plots/Relative_Gain_Plot2_{}.png".format(title)

	import csv
	with open("Plots/Relative_Gains_CSV2_{}.csv".format(title), "w") as csvFile:
		writer = csv.writer(csvFile)
		writer.writerow(CSV2_thresholds)
		writer.writerow(CSV2_ANN_withPT_gains)
		writer.writerow(CSV2_ANN_withPT_errors)
	csvFile.close()
	print "saved gains wrt CSV2 as Plots/Relative_Gains_CSV2_{}.csv".format(title)
	with open("Plots/Relative_Gains_deepCSV_{}.csv".format(title), "w") as csvFile:
		writer = csv.writer(csvFile)
		writer.writerow(deepCSV_thresholds)
		writer.writerow(deepCSV_ANN_withPT_gains)
		writer.writerow(deepCSV_ANN_withPT_errors)
	csvFile.close()
	print "saved gains wrt deepCSV as Plots/Relative_Gains_deepCSV_{}.csv".format(title)

def Efficiency_vs_pT(title,histlist, hist_all_jets,y_max,legend_shift=False,BG=False, LargeLegend=False):
        """plots for each histogram of tagged jets given in a list of tuples (histogram, title, colorindex(optional)) the efficiency for each bin, where the x-axis corresponds to the feature given as string (see FeatureDict)."""
	rt.gROOT.SetBatch(True)
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
		graphlist[n].SetName(histlist[n][1])
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
	for graph in graphlist:
		graph.Write()
	f2.Close()
	print "saved as Plots/"+title+"_vs_pT.root"

def DrawROCs(title, signal_MANtag, signal_CSV, signal_deepCSV, signal_pT, bg_MANtag, bg_CSV, bg_deepCSV, bg_pT, pT_range):
	from sklearn.metrics import roc_curve
	labels = np.concatenate((np.ones(signal_pT.size), np.zeros(bg_pT.size)))
	condition = np.logical_and((np.concatenate((signal_pT,bg_pT))>=pT_range[0]), (np.concatenate((signal_pT,bg_pT))<pT_range[1]))
	MANtag_fpr, MANtag_tpr, MANtag_thr = roc_curve(labels[condition], np.concatenate((signal_MANtag, bg_MANtag))[condition])
	CSV_fpr, CSV_tpr, CSV_thr = roc_curve(labels[condition], np.concatenate((signal_CSV, bg_CSV))[condition])
	deepCSV_fpr, deepCSV_tpr, deepCSV_thr = roc_curve(labels[condition], np.concatenate((signal_deepCSV, bg_deepCSV))[condition])

	plt.figure()
        plt.semilogy(MANtag_tpr, MANtag_fpr, color='green', label=r"MANtag")
        plt.semilogy(CSV_tpr, CSV_fpr, color='blue', label=r"CSVv2")
        plt.semilogy(deepCSV_tpr, deepCSV_fpr, color='magenta', label=r"deepCSV")
	plt.semilogy([0,1],[0.1,0.1],'k:')
        plt.xlim(0,1)
        plt.ylim(1e-3,1.)
        plt.xlabel(r"signal efficiency")
        plt.ylabel(r"mistag rate")
        plt.legend(loc=4)
	plt.figtext(0.16,0.83,r'{}<$p_T$<{}'.format(*pT_range))
        plt.savefig("Plots/ROCS_{}_pT{}to{}.png".format(title,pT_range[0],pT_range[1]))
        print "saved figure as Plots/ROCS_{}_pT{}to{}.png".format(title,pT_range[0],pT_range[1])

def ANN_bin_selection(pT,bins):
        """numpy array version of the bin_selection() function: takes an array of all pT values and the pT-bins and returns an array of same length as pT with the corresponing bins index at each entry. pT values outside the bins are labeled with -100"""
        bin_numbers = np.zeros(len(pT))
        for n in range(len(bins)-1):
                bin_numbers += (n+100)*(pT>bins[n])*(pT<bins[n+1])
        bin_numbers -=100
        return bin_numbers.astype(int)

def Efficiency_vs_PU(output_file, MANtag, CSV2, deepCSV, pT, PV, ANN_Cuts, CSV2_Cuts, deepCSV_Cuts, bins, y_max, pT_Cut=200, BG=False):
        assert len(MANtag) == len(pT) == len(CSV2) == len(deepCSV), "data inputs need to have the same length"

        ran = (0,80)
        nbins = 80
        import array
        if BG:
                bins_ = array.array('d',[0.0, 11.0]+range(19,41,8)+[42.0,  52.0, 80])
        else:
                bins_ = array.array('d',[0.0, 11.0]+range(15,41,4)+[42.0, 52.0, 58.0, 65.0, 80])

        if pT_Cut >= 1200:
                bins_ = array.array('d',[0.0, 20.0, 40.0, 80.0])

        AllJets_Hist = rt.TH1D("AllJets","AllJets",nbins,ran[0],ran[1])
        ANN_Hist = rt.TH1D("MANtag","MANtag",nbins,ran[0],ran[1])
        CSV2_Hist = rt.TH1D("CSVv2","CSVv2",nbins,ran[0],ran[1])
	deepCSV_Hist = rt.TH1D("deepCSV","deepCSV",nbins,ran[0],ran[1])
        AllJets_Hist = AllJets_Hist.Rebin(len(bins_)-1,"AllJets",bins_)
        ANN_Hist = ANN_Hist.Rebin(len(bins_)-1,"MANtag",bins_)
        CSV2_Hist = CSV2_Hist.Rebin(len(bins_)-1,"CSVv2",bins_)
        deepCSV_Hist = deepCSV_Hist.Rebin(len(bins_)-1,"deepCSV",bins_)

        bin_numbers = ANN_bin_selection(pT,bins)

        for i,pT_value in enumerate(pT):
                        if pT_value < pT_Cut: continue
                        if bin_numbers[i] == -100: continue
                        AllJets_Hist.Fill(PV[i])
                        if CSV2[i] >= CSV2_Cuts[bin_numbers[i]]: CSV2_Hist.Fill(PV[i])
			if deepCSV[i] >= deepCSV_Cuts[bin_numbers[i]]: deepCSV_Hist.Fill(PV[i])
                        if MANtag[i] >= ANN_Cuts[bin_numbers[i]]: ANN_Hist.Fill(PV[i])

        canvas = rt.TCanvas('canvas','canvas',600,600)
        rt.gStyle.SetOptTitle(0)
        legend = rt.TLegend(0.1,0.9,0.35,0.75)
        ANN_Graph = rt.TGraphAsymmErrors()
        CSV2_Graph = rt.TGraphAsymmErrors()
	deepCSV_Graph = rt.TGraphAsymmErrors()
        ANN_Graph.Divide(ANN_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
        CSV2_Graph.Divide(CSV2_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
        deepCSV_Graph.Divide(deepCSV_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
	ANN_Graph.SetLineColor(8)
        CSV2_Graph.SetLineColor(4)
	deepCSV_Graph.SetLineColor(6)
        legend.AddEntry(ANN_Graph, "MANtag", "LEP")
        legend.AddEntry(CSV2_Graph, "CSVv2", "LEP")
	legend.AddEntry(deepCSV_Graph, "deepCSV", "LEP")
        if BG:
                ANN_Graph.GetYaxis().SetTitle('mistag rate')
        else:
                ANN_Graph.GetYaxis().SetTitle('efficiency')
	ANN_Graph.GetXaxis().SetTitle('nPV')
        ANN_Graph.GetYaxis().SetTitleOffset(1.5)
        ANN_Graph.SetMinimum(0.)
        ANN_Graph.SetMaximum(y_max)
        ANN_Graph.Draw()
        CSV2_Graph.Draw("SAME")
        deepCSV_Graph.Draw("SAME")
	legend.Draw()
        canvas.SaveAs(output_file)
        
def RebinHist(hist,name,bins):
        """Rebins an efficiency-vs-pT histogram such that the lower statistics at high pT are taken into account"""
        import array
        bins_ = array.array('d',bins)
        return  hist.Rebin(len(bins_)-1,"rebinned_"+name,bins_)

def AnalysisStep0_LoadData(signal_files, bg_files, ZPrime_matching=False):
	signal_MANtag_data, signal_pT_data, signal_CSV_data, signal_deep_CSV_data, signal_PV = MANtag_pT_CSV_from_multiple_files(signal_files, MC=True, BG=False, ZPrime_matching=ZPrime_matching)
        bg_MANtag_data, bg_pT_data, bg_CSV_data, bg_deep_CSV_data, bg_PV = MANtag_pT_CSV_from_multiple_files(bg_files, MC=True, BG=True)
	return (signal_MANtag_data, signal_pT_data, signal_CSV_data, signal_deep_CSV_data, signal_PV), (bg_MANtag_data, bg_pT_data, bg_CSV_data, bg_deep_CSV_data, bg_PV)

def AnalysisStep1_DeriveCuts(cut_bins, bg_data):
        bg_MANtag_data, bg_pT_data, bg_CSV_data, bg_deep_CSV_data, PV = bg_data

        ANN_Make_Binned_ROC_histograms("BG_large", bg_MANtag_data, bg_pT_data, bg_CSV_data, bg_deep_CSV_data, cut_bins)
        MANtag_cuts, CSV_cuts, deep_CSV_cuts = find_cuts('ROC_histograms/BG_large_ROC_histograms.root',cut_bins)
	return MANtag_cuts, CSV_cuts, deep_CSV_cuts

def AnalysisStep2_GeneratePlots(title, cut_bins, signal_data, bg_data, MANtag_cuts, CSV_cuts, deep_CSV_cuts):
	signal_MANtag_data, signal_pT_data, signal_CSV_data, signal_deep_CSV_data, signal_PV = signal_data
	bg_MANtag_data, bg_pT_data, bg_CSV_data, bg_deep_CSV_data, bg_PV = bg_data
	ANN_binned_tagged_jets_hist([(signal_MANtag_data, signal_pT_data, signal_CSV_data, signal_deep_CSV_data, "signal_"+title, (200,2500))], MANtag_cuts, CSV_cuts, deep_CSV_cuts, cut_bins, 60)
	ANN_binned_tagged_jets_hist([(bg_MANtag_data, bg_pT_data, bg_CSV_data, bg_deep_CSV_data, "BG_"+title, (200,2500))], MANtag_cuts, CSV_cuts, deep_CSV_cuts, cut_bins, 60)
	
	plot_bins = [0, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 5000]
	#plot_bins = [0,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,1650,1700,1750,1800,1850,1900,1950,2000,2050,2100,2150,2200,2250,2300,3000] 
	#plot_bins = [0]
	#for num in range(500,3001, 5):
	#	plot_bins.append(num)
	signal_tagged_jets_file =   rt.TFile.Open("tagged_jets_histograms/tagged_jets_vs_pT_signal_{}ANN.root".format(title))
        signal_AllJets_hist = signal_tagged_jets_file.Get('signal_{}_AllJets'.format(title))
	signal_ANN_hist = RebinHist(signal_tagged_jets_file.Get('signal_{}_Discriminant'.format(title)), "ANN", plot_bins)
	signal_CSV_hist = RebinHist(signal_tagged_jets_file.Get('signal_{}_CSV'.format(title)), "CSV", plot_bins)
	signal_deep_CSV_hist = RebinHist(signal_tagged_jets_file.Get('signal_{}_deep_CSV'.format(title)), "deepCSV", plot_bins)
	BG_tagged_jets_file =   rt.TFile.Open("tagged_jets_histograms/tagged_jets_vs_pT_BG_{}ANN.root".format(title))
	BG_AllJets_hist = RebinHist(BG_tagged_jets_file.Get('BG_{}_AllJets'.format(title)), "AllJets", plot_bins)
	BG_ANN_hist = RebinHist(BG_tagged_jets_file.Get('BG_{}_Discriminant'.format(title)), "ANN", plot_bins)
        BG_CSV_hist = RebinHist(BG_tagged_jets_file.Get('BG_{}_CSV'.format(title)), "CSV", plot_bins)
	BG_deep_CSV_hist = RebinHist(BG_tagged_jets_file.Get('BG_{}_deep_CSV'.format(title)), "deepCSV", plot_bins)
	
	ANN_exclusive_tagged_jets_hist("signal_{}_CSV2".format(title), signal_MANtag_data, signal_pT_data, signal_CSV_data, MANtag_cuts, CSV_cuts, cut_bins, (200,2500), 60, signal_AllJets_hist)
	ANN_exclusive_tagged_jets_hist("signal_{}_deepCSV".format(title), signal_MANtag_data, signal_pT_data, signal_deep_CSV_data, MANtag_cuts, deep_CSV_cuts, cut_bins, (200,2500), 60, signal_AllJets_hist)

	#Relative_Gain_Plots("exclusive_tagged_histograms/signal_large_CSV2_tagged_jets_vs_pT_exclusive.root", "signal_large_CSV2", (200,2000))
	#Relative_Gain_Plots("exclusive_tagged_histograms/signal_large_deepCSV_tagged_jets_vs_pT_exclusive.root", "signal_large_deepCSV", (200,2000))
	#Relative_Gain_Plots2(title, "exclusive_tagged_histograms/signal_{}_CSV2_tagged_jets_vs_pT_exclusive.root".format(title), "exclusive_tagged_histograms/signal_{}_deepCSV_tagged_jets_vs_pT_exclusive.root".format(title), "signal_{}_CSV2".format(title), "signal_{}_deepCSV".format(title), (200,2000))

	signal_AllJets_hist = RebinHist(signal_tagged_jets_file.Get('signal_{}_AllJets'.format(title)), "AllJets", plot_bins)

	Efficiency_vs_pT("signal_"+title,[(signal_ANN_hist, "MANtag", 8), (signal_CSV_hist, "CSVv2", 4), (signal_deep_CSV_hist, "deepCSV", 6)], signal_AllJets_hist,0.7, legend_shift=True, BG=False)
	Efficiency_vs_pT("BG_"+title,[(BG_ANN_hist, "MANtag", 8), (BG_CSV_hist, "CSVv2", 4), (BG_deep_CSV_hist, "deepCSV", 6)], BG_AllJets_hist,0.3,BG=True)

	Efficiency_vs_pT("onlyCSV_signal_"+title,[(signal_CSV_hist, "CSVv2", 4), (signal_deep_CSV_hist, "deepCSV", 6)], signal_AllJets_hist,0.7, legend_shift=True, BG=False)
        Efficiency_vs_pT("onlyCSV_BG_"+title,[(BG_CSV_hist, "CSVv2", 4), (BG_deep_CSV_hist, "deepCSV", 6)], BG_AllJets_hist,0.3,BG=True)

	Efficiency_vs_PU("Plots/"+title+"eff_vs_PU.png", signal_MANtag_data, signal_CSV_data, signal_deep_CSV_data, signal_pT_data, signal_PV, MANtag_cuts, CSV_cuts, deep_CSV_cuts, cut_bins, 0.7, pT_Cut=200, BG=False)
	Efficiency_vs_PU("Plots/"+title+"mistag_vs_PU.png", bg_MANtag_data, bg_CSV_data, bg_deep_CSV_data, bg_pT_data, bg_PV, MANtag_cuts, CSV_cuts, deep_CSV_cuts, cut_bins, 0.3, pT_Cut=200, BG=True)

if __name__ == "__main__":

	#model_path = "ANN/test/model_large.h5"
	#from AddMANtag_list import AddMANtag_to_dataset(model_path)
	
	signal_path = "/afs/cern.ch/work/m/msommerh/public/HitAnalyzer_updated_ptcut/{}/M{}/flatTuple_{}.root"
        bg_path = "/afs/cern.ch/work/m/msommerh/public/HitAnalyzer_updated_ptcut/QCD/{}/flatTuple_{}.root"

        M0_list = ['1000', '1200', '1400', '1600', '1800', '2000', '2500', '3000', '3500', '4000', '4500', '5000', '5500', '6000']
        bin_list = [('170','300'), ('300','470'), ('470','600'), ('600','800'), ('800','1000'), ('1000','1400'), ('1400','1800'), ('1800','2400'), ('2400','3200'), ('3200', 'Inf')]

        signal_files = []
        for M0 in M0_list[5:]:
                for i in range(41,61):
                        signal_files.append(signal_path.format('2017',M0,i))
	for M0 in M0_list[5:]:
                for i in range(41,61):
                        signal_files.append(signal_path.format('2018',M0,i))
        bg_files = []
        for bin_ in bin_list:
                for i in range(61,168):
                        bg_files.append(bg_path.format(bin_[0]+"to"+bin_[1],i))

	cut_bins = [0, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250, 3500, 3750, 4000, 5000]
	
	signal_data, bg_data = AnalysisStep0_LoadData(signal_files, bg_files, ZPrime_matching=False)

	#signal_MANtag, signal_pT, signal_CSV, signal_deepCSV = signal_data
	#bg_MANtag, bg_pT, bg_CSV, bg_deepCSV = bg_data
	#DrawROCs("summary_ZPrime-matching", signal_MANtag, signal_CSV, signal_deepCSV, signal_pT, bg_MANtag, bg_CSV, bg_deepCSV, bg_pT, (0,2500))
	#DrawROCs("summary_ZPrime-matching", signal_MANtag, signal_CSV, signal_deepCSV, signal_pT, bg_MANtag, bg_CSV, bg_deepCSV, bg_pT, (200,1200))
	#DrawROCs("summary_ZPrime-matching", signal_MANtag, signal_CSV, signal_deepCSV, signal_pT, bg_MANtag, bg_CSV, bg_deepCSV, bg_pT, (1200,2500))

	MANtag_cuts, CSV_cuts, deep_CSV_cuts = AnalysisStep1_DeriveCuts(cut_bins, bg_data)

	AnalysisStep2_GeneratePlots("updated_ptcut", cut_bins, signal_data, bg_data, MANtag_cuts, CSV_cuts, deep_CSV_cuts)


