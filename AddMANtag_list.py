from AddMANtag import AddMANtag
import keras as kr
import numpy as np
from sys import argv

def AddMANtag_to_dataset(model_path):
	model = kr.models.load_model(model_path)
	np.seterr(divide='ignore', invalid='ignore')

	input_signal_path = "/eos/user/m/msommerh/HitAnalyzer4_collected/{}/M{}/flatTuple_{}.root"
	input_bg_path = "/eos/user/m/msommerh/HitAnalyzer4_collected/QCD/{}/flatTuple_{}.root"
	output_signal_path = "/afs/cern.ch/work/m/msommerh/public/HitAnalyzer4/{}/M{}/flatTuple_{}.root"
	output_bg_path = "/afs/cern.ch/work/m/msommerh/public/HitAnalyzer4/QCD/{}/flatTuple_{}.root"

	M0_list = ['1000', '1200', '1400', '1600', '1800', '2000', '2500', '3000', '3500', '4000', '4500', '5000', '5500', '6000']
	bin_list = [('170','300'), ('300','470'), ('470','600'), ('600','800'), ('800','1000'),('1000','1400'), ('1400','1800'), ('1800','2400'), ('2400','3200'), ('3200', 'Inf')]
	
	for year in ['2017', '2018']:
		for M0 in M0_list:
			for i in range(1,61):
				try:
					AddMANtag(model, input_signal_path.format(year, M0, i), output_signal_path.format(year, M0, i))
				except ReferenceError:
					print "skipping file"
					#print input_signal_path.format(year, M0, i)

	for bin_ in bin_list:
		for i in range(1,168):
			try:
				AddMANtag(model, input_bg_path.format(bin_[0]+"to"+bin_[1],i), output_bg_path.format(bin_[0]+"to"+bin_[1],i))
			except ReferenceError:
				print "skipping file"
				#print input_bg_path.format(bin_[0]+"to"+bin_[1],i)

if __name__ == "__main__":
	if len(argv) > 1:
		model_path = argv[1]
		AddMANtag_to_dataset(model_path)
	else:
		print "need to give model path as argument, e.g. ANN/test/model_large.h5"
