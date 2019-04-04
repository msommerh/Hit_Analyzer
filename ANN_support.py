import numpy as np
import keras as kr


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
	"""nClusers is an numpy array of 20 entries corresponding to the number of clusters in each of the 4 layers in each of the 5 cone sizes"""
	nClusters = np.vstack((np.ndarray((0,20)),nClusters))
	pT = np.array([pt])
	X = ANN_functional_shape(nClusters)
        X.append(pT)
	return model.predict(X).item()

if __name__ == "__main__":
	
	#input_data = np.ndarray((0,20) , dtype=np.int32)
	#input_data = np.vstack((input_data,[ 1,  1,  2,  2,  7,  0,  0,  1,  3,  8,  2,  2,  3,  7, 16,  0,  1,  4,  7, 16]))
	input_data = np.array([ 1,  1,  2,  2,  7,  0,  0,  1,  3,  8,  2,  2,  3,  7, 16,  0,  1,  4,  7, 16])
	#input_data = np.vstack((input_data, np.ndarray((0,20))))
	model = kr.models.load_model("ANN/model.h5")
	print Predict_from_best_model(input_data, 58.4362373352, model)
