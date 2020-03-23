import numpy as np
import math
from collections import Counter
from collections import defaultdict
import scipy as sp

data = sp.genfromtxt("training1.tsv", delimiter="\t")
outcome = data[:,0]
training = data[:,1:]

data1 = sp.genfromtxt("test1.tsv", delimiter="\t")
testing = data[6:8,1:]

def LogPrior(outcome):
    no_of_samples = len(outcome)
    count = 0
    for i in outcome:
        if i == 1:
            count += 1
    P_C = count/no_of_samples
    P_C_bar = (no_of_samples-count)/no_of_samples
    LogPrior = math.log10(P_C/P_C_bar)
    print ("Probability that the proteins form a complex is: ",P_C)
    print ("Probability that there is no protein interaction is: ",P_C_bar)
    print ("Log prior odd is: ", LogPrior)
    return LogPrior


def occurrence(outcome):
    no_of_samples = len(outcome)
    prob = dict(Counter(outcome))
    for key in prob.keys():
        prob[key] = (prob[key]/ float(no_of_samples))
    return prob

def NB(training, outcome, testing):
    LogPrior(outcome)
    classes = np.unique(outcome)
    likelihoods = {}
    for cls in classes:
        likelihoods[cls] = defaultdict(list)
    
    for cls in classes:
        row_indices = np.where(outcome == cls)[0]
        subset = training[row_indices, :]
        
        r,c = np.shape(subset)
        for j in range(0, c):
            likelihoods[cls][j] += list(subset[:,j])
     
    for cls in classes:
        for j in range(0,c):
            likelihoods[cls][j]= occurrence(likelihoods[cls][j])     
       
    rows, cols = np.shape(testing)    
    
    class_probability = {}
    for cls in classes:
        class_probability[cls] = defaultdict(list)
        for i in range(cols):
            for j in range(rows):
                if new_sample[j][i] in list(likelihoods[cls][i]):
                    class_probability[cls][i] = math.log10(likelihoods[cls][i][testing[j][i]])
                else:
                    class_probability[cls][i] = 0
    
    likelihoods_ratio = {}
    abs_likelihoods_ratio = {}
    for i in range(cols):
        likelihoods_ratio[i] = (class_probability[1][i] - class_probability[0][i])
        abs_likelihoods_ratio[i] = abs(class_probability[1][i] - class_probability[0][i])
    print (likelihoods_ratio)
    print ("The ten features with the highest absolute log-likelihood ratios are: ")
    for key in dict(Counter(abs_likelihoods_ratio).most_common(10)).keys():
        print ("Feature ", key,": ",likelihoods_ratio[key])
'''
    LogPosterior = {}
    for i in likelihoods_ratio.keys():
        LogPosterior = likelihoods_ratio + LogPrior(outcome)
    
    prediction = []    
    if LogPosterior < 0:
        prediction.append("0")
    else:
        prediction.append("1")           
    prediction = np.asarray(prediction)
    print (prediction)

def accuracy(testing, prediction):
    correct = 0
    for i in testing[:,0]:
        if i == prediction[i]:
            correct += 1
    return (correct / len(testing)) * 100.0
    
'''
if __name__ == "__main__":
    training   = np.asarray(training)
    outcome    = np.asarray(outcome)
    new_sample = np.asarray(testing)
    NB(training, outcome, new_sample) 
