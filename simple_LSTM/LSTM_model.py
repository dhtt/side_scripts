#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 21:11:51 2019
"""

"""Import required packages"""
import pandas as pd
import numpy as np
import tensorflow as tf
import pydot, graphviz, pydot_ng
import matplotlib.pyplot as plt
from tensorflow.python import keras
from tensorflow.python.keras.models import Sequential
from tensorflow.python.keras.layers import LSTM
from tensorflow.python.keras.utils import plot_model
from tensorflow.python.keras.callbacks import EarlyStopping, TensorBoard, ModelCheckpoint
#Ignore warning messages while training
import os
import warnings
warnings.filterwarnings("ignore")
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
#Report running time
from time import time
start_time = time()


"""Dataset specification (as described in dataset_readme.txt)"""

peptide_no = 700 
feature_no = 57
aminoacid_no = 22 #Aminoacid: 0-22
secondarystructure_no = 9 #Secondary structure: 22-31
terminal_no = 2 #Terminal C or N: 31-33
solventaccess_no = 2 #Relative and absolute solvent accessibility: 33-35
sequence_no = 22 #Sequence profile: 35-57


"""Define load_dataset(), get_traintestval_data(), build_model(), fit_model()"""

def load_dataset(path_to_data, peptide_no, feature_no):
    pdb_data = np.load(path_to_data) #Loaded data has dimension: 5926 proteins * 39900 features
    reshaped_data = np.reshape(pdb_data, (pdb_data.shape[0], peptide_no, feature_no)) #Reshape to get dimension: 5926 proteins * 700 aminoacids * 57 features
    return(reshaped_data)

def get_traintestval_data(reshaped_data):
    train_set = reshaped_data[0:5430, :, 0:(aminoacid_no + secondarystructure_no)]
    X_train = train_set[:,:, :aminoacid_no]
    Y_train = train_set[:,:, aminoacid_no:(aminoacid_no + secondarystructure_no)]
    
    test_set = reshaped_data[5435:5690, :, 0:(aminoacid_no + secondarystructure_no)]
    X_test = test_set[:,:, :aminoacid_no]
    Y_test = test_set[:,:, aminoacid_no:(aminoacid_no + secondarystructure_no)]
    
    validation_set = reshaped_data[5690:5926, :, 0:(aminoacid_no + secondarystructure_no)]
    X_val = validation_set[:,:, :aminoacid_no]
    Y_val = validation_set[:,:, aminoacid_no:(aminoacid_no + secondarystructure_no)]
    return(X_train, Y_train, X_test, Y_test, X_val, Y_val)
    
    
def build_model(peptide_no, aminoacid_no, secondarystructure_no):
    #Create a similar to model "Stacked LSTM for sequence classification"
    #Link: https://keras.io/getting-started/sequential-model-guide/#sequence-classification-with-lstm
    
    ###You can add more LSTM 
    model = Sequential()
    model.add(LSTM(units = 16, ###Dimension of output: Positive integer
                   activation = 'tanh', ###Default activation function: 'tanh'. Other function: 'linear'
                   recurrent_activation ='hard_sigmoid', ###Default recurrent activation function: 'hard_sigmoid'. Other function: linear
                   recurrent_dropout = 0.5, ###Dropout rate for each recurrence: between 0 and 1
                   return_sequences=True,
                   input_shape =(peptide_no, aminoacid_no))) #The first layer of a network ALWAYS contains input_shape, which is dimension of X_train/X_test
#   
#    model.add(LSTM(units = 32, ###
#                   activation = 'tanh', ###
#                   recurrent_activation ='hard_sigmoid', ###
#                   recurrent_dropout = 0.5, ###
#                   return_sequences=True))
     
    model.add(LSTM(units = secondarystructure_no,
                   activation = 'tanh', ###
                   recurrent_activation ='hard_sigmoid', ###
                   recurrent_dropout = 0.5, ###
                   return_sequences=True))
    
    model.compile(loss='binary_crossentropy',
                optimizer='adam', ###
                metrics=['accuracy'])
    
    return model

def fit_model(built_model, x_train, y_train, x_val, y_val, epochs, batch_size):
    model = built_model
    #Define callbacks
    early_stopping = EarlyStopping(monitor='val_loss', patience=3)
    tensorboard = TensorBoard(log_dir='./logs', histogram_freq=10, write_graph=True, write_grads=False, write_images=True, embeddings_freq=0, embeddings_layer_names=None, embeddings_metadata=None, embeddings_data=None, update_freq='epoch')
    filepath="weights.{epoch:02d}-{val_loss:.2f}.hdf5"
    checkpoint = ModelCheckpoint(filepath, monitor='val_acc', verbose=0, save_weights_only=False, save_best_only=False, mode='max')
    
    #Fit model
    history = model.fit(x = x_train, y = y_train, 
              validation_data = (x_val, y_val),
              epochs = epochs, 
              batch_size = batch_size, 
              callbacks=[early_stopping, tensorboard, checkpoint]) ###if you don't need a callback, just delete it from this list
    return(history)
    
if __name__ == '__main__':
    #Load dataset
    pbd = load_dataset(path_to_data = "data/cullpdb+profile_5926.npy", peptide_no = peptide_no, feature_no = feature_no)
    
    #Divide dataset into train, test, validation sets:
    whole_dataset = get_traintestval_data(pbd)
    X_train, Y_train, = whole_dataset[0], whole_dataset[1]
    X_test, Y_test = whole_dataset[2], whole_dataset[3]
    X_val, Y_val = whole_dataset[4], whole_dataset[5]
    
    #Build model with plot and summary
    LSTM_model = build_model(peptide_no, aminoacid_no, secondarystructure_no)
    plot_model(LSTM_model, show_shapes=True, to_file='LSTM.png')
    LSTM_model.summary()
    
    #Fit model
    LSTM_model_history = fit_model(built_model = LSTM_model, 
                                   x_train = X_train, y_train = Y_train, 
                                   x_val = X_val, y_val = Y_val, 
                                   epochs = 1, batch_size = 100 ###Specify your prefered epochs and batch_size here    
                                   ).history             
    print(LSTM_model_history)
    
    #Plot training & validation accuracy values
    #Link: https://keras.io/visualization/#training-history-visualization
    plt.plot(LSTM_model_history['accuracy'])
    plt.plot(LSTM_model_history['val_accuracy'])
    plt.title('Model accuracy')
    plt.ylabel('Accuracy')
    plt.xlabel('Epoch')
    plt.legend(['Train', 'Test'], loc='upper left')
    plt.show()
    
    # Plot training & validation loss values
    plt.plot(LSTM_model_history['loss'])
    plt.plot(LSTM_model_history['val_loss'])
    plt.title('Model loss')
    plt.ylabel('Loss')
    plt.xlabel('Epoch')
    plt.legend(['Train', 'Test'], loc='upper left')
    plt.show()
    
    Y_pred = LSTM_model.predict(X_test)
    print("--- %s seconds ---" % (time() - start_time))
