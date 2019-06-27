# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 09:08:49 2019

@author: jimmyli2002
"""

import pandas as pd
import numpy as np
import glob
import re

from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split

import torch
import torch.utils.data
import torch.nn as nn

#DEBUG
import pdb

"""
SECTION
Name: Definitions
Description: Global variables and hyperparameters
"""

# Device configuration
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# Fix random seed for reproducibility
seed = 7
np.random.seed(seed)

num_epochs = 100
num_classes = 2
batch_size = 32
learning_rate = 0.001

column_names = [
        "mol_name",
        "res_name",
        "res_num",
        "atom_x",
        "atom_y",
        "atom_z",
        "occupancy",
        "temp_fac",
        "label"
        ]

labels = [
        "MAN",
        "Other"
        ]

# Hyper parameters to adjust
segment_time_size = 200 # Sliding window size, must include important atoms
time_step = 100
num_features = 4
hidden_size = 60
num_lstm_layers = 2
emb_dropout= 0.2
label_threshold = 0.20

"""
SECTION
Name: Pre-processing
Description: Pre-process the data and load the information
"""

list_pdbID = []
pdbFileMatch = re.compile(r"\\(.*).pdb")
for file in glob.glob("./pdbdata/*"):
    list_pdbID.append(pdbFileMatch.findall(file)[0])
    
thead = pd.read_csv("./training_data/" + str(list_pdbID[0]) + "_data.csv", nrows=5) # just read in a few lines to get the column headers from first file
dtypes = dict(zip(thead.columns.values, ["str", "str", "int32", "float32", "float32", "float32", "float32", "float32", "str"]))   # datatypes as given by the data page

data_convoluted = []
labels = []

# Need to create label encoder which takes in all the information

molNameEncoder = LabelEncoder()

atomNames = []

for pdbID in list_pdbID:
    
    data = pd.read_csv("./training_data/" + str(pdbID) + "_data.csv", header=0, nrows=100000, dtype=dtypes, names=column_names).dropna()
    
    for i in range(len(data["mol_name"].values)):
        
        if data["mol_name"].values[i][0] not in atomNames:
            
            atomNames.append(data["mol_name"].values[i][0])
            
molNameEncoder.fit(atomNames)

print("Molecule name encoder has been created and fitted")

for pdbID in list_pdbID:
    
    # For each pdb CSV, we need to load it into the data_convoluted array
    
    data = pd.read_csv("./training_data/" + str(pdbID) + "_data.csv", header=0, nrows=100000, dtype=dtypes, names=column_names).dropna()

    # Slide a "SEGMENT_TIME_SIZE" wide window with a step size of "TIME_STEP"
    for i in range(0, len(data) - segment_time_size, time_step):
        
        completeMolName = data["mol_name"].values[i: i+ segment_time_size]
        molName = list(map(lambda mol: mol[0], completeMolName))
        molName = molNameEncoder.transform(molName) # Label encode the information
        
        x = data["atom_x"].values[i: i + segment_time_size]
        y = data["atom_y"].values[i: i + segment_time_size]
        z = data["atom_z"].values[i: i + segment_time_size]
        
        data_convoluted.append([molName, x, y, z])
    
        # Label is the most common label which is not other and passes threshold
        
        numUnique = np.unique(data["label"][i: i + segment_time_size], return_counts=True)
        combinedInfoTupples = list(zip(numUnique[0], numUnique[1])) # Convert the information into tupples
        sortedInfo = sorted(combinedInfoTupples, key=lambda x: (-x[1])) # Sort in descending order
        
        label = "Other"
        
        for i in range(len(sortedInfo)):
            
            if sortedInfo[i][0] == "Other": continue # If the label is other than ignore this loop
        
            if sortedInfo[i][1] / segment_time_size > label_threshold:
                
                label = sortedInfo[i][0]
        
            break
        
        labels.append(label)

# Convert to numpy
data_convoluted = np.asarray(data_convoluted).transpose(0, 2, 1)

# Encode labels (cross entropy loss takes labels, not one-hot [mathematically same])
encoded_labels = LabelEncoder().fit_transform(labels)

# Store the embedding data
# Create emb_dim for mol_name
vocab_size = int(data["mol_name"].nunique())
emb_dims = [(vocab_size, min(50, (vocab_size + 1) // 2))]

"""
SECTION
Name: Dataloader and Model definition
"""

class proteinDataset(object):
    
    def __init__(self, features, labels, transform=None):
        
        self.features = features
        self.labels = labels
        self.transform = transform
        
    def __len__(self):
        
        return len(self.features)
    
    def __getitem__(self, index):
        
        features = self.features[i]
        label = self.labels[i]
        
        if self.transform is not None:
            features = self.transform(features)
            
        return torch.Tensor(features), label
    
# Recurrent neural network (many-to-one)
class RNNWithEmbed(nn.Module):
    
    def __init__(self, input_size, emb_dims, hidden_size, num_lstm_layers, num_classes):
        
        super(RNNWithEmbed, self).__init__()
        
        # Make a list of embedding layers per categorical feature [emb_dims is tupple (vocab, emb_dim)]
        self.emb_layers = nn.ModuleList([nn.Embedding(x, y) for x, y in emb_dims])
        
        # Dropout layer for embedding layers
        self.emb_dropout_layer = nn.Dropout(emb_dropout)
        
        # Begin the LSTM layers
        self.hidden_size = hidden_size
        self.num_lstm_layers = num_lstm_layers
        
        self.lstm = nn.LSTM(input_size, hidden_size, num_lstm_layers, batch_first=True)
        self.fc = nn.Linear(hidden_size, num_classes)
        
    def forward(self, categorical_data, continuous_data):
        
        # Categorical features are put in their own embedding layers which are concat at w/ continuous features
        
        embedData = [emb_layer(categorical_data[i, :]) for i, emb_layer in enumerate(self.emb_layers)]
        embedData = torch.cat(embedData, 1)
        embedData = self.emb_dropout_layer(embedData)
        pdb.set_trace()
        x = torch.cat([embedData, continuous_data], 1)
        
        # Set initial hidden and cell states 
        h0 = torch.zeros(self.num_lstm_layers, x.size(0), self.hidden_size).to(device) 
        c0 = torch.zeros(self.num_lstm_layers, x.size(0), self.hidden_size).to(device)
        
        # Forward propagate LSTM
        out, _ = self.lstm(x, (h0, c0)) # Output format: tensor (batch_size, seeq_len, hidden_size)
        
        # Generate meaning from last time step
        out = self.fc(out[:, -1, :])
        
        return out

"""
SECTION
Name: Model creation and training
"""

model = RNNWithEmbed(num_features, emb_dims, hidden_size, num_lstm_layers, num_classes)

# Split data into training and validation sets

X_train, X_test, y_train, y_test = train_test_split(data_convoluted, encoded_labels, test_size=0.1, random_state=seed)
print("X train size: ", len(X_train))
print("X test size: ", len(X_test))
print("y train size: ", len(y_train))
print("y test size: ", len(y_test))

train_dataset = proteinDataset(X_train, y_train, transform=None)
train_loader = torch.utils.data.DataLoader(dataset=train_dataset, batch_size=batch_size, shuffle=True)

test_dataset = proteinDataset(X_test, y_test, transform=None)
test_loader = torch.utils.data.DataLoader(dataset=test_dataset, batch_size=batch_size, shuffle=True)

# Loss and optimizer definition

loss_fn = nn.CrossEntropyLoss()
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

# Train the model

total_step = len(train_loader)

for epoch in range(num_epochs):
    
    for i, (features, labels) in enumerate(train_loader):
        
        features = features.to(device)
        labels = labels.type(torch.LongTensor)
        labels = labels.to(device)
        
        # Forward pass
        outputs = model(features[:, :, 0].long(), features[:, :, 1:].float())
        loss = loss_fn(outputs. labels)
        
        # Backward pass and optimize
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        if (i + 1) % 5 == 0:
            
            print ("Epoch [{}/{}], Step [{}/{}], Loss: {:.4f}" .format(epoch+1, num_epochs, i+1, total_step, loss.item()))
