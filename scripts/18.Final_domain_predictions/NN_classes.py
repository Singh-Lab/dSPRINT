import numpy as np
import pandas as pd

# Neural Net imports
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.optim.lr_scheduler import LambdaLR, ReduceLROnPlateau
from sklearn.model_selection import RepeatedStratifiedKFold

#Checking which device NN should run on
curr_device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
#====================================================================================================================#

#Define the Network
#Tutorial for Neural Net Architecture: https://pytorch.org/tutorials/beginner/blitz/neural_networks_tutorial.html

#Utilize batch normalization, as explained here: https://www.youtube.com/watch?v=fv1Luwd-LOI&index=69&list=PLBAGcD3siRDguyYYzhVwZ3tLvOyyG5k6K

# The network without early stopping 
class Net(nn.Module):
    def __init__(self, dropout_parameter = 0.5, hidden_units_1 = 200, 
                 hidden_units_2 = 400, batch_size = 75, 
                 learning_rate = 1e-5, beta = 0.9, 
                 weight_decay = 1e-4, epoch_count = 15, weight="balanced", input_size = 750, rseed=0):
        
        torch.manual_seed(rseed)
        super(Net, self).__init__()   
        self.input = nn.Linear(input_size, hidden_units_1) 
        self.hidden1 = nn.Linear(hidden_units_1, hidden_units_2)
        self.hidden1_bn = nn.BatchNorm1d(hidden_units_2)
        self.hidden2 = nn.Linear(hidden_units_2, hidden_units_2)
        self.hidden2_bn = nn.BatchNorm1d(hidden_units_2)
        self.hidden3 = nn.Linear(hidden_units_2, hidden_units_1)
        self.hidden3_bn = nn.BatchNorm1d(hidden_units_1)
        self.dropout = nn.Dropout(p = dropout_parameter)
        self.output = nn.Linear(hidden_units_1,2)
        self.learning_rate = learning_rate
        self.beta = beta
        self.batch_size = batch_size
        self.weight_decay = weight_decay
        self.epoch_count = epoch_count
        self.weight = weight
        self.rseed=rseed
        
    def forward(self, x):
        sf = nn.Softmax()
        x = F.rrelu(self.input(x))
        x = self.dropout(F.rrelu(self.hidden1_bn(self.hidden1(x))))
        x = self.dropout(F.rrelu(self.hidden2_bn(self.hidden2(x))))
        x = self.dropout(F.rrelu(self.hidden3_bn(self.hidden3(x))))
        x = self.output(x)
        return x
    
    def fit(self, train_valid_data, train_valid_labels, weight):
        # set in training mode
        self.train()
        
        # set random seed
        torch.manual_seed(self.rseed)
          
        trainset = pd.concat([train_valid_data,train_valid_labels],axis=1)
        trainset = shuffle(trainset, random_state=self.rseed)

        train_valid_data = trainset.iloc[:,:trainset.shape[1]-1]
        train_valid_labels = trainset.iloc[:,trainset.shape[1]-1]

        # create loss function
        loss = nn.BCEWithLogitsLoss(weight = weight)
        # mini-batching
        batch_size = self.batch_size
        
        BETA_2 = 0.999        
        no_batch_minus_1 = train_valid_data.shape[0] / batch_size 

        skf_2 = RepeatedStratifiedKFold(n_splits=no_batch_minus_1,n_repeats = self.epoch_count,random_state=0)

        # create adam optimizer for Phase 2
        optimizer_2 = optim.Adam(self.parameters(), lr=self.learning_rate,betas = (self.beta,BETA_2), 
                                 weight_decay = self.weight_decay)
        
        lambda1 = lambda epoch_count: 0.99 ** epoch_count 
        scheduler = LambdaLR(optimizer_2, lr_lambda=lambda1)
        
        count = 0
        epoch_count = 0
        
        for train,test in skf_2.split(train_valid_data,train_valid_labels):
            data = train_valid_data.iloc[test,:]
            data = torch.Tensor(data.values.astype(np.float32)).to(device=curr_device)
            # forward pass          
            output = self.forward(data)
            output.data = output.data.view(data.shape[0],2)

            labels = train_valid_labels[test]
            labels = labels.astype(int)
            labels = torch.Tensor(np.eye(2)[labels]).to(device=curr_device)
            labels = torch.autograd.Variable(labels, requires_grad = False)

            # zero the gradient buffers
            optimizer_2.zero_grad()
            # compute loss and gradients
            loss_output = loss(output,labels)
            loss_output.backward()
            # Does the update
            optimizer_2.step()
            
            count = count + 1

            # Early Stopping
            if count == no_batch_minus_1 + 1:
                count = 0
                epoch_count = epoch_count + 1
                scheduler.step()
            
    #prediction probabilities array
    def predict_proba(self, X_test):
        self.eval()
        #forward pass
        test = torch.Tensor(X_test.values.astype(np.float32)).to(device=curr_device)
        output = self.forward(test)
        sf = nn.Softmax()
        probs = sf(output.data)
        probs_list = []
        for i in range(len(probs)):
            probs_list.append(probs[i][1].item())          
        return probs_list
#====================================================================================================================#

# define the network with early stopping
class Net_tune(nn.Module):
    def __init__(self, dropout_parameter = 0.5, hidden_units_1 = 200, 
                 hidden_units_2 = 400, batch_size = 75, 
                 learning_rate = 1e-5, beta = 0.9, weight_decay = 1e-4, input_size = 750):
        torch.manual_seed(0)
        super(Net_tune, self).__init__()
        self.input = nn.Linear(input_size, hidden_units_1) # read input size from the .shape of data table
        self.hidden1 = nn.Linear(hidden_units_1, hidden_units_2)
        self.hidden1_bn = nn.BatchNorm1d(hidden_units_2)
        self.hidden2 = nn.Linear(hidden_units_2, hidden_units_2)
        self.hidden2_bn = nn.BatchNorm1d(hidden_units_2)
        self.hidden3 = nn.Linear(hidden_units_2, hidden_units_1)
        self.hidden3_bn = nn.BatchNorm1d(hidden_units_1)
        self.dropout = nn.Dropout(p = dropout_parameter)
        self.output = nn.Linear(hidden_units_1,2)
        self.learning_rate = learning_rate
        self.beta = beta
        self.batch_size = batch_size
        self.weight_decay = weight_decay
  
    def forward(self, x):
        x = F.rrelu(self.input(x))
        x = self.dropout(F.rrelu(self.hidden1_bn(self.hidden1(x))))
        x = self.dropout(F.rrelu(self.hidden2_bn(self.hidden2(x))))
        x = self.dropout(F.rrelu(self.hidden3_bn(self.hidden3(x))))
        x = self.output(x)
        return x
    
    def fit(self, X_train, y_train_label, X_valid, y_valid, weight):
            # sets model in training mode because batch normalization behavior in training and testing modes are different
            self.train()
            # set random seed for weights and biases
            torch.manual_seed(0)

            # dataset
            dataset = pd.concat([X_train,y_train_label],axis=1)
            dataset = shuffle(dataset, random_state = 0)

            X_train = dataset.iloc[:,:dataset.shape[1]-1]
            y_train_label = dataset.iloc[:,dataset.shape[1]-1]

            # create loss function
            loss = nn.BCEWithLogitsLoss(weight = weight)
            # mini-batching
            batch_size = self.batch_size

            BETA_2 = 0.999
            TOTAL_EPOCHS_TRAINED = 10**4

            # create adam optimizer for Phase 1
            optimizer_1 = optim.Adam(self.parameters(), lr=self.learning_rate,betas=(self.beta,BETA_2), 
                                     weight_decay = self.weight_decay)

            lambda1 = lambda epoch_count: 0.99 ** epoch_count 
            scheduler = LambdaLR(optimizer_1, lr_lambda=lambda1)
            no_batch_minus_1 = X_train.shape[0] / batch_size 

            # Repeated Stratified K Fold to ensure positives are evenly distributed across batches
            skf_1 = RepeatedStratifiedKFold(n_splits=no_batch_minus_1,n_repeats=TOTAL_EPOCHS_TRAINED,random_state=0)

            INITIAL_PATIENCE = 100
            count = 0
            epoch_count = 0
            max_auprc = 0
            ideal_epoch_count = 0 
            patience = INITIAL_PATIENCE
            #print "initial patience = "+str(patience)
            patience_j = 0

            for train,test in skf_1.split(X_train,y_train_label):
                data = X_train.iloc[test,:]
                data = torch.Tensor(data.values.astype(np.float32)).to(device=curr_device)
                 # forward pass
                output = self.forward(data)
                output.data = output.data.view(data.shape[0],2)

                labels = y_train_label[test]
                labels = labels.astype(int)
                labels = torch.Tensor(np.eye(2)[labels]).to(device=curr_device)
                labels = torch.autograd.Variable(labels, requires_grad = False)

                # zero the gradient buffers
                optimizer_1.zero_grad()
                # compute loss and gradients
                loss_output = loss(output,labels)
                loss_output.backward()
                # Does the update
                optimizer_1.step()

                count = count + 1

                # Early Stopping
                if count == no_batch_minus_1 + 1:
                    count = 0
                    epoch_count = epoch_count + 1
                    scheduler.step()
                    probs_valid = self.predict_proba(X_valid)
                    precision, recall, _ = precision_recall_curve(y_valid, probs_valid)
                    auprc = auc(recall, precision)
                    #print auprc
                    if auprc > max_auprc:
                        max_auprc = auprc
                        #print "max_auprc = "+str(max_auprc)
                        ideal_epoch_count = epoch_count
                        patience = INITIAL_PATIENCE + epoch_count
                        #print "Updating patience to "+str(patience)
                        patience_j = 0
                    else:
                        #print "patience_j="+str(patience_j)
                        patience_j = patience_j + 1 
                        if patience_j == patience: break
            
                self.train()
            return max_auprc, ideal_epoch_count

        
    #prediction probabilities array
    def predict_proba(self, X_test):
        self.eval()
        #forward pass
        test = torch.Tensor(X_test.values.astype(np.float32)).to(device=curr_device)
        output = self.forward(test)
        sf = nn.Softmax()
        probs = sf(output.data)
        return probs[:,1]
    
    def evaluate_loss(self,X_valid, y_valid, weight):
        data = torch.Tensor(X_valid.values.astype(np.float32)).to(device=curr_device)
        loss = nn.BCEWithLogitsLoss(weight = weight)
         # forward pass
        output = self.forward(data)
        output.data = output.data.view(data.shape[0],2)

        labels = y_valid
        labels = labels.astype(int)
        labels = torch.Tensor(np.eye(2)[labels]).to(device=curr_device)
        labels = torch.autograd.Variable(labels, requires_grad = False)

        # compute loss and gradients
        loss_output = loss(output,labels)

        return loss_output        
#====================================================================================================================#