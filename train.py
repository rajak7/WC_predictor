import torch 
from torch import nn
from torch import optim
import numpy as np
from dataset import construct_dataset 
from model import WC_Dataset,save
from torch.utils import data
from torch.utils.data import DataLoader
from test import eval_performance

def train(model,atype_dict,n_run,model_name,train_path,train_frame,test_path,test_frame):
    
    #read training and test data set
    X_coor_trj,X_type_trj,Y_corr_trj,trans_vector_trj,trans_mat_trj,water_real_trj,water_mol_trj = construct_dataset(train_frame,train_path,atype_dict)
    Xtest_coor_trj,Xtest_type_trj,Ytest_corr_trj,trans_test_vector_trj,trans_test_mat_trj,water_test_real_trj,water_test_mol_trj = construct_dataset(test_frame,test_path,atype_dict)
    print("Training data info:")
    print(X_coor_trj.shape,X_type_trj.shape,Y_corr_trj.shape)
    print(trans_vector_trj.shape,trans_mat_trj.shape,water_real_trj.shape,water_mol_trj.shape)
    print("Test data info:")
    print(Xtest_coor_trj.shape,Xtest_type_trj.shape,Ytest_corr_trj.shape)
    print(trans_test_vector_trj.shape,trans_test_mat_trj.shape,water_test_real_trj.shape,water_test_mol_trj.shape)
    #training set
    train_data=WC_Dataset(X_coor_trj,X_type_trj,Y_corr_trj,trans_vector_trj,trans_mat_trj,
                      water_real_trj,water_mol_trj)
    train_X = DataLoader(train_data, batch_size=64, shuffle=True)
    #test set
    test_data=WC_Dataset(Xtest_coor_trj,Xtest_type_trj,Ytest_corr_trj,trans_test_vector_trj,trans_test_mat_trj,
                    water_test_real_trj,water_test_mol_trj)
    test_X = DataLoader(train_data, batch_size=64, shuffle=True)

    # construct optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)
    # strat training
    step = 0
    tot_loss_trj = []
    for epoch in range(n_run): 
        model.train()
        for data in train_X:
            step+=1
            optimizer.zero_grad()
            X=data['feat_vec'].reshape(data['feat_vec'].size()[0],-1)
            Y=data['wc'].reshape(data['wc'].size()[0],-1)
            output=model.forward(X)
            loss=model.compute_loss(output,Y)
            loss=loss.mean()
            tot_loss_trj.append(loss.item())
            if step % 100 == 0:
               print("epoch: ",epoch,"step: ",step,"loss: ",loss.item())
            loss.backward()
            optimizer.step()

    # save trained model
    save(model, model_name)
    print("Model Evaluation")
    train_loss = eval_performance(model,train_X)
    test_loss = eval_performance(model,test_X)
    train_loss=train_loss.reshape(-1).numpy()
    test_loss=test_loss.reshape(-1).numpy()
    print("Train MSE and STD:",train_loss.mean(),train_loss.std())
    print("Test MSE and STD:",test_loss.mean(),test_loss.std())
    print("Max Error in Training:",np.max(np.sqrt(train_loss)))
    print("Max Error in Test:",np.max(np.sqrt(test_loss)))
    print("Min Error in Training:",np.min(np.sqrt(train_loss)))
    print("Min Error in Test:",np.min(np.sqrt(test_loss)))
    print("Train RMSE: ",np.sqrt(train_loss.mean()))
    print("Test RMSE: ",np.sqrt(test_loss.mean()))
    return
