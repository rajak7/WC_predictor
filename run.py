import sys
import torch 
from torch import nn
from torch.nn import functional as F
from torch import optim
import numpy as np
from dataset import read_frame 
from model import mlp,save,resume
from train import train
from test import make_prediction

run_type = sys.argv[1]

if run_type not in ['train','test']:
    raise Exception('Invalid run type. only train and test allowed')
else:
    print("run type: ",run_type)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

#parameters
input_dim = 22*4
output_dim = 12
n_run = 25
model_name = 'wc_predicter.pt'
train_path = 'data/Train_data/'
train_frame = 200
test_path = 'data/Test_data/'
predict_path = 'data/Test_ML_MD/'
predict_frame = 20
test_frame = 50
atype_dict={'O':0,'H1':1,'H2':2,'WC0':3,'WC1':4,'WC2':5,'WC3':6}
atype_inv={0:'O',1:'H1',2:'H2'}

#create model
model = mlp(input_dim=122*4,output_dim=12).to(device)

if run_type == 'train':
    train(model,atype_dict,n_run,model_name,train_path,train_frame,test_path,test_frame)
else:
    resume(model_name, model)
    model.eval()
    make_prediction(model,predict_path,predict_frame,atype_dict,savedir='predicted_temp')

