import torch 
from torch import nn
from torch.nn import functional as F
from torch.utils import data
from torch.utils.data import DataLoader

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

#construct a data loader
class WC_Dataset(data.Dataset):
    def __init__(self,feature_vector,atype,output_wc,t_vector,t_matrix,water_real,water_mol):
        self.feature_vector=torch.tensor(feature_vector,dtype=torch.float)
        self.atype=torch.tensor(atype)
        self.output_wc=torch.tensor(output_wc,dtype=torch.float)
        self.t_vector = torch.tensor(t_vector,dtype=torch.float)
        self.t_matrix = torch.tensor(t_matrix,dtype=torch.float)
        self.water_real = torch.tensor(water_real,dtype=torch.float)
        self.water_mol = torch.tensor(water_mol,dtype=torch.float)
        
    def __len__(self):
         return self.feature_vector.size(0)
        
    def __getitem__(self,index):
        f_ = self.feature_vector[index]
        a_ = self.atype[index]
        w_ = self.output_wc[index]
        t_v = self.t_vector[index]
        t_m = self.t_matrix[index]
        wat_r = self.water_real[index]
        wat_m = self.water_mol[index]
        return {'feat_vec':f_,'atype':a_,'wc':w_,'t_vec':t_v,'t_mat':t_m,'wat_r':wat_r,'wat_m':wat_m}

# nn model
class mlp(torch.nn.Module):
    def __init__(self,input_dim=400,output_dim=12):
        super(mlp, self).__init__()
        self.fnn1 = nn.Linear(input_dim,512)
        self.fnn2 = nn.Linear(512,256)
        self.fnn3 = nn.Linear(256,128)
        self.fnn4 = nn.Linear(128,64)
        self.outpur = nn.Linear(64,output_dim)
        self.mseloss = torch.nn.MSELoss(reduction='none')
    def forward(self,X):
        hx0 = torch.tanh(self.fnn1(X))  #F.softplus, torch.tanh,F.relu
        hx0 = F.dropout(hx0,p=0.5)
        hx1 = torch.tanh(self.fnn2(hx0))
        hx1 = F.dropout(hx1,p=0.5)
        hx2 = torch.tanh(self.fnn3(hx1))
        hx2 = F.dropout(hx2,p=0.5)
        hx3 = torch.tanh(self.fnn4(hx2))
        hx3 = F.dropout(hx3,p=0.5)
        y = self.outpur(hx3)
        return y
    def compute_loss(self,predicted,target):
        loss = self.mseloss(input=predicted,target=target)
        return loss

#save nn model
def save(model, m_name='wc_predicter.pt'):
    torch.save(model.state_dict(), m_name)

#load saved model
def resume(m_name, model):
    model.load_state_dict(torch.load(m_name,map_location=device))
