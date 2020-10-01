import torch 
from torch import nn
from dataset import read_frame 

def eval_performance(model,input_X):
    #evaluate performance of the model
    step = 0
    model.eval()
    tot_loss = []
    with torch.no_grad(): 
        for data in input_X:
            step+=1
            X=data['feat_vec'].reshape(data['feat_vec'].size()[0],-1)
            Y=data['wc'].reshape(data['wc'].size()[0],-1)
            output=model.forward(X)
            loss=model.compute_loss(output,Y)
            tot_loss.append(loss)
    tot_loss=torch.cat(tot_loss)
    return tot_loss

# make prediction on the MD frame
def make_prediction(model,path,n_frame,atype_dict,savedir='predicted'):
    cal_inv = lambda x,y,translation : (x*y).sum(-1) + translation
    model.eval()
    for i in range(n_frame):
        f_name=path+str(i+1)+'.ft'
        o_name=path+str(i+1)+'.wc'
        X_coor_, X_type_,Y_corr_,_,trans_vector,trans_mat,water_real,water_mol = read_frame(f_name,o_name,atype_dict)
        #make prediction
        feat_vect = torch.tensor(X_coor_,dtype=torch.float)
        X = feat_vect.reshape(feat_vect.size()[0],-1)
        with torch.no_grad(): 
            output=model.forward(X)
        #write to a file
        tot_atoms = X.size()[0]
        f_name =savedir + '/'+str(i+1)+'.xyz'
        print("writing file: ",f_name)
        with open(f_name,'w') as outfile:
            outfile.write('{} \n'.format(tot_atoms))
            for index in range(tot_atoms):
                t_vector = trans_vector[index,0]
                outfile.write('Mol {0:6d}, Tranalation vector: {1:12.6f} {2:12.6f} {3:12.6f} \n'.format(index+1,
                                                                        t_vector[0],t_vector[1],t_vector[2]))
                Y_predicted = output[index].reshape(-1,3)
                # wrire cooridnate of the water model and its WC center
                wat_r = water_real[index]
                outfile.write('O  {0:12.6f} {1:12.6f} {2:12.6f} \n'.format(wat_r[0,0],wat_r[0,1],wat_r[0,2]))
                outfile.write('H1 {0:12.6f} {1:12.6f} {2:12.6f} \n'.format(wat_r[1,0],wat_r[1,1],wat_r[1,2]))
                outfile.write('H2 {0:12.6f} {1:12.6f} {2:12.6f} \n'.format(wat_r[2,0],wat_r[2,1],wat_r[2,2]))
                for i in range(4):
                    WC_pred = cal_inv(Y_predicted[i],trans_mat[index],trans_vector[index])[0].numpy()
                    outfile.write('W  {0:12.6f} {1:12.6f} {2:12.6f} \n'.format(WC_pred[0],WC_pred[1],WC_pred[2]))
