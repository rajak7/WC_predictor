import numpy as np 

def read_frame(f_name,o_name,atype_dict,n_feat=122):
    X_coor_trj=[]
    X_type_trj=[]
    #read feature vector of molecules
    with open(f_name) as in_file:
        val = in_file.readline()  #oxygen 1 coordinae
        _ = in_file.readline()   #inverse matrix
        _ = in_file.readline()   #inverse matrix
        _ = in_file.readline()   #inverse matrix
        _ = in_file.readline() # atom info
        count = 0
        mol=0
        feat_coor=[]
        feat_type=[]
        for val in in_file:
            if count == n_feat:
                mol+=1
                #print("reading molecule",mol,val.strip())
                val = in_file.readline()
                val = in_file.readline()
                val = in_file.readline()
                val = in_file.readline()
                count = 0
                X_coor_trj.append(feat_coor)
                X_type_trj.append(feat_type)
                feat_coor=[]
                feat_type=[]
            else:
                count+=1
                val=val.strip().split()
                feat_coor.append([float(val[1]),float(val[2]),float(val[3]),float(val[4])])
                feat_type.append(atype_dict[val[0]])
        if count != 0:
            X_coor_trj.append(feat_coor)
            X_type_trj.append(feat_type)  
    X_coor_trj=np.asarray(X_coor_trj)
    X_type_trj=np.asarray(X_type_trj)
    #read wc and atomic cooridante of the molecules
    trans_vector_trj =[]
    trans_mat_trj =[]
    water_real_trj =[]
    water_mol_trj = []
    WC_corr_trj=[]
    WC_type_trj=[]
    
    with open(o_name) as in_file:
        for val in in_file:
            val = val.strip().split()
            #print('reading mol',val)
            translation_vector = []
            translation_vector.append([float(val[6]),float(val[7]),float(val[8])])
            trans_vector_trj.append(translation_vector)

            transformaton_matrix = []
            val = in_file.readline().strip().split()   #matrix - 00
            transformaton_matrix.append([float(val[0]),float(val[1]),float(val[2])])
            val = in_file.readline().strip().split()   #matrix - 01
            transformaton_matrix.append([float(val[0]),float(val[1]),float(val[2])])
            val = in_file.readline().strip().split()   #matrix - 02
            transformaton_matrix.append([float(val[0]),float(val[1]),float(val[2])])
            trans_mat_trj.append(transformaton_matrix)

            water_real_coor = []
            val = in_file.readline().strip().split()   #O read
            water_real_coor.append([float(val[1]),float(val[2]),float(val[3])])
            val = in_file.readline().strip().split()   #H1 real
            water_real_coor.append([float(val[1]),float(val[2]),float(val[3])])
            val = in_file.readline().strip().split()   #H2 real
            water_real_coor.append([float(val[1]),float(val[2]),float(val[3])])
            water_real_trj.append(water_real_coor)

            _ = in_file.readline() 
            _ = in_file.readline()
            _ = in_file.readline()

            water_mol_coor = []
            val = in_file.readline().strip().split()   #O mol
            water_mol_coor.append([float(val[1]),float(val[2]),float(val[3])])
            val = in_file.readline().strip().split()   #H1 mol
            water_mol_coor.append([float(val[1]),float(val[2]),float(val[3])])
            val = in_file.readline().strip().split()   #H2 mol
            water_mol_coor.append([float(val[1]),float(val[2]),float(val[3])])
            water_mol_trj.append(water_mol_coor)

            #reading WC 
            w_corr=[]
            w_type=[]
            val = in_file.readline().strip().split() 
            w_corr.append([float(val[1]),float(val[2]),float(val[3])])
            w_type.append(atype_dict[val[0]])

            val = in_file.readline().strip().split()
            w_corr.append([float(val[1]),float(val[2]),float(val[3])])
            w_type.append(atype_dict[val[0]])

            val = in_file.readline().strip().split()
            w_corr.append([float(val[1]),float(val[2]),float(val[3])])
            w_type.append(atype_dict[val[0]])

            val = in_file.readline().strip().split()
            w_corr.append([float(val[1]),float(val[2]),float(val[3])])
            w_type.append(atype_dict[val[0]])

            WC_corr_trj.append(w_corr)
            WC_type_trj.append(w_type)

    WC_corr_trj=np.asarray(WC_corr_trj)
    WC_type_trj=np.asarray(WC_type_trj)
    trans_vector_trj = np.asarray(trans_vector_trj)
    trans_mat_trj = np.asarray(trans_mat_trj)
    water_real_trj = np.asarray(water_real_trj)
    water_mol_trj = np.asarray(water_mol_trj)
    return X_coor_trj,X_type_trj,WC_corr_trj,WC_type_trj,trans_vector_trj,trans_mat_trj,water_real_trj,water_mol_trj
        

def construct_dataset(nframe,path,atype_dict):
    for i in range(nframe):
        f_name=path+str(i+1)+'.ft'
        o_name=path+str(i+1)+'.wc'
        X_coor_, X_type_,Y_corr_,_,trans_vector,trans_mat,water_real,water_mol = read_frame(f_name,o_name,atype_dict)
        if i==0:
            X_coor_trj=X_coor_
            X_type_trj=X_type_
            Y_corr_trj=Y_corr_
            trans_vector_trj = trans_vector
            trans_mat_trj = trans_mat
            water_real_trj = water_real
            water_mol_trj = water_mol
        else:
            X_coor_trj=np.concatenate((X_coor_trj, X_coor_), axis=0)
            X_type_trj=np.concatenate((X_type_trj, X_type_), axis=0)
            Y_corr_trj=np.concatenate((Y_corr_trj, Y_corr_), axis=0)
            trans_vector_trj=np.concatenate((trans_vector_trj, trans_vector), axis=0)
            trans_mat_trj=np.concatenate((trans_mat_trj, trans_mat), axis=0)
            water_real_trj=np.concatenate((water_real_trj, water_real), axis=0)
            water_mol_trj=np.concatenate((water_mol_trj, water_mol), axis=0)
    return X_coor_trj,X_type_trj,Y_corr_trj,trans_vector_trj,trans_mat_trj,water_real_trj,water_mol_trj