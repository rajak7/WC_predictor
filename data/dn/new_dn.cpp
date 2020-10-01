// Last updated 09/27/2020
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>
#include <tuple>
#include <functional>

#include <cstdlib>
#include <cmath>
#include <cassert>

#include "Eigen/Eigen"

#define MAX_NEIGHBOR_MOLS 40

struct MDFrame
{
		std::string filename; 

    int natoms;
    double lattice[6];

    std::vector<double> x, y, z;
    std::vector<double> vx, vy, vz;
		std::vector<int> mol_id;  // molecular Id
    std::vector<std::string> name;

    double apply_pbc(double x, double lattice)
    {
        if(x>=0.5*lattice) x-=lattice;
        if(x<-0.5*lattice) x+=lattice;
        return x;
    };

		void print()
		{
    	std::cout << " # of Atoms : " << natoms << std::endl;
	    std::cout << " Lattice Consts : " <<
              lattice[0] << " " << lattice[1] << " " << lattice[2] << " " << 
              lattice[3] << " " << lattice[4] << " " << lattice[5] << std::endl;
		}
};

MDFrame read_single_mdframe(std::ifstream &in, std::string _filename="NA")
{
    MDFrame mdframe;
    std::string str;

    std::getline(in,str);
		mdframe.filename = _filename;
    mdframe.natoms = std::atoi(str.c_str());

    double dummy;
    std::stringstream ss;
    std::getline(in,str);
    ss << str;
		ss >> str >> mdframe.lattice[0] >> dummy >> dummy >>
			dummy >> mdframe.lattice[1] >> dummy >>
			dummy >> dummy >> mdframe.lattice[2];
    //ss >> mdframe.lattice[0] >> mdframe.lattice[1] >> mdframe.lattice[2];

		mdframe.name.resize(mdframe.natoms);
		mdframe.x.resize(mdframe.natoms);
		mdframe.y.resize(mdframe.natoms);
		mdframe.z.resize(mdframe.natoms);
		mdframe.mol_id.resize(mdframe.natoms);
		std::cout << mdframe.natoms << std::endl;

    for (int i=0; i<mdframe.natoms; i++)
    {
        std::string name;
        float x,y,z,vx,vy,vz;
				int id; 

        std::stringstream ss;
        std::getline(in,str);
        ss << str;

        ss >> name >> x >> y >> z >> vx >> vy >> vz;
/*
        ss >> name >> x >> y >> z >> id; 
        mdframe.name[id-1] = name;
        mdframe.x[id-1] = x;
        mdframe.y[id-1] = y;
        mdframe.z[id-1] = z;
        mdframe.mol_id[id-1] = id;
*/

        ss >> name >> x >> y >> z >> id; 
        mdframe.name[i] = name;
        mdframe.x[i] = x;
        mdframe.y[i] = y;
        mdframe.z[i] = z;
        mdframe.mol_id[i] = id;
        //std::cout << mdframe.name[i] << " " << mdframe.x[i] << " " << mdframe.y[i] << " " << mdframe.z[i] << std::endl; 
        //std::cout << name << " " <<x << " " << y << " " << z << " " << vx << " " << vy << " " << vz << std::endl;
    }

    return mdframe;
};

typedef std::tuple<double, int, std::string> tuple_t;
typedef std::vector<std::vector<tuple_t>> nbrlist_t;

struct NeighborList
{
	int max_neighbors;

	nbrlist_t nbrlist, nbrlist_wc;

	NeighborList(MDFrame & mdframe, int _max=MAX_NEIGHBOR_MOLS) : max_neighbors(_max)
	{
		nbrlist.resize(mdframe.natoms);	
		nbrlist_wc.resize(mdframe.natoms);	

		for (int i=0; i<mdframe.natoms; i++)
		{
			//std::cout << mdframe.name[i] << " " << mdframe.x[i] << " " << mdframe.y[i] << " " << mdframe.z[i] << std::endl; 

			// construct nbrlist from only O
			if (mdframe.name[i].find("O") == std::string::npos) continue; 

			for (int j=0; j<mdframe.natoms; j++)
			{
				double dx = mdframe.x[j] - mdframe.x[i];
				double dy = mdframe.y[j] - mdframe.y[i];
				double dz = mdframe.z[j] - mdframe.z[i];

				dx = mdframe.apply_pbc(dx, mdframe.lattice[0]);
				dy = mdframe.apply_pbc(dy, mdframe.lattice[1]);
				dz = mdframe.apply_pbc(dz, mdframe.lattice[2]);

				double dr = sqrt(dx*dx + dy*dy + dz*dz); 
				if (dr == 0.0) continue;

				// use sparate nbrlist for WCs
				if (mdframe.name[j].find("X") == std::string::npos) 
					nbrlist[i].push_back(std::make_tuple(dr, j, mdframe.name[j]));

			}

			// use sparate nbrlist for WCs
			for (int j=i+1; j<mdframe.natoms; j++)
			{
				if (mdframe.name[j].find("X") != std::string::npos) 
				{
					double dx = mdframe.x[j] - mdframe.x[i];
					double dy = mdframe.y[j] - mdframe.y[i];
					double dz = mdframe.z[j] - mdframe.z[i];

					dx = mdframe.apply_pbc(dx, mdframe.lattice[0]);
					dy = mdframe.apply_pbc(dy, mdframe.lattice[1]);
					dz = mdframe.apply_pbc(dz, mdframe.lattice[2]);

					double dr = sqrt(dx*dx + dy*dy + dz*dz); 
					nbrlist_wc[i].push_back(std::make_tuple(dr, j, mdframe.name[j]));
				}
			}


		}

		for (std::vector<tuple_t> & n : nbrlist) 
		{
			// erase non-atom data
			n.erase(std::remove_if(n.begin(), n.end(), [](auto const &t) {return std::get<0>(t) == 0.0;}), n.end());

			// sort by distance
			std::sort(begin(n), end(n), [](auto const & t1, auto const & t2) { return std::get<0>(t1) < std::get<0>(t2);} );

/*
			// cap at max_neighbors
			if(n.size() > max_neighbors) n.resize(max_neighbors);
*/
		}

		
		for (auto & nbr_wc : nbrlist_wc) 
		{
			if(nbr_wc.size()>0)
			{
				assert(nbr_wc.size()>=4);
				nbr_wc.resize(4);
			}
		}

	}

	void print()
	{
		for (int i=0; i<nbrlist.size(); i++)
		{
			std::cout << std::endl << i << "th neighbr list size : " << nbrlist[i].size() << std::endl << std::endl;
			for (const auto & l : nbrlist[i])
				std::cout << i << " " << std::get<0>(l) << " " << std::get<1>(l) << " " << std::get<2>(l) << std::endl; 
				std::cout << std::endl; 
		}
	}

};


Eigen::Matrix3d get_transformation_matrix(MDFrame & mdframe, nbrlist_t const& nbrlist, const int i)
{
	// make the oxygen as reference 
	Eigen::Vector3d ref(mdframe.x[i], mdframe.y[i], mdframe.z[i]);

	auto const & h1_atom = nbrlist[i][0];

	// closest atom must be H1, but this condition might be relaxed for non-water system later
	//std::cout << i << " std::get<2>(h1_atom)  " << std::get<2>(h1_atom) << " " << std::get<0>(h1_atom) << std::endl;
	assert(std::get<2>(h1_atom).find("H1") != std::string::npos); 

	int const & idx = std::get<1>(h1_atom); 

	Eigen::Vector3d r_h1(mdframe.x[idx], mdframe.y[idx], mdframe.z[idx]);

	r_h1 -= ref;
	for (int ia=0; ia<3; ia++) r_h1(ia) = mdframe.apply_pbc(r_h1(ia), mdframe.lattice[ia]);
	auto const r1 = r_h1.normalized();

	//std::cout << "r1 " << r1(0) << " " << r1(1) << " " << r1(2) << " " << r1.norm() << std::endl;

	auto const & h2_atom = nbrlist[i][1];

	// 2nd closest atom must be H2, but this condition might be relaxed for non-water system later
	//std::cout << i << " std::get<2>(h2_atom)  " << std::get<2>(h2_atom) << " " << std::get<0>(h2_atom) << std::endl;
	assert(std::get<2>(h2_atom).find("H2") != std::string::npos); 

	int const & idx2 = std::get<1>(h2_atom); 
	Eigen::Vector3d r_h2(mdframe.x[idx2], mdframe.y[idx2], mdframe.z[idx2]);

	r_h2 -= ref;
	for (int ia=0; ia<3; ia++) r_h2(ia) = mdframe.apply_pbc(r_h2(ia), mdframe.lattice[ia]);
	auto r2 = r_h2.normalized();

	std::cout << "H1 " << std::get<2>(nbrlist[i][0]) << " " << mdframe.x[idx] << " " << mdframe.y[idx] << " " << mdframe.z[idx] << std::endl;
	std::cout << "H2 " << std::get<2>(nbrlist[i][1]) << " " << mdframe.x[idx2] << " " << mdframe.y[idx2] << " " << mdframe.z[idx2] << std::endl;

	//std::cout << "r2 " << r2(0) << " " << r2(1) << " " << r2(2) << " " << r2.norm() << std::endl;

	auto n1 = r1.cross(r2); n1.normalize();
	//std::cout << "n1 " << n1(0) << " " << n1(1) << " " << n1(2) << " " << n1.norm() << std::endl;

	auto n2 = r1.cross(n1); n2.normalize();
	//std::cout << "n2 " << n2(0) << " " << n2(1) << " " << n2(2) << " " << n2.norm() << std::endl;

	Eigen::Matrix3d M;

	// reverse rotation matrix. from molecule to lab coordinate.
	M << r1(0), n1(0), n2(0),
			r1(1), n1(1), n2(1),
			r1(2), n1(2), n2(2);

	return M;
}

void get_dp_features(MDFrame & mdframe, NeighborList const& Nbr)
{
	std::string labframe = mdframe.filename + "_labframe.xyz";
	std::string molframe = mdframe.filename + "_molframe.xyz"; 
	std::string nbratoms = mdframe.filename + "_nbratoms.xyz"; 
	std::string feature_file = mdframe.filename +".ft";
	std::string groundtruth_file = mdframe.filename +".wc";

	std::ofstream xyz_lab(labframe);
	std::ofstream xyz_mol(molframe);
	std::ofstream xyz_nbr(nbratoms);
	std::ofstream ft_out(feature_file);
	std::ofstream wc_out(groundtruth_file);

	nbrlist_t const& nbrlist = Nbr.nbrlist;
	nbrlist_t const& nbrlist_wc = Nbr.nbrlist_wc;

	for (int i=0; i<mdframe.natoms; i++)
	{
		// return if not oxygen 
		if (mdframe.name[i].find("O") == std::string::npos) continue;

		// rotaiton matrix. from lab to molecule coordinate.
		Eigen::Matrix3d M = get_transformation_matrix(mdframe, nbrlist, i);
		Eigen::Matrix3d Minv = M.inverse();

		std::cout << "-------------------" << std::endl;
		std::cout << M << std::endl;
		std::cout << "-------------------" << std::endl;
		std::cout << Minv << std::endl;
		//std::cout << Minv*M << std::endl;
		//assert(Minv*M == Eigen::Matrix3d::Identity());
		std::cout << "-------------------" << std::endl;

		std::vector<int>  nbridx_o, nbridx_h1, nbridx_h2;

		auto const & idx_h1 = std::get<1>(nbrlist[i][0]);
		auto const & idx_h2 = std::get<1>(nbrlist[i][1]);

		nbridx_o.push_back(i);
		nbridx_h1.push_back(idx_h1);
		nbridx_h2.push_back(idx_h2);

		for (const auto & nbr : nbrlist[i])
		{
			int const & idx_o = std::get<1>(nbr); 

			if (mdframe.name[idx_o].find("O") != std::string::npos) 
			{
				nbridx_o.push_back(idx_o);
				auto const & nbr_o = nbrlist[idx_o];

				int const & idx_h1 = std::get<1>(nbr_o[0]); 
				int const & idx_h2 = std::get<1>(nbr_o[1]); 

				nbridx_h1.push_back(idx_h1);
				nbridx_h2.push_back(idx_h2);

/*
				std::cout << "idex_[o,h1,h2] " << 
					nbridx_o.size() << " " << nbridx_o.size() << " " << nbridx_o.size() << "\t\t" << 
					mdframe.name[idx_o] << " " << idx_o << "\t" << 
					mdframe.name[idx_h1] << " " << idx_h1 << "\t" << 
					mdframe.name[idx_h2] << " " << idx_h2 << std::endl;
*/
			}
		}

		std::cout << "num_nbr_atoms " << nbrlist[i].size() << ", num_wc_atoms " << nbrlist_wc[i].size() << std::endl;
		std::cout << "O " << nbridx_o.size() << ", H1 " << nbridx_h1.size() << ", H2 " << nbridx_h2.size() << std::endl;
		if(nbridx_o.size() > MAX_NEIGHBOR_MOLS+1) nbridx_o.resize(MAX_NEIGHBOR_MOLS+1);
		if(nbridx_h1.size() > MAX_NEIGHBOR_MOLS+1) nbridx_h1.resize(MAX_NEIGHBOR_MOLS+1);
		if(nbridx_h2.size() > MAX_NEIGHBOR_MOLS+1) nbridx_h2.resize(MAX_NEIGHBOR_MOLS+1);

		std::vector<std::vector<int>> nbridx_lists = {nbridx_o, nbridx_h1, nbridx_h2};

		// ref: center oxygen position
		Eigen::Vector3d ref(mdframe.x[i], mdframe.y[i], mdframe.z[i]);


//FIXME the boilerplate
//=========================================================================================================
		wc_out << "mol_id " << mdframe.mol_id[i] << " , oxygen in labframe " 
			<< " " << ref(0) << " " << ref(1) << " " << ref(2) << std::endl;
		wc_out << M << std::endl; 

		// real coordinates
		for(int k=i; k<i+3; k++)
			wc_out <<mdframe.name[i]<<"(real)" << " "<< mdframe.x[k] <<  " " << mdframe.y[k] <<  " " << mdframe.z[k] << std::endl;

		// lab coordnates, shifted to the center O
		for(int k=i; k<i+3; k++)
		{
			Eigen::Vector3d r_lab(mdframe.x[k], mdframe.y[k], mdframe.z[k]);
			r_lab -= ref; 
			wc_out <<mdframe.name[i]<<"(lab) " << " "<< r_lab[0] <<  " " << r_lab[1] <<  " " << r_lab[2] << std::endl;
		}
		// mol coordinates, shifted to the center O, then rotated
		for(int k=i; k<i+3; k++)
		{
			Eigen::Vector3d r_lab(mdframe.x[k], mdframe.y[k], mdframe.z[k]);
			r_lab -= ref; 
			for (int ia=0; ia<3; ia++) r_lab(ia) = mdframe.apply_pbc(r_lab(ia), mdframe.lattice[ia]);
			auto const r_mol = Minv*r_lab;
			wc_out <<mdframe.name[i]<<"(mol) " << " "<< r_mol[0] <<  " " << r_mol[1] <<  " " << r_mol[2] << std::endl;
		}

		if(nbrlist_wc[i].size()>0)
		{
			int counter = 0; 
			for(auto l : nbrlist_wc[i])
			{
				int iw = std::get<1>(l);
				wc_out << "WC" << std::to_string(counter++) << " " << mdframe.x[iw] <<  " " << mdframe.y[iw] <<  " " << mdframe.z[iw] << std::endl;
			}
		} else {
			for(int iw=0; iw<4; iw++) wc_out << "WC" << std::to_string(iw) << " 0.0 0.0 0.0" << std::endl;
		}
//=========================================================================================================
		xyz_lab << nbridx_o.size() + nbridx_h1.size() + nbridx_h2.size() + nbrlist_wc[i].size() + 4 << std::endl;
		xyz_lab << "atom # " << i << " lattice " << 
			mdframe.lattice[0] << " " << mdframe.lattice[1] << " " << mdframe.lattice[2] << std::endl;
		xyz_lab << "N 0.0 0.0 0.0  -1" << std::endl;
		xyz_lab << "S " << M(0,0) << " " << M(1,0) << " " << M(2,0) << " " << -1 << std::endl;
		xyz_lab << "S " << M(0,1) << " " << M(1,1) << " " << M(2,1) << " " << -1 << std::endl;
		xyz_lab << "S " << M(0,2) << " " << M(1,2) << " " << M(2,2) << " " << -1 << std::endl;

		for(auto l : nbrlist_wc[i])
		{
			int iw = std::get<1>(l);
			Eigen::Vector3d wc(mdframe.x[iw],mdframe.y[iw],mdframe.z[iw]);
			wc -= ref;
			xyz_lab << "WC " << wc(0) << " " << wc(1) << " " << wc(2) << " " << iw << std::endl;
		}

		for (const auto nbridx : nbridx_lists)
		{
			for (const int j : nbridx)
			{
				Eigen::Vector3d r_lab(mdframe.x[j], mdframe.y[j], mdframe.z[j]);
				r_lab -= ref; 
				xyz_lab << mdframe.name[j] << " " << r_lab(0) << " " << r_lab(1) << " " << r_lab(2) << " " << j << std::endl;
			}
		}
//=========================================================================================================
		xyz_mol << nbridx_o.size() + nbridx_h1.size() + nbridx_h2.size() + nbrlist_wc[i].size() + 4 << std::endl;
		xyz_mol << "atom # " << i << " lattice " << 
			mdframe.lattice[0] << " " << mdframe.lattice[1] << " " << mdframe.lattice[2] << std::endl;
		xyz_mol << "N 0.0 0.0 0.0 -1" << std::endl; 
		xyz_mol << "S 1.0 0.0 0.0 -1" << std::endl;
		xyz_mol << "S 0.0 1.0 0.0 -1" << std::endl;
		xyz_mol << "S 0.0 0.0 1.0 -1" << std::endl;

		for(auto l : nbrlist_wc[i])
		{
			int iw = std::get<1>(l);
			Eigen::Vector3d wc(mdframe.x[iw],mdframe.y[iw],mdframe.z[iw]);
			wc -= ref;
			for (int ia=0; ia<3; ia++) wc(ia) = mdframe.apply_pbc(wc(ia), mdframe.lattice[ia]);
			auto const wc_mol = Minv*wc;
			xyz_mol<< "WC " << wc_mol(0) << " " << wc_mol(1) << " " << wc_mol(2) << " " << iw << std::endl;
		}

		for (const auto nbridx : nbridx_lists)
		{
			for (const int j : nbridx)
			{
				Eigen::Vector3d r_lab(mdframe.x[j], mdframe.y[j], mdframe.z[j]);
				r_lab -= ref; 
				for (int ia=0; ia<3; ia++) r_lab(ia) = mdframe.apply_pbc(r_lab(ia), mdframe.lattice[ia]);
				auto const r_mol = Minv*r_lab;
				xyz_mol<< mdframe.name[j] << " " << r_mol(0) << " " << r_mol(1) << " " << r_mol(2) << " " << j << std::endl;
			}
		}
//=========================================================================================================
		xyz_nbr << nbridx_o.size() + nbridx_h1.size() + nbridx_h2.size()  + 1 << std::endl;
		xyz_nbr << "atom # " << i << std::endl;
		xyz_nbr << "N" << " " << ref(0) << " " << ref(1) << " " << ref(2) << std::endl;

		for (const auto nbridx : nbridx_lists)
		{
			for (const int j : nbridx)
			{
				Eigen::Vector3d r_lab(mdframe.x[j], mdframe.y[j], mdframe.z[j]);
				xyz_nbr << mdframe.name[j] << " " << r_lab(0) << " " << r_lab(1) << " " << r_lab(2) << std::endl;
			}
		}
//=========================================================================================================
		ft_out << "mol_id " << mdframe.mol_id[i] << " , oxygen in labframe " 
			<< " " << ref(0) << " " << ref(1) << " " << ref(2) << std::endl;
		ft_out << M << std::endl;
		ft_out << "O " << nbridx_o.size()-1 << ", H1 " << nbridx_h1.size() << ", H2 "  << nbridx_h2.size() << std::endl;

		for (const auto nbridx : nbridx_lists)
		{
			for (const int j : nbridx)
			{
				Eigen::Vector3d r_lab(mdframe.x[j], mdframe.y[j], mdframe.z[j]);
				r_lab -= ref; 
				if(r_lab.norm() == 0) continue; // must be the center O.
				for (int ia=0; ia<3; ia++) r_lab(ia) = mdframe.apply_pbc(r_lab(ia), mdframe.lattice[ia]);

				auto r_mol = Minv*r_lab;

				const double r_mol_inv = 1.0/r_mol.norm();
				const double r_mol_inv_sq = r_mol_inv*r_mol_inv; 

				const auto ft_vec = r_mol*r_mol_inv_sq; 
				ft_out << mdframe.name[j] << " " << r_mol_inv << " " << ft_vec(0) << " " << ft_vec(1) << " " << ft_vec(2) << std::endl;
			}
		}
//=========================================================================================================
	}

	xyz_lab.close();
	xyz_mol.close();
	ft_out.close();
	wc_out.close();
}

int main(int argc, char *argv[])
{
    std::vector<MDFrame> mdframes;

		if ( std::string("-h").compare(std::string(argv[1])) == 0)
		{
			std::cout << "help message here" << std::endl;
			exit(0);
		}

		std::string filename = argv[1];
    std::ifstream ifile(filename);

		std::cout << "loading " << filename << std::endl;

		MDFrame single_frame = read_single_mdframe(ifile, filename); 
		single_frame.print();

		get_dp_features(single_frame, NeighborList(single_frame));

		std::cout << "\nsuccessfully finished\n";
		return 0;
}
