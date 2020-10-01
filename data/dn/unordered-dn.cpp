// Last updated 09/10/2020
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

#define MAX_NEIGHBORS 100

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

    for (int i=0; i<mdframe.natoms; i++)
    {
        std::string name;
        float x,y,z,vx,vy,vz;
				int id; 

        std::stringstream ss;
        std::getline(in,str);
        ss << str;
        //ss >> name >> x >> y >> z >> vx >> vy >> vz;
        ss >> name >> x >> y >> z >> id; 

        //std::cout << x << " " << y << " " << z << " " << vx << " " << vy << " " << vz << std::endl;

        mdframe.name.push_back(name);
        mdframe.x.push_back(x);
        mdframe.y.push_back(y);
        mdframe.z.push_back(z);
        mdframe.mol_id.push_back(id);
        //mdframe.vx.push_back(vx);
        //mdframe.vy.push_back(vy);
        //mdframe.vz.push_back(vz);
    }

    return mdframe;
};

typedef std::tuple<double, int, std::string> tuple_t;
typedef std::vector<std::vector<tuple_t>> nbrlist_t;

struct NeighborList
{
	int max_neighbors;

	nbrlist_t nbrlist;

	NeighborList(MDFrame & mdframe, int _max=MAX_NEIGHBORS) : max_neighbors(_max)
	{
		nbrlist.resize(mdframe.natoms);	

		for (int i=0; i<nbrlist.size(); i++) 
			nbrlist[i].resize(mdframe.natoms);

		for (int i=0; i<mdframe.natoms; i++)
		{
			if (mdframe.name[i].find("X") != std::string::npos) continue; // skip WC during nbrlist construction
			for (int j=0; j<mdframe.natoms; j++)
			{
				if (mdframe.name[j].find("X") != std::string::npos) continue; // skip WC during nbrlist construction
				double dx = mdframe.x[j] - mdframe.x[i];
				double dy = mdframe.y[j] - mdframe.y[i];
				double dz = mdframe.z[j] - mdframe.z[i];

				dx = mdframe.apply_pbc(dx, mdframe.lattice[0]);
				dy = mdframe.apply_pbc(dy, mdframe.lattice[1]);
				dz = mdframe.apply_pbc(dz, mdframe.lattice[2]);

				double dr = sqrt(dx*dx + dy*dy + dz*dz); 
				if (dr == 0.0) continue;

				nbrlist[i].push_back(std::make_tuple(dr, j, mdframe.name[j]));
			}
		}

		for (std::vector<tuple_t> & n : nbrlist) 
		{
			// erase non-atom data
			n.erase(std::remove_if(n.begin(), n.end(), [](auto const &t) {return std::get<0>(t) == 0.0;}), n.end());

			// sort by distance
			std::sort(begin(n), end(n), [](auto const & t1, auto const & t2) { return std::get<0>(t1) < std::get<0>(t2);} );

			// cap at max_neighbors
			if(n.size() > max_neighbors) n.resize(max_neighbors);
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

std::vector<Eigen::Vector3d> get_wc_and_atom_positions_in_molframe(MDFrame const& mdframe, int oid, Eigen::Matrix3d const& M)
{
	std::vector<Eigen::Vector3d> wcs; 
	
	Eigen::Vector3d ref(mdframe.x[oid],mdframe.y[oid],mdframe.z[oid]);

// real :  original coordinates from MD
	for (int i=0; i<3; i++)
	{
		Eigen::Vector3d h2o_real(mdframe.x[oid+i],mdframe.y[oid+i],mdframe.z[oid+i]);
		wcs.push_back(h2o_real);
	}

// lab : oxygen-centered real coordinates
	for (int i=0; i<3; i++)
	{
		Eigen::Vector3d h2o_lab(mdframe.x[oid+i],mdframe.y[oid+i],mdframe.z[oid+i]);
		h2o_lab -= ref; 
		wcs.push_back(h2o_lab);
	}

// mol : H2O atoms in molecular coordinates
	for (int i=0; i<3; i++)
	{
		Eigen::Vector3d h2o_mol(mdframe.x[oid+i],mdframe.y[oid+i],mdframe.z[oid+i]);
		h2o_mol -= ref; 
		h2o_mol = M*h2o_mol;
		wcs.push_back(h2o_mol);
	}

// mol : WC in molecular coordinates
	for (int i=3; i<7; i++)
	{
		//std::cout << "oid " << oid << " wc_id " << oid + i << std::endl;
		Eigen::Vector3d wc(mdframe.x[oid+i],mdframe.y[oid+i],mdframe.z[oid+i]);
		wc -= ref; 
		wc = M*wc;

		wcs.push_back(wc);
	}

	return wcs;
}

void get_dp_features(MDFrame & mdframe, nbrlist_t const& nbrlist)
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

	for (int i=0; i<mdframe.natoms; i++)
	{
		// return if not oxygen 
		if (mdframe.name[i].find("O") == std::string::npos) continue;

		// make the oxygen as reference 
		Eigen::Vector3d ref(mdframe.x[i], mdframe.y[i], mdframe.z[i]);

		auto const & first = nbrlist[i][0];

		// closest atom must be H1, but this condition might be relaxed for non-water system later
		assert(std::get<2>(first).find("H1") != std::string::npos); 

		int const & idx = std::get<1>(first); 

		Eigen::Vector3d r_h1(mdframe.x[idx], mdframe.y[idx], mdframe.z[idx]);

		r_h1 -= ref;
		for (int ia=0; ia<3; ia++) r_h1(ia) = mdframe.apply_pbc(r_h1(ia), mdframe.lattice[ia]);
		auto const r1 = r_h1.normalized();

		//std::cout << "r1 " << r1(0) << " " << r1(1) << " " << r1(2) << " " << r1.norm() << std::endl;

		auto const & second = nbrlist[i][1];

		// 2nd closest atom must be H2, but this condition might be relaxed for non-water system later
		assert(std::get<2>(second).find("H2") != std::string::npos); 

		int const & idx2 = std::get<1>(second); 
		Eigen::Vector3d r_h2(mdframe.x[idx2], mdframe.y[idx2], mdframe.z[idx2]);

		r_h2 -= ref;
		for (int ia=0; ia<3; ia++) r_h2(ia) = mdframe.apply_pbc(r_h2(ia), mdframe.lattice[ia]);
		auto r2 = r_h2.normalized();

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

		// rotaiton matrix. from lab to molecule coordinate.
		Eigen::Matrix3d Minv = M.inverse();

/*
		std::cout << "-------------------" << std::endl;
		std::cout << M << std::endl;
		std::cout << Minv << std::endl;
		std::cout << Minv*M << std::endl;
		//assert(Minv*M == Eigen::Matrix3d::Identity());
		std::cout << "-------------------" << std::endl;
*/
		auto const wc_p = get_wc_and_atom_positions_in_molframe(mdframe, i, Minv);

		wc_out << "mol_id " << mdframe.mol_id[i] << " , oxygen in labframe " 
			<< " " << ref(0) << " " << ref(1) << " " << ref(2) << std::endl;
		wc_out << M << std::endl; 

		std::vector<std::string> frames = {"real","lab","mol"};
		int iw = 0; 
		for (auto str : frames)
		{
			wc_out <<"O("<<str<<") "<< wc_p[iw][0] <<  " " << wc_p[iw][1] <<  " " << wc_p[iw][2] << std::endl;
			wc_out <<"H1("<<str<<") "<< wc_p[iw+1][0] <<  " " << wc_p[iw+1][1] <<  " " << wc_p[iw+1][2] << std::endl;
			wc_out <<"H2("<<str<<") "<< wc_p[iw+2][0] <<  " " << wc_p[iw+2][1] <<  " " << wc_p[iw+2][2] << std::endl;

			iw+=3;
		}

		for (int iw = 9; iw<13; iw++) 
				wc_out << "WC" << std::to_string(iw) << " " << 
				wc_p[iw][0] <<  " " << wc_p[iw][1] <<  " " << wc_p[iw][2] << std::endl;

		xyz_lab << nbrlist[i].size() + 4 + 4 << std::endl;
		xyz_lab << "atom # " << i << " lattice " << 
			mdframe.lattice[0] << " " << mdframe.lattice[1] << " " << mdframe.lattice[2] << std::endl;
		xyz_lab << "N 0.0 0.0 0.0  -1" << std::endl;
		xyz_lab << "S " << r1(0) << " " << r1(1) << " " << r1(2) << " " << -1 << std::endl;
		xyz_lab << "S " << n1(0) << " " << n1(1) << " " << n1(2) << " " << -1 << std::endl;
		xyz_lab << "S " << n2(0) << " " << n2(1) << " " << n2(2) << " " << -1 << std::endl;

		for (int iw=i+3; iw<i+7; iw++)
		{
			Eigen::Vector3d wc(mdframe.x[iw],mdframe.y[iw],mdframe.z[iw]);
			wc -= ref;
			xyz_lab << "WC " << wc(0) << " " << wc(1) << " " << wc(2) << " " << mdframe.mol_id[i] << std::endl;
		}

		xyz_mol << nbrlist[i].size() + 4 + 4 << std::endl;
		xyz_mol << "atom # " << i << " lattice " << 
			mdframe.lattice[0] << " " << mdframe.lattice[1] << " " << mdframe.lattice[2] << std::endl;
		xyz_mol << "N 0.0 0.0 0.0 -1" << std::endl; 
		xyz_mol << "S 1.0 0.0 0.0 -1" << std::endl;
		xyz_mol << "S 0.0 1.0 0.0 -1" << std::endl;
		xyz_mol << "S 0.0 0.0 1.0 -1" << std::endl;
		for (auto const& w : wc_p)
			xyz_mol << "WC " << w(0) << " " << w(1) << " " << w(2) << " " << mdframe.mol_id[i] << std::endl;

		xyz_nbr << nbrlist[i].size() + 1 << std::endl;
		xyz_nbr << "atom # " << i << std::endl;
		xyz_nbr << "N" << " " << ref(0) << " " << ref(1) << " " << ref(2) << std::endl;

		ft_out << "mol_id " << mdframe.mol_id[i] << " , oxygen in labframe " 
			<< " " << ref(0) << " " << ref(1) << " " << ref(2) << std::endl;
		ft_out << M;
		ft_out << nbrlist[i].size() << std::endl;

			//std::cout << w(0) << " " << w(1) << " " << w(2) << " " << mdframe.mol_id[i] << std::endl;

		for (const auto & l : nbrlist[i])
		{
			int j = std::get<1>(l);

			Eigen::Vector3d r_lab(mdframe.x[j], mdframe.y[j], mdframe.z[j]);
			xyz_nbr << std::get<2>(l) << " " << r_lab(0) << " " << r_lab(1) << " " << r_lab(2) << std::endl;

			r_lab -= ref; 
			for (int ia=0; ia<3; ia++) r_lab(ia) = mdframe.apply_pbc(r_lab(ia), mdframe.lattice[ia]);

			auto r_mol = Minv*r_lab;
			auto r_mol2lab = M*r_mol;

			//std::cout << "  r_lab " << r_lab(0) << " " << r_lab(1) << " " << r_lab(2) << " " << r_lab.norm() << std::endl;
			//std::cout << "  r_mol " << r_mol(0) << " " << r_mol(1) << " " << r_mol(2) << " " << r_mol.norm() << std::endl;

/* How to revert molecule frame to lab frame
			auto mol2lab_x = M(0,0)*r_mol(0)+M(0,1)*r_mol(1)+M(0,2)*r_mol(2);
			auto mol2lab_y = M(1,0)*r_mol(0)+M(1,1)*r_mol(1)+M(1,2)*r_mol(2);
			auto mol2lab_z = M(2,0)*r_mol(0)+M(1,2)*r_mol(1)+M(2,2)*r_mol(2);

			std::cout << "  r_lab " << r_lab(0) << " " << r_lab(1) << " " << r_lab(2) << " " << r_lab.norm() << std::endl;
			std::cout << "  r_mol2lab " << mol2lab_x << " " << mol2lab_y << " " << mol2lab_z << " " << r_mol2lab.norm() << std::endl << std::endl;;
*/

			xyz_lab << std::get<2>(l) << " " << r_lab(0) << " " << r_lab(1) << " " << r_lab(2) << " " << mdframe.mol_id[j] << std::endl;
			xyz_mol << std::get<2>(l) << " " << r_mol(0) << " " << r_mol(1) << " " << r_mol(2) << " " << mdframe.mol_id[j] << std::endl;

			const double r_mol_inv = 1.0/r_mol.norm();
			const double r_mol_inv_sq = r_mol_inv*r_mol_inv; 

			const auto ft_vec = r_mol*r_mol_inv_sq; 
			ft_out << std::get<2>(l) << " " << r_mol_inv << " " << ft_vec(0) << " " << ft_vec(1) << " " << ft_vec(2) << std::endl;
		}

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
		//single_frame.print();

		get_dp_features(single_frame, NeighborList(single_frame).nbrlist);

		return 0;
}
