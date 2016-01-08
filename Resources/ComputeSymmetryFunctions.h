#ifndef COMPUTE_SYMMETRY_FUNCTIONS

#define COMPUTE_SYMMETRY_FUNCTIONS

#include <list>
#include <vector>
#include "../Resources/std_include.h"
#include "../State/StaticState.h"
//#include "../Resources/Settings.h"

namespace LiuJamming {

template <int Dim>
void ComputeSymmetryFunctions(CStaticState<Dim> &State,CSimpleGrid<Dim> *Grid,CSimpleGrid<Dim> *Angular_Grid, Settings &settings, list<int> Particles, list<vector<double> > &Symmetry)
{
	//Load the symmetry function
	vector<double> Mu = settings["Radial_Mu"];
	double Sigma = settings["Radial_Sigma"];
	vector<double> Eta = settings["Angular_Eta"];
	vector<double> Zeta = settings["Zeta"];
	vector<double> Lambda = settings["Lambda"]; 

	double Cutoff = settings["Cutoff"];
	double AngularCutoff = settings["Angular_Cutoff"];

	int Number_A = settings["Number_A"];

	Eigen::Matrix<double,Dim,1> Displacement_JK;

	for(list<int>::iterator it = Particles.begin() ; it!=Particles.end() ; it++)
	{
	//	cout << "Computing symmetry functions for particle " << *it << endl;
		
		vector<double> Symmetry_Functions;

		//compute radial functions
		for(vector<double>::iterator mu_it = Mu.begin() ; mu_it!= Mu.end() ; mu_it++)
		{	
			double G1XA = 0.0;
			double G1XB = 0.0;

			const list<CSimpleNeighbor<Dim> > &neighbors = Grid->GetNeighbors(*it);

			for(typename list<CSimpleNeighbor<Dim> >::const_iterator neigh_it = neighbors.begin() ; neigh_it != neighbors.end() ; neigh_it++)
			{
				if((*neigh_it).j < Number_A)
				{
					double norm = (*neigh_it).Displacement.norm();
					G1XA += exp(-1.0/2.0/Sigma/Sigma*(norm - *mu_it)*(norm - *mu_it));
				}else{
					double norm = (*neigh_it).Displacement.norm();
					G1XB += exp(-1.0/2.0/Sigma/Sigma*(norm - *mu_it)*(norm - *mu_it));
				}
			}
			Symmetry_Functions.push_back(G1XA);
			Symmetry_Functions.push_back(G1XB);
		}

		//compute angular functions
		for(int a_it = 0 ; a_it < Eta.size() ; a_it++)
		{
	//		cout << "Computing angular function " << a_it << endl;
			double G2XAA = 0.0;
			double G2XAB = 0.0;
			double G2XBB = 0.0;

			double eta = Eta[a_it];
			double zeta = Zeta[a_it];
			double lambda = Lambda[a_it];

			list<CSimpleNeighbor<Dim> > neighbors = Angular_Grid->GetNeighbors(*it);
			
			for(typename list<CSimpleNeighbor<Dim> >::iterator neighbor_J = neighbors.begin() ; neighbor_J != neighbors.end() ; neighbor_J++)
			{
				for(typename list<CSimpleNeighbor<Dim> >::iterator neighbor_K = neighbors.begin() ; neighbor_K != neighbors.end() ; neighbor_K++)
             		        {
					if((*neighbor_K).j != (*neighbor_J).j)
					{
						State.GetDisplacement((*neighbor_J).j,(*neighbor_K).j,Displacement_JK);
									
						double norm_JK = Displacement_JK.squaredNorm();
						if(norm_JK < AngularCutoff*AngularCutoff)
						{
							norm_JK = sqrt(norm_JK);
							double norm_IJ = (*neighbor_J).Displacement.norm();
                                                	double norm_IK = (*neighbor_K).Displacement.norm();
	
							double ctheta = (*neighbor_J).Displacement.dot((*neighbor_K).Displacement)/norm_IJ/norm_IK;											
							if((*neighbor_J).j < Number_A)
							{
								/*
								if((*neighbor_K).j < Number_A)
									G2XAA += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
								else
									G2XAB += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
								*/
								
								if((*neighbor_K).j< Number_A)
                                                                        G2XAA += 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/AngularCutoff)+1)*(cos(PI*norm_IK/AngularCutoff)+1)*(cos(PI*norm_JK/AngularCutoff)+1);
                                                                else
                                                                        G2XAB += 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/AngularCutoff)+1)*(cos(PI*norm_IK/AngularCutoff)+1)*(cos(PI*norm_JK/AngularCutoff)+1);
							}else{
								if((*neighbor_K).j > Number_A)
									G2XBB += 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/AngularCutoff)+1)*(cos(PI*norm_IK/AngularCutoff)+1)*(cos(PI*norm_JK/AngularCutoff)+1);
								/*
								if((*neighbor_K).j < Number_A)
                                                                        G2XAB += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
                                                                else
                                                      	                G2XBB += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
								*/
							}
						}							
					}
				}
			}
			
			Symmetry_Functions.push_back(G2XAA);
			Symmetry_Functions.push_back(G2XAB);
			Symmetry_Functions.push_back(G2XBB);
		}	
		
		Symmetry.push_back(Symmetry_Functions);
	}
}


template <int Dim>
void ComputeSymmetryFunctions(CStaticState<Dim> &State,list<list<CSimpleNeighbor<Dim> > > &NeighborList, Settings &settings, list<int> Particles, list<vector<double> > &Symmetry)
{
	//Load the symmetry function
	vector<double> Mu = settings["Radial_Mu"];
	double Sigma = settings["Radial_Sigma"];
	vector<double> Eta = settings["Angular_Eta"];
	vector<double> Zeta = settings["Zeta"];
	vector<double> Lambda = settings["Lambda"]; 

	double Cutoff = settings["Cutoff"];
	double AngularCutoff = settings["Angular_Cutoff"];

	int Number_A = settings["Number_A"];

	Eigen::Matrix<double,Dim,1> Displacement_JK;

	typename list<list<CSimpleNeighbor<Dim> > >::iterator neighbor_it = NeighborList.begin();

	for(list<int>::iterator it = Particles.begin() ; it!=Particles.end() ; it++)
	{
		cout << "Computing symmetry functions for particle " << *it << endl;
		
		vector<double> Symmetry_Functions;

		//compute radial functions
		for(vector<double>::iterator mu_it = Mu.begin() ; mu_it!= Mu.end() ; mu_it++)
		{	
			double G1XA = 0.0;
			double G1XB = 0.0;

			list<CSimpleNeighbor<Dim> > &neighbors = *neighbor_it;

			for(typename list<CSimpleNeighbor<Dim> >::iterator neigh_it = neighbors.begin() ; neigh_it != neighbors.end() ; neigh_it++)
			{
				if((*neigh_it).j < Number_A)
				{
					double norm = (*neigh_it).Displacement.norm();
					G1XA += exp(-1.0/2.0/Sigma/Sigma*(norm - *mu_it)*(norm - *mu_it));
				}else{
					double norm = (*neigh_it).Displacement.norm();
					G1XB += exp(-1.0/2.0/Sigma/Sigma*(norm - *mu_it)*(norm - *mu_it));
				}
			}
			Symmetry_Functions.push_back(G1XA);
			Symmetry_Functions.push_back(G1XB);
		}

		//compute angular functions
		for(int a_it = 0 ; a_it < Eta.size() ; a_it++)
		{
			cout << "Computing angular function " << a_it << endl;
			double G2XAA = 0.0;
			double G2XAB = 0.0;
			double G2XBB = 0.0;

			double eta = Eta[a_it];
			double zeta = Zeta[a_it];
			double lambda = Lambda[a_it];

			list<CSimpleNeighbor<Dim> > &neighbors = *neighbor_it;	
		
			for(typename list<CSimpleNeighbor<Dim> >::iterator neighbor_J = neighbors.begin() ; neighbor_J != neighbors.end() ; neighbor_J++)
			{
				for(typename list<CSimpleNeighbor<Dim> >::iterator neighbor_K = neighbors.begin() ; neighbor_K != neighbors.end() ; neighbor_K++)
             		        {
					if((*neighbor_K).j != (*neighbor_J).j)
					{
						State.GetDisplacement((*neighbor_J).j,(*neighbor_K).j,Displacement_JK);
									
						double norm_IJ = (*neighbor_J).Displacement.norm();
						double norm_IK = (*neighbor_K).Displacement.norm();
						double norm_JK = Displacement_JK.norm();
						if(norm_JK < AngularCutoff && norm_IJ < AngularCutoff && norm_IK < AngularCutoff)
						{
							double ctheta = (*neighbor_J).Displacement.dot((*neighbor_K).Displacement)/norm_IJ/norm_IK;											
							if((*neighbor_J).j < Number_A)
							{
								/*
								if((*neighbor_K).j < Number_A)
									G2XAA += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
								else
									G2XAB += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
								*/
								if((*neighbor_K).j< Number_A)
                                                                        G2XAA += 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/AngularCutoff)+1)*(cos(PI*norm_IK/AngularCutoff)+1)*(cos(PI*norm_JK/AngularCutoff)+1);
                                                                else
                                                                        G2XAB += 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/AngularCutoff)+1)*(cos(PI*norm_IK/AngularCutoff)+1)*(cos(PI*norm_JK/AngularCutoff)+1);
                                                        }else{
                                                                if((*neighbor_K).j > Number_A)
                                                                        G2XBB += 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/AngularCutoff)+1)*(cos(PI*norm_IK/AngularCutoff)+1)*(cos(PI*norm_JK/AngularCutoff)+1);/*
							}else{
								if((*neighbor_K).j < Number_A)
                                                                        G2XAB += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
                                                                else
                                                      	                G2XBB += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
							*/
							}
						}							
					}
				}
			}
			
			Symmetry_Functions.push_back(G2XAA);
			Symmetry_Functions.push_back(G2XAB);
			Symmetry_Functions.push_back(G2XBB);
		}	
		*neighbor_it++;
		Symmetry.push_back(Symmetry_Functions);
	}
}


}
#endif
