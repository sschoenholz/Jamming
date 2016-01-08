#include <iostream>
#include "State/StaticState.h"
#include "Potentials/Potentials.h"
//#include "Potentials/RegisterPotentials.h"
#include "Boundaries/Boxes.h"
#include "Resources/Exception.h"
#include "Resources/Settings.h"
#include "Computers/StaticComputer.h"
#include "Computers/WallComputer.h"
#include "Computers/MatrixInterface.h"
#include "Minimization/minimizer.h"
#include "Database/Database.h"
#include "Computers/SoftSpotComputer.h"
#include "Computers/BarrierHeight.h"
#include "Resources/OutputLAMMPS.h"

using namespace LiuJamming;
#define DIM 2
using namespace std;


int main(int argc, char* argv[])
{
	//Disable harsh errors.
	NcError nc_err(NcError::silent_nonfatal);

	//the number of particles
	int number = 1024;
	double phi = 0.9;
	int seed = 0;
	double maxStrain = 0.5;
	double deltaStrain = 0.001;
	int strainOutput = 1;

	//get the inputs
	for(int i = 1 ; i < argc ; i++)
	{
		if(strcmp(argv[i],"-N")==0)
			number = atoi(argv[i+1]);
		if(strcmp(argv[i],"-P")==0)
			phi = atof(argv[i+1]);
		if(strcmp(argv[i],"-S")==0)
			seed = atoi(argv[i+1]);
		if(strcmp(argv[i],"-sM")==0)
			maxStrain = atof(argv[i+1]);
		if(strcmp(argv[i],"-sD")==0)
			deltaStrain = atof(argv[i+1]);
		if(strcmp(argv[i],"-sO")==0)
			strainOutput = atoi(argv[i+1]);
		if(strcmp(argv[i],"-h")==0)
		{
			cout << "This program generates states in strain intervals ds.\n";
			cout << "Commands:\n\t-N\tNumber of Particles\n\t-P\tPacking Fraction\n\t-S\tSeed\n\t-sM\tMaximum Strain\n\t-sD\tChange in Strain per Step.\n";
			return 0;
		}
	}

	//create a database to output the system
	char filename[256];
	sprintf(filename,"/data0/home/schsam/BaseCode_Amy/data/strain_data/strain_N%i_S%i.nc",number,seed);
	CStaticDatabase<DIM> db(number,filename,NcFile::Replace);
        
	//create the system
	CStaticState<DIM> System(number);

	//sets the particle radii to be bidisperse 50:50 with size ratio 1:1.4
	System.SetRadiiBi();	

	//sets the packing fraction to the desired packing fraction
	System.SetPackingFraction(phi);

	System.SetPotentialHertzian();
	
	//randomize the positions to an infinite temperature configuration
	System.RandomizePositions(seed);

	//create a computer (an object that computes various properties of the system)
	CStaticComputer<DIM> Computer(System);

	//create a minimizer (an object that will minimize the system to its inherent structure)
	CSimpleMinimizer<DIM> Minimizer(Computer);

	cout << "Beginning strain in steps of size " << deltaStrain << " to a maximum strain of " << maxStrain << endl;
	cout << "step\tstrain\tstress\tenergy\n";
	cout << "--------------------------------------------------\n";
	int step = 0 ;
	for(double strain = 0.0 ; strain < maxStrain ; strain += deltaStrain) 
	{
		Eigen::Matrix<double,DIM,DIM> box;
		System.GetBox()->GetTransformation(box);
		box(1,0) = strain*box(0,0);
		System.GetBox()->SetTransformation(box);

		//perform the minimization
		Minimizer.minimizeFIRE();

		//calculate and print out standard information about the system including energy / moduli etc
		Computer.StdPrepareSystem();
		Computer.CalculateStdData();

		cout << step << "\t" << strain << "\t" << Computer.Data.Stress(1,0) << "\t" << Computer.Data.Energy << endl;
	
		if(step % strainOutput == 0)
			db.WriteState(System,Computer.Data.Stress);
		step++;
	}
		

	//Computer.Data.Print();
	
	//write the state
	
	return 0;
}
