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
	int Number = 100;
	double phi = 0.9;

	//create a database to output the system
	CStaticDatabase<DIM> db(Number,"/data0/home/schsam/BaseCode_Amy/sdata/test_out.nc",NcFile::Replace);
        
	//create the system
	CStaticState<DIM> System(Number);

	//sets the particle radii to be bidisperse 50:50 with size ratio 1:1.4
	System.SetRadiiBi();	

	//sets the packing fraction to the desired packing fraction
	System.SetPackingFraction(phi);

	//randomize the positions to an infinite temperature configuration
	System.RandomizePositions();

	//create a computer (an object that computes various properties of the system)
	CStaticComputer<DIM> Computer(System);

	//create a minimizer (an object that will minimize the system to its inherent structure)
	CSimpleMinimizer<DIM> Minimizer(Computer);

	//perform the minimization
	Minimizer.minimizeFIRE();

	//calculate and print out standard information about the system including energy / moduli etc
	Computer.StdPrepareSystem();
	Computer.CalculateStdData();
	Computer.Data.Print();
	
	//write the state
	db.WriteState(System);
	
	return 0;
}
