# Jamming

Codebase used to create and analyze amorphous / crystalline solids. The main piece of code is the CStaticState class which holds particle position, radius, potential, and boundary condition information. Particle potentials are defined by overloading the CPotential class, likewise different boundary conditions can be implemented by overloading the CBox class. The CStaticDatabase class handles output / input of states to/from netcdf files. Minimization to inherent structures is implemented by default using the FIRE and LBFGS algorithms. However, different minimizers can easily be added. Finally, information (for example energies, gradients, and hessians) can be computed using CStaticComputer. Different kinds of minimization (for example on spring networks instead of collections of particles) can be implemented by writing a new CStaticState class and then overloading CBaseComputer. This allows seamless computation of elastic constants and minimization over relatively disperate systems.

A simple example to create and quench an amorphous configuration to its inherent structure and output this structure do a database would be,
```
char filename[256];
int number = 1024;
sprintf(filename,"/data0/home/schsam/BaseCode_Amy/data/strain_data/strain_N%i_S%i.nc",number,seed);
CStaticDatabase<2> db(number,filename,NcFile::Replace);
        
//create the system
CStaticState<2> System(number);

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
	
Minimizer.minimizeFIRE();
db.WriteState(System);
```

Here is an image of such a configuration.
![](https://samschoenholz.files.wordpress.com/2015/10/jammed.png)
