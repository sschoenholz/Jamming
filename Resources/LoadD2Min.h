#ifndef LOAD_D2MIN

#define LOAD_D2MIN

#include "std_include.h"


void LoadD2Min(string Filename, Eigen::VectorXd &D2Min)
{
	ifstream D2Min_In(Filename.c_str());

	int i;

	//while(!D2Min_In.eof())
	for(int t = 0 ; t < D2Min.rows() ; t++)
	{
		D2Min_In >> i;
		D2Min_In >> D2Min(i);	
	}
	D2Min_In.close();
}


#endif
