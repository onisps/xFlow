#include <omp.h>
#include <csignal>
//#include "easylogging++.h"
#include "XProblem.h"
#include <iostream>
#include <math.h>
//_INITIALIZE_EASYLOGGINGPP

int main(int argc, const char* argv[])
{	
    //el::Configurations conf("log.conf");
    //el::Loggers::reconfigureLogger("default", conf);
    //el::Loggers::reconfigureAllLoggers(conf);

    //LOG(INFO) << "Start calculation";

    try
    {
		
		double Tend = 12;
		double Tstep = 0.001;
		double Stiffness = 1200;
		double Radius = 0.5;
		double Lx = 6.183800;
		double Ly = 2;
		double Lz = 2;
		double Pin = 6.;
		double Pout = 0;
		double SpaceStep = 0.05;
		
        TProblem* problem;

		TProblem::GlobalTurbulenceMode=false;
		TProblem::GlobalUsePackOperator = false;
		TProblem::GlobalTurbulenceParameter=1.7;
		TProblem::GlobalCatalog = "../Res/";
		TProblem::GlobalLength = 1;
		TProblem::GlobalVelocity = 1;
		TProblem::GlobalPressure = 1;
		TProblem::GlobalSaveStep = 50;
		/*
		//in meters
		double nu1= 1.2*1e-6;
		double nu2= 3.6*1e-6;
		double rho1 = 1030;
		double rho2 = 1050;
		//*/
		//*
		//in santimeters
		double nu1= 1.2*1e-3;
		double nu2= 3.6*1e-3;
		double rho1 = 1030*1e-3;
		double rho2 = 1050*1e-3;
		//*/


		double mu1= nu1;//*rho1;
		double mu2= nu2;//*rho2;

		//problem = TCanalIBMWithElasticBoundary::CreateInstance(12, 0.01, 400, 1200, 0.13, 3, 1, 1, 6.0, 1.0, 0.025);
		// TimeEndValue, TimeStepValue, Re, Stiffness, Radius, LengthX, LengthY, LengthZ, PressureInput, PressureOutput,SpaceStepValue

		//problem = TCanalIBMWithElasticBoundary::CreateInstance(12, 0.001, mu1, mu2, rho1, rho2, 1200, 0.125, 3, 1, 1, 6.0, 1.0, 0.025);
		problem = TCanalIBMWithElasticBoundary::CreateInstance(Tend, Tstep, mu1, mu2, rho1, rho2, Stiffness, Radius, Lx, Ly, Lz, Pin, Pout, SpaceStep);	
		//TimeEndValue, TimeStepValue, mu1, mu2, rho1, rho2, Stiffness,Radius,LX,LY,LZ,PIn,POut,SpaceStepValue 
		//mu1 and rho1 - for plasma
		//mu2 and rho2 - for blood particles
		problem->Solve();


        delete problem;		
    }
    catch (const char* s)
    {
        printf("\n");
        printf("Error: ");
        printf(s);
        printf("\n");
    }

	printf("\n");
	printf("Press any key to exit...");
	getchar();

	return 0;
}
