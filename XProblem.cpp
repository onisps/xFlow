#pragma once
#pragma hdrstop


#include "XProblem.h"
#include "XSLAE.h"
#include "XUtils.h"
#include <string.h>
#include <csignal>
#include <ctime>
#include <cmath>
#include "XImmersedBoundary.h"

#if CNPY_VISUALIZATION
	#include"cnpy.h"
#endif

#include <omp.h>
#include <functional>
#include <stdexcept>
#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>

//#include "easylogging++.h"

#define RANGE 3
#define ISNAN(value) value != value
//#define M_PI 3.14159265358979323846
//#define M_PI_4 M_PI/4.
#define CNPY_VISUALIZATION 0

//using namespace std;

using std::vector;

//#################################################################


//Temporary step (Begin)
bool TProblem::GlobalTurbulenceMode=false;
double TProblem::GlobalTurbulenceParameter=0;
bool TProblem::GlobalUsePackOperator=true;
char* TProblem::GlobalCatalog = NULL;
double TProblem::GlobalLength=0;
double TProblem::GlobalVelocity=0;
double TProblem::GlobalPressure=0;
int TProblem::GlobalSaveStep = 0;


int iter = 0;


//Temporary step (End)


//#################################################################

//*************************************
T1dPoissonProblem::T1dPoissonProblem(double lH):fH(lH)
{
   printf("The problem is preparing data...\n");

   double h = fH;

   TSeparator *sep1 = new TUniformSeparator(0,1,h);

   int N = sep1->EndIndex;

   T1dNumberMask &mask = *( new T1dNumberMask(N+1));
   T1dVector** normals = new T1dVector*[N+1];

   for (int i=0;i<N+1;i++)
   {
         T1dVector* normal = NULL;
         mask[i]=TFictivePoint;

		 if ( (i!=0) && (i!=N) )
         {
             mask[i]=TActualPoint;
         }
         else
         {
            if ( (i==0) )
            {
				mask[i]=TDefinedBorderPoint;
				normal = new T1dVector;
				normal->X = -1;
			}
            else if ( (i==N) )
			{
				mask[i]=TDefinedBorderPoint;
				//mask[i]=TNormalBorderPoint;
				normal = new T1dVector;
				normal->X = 1;
			}
         }
		 normals[i] = normal;
    }

    fGrid = new T1dNormalGrid(sep1,mask,normals);

    delete &mask;

	for (int i=0;i<N+1;i++)
    {
			if ( (normals[i])!=NULL ) {delete normals[i];}
	}
	delete[] normals;// normals = NULL;
}
//*************************************
T1dPoissonProblem::~T1dPoissonProblem()
{
   delete fGrid;
   printf("The problem has been completed successfully.\n");
}
//*************************************
double T1dPoissonProblem::SolutionFunction(double x)
{
	return sin(x);
}
//*************************************
double T1dPoissonProblem::SolutionFunctionDX(double x)
{
	return cos(x);
}
//*************************************
double T1dPoissonProblem::RightPartFunction(double x)
{
   return sin(x);
}
//*************************************
void T1dPoissonProblem::Solve()
{
   printf("The problem is starting ...\n");

   double h = fH;
   T1dNormalGrid& grid = *fGrid;

   int N = grid.GetSeparator1().EndIndex;
   const double* x = grid.GetSeparator1().Dimension;

   const T1dNumberMask& mask = grid.Mask;

   const T1dVector* const* normals = grid.Normals;


   T1dNumberMask& spaceMask = mask.Copy();
   for (int i=0;i<N+1;i++)
   {
	  if (spaceMask[i]!=TDefinedBorderPoint)
	  {
	      spaceMask[i] = TActualPoint;
	  }
	  else
	  {
		  spaceMask[i] = TFictivePoint;
	  }
   }


   TMask** spaceMasks = new TMask*[1];
   spaceMasks[0] = &spaceMask;

   TRnRelease1dSpaceModel& spaceModel = *(new TRnRelease1dSpaceModel(1,spaceMasks));

   int spaceModelName = 5;
   TRnSpace::AddModel(spaceModelName,&spaceModel,false);

   TRnRelease1dSpace& solution = *(new TRnRelease1dSpace(spaceModelName));
   TRnRelease1dSpace& rightPart = *(new TRnRelease1dSpace(spaceModelName));

   TRnRelease1dSpace& approximation = *(new TRnRelease1dSpace(spaceModelName));
   TRnRelease1dSpace& zeroHelp = *(new TRnRelease1dSpace(spaceModelName));
   TRnRelease1dSpace& rightPartHelp = *(new TRnRelease1dSpace(spaceModelName));

   TRnRelease1dSpace& error = *(new TRnRelease1dSpace(spaceModelName));

   for (int i=0;i<N+1;i++)
   {
	   error[i] = 0;
	   solution[i] = SolutionFunction(x[i]);
	   if (mask[i] == TFictivePoint)
	   {
		  solution[i] = 0;
		  rightPart[i] = 0;
		  approximation[i] = 0;
	   }
	   else if (mask[i] == TActualPoint)
	   {
		  rightPart[i] = RightPartFunction(x[i]);
		  approximation[i] = 0;
	   }
	   else if (mask[i]==TDefinedBorderPoint)
	   {
		  rightPart[i] = 0;
		  approximation[i] = solution[i];
	   }
	   else if (mask[i]==TNormalBorderPoint)
	   {
		  rightPart[i] = (SolutionFunctionDX(x[i]))*(normals[i]->X);
		  rightPart[i] = SolutionFunction(x[i])/2 + (rightPart[i])/h;
		  approximation[i] = 0;
	   }
	   else {throw "Error in T1dPoissonProblem::Solve: Unknown border point type.";  }
   }


   TNodesGrid** grids = new TNodesGrid*[1];
   grids[0] = &grid;

   T1dLaplaceOperator& A = *(new T1dLaplaceOperator(spaceModelName,1,grids));
   delete[] grids;// grids = NULL;


   //Pack operator (Begin)
   zeroHelp.Assign(approximation);
   zeroHelp|=0;
   rightPartHelp|= ((TRnLinearOperator&)A)(zeroHelp) - rightPart;
   rightPartHelp|=(-1)*rightPartHelp;



   int L = approximation.EndIndex;
   double* helpVector = new double[L];
   for (int l=0;l<L;l++) {helpVector[l]=0;}

   T1dNumberMask& spaceMaskRn = *(new T1dNumberMask(L));
   for (int l=0;l<L;l++) {spaceMaskRn[l] = TActualPoint;}
   TMask** spaceMasksRn = new TMask*[1];
   spaceMasksRn[0] = &spaceMaskRn;
   TRnRelease1dSpaceModel& spaceModelRn = *(new TRnRelease1dSpaceModel(1,spaceMasksRn));
   int spaceModelNameRn = 55;
   TRnSpace::AddModel(spaceModelNameRn,&spaceModelRn,false);



   TRnRelease1dSpace& solutionRn = *(new TRnRelease1dSpace(spaceModelNameRn));
   TRnRelease1dSpace& rightPartRn = *(new TRnRelease1dSpace(spaceModelNameRn));
   TRnRelease1dSpace& approximationRn = *(new TRnRelease1dSpace(spaceModelNameRn));
   TRnRelease1dSpace& errorRn = *(new TRnRelease1dSpace(spaceModelNameRn));

   for (int l=1;l<L+1;l++)
   {
	   solutionRn[l-1] = ((TRnSpace&)solution)[l];

	   rightPartRn[l-1] = ((TRnSpace&)rightPartHelp)[l];
	   helpVector[l-1] = rightPartRn[l-1];

	   approximationRn[l-1] = ((TRnSpace&)approximation)[l];

       errorRn[l-1] = ((TRnSpace&)error)[l];
   }


   printf("Preparing pack operator: Starting...\n");
   bool flagTransform = true;
   T1dPackOperator* ARnPtr;
   if (flagTransform == false)
   {
	   double* nullVector = NULL;
       ARnPtr = new T1dPackOperator(spaceModelNameRn,A,&nullVector);
   }
   else
   {
	   ARnPtr = new T1dPackOperator(spaceModelNameRn,A,&helpVector);
	   for (int l=0;l<L;l++)
       {
         rightPartRn[l] = helpVector[l];
       }
   }
   T1dPackOperator& ARn = *ARnPtr;
   printf("Preparing pack operator: OK.\n");

   //Temp (Begin)
   /*
   //error|= ((TRnLinearOperator&)A)(approximation) - ((TRnSpace&)rightPart);
   //errorRn|= ((TRnLinearOperator&)ARn)(approximationRn) - ((TRnSpace&)rightPartRn);
   //printf("\n");
   //printf("%.10f",error.GetAbsMaxElement());
   //printf("\n");
   //printf("%.10f",errorRn.GetAbsMaxElement());
   //printf("\n");
   //return;

   A.ConstructMatrix();
   ARn.ConstructMatrix();

   const double* const* matrTemp = A.Matrix;
   const double* const* matrRn = ARn.Matrix;

   double** matr = new double*[L];
   for (int l=0;l<L;l++) {matr[l] = new double[L];}
   for (int i=0;i<L;i++)
   {
	   for (int j=0;j<L;j++)
	   {
		   double v = 0;
		   for  (int l=0;l<L;l++)
		   {
			  v = v + matrTemp[l][i]*matrTemp[l][j];
		   }
		   matr[i][j] = v;
	   }
   }


   for (int i=0;i<L;i++)
   {
	   for (int j=0;j<L;j++)
	   {
		   if (matr[i][j]!=matrRn[i][j]) {throw "Matrixes are not equal.";}
	   }
   }
   //throw "User stop.";
   */
   //Temp (End)

   TSLAE& slae = *(new TSLAE(ARn,approximationRn,rightPartRn));
   //slae.SolveByGauss();
   slae.SolveByConjugateGradients();
   //slae.SolveByMinimalResiduals();
   approximationRn|=slae.Result;
   delete &ARn;
   delete &slae;

   errorRn|=solutionRn-approximationRn;
   double errorNorma = errorRn.GetAbsMaxElement();

   delete &spaceMaskRn;
   delete[] spaceMasksRn;//  spaceMaskRn = NULL;
   delete &spaceModelRn;

   delete &solutionRn;
   delete &rightPartRn;
   delete &approximationRn;
   delete &errorRn;

   delete helpVector;
   //Pack operator (End)


   /*
   //Ordinary operator (Begin)
   TSLAE& slae = *(new TSLAE(A,approximation,rightPart));
   slae.SolveByConjugateGradients();
   //slae.SolveByMinimalResiduals();
   approximation|=slae.Result;
   delete &A;
   delete &slae;

   error|=solution-approximation;
   double errorNorma = error.GetAbsMaxElement();
   //Ordinary operator (End)
   */

   printf("\n");
   printf("The norma is ");
   printf("%.10f",errorNorma);
   printf("\n");

   delete &spaceMask;
   delete[] spaceMasks;//  spaceMask = NULL;
   delete &spaceModel;

   delete &solution;
   delete &rightPart;
   delete &approximation;
   delete &zeroHelp;
   delete &rightPartHelp;
   delete &error;

   TRnSpace::RemoveAllModels();
}
//*************************************



//#################################################################



//*************************************
T3dPoissonProblem::T3dPoissonProblem(double lH):fH(lH)
{

   printf("The problem is preparing data...\n");

   double h = fH;

  TSeparator* sep1F = new TUniformSeparator(-h/2,1+h/2,h);
  TSeparator* sep2F = new TUniformSeparator(-h/2,1+h/2,h);
  TSeparator* sep3F = new TUniformSeparator(-h/2,1+h/2,h);
  T3dNumberMask* maskF = new T3dNumberMask((sep1F->EndIndex)+1,(sep2F->EndIndex)+1,(sep3F->EndIndex)+1);



  TSeparator *sep1G = new TUniformSeparator(0,1,h);
  TSeparator *sep2G = new TUniformSeparator(0,1,h);
  TSeparator *sep3G = new TUniformSeparator(0,1,h);
  T3dNumberMask* maskG = new T3dNumberMask((sep1G->EndIndex)+1,(sep2G->EndIndex)+1,(sep3G->EndIndex)+1);



  for (int gr=1;gr<=2;gr++)
  {
   TSeparator* sep1;
   TSeparator* sep2;
   TSeparator* sep3;
   T3dNumberMask* maskPtr;
   if (gr==1)
   {
	   sep1 = sep1F;
	   sep2 = sep2F;
	   sep3 = sep3F;
	   maskPtr = maskF;
   }
   else
   {
	   sep1 = sep1G;
	   sep2 = sep2G;
	   sep3 = sep3G;
	   maskPtr = maskG;
   }

   int N = sep1->EndIndex;
   int L = sep2->EndIndex;
   int M = sep3->EndIndex;

   T3dVector**** normals = new T3dVector***[N+1];
   T3dNumberMask& mask = *maskPtr;

   for (int i=0;i<N+1;i++)
   {
     normals[i] = new T3dVector**[L+1];
     for (int j=0;j<L+1;j++)
     {
       normals[i][j] = new T3dVector*[M+1];
       for (int k=0;k<M+1;k++)
       {
         T3dVector* normal = NULL;
		 mask[i][j][k]=TFictivePoint;

		 if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
         {
             mask[i][j][k]=TActualPoint;
         }
         else
         {
            if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==0) )
			{
				if ( (i!=1) && (i!=N-1) && (j!=1) && (j!=L-1) )
				{
				  mask[i][j][k]=TPreDefinedBorderPoint;
				}
				else
				{
				  mask[i][j][k]=TPreNormalBorderPoint;
				}
				normal = new T3dVector;
				normal->X = 0;
				normal->Y = 0;
				normal->Z = -1;
			}
            else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==M) )
            {
				if ( (i!=1) && (i!=N-1) && (j!=1) && (j!=L-1) )
				{
				  mask[i][j][k]=TPreDefinedBorderPoint;
				}
				else
				{
				  mask[i][j][k]=TPreNormalBorderPoint;
				}
				normal = new T3dVector;
				normal->X = 0;
				normal->Y = 0;
				normal->Z = 1;
			}

            else if ( (i==0) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
            {
				mask[i][j][k]=TPreNormalBorderPoint;
				normal = new T3dVector;
				normal->X = -1;
				normal->Y = 0;
				normal->Z = 0;
			}
            else if ( (i==N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreNormalBorderPoint;
				normal = new T3dVector;
				normal->X = 1;
				normal->Y = 0;
				normal->Z = 0;
			}

			else if ( (i!=0) && (i!=N) && (j==0) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreNormalBorderPoint;
				normal = new T3dVector;
				normal->X = 0;
				normal->Y = -1;
				normal->Z = 0;
			}
			else if ( (i!=0) && (i!=N) && (j==L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreNormalBorderPoint;
				normal = new T3dVector;
				normal->X = 0;
				normal->Y = 1;
				normal->Z = 0;
			}
         }
		 normals[i][j][k] = normal;
       }
     }
   }

   if (gr==1)
   {
     fGrid = new T3dNormalGrid(sep1,sep2,sep3,mask,normals);
   }
   else
   {
     fGridGeometry = new T3dNormalGrid(sep1,sep2,sep3,mask,normals);
   }

   delete &mask;

   for (int i=0;i<N+1;i++)
   {
     for (int j=0;j<L+1;j++)
	 {
       for (int k=0;k<M+1;k++)
	   {
			if ( (normals[i][j][k])!=NULL ) {delete normals[i][j][k];}
	   }
	   delete[] normals[i][j];
	 }
	 delete[] normals[i];
   }
   delete[] normals;// normals = NULL;
  } // by grids
}
//*************************************
T3dPoissonProblem::~T3dPoissonProblem()
{
   delete fGrid;
   printf("The problem has been completed successfully.\n");
}
//*************************************
double T3dPoissonProblem::SolutionFunction(double x, double y, double z)
{
	return sin(x)*sin(y)*sin(z);
}
//*************************************
double T3dPoissonProblem::SolutionFunctionDX(double x, double y, double z)
{
	return cos(x)*sin(y)*sin(z);
}
//*************************************
double T3dPoissonProblem::SolutionFunctionDY(double x, double y, double z)
{
	return sin(x)*cos(y)*sin(z);
}
//*************************************
double T3dPoissonProblem::SolutionFunctionDZ(double x, double y, double z)
{
	return sin(x)*sin(y)*cos(z);
}
//*************************************
double T3dPoissonProblem::RightPartFunction(double x, double y, double z)
{
   return 3*sin(x)*sin(y)*sin(z);
}
//*************************************
void T3dPoissonProblem::Solve()
{
   printf("The problem is starting ...\n");

   double h = fH;
   T3dNormalGrid& grid = *fGrid;

   int N = grid.GetSeparator1().EndIndex;
   int L = grid.GetSeparator2().EndIndex;
   int M = grid.GetSeparator3().EndIndex;

   const double* x = grid.GetSeparator1().Dimension;
   const double* y = grid.GetSeparator2().Dimension;
   const double* z = grid.GetSeparator3().Dimension;

   int NG = fGridGeometry->GetSeparator1().EndIndex;
   int LG = fGridGeometry->GetSeparator2().EndIndex;
   int MG = fGridGeometry->GetSeparator3().EndIndex;
   const double* xG = fGridGeometry->GetSeparator1().Dimension;
   const double* yG = fGridGeometry->GetSeparator2().Dimension;
   const double* zG = fGridGeometry->GetSeparator3().Dimension;

   const T3dNumberMask& mask = grid.Mask;

   const T3dVector* const* const* const* normals = grid.Normals;


   T3dNumberMask& spaceMask = mask.Copy();
   for (int i=0;i<N+1;i++)
   {
	   for (int j=0;j<L+1;j++)
	   {
		   for (int k=0;k<M+1;k++)
		   {
			   if (spaceMask[i][j][k]!=TFictivePoint)
			   {
					   spaceMask[i][j][k] = TActualPoint;
			   }
		   }
	   }
   }


   TMask** spaceMasks = new TMask*[1];
   spaceMasks[0] = &spaceMask;

   TRnRelease3dSpaceModel& spaceModel = *(new TRnRelease3dSpaceModel(1,spaceMasks));

   int spaceModelName = 555;
   TRnSpace::AddModel(spaceModelName,&spaceModel,false);

   TRnRelease3dSpace& rightPart = *(new TRnRelease3dSpace(spaceModelName));

   TRnRelease3dSpace& solution = *(new TRnRelease3dSpace(spaceModelName));

   TRnRelease3dSpace& approximation = *(new TRnRelease3dSpace(spaceModelName));

   TRnRelease3dSpace& error = *(new TRnRelease3dSpace(spaceModelName));

   for (int i=0;i<N+1;i++)
   {
	   for (int j=0;j<L+1;j++)
	   {
		   for (int k=0;k<M+1;k++)
		   {
			   error[i][j][k] = 0;
			   solution[i][j][k] = SolutionFunction(x[i],y[j],z[k]);
			   if (mask[i][j][k] == TFictivePoint)
			   {
			       rightPart[i][j][k] = 0;
			       approximation[i][j][k] = 0;
				   solution[i][j][k] = 0;
			   }
			   else if (mask[i][j][k] == TActualPoint)
			   {
			     rightPart[i][j][k] = RightPartFunction(x[i],y[j],z[k]);
			     approximation[i][j][k] = 0;
			   }
			   else if (mask[i][j][k]==TPreDefinedBorderPoint)
			   {
				   if ( i==0 ) {rightPart[i][j][k] = SolutionFunction(xG[i],y[j],z[k]);}
				   else if ( i==N ) {rightPart[i][j][k] = SolutionFunction(xG[i-1],y[j],z[k]);}
				   else if ( j==0 ) {rightPart[i][j][k] = SolutionFunction(x[i],yG[j],z[k]);}
				   else if ( j==L ) {rightPart[i][j][k] = SolutionFunction(x[i],yG[j-1],z[k]);}
				   else if ( k==0 ) {rightPart[i][j][k] = SolutionFunction(x[i],y[j],zG[k]);}
				   else if ( k==M ) {rightPart[i][j][k] = SolutionFunction(x[i],y[j],zG[k-1]);}
			       approximation[i][j][k] = 0;
			   }
			   else if (mask[i][j][k]==TPreNormalBorderPoint)
			   {
				   double nx = normals[i][j][k]->X;
				   double ny = normals[i][j][k]->Y;
				   double nz = normals[i][j][k]->Z;
				   if ( i==0 ) {rightPart[i][j][k] = nx*SolutionFunctionDX(xG[i],y[j],z[k]);}
				   else if ( i==N ) {rightPart[i][j][k] = nx*SolutionFunctionDX(xG[i-1],y[j],z[k]);}
				   else if ( j==0 ) {rightPart[i][j][k] = ny*SolutionFunctionDY(x[i],yG[j],z[k]);}
				   else if ( j==L ) {rightPart[i][j][k] = ny*SolutionFunctionDY(x[i],yG[j-1],z[k]);}
				   else if ( k==0 ) {rightPart[i][j][k] = nz*SolutionFunctionDZ(x[i],y[j],zG[k]);}
				   else if ( k==M ) {rightPart[i][j][k] = nz*SolutionFunctionDZ(x[i],y[j],zG[k-1]);}
			       approximation[i][j][k] = 0;
			   }
			   else
			   {
				   throw "Unknown point.";
			   }
		   }
	   }
   }

   TNodesGrid** grids = new TNodesGrid*[1];
   grids[0] = &grid;

   T3dLaplaceOperator& A = *(new T3dLaplaceOperator(spaceModelName,1,grids));
   delete[] grids;// grids = NULL;

   error|= ((TRnLinearOperator&)A)(solution) - rightPart;
   printf("\n");
   printf("The approximation norma is ");
   printf("%.20f",error.GetAbsMaxElement());
   printf("\n");


   //Pack operator (Begin)
   TRnRelease3dSpace& zeroHelp = *(new TRnRelease3dSpace(spaceModelName));
   TRnRelease3dSpace& rightPartHelp = *(new TRnRelease3dSpace(spaceModelName));

   zeroHelp.Assign(approximation);
   zeroHelp|=0;
   rightPartHelp|= ((TRnLinearOperator&)A)(zeroHelp) - rightPart;
   rightPartHelp|=(-1)*rightPartHelp;

   delete &zeroHelp;


   int RN = approximation.EndIndex;
   double* helpVector = new double[RN];
   for (int rn=0;rn<RN;rn++) {helpVector[rn]=0;}

   T1dNumberMask& spaceMaskRn = *(new T1dNumberMask(RN));
   for (int rn=0;rn<RN;rn++) {spaceMaskRn[rn] = TActualPoint;}
   TMask** spaceMasksRn = new TMask*[1];
   spaceMasksRn[0] = &spaceMaskRn;
   TRnRelease1dSpaceModel& spaceModelRn = *(new TRnRelease1dSpaceModel(1,spaceMasksRn));
   int spaceModelNameRn = 777;
   TRnSpace::AddModel(spaceModelNameRn,&spaceModelRn,false);


   TRnRelease1dSpace& rightPartRn = *(new TRnRelease1dSpace(spaceModelNameRn));
   TRnRelease1dSpace& approximationRn = *(new TRnRelease1dSpace(spaceModelNameRn));

   for (int rn=1;rn<RN+1;rn++)
   {
	   rightPartRn[rn-1] = ((TRnSpace&)rightPartHelp)[rn];
	   helpVector[rn-1] = rightPartRn[rn-1];

	   approximationRn[rn-1] = ((TRnSpace&)approximation)[rn];
   }


   delete &rightPartHelp;

   printf("Preparing pack operator: Starting...\n");
   bool flagTransform = true;
   T1dPackOperator* ARnPtr;
   if (flagTransform == false)
   {
	   double* nullVector = NULL;
       ARnPtr = new T1dPackOperator(spaceModelNameRn,A,&nullVector);
   }
   else
   {
	   ARnPtr = new T1dPackOperator(spaceModelNameRn,A,&helpVector);
	   for (int rn=0;rn<RN;rn++)
       {
         rightPartRn[rn] = helpVector[rn];
       }
   }

   delete[] helpVector;// helpVector = NULL;

   T1dPackOperator& ARn = *ARnPtr;
   printf("Preparing pack operator: OK.\n");


   TSLAE& slae = *(new TSLAE(ARn,approximationRn,rightPartRn));
   //slae.SolveByMinimalResiduals();
   slae.SolveByConjugateGradients();
   approximationRn|=slae.Result;
   delete &ARn;
   delete &slae;

   delete &spaceMaskRn;
   delete[] spaceMasksRn;//  spaceMaskRn = NULL;
   delete &spaceModelRn;


   for (int rn=1;rn<RN+1;rn++)
   {

	   ((TRnSpace&)approximation)[rn] = approximationRn[rn-1];
   }

   delete &rightPartRn;
   delete &approximationRn;
   //Pack operator (End)


   /*
   //Ordinary operator (Begin)
   TSLAE& slae = *(new TSLAE(A,approximation,rightPart));
   //A.ConstructMatrix();
   //slae.SolveByGauss();
   //slae.SolveByConjugateGradients();
   slae.SolveByMinimalResiduals();
   approximation|=slae.Result;
   delete &A;
   delete &slae;
   //Ordinary operator (Begin)
   */

   error|=solution-approximation;
   printf("\n");
   printf("The error norma is ");
   printf("%.20f",error.GetAbsMaxElement());
   printf("\n");
   printf("\n");


   delete &spaceMask;
   delete[] spaceMasks;
   delete &spaceModel;

   delete &rightPart;
   delete &solution;
   delete &approximation;
   delete &error;

   TRnSpace::RemoveAllModels();
}
//*************************************



//#################################################################


void THydrodynamicsProblem::SolveTransferEquation(TRnRelease3dSpace& lVecCur, TRnRelease3dSpace& lVecPrev,
		                               TRnRelease3dSpace&  lU, TRnRelease3dSpace& lV, TRnRelease3dSpace& lW, TRnRelease3dSpace& lRP,
								       const TNodesGrid& lGrid)
{
    SolveTransferEquation(lVecCur,lVecPrev,lU,lV,lW,lRP,lGrid,1/fRe);
}
//*************************************
void THydrodynamicsProblem::SolveTransferEquation(TRnRelease3dSpace& lVecCur, TRnRelease3dSpace& lVecPrev,
		                               TRnRelease3dSpace&  lU, TRnRelease3dSpace& lV, TRnRelease3dSpace& lW, TRnRelease3dSpace& lRP,
								       const TNodesGrid& lGrid, double lDiffusionCoefficient)
{	
    //Preparing variables (Begin)
	TRnRelease3dSpace& f = lVecCur;
    TRnRelease3dSpace& f_1 = lVecPrev;
    TRnRelease3dSpace& u = lU;
	TRnRelease3dSpace& v = lV;
	TRnRelease3dSpace& w = lW;
	TRnRelease3dSpace& rp = lRP;

	TRnRelease3dSpace& f_1_3 = dynamic_cast<TRnRelease3dSpace&>(f_1.Copy());
	TRnRelease3dSpace& f_2_3 = dynamic_cast<TRnRelease3dSpace&>(f_1.Copy());

	const T3dNormalGrid& grid = dynamic_cast<const T3dNormalGrid&>(lGrid);

	const T3dNumberMask& mask = grid.Mask;

	const TSeparator& sepX = grid.Separator1;
	const TSeparator& sepY = grid.Separator2;
	const TSeparator& sepZ = grid.Separator3;

	int N = sepX.EndIndex;
    int L = sepY.EndIndex;
    int M = sepZ.EndIndex;

	const double* hx = sepX.Separation;
	const double* hy = sepY.Separation;
	const double* hz = sepZ.Separation;

	double tau = fTimeSeparator->SeparationValue;

	double diffCoef = lDiffusionCoefficient;

	//Preparing variables (End)

	//Direction order (Begin)
	int directionOrder = 1;
	int directionX;
	int directionY;
	int directionZ;
	if (directionOrder == 1) {directionX=1;directionY=2;directionZ=3;}
	else if (directionOrder == 2) {directionX=1;directionY=3;directionZ=2;}
	else if (directionOrder == 3) {directionX=2;directionY=1;directionZ=3;}
	else if (directionOrder == 4) {directionX=2;directionY=3;directionZ=1;}
	else if (directionOrder == 5) {directionX=3;directionY=1;directionZ=2;}
	else if (directionOrder == 6) {directionX=3;directionY=2;directionZ=1;}
	else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction order.";}
	//Direction order (End)

	//Terms using (Begin)
	int flagConvectiveX = 1; int flagDiffusionX = 1;
    int flagConvectiveY = 1; int flagDiffusionY = 1;
    int flagConvectiveZ = 1; int flagDiffusionZ = 1;
	//Terms using (Begin)


	for (int directionGlobal=1;directionGlobal<4;directionGlobal++) //splitting scheame steps
	{
	for (int directionFake=1;directionFake<4;directionFake++) //for border condtions recalculating
    {
	   int direction;
	   if (directionFake==1) {direction = directionGlobal;}
	   else
	   {
		   if (directionGlobal==1) {direction=directionFake;}
		   else if (directionGlobal==2)
		   {
			   if (directionFake==2) {direction=1;}
			   else if (directionFake==3) {direction=3;}
			   else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown fake direction.";}
		   }
		   else if (directionGlobal==3) {direction=directionFake-1;}
		   else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown global direction.";}
	   }


	   int dir1EndIndex;
	   int dir2EndIndex;
	   int dirActEndIndex;

	   if (direction==directionX)
	   {
		 dir1EndIndex = M;
	     dir2EndIndex = L;
	     dirActEndIndex = N;
	   }
	   else if (direction==directionY)
	   {
         dir1EndIndex = M;
	     dir2EndIndex = N;
	     dirActEndIndex = L;
	   }
	   else if (direction==directionZ)
	   {
	     dir1EndIndex = L;
	     dir2EndIndex = N;
	     dirActEndIndex = M;
	   }
	   else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}

	   int flagConvective = 0;
	   if (
		    ( (direction==directionX) && (flagConvectiveX==1) ) ||
		    ( (direction==directionY) && (flagConvectiveY==1) ) ||
			( (direction==directionZ) && (flagConvectiveZ==1) )
		  ) {flagConvective=1;}

	   int flagDiffusion = 0;
	   if (
		    ( (direction==directionX) && (flagDiffusionX==1) ) ||
		    ( (direction==directionY) && (flagDiffusionY==1) ) ||
			( (direction==directionZ) && (flagDiffusionZ==1) )
		  ) {flagDiffusion=1;}


	   for (int dir1=0;dir1<dir1EndIndex+1;dir1++)
	   {
         for (int dir2=0;dir2<dir2EndIndex+1;dir2++)
         {
           int dirActSecBeg = -1;
		   int maskActSecBeg;

           int dirActSecEnd;
		   int maskActSecEnd;

           for (int dirActSec=0; dirActSec<dirActEndIndex+1; dirActSec++)
           {
			     int maskActSec;
                 if (direction==directionX) {maskActSec = mask[dirActSec][dir2][dir1];}
	             else if (direction==directionY) {maskActSec = mask[dir2][dirActSec][dir1];}
	             else if (direction==directionZ) {maskActSec = mask[dir2][dir1][dirActSec];}
				 else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}

                 if (dirActSecBeg<0)
                 {
                   if (maskActSec==TActualPoint)
                   {
                      if (dirActSec==0)
                      {
                        throw "Error in THydrodynamicsProblem::SolveTransferEquation: Initial border not found.";
                      }
                      else
                      {
                        dirActSecBeg = dirActSec-1;
						if (direction==directionX) {maskActSecBeg = mask[dirActSecBeg][dir2][dir1];}
	                    else if (direction==directionY) {maskActSecBeg = mask[dir2][dirActSecBeg][dir1];}
	                    else if (direction==directionZ) {maskActSecBeg = mask[dir2][dir1][dirActSecBeg];}
				        else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
                      }
                   }
                 }
                 else if (maskActSec!=TActualPoint)
                 {
                    dirActSecEnd=dirActSec;
                    maskActSecEnd = maskActSec;

					//Сhecking sector (Begin)
					if ( (maskActSecBeg == TActualPoint) || (maskActSecBeg == TFictivePoint) || (maskActSecEnd == TActualPoint) || (maskActSecEnd == TFictivePoint) )
					{
                        throw "Error in THydrodynamicsProblem::SolveTransferEquation: Incorrect sector border.";
					}

					if ( (dirActSecBeg<0) || (dirActSecEnd-dirActSecBeg)<2 )
					{
                        throw "Error in THydrodynamicsProblem::SolveTransferEquation: Incorrect sector size.";
					}
					//Сhecking sector (End)



					//Only borders recalculation and continue (Begin)
					if (direction!=directionGlobal)
					{
                      TRnRelease3dSpace* f_AssignPtr;
					  if (directionGlobal==1) {f_AssignPtr = &f_1_3;}
				      else if (directionGlobal==2) {f_AssignPtr = &f_2_3;}
				      else if (directionGlobal==3) {f_AssignPtr = &f;}
					  else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown global direction.";}
					  TRnRelease3dSpace& f_Assign = *f_AssignPtr;

					  //Sector begin (Begin)
					  if  ( (maskActSecBeg == TDefinedBorderPoint) || (maskActSecBeg == TEquationBorderPoint) )
				   	  {
					        if (direction==directionX) {f_Assign[dirActSecBeg][dir2][dir1] = rp[dirActSecBeg][dir2][dir1];}
	                        else if (direction==directionY) {f_Assign[dir2][dirActSecBeg][dir1] = rp[dir2][dirActSecBeg][dir1];}
	                        else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecBeg] = rp[dir2][dir1][dirActSecBeg];}
				            else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
					  else if  (maskActSecBeg == TPreDefinedBorderPoint)
				   	  {
					        if (direction==directionX) {f_Assign[dirActSecBeg][dir2][dir1] = 2*rp[dirActSecBeg][dir2][dir1] - f_Assign[dirActSecBeg+1][dir2][dir1];}
	                        else if (direction==directionY) {f_Assign[dir2][dirActSecBeg][dir1] = 2*rp[dir2][dirActSecBeg][dir1] - f_Assign[dir2][dirActSecBeg+1][dir1];}
	                        else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecBeg] = 2*rp[dir2][dir1][dirActSecBeg] - f_Assign[dir2][dir1][dirActSecBeg+1];}
				            else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
					  else if ( (maskActSecBeg == TNormalBorderPoint) || (maskActSecBeg == TPreNormalBorderPoint) )
					  {
                       //! need to check acceptable normal direction

					   if (direction==directionX) {f_Assign[dirActSecBeg][dir2][dir1] = -(hx[dirActSecBeg])*rp[dirActSecBeg][dir2][dir1] + f_Assign[dirActSecBeg+1][dir2][dir1];}
	                   else if (direction==directionY) {f_Assign[dir2][dirActSecBeg][dir1] = -(hy[dirActSecBeg])*rp[dir2][dirActSecBeg][dir1] + f_Assign[dir2][dirActSecBeg+1][dir1];}
	                   else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecBeg] = -(hz[dirActSecBeg])*rp[dir2][dir1][dirActSecBeg] + f_Assign[dir2][dir1][dirActSecBeg+1];}
				       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
					  else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of begin sector.";}
                      //Sector begin (End)


					  //Sector end (Begin)
                      if ( (maskActSecEnd == TDefinedBorderPoint) || (maskActSecEnd == TEquationBorderPoint) )
					  {
					        if (direction==directionX) {f_Assign[dirActSecEnd][dir2][dir1] = rp[dirActSecEnd][dir2][dir1];}
					        else if (direction==directionY) {f_Assign[dir2][dirActSecEnd][dir1] = rp[dir2][dirActSecEnd][dir1];}
					        else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecEnd] = rp[dir2][dir1][dirActSecEnd];}
				            else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
                      else if (maskActSecEnd == TPreDefinedBorderPoint)
					  {
					        if (direction==directionX) {f_Assign[dirActSecEnd][dir2][dir1] = 2*rp[dirActSecEnd][dir2][dir1] - f_Assign[dirActSecEnd-1][dir2][dir1];}
					        else if (direction==directionY) {f_Assign[dir2][dirActSecEnd][dir1] = 2*rp[dir2][dirActSecEnd][dir1] - f_Assign[dir2][dirActSecEnd-1][dir1];}
					        else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecEnd] = 2*rp[dir2][dir1][dirActSecEnd] - f_Assign[dir2][dir1][dirActSecEnd-1];}
				            else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
					  else if ( (maskActSecEnd == TNormalBorderPoint) || (maskActSecEnd == TPreNormalBorderPoint) )
					  {
						    //! need to check acceptable normal direction
					        if (direction==directionX) {f_Assign[dirActSecEnd][dir2][dir1] = (hx[dirActSecEnd-1])*rp[dirActSecEnd][dir2][dir1] + f_Assign[dirActSecEnd-1][dir2][dir1];}
					        else if (direction==directionY) {f_Assign[dir2][dirActSecEnd][dir1] = (hy[dirActSecEnd-1])*rp[dir2][dirActSecEnd][dir1] + f_Assign[dir2][dirActSecEnd-1][dir1];}
					        else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecEnd] = (hz[dirActSecEnd-1])*rp[dir2][dir1][dirActSecEnd] + f_Assign[dir2][dir1][dirActSecEnd-1];}
				            else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
					  else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of end sector.";}
					  //Sector end (End)

					  dirActSecBeg=-1; //sector has been calculated
					  continue;
					} //if it is fake direction
					//Only borders recalculation and continue (End)





					//Need for fictive points (Begin)
					int fictiveLeft = 0;
					if ( (maskActSecBeg!=TDefinedBorderPoint) && (maskActSecBeg!=TEquationBorderPoint) )  {fictiveLeft=1;}

					int fictiveRight = 0;
					if ( (maskActSecEnd!=TDefinedBorderPoint) && (maskActSecEnd!=TEquationBorderPoint) )  {fictiveRight=1;}

					int dimActSec = dirActSecEnd-dirActSecBeg+1+fictiveLeft+fictiveRight;
					//Need for fictive points (End)


				    //Preparing scalar objects (Begin)
					int scalDim = dimActSec-1;
					double* scalVar = new double[scalDim+1];
					double* scalA = new double[scalDim+1];
					double* scalB = new double[scalDim+1];
					double* scalC = new double[scalDim+1];
					double* scalRP = new double[scalDim+1];

					scalA[0]  = 0;  scalA[scalDim]  = 0;
					scalB[0]  = 0;  scalB[scalDim]  = 0;
					scalC[0]  = 0;  scalC[scalDim]  = 0;
					scalRP[0] = 0;  scalRP[scalDim] = 0;


					//Borders points (Begin)
					if ( (maskActSecBeg == TDefinedBorderPoint) || (maskActSecBeg == TEquationBorderPoint) )
					{
						if (direction==directionX) {scalVar[0] = rp[dirActSecBeg][dir2][dir1];}
	                    else if (direction==directionY) {scalVar[0] = rp[dir2][dirActSecBeg][dir1];}
	                    else if (direction==directionZ) {scalVar[0] = rp[dir2][dir1][dirActSecBeg];}
				        else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else if  (maskActSecBeg == TPreDefinedBorderPoint)
					{
					   scalVar[0] = 0;
					   scalA[1] = 1;
					   scalB[1] = 0.5;
					   scalC[1] = 0.5;
					   if (direction==directionX) {scalRP[1] = rp[dirActSecBeg][dir2][dir1];}
	                   else if (direction==directionY) {scalRP[1] = rp[dir2][dirActSecBeg][dir1];}
	                   else if (direction==directionZ) {scalRP[1] = rp[dir2][dir1][dirActSecBeg];}
				       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else if  ( (maskActSecBeg == TNormalBorderPoint) ||  (maskActSecBeg == TPreNormalBorderPoint) )
					{
                       //! need to check acceptable normal direction
					   scalVar[0] = 0;
					   scalA[1] = 1;
					   if (direction==directionX)
					   {
						   scalB[1] = -1/hx[dirActSecBeg];
					       scalC[1] = 1/hx[dirActSecBeg];
						   scalRP[1] = rp[dirActSecBeg][dir2][dir1];
					   }
	                   else if (direction==directionY)
					   {
						   scalB[1] = -1/hy[dirActSecBeg];
					       scalC[1] = 1/hy[dirActSecBeg];
						   scalRP[1] = rp[dir2][dirActSecBeg][dir1];
					   }
	                   else if (direction==directionZ)
					   {
						   scalB[1] = -1/hz[dirActSecBeg];
					       scalC[1] = 1/hz[dirActSecBeg];
						   scalRP[1] = rp[dir2][dir1][dirActSecBeg];
					   }
				       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of begin sector.";}



					if ( (maskActSecEnd == TDefinedBorderPoint) || (maskActSecEnd == TEquationBorderPoint) )
					{
						if (direction==directionX) {scalVar[scalDim] = rp[dirActSecEnd][dir2][dir1];}
	                    else if (direction==directionY) {scalVar[scalDim] = rp[dir2][dirActSecEnd][dir1];}
	                    else if (direction==directionZ) {scalVar[scalDim] = rp[dir2][dir1][dirActSecEnd];}
				        else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else if (maskActSecEnd == TPreDefinedBorderPoint)
					{
					   scalVar[scalDim] = 0;
					   scalC[scalDim-1] = 1;
					   scalA[scalDim-1] = 0.5;
					   scalB[scalDim-1] = 0.5;
					   if (direction==directionX) {scalRP[scalDim-1] = rp[dirActSecEnd][dir2][dir1];}
					   else if (direction==directionY) {scalRP[scalDim-1] = rp[dir2][dirActSecEnd][dir1];}
					   else if (direction==directionZ) {scalRP[scalDim-1] = rp[dir2][dir1][dirActSecEnd];}
				       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else if  ( (maskActSecEnd == TNormalBorderPoint) || (maskActSecEnd == TPreNormalBorderPoint) )
					{
                       scalVar[scalDim] = 0;
					   scalC[scalDim-1] = 1;
					   if (direction==directionX)
					   {
						   scalA[scalDim-1] = -1/hx[dirActSecEnd-1];
					       scalB[scalDim-1] = 1/hx[dirActSecEnd-1];
						   scalRP[scalDim-1] = rp[dirActSecEnd][dir2][dir1];
					   }
	                   else if (direction==directionY)
					   {
						   scalA[scalDim-1] = -1/hy[dirActSecEnd-1];
					       scalB[scalDim-1] = 1/hy[dirActSecEnd-1];
						   scalRP[scalDim-1] = rp[dir2][dirActSecEnd][dir1];
					   }
	                   else if (direction==directionZ)
					   {
						   scalA[scalDim-1] = -1/hz[dirActSecEnd-1];
					       scalB[scalDim-1] = 1/hz[dirActSecEnd-1];
						   scalRP[scalDim-1] = rp[dir2][dir1][dirActSecEnd];
					   }
				       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of end sector.";}
					//Borders points (End)


					//Actual points (Begin)
					for (int scalLoop = 1+fictiveLeft;scalLoop<scalDim-fictiveRight;scalLoop++)
					{
						int dirAct = dirActSecBeg + scalLoop - fictiveLeft;

						int i; int j; int k;
						TRnRelease3dSpace* vecTrPtr;
						if (direction==directionX) {i = dirAct; j = dir2; k = dir1; vecTrPtr = &u;}
						else if (direction==directionY) {i = dir2; j = dirAct; k = dir1; vecTrPtr = &v;}
						else if (direction==directionZ) {i = dir2; j = dir1; k = dirAct; vecTrPtr = &w;}
						else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
						TRnRelease3dSpace& vecTr = *vecTrPtr;

						if (mask[i][j][k]!=TActualPoint) {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Incorrect sector innerr point.";}

						double h_1X; double hX; double hx2X; double hxhxhX;
						h_1X = hx[i-1];
					    hX = hx[i];
						hx2X = hx[i-1]+hx[i];
						hxhxhX = hx[i-1]*hx[i]*(hx[i-1]+hx[i])/2;

						double h_1Y; double hY; double hx2Y; double hxhxhY;
						h_1Y = hy[j-1];
					    hY = hy[j];
					    hx2Y = hy[j-1]+hy[j];
						hxhxhY = hy[j-1]*hy[j]*(hy[j-1]+hy[j])/2;

						double h_1Z; double hZ; double hx2Z; double hxhxhZ;
						h_1Z = hz[k-1];
						hZ = hz[k];
						hx2Z = hz[k-1]+hz[k];
						hxhxhZ = hz[k-1]*hz[k]*(hz[k-1]+hz[k])/2;

						double h_1; double h; double hx2; double hxhxh;
						if (direction==directionX) {h_1 = h_1X; h = hX; hx2 = hx2X; hxhxh = hxhxhX;}
						else if (direction==directionY) {h_1 = h_1Y; h = hY; hx2 = hx2Y; hxhxh = hxhxhY;}
						else if (direction==directionZ) {h_1 = h_1Z; h = hZ; hx2 = hx2Z; hxhxh = hxhxhZ;}
					    else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}


						scalA[scalLoop] = flagConvective*( 0.5*(-vecTr[i][j][k]/hx2) ) + flagDiffusion*( 0.5*(-(diffCoef)*h/hxhxh) );
						scalB[scalLoop] = 1/tau + flagDiffusion*( 0.5*(diffCoef)*hx2/hxhxh );
						scalC[scalLoop] = flagConvective*( 0.5*(vecTr[i][j][k]/hx2) ) + flagDiffusion*( 0.5*(-(diffCoef)*h_1/hxhxh) );

						double F = (1/tau)*f_1[i][j][k] + rp[i][j][k];

						//Help variables for simple writing F (Begin)
						double FC;
						double FD;
						TRnRelease3dSpace* fDir=NULL;
						//Help variables for simple writing F (End)

						if (direction==directionX)
					    {
						   F = F + flagConvectiveX*( 0.5*(-u[i][j][k])*(f_1[i+1][j][k]-f_1[i-1][j][k])/hx2X );
						   F = F + flagDiffusionX*(  0.5*(diffCoef)*( (hX*f_1[i-1][j][k] - hx2X*f_1[i][j][k] + h_1X*f_1[i+1][j][k])/hxhxhX )  );

						   //Y (Begin)
						   if (directionY>direction) {fDir = &f_1;}
						   else if (directionY==1) {fDir = &f_1_3;}
						   else if (directionY==2) {fDir = &f_2_3;}

						   FC = 0; FD = 0;
						   for (int item=1;item<=2;item++)
						   {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   FC = FC + 0.5*(-v[i][j][k])*(f_[i][j+1][k]-f_[i][j-1][k])/hx2Y;
						       FD = FD + 0.5*(diffCoef)*( (hY*f_[i][j-1][k] - hx2Y*f_[i][j][k] + h_1Y*f_[i][j+1][k])/hxhxhY );
						   }
						   F = F + flagConvectiveY*FC + flagDiffusionY*FD;
						   //Y (End)

						   //Z (Begin)
						   if (directionZ>direction) {fDir = &f_1;}
						   else if (directionZ==1) {fDir = &f_1_3;}
						   else if (directionZ==2) {fDir = &f_2_3;}

						   FC = 0; FD = 0;
						   for (int item=1;item<=2;item++)
						   {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k+1]-f_[i][j][k-1])/hx2Z;
						       FD = FD + 0.5*(diffCoef)*( (hZ*f_[i][j][k-1] - hx2Z*f_[i][j][k] + h_1Z*f_[i][j][k+1])/hxhxhZ );
						   }
						   F = F + flagConvectiveZ*FC + flagDiffusionZ*FD;
						   //Z (End)
						}
						else if (direction==directionY)
					    {
						    F = F + flagConvectiveY*( 0.5*(-v[i][j][k])*(f_1[i][j+1][k]-f_1[i][j-1][k])/hx2Y );
						    F = F + flagDiffusionY*(  0.5*(diffCoef)*( (hY*f_1[i][j-1][k] - hx2Y*f_1[i][j][k] + h_1Y*f_1[i][j+1][k])/hxhxhY )  );

                            //X (Begin)
						    if (directionX>direction) {fDir = &f_1;}
						    else if (directionX==1) {fDir = &f_1_3;}
						    else if (directionX==2) {fDir = &f_2_3;}

						    FC = 0; FD = 0;
						    for (int item=1;item<=2;item++)
						    {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   FC = FC + 0.5*(-u[i][j][k])*(f_[i+1][j][k]-f_[i-1][j][k])/hx2X;
						       FD = FD + 0.5*(diffCoef)*( (hX*f_[i-1][j][k] - hx2X*f_[i][j][k] + h_1X*f_[i+1][j][k])/hxhxhX );
						    }
						    F = F + flagConvectiveX*FC + flagDiffusionX*FD;
						    //X (End)


							//Z (Begin)
						    if (directionZ>direction) {fDir = &f_1;}
						    else if (directionZ==1) {fDir = &f_1_3;}
						    else if (directionZ==2) {fDir = &f_2_3;}

						    FC = 0; FD = 0;
						    for (int item=1;item<=2;item++)
						    {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k+1]-f_[i][j][k-1])/hx2Z;
						       FD = FD + 0.5*(diffCoef)*( (hZ*f_[i][j][k-1] - hx2Z*f_[i][j][k] + h_1Z*f_[i][j][k+1])/hxhxhZ );
						    }
						    F = F + flagConvectiveZ*FC + flagDiffusionZ*FD;
						    //Z (End)
					    }
	                    else if (direction==directionZ)
					    {
						    F = F + flagConvectiveZ*( 0.5*(-w[i][j][k])*(f_1[i][j][k+1]-f_1[i][j][k-1])/hx2Z );
						    F = F + flagDiffusionZ*(  0.5*(diffCoef)*( (hZ*f_1[i][j][k-1] - hx2Z*f_1[i][j][k] + h_1Z*f_1[i][j][k+1])/hxhxhZ )  );

						    //X (Begin)
						    if (directionX>direction) {fDir = &f_1;}
						    else if (directionX==1) {fDir = &f_1_3;}
						    else if (directionX==2) {fDir = &f_2_3;}

						    FC = 0; FD = 0;
						    for (int item=1;item<=2;item++)
						    {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   FC = FC + 0.5*(-u[i][j][k])*(f_[i+1][j][k]-f_[i-1][j][k])/hx2X;
						       FD = FD + 0.5*(diffCoef)*( (hX*f_[i-1][j][k] - hx2X*f_[i][j][k] + h_1X*f_[i+1][j][k])/hxhxhX );
						    }
						    F = F + flagConvectiveX*FC + flagDiffusionX*FD;
						    //X (End)

							//Y (Begin)
						    if (directionY>direction) {fDir = &f_1;}
						    else if (directionY==1) {fDir = &f_1_3;}
						    else if (directionY==2) {fDir = &f_2_3;}

						    FC = 0; FD = 0;
						    for (int item=1;item<=2;item++)
						    {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   FC = FC + 0.5*(-v[i][j][k])*(f_[i][j+1][k]-f_[i][j-1][k])/hx2Y;
						       FD = FD + 0.5*(diffCoef)*( (hY*f_[i][j-1][k] - hx2Y*f_[i][j][k] + h_1Y*f_[i][j+1][k])/hxhxhY );
						    }
						    F = F + flagConvectiveY*FC + flagDiffusionY*FD;
						    //Y (End)
					    }
				        else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
						scalRP[scalLoop] = F;
					} // by act direction
					//Actual points (End)

					//Preparing scalar objects (End)


					//Scalar running (Begin)
					TSLAE::SolveByScalarRunning(scalDim,scalVar,scalA,scalB,scalC,scalRP);
					//Scalar running (End)


					//Assigninig (Begin)
					TRnRelease3dSpace* f_AssignPtr;
					if (direction==1) {f_AssignPtr = &f_1_3;}
				    else if (direction==2) {f_AssignPtr = &f_2_3;}
				    else if (direction==3) {f_AssignPtr = &f;}
					else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					TRnRelease3dSpace& f_Assign = *f_AssignPtr;

					for (int scalLoop = 0+fictiveLeft;scalLoop<scalDim+1-fictiveRight;scalLoop++)
					{
						int dirAct = dirActSecBeg + scalLoop - fictiveLeft;
						int i; int j; int k;
						if (direction==directionX) {i = dirAct; j = dir2; k = dir1;}
						else if (direction==directionY) {i = dir2; j = dirAct; k = dir1;}
						else if (direction==directionZ) {i = dir2; j = dir1; k = dirAct;}
						else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}

						f_Assign[i][j][k] = scalVar[scalLoop];
					}
					//Assigninig (End)

					//Deleting scalar objects (Begin)
					delete[] scalVar;// scalVar = NULL;
					delete[] scalA;// scalA = NULL;
					delete[] scalB;// scalB = NULL;
					delete[] scalC;// scalC = NULL;
					delete[] scalRP;// scalRP = NULL;
					//Deleting scalar objects (End)

                    dirActSecBeg=-1; //sector has been calculated
                 } // sector calculating
           } // by sectors of act direction
         } // by direction 2
       } // by direction 1
    } // by fake direction (x,y,z) - for border conditions recalculating
	} // by global direction (x,y,z)
	delete &f_1_3;
	delete &f_2_3;
}
//*************************************
void THydrodynamicsProblem::SolveTransferEquation_(TRnRelease3dSpace& lVecCur, TRnRelease3dSpace& lVecPrev,
		                               TRnRelease3dSpace&  lU, TRnRelease3dSpace& lV, TRnRelease3dSpace& lW, TRnRelease3dSpace& lRP,
								       const TNodesGrid& lGrid, double lDiffusionCoefficient)
{	
    //Preparing variables (Begin)
	TRnRelease3dSpace& f = lVecCur;
    TRnRelease3dSpace& f_1 = lVecPrev;
    TRnRelease3dSpace& u = lU;
	TRnRelease3dSpace& v = lV;
	TRnRelease3dSpace& w = lW;
	TRnRelease3dSpace& rp = lRP;

	TRnRelease3dSpace& f_1_3 = dynamic_cast<TRnRelease3dSpace&>(f_1.Copy());
	TRnRelease3dSpace& f_2_3 = dynamic_cast<TRnRelease3dSpace&>(f_1.Copy());

	const T3dNormalGrid& grid = dynamic_cast<const T3dNormalGrid&>(lGrid);

	const T3dNumberMask& mask = grid.Mask;

	const TSeparator& sepX = grid.Separator1;
	const TSeparator& sepY = grid.Separator2;
	const TSeparator& sepZ = grid.Separator3;

	int N = sepX.EndIndex;
    int L = sepY.EndIndex;
    int M = sepZ.EndIndex;

	const double* hx = sepX.Separation;
	const double* hy = sepY.Separation;
	const double* hz = sepZ.Separation;

	double tau = fTimeSeparator->SeparationValue;

	double diffCoef = lDiffusionCoefficient;

	//Preparing variables (End)

	//Direction order (Begin)
	int directionOrder = 1;
	int directionX;
	int directionY;
	int directionZ;
	if (directionOrder == 1) {directionX=1;directionY=2;directionZ=3;}
	else if (directionOrder == 2) {directionX=1;directionY=3;directionZ=2;}
	else if (directionOrder == 3) {directionX=2;directionY=1;directionZ=3;}
	else if (directionOrder == 4) {directionX=2;directionY=3;directionZ=1;}
	else if (directionOrder == 5) {directionX=3;directionY=1;directionZ=2;}
	else if (directionOrder == 6) {directionX=3;directionY=2;directionZ=1;}
	else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction order.";}
	//Direction order (End)

	//Terms using (Begin)
	int flagConvectiveX = 1; int flagDiffusionX = 1;
    int flagConvectiveY = 1; int flagDiffusionY = 1;
    int flagConvectiveZ = 1; int flagDiffusionZ = 1;
	//Terms using (Begin)

	//Scheme type (Begin)
	 //0 - central;
	 //1 - upstream;
	 int schemeType = 1;
	//Scheme type (End)

	for (int directionGlobal=1;directionGlobal<4;directionGlobal++) //splitting scheame steps
	{
	for (int directionFake=1;directionFake<4;directionFake++) //for border condtions recalculating
    {
	   int direction;
	   if (directionFake==1) {direction = directionGlobal;}
	   else
	   {
		   if (directionGlobal==1) {direction=directionFake;}
		   else if (directionGlobal==2)
		   {
			   if (directionFake==2) {direction=1;}
			   else if (directionFake==3) {direction=3;}
			   else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown fake direction.";}
		   }
		   else if (directionGlobal==3) {direction=directionFake-1;}
		   else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown global direction.";}
	   }


	   int dir1EndIndex;
	   int dir2EndIndex;
	   int dirActEndIndex;

	   if (direction==directionX)
	   {
		 dir1EndIndex = M;
	     dir2EndIndex = L;
	     dirActEndIndex = N;
	   }
	   else if (direction==directionY)
	   {
         dir1EndIndex = M;
	     dir2EndIndex = N;
	     dirActEndIndex = L;
	   }
	   else if (direction==directionZ)
	   {
	     dir1EndIndex = L;
	     dir2EndIndex = N;
	     dirActEndIndex = M;
	   }
	   else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}

	   int flagConvective = 0;
	   if (
		    ( (direction==directionX) && (flagConvectiveX==1) ) ||
		    ( (direction==directionY) && (flagConvectiveY==1) ) ||
			( (direction==directionZ) && (flagConvectiveZ==1) )
		  ) {flagConvective=1;}

	   int flagDiffusion = 0;
	   if (
		    ( (direction==directionX) && (flagDiffusionX==1) ) ||
		    ( (direction==directionY) && (flagDiffusionY==1) ) ||
			( (direction==directionZ) && (flagDiffusionZ==1) )
		  ) {flagDiffusion=1;}


	   for (int dir1=0;dir1<dir1EndIndex+1;dir1++)
	   {
         for (int dir2=0;dir2<dir2EndIndex+1;dir2++)
         {
           int dirActSecBeg = -1;
		   int maskActSecBeg;

           int dirActSecEnd;
		   int maskActSecEnd;

           for (int dirActSec=0; dirActSec<dirActEndIndex+1; dirActSec++)
           {
			     int maskActSec;
                 if (direction==directionX) {maskActSec = mask[dirActSec][dir2][dir1];}
	             else if (direction==directionY) {maskActSec = mask[dir2][dirActSec][dir1];}
	             else if (direction==directionZ) {maskActSec = mask[dir2][dir1][dirActSec];}
				 else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}

                 if (dirActSecBeg<0)
                 {
                   if (maskActSec==TActualPoint)
                   {
                      if (dirActSec==0)
                      {
                        throw "Error in THydrodynamicsProblem::SolveTransferEquation: Initial border not found.";
                      }
                      else
                      {
                        dirActSecBeg = dirActSec-1;
						if (direction==directionX) {maskActSecBeg = mask[dirActSecBeg][dir2][dir1];}
	                    else if (direction==directionY) {maskActSecBeg = mask[dir2][dirActSecBeg][dir1];}
	                    else if (direction==directionZ) {maskActSecBeg = mask[dir2][dir1][dirActSecBeg];}
				        else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
                      }
                   }
                 }
                 else if (maskActSec!=TActualPoint)
                 {
                    dirActSecEnd=dirActSec;
                    maskActSecEnd = maskActSec;

					//Сhecking sector (Begin)
					if ( (maskActSecBeg == TActualPoint) || (maskActSecBeg == TFictivePoint) || (maskActSecEnd == TActualPoint) || (maskActSecEnd == TFictivePoint) )
					{
                        throw "Error in THydrodynamicsProblem::SolveTransferEquation: Incorrect sector border.";
					}

					if ( (dirActSecBeg<0) || (dirActSecEnd-dirActSecBeg)<2 )
					{
                        throw "Error in THydrodynamicsProblem::SolveTransferEquation: Incorrect sector size.";
					}
					//Сhecking sector (End)



					//Only borders recalculation and continue (Begin)
					if (direction!=directionGlobal)
					{
                      TRnRelease3dSpace* f_AssignPtr;
					  if (directionGlobal==1) {f_AssignPtr = &f_1_3;}
				      else if (directionGlobal==2) {f_AssignPtr = &f_2_3;}
				      else if (directionGlobal==3) {f_AssignPtr = &f;}
					  else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown global direction.";}
					  TRnRelease3dSpace& f_Assign = *f_AssignPtr;

					  //Sector begin (Begin)
					  if  ( (maskActSecBeg == TDefinedBorderPoint) || (maskActSecBeg == TEquationBorderPoint) )
				   	  {
					        if (direction==directionX) {f_Assign[dirActSecBeg][dir2][dir1] = rp[dirActSecBeg][dir2][dir1];}
	                        else if (direction==directionY) {f_Assign[dir2][dirActSecBeg][dir1] = rp[dir2][dirActSecBeg][dir1];}
	                        else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecBeg] = rp[dir2][dir1][dirActSecBeg];}
				            else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
					  else if  (maskActSecBeg == TPreDefinedBorderPoint)
				   	  {
					        if (direction==directionX) {f_Assign[dirActSecBeg][dir2][dir1] = 2*rp[dirActSecBeg][dir2][dir1] - f_Assign[dirActSecBeg+1][dir2][dir1];}
	                        else if (direction==directionY) {f_Assign[dir2][dirActSecBeg][dir1] = 2*rp[dir2][dirActSecBeg][dir1] - f_Assign[dir2][dirActSecBeg+1][dir1];}
	                        else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecBeg] = 2*rp[dir2][dir1][dirActSecBeg] - f_Assign[dir2][dir1][dirActSecBeg+1];}
				            else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
					  else if ( (maskActSecBeg == TNormalBorderPoint) || (maskActSecBeg == TPreNormalBorderPoint) )
					  {
                       //! need to check acceptable normal direction

					   if (direction==directionX) {f_Assign[dirActSecBeg][dir2][dir1] = -(hx[dirActSecBeg])*rp[dirActSecBeg][dir2][dir1] + f_Assign[dirActSecBeg+1][dir2][dir1];}
	                   else if (direction==directionY) {f_Assign[dir2][dirActSecBeg][dir1] = -(hy[dirActSecBeg])*rp[dir2][dirActSecBeg][dir1] + f_Assign[dir2][dirActSecBeg+1][dir1];}
	                   else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecBeg] = -(hz[dirActSecBeg])*rp[dir2][dir1][dirActSecBeg] + f_Assign[dir2][dir1][dirActSecBeg+1];}
				       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
					  else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of begin sector.";}
                      //Sector begin (End)


					  //Sector end (Begin)
                      if ( (maskActSecEnd == TDefinedBorderPoint) || (maskActSecEnd == TEquationBorderPoint) )
					  {
					        if (direction==directionX) {f_Assign[dirActSecEnd][dir2][dir1] = rp[dirActSecEnd][dir2][dir1];}
					        else if (direction==directionY) {f_Assign[dir2][dirActSecEnd][dir1] = rp[dir2][dirActSecEnd][dir1];}
					        else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecEnd] = rp[dir2][dir1][dirActSecEnd];}
				            else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
                      else if (maskActSecEnd == TPreDefinedBorderPoint)
					  {
					        if (direction==directionX) {f_Assign[dirActSecEnd][dir2][dir1] = 2*rp[dirActSecEnd][dir2][dir1] - f_Assign[dirActSecEnd-1][dir2][dir1];}
					        else if (direction==directionY) {f_Assign[dir2][dirActSecEnd][dir1] = 2*rp[dir2][dirActSecEnd][dir1] - f_Assign[dir2][dirActSecEnd-1][dir1];}
					        else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecEnd] = 2*rp[dir2][dir1][dirActSecEnd] - f_Assign[dir2][dir1][dirActSecEnd-1];}
				            else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
					  else if ( (maskActSecEnd == TNormalBorderPoint) || (maskActSecEnd == TPreNormalBorderPoint) )
					  {
						    //! need to check acceptable normal direction
					        if (direction==directionX) {f_Assign[dirActSecEnd][dir2][dir1] = (hx[dirActSecEnd-1])*rp[dirActSecEnd][dir2][dir1] + f_Assign[dirActSecEnd-1][dir2][dir1];}
					        else if (direction==directionY) {f_Assign[dir2][dirActSecEnd][dir1] = (hy[dirActSecEnd-1])*rp[dir2][dirActSecEnd][dir1] + f_Assign[dir2][dirActSecEnd-1][dir1];}
					        else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecEnd] = (hz[dirActSecEnd-1])*rp[dir2][dir1][dirActSecEnd] + f_Assign[dir2][dir1][dirActSecEnd-1];}
				            else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
					  else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of end sector.";}
					  //Sector end (End)

					  dirActSecBeg=-1; //sector has been calculated
					  continue;
					} //if it is fake direction
					//Only borders recalculation and continue (End)





					//Need for fictive points (Begin)
					int fictiveLeft = 0;
					if ( (maskActSecBeg!=TDefinedBorderPoint) && (maskActSecBeg!=TEquationBorderPoint) )  {fictiveLeft=1;}

					int fictiveRight = 0;
					if ( (maskActSecEnd!=TDefinedBorderPoint) && (maskActSecEnd!=TEquationBorderPoint) )  {fictiveRight=1;}

					int dimActSec = dirActSecEnd-dirActSecBeg+1+fictiveLeft+fictiveRight;
					//Need for fictive points (End)


				    //Preparing scalar objects (Begin)
					int scalDim = dimActSec-1;
					double* scalVar = new double[scalDim+1];
					double* scalA = new double[scalDim+1];
					double* scalB = new double[scalDim+1];
					double* scalC = new double[scalDim+1];
					double* scalRP = new double[scalDim+1];

					scalA[0]  = 0;  scalA[scalDim]  = 0;
					scalB[0]  = 0;  scalB[scalDim]  = 0;
					scalC[0]  = 0;  scalC[scalDim]  = 0;
					scalRP[0] = 0;  scalRP[scalDim] = 0;


					//Borders points (Begin)
					if ( (maskActSecBeg == TDefinedBorderPoint) || (maskActSecBeg == TEquationBorderPoint) )
					{
						if (direction==directionX) {scalVar[0] = rp[dirActSecBeg][dir2][dir1];}
	                    else if (direction==directionY) {scalVar[0] = rp[dir2][dirActSecBeg][dir1];}
	                    else if (direction==directionZ) {scalVar[0] = rp[dir2][dir1][dirActSecBeg];}
				        else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else if  (maskActSecBeg == TPreDefinedBorderPoint)
					{
					   scalVar[0] = 0;
					   scalA[1] = 1;
					   scalB[1] = 0.5;
					   scalC[1] = 0.5;
					   if (direction==directionX) {scalRP[1] = rp[dirActSecBeg][dir2][dir1];}
	                   else if (direction==directionY) {scalRP[1] = rp[dir2][dirActSecBeg][dir1];}
	                   else if (direction==directionZ) {scalRP[1] = rp[dir2][dir1][dirActSecBeg];}
				       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else if  ( (maskActSecBeg == TNormalBorderPoint) ||  (maskActSecBeg == TPreNormalBorderPoint) )
					{
                       //! need to check acceptable normal direction
					   scalVar[0] = 0;
					   scalA[1] = 1;
					   if (direction==directionX)
					   {
						   scalB[1] = -1/hx[dirActSecBeg];
					       scalC[1] = 1/hx[dirActSecBeg];
						   scalRP[1] = rp[dirActSecBeg][dir2][dir1];
					   }
	                   else if (direction==directionY)
					   {
						   scalB[1] = -1/hy[dirActSecBeg];
					       scalC[1] = 1/hy[dirActSecBeg];
						   scalRP[1] = rp[dir2][dirActSecBeg][dir1];
					   }
	                   else if (direction==directionZ)
					   {
						   scalB[1] = -1/hz[dirActSecBeg];
					       scalC[1] = 1/hz[dirActSecBeg];
						   scalRP[1] = rp[dir2][dir1][dirActSecBeg];
					   }
				       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of begin sector.";}



					if ( (maskActSecEnd == TDefinedBorderPoint) || (maskActSecEnd == TEquationBorderPoint) )
					{
						if (direction==directionX) {scalVar[scalDim] = rp[dirActSecEnd][dir2][dir1];}
	                    else if (direction==directionY) {scalVar[scalDim] = rp[dir2][dirActSecEnd][dir1];}
	                    else if (direction==directionZ) {scalVar[scalDim] = rp[dir2][dir1][dirActSecEnd];}
				        else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else if (maskActSecEnd == TPreDefinedBorderPoint)
					{
					   scalVar[scalDim] = 0;
					   scalC[scalDim-1] = 1;
					   scalA[scalDim-1] = 0.5;
					   scalB[scalDim-1] = 0.5;
					   if (direction==directionX) {scalRP[scalDim-1] = rp[dirActSecEnd][dir2][dir1];}
					   else if (direction==directionY) {scalRP[scalDim-1] = rp[dir2][dirActSecEnd][dir1];}
					   else if (direction==directionZ) {scalRP[scalDim-1] = rp[dir2][dir1][dirActSecEnd];}
				       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else if  ( (maskActSecEnd == TNormalBorderPoint) || (maskActSecEnd == TPreNormalBorderPoint) )
					{
                       scalVar[scalDim] = 0;
					   scalC[scalDim-1] = 1;
					   if (direction==directionX)
					   {
						   scalA[scalDim-1] = -1/hx[dirActSecEnd-1];
					       scalB[scalDim-1] = 1/hx[dirActSecEnd-1];
						   scalRP[scalDim-1] = rp[dirActSecEnd][dir2][dir1];
					   }
	                   else if (direction==directionY)
					   {
						   scalA[scalDim-1] = -1/hy[dirActSecEnd-1];
					       scalB[scalDim-1] = 1/hy[dirActSecEnd-1];
						   scalRP[scalDim-1] = rp[dir2][dirActSecEnd][dir1];
					   }
	                   else if (direction==directionZ)
					   {
						   scalA[scalDim-1] = -1/hz[dirActSecEnd-1];
					       scalB[scalDim-1] = 1/hz[dirActSecEnd-1];
						   scalRP[scalDim-1] = rp[dir2][dir1][dirActSecEnd];
					   }
				       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of end sector.";}
					//Borders points (End)


					//Actual points (Begin)
					for (int scalLoop = 1+fictiveLeft;scalLoop<scalDim-fictiveRight;scalLoop++)
					{
						int dirAct = dirActSecBeg + scalLoop - fictiveLeft;

						int i; int j; int k;
						TRnRelease3dSpace* vecTrPtr;
						if (direction==directionX) {i = dirAct; j = dir2; k = dir1; vecTrPtr = &u;}
						else if (direction==directionY) {i = dir2; j = dirAct; k = dir1; vecTrPtr = &v;}
						else if (direction==directionZ) {i = dir2; j = dir1; k = dirAct; vecTrPtr = &w;}
						else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
						TRnRelease3dSpace& vecTr = *vecTrPtr;

						if (mask[i][j][k]!=TActualPoint) {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Incorrect sector innerr point.";}

						double h_1X; double hX; double hx2X; double hxhxhX;
						h_1X = hx[i-1];
					    hX = hx[i];
						hx2X = hx[i-1]+hx[i];
						hxhxhX = hx[i-1]*hx[i]*(hx[i-1]+hx[i])/2;

						double h_1Y; double hY; double hx2Y; double hxhxhY;
						h_1Y = hy[j-1];
					    hY = hy[j];
					    hx2Y = hy[j-1]+hy[j];
						hxhxhY = hy[j-1]*hy[j]*(hy[j-1]+hy[j])/2;

						double h_1Z; double hZ; double hx2Z; double hxhxhZ;
						h_1Z = hz[k-1];
						hZ = hz[k];
						hx2Z = hz[k-1]+hz[k];
						hxhxhZ = hz[k-1]*hz[k]*(hz[k-1]+hz[k])/2;

						double h_1; double h; double hx2; double hxhxh;
						if (direction==directionX) {h_1 = h_1X; h = hX; hx2 = hx2X; hxhxh = hxhxhX;}
						else if (direction==directionY) {h_1 = h_1Y; h = hY; hx2 = hx2Y; hxhxh = hxhxhY;}
						else if (direction==directionZ) {h_1 = h_1Z; h = hZ; hx2 = hx2Z; hxhxh = hxhxhZ;}
					    else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}


						if (schemeType == 0)
						{
						  scalA[scalLoop] = flagConvective*( 0.5*(-vecTr[i][j][k]/hx2) ) + flagDiffusion*( 0.5*(-(diffCoef)*h/hxhxh) );
						  scalB[scalLoop] = 1/tau + flagDiffusion*( 0.5*(diffCoef)*hx2/hxhxh );
						  scalC[scalLoop] = flagConvective*( 0.5*(vecTr[i][j][k]/hx2) ) + flagDiffusion*( 0.5*(-(diffCoef)*h_1/hxhxh) );
						}
                        else if (schemeType == 1)
						{
							double vecFlow = vecTr[i][j][k];
							if (vecFlow<0)
							{
							   scalA[scalLoop] = flagDiffusion*( 0.5*(-(diffCoef)*h/hxhxh) );
						       scalB[scalLoop] = 1/tau + flagConvective*( 0.5*(-vecTr[i][j][k]/h) ) + flagDiffusion*( 0.5*(diffCoef)*hx2/hxhxh );
							   scalC[scalLoop] = flagConvective*( 0.5*(vecTr[i][j][k]/h) ) + flagDiffusion*( 0.5*(-(diffCoef)*h_1/hxhxh) );
							}
							else
							{
								scalA[scalLoop] = flagConvective*( 0.5*(-vecTr[i][j][k]/h) ) + flagDiffusion*( 0.5*(-(diffCoef)*h/hxhxh) );
						        scalB[scalLoop] = 1/tau + flagConvective*( 0.5*(vecTr[i][j][k]/h) ) + flagDiffusion*( 0.5*(diffCoef)*hx2/hxhxh );
							    scalC[scalLoop] = flagDiffusion*( 0.5*(-(diffCoef)*h_1/hxhxh) );
							}
						}
						else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						double F = (1/tau)*f_1[i][j][k] + rp[i][j][k];

						//Help variables for simple writing F (Begin)
						double FC;
						double FD;
						TRnRelease3dSpace* fDir=NULL;
						//Help variables for simple writing F (End)

						if (direction==directionX)
					    {
							if (schemeType == 0) { F = F + flagConvectiveX*( 0.5*(-u[i][j][k])*(f_1[i+1][j][k]-f_1[i-1][j][k])/hx2X ); }
						    else if (schemeType == 1)
						    {
							   double vecFlow = u[i][j][k];
							   if (vecFlow<0) { F = F + flagConvectiveX*( 0.5*(-u[i][j][k])*(f_1[i+1][j][k]-f_1[i][j][k])/hX ); }
							   else           { F = F + flagConvectiveX*( 0.5*(-u[i][j][k])*(f_1[i][j][k]-f_1[i-1][j][k])/hX ); }
						   }
						   else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						   F = F + flagDiffusionX*(  0.5*(diffCoef)*( (hX*f_1[i-1][j][k] - hx2X*f_1[i][j][k] + h_1X*f_1[i+1][j][k])/hxhxhX )  );


						   //Y (Begin)
						   if (directionY>direction) {fDir = &f_1;}
						   else if (directionY==1) {fDir = &f_1_3;}
						   else if (directionY==2) {fDir = &f_2_3;}

						   FC = 0; FD = 0;
						   for (int item=1;item<=2;item++)
						   {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   if (schemeType == 0) { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j+1][k]-f_[i][j-1][k])/hx2Y; }
							   else if (schemeType == 1)
						       {
							      double vecFlow = v[i][j][k];
							      if (vecFlow<0) { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j+1][k]-f_[i][j][k])/hY; }
							      else           { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j][k]-f_[i][j-1][k])/hY; }
						       }
						       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						       FD = FD + 0.5*(diffCoef)*( (hY*f_[i][j-1][k] - hx2Y*f_[i][j][k] + h_1Y*f_[i][j+1][k])/hxhxhY );
						   }
						   F = F + flagConvectiveY*FC + flagDiffusionY*FD;
						   //Y (End)

						   //Z (Begin)
						   if (directionZ>direction) {fDir = &f_1;}
						   else if (directionZ==1) {fDir = &f_1_3;}
						   else if (directionZ==2) {fDir = &f_2_3;}

						   FC = 0; FD = 0;
						   for (int item=1;item<=2;item++)
						   {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   if (schemeType == 0) { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k+1]-f_[i][j][k-1])/hx2Z; }
							   else if (schemeType == 1)
						       {
							      double vecFlow = w[i][j][k];
							      if (vecFlow<0) { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k+1]-f_[i][j][k])/hZ; }
							      else           { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k]-f_[i][j][k-1])/hZ; }
						       }
						       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						       FD = FD + 0.5*(diffCoef)*( (hZ*f_[i][j][k-1] - hx2Z*f_[i][j][k] + h_1Z*f_[i][j][k+1])/hxhxhZ );
						   }
						   F = F + flagConvectiveZ*FC + flagDiffusionZ*FD;
						   //Z (End)
						}
						else if (direction==directionY)
					    {
							if (schemeType == 0) { F = F + flagConvectiveY*( 0.5*(-v[i][j][k])*(f_1[i][j+1][k]-f_1[i][j-1][k])/hx2Y ); }
						    else if (schemeType == 1)
						    {
							   double vecFlow = v[i][j][k];
							   if (vecFlow<0) { F = F + flagConvectiveY*( 0.5*(-v[i][j][k])*(f_1[i][j+1][k]-f_1[i][j][k])/hY ); }
							   else           { F = F + flagConvectiveY*( 0.5*(-v[i][j][k])*(f_1[i][j][k]-f_1[i][j-1][k])/hY ); }
						    }
						    else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}


						    F = F + flagDiffusionY*(  0.5*(diffCoef)*( (hY*f_1[i][j-1][k] - hx2Y*f_1[i][j][k] + h_1Y*f_1[i][j+1][k])/hxhxhY )  );

                            //X (Begin)
						    if (directionX>direction) {fDir = &f_1;}
						    else if (directionX==1) {fDir = &f_1_3;}
						    else if (directionX==2) {fDir = &f_2_3;}

						    FC = 0; FD = 0;
						    for (int item=1;item<=2;item++)
						    {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   if (schemeType == 0) { FC = FC + 0.5*(-u[i][j][k])*(f_[i+1][j][k]-f_[i-1][j][k])/hx2X; }
							   else if (schemeType == 1)
						       {
							      double vecFlow = u[i][j][k];
							      if (vecFlow<0) { FC = FC + 0.5*(-u[i][j][k])*(f_[i+1][j][k]-f_[i][j][k])/hX; }
							      else           { FC = FC + 0.5*(-u[i][j][k])*(f_[i][j][k]-f_[i-1][j][k])/hX; }
						       }
						       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						       FD = FD + 0.5*(diffCoef)*( (hX*f_[i-1][j][k] - hx2X*f_[i][j][k] + h_1X*f_[i+1][j][k])/hxhxhX );
						    }
						    F = F + flagConvectiveX*FC + flagDiffusionX*FD;
						    //X (End)


							//Z (Begin)
						    if (directionZ>direction) {fDir = &f_1;}
						    else if (directionZ==1) {fDir = &f_1_3;}
						    else if (directionZ==2) {fDir = &f_2_3;}

						    FC = 0; FD = 0;
						    for (int item=1;item<=2;item++)
						    {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   if (schemeType == 0) { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k+1]-f_[i][j][k-1])/hx2Z; }
							   else if (schemeType == 1)
						       {
							      double vecFlow = w[i][j][k];
							      if (vecFlow<0) { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k+1]-f_[i][j][k])/hZ; }
							      else           { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k]-f_[i][j][k-1])/hZ; }
						       }
						       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						       FD = FD + 0.5*(diffCoef)*( (hZ*f_[i][j][k-1] - hx2Z*f_[i][j][k] + h_1Z*f_[i][j][k+1])/hxhxhZ );
						    }
						    F = F + flagConvectiveZ*FC + flagDiffusionZ*FD;
						    //Z (End)
					    }
	                    else if (direction==directionZ)
					    {
							if (schemeType == 0) { F = F + flagConvectiveZ*( 0.5*(-w[i][j][k])*(f_1[i][j][k+1]-f_1[i][j][k-1])/hx2Z ); }
						    else if (schemeType == 1)
						    {
							   double vecFlow = w[i][j][k];
							   if (vecFlow<0) { F = F + flagConvectiveZ*( 0.5*(-w[i][j][k])*(f_1[i][j][k+1]-f_1[i][j][k])/hZ ); }
							   else           { F = F + flagConvectiveZ*( 0.5*(-w[i][j][k])*(f_1[i][j][k]-f_1[i][j][k-1])/hZ ); }
						    }
						    else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						    F = F + flagDiffusionZ*(  0.5*(diffCoef)*( (hZ*f_1[i][j][k-1] - hx2Z*f_1[i][j][k] + h_1Z*f_1[i][j][k+1])/hxhxhZ )  );

						    //X (Begin)
						    if (directionX>direction) {fDir = &f_1;}
						    else if (directionX==1) {fDir = &f_1_3;}
						    else if (directionX==2) {fDir = &f_2_3;}

						    FC = 0; FD = 0;
						    for (int item=1;item<=2;item++)
						    {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   if (schemeType == 0) { FC = FC + 0.5*(-u[i][j][k])*(f_[i+1][j][k]-f_[i-1][j][k])/hx2X; }
							   else if (schemeType == 1)
						       {
							      double vecFlow = u[i][j][k];
							      if (vecFlow<0) { FC = FC + 0.5*(-u[i][j][k])*(f_[i+1][j][k]-f_[i][j][k])/hX; }
							      else           { FC = FC + 0.5*(-u[i][j][k])*(f_[i][j][k]-f_[i-1][j][k])/hX; }
						       }
						       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						       FD = FD + 0.5*(diffCoef)*( (hX*f_[i-1][j][k] - hx2X*f_[i][j][k] + h_1X*f_[i+1][j][k])/hxhxhX );
						    }
						    F = F + flagConvectiveX*FC + flagDiffusionX*FD;
						    //X (End)

							//Y (Begin)
						    if (directionY>direction) {fDir = &f_1;}
						    else if (directionY==1) {fDir = &f_1_3;}
						    else if (directionY==2) {fDir = &f_2_3;}

						    FC = 0; FD = 0;
						    for (int item=1;item<=2;item++)
						    {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   if (schemeType == 0) { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j+1][k]-f_[i][j-1][k])/hx2Y; }
							   else if (schemeType == 1)
						       {
							      double vecFlow = v[i][j][k];
							      if (vecFlow<0) { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j+1][k]-f_[i][j][k])/hY; }
							      else           { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j][k]-f_[i][j-1][k])/hY; }
						       }
						       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						       FD = FD + 0.5*(diffCoef)*( (hY*f_[i][j-1][k] - hx2Y*f_[i][j][k] + h_1Y*f_[i][j+1][k])/hxhxhY );
						    }
						    F = F + flagConvectiveY*FC + flagDiffusionY*FD;
						    //Y (End)
					    }
				        else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
						scalRP[scalLoop] = F;
					} // by act direction
					//Actual points (End)

					//Preparing scalar objects (End)


					//Scalar running (Begin)
					TSLAE::SolveByScalarRunning(scalDim,scalVar,scalA,scalB,scalC,scalRP);
					//Scalar running (End)


					//Assigninig (Begin)
					TRnRelease3dSpace* f_AssignPtr;
					if (direction==1) {f_AssignPtr = &f_1_3;}
				    else if (direction==2) {f_AssignPtr = &f_2_3;}
				    else if (direction==3) {f_AssignPtr = &f;}
					else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					TRnRelease3dSpace& f_Assign = *f_AssignPtr;

					for (int scalLoop = 0+fictiveLeft;scalLoop<scalDim+1-fictiveRight;scalLoop++)
					{
						int dirAct = dirActSecBeg + scalLoop - fictiveLeft;
						int i; int j; int k;
						if (direction==directionX) {i = dirAct; j = dir2; k = dir1;}
						else if (direction==directionY) {i = dir2; j = dirAct; k = dir1;}
						else if (direction==directionZ) {i = dir2; j = dir1; k = dirAct;}
						else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}

						f_Assign[i][j][k] = scalVar[scalLoop];
					}
					//Assigninig (End)

					//Deleting scalar objects (Begin)
					delete[] scalVar;// scalVar = NULL;
					delete[] scalA;// scalA = NULL;
					delete[] scalB;// scalB = NULL;
					delete[] scalC;// scalC = NULL;
					delete[] scalRP;// scalRP = NULL;
					//Deleting scalar objects (End)

                    dirActSecBeg=-1; //sector has been calculated
                 } // sector calculating
           } // by sectors of act direction
         } // by direction 2
       } // by direction 1
    } // by fake direction (x,y,z) - for border conditions recalculating
	} // by global direction (x,y,z)
	delete &f_1_3;
	delete &f_2_3;
}
//*************************************
// old one
/*
void THydrodynamicsProblem::SolveTransferEquation__(TRnRelease3dSpace& lVecCur, TRnRelease3dSpace& lVecPrev,
		                               TRnRelease3dSpace&  lU, TRnRelease3dSpace& lV, TRnRelease3dSpace& lW, TRnRelease3dSpace& lRP,
								       const TNodesGrid& lGrid, TRnRelease3dSpace& lDiffusionCoefficient, int lSchemeType)
{
    //Preparing variables (Begin)
	TRnRelease3dSpace& f = lVecCur;
    TRnRelease3dSpace& f_1 = lVecPrev;
    TRnRelease3dSpace& u = lU;
	TRnRelease3dSpace& v = lV;
	TRnRelease3dSpace& w = lW;
	TRnRelease3dSpace& rp = lRP;

	TRnRelease3dSpace& f_1_3 = dynamic_cast<TRnRelease3dSpace&>(f_1.Copy());
	TRnRelease3dSpace& f_2_3 = dynamic_cast<TRnRelease3dSpace&>(f_1.Copy());

	const T3dNormalGrid& grid = dynamic_cast<const T3dNormalGrid&>(lGrid);

	const T3dNumberMask& mask = grid.Mask;

	const TSeparator& sepX = grid.Separator1;
	const TSeparator& sepY = grid.Separator2;
	const TSeparator& sepZ = grid.Separator3;

	int N = sepX.EndIndex;
    int L = sepY.EndIndex;
    int M = sepZ.EndIndex;

	const double* hx = sepX.GetSeparation();
	const double* hy = sepY.GetSeparation();
	const double* hz = sepZ.GetSeparation();

	double tau = fTimeSeparator->SeparationValue;

	TRnRelease3dSpace& diffCoefVar = lDiffusionCoefficient;

	//Preparing variables (End)

	//Direction order (Begin)
	int directionOrder = 1;
	int directionX;
	int directionY;
	int directionZ;
	if (directionOrder == 1) {directionX=1;directionY=2;directionZ=3;}
	else if (directionOrder == 2) {directionX=1;directionY=3;directionZ=2;}
	else if (directionOrder == 3) {directionX=2;directionY=1;directionZ=3;}
	else if (directionOrder == 4) {directionX=2;directionY=3;directionZ=1;}
	else if (directionOrder == 5) {directionX=3;directionY=1;directionZ=2;}
	else if (directionOrder == 6) {directionX=3;directionY=2;directionZ=1;}
	else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction order.";}
	//Direction order (End)

	//Terms using (Begin)
	int flagConvectiveX = 1; int flagDiffusionX = 1;
    int flagConvectiveY = 1; int flagDiffusionY = 1;
    int flagConvectiveZ = 1; int flagDiffusionZ = 1;
	//Terms using (Begin)

	//Scheme type (Begin)
	 //0 - central;
	 //1 - upstream;
	 int schemeType = lSchemeType;
	//Scheme type (End)

	for (int directionGlobal=1;directionGlobal<4;directionGlobal++) //splitting scheame steps
	{
	for (int directionFake=1;directionFake<4;directionFake++) //for border condtions recalculating
    {
	   int direction;
	   if (directionFake==1) {direction = directionGlobal;}
	   else
	   {
		   if (directionGlobal==1) {direction=directionFake;}
		   else if (directionGlobal==2)
		   {
			   if (directionFake==2) {direction=1;}
			   else if (directionFake==3) {direction=3;}
			   else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown fake direction.";}
		   }
		   else if (directionGlobal==3) {direction=directionFake-1;}
		   else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown global direction.";}
	   }


	   int dir1EndIndex;
	   int dir2EndIndex;
	   int dirActEndIndex;

	   if (direction==directionX)
	   {
		 dir1EndIndex = M;
	     dir2EndIndex = L;
	     dirActEndIndex = N;
	   }
	   else if (direction==directionY)
	   {
         dir1EndIndex = M;
	     dir2EndIndex = N;
	     dirActEndIndex = L;
	   }
	   else if (direction==directionZ)
	   {
	     dir1EndIndex = L;
	     dir2EndIndex = N;
	     dirActEndIndex = M;
	   }
	   else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}

	   int flagConvective = 0;
	   if (
		    ( (direction==directionX) && (flagConvectiveX==1) ) ||
		    ( (direction==directionY) && (flagConvectiveY==1) ) ||
			( (direction==directionZ) && (flagConvectiveZ==1) )
		  ) {flagConvective=1;}

	   int flagDiffusion = 0;
	   if (
		    ( (direction==directionX) && (flagDiffusionX==1) ) ||
		    ( (direction==directionY) && (flagDiffusionY==1) ) ||
			( (direction==directionZ) && (flagDiffusionZ==1) )
		  ) {flagDiffusion=1;}


	   for (int dir1=0;dir1<dir1EndIndex+1;dir1++)
	   {
         for (int dir2=0;dir2<dir2EndIndex+1;dir2++)
         {
           int dirActSecBeg = -1;
		   int maskActSecBeg;

           int dirActSecEnd;
		   int maskActSecEnd;

           for (int dirActSec=0; dirActSec<dirActEndIndex+1; dirActSec++)
           {
			     int maskActSec;
                 if (direction==directionX) {maskActSec = mask[dirActSec][dir2][dir1];}
	             else if (direction==directionY) {maskActSec = mask[dir2][dirActSec][dir1];}
	             else if (direction==directionZ) {maskActSec = mask[dir2][dir1][dirActSec];}
				 else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}

                 if (dirActSecBeg<0)
                 {
                   if (maskActSec==TActualPoint)
                   {
                      if (dirActSec==0)
                      {
                        throw "Error in THydrodynamicsProblem::SolveTransferEquation: Initial border not found.";
                      }
                      else
                      {
                        dirActSecBeg = dirActSec-1;
						if (direction==directionX) {maskActSecBeg = mask[dirActSecBeg][dir2][dir1];}
	                    else if (direction==directionY) {maskActSecBeg = mask[dir2][dirActSecBeg][dir1];}
	                    else if (direction==directionZ) {maskActSecBeg = mask[dir2][dir1][dirActSecBeg];}
				        else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
                      }
                   }
                 }
                 else if (maskActSec!=TActualPoint)
                 {
                    dirActSecEnd=dirActSec;
                    maskActSecEnd = maskActSec;

					//Сhecking sector (Begin)
					if ( (maskActSecBeg == TActualPoint) || (maskActSecBeg == TFictivePoint) || (maskActSecEnd == TActualPoint) || (maskActSecEnd == TFictivePoint) )
					{
                        throw "Error in THydrodynamicsProblem::SolveTransferEquation: Incorrect sector border.";
					}

					if ( (dirActSecBeg<0) || (dirActSecEnd-dirActSecBeg)<2 )
					{
                        throw "Error in THydrodynamicsProblem::SolveTransferEquation: Incorrect sector size.";
					}
					//Сhecking sector (End)



					//Only borders recalculation and continue (Begin)
					if (direction!=directionGlobal)
					{
                      TRnRelease3dSpace* f_AssignPtr;
					  if (directionGlobal==1) {f_AssignPtr = &f_1_3;}
				      else if (directionGlobal==2) {f_AssignPtr = &f_2_3;}
				      else if (directionGlobal==3) {f_AssignPtr = &f;}
					  else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown global direction.";}
					  TRnRelease3dSpace& f_Assign = *f_AssignPtr;

					  //Sector begin (Begin)
					  if  ( (maskActSecBeg == TDefinedBorderPoint) || (maskActSecBeg == TEquationBorderPoint) )
				   	  {
					        if (direction==directionX) {f_Assign[dirActSecBeg][dir2][dir1] = rp[dirActSecBeg][dir2][dir1];}
	                        else if (direction==directionY) {f_Assign[dir2][dirActSecBeg][dir1] = rp[dir2][dirActSecBeg][dir1];}
	                        else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecBeg] = rp[dir2][dir1][dirActSecBeg];}
				            else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
					  else if  (maskActSecBeg == TPreDefinedBorderPoint)
				   	  {
					        if (direction==directionX) {f_Assign[dirActSecBeg][dir2][dir1] = 2*rp[dirActSecBeg][dir2][dir1] - f_Assign[dirActSecBeg+1][dir2][dir1];}
	                        else if (direction==directionY) {f_Assign[dir2][dirActSecBeg][dir1] = 2*rp[dir2][dirActSecBeg][dir1] - f_Assign[dir2][dirActSecBeg+1][dir1];}
	                        else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecBeg] = 2*rp[dir2][dir1][dirActSecBeg] - f_Assign[dir2][dir1][dirActSecBeg+1];}
				            else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
					  else if ( (maskActSecBeg == TNormalBorderPoint) || (maskActSecBeg == TPreNormalBorderPoint) )
					  {
                       //! need to check acceptable normal direction

					   if (direction==directionX) {f_Assign[dirActSecBeg][dir2][dir1] = -(hx[dirActSecBeg])*rp[dirActSecBeg][dir2][dir1] + f_Assign[dirActSecBeg+1][dir2][dir1];}
	                   else if (direction==directionY) {f_Assign[dir2][dirActSecBeg][dir1] = -(hy[dirActSecBeg])*rp[dir2][dirActSecBeg][dir1] + f_Assign[dir2][dirActSecBeg+1][dir1];}
	                   else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecBeg] = -(hz[dirActSecBeg])*rp[dir2][dir1][dirActSecBeg] + f_Assign[dir2][dir1][dirActSecBeg+1];}
				       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
					  else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of begin sector.";}
                      //Sector begin (End)


					  //Sector end (Begin)
                      if ( (maskActSecEnd == TDefinedBorderPoint) || (maskActSecEnd == TEquationBorderPoint) )
					  {
					        if (direction==directionX) {f_Assign[dirActSecEnd][dir2][dir1] = rp[dirActSecEnd][dir2][dir1];}
					        else if (direction==directionY) {f_Assign[dir2][dirActSecEnd][dir1] = rp[dir2][dirActSecEnd][dir1];}
					        else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecEnd] = rp[dir2][dir1][dirActSecEnd];}
				            else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
                      else if (maskActSecEnd == TPreDefinedBorderPoint)
					  {
					        if (direction==directionX) {f_Assign[dirActSecEnd][dir2][dir1] = 2*rp[dirActSecEnd][dir2][dir1] - f_Assign[dirActSecEnd-1][dir2][dir1];}
					        else if (direction==directionY) {f_Assign[dir2][dirActSecEnd][dir1] = 2*rp[dir2][dirActSecEnd][dir1] - f_Assign[dir2][dirActSecEnd-1][dir1];}
					        else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecEnd] = 2*rp[dir2][dir1][dirActSecEnd] - f_Assign[dir2][dir1][dirActSecEnd-1];}
				            else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
					  else if ( (maskActSecEnd == TNormalBorderPoint) || (maskActSecEnd == TPreNormalBorderPoint) )
					  {
						    //! need to check acceptable normal direction

					        if (direction==directionX) {f_Assign[dirActSecEnd][dir2][dir1] = (hx[dirActSecEnd-1])*rp[dirActSecEnd][dir2][dir1] + f_Assign[dirActSecEnd-1][dir2][dir1];}
					        else if (direction==directionY) {f_Assign[dir2][dirActSecEnd][dir1] = (hy[dirActSecEnd-1])*rp[dir2][dirActSecEnd][dir1] + f_Assign[dir2][dirActSecEnd-1][dir1];}
					        else if (direction==directionZ) {f_Assign[dir2][dir1][dirActSecEnd] = (hz[dirActSecEnd-1])*rp[dir2][dir1][dirActSecEnd] + f_Assign[dir2][dir1][dirActSecEnd-1];}
				            else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					  }
					  else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of end sector.";}
					  //Sector end (End)

					  dirActSecBeg=-1; //sector has been calculated
					  continue;
					} //if it is fake direction
					//Only borders recalculation and continue (End)





					//Need for fictive points (Begin)
					int fictiveLeft = 0;
					if ( (maskActSecBeg!=TDefinedBorderPoint) && (maskActSecBeg!=TEquationBorderPoint) )  {fictiveLeft=1;}

					int fictiveRight = 0;
					if ( (maskActSecEnd!=TDefinedBorderPoint) && (maskActSecEnd!=TEquationBorderPoint) )  {fictiveRight=1;}

					int dimActSec = dirActSecEnd-dirActSecBeg+1+fictiveLeft+fictiveRight;
					//Need for fictive points (End)


				    //Preparing scalar objects (Begin)
					int scalDim = dimActSec-1;
					double* scalVar = new double[scalDim+1];
					double* scalA = new double[scalDim+1];
					double* scalB = new double[scalDim+1];
					double* scalC = new double[scalDim+1];
					double* scalRP = new double[scalDim+1];

					scalA[0]  = 0;  scalA[scalDim]  = 0;
					scalB[0]  = 0;  scalB[scalDim]  = 0;
					scalC[0]  = 0;  scalC[scalDim]  = 0;
					scalRP[0] = 0;  scalRP[scalDim] = 0;


					//Borders points (Begin)
					if ( (maskActSecBeg == TDefinedBorderPoint) || (maskActSecBeg == TEquationBorderPoint) )
					{
						if (direction==directionX) {scalVar[0] = rp[dirActSecBeg][dir2][dir1];}
	                    else if (direction==directionY) {scalVar[0] = rp[dir2][dirActSecBeg][dir1];}
	                    else if (direction==directionZ) {scalVar[0] = rp[dir2][dir1][dirActSecBeg];}
				        else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else if  (maskActSecBeg == TPreDefinedBorderPoint)
					{
					   scalVar[0] = 0;
					   scalA[1] = 1;
					   scalB[1] = 0.5;
					   scalC[1] = 0.5;
					   if (direction==directionX) {scalRP[1] = rp[dirActSecBeg][dir2][dir1];}
	                   else if (direction==directionY) {scalRP[1] = rp[dir2][dirActSecBeg][dir1];}
	                   else if (direction==directionZ) {scalRP[1] = rp[dir2][dir1][dirActSecBeg];}
				       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else if  ( (maskActSecBeg == TNormalBorderPoint) ||  (maskActSecBeg == TPreNormalBorderPoint) )
					{
                       //! need to check acceptable normal direction
					   scalVar[0] = 0;
					   scalA[1] = 1;
					   if (direction==directionX)
					   {
						   scalB[1] = -1/hx[dirActSecBeg];
					       scalC[1] = 1/hx[dirActSecBeg];
						   scalRP[1] = rp[dirActSecBeg][dir2][dir1];
					   }
	                   else if (direction==directionY)
					   {
						   scalB[1] = -1/hy[dirActSecBeg];
					       scalC[1] = 1/hy[dirActSecBeg];
						   scalRP[1] = rp[dir2][dirActSecBeg][dir1];
					   }
	                   else if (direction==directionZ)
					   {
						   scalB[1] = -1/hz[dirActSecBeg];
					       scalC[1] = 1/hz[dirActSecBeg];
						   scalRP[1] = rp[dir2][dir1][dirActSecBeg];
					   }
				       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of begin sector.";}



					if ( (maskActSecEnd == TDefinedBorderPoint) || (maskActSecEnd == TEquationBorderPoint) )
					{
						if (direction==directionX) {scalVar[scalDim] = rp[dirActSecEnd][dir2][dir1];}
	                    else if (direction==directionY) {scalVar[scalDim] = rp[dir2][dirActSecEnd][dir1];}
	                    else if (direction==directionZ) {scalVar[scalDim] = rp[dir2][dir1][dirActSecEnd];}
				        else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else if (maskActSecEnd == TPreDefinedBorderPoint)
					{
					   scalVar[scalDim] = 0;
					   scalC[scalDim-1] = 1;
					   scalA[scalDim-1] = 0.5;
					   scalB[scalDim-1] = 0.5;
					   if (direction==directionX) {scalRP[scalDim-1] = rp[dirActSecEnd][dir2][dir1];}
					   else if (direction==directionY) {scalRP[scalDim-1] = rp[dir2][dirActSecEnd][dir1];}
					   else if (direction==directionZ) {scalRP[scalDim-1] = rp[dir2][dir1][dirActSecEnd];}
				       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else if  ( (maskActSecEnd == TNormalBorderPoint) || (maskActSecEnd == TPreNormalBorderPoint) )
					{
                       scalVar[scalDim] = 0;
					   scalC[scalDim-1] = 1;
					   if (direction==directionX)
					   {
						   scalA[scalDim-1] = -1/hx[dirActSecEnd-1];
					       scalB[scalDim-1] = 1/hx[dirActSecEnd-1];
						   scalRP[scalDim-1] = rp[dirActSecEnd][dir2][dir1];
					   }
	                   else if (direction==directionY)
					   {
						   scalA[scalDim-1] = -1/hy[dirActSecEnd-1];
					       scalB[scalDim-1] = 1/hy[dirActSecEnd-1];
						   scalRP[scalDim-1] = rp[dir2][dirActSecEnd][dir1];
					   }
	                   else if (direction==directionZ)
					   {
						   scalA[scalDim-1] = -1/hz[dirActSecEnd-1];
					       scalB[scalDim-1] = 1/hz[dirActSecEnd-1];
						   scalRP[scalDim-1] = rp[dir2][dir1][dirActSecEnd];
					   }
				       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					}
					else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of end sector.";}
					//Borders points (End)


					//Actual points (Begin)
					for (int scalLoop = 1+fictiveLeft;scalLoop<scalDim-fictiveRight;scalLoop++)
					{
						int dirAct = dirActSecBeg + scalLoop - fictiveLeft;

						int i; int j; int k;
						TRnRelease3dSpace* vecTrPtr;
						if (direction==directionX) {i = dirAct; j = dir2; k = dir1; vecTrPtr = &u;}
						else if (direction==directionY) {i = dir2; j = dirAct; k = dir1; vecTrPtr = &v;}
						else if (direction==directionZ) {i = dir2; j = dir1; k = dirAct; vecTrPtr = &w;}
						else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
						TRnRelease3dSpace& vecTr = *vecTrPtr;

						if (mask[i][j][k]!=TActualPoint) {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Incorrect sector innerr point.";}

						double h_1X; double hX; double hx2X; double hxhxhX;
						h_1X = hx[i-1];
					    hX = hx[i];
						hx2X = hx[i-1]+hx[i];
						hxhxhX = hx[i-1]*hx[i]*(hx[i-1]+hx[i])/2;

						double h_1Y; double hY; double hx2Y; double hxhxhY;
						h_1Y = hy[j-1];
					    hY = hy[j];
					    hx2Y = hy[j-1]+hy[j];
						hxhxhY = hy[j-1]*hy[j]*(hy[j-1]+hy[j])/2;

						double h_1Z; double hZ; double hx2Z; double hxhxhZ;
						h_1Z = hz[k-1];
						hZ = hz[k];
						hx2Z = hz[k-1]+hz[k];
						hxhxhZ = hz[k-1]*hz[k]*(hz[k-1]+hz[k])/2;

						double h_1; double h; double hx2; double hxhxh;
						if (direction==directionX) {h_1 = h_1X; h = hX; hx2 = hx2X; hxhxh = hxhxhX;}
						else if (direction==directionY) {h_1 = h_1Y; h = hY; hx2 = hx2Y; hxhxh = hxhxhY;}
						else if (direction==directionZ) {h_1 = h_1Z; h = hZ; hx2 = hx2Z; hxhxh = hxhxhZ;}
					    else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}


						double diffCoef = diffCoefVar[i][j][k];

						if (schemeType == 0)
						{
						  scalA[scalLoop] = flagConvective*( 0.5*(-vecTr[i][j][k]/hx2) ) + flagDiffusion*( 0.5*(-(diffCoef)*h/hxhxh) );
						  scalB[scalLoop] = 1/tau + flagDiffusion*( 0.5*(diffCoef)*hx2/hxhxh );
						  scalC[scalLoop] = flagConvective*( 0.5*(vecTr[i][j][k]/hx2) ) + flagDiffusion*( 0.5*(-(diffCoef)*h_1/hxhxh) );
						}
                        else if (schemeType == 1)
						{
							double vecFlow = vecTr[i][j][k];
							if (vecFlow<0)
							{
							   scalA[scalLoop] = flagDiffusion*( 0.5*(-(diffCoef)*h/hxhxh) );
						       scalB[scalLoop] = 1/tau + flagConvective*( 0.5*(-vecTr[i][j][k]/h) ) + flagDiffusion*( 0.5*(diffCoef)*hx2/hxhxh );
							   scalC[scalLoop] = flagConvective*( 0.5*(vecTr[i][j][k]/h) ) + flagDiffusion*( 0.5*(-(diffCoef)*h_1/hxhxh) );
							}
							else
							{
								scalA[scalLoop] = flagConvective*( 0.5*(-vecTr[i][j][k]/h) ) + flagDiffusion*( 0.5*(-(diffCoef)*h/hxhxh) );
						        scalB[scalLoop] = 1/tau + flagConvective*( 0.5*(vecTr[i][j][k]/h) ) + flagDiffusion*( 0.5*(diffCoef)*hx2/hxhxh );
							    scalC[scalLoop] = flagDiffusion*( 0.5*(-(diffCoef)*h_1/hxhxh) );
							}
						}
						else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						double F = (1/tau)*f_1[i][j][k] + rp[i][j][k];

						//Help variables for simple writing F (Begin)
						double FC;
						double FD;
						TRnRelease3dSpace* fDir=NULL;
						//Help variables for simple writing F (End)

						if (direction==directionX)
					    {
							if (schemeType == 0) { F = F + flagConvectiveX*( 0.5*(-u[i][j][k])*(f_1[i+1][j][k]-f_1[i-1][j][k])/hx2X ); }
						    else if (schemeType == 1)
						    {
							   double vecFlow = u[i][j][k];
							   if (vecFlow<0) { F = F + flagConvectiveX*( 0.5*(-u[i][j][k])*(f_1[i+1][j][k]-f_1[i][j][k])/hX ); }
							   else           { F = F + flagConvectiveX*( 0.5*(-u[i][j][k])*(f_1[i][j][k]-f_1[i-1][j][k])/hX ); }
						   }
						   else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						   F = F + flagDiffusionX*(  0.5*(diffCoef)*( (hX*f_1[i-1][j][k] - hx2X*f_1[i][j][k] + h_1X*f_1[i+1][j][k])/hxhxhX )  );


						   //Y (Begin)
						   if (directionY>direction) {fDir = &f_1;}
						   else if (directionY==1) {fDir = &f_1_3;}
						   else if (directionY==2) {fDir = &f_2_3;}

						   FC = 0; FD = 0;
						   for (int item=1;item<=2;item++)
						   {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   if (schemeType == 0) { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j+1][k]-f_[i][j-1][k])/hx2Y; }
							   else if (schemeType == 1)
						       {
							      double vecFlow = v[i][j][k];
							      if (vecFlow<0) { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j+1][k]-f_[i][j][k])/hY; }
							      else           { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j][k]-f_[i][j-1][k])/hY; }
						       }
						       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						       FD = FD + 0.5*(diffCoef)*( (hY*f_[i][j-1][k] - hx2Y*f_[i][j][k] + h_1Y*f_[i][j+1][k])/hxhxhY );
						   }
						   F = F + flagConvectiveY*FC + flagDiffusionY*FD;
						   //Y (End)

						   //Z (Begin)
						   if (directionZ>direction) {fDir = &f_1;}
						   else if (directionZ==1) {fDir = &f_1_3;}
						   else if (directionZ==2) {fDir = &f_2_3;}

						   FC = 0; FD = 0;
						   for (int item=1;item<=2;item++)
						   {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   if (schemeType == 0) { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k+1]-f_[i][j][k-1])/hx2Z; }
							   else if (schemeType == 1)
						       {
							      double vecFlow = w[i][j][k];
							      if (vecFlow<0) { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k+1]-f_[i][j][k])/hZ; }
							      else           { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k]-f_[i][j][k-1])/hZ; }
						       }
						       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						       FD = FD + 0.5*(diffCoef)*( (hZ*f_[i][j][k-1] - hx2Z*f_[i][j][k] + h_1Z*f_[i][j][k+1])/hxhxhZ );
						   }
						   F = F + flagConvectiveZ*FC + flagDiffusionZ*FD;
						   //Z (End)
						}
						else if (direction==directionY)
					    {
							if (schemeType == 0) { F = F + flagConvectiveY*( 0.5*(-v[i][j][k])*(f_1[i][j+1][k]-f_1[i][j-1][k])/hx2Y ); }
						    else if (schemeType == 1)
						    {
							   double vecFlow = v[i][j][k];
							   if (vecFlow<0) { F = F + flagConvectiveY*( 0.5*(-v[i][j][k])*(f_1[i][j+1][k]-f_1[i][j][k])/hY ); }
							   else           { F = F + flagConvectiveY*( 0.5*(-v[i][j][k])*(f_1[i][j][k]-f_1[i][j-1][k])/hY ); }
						    }
						    else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}


						    F = F + flagDiffusionY*(  0.5*(diffCoef)*( (hY*f_1[i][j-1][k] - hx2Y*f_1[i][j][k] + h_1Y*f_1[i][j+1][k])/hxhxhY )  );

                            //X (Begin)
						    if (directionX>direction) {fDir = &f_1;}
						    else if (directionX==1) {fDir = &f_1_3;}
						    else if (directionX==2) {fDir = &f_2_3;}

						    FC = 0; FD = 0;
						    for (int item=1;item<=2;item++)
						    {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   if (schemeType == 0) { FC = FC + 0.5*(-u[i][j][k])*(f_[i+1][j][k]-f_[i-1][j][k])/hx2X; }
							   else if (schemeType == 1)
						       {
							      double vecFlow = u[i][j][k];
							      if (vecFlow<0) { FC = FC + 0.5*(-u[i][j][k])*(f_[i+1][j][k]-f_[i][j][k])/hX; }
							      else           { FC = FC + 0.5*(-u[i][j][k])*(f_[i][j][k]-f_[i-1][j][k])/hX; }
						       }
						       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						       FD = FD + 0.5*(diffCoef)*( (hX*f_[i-1][j][k] - hx2X*f_[i][j][k] + h_1X*f_[i+1][j][k])/hxhxhX );
						    }
						    F = F + flagConvectiveX*FC + flagDiffusionX*FD;
						    //X (End)


							//Z (Begin)
						    if (directionZ>direction) {fDir = &f_1;}
						    else if (directionZ==1) {fDir = &f_1_3;}
						    else if (directionZ==2) {fDir = &f_2_3;}

						    FC = 0; FD = 0;
						    for (int item=1;item<=2;item++)
						    {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   if (schemeType == 0) { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k+1]-f_[i][j][k-1])/hx2Z; }
							   else if (schemeType == 1)
						       {
							      double vecFlow = w[i][j][k];
							      if (vecFlow<0) { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k+1]-f_[i][j][k])/hZ; }
							      else           { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k]-f_[i][j][k-1])/hZ; }
						       }
						       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						       FD = FD + 0.5*(diffCoef)*( (hZ*f_[i][j][k-1] - hx2Z*f_[i][j][k] + h_1Z*f_[i][j][k+1])/hxhxhZ );
						    }
						    F = F + flagConvectiveZ*FC + flagDiffusionZ*FD;
						    //Z (End)
					    }
	                    else if (direction==directionZ)
					    {
							if (schemeType == 0) { F = F + flagConvectiveZ*( 0.5*(-w[i][j][k])*(f_1[i][j][k+1]-f_1[i][j][k-1])/hx2Z ); }
						    else if (schemeType == 1)
						    {
							   double vecFlow = w[i][j][k];
							   if (vecFlow<0) { F = F + flagConvectiveZ*( 0.5*(-w[i][j][k])*(f_1[i][j][k+1]-f_1[i][j][k])/hZ ); }
							   else           { F = F + flagConvectiveZ*( 0.5*(-w[i][j][k])*(f_1[i][j][k]-f_1[i][j][k-1])/hZ ); }
						    }
						    else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						    F = F + flagDiffusionZ*(  0.5*(diffCoef)*( (hZ*f_1[i][j][k-1] - hx2Z*f_1[i][j][k] + h_1Z*f_1[i][j][k+1])/hxhxhZ )  );

						    //X (Begin)
						    if (directionX>direction) {fDir = &f_1;}
						    else if (directionX==1) {fDir = &f_1_3;}
						    else if (directionX==2) {fDir = &f_2_3;}

						    FC = 0; FD = 0;
						    for (int item=1;item<=2;item++)
						    {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   if (schemeType == 0) { FC = FC + 0.5*(-u[i][j][k])*(f_[i+1][j][k]-f_[i-1][j][k])/hx2X; }
							   else if (schemeType == 1)
						       {
							      double vecFlow = u[i][j][k];
							      if (vecFlow<0) { FC = FC + 0.5*(-u[i][j][k])*(f_[i+1][j][k]-f_[i][j][k])/hX; }
							      else           { FC = FC + 0.5*(-u[i][j][k])*(f_[i][j][k]-f_[i-1][j][k])/hX; }
						       }
						       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						       FD = FD + 0.5*(diffCoef)*( (hX*f_[i-1][j][k] - hx2X*f_[i][j][k] + h_1X*f_[i+1][j][k])/hxhxhX );
						    }
						    F = F + flagConvectiveX*FC + flagDiffusionX*FD;
						    //X (End)

							//Y (Begin)
						    if (directionY>direction) {fDir = &f_1;}
						    else if (directionY==1) {fDir = &f_1_3;}
						    else if (directionY==2) {fDir = &f_2_3;}

						    FC = 0; FD = 0;
						    for (int item=1;item<=2;item++)
						    {
							   TRnRelease3dSpace* f_Ptr = fDir; if (item==2) {f_Ptr = &f_1;}
							   TRnRelease3dSpace& f_ = *f_Ptr;
							   if (schemeType == 0) { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j+1][k]-f_[i][j-1][k])/hx2Y; }
							   else if (schemeType == 1)
						       {
							      double vecFlow = v[i][j][k];
							      if (vecFlow<0) { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j+1][k]-f_[i][j][k])/hY; }
							      else           { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j][k]-f_[i][j-1][k])/hY; }
						       }
						       else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type.";}

						       FD = FD + 0.5*(diffCoef)*( (hY*f_[i][j-1][k] - hx2Y*f_[i][j][k] + h_1Y*f_[i][j+1][k])/hxhxhY );
						    }
						    F = F + flagConvectiveY*FC + flagDiffusionY*FD;
						    //Y (End)
					    }
				        else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
						scalRP[scalLoop] = F;
					} // by act direction
					//Actual points (End)

					//Preparing scalar objects (End)


					//Scalar running (Begin)
					TSLAE::SolveByScalarRunning(scalDim,scalVar,scalA,scalB,scalC,scalRP);
					//Scalar running (End)


					//Assigninig (Begin)
					TRnRelease3dSpace* f_AssignPtr;
					if (direction==1) {f_AssignPtr = &f_1_3;}
				    else if (direction==2) {f_AssignPtr = &f_2_3;}
				    else if (direction==3) {f_AssignPtr = &f;}
					else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}
					TRnRelease3dSpace& f_Assign = *f_AssignPtr;

					for (int scalLoop = 0+fictiveLeft;scalLoop<scalDim+1-fictiveRight;scalLoop++)
					{
						int dirAct = dirActSecBeg + scalLoop - fictiveLeft;
						int i; int j; int k;
						if (direction==directionX) {i = dirAct; j = dir2; k = dir1;}
						else if (direction==directionY) {i = dir2; j = dirAct; k = dir1;}
						else if (direction==directionZ) {i = dir2; j = dir1; k = dirAct;}
						else {throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction.";}

						f_Assign[i][j][k] = scalVar[scalLoop];
					}
					//Assigninig (End)

					//Deleting scalar objects (Begin)
					delete[] scalVar; scalVar = NULL;
					delete[] scalA; scalA = NULL;
					delete[] scalB; scalB = NULL;
					delete[] scalC; scalC = NULL;
					delete[] scalRP; scalRP = NULL;
					//Deleting scalar objects (End)

                    dirActSecBeg=-1; //sector has been calculated
                 } // sector calculating
           } // by sectors of act direction
         } // by direction 2
       } // by direction 1
    } // by fake direction (x,y,z) - for border conditions recalculating
	} // by global direction (x,y,z)
	delete &f_1_3;
	delete &f_2_3;
}
//*************************************
//*/


// new one

//*************************************
void THydrodynamicsProblem::SolveTransferEquation__(TRnRelease3dSpace& lVecCur, TRnRelease3dSpace& lVecPrev,
	TRnRelease3dSpace&  lU, TRnRelease3dSpace& lV, TRnRelease3dSpace& lW, TRnRelease3dSpace& lRP,
	const TNodesGrid& lGrid, TRnRelease3dSpace& lDiffusionCoefficient, int lSchemeType)
{	
	//Preparing variables (Begin)
	TRnRelease3dSpace& f = lVecCur;
	TRnRelease3dSpace& f_1 = lVecPrev;
	TRnRelease3dSpace& u = lU;
	TRnRelease3dSpace& v = lV;
	TRnRelease3dSpace& w = lW;
	TRnRelease3dSpace& rp = lRP;

	TRnRelease3dSpace& f_1_3 = dynamic_cast<TRnRelease3dSpace&>(f_1.Copy());
	TRnRelease3dSpace& f_2_3 = dynamic_cast<TRnRelease3dSpace&>(f_1.Copy());

	const T3dNormalGrid& grid = dynamic_cast<const T3dNormalGrid&>(lGrid);

	const T3dNumberMask& mask = grid.Mask;

	const TSeparator& sepX = grid.Separator1;
	const TSeparator& sepY = grid.Separator2;
	const TSeparator& sepZ = grid.Separator3;

	int N = sepX.EndIndex;
	int L = sepY.EndIndex;
	int M = sepZ.EndIndex;

	const double* hx = sepX.GetSeparation();
	const double* hy = sepY.GetSeparation();
	const double* hz = sepZ.GetSeparation();

	double tau = fTimeSeparator->SeparationValue;

	TRnRelease3dSpace& diffCoefVar = lDiffusionCoefficient;

	//Preparing variables (End)

	//Direction order (Begin)
	int directionOrder = 1;
	int directionX;
	int directionY;
	int directionZ;
	if (directionOrder == 1) { directionX = 1; directionY = 2; directionZ = 3; }
	else if (directionOrder == 2) { directionX = 1; directionY = 3; directionZ = 2; }
	else if (directionOrder == 3) { directionX = 2; directionY = 1; directionZ = 3; }
	else if (directionOrder == 4) { directionX = 2; directionY = 3; directionZ = 1; }
	else if (directionOrder == 5) { directionX = 3; directionY = 1; directionZ = 2; }
	else if (directionOrder == 6) { directionX = 3; directionY = 2; directionZ = 1; }
	else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction order."; }
	//Direction order (End)

	//Terms using (Begin)
	int flagConvectiveX = 1; int flagDiffusionX = 1;
	int flagConvectiveY = 1; int flagDiffusionY = 1;
	int flagConvectiveZ = 1; int flagDiffusionZ = 1;
	//Terms using (Begin)

	//Scheme type (Begin)
	//0 - central;
	//1 - upstream;
	int schemeType = lSchemeType;
	//Scheme type (End)

	for (int directionGlobal = 1; directionGlobal<4; directionGlobal++) //splitting scheame steps
	{
		for (int directionFake = 1; directionFake<4; directionFake++) //for border condtions recalculating
		{
			int direction;
			if (directionFake == 1) { direction = directionGlobal; }
			else
			{
				if (directionGlobal == 1) { direction = directionFake; }
				else if (directionGlobal == 2)
				{
					if (directionFake == 2) { direction = 1; }
					else if (directionFake == 3) { direction = 3; }
					else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown fake direction."; }
				}
				else if (directionGlobal == 3) { direction = directionFake - 1; }
				else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown global direction."; }
			}


			int dir1EndIndex;
			int dir2EndIndex;
			int dirActEndIndex;

			if (direction == directionX)
			{
				dir1EndIndex = M;
				dir2EndIndex = L;
				dirActEndIndex = N;
			}
			else if (direction == directionY)
			{
				dir1EndIndex = M;
				dir2EndIndex = N;
				dirActEndIndex = L;
			}
			else if (direction == directionZ)
			{
				dir1EndIndex = L;
				dir2EndIndex = N;
				dirActEndIndex = M;
			}
			else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }

			int flagConvective = 0;
			if (
				((direction == directionX) && (flagConvectiveX == 1)) ||
				((direction == directionY) && (flagConvectiveY == 1)) ||
				((direction == directionZ) && (flagConvectiveZ == 1))
				) {
				flagConvective = 1;
			}

			int flagDiffusion = 0;
			if (
				((direction == directionX) && (flagDiffusionX == 1)) ||
				((direction == directionY) && (flagDiffusionY == 1)) ||
				((direction == directionZ) && (flagDiffusionZ == 1))
				) {
				flagDiffusion = 1;
			}


			for (int dir1 = 0; dir1<dir1EndIndex + 1; dir1++)
			{
				for (int dir2 = 0; dir2<dir2EndIndex + 1; dir2++)
				{
					int dirActSecBeg = -1;
					int maskActSecBeg;

					int dirActSecEnd;
					int maskActSecEnd;

					for (int dirActSec = 0; dirActSec<dirActEndIndex + 1; dirActSec++)
					{
						int maskActSec;
						if (direction == directionX) { maskActSec = mask[dirActSec][dir2][dir1]; }
						else if (direction == directionY) { maskActSec = mask[dir2][dirActSec][dir1]; }
						else if (direction == directionZ) { maskActSec = mask[dir2][dir1][dirActSec]; }
						else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }

						if (dirActSecBeg<0)
						{
							if (maskActSec == TActualPoint)
							{
								if (dirActSec == 0)
								{
									throw "Error in THydrodynamicsProblem::SolveTransferEquation: Initial border not found.";
								}
								else
								{
									dirActSecBeg = dirActSec - 1;
									if (direction == directionX) { maskActSecBeg = mask[dirActSecBeg][dir2][dir1]; }
									else if (direction == directionY) { maskActSecBeg = mask[dir2][dirActSecBeg][dir1]; }
									else if (direction == directionZ) { maskActSecBeg = mask[dir2][dir1][dirActSecBeg]; }
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								}
							}
						}
						else if (maskActSec != TActualPoint)
						{
							dirActSecEnd = dirActSec;
							maskActSecEnd = maskActSec;

							//Сhecking sector (Begin)
							if ((maskActSecBeg == TActualPoint) || (maskActSecBeg == TFictivePoint) || (maskActSecEnd == TActualPoint) || (maskActSecEnd == TFictivePoint))
							{
								throw "Error in THydrodynamicsProblem::SolveTransferEquation: Incorrect sector border.";
							}

							if ((dirActSecBeg<0) || (dirActSecEnd - dirActSecBeg)<2)
							{
								throw "Error in THydrodynamicsProblem::SolveTransferEquation: Incorrect sector size.";
							}
							//Сhecking sector (End)



							//Only borders recalculation and continue (Begin)
							if (direction != directionGlobal)
							{
								TRnRelease3dSpace* f_AssignPtr;
								if (directionGlobal == 1) { f_AssignPtr = &f_1_3; }
								else if (directionGlobal == 2) { f_AssignPtr = &f_2_3; }
								else if (directionGlobal == 3) { f_AssignPtr = &f; }
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown global direction."; }
								TRnRelease3dSpace& f_Assign = *f_AssignPtr;

								//Sector begin (Begin)
								if ((maskActSecBeg == TDefinedBorderPoint) || (maskActSecBeg == TEquationBorderPoint))
								{
									if (direction == directionX) { f_Assign[dirActSecBeg][dir2][dir1] = rp[dirActSecBeg][dir2][dir1]; }
									else if (direction == directionY) { f_Assign[dir2][dirActSecBeg][dir1] = rp[dir2][dirActSecBeg][dir1]; }
									else if (direction == directionZ) { f_Assign[dir2][dir1][dirActSecBeg] = rp[dir2][dir1][dirActSecBeg]; }
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								}
								else if (maskActSecBeg == TPreDefinedBorderPoint)
								{
									if (direction == directionX) { f_Assign[dirActSecBeg][dir2][dir1] = 2 * rp[dirActSecBeg][dir2][dir1] - f_Assign[dirActSecBeg + 1][dir2][dir1]; }
									else if (direction == directionY) { f_Assign[dir2][dirActSecBeg][dir1] = 2 * rp[dir2][dirActSecBeg][dir1] - f_Assign[dir2][dirActSecBeg + 1][dir1]; }
									else if (direction == directionZ) { f_Assign[dir2][dir1][dirActSecBeg] = 2 * rp[dir2][dir1][dirActSecBeg] - f_Assign[dir2][dir1][dirActSecBeg + 1]; }
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								}
								else if ((maskActSecBeg == TNormalBorderPoint) || (maskActSecBeg == TPreNormalBorderPoint))
								{
									//! need to check acceptable normal direction

									if (direction == directionX) { f_Assign[dirActSecBeg][dir2][dir1] = -(hx[dirActSecBeg])*rp[dirActSecBeg][dir2][dir1] + f_Assign[dirActSecBeg + 1][dir2][dir1]; }
									else if (direction == directionY) { f_Assign[dir2][dirActSecBeg][dir1] = -(hy[dirActSecBeg])*rp[dir2][dirActSecBeg][dir1] + f_Assign[dir2][dirActSecBeg + 1][dir1]; }
									else if (direction == directionZ) { f_Assign[dir2][dir1][dirActSecBeg] = -(hz[dirActSecBeg])*rp[dir2][dir1][dirActSecBeg] + f_Assign[dir2][dir1][dirActSecBeg + 1]; }
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								}
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of begin sector."; }
								//Sector begin (End)


								//Sector end (Begin)
								if ((maskActSecEnd == TDefinedBorderPoint) || (maskActSecEnd == TEquationBorderPoint))
								{
									if (direction == directionX) { f_Assign[dirActSecEnd][dir2][dir1] = rp[dirActSecEnd][dir2][dir1]; }
									else if (direction == directionY) { f_Assign[dir2][dirActSecEnd][dir1] = rp[dir2][dirActSecEnd][dir1]; }
									else if (direction == directionZ) { f_Assign[dir2][dir1][dirActSecEnd] = rp[dir2][dir1][dirActSecEnd]; }
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								}
								else if (maskActSecEnd == TPreDefinedBorderPoint)
								{
									if (direction == directionX) { f_Assign[dirActSecEnd][dir2][dir1] = 2 * rp[dirActSecEnd][dir2][dir1] - f_Assign[dirActSecEnd - 1][dir2][dir1]; }
									else if (direction == directionY) { f_Assign[dir2][dirActSecEnd][dir1] = 2 * rp[dir2][dirActSecEnd][dir1] - f_Assign[dir2][dirActSecEnd - 1][dir1]; }
									else if (direction == directionZ) { f_Assign[dir2][dir1][dirActSecEnd] = 2 * rp[dir2][dir1][dirActSecEnd] - f_Assign[dir2][dir1][dirActSecEnd - 1]; }
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								}
								else if ((maskActSecEnd == TNormalBorderPoint) || (maskActSecEnd == TPreNormalBorderPoint))
								{
									//! need to check acceptable normal direction

									if (direction == directionX) { f_Assign[dirActSecEnd][dir2][dir1] = (hx[dirActSecEnd - 1])*rp[dirActSecEnd][dir2][dir1] + f_Assign[dirActSecEnd - 1][dir2][dir1]; }
									else if (direction == directionY) { f_Assign[dir2][dirActSecEnd][dir1] = (hy[dirActSecEnd - 1])*rp[dir2][dirActSecEnd][dir1] + f_Assign[dir2][dirActSecEnd - 1][dir1]; }
									else if (direction == directionZ) { f_Assign[dir2][dir1][dirActSecEnd] = (hz[dirActSecEnd - 1])*rp[dir2][dir1][dirActSecEnd] + f_Assign[dir2][dir1][dirActSecEnd - 1]; }
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								}
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of end sector."; }
								//Sector end (End)

								dirActSecBeg = -1; //sector has been calculated
								continue;
							} //if it is fake direction
							  //Only borders recalculation and continue (End)





							  //Need for fictive points (Begin)
							int fictiveLeft = 0;
							if ((maskActSecBeg != TDefinedBorderPoint) && (maskActSecBeg != TEquationBorderPoint)) { fictiveLeft = 1; }

							int fictiveRight = 0;
							if ((maskActSecEnd != TDefinedBorderPoint) && (maskActSecEnd != TEquationBorderPoint)) { fictiveRight = 1; }

							int dimActSec = dirActSecEnd - dirActSecBeg + 1 + fictiveLeft + fictiveRight;
							//Need for fictive points (End)


							//Preparing scalar objects (Begin)
							int scalDim = dimActSec - 1;
							double* scalVar = new double[scalDim + 1];
							double* scalA = new double[scalDim + 1];
							double* scalB = new double[scalDim + 1];
							double* scalC = new double[scalDim + 1];
							double* scalRP = new double[scalDim + 1];

							scalA[0] = 0;  scalA[scalDim] = 0;
							scalB[0] = 0;  scalB[scalDim] = 0;
							scalC[0] = 0;  scalC[scalDim] = 0;
							scalRP[0] = 0;  scalRP[scalDim] = 0;


							//Borders points (Begin)
							if ((maskActSecBeg == TDefinedBorderPoint) || (maskActSecBeg == TEquationBorderPoint))
							{
								if (direction == directionX) { scalVar[0] = rp[dirActSecBeg][dir2][dir1]; }
								else if (direction == directionY) { scalVar[0] = rp[dir2][dirActSecBeg][dir1]; }
								else if (direction == directionZ) { scalVar[0] = rp[dir2][dir1][dirActSecBeg]; }
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
							}
							else if (maskActSecBeg == TPreDefinedBorderPoint)
							{
								scalVar[0] = 0;
								scalA[1] = 1;
								scalB[1] = 0.5;
								scalC[1] = 0.5;
								if (direction == directionX) { scalRP[1] = rp[dirActSecBeg][dir2][dir1]; }
								else if (direction == directionY) { scalRP[1] = rp[dir2][dirActSecBeg][dir1]; }
								else if (direction == directionZ) { scalRP[1] = rp[dir2][dir1][dirActSecBeg]; }
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
							}
							else if ((maskActSecBeg == TNormalBorderPoint) || (maskActSecBeg == TPreNormalBorderPoint))
							{
								//! need to check acceptable normal direction
								scalVar[0] = 0;
								scalA[1] = 1;
								if (direction == directionX)
								{
									scalB[1] = -1 / hx[dirActSecBeg];
									scalC[1] = 1 / hx[dirActSecBeg];
									scalRP[1] = rp[dirActSecBeg][dir2][dir1];
								}
								else if (direction == directionY)
								{
									scalB[1] = -1 / hy[dirActSecBeg];
									scalC[1] = 1 / hy[dirActSecBeg];
									scalRP[1] = rp[dir2][dirActSecBeg][dir1];
								}
								else if (direction == directionZ)
								{
									scalB[1] = -1 / hz[dirActSecBeg];
									scalC[1] = 1 / hz[dirActSecBeg];
									scalRP[1] = rp[dir2][dir1][dirActSecBeg];
								}
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
							}
							else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of begin sector."; }



							if ((maskActSecEnd == TDefinedBorderPoint) || (maskActSecEnd == TEquationBorderPoint))
							{
								if (direction == directionX) { scalVar[scalDim] = rp[dirActSecEnd][dir2][dir1]; }
								else if (direction == directionY) { scalVar[scalDim] = rp[dir2][dirActSecEnd][dir1]; }
								else if (direction == directionZ) { scalVar[scalDim] = rp[dir2][dir1][dirActSecEnd]; }
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
							}
							else if (maskActSecEnd == TPreDefinedBorderPoint)
							{
								scalVar[scalDim] = 0;
								scalC[scalDim - 1] = 1;
								scalA[scalDim - 1] = 0.5;
								scalB[scalDim - 1] = 0.5;
								if (direction == directionX) { scalRP[scalDim - 1] = rp[dirActSecEnd][dir2][dir1]; }
								else if (direction == directionY) { scalRP[scalDim - 1] = rp[dir2][dirActSecEnd][dir1]; }
								else if (direction == directionZ) { scalRP[scalDim - 1] = rp[dir2][dir1][dirActSecEnd]; }
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
							}
							else if ((maskActSecEnd == TNormalBorderPoint) || (maskActSecEnd == TPreNormalBorderPoint))
							{
								scalVar[scalDim] = 0;
								scalC[scalDim - 1] = 1;
								if (direction == directionX)
								{
									scalA[scalDim - 1] = -1 / hx[dirActSecEnd - 1];
									scalB[scalDim - 1] = 1 / hx[dirActSecEnd - 1];
									scalRP[scalDim - 1] = rp[dirActSecEnd][dir2][dir1];
								}
								else if (direction == directionY)
								{
									scalA[scalDim - 1] = -1 / hy[dirActSecEnd - 1];
									scalB[scalDim - 1] = 1 / hy[dirActSecEnd - 1];
									scalRP[scalDim - 1] = rp[dir2][dirActSecEnd][dir1];
								}
								else if (direction == directionZ)
								{
									scalA[scalDim - 1] = -1 / hz[dirActSecEnd - 1];
									scalB[scalDim - 1] = 1 / hz[dirActSecEnd - 1];
									scalRP[scalDim - 1] = rp[dir2][dir1][dirActSecEnd];
								}
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
							}
							else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of end sector."; }
							//Borders points (End)


							//Actual points (Begin)
							for (int scalLoop = 1 + fictiveLeft; scalLoop<scalDim - fictiveRight; scalLoop++)
							{
								int dirAct = dirActSecBeg + scalLoop - fictiveLeft;

								int i; int j; int k;
								TRnRelease3dSpace* vecTrPtr;
								if (direction == directionX) { i = dirAct; j = dir2; k = dir1; vecTrPtr = &u; }
								else if (direction == directionY) { i = dir2; j = dirAct; k = dir1; vecTrPtr = &v; }
								else if (direction == directionZ) { i = dir2; j = dir1; k = dirAct; vecTrPtr = &w; }
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								TRnRelease3dSpace& vecTr = *vecTrPtr;

								if (mask[i][j][k] != TActualPoint) { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Incorrect sector innerr point."; }

								double h_1X; double hX; double hx2X; double hxhxhX;
								h_1X = hx[i - 1];
								hX = hx[i];
								hx2X = hx[i - 1] + hx[i];
								hxhxhX = hx[i - 1] * hx[i] * (hx[i - 1] + hx[i]) / 2;

								double h_1Y; double hY; double hx2Y; double hxhxhY;
								h_1Y = hy[j - 1];
								hY = hy[j];
								hx2Y = hy[j - 1] + hy[j];
								hxhxhY = hy[j - 1] * hy[j] * (hy[j - 1] + hy[j]) / 2;

								double h_1Z; double hZ; double hx2Z; double hxhxhZ;
								h_1Z = hz[k - 1];
								hZ = hz[k];
								hx2Z = hz[k - 1] + hz[k];
								hxhxhZ = hz[k - 1] * hz[k] * (hz[k - 1] + hz[k]) / 2;

								double h_1; double h; double hx2; double hxhxh;
								if (direction == directionX) { h_1 = h_1X; h = hX; hx2 = hx2X; hxhxh = hxhxhX; }
								else if (direction == directionY) { h_1 = h_1Y; h = hY; hx2 = hx2Y; hxhxh = hxhxhY; }
								else if (direction == directionZ) { h_1 = h_1Z; h = hZ; hx2 = hx2Z; hxhxh = hxhxhZ; }
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }


								double diffCoef = diffCoefVar[i][j][k];

								if (schemeType == 0)
								{
									scalA[scalLoop] = flagConvective*(0.5*(-vecTr[i][j][k] / hx2)) + flagDiffusion*(0.5*(-(diffCoef)*h / hxhxh));
									scalB[scalLoop] = 1 / tau + flagDiffusion*(0.5*(diffCoef)*hx2 / hxhxh);
									scalC[scalLoop] = flagConvective*(0.5*(vecTr[i][j][k] / hx2)) + flagDiffusion*(0.5*(-(diffCoef)*h_1 / hxhxh));
								}
								else if (schemeType == 1)
								{
									double vecFlow = vecTr[i][j][k];
									if (vecFlow<0)
									{
										scalA[scalLoop] = flagDiffusion*(0.5*(-(diffCoef)*h / hxhxh));
										scalB[scalLoop] = 1 / tau + flagConvective*(0.5*(-vecTr[i][j][k] / h)) + flagDiffusion*(0.5*(diffCoef)*hx2 / hxhxh);
										scalC[scalLoop] = flagConvective*(0.5*(vecTr[i][j][k] / h)) + flagDiffusion*(0.5*(-(diffCoef)*h_1 / hxhxh));
									}
									else
									{
										scalA[scalLoop] = flagConvective*(0.5*(-vecTr[i][j][k] / h)) + flagDiffusion*(0.5*(-(diffCoef)*h / hxhxh));
										scalB[scalLoop] = 1 / tau + flagConvective*(0.5*(vecTr[i][j][k] / h)) + flagDiffusion*(0.5*(diffCoef)*hx2 / hxhxh);
										scalC[scalLoop] = flagDiffusion*(0.5*(-(diffCoef)*h_1 / hxhxh));
									}
								}
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

								double F = (1 / tau)*f_1[i][j][k] + rp[i][j][k];

								//Help variables for simple writing F (Begin)
								double FC;
								double FD;
								TRnRelease3dSpace* fDir = nullptr;
								//Help variables for simple writing F (End)

								if (direction == directionX)
								{
									if (schemeType == 0) { F = F + flagConvectiveX*(0.5*(-u[i][j][k])*(f_1[i + 1][j][k] - f_1[i - 1][j][k]) / hx2X); }
									else if (schemeType == 1)
									{
										double vecFlow = u[i][j][k];
										if (vecFlow<0) { F = F + flagConvectiveX*(0.5*(-u[i][j][k])*(f_1[i + 1][j][k] - f_1[i][j][k]) / hX); }
										else { F = F + flagConvectiveX*(0.5*(-u[i][j][k])*(f_1[i][j][k] - f_1[i - 1][j][k]) / hX); }
									}
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

									F = F + flagDiffusionX*(0.5*(diffCoef)*((hX*f_1[i - 1][j][k] - hx2X*f_1[i][j][k] + h_1X*f_1[i + 1][j][k]) / hxhxhX));


									//Y (Begin)
									if (directionY>direction) { fDir = &f_1; }
									else if (directionY == 1) { fDir = &f_1_3; }
									else if (directionY == 2) { fDir = &f_2_3; }

									FC = 0; FD = 0;
									for (int item = 1; item <= 2; item++)
									{
										TRnRelease3dSpace* f_Ptr = fDir; if (item == 2) { f_Ptr = &f_1; }
										TRnRelease3dSpace& f_ = *f_Ptr;
										if (schemeType == 0) { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j + 1][k] - f_[i][j - 1][k]) / hx2Y; }
										else if (schemeType == 1)
										{
											double vecFlow = v[i][j][k];
											if (vecFlow<0) { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j + 1][k] - f_[i][j][k]) / hY; }
											else { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j][k] - f_[i][j - 1][k]) / hY; }
										}
										else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

										FD = FD + 0.5*(diffCoef)*((hY*f_[i][j - 1][k] - hx2Y*f_[i][j][k] + h_1Y*f_[i][j + 1][k]) / hxhxhY);
									}
									F = F + flagConvectiveY*FC + flagDiffusionY*FD;
									//Y (End)

									//Z (Begin)
									if (directionZ>direction) { fDir = &f_1; }
									else if (directionZ == 1) { fDir = &f_1_3; }
									else if (directionZ == 2) { fDir = &f_2_3; }

									FC = 0; FD = 0;
									for (int item = 1; item <= 2; item++)
									{
										TRnRelease3dSpace* f_Ptr = fDir; if (item == 2) { f_Ptr = &f_1; }
										TRnRelease3dSpace& f_ = *f_Ptr;
										if (schemeType == 0) { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k + 1] - f_[i][j][k - 1]) / hx2Z; }
										else if (schemeType == 1)
										{
											double vecFlow = w[i][j][k];
											if (vecFlow<0) { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k + 1] - f_[i][j][k]) / hZ; }
											else { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k] - f_[i][j][k - 1]) / hZ; }
										}
										else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

										FD = FD + 0.5*(diffCoef)*((hZ*f_[i][j][k - 1] - hx2Z*f_[i][j][k] + h_1Z*f_[i][j][k + 1]) / hxhxhZ);
									}
									F = F + flagConvectiveZ*FC + flagDiffusionZ*FD;
									//Z (End)
								}
								else if (direction == directionY)
								{
									if (schemeType == 0) { F = F + flagConvectiveY*(0.5*(-v[i][j][k])*(f_1[i][j + 1][k] - f_1[i][j - 1][k]) / hx2Y); }
									else if (schemeType == 1)
									{
										double vecFlow = v[i][j][k];
										if (vecFlow<0) { F = F + flagConvectiveY*(0.5*(-v[i][j][k])*(f_1[i][j + 1][k] - f_1[i][j][k]) / hY); }
										else { F = F + flagConvectiveY*(0.5*(-v[i][j][k])*(f_1[i][j][k] - f_1[i][j - 1][k]) / hY); }
									}
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }


									F = F + flagDiffusionY*(0.5*(diffCoef)*((hY*f_1[i][j - 1][k] - hx2Y*f_1[i][j][k] + h_1Y*f_1[i][j + 1][k]) / hxhxhY));

									//X (Begin)
									if (directionX>direction) { fDir = &f_1; }
									else if (directionX == 1) { fDir = &f_1_3; }
									else if (directionX == 2) { fDir = &f_2_3; }

									FC = 0; FD = 0;
									for (int item = 1; item <= 2; item++)
									{
										TRnRelease3dSpace* f_Ptr = fDir; if (item == 2) { f_Ptr = &f_1; }
										TRnRelease3dSpace& f_ = *f_Ptr;
										if (schemeType == 0) { FC = FC + 0.5*(-u[i][j][k])*(f_[i + 1][j][k] - f_[i - 1][j][k]) / hx2X; }
										else if (schemeType == 1)
										{
											double vecFlow = u[i][j][k];
											if (vecFlow<0) { FC = FC + 0.5*(-u[i][j][k])*(f_[i + 1][j][k] - f_[i][j][k]) / hX; }
											else { FC = FC + 0.5*(-u[i][j][k])*(f_[i][j][k] - f_[i - 1][j][k]) / hX; }
										}
										else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

										FD = FD + 0.5*(diffCoef)*((hX*f_[i - 1][j][k] - hx2X*f_[i][j][k] + h_1X*f_[i + 1][j][k]) / hxhxhX);
									}
									F = F + flagConvectiveX*FC + flagDiffusionX*FD;
									//X (End)


									//Z (Begin)
									if (directionZ>direction) { fDir = &f_1; }
									else if (directionZ == 1) { fDir = &f_1_3; }
									else if (directionZ == 2) { fDir = &f_2_3; }

									FC = 0; FD = 0;
									for (int item = 1; item <= 2; item++)
									{
										TRnRelease3dSpace* f_Ptr = fDir; if (item == 2) { f_Ptr = &f_1; }
										TRnRelease3dSpace& f_ = *f_Ptr;
										if (schemeType == 0) { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k + 1] - f_[i][j][k - 1]) / hx2Z; }
										else if (schemeType == 1)
										{
											double vecFlow = w[i][j][k];
											if (vecFlow<0) { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k + 1] - f_[i][j][k]) / hZ; }
											else { FC = FC + 0.5*(-w[i][j][k])*(f_[i][j][k] - f_[i][j][k - 1]) / hZ; }
										}
										else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

										FD = FD + 0.5*(diffCoef)*((hZ*f_[i][j][k - 1] - hx2Z*f_[i][j][k] + h_1Z*f_[i][j][k + 1]) / hxhxhZ);
									}
									F = F + flagConvectiveZ*FC + flagDiffusionZ*FD;
									//Z (End)
								}
								else if (direction == directionZ)
								{
									if (schemeType == 0) { F = F + flagConvectiveZ*(0.5*(-w[i][j][k])*(f_1[i][j][k + 1] - f_1[i][j][k - 1]) / hx2Z); }
									else if (schemeType == 1)
									{
										double vecFlow = w[i][j][k];
										if (vecFlow<0) { F = F + flagConvectiveZ*(0.5*(-w[i][j][k])*(f_1[i][j][k + 1] - f_1[i][j][k]) / hZ); }
										else { F = F + flagConvectiveZ*(0.5*(-w[i][j][k])*(f_1[i][j][k] - f_1[i][j][k - 1]) / hZ); }
									}
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

									F = F + flagDiffusionZ*(0.5*(diffCoef)*((hZ*f_1[i][j][k - 1] - hx2Z*f_1[i][j][k] + h_1Z*f_1[i][j][k + 1]) / hxhxhZ));

									//X (Begin)
									if (directionX>direction) { fDir = &f_1; }
									else if (directionX == 1) { fDir = &f_1_3; }
									else if (directionX == 2) { fDir = &f_2_3; }

									FC = 0; FD = 0;
									for (int item = 1; item <= 2; item++)
									{
										TRnRelease3dSpace* f_Ptr = fDir; if (item == 2) { f_Ptr = &f_1; }
										TRnRelease3dSpace& f_ = *f_Ptr;
										if (schemeType == 0) { FC = FC + 0.5*(-u[i][j][k])*(f_[i + 1][j][k] - f_[i - 1][j][k]) / hx2X; }
										else if (schemeType == 1)
										{
											double vecFlow = u[i][j][k];
											if (vecFlow<0) { FC = FC + 0.5*(-u[i][j][k])*(f_[i + 1][j][k] - f_[i][j][k]) / hX; }
											else { FC = FC + 0.5*(-u[i][j][k])*(f_[i][j][k] - f_[i - 1][j][k]) / hX; }
										}
										else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

										FD = FD + 0.5*(diffCoef)*((hX*f_[i - 1][j][k] - hx2X*f_[i][j][k] + h_1X*f_[i + 1][j][k]) / hxhxhX);
									}
									F = F + flagConvectiveX*FC + flagDiffusionX*FD;
									//X (End)

									//Y (Begin)
									if (directionY>direction) { fDir = &f_1; }
									else if (directionY == 1) { fDir = &f_1_3; }
									else if (directionY == 2) { fDir = &f_2_3; }

									FC = 0; FD = 0;
									for (int item = 1; item <= 2; item++)
									{
										TRnRelease3dSpace* f_Ptr = fDir; if (item == 2) { f_Ptr = &f_1; }
										TRnRelease3dSpace& f_ = *f_Ptr;
										if (schemeType == 0) { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j + 1][k] - f_[i][j - 1][k]) / hx2Y; }
										else if (schemeType == 1)
										{
											double vecFlow = v[i][j][k];
											if (vecFlow<0) { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j + 1][k] - f_[i][j][k]) / hY; }
											else { FC = FC + 0.5*(-v[i][j][k])*(f_[i][j][k] - f_[i][j - 1][k]) / hY; }
										}
										else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

										FD = FD + 0.5*(diffCoef)*((hY*f_[i][j - 1][k] - hx2Y*f_[i][j][k] + h_1Y*f_[i][j + 1][k]) / hxhxhY);
									}
									F = F + flagConvectiveY*FC + flagDiffusionY*FD;
									//Y (End)
								}
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								scalRP[scalLoop] = F;
							} // by act direction
							  //Actual points (End)

							  //Preparing scalar objects (End)


							//Scalar running (Begin)
							TSLAE::SolveByScalarRunning(scalDim, scalVar, scalA, scalB, scalC, scalRP);
							//Scalar running (End)


							//Assigninig (Begin)
							TRnRelease3dSpace* f_AssignPtr;
							if (direction == 1) { f_AssignPtr = &f_1_3; }
							else if (direction == 2) { f_AssignPtr = &f_2_3; }
							else if (direction == 3) { f_AssignPtr = &f; }
							else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
							TRnRelease3dSpace& f_Assign = *f_AssignPtr;

							for (int scalLoop = 0 + fictiveLeft; scalLoop<scalDim + 1 - fictiveRight; scalLoop++)
							{
								int dirAct = dirActSecBeg + scalLoop - fictiveLeft;
								int i; int j; int k;
								if (direction == directionX) { i = dirAct; j = dir2; k = dir1; }
								else if (direction == directionY) { i = dir2; j = dirAct; k = dir1; }
								else if (direction == directionZ) { i = dir2; j = dir1; k = dirAct; }
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }

								f_Assign[i][j][k] = scalVar[scalLoop];
							}
							//Assigninig (End)

							//Deleting scalar objects (Begin)
							delete[] scalVar;
							delete[] scalA;
							delete[] scalB;
							delete[] scalC;
							delete[] scalRP;
							//Deleting scalar objects (End)

							dirActSecBeg = -1; //sector has been calculated
						} // sector calculating
					} // by sectors of act direction
				} // by direction 2
			} // by direction 1
		} // by fake direction (x,y,z) - for border conditions recalculating
	} // by global direction (x,y,z)
	delete &f_1_3;
	delete &f_2_3;
}
//*************************************
//#################################################################


//*************************************
TNodesGrid* T3dTranferProblem::BuildNormalGrid(TSeparator& lSep1, TSeparator& lSep2, TSeparator& lSep3, const TMask& lMask)
{	
   int N = lSep1.EndIndex;
   int L = lSep2.EndIndex;
   int M = lSep3.EndIndex;

   T3dVector**** normals = new T3dVector***[N+1];

   for (int i=0;i<N+1;i++)
   {
     normals[i] = new T3dVector**[L+1];
     for (int j=0;j<L+1;j++)
     {
       normals[i][j] = new T3dVector*[M+1];
       for (int k=0;k<M+1;k++)
       {
         T3dVector* normal = NULL;
		 if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) ) {normal = NULL;}
         else
         {
            if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==0) )
			{
				normal = new T3dVector;
				normal->X = 0;
				normal->Y = 0;
				normal->Z = -1;
			}
            else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==M) )
            {
				normal = new T3dVector;
				normal->X = 0;
				normal->Y = 0;
				normal->Z = 1;
			}

            else if ( (i==0) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
            {
				normal = new T3dVector;
				normal->X = -1;
				normal->Y = 0;
				normal->Z = 0;
			}
            else if ( (i==N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
			{
				normal = new T3dVector;
				normal->X = 1;
				normal->Y = 0;
				normal->Z = 0;
			}

			else if ( (i!=0) && (i!=N) && (j==0) && (k!=0) && (k!=M) )
			{
				normal = new T3dVector;
				normal->X = 0;
				normal->Y = -1;
				normal->Z = 0;
			}
			else if ( (i!=0) && (i!=N) && (j==L) && (k!=0) && (k!=M) )
			{
				normal = new T3dVector;
				normal->X = 0;
				normal->Y = 1;
				normal->Z = 0;
			}
         }
		 normals[i][j][k] = normal;
        }
       }
      }

      TNodesGrid* lGrid = new T3dNormalGrid(&lSep1,&lSep2,&lSep3,lMask,normals);

	  for (int i=0;i<N+1;i++)
      {
       for (int j=0;j<L+1;j++)
	   {
        for (int k=0;k<M+1;k++)
		{
			if ( (normals[i][j][k])!=NULL ) {delete normals[i][j][k];}
		}
		delete[] normals[i][j];
	   }
	   delete[] normals[i];
	  }
	  delete[] normals;// normals = NULL;

	  return lGrid;

}
//*************************************
TNodesGrid* T3dTranferProblem::BuildGrid()
{
   double h = fStep;
   double b = 1;
   TSeparator& sepU1 = *(new TUniformSeparator(0,b,h));
   TSeparator& sepU2 = *(new TUniformSeparator(0,b,h));
   TSeparator& sepU3 = *(new TUniformSeparator(0,b,h));

   TSeparator& sep1 = *(new TRandomSeparator(sepU1.CountNodes,sepU1.Dimension));
   delete &sepU1;
   TSeparator& sep2 = *(new TRandomSeparator(sepU2.CountNodes,sepU2.Dimension));
   delete &sepU2;
   TSeparator& sep3 = *(new TRandomSeparator(sepU3.CountNodes,sepU3.Dimension));
   delete &sepU3;

   int N = sep1.EndIndex;
   int L = sep2.EndIndex;
   int M = sep3.EndIndex;

   T3dNumberMask &mask = *( new T3dNumberMask(N+1,L+1,M+1));

   for (int i=0;i<N+1;i++)
   {
     for (int j=0;j<L+1;j++)
     {
       for (int k=0;k<M+1;k++)
       {
		 mask[i][j][k]=TFictivePoint;
		 if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
         {
             mask[i][j][k]=TActualPoint;
         }
         else
         {
            if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==0) )
			{
				mask[i][j][k]=TBorderPoint;
			}
            else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==M) )
            {
				mask[i][j][k]=TBorderPoint;
			}

            else if ( (i==0) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
            {
				mask[i][j][k]=TBorderPoint;
			}
            else if ( (i==N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TBorderPoint;
			}

			else if ( (i!=0) && (i!=N) && (j==0) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TBorderPoint;
			}
         }
       }
	 }
   }

   TNodesGrid* grid = BuildNormalGrid(sep1,sep2,sep3,mask);

   delete &mask;

   return grid;
}
//*************************************
TNodesGrid* T3dTranferProblem::BuildGridSolution()
{
   double h = fStep;
   double b = 1;
   TSeparator& sepU1 = *(new TUniformSeparator(0,b,h));
   TSeparator& sepU2 = *(new TUniformSeparator(0,b,h));
   TSeparator& sepU3 = *(new TUniformSeparator(0,b,h));

   TSeparator& sep1 = *(new TRandomSeparator(sepU1.CountNodes,sepU1.Dimension));
   TSeparator& sep2 = *(new TRandomSeparator(sepU2.CountNodes,sepU2.Dimension));
   TSeparator& sep3 = *(new TRandomSeparator(sepU3.CountNodes,sepU3.Dimension));

   delete &sepU1;
   delete &sepU2;
   delete &sepU3;

   int N = sep1.EndIndex;
   int L = sep2.EndIndex;
   int M = sep3.EndIndex;

   T3dNumberMask &mask = *( new T3dNumberMask(N+1,L+1,M+1));

   for (int i=0;i<N+1;i++)
   {
     for (int j=0;j<L+1;j++)
     {
       for (int k=0;k<M+1;k++)
       {
		 mask[i][j][k]=TFictivePoint;
		 if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
         {
             mask[i][j][k]=TActualPoint;
         }
         else
         {
            if ( (i==0) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
            {
				mask[i][j][k]=TDefinedBorderPoint;
			}
            else if ( (i==N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==0) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==0) )
			{
				mask[i][j][k]=TDefinedBorderPoint;
			}
            else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==M) )
            {
				mask[i][j][k]=TDefinedBorderPoint;
			}
         }
       }
     }
   }

   TNodesGrid* grid = BuildNormalGrid(sep1,sep2,sep3,mask);

   delete &mask;

   return grid;
}
//*************************************
void T3dTranferProblem::PrepareVariables()
{
  printf("Preparing variables: Starting....");
  printf("\n");

  T3dNormalGrid* gridPtr = dynamic_cast<T3dNormalGrid*>(fGridSolution);
  const T3dNormalGrid& grid = *gridPtr;
  int N = grid.GetSeparator1().EndIndex;
  int L = grid.GetSeparator2().EndIndex;
  int M = grid.GetSeparator3().EndIndex;
  const T3dNumberMask& mask = grid.Mask;

  const T3dVector* const* const* const* normals = grid.Normals;

  T3dNumberMask& spaceMask = mask.Copy();
  for (int i=0;i<N+1;i++)
  {
	for (int j=0;j<L+1;j++)
	{
	  for (int k=0;k<M+1;k++)
	  {
			   int mijk = spaceMask[i][j][k];
			   if ( (mijk==TActualPoint) || (mijk==TFictivePoint) ) {spaceMask[i][j][k] = mijk;}
			   else if (mijk==TDefinedBorderPoint)
			   {
				 spaceMask[i][j][k] = TFictivePoint;
			   }
			   else if (mijk==TEquationBorderPoint)
			   {
				 spaceMask[i][j][k] = TFictivePoint;
			   }
			   else if ( (mijk==TPreDefinedBorderPoint) || (mijk==TNormalBorderPoint) || (mijk==TPreNormalBorderPoint) || (mijk==TPreEquationBorderPoint) )
			   {				 spaceMask[i][j][k] = TActualPoint;
			   }
			   else
			   {
                 throw "Error in T3dTranferProblem::PrepareVariables: Unknown point mask type.";
			   }
	   }
	 }
    }

    TMask** spaceMasks = new TMask*[1];
    spaceMasks[0] = &spaceMask;

    TRnRelease3dSpaceModel& spaceModel = *(new TRnRelease3dSpaceModel(1,spaceMasks));
	int spaceModelName = 1;
    TRnSpace::AddModel(spaceModelName,&spaceModel,false);

    fSolution = new TRnRelease3dSpace(spaceModelName);


	delete &spaceMask;
	delete[] spaceMasks;
	delete &spaceModel;

    printf("Preparing variables: OK.");
    printf("\n");
}
//*************************************
T3dTranferProblem::~T3dTranferProblem()
{
	printf("Deleting transfer objects: Starting....");
    printf("\n");

    delete fGridSolution;

	delete fSolution;

	TRnSpace::RemoveAllModels();

	printf("Deleting transfer objects: OK.");
}
//*************************************
void T3dTranferProblem::Solve()
{	
  printf("Solver: Starting....");
  printf("\n");

  //Grid (Begin)
  T3dNormalGrid& grid = dynamic_cast<T3dNormalGrid&>(*fGridSolution);
  int N = grid.GetSeparator1().EndIndex;
  int L = grid.GetSeparator2().EndIndex;
  int M = grid.GetSeparator3().EndIndex;
  const T3dNumberMask& mask = grid.Mask;

  const double* x = grid.GetSeparator1().Dimension;
  const double* y = grid.GetSeparator2().Dimension;
  const double* z = grid.GetSeparator3().Dimension;
  //Grid (End)


  //Vectors defining (Begin)

  //Basic - current time step (Begin)
  TRnRelease3dSpace &f = *fSolution;

  //Initial conditions (Begin)
  for (int i=0;i<N+1;i++)
  {
   for (int j=0;j<L+1;j++)
   {
	for (int k=0;k<M+1;k++)
	{
		  int mijk = mask[i][j][k];
	      if (mijk==TFictivePoint)
	      {
		    f[i][j][k] = 0;
	      }
	      else if (mijk==TActualPoint)
	      {
		    f[i][j][k] = SolutionFunction(x[i],y[j],z[k],0);
	      }
	      else if ( (mijk==TDefinedBorderPoint) )
	      {
          	f[i][j][k] = SolutionFunction(x[i],y[j],z[k],0);
	      }
	      else if ( (mijk==TPreDefinedBorderPoint) || (mijk==TNormalBorderPoint) || (mijk==TPreNormalBorderPoint) || (mijk==TEquationBorderPoint) || (mijk==TPreEquationBorderPoint) )
	      {
            throw "Error in T3dTranferProblem::Solve: Method is not implemented for this point border type.";
	      }
	      else
	      {
            throw "Error in T3dTranferProblem::Solve: Unknown point mask type.";
	      }
	}
   }
  }
  //Initial conditions (End)

  //Basic - current time step (End)

  //Previous time step (Begin)
  TRnRelease3dSpace &f_1 = dynamic_cast<TRnRelease3dSpace &>(f.Copy());
  //Previous time step (End)

  //Transfer terms (Begin)
  TRnRelease3dSpace& uTr = dynamic_cast<TRnRelease3dSpace &>(f.CreateInstance());
  TRnRelease3dSpace& vTr = dynamic_cast<TRnRelease3dSpace &>(f.CreateInstance());
  TRnRelease3dSpace& wTr = dynamic_cast<TRnRelease3dSpace &>(f.CreateInstance());
  TRnRelease3dSpace& rpTr = dynamic_cast<TRnRelease3dSpace &>(f.CreateInstance());
  //Transfer terms (End)

  //Explicit solution (Begin)
  TRnRelease3dSpace& solution = dynamic_cast<TRnRelease3dSpace &>(f.CreateInstance());
  //Explicit solution, error (End)

  //Vectors defining (End)


  //Error (Begin)
  printf("Time step: ");
  printf("%d",0);
  printf(", ");
  {
  double errorNorma = 0;
  for (int i=0;i<N+1;i++)
  {
	 for (int j=0;j<L+1;j++)
	 {
      for (int k=0;k<M+1;k++)
      {
         double er = fabs(solution[i][j][k] - f[i][j][k]);
		 if (er>errorNorma) {errorNorma = er;}
	  }
	 }
  }
  printf("Error: ");
  printf("%.20f",errorNorma);
  printf("\n");
  }
  //Error (End)


  const double* t = fTimeSeparator->Dimension;
  double tau = fTimeSeparator->SeparationValue;
  int timeStepsCount = fTimeSeparator->EndIndex;


  for (int timeStepNumber=1;timeStepNumber<timeStepsCount+1;timeStepNumber++)
  {
    printf("Time step: ");
	printf("%d",timeStepNumber);
	printf(", ");

	//Transfer equation (Begin)

	//Convective items, right part, border conditions (Begin)
	int p = timeStepNumber;
    for (int i=0;i<N+1;i++)
    {
	 for (int j=0;j<L+1;j++)
	 {
      for (int k=0;k<M+1;k++)
      {
	   int mijk = mask[i][j][k];
	   if (mijk==TFictivePoint)
	   {
           uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0; rpTr[i][j][k] = 0;
		   f[i][j][k] = 0;
           solution[i][j][k] = 0;
	   }
	   else if (mijk==TActualPoint)
	   {
		   uTr[i][j][k]  = U(x[i],y[j],z[k],t[p]);
		   vTr[i][j][k]  = V(x[i],y[j],z[k],t[p]);
		   wTr[i][j][k]  = W(x[i],y[j],z[k],t[p]);
		   rpTr[i][j][k] = RightPartFunction(x[i],y[j],z[k],t[p]);

		   f[i][j][k] = 0;
           solution[i][j][k] = SolutionFunction(x[i],y[j],z[k],t[p]);
	   }
	   else if ( (mijk==TDefinedBorderPoint) )
	   {
          uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0;
		  rpTr[i][j][k] = SolutionFunction(x[i],y[j],z[k],t[p]);
		  //f[i][j][k] = SolutionFunction(x[i],y[j],z[k],t[p]);
		  f[i][j][k] = 0;
          solution[i][j][k] = SolutionFunction(x[i],y[j],z[k],t[p]);
	   }
	   else if ( (mijk==TPreDefinedBorderPoint) || (mijk==TNormalBorderPoint) || (mijk==TPreNormalBorderPoint) || (mijk==TEquationBorderPoint) || (mijk==TPreEquationBorderPoint) )
	   {
           throw "Error in T3dTranferProblem::Solve: Method is not implemented for this point border type.";
	   }
	   else
	   {
          throw "Error in T3dTranferProblem::Solve: Unknown point mask type.";
	   }
      } //by k
     } //by j
    } //by i
    //Convective items, right part, border conditions (End)

	SolveTransferEquation(f,f_1,uTr,vTr,wTr,rpTr,grid);
    //Transfer equation (End)

	//Error (Begin)
	double errorNorma = 0;
	for (int i=0;i<N+1;i++)
    {
	 for (int j=0;j<L+1;j++)
	 {
      for (int k=0;k<M+1;k++)
      {
         double er = fabs(solution[i][j][k] - f[i][j][k]);
		 if (er>errorNorma) {errorNorma = er;}
	  }
	 }
	}
	printf("Error: ");
	printf("%.20f",errorNorma);
	printf("\n");
	//Error (End)

	//Reassigning (Begin)
	f_1.Assign(f);
	//Reassigning (End)

  } // by time

  delete &f_1;

  delete &uTr;
  delete &vTr;
  delete &wTr;
  delete &rpTr;

  delete &solution;

  printf("Solver: OK.");
  printf("\n");
}
//*************************************



//#################################################################



//*************************************
TNodesGrid* TVPProblem::BuildLocalGrid(
		                                 double *lX, double *lY, double *lZ,
										 int lNx, int lNy, int lNz,
										 glTMaskWithNormals ***lMask
									  )
{
    TSeparator& sep1 = *(new TRandomSeparator(lNx,lX));
    TSeparator& sep2 = *(new TRandomSeparator(lNy,lY));
    TSeparator& sep3 = *(new TRandomSeparator(lNz,lZ));

    int N = sep1.EndIndex;
    int L = sep2.EndIndex;
    int M = sep3.EndIndex;

    T3dNumberMask &mask = *( new T3dNumberMask(N+1,L+1,M+1));

	T3dVector**** normals = new T3dVector***[N+1];

    for (int i=0;i<N+1;i++)
    {
      normals[i] = new T3dVector**[L+1];
      for (int j=0;j<L+1;j++)
      {
        normals[i][j] = new T3dVector*[M+1];
        for (int k=0;k<M+1;k++)
        {
          T3dVector* normal = NULL;
		  glTVector* n = (lMask[i][j][k]).normal;
		  if (n!=NULL)
		  {
		   normal = new T3dVector;
		   normal->X = n->x;
		   normal->Y = n->y;
		   normal->Z = n->z;
		  }
		  normals[i][j][k] = normal;
		  mask[i][j][k] = (lMask[i][j][k]).mask;
        }
      }
    }

    TNodesGrid* grid = new T3dNormalGrid(&sep1,&sep2,&sep3,mask,normals);

	for (int i=0;i<N+1;i++)
    {
      for (int j=0;j<L+1;j++)
	  {
        for (int k=0;k<M+1;k++)
		{
			if ( (normals[i][j][k])!=NULL ) {delete normals[i][j][k];}
		}
		delete[] normals[i][j];
	  }
	  delete[] normals[i];
	 }
	delete[] normals;// normals = NULL;

     delete &mask;

     return grid;
}
//*************************************
TNodesGrid* TVPProblem::BuildGrid()
{
	if (fStatement == NULL) {
        throw std::runtime_error("TVPProblem::BuildGrid: No statement.");
    }
	TMaskGenerator& state = *fStatement;
	return BuildLocalGrid(state.X,state.Y,state.Z,state.Nx,state.Ny,state.Nz,state.Mask);
}
//*************************************
TNodesGrid* TVPProblem::BuildGridU()
{
	if (fStatement == NULL) {throw "TVPProblem::BuildGridU: No statement.";}
	TMaskGenerator& state = *fStatement;
	return BuildLocalGrid(state.XU,state.YU,state.ZU,state.NxU,state.NyU,state.NzU,state.MaskU);
}
//*************************************
double TVPProblem::GetBorderConditionU(int lI, int lJ, int lK, int lN)
{
	if (fStatement == NULL) {throw "TVPProblem::GetBorderConditionU: No statement.";}
	TMaskGenerator& state = *fStatement;
    return state.getUCondition(lI,lJ,lK,lN);
}
//*************************************
TNodesGrid* TVPProblem::BuildGridV()
{
	if (fStatement == NULL) {throw "TVPProblem::BuildGridV: No statement.";}
	TMaskGenerator& state = *fStatement;
	return BuildLocalGrid(state.XV,state.YV,state.ZV,state.NxV,state.NyV,state.NzV,state.MaskV);
}
//*************************************
double TVPProblem::GetBorderConditionV(int lI, int lJ, int lK, int lN)
{
	if (fStatement == NULL) {throw "TVPProblem::GetBorderConditionV: No statement.";}
	TMaskGenerator& state = *fStatement;
	return state.getVCondition(lI,lJ,lK,lN);
}
//*************************************
TNodesGrid* TVPProblem::BuildGridW()
{
	if (fStatement == NULL) {throw "TVPProblem::BuildGridW: No statement.";}
	TMaskGenerator& state = *fStatement;
	return BuildLocalGrid(state.XW,state.YW,state.ZW,state.NxW,state.NyW,state.NzW,state.MaskW);
}
//*************************************
double TVPProblem::GetBorderConditionW(int lI, int lJ, int lK, int lN)
{
	if (fStatement == NULL) {throw "TVPProblem::GetBorderConditionW: No statement.";}
	TMaskGenerator& state = *fStatement;
	return state.getWCondition(lI,lJ,lK,lN);
}
//*************************************
TNodesGrid* TVPProblem::BuildGridP()
{
	if (fStatement == NULL) {throw "TVPProblem::BuildGridP: No statement.";}
	TMaskGenerator& state = *fStatement;
    return BuildLocalGrid(state.XP,state.YP,state.ZP,state.NxP,state.NyP,state.NzP,state.MaskP);
}
//*************************************
double TVPProblem::GetBorderConditionP(int lI, int lJ, int lK, int lN)
{
	if (fStatement == NULL) {throw "TVPProblem::GetBorderConditionP: No statement.";}
	TMaskGenerator& state = *fStatement;
	return state.getPCondition(lI,lJ,lK,lN);
}
//*************************************
void TVPProblem::LocalPrepareVariables()
{
  printf("Preparing variables: Starting....");
  printf("\n");

  int spaceModelsCount = 4;

  if (fSpaceModelNameU==5)
  {
   fSpaceModelNameU = fSpaceModelNameU - spaceModelsCount;
   fSpaceModelNameV = fSpaceModelNameV - spaceModelsCount;
   fSpaceModelNameW = fSpaceModelNameW - spaceModelsCount;
   fSpaceModelNameP = fSpaceModelNameP - spaceModelsCount;
  }
  else
  {
   fSpaceModelNameU = fSpaceModelNameU + spaceModelsCount;
   fSpaceModelNameV = fSpaceModelNameV + spaceModelsCount;
   fSpaceModelNameW = fSpaceModelNameW + spaceModelsCount;
   fSpaceModelNameP = fSpaceModelNameP + spaceModelsCount;
  }

  for (int spaceModelName=fSpaceModelNameU;spaceModelName<fSpaceModelNameU+spaceModelsCount;spaceModelName++)
  {
	T3dNormalGrid* gridPtr;
	if (spaceModelName==fSpaceModelNameU) {gridPtr = dynamic_cast<T3dNormalGrid*>(fGridU);}
	else if (spaceModelName==fSpaceModelNameV) {gridPtr = dynamic_cast<T3dNormalGrid*>(fGridV);}
	else if (spaceModelName==fSpaceModelNameW) {gridPtr = dynamic_cast<T3dNormalGrid*>(fGridW);}
	else if (spaceModelName==fSpaceModelNameP) {gridPtr = dynamic_cast<T3dNormalGrid*>(fGridP);}
	else {throw "Error in TVPProblem::PrepareVariables: Unknown space model name.";  }

    const T3dNormalGrid& grid = *gridPtr;
    int N = grid.GetSeparator1().EndIndex;
    int L = grid.GetSeparator2().EndIndex;
    int M = grid.GetSeparator3().EndIndex;
    const T3dNumberMask& mask = grid.Mask;

    const T3dVector* const* const* const* normals = grid.Normals;

    T3dNumberMask& spaceMask = mask.Copy();
    for (int i=0;i<N+1;i++)
    {
	  for (int j=0;j<L+1;j++)
	  {
		   for (int k=0;k<M+1;k++)
		   {
			   int mijk = spaceMask[i][j][k];
			   if ( (mijk==TActualPoint) || (mijk==TFictivePoint) ) {spaceMask[i][j][k] = mijk;}
			   else if (mijk==TDefinedBorderPoint)
			   {
				 spaceMask[i][j][k] = TFictivePoint;
			   }
			   else if (mijk==TEquationBorderPoint)
			   {
				 spaceMask[i][j][k] = TFictivePoint;
			   }
			   else if ( (mijk==TPreDefinedBorderPoint) || (mijk==TNormalBorderPoint) || (mijk==TPreNormalBorderPoint) || (mijk==TPreEquationBorderPoint) )
			   {
				 spaceMask[i][j][k] = TActualPoint;
			   }
			   else
			   {
                 throw "Error in TVPProblem::PrepareVariables: Unknown point mask type.";
			   }
		   }
	   }
    }

    TMask** spaceMasks = new TMask*[1];
    spaceMasks[0] = &spaceMask;

    TRnRelease3dSpaceModel& spaceModel = *(new TRnRelease3dSpaceModel(1,spaceMasks));
    TRnSpace::AddModel(spaceModelName,&spaceModel,false);

    if (spaceModelName==fSpaceModelNameU) {fU = new TRnRelease3dSpace(spaceModelName);}
	else if (spaceModelName==fSpaceModelNameV) {fV = new TRnRelease3dSpace(spaceModelName);}
	else if (spaceModelName==fSpaceModelNameW) {fW = new TRnRelease3dSpace(spaceModelName);}
	else if (spaceModelName==fSpaceModelNameP) {fP = new TRnRelease3dSpace(spaceModelName);}
	else {throw "Error in TVPProblem::PrepareVariables: Unknown space model name.";}

	delete &spaceMask;
	delete[] spaceMasks;
	delete &spaceModel;

	printf("Preparing variables [");
	printf("%d",spaceModelName);
	printf("]: OK.");
	printf("\n");
 }
 printf("Preparing variables: OK.");
  printf("\n");
}
//*************************************
void TVPProblem::PrepareVariables()
{
  fSpaceModelNameU = -3;
  fSpaceModelNameV = -2;
  fSpaceModelNameW = -1;
  fSpaceModelNameP =  0;

  LocalPrepareVariables();
}
//*************************************
TVPProblem::~TVPProblem()
{
	printf("\n\n\n");
	printf("Deleting VP objects: Starting....\n");



	delete fGridU;
	delete fGridV;
	delete fGridW;
	delete fGridP;

	delete fU;
	delete fV;
	delete fW;
	delete fP;

	TRnSpace::RemoveAllModels();

	if (fStatement!=NULL) {delete fStatement;}


	printf("Deleting VP objects: OK.\n");
}

double CountDistribution(double r)
{
	double const_to=0.3;

	if (r >= const_to)
		return 1.0;
	else
		return r / const_to;

}
//*************************************
void TVPProblem::SolveBySplitting()
{
	printf("\n\n\n");
	printf("Splitting method: Starting....\n");

	//Global parameters (Begin)
	bool turbulenceMode = TProblem::GlobalTurbulenceMode;
	bool usePackOperator = TProblem::GlobalUsePackOperator;
	//Global parameters (End)

	//Vectors and operators global pointers defining (Begin)

	//Current time step (Begin)
	//Already defined as class fields
	//Current time step (End)

	//Previous time step (Begin)
	TRnRelease3dSpace* u_1Ptr = NULL;
	TRnRelease3dSpace* v_1Ptr = NULL;
	TRnRelease3dSpace* w_1Ptr = NULL;
	//Previous time step (End)

	//Intermediate (Begin)
	TRnRelease3dSpace* u_Ptr = NULL;
	TRnRelease3dSpace* v_Ptr = NULL;
	TRnRelease3dSpace* w_Ptr = NULL;
	//Intermediate (End)

	//For control stating (Begin)
	TRnRelease3dSpace* uSPtr = NULL;
	TRnRelease3dSpace* vSPtr = NULL;
	TRnRelease3dSpace* wSPtr = NULL;
	//For control stating (End)

	//Poisson equation right part (Begin)
	TRnRelease3dSpace* rightPartPtr = NULL;
	//Poisson equation right part (End)

	//Concentration (Begin)
	TRnRelease3dSpace* C_1Ptr = NULL;
	TRnRelease3dSpace* CPtr = NULL;
	T3dNormalGrid* gridCPtr = NULL;
	double visValue = 1.5 / fRe;
	//Concentration (End)

	TRnRelease3dSpace* densityPtr = NULL;
	T3dNormalGrid* gridDensityPtr = NULL;
	double soluteDensityValue = 1.0;
	double envDensityValue = 1.0;

	double soluteConcentratonValue = 0.0;
	//double envDensityValue = 1.0;

	//Viscosity (Begin)
	T3dNormalGrid* gridVisPtr = NULL;
	TRnRelease3dSpace* visPtr = NULL;
	//Viscosity (End)
	
	T3dVariableDivergentedOperator* APtr = NULL;


	int spaceModelNameRn = 1054;
	T1dPackOperator* ARnPtr = NULL;
	TSLAE* slae = NULL;
	//Vectors and operators global pointers defining (End)




	//Initial conditions (Begin)
	(*fU) |= 0;
	(*fV) |= 0;
	(*fW) |= 0;
	(*fP) |= 0;
	//Initial conditions (End)


	double tau = fTimeSeparator->SeparationValue;
	int timeStepsCount = fTimeSeparator->EndIndex;


	bool flagRestruct = true;


	clock_t start = clock();
	clock_t end = clock();
	double diff = 0;


	for (int timeStepNumber = 1; timeStepNumber<timeStepsCount + 1; timeStepNumber++)
	{

		//ChangeStifness(timeStepNumber);
		printf("\n\n");
		printf("*****************************************************\n");
		printf("Time step = %d: Starting...\n", timeStepNumber);

		//Restructing (Begin)

			if (flagRestruct == true)
			{
				printf("\n");
				printf("####################################\n");
				printf("Restructing: Starting...\n");


				//Class fields (Begin)
				int spaceModelNameU_Old = fSpaceModelNameU;
				int spaceModelNameV_Old = fSpaceModelNameV;
				int spaceModelNameW_Old = fSpaceModelNameW;
				int spaceModelNameP_Old = fSpaceModelNameP;

				TRnRelease3dSpace* u_Old = fU;
				TRnRelease3dSpace* v_Old = fV;
				TRnRelease3dSpace* w_Old = fW;
				TRnRelease3dSpace* p_Old = fP;


				printf("ReBuilding VP objects: Starting....");
				printf("\n");
				ReBuildGrids();
				printf("ReBuilding VP objects: OK.\n");

				LocalPrepareVariables();


				//ReAssigning (Begin)
				for (int vp = 1; vp<5; vp++)
				{
					TNodesGrid* gridPtr;
					if (vp == 1) { gridPtr = fGridU; }
					else if (vp == 2) { gridPtr = fGridV; }
					else if (vp == 3) { gridPtr = fGridW; }
					else if (vp == 4) { gridPtr = fGridP; }
					else { throw "Error in TVPProblem::SolveBySplitting: Unknown VP-object number."; }

					T3dNormalGrid& grid = dynamic_cast<T3dNormalGrid&>(*gridPtr);
					int N = grid.GetSeparator1().EndIndex;
					int L = grid.GetSeparator2().EndIndex;
					int M = grid.GetSeparator3().EndIndex;
					const T3dNumberMask& mask = grid.Mask;

					TRnRelease3dSpace* vpPrevPtr;
					TRnRelease3dSpace* vpCurPtr;
					if (vp == 1) { vpPrevPtr = u_Old; vpCurPtr = fU; }
					else if (vp == 2) { vpPrevPtr = v_Old; vpCurPtr = fV; }
					else if (vp == 3) { vpPrevPtr = w_Old; vpCurPtr = fW; }
					else if (vp == 4) { vpPrevPtr = p_Old; vpCurPtr = fP; }
					else { throw "Error in TVPProblem::SolveBySplitting: Unknown VP-object number."; }
					TRnRelease3dSpace& vpPrev = *vpPrevPtr;
					TRnRelease3dSpace& vpCur = *vpCurPtr;

					for (int i = 0; i<N + 1; i++)
					{
						for (int j = 0; j<L + 1; j++)
						{
							for (int k = 0; k<M + 1; k++)
							{
								int mijk = mask[i][j][k];
								if (mijk == TFictivePoint)
								{
									vpCur[i][j][k] = 0;
								}
								else
								{
									vpCur[i][j][k] = vpPrev[i][j][k];
								}
							} //by k
						} //by j
					} //by i
				} // by vp-object number
				//ReAssigning (End)

				//Class fields (End)


				//Previous time step (Begin)
				if (u_1Ptr != NULL) { delete u_1Ptr; }
				u_1Ptr = dynamic_cast<TRnRelease3dSpace*>(&(fU->Copy()));
				if (v_1Ptr != NULL) { delete v_1Ptr; }
				v_1Ptr = dynamic_cast<TRnRelease3dSpace*>(&(fV->Copy()));
				if (w_1Ptr != NULL) { delete w_1Ptr; }
				w_1Ptr = dynamic_cast<TRnRelease3dSpace*>(&(fW->Copy()));
				//Previous time step (End)

				//Intermediate (Begin)
				if (u_Ptr != NULL) { delete u_Ptr; }
				u_Ptr = dynamic_cast<TRnRelease3dSpace*>(&(fU->Copy()));
				if (v_Ptr != NULL) { delete v_Ptr; }
				v_Ptr = dynamic_cast<TRnRelease3dSpace*>(&(fV->Copy()));
				if (w_Ptr != NULL) { delete w_Ptr; }
				w_Ptr = dynamic_cast<TRnRelease3dSpace*>(&(fW->Copy()));
				//Intermediate (End)

				//For control stating (Begin)
				if (uSPtr != NULL) { delete uSPtr; }
				uSPtr = dynamic_cast<TRnRelease3dSpace*>(&(fU->Copy()));
				(*uSPtr) |= 0;
				if (vSPtr != NULL) { delete vSPtr; }
				vSPtr = dynamic_cast<TRnRelease3dSpace*>(&(fV->Copy()));
				(*vSPtr) |= 0;
				if (wSPtr != NULL) { delete wSPtr; }
				wSPtr = dynamic_cast<TRnRelease3dSpace*>(&(fW->Copy()));
				(*wSPtr) |= 0;
				//For control stating (End)

				//Poisson equation right part (Begin)
				if (rightPartPtr != NULL) { delete rightPartPtr; }
				rightPartPtr = dynamic_cast<TRnRelease3dSpace*>(&(fP->Copy()));
				(*rightPartPtr) |= 0;
				//Poisson equation right part (End)


				//Constructing pressure operator (Begin)
				{
					TNodesGrid** grids = new TNodesGrid*[1];
					grids[0] = fGridP;

					if (APtr != NULL) { delete APtr; }
					APtr = new T3dVariableDivergentedOperator(fSpaceModelNameP, 1, grids);

					delete[] grids; grids = NULL;
				}



				if (usePackOperator == true)
				{
					int RN = fP->EndIndex;

					//Model Rn (Begin)
					T1dNumberMask& spaceMaskRn = *(new T1dNumberMask(RN));
					for (int rn = 0; rn<RN; rn++) { spaceMaskRn[rn] = TActualPoint; }
					TMask** spaceMasksRn = new TMask*[1];
					spaceMasksRn[0] = &spaceMaskRn;
					TRnRelease1dSpaceModel& spaceModelRn = *(new TRnRelease1dSpaceModel(1, spaceMasksRn));


					int spaceModelNameRn_Old = spaceModelNameRn;
					if (spaceModelNameRn == 1056) { spaceModelNameRn = spaceModelNameRn - 1; }
					else { spaceModelNameRn = spaceModelNameRn + 1; }

					TRnSpace::RemoveModel(spaceModelNameRn_Old);

					TRnSpace::AddModel(spaceModelNameRn, &spaceModelRn, false);

					delete &spaceMaskRn;
					delete[] spaceMasksRn;
					delete &spaceModelRn;
					//Model Rn (End)

					//Pack operator (Begin)
					printf("Constructing pack operator: Starting...\n");
					double* nullVector = NULL;

					if (ARnPtr != NULL) { delete ARnPtr; }
					ARnPtr = new T1dPackOperator(spaceModelNameRn, *APtr, &nullVector);

					printf("Constructing pack operator: OK.\n");
					//Pack operator (End)
				}
				else
				{
					//printf("Constracting matrix: Starting...");
					//printf("\n");

					//APtr->ConstructMatrix();

					//printf("Constracting matrix: OK.");
					//printf("\n");
				}
				//Constructing pressure operator (End)


				printf("Deleting old VP objects: Starting....");
				printf("\n");
				delete u_Old;
				delete v_Old;
				delete w_Old;
				delete p_Old;

				TRnSpace::RemoveModel(spaceModelNameU_Old);
				TRnSpace::RemoveModel(spaceModelNameV_Old);
				TRnSpace::RemoveModel(spaceModelNameW_Old);
				TRnSpace::RemoveModel(spaceModelNameP_Old);
				printf("Deleting old VP objects: OK.\n");

				//Concentration (Begin)
				{
					if ((C_1Ptr != NULL) || (CPtr != NULL))
					{
						throw "Error in TVPProblem::SolveBySplitting: Concentration does not support restructing.";
					}

					T3dNormalGrid& grid = dynamic_cast<T3dNormalGrid&>(*fGridP);
					TSeparator* sep1 = const_cast<TSeparator*>(&(grid.GetSeparator1()));
					TSeparator* sep2 = const_cast<TSeparator*>(&(grid.GetSeparator2()));
					TSeparator* sep3 = const_cast<TSeparator*>(&(grid.GetSeparator3()));
					T3dVector *const *const *const * normals = const_cast<T3dVector *const *const *const * >(grid.GetNormals());

					T3dNumberMask& mask = grid.CopyMask();
					//Mask modification for grid (Begin)
					int N = sep1->EndIndex;
					int L = sep2->EndIndex;
					int M = sep3->EndIndex;
					for (int i = 0; i<N + 1; i++)
					{
						for (int j = 0; j<L + 1; j++)
						{
							for (int k = 0; k<M + 1; k++)
							{
								if (i == 0 && ConcentrationInletMask(i, j, k, N, L, M))
								{
									mask[i][j][k] = TPreDefinedBorderPoint;
									continue;
								}
								else if (!ConcentrationInletMask(i, j, k, N, L, M))
								{
									mask[i][j][k] = TPreNormalBorderPoint;
									continue;
								}

								int mijk = mask[i][j][k];
								if (mijk == TActualPoint) { mask[i][j][k] = mijk; }
								else if (mijk == TPreDefinedBorderPoint)
								{
									if (i == 0 || j == 0 || k == 0) { mask[i][j][k] = mijk; }
									else if (i == N || j == L || k == M) { mask[i][j][k] = TPreNormalBorderPoint; }
									else
									{
										printf("N = %d, L = %d, M = %d\n", N, L, M);
										printf("i = %d, j = %d, k = %d\n", i, j, k);
										throw "Error in TVPProblem::SolveBySplitting: Unknown concentration PreDefined point mask type.";
									}
								}
								else if (mijk == TFictivePoint)
								{
									if (((j == 0) || (j == L)) && (i != 0) && (i != N) && (k != 0) && (k != M))
									{
										mask[i][j][k] = TPreNormalBorderPoint;
									}
									else if (((k == 0) || (k == M)) && (i != 0) && (i != N) && (j != 0) && (j != L))
									{
										mask[i][j][k] = TPreNormalBorderPoint;
									}
									else if (((i == 0) || (i == N)) && (k != 0) && (k != M) && (j != 0) && (j != L))
									{
										mask[i][j][k] = TPreNormalBorderPoint;
									}
									else { mask[i][j][k] = mijk; }
								}
								else
								{
									throw "Error in TVPProblem::SolveBySplitting: Unknown concentration point mask type.";
								}
							} //by k
						} //by j
					} //by i
					  //Mask modification for grid (End)
					gridCPtr = new T3dNormalGrid(sep1, sep2, sep3, mask, normals, false);

					//Mask modification for space (Begin)
					for (int i = 0; i<N + 1; i++)
					{
						for (int j = 0; j<L + 1; j++)
						{
							for (int k = 0; k<M + 1; k++)
							{
								int mijk = mask[i][j][k];
								if ((mijk == TActualPoint) || (mijk == TFictivePoint)) { mask[i][j][k] = mijk; }
								else if (mijk == TPreDefinedBorderPoint)
								{
									mask[i][j][k] = TFictivePoint;
								}
								else if (mijk == TPreNormalBorderPoint)
								{
									mask[i][j][k] = TActualPoint;
								}
								else
								{
									throw "Error in TVPProblem::SolveBySplitting: Unknown concentration point mask type.";
								}
							}
						}
					}
					//Mask modification for space (End)

					TMask** masks = new TMask*[1];
					masks[0] = &mask;

					int spaceModelName = 55555;

					TRnRelease3dSpaceModel& spaceModel = *(new TRnRelease3dSpaceModel(1, masks));
					TRnSpace::AddModel(spaceModelName, &spaceModel, false);

					delete &mask;
					delete[] masks;
					delete &spaceModel;

					C_1Ptr = new TRnRelease3dSpace(spaceModelName);
					CPtr = new TRnRelease3dSpace(spaceModelName);

					(*C_1Ptr) |= 0;
					(*CPtr) |= 0;

					//Initial conditions (Begin)
					for (int i = 0; i<N + 1; i++)
					{
						for (int j = 0; j<L + 1; j++)
						{
							for (int k = 0; k<M + 1; k++)
							{
								if (InitialConcentrationDistribution(i, j, k, N, L, M))
								{
									(*C_1Ptr)[i][j][k] = 0.45;
									(*CPtr)[i][j][k] = 0.45;
								}
							}
						}
					}

					//Initial conditions (End)
				}
				//Concentration (End)

				//Dencity (Begin)
				{
					if (densityPtr != NULL)
					{
						throw "Error in TVPProblem::SolveBySplitting: Density does not support restructing.";
					}

					T3dNormalGrid& grid = dynamic_cast<T3dNormalGrid&>(*fGridP);
					TSeparator* sep1 = const_cast<TSeparator*>(&(grid.GetSeparator1()));
					TSeparator* sep2 = const_cast<TSeparator*>(&(grid.GetSeparator2()));
					TSeparator* sep3 = const_cast<TSeparator*>(&(grid.GetSeparator3()));
					T3dVector *const *const *const * normals = const_cast<T3dVector *const *const *const * >(grid.GetNormals());

					T3dNumberMask& mask = grid.CopyMask();

					//Mask modification for grid (Begin)
					int N = sep1->EndIndex;
					int L = sep2->EndIndex;
					int M = sep3->EndIndex;
					printf("\nDENSITY:   N=%d, L=%d, M=%d\n", N, L, M);

					for (int i = 0; i < N + 1; i++)
					{
						for (int j = 0; j < L + 1; j++)
						{
							for (int k = 0; k < M + 1; k++)
							{
								int mijk = mask[i][j][k];
								if (mijk == TActualPoint)
								{
									mask[i][j][k] = mijk;
								}
								else if (mijk == TPreDefinedBorderPoint)
								{
									if (i == 0 || j == 0 || k == 0)
									{
										mask[i][j][k] = mijk;
									}
									else if (i == N || j == L || k == M)
									{
										mask[i][j][k] = TPreNormalBorderPoint;
									}
									else
									{
										printf("N = %d, L = %d, M = %d\n", N, L, M);
										printf("i = %d, j = %d, k = %d\n", i, j, k);
										throw "Error in TVPProblem::SolveBySplitting: Unknown concentration PreDefined point mask type.";
									}
								}
								else if (mijk == TFictivePoint)
								{
									if (((j == 0) || (j == L)) && (i != 0) && (i != N) && (k != 0) && (k != M))
									{
										mask[i][j][k] = TPreNormalBorderPoint;
									}
									else if (((k == 0) || (k == M)) && (i != 0) && (i != N) && (j != 0) && (j != L))
									{
										mask[i][j][k] = TPreNormalBorderPoint;
									}
									else if (((i == 0) || (i == N)) && (k != 0) && (k != M) && (j != 0) && (j != L))
									{
										mask[i][j][k] = TPreNormalBorderPoint;
									}
									else
									{
										mask[i][j][k] = mijk;
									}
								}
								else
								{
									throw "Error in TVPProblem::SolveBySplitting: Unknown density point mask type.";
								}
							}
						}
					}

					gridDensityPtr = new T3dNormalGrid(sep1, sep2, sep3, mask, normals, false);

					for (int i = 0; i < N + 1; i++)
					{
						for (int j = 0; j < L + 1; j++)
						{
							for (int k = 0; k < M + 1; k++)
							{
								int mijk = mask[i][j][k];
								if ((mijk == TActualPoint) || (mijk == TFictivePoint))
								{
									mask[i][j][k] = mijk;
								}
								else if (mijk == TPreDefinedBorderPoint)
								{
									mask[i][j][k] = TFictivePoint;
								}
								else if (mijk == TPreNormalBorderPoint)
								{
									mask[i][j][k] = TActualPoint;
								}
								else
								{
									throw "Error in TVPProblem::SolveBySplitting: Unknown density point mask type.";
								}
							}
						}
					}

					TMask** masks = new TMask*[1];
					masks[0] = &mask;

					int spaceModelName = 88888;

					TRnRelease3dSpaceModel& spaceModel = *(new TRnRelease3dSpaceModel(1, masks));
					TRnSpace::AddModel(spaceModelName, &spaceModel, false);

					delete &mask;
					delete[] masks;
					delete &spaceModel;

					densityPtr = new TRnRelease3dSpace(spaceModelName);

					(*densityPtr) |= 0;


					//Initial conditions (Begin)
					for (int i = 0; i < N + 1; i++)
					{
						for (int j = 0; j < L + 1; j++)
						{
							for (int k = 0; k < M + 1; k++)
							{
								//int mijk = (*gridDensityPtr).GetMask()[i][j][k];
								//if (mijk == TActualPoint) { (*densityPtr)[i][j][k] = (*C_1Ptr)[i][j][k]*(fRho2 - fRho1) + fRho1; }
								//else if (mijk == TFictivePoint) { (*densityPtr)[i][j][k] = 0; }
								//else { throw "Error in TVPProblem::SolveBySplitting: Unknown density point mask type."; }
								(*densityPtr)[i][j][k] = (*CPtr)[i][j][k]*(fRho2 - fRho1) + fRho1;
							}
						}
					}
					/*
					for (int i = 0; i < N + 1; i++)
					{
						for (int j = 0; j < L; j++)
						{
							for (int k = 0; k < M; k++)
							{
								if (InitialDensityDistribution(i, j, k, N, L, M))
								{
									//r_ = GetDistribution(i, j, k, hX_, hY_, hZ_);
									//distr = CountDistribution(r_);
									(*densityPtr)[i][j][k] = soluteDensityValue; // +(1.00 - distr)*envDensityValue;
								}
							}
						}
					}
					*/
				}
				//Dencity (End)
				
				//Viscosity (Begin)
				{
					if ((visPtr != NULL))
					{
						throw "Error in TVPProblem::SolveBySplitting: Viscosity does not support restructing.";
					}

					T3dNormalGrid& grid = dynamic_cast<T3dNormalGrid&>(*fGridP);
					TSeparator* sep1 = const_cast<TSeparator*>(&(grid.GetSeparator1()));
					TSeparator* sep2 = const_cast<TSeparator*>(&(grid.GetSeparator2()));
					TSeparator* sep3 = const_cast<TSeparator*>(&(grid.GetSeparator3()));
					T3dVector *const *const *const * normals = const_cast<T3dVector *const *const *const * >(grid.GetNormals());

					T3dNumberMask& mask = grid.CopyMask();
					//Mask modification for grid (Begin)
					int N = sep1->EndIndex;
					int L = sep2->EndIndex;
					int M = sep3->EndIndex;
					for (int i = 0; i<N + 1; i++)
					{
						for (int j = 0; j<L + 1; j++)
						{
							for (int k = 0; k<M + 1; k++)
							{
								int mijk = mask[i][j][k];
								if (mijk == TActualPoint) { mask[i][j][k] = mijk; }
								else if (mijk == TPreDefinedBorderPoint) { mask[i][j][k] = TFictivePoint; }
								else if (mijk == TFictivePoint) { mask[i][j][k] = mijk; }
								else
								{
									throw "Error in TVPProblem::SolveBySplitting: Unknown viscosity point mask type.";
								}
							} //by k
						} //by j
					} //by i
					  //Mask modification for grid (End)
					gridVisPtr = new T3dNormalGrid(sep1, sep2, sep3, mask, normals, false);

					//Mask modification for space (Begin)
					for (int i = 0; i<N + 1; i++)
					{
						for (int j = 0; j<L + 1; j++)
						{
							for (int k = 0; k<M + 1; k++)
							{
								int mijk = mask[i][j][k];
								if ((mijk == TActualPoint) || (mijk == TFictivePoint)) { mask[i][j][k] = mijk; }
								else
								{
									throw "Error in TVPProblem::SolveBySplitting: Unknown viscosity point mask type.";
								}
							}
						}
					}
					//Mask modification for space (End)

					TMask** masks = new TMask*[1];
					masks[0] = &mask;

					int spaceModelName = 77777;

					TRnRelease3dSpaceModel& spaceModel = *(new TRnRelease3dSpaceModel(1, masks));
					TRnSpace::AddModel(spaceModelName, &spaceModel, false);

					delete &mask;
					delete[] masks;
					delete &spaceModel;

					visPtr = new TRnRelease3dSpace(spaceModelName);

					(*visPtr) |= 0;

					//Initial conditions (Begin)
					for (int i = 0; i<N + 1; i++)
					{
						for (int j = 0; j<L + 1; j++)
						{
							for (int k = 0; k<M + 1; k++)
							{
								//int mijk = (*gridVisPtr).GetMask()[i][j][k];
								//if (mijk == TActualPoint) { (*visPtr)[i][j][k] = (*C_1Ptr)[i][j][k]*(fMu2 - fMu1) + fMu1; }
								//else if (mijk == TFictivePoint) { (*visPtr)[i][j][k] = 0; }
								//else { throw "Error in TVPProblem::SolveBySplitting: Unknown viscosity point mask type."; }
								(*visPtr)[i][j][k] = (*CPtr)[i][j][k]*(fMu2 - fMu1) + fMu1;
							}
						}
					}

					const double* hx = grid.GetSeparator1().Separation;
					const double* hy = grid.GetSeparator2().Separation;
					const double* hz = grid.GetSeparator2().Separation;

					double hX_ = hx[1], hY_ = hy[1], hZ_ = hz[1];
					double r_ = 0, distr = 0.;
					/*
					for (int i = 0; i<N + 1; i++)
					{
						for (int j = 0; j < L; j++)
						{
							for (int k = 0; k < M; k++)
							{
								if (InitialDensityDistribution(i, j, k, N, L, M))
								{
									(*visPtr)[i][j][k] = 1 / fRe;
								}
							}
						}
					}
					*/
					//Initial conditions (End)
				}
				//Viscosity (End)			

				
				printf("Restructing: OK.\n");
				printf("####################################\n");
			}//need to restruct
			 //Restructing (End)



			 //Vectors defining (Begin)

			 //Basic - current time step (Begin)
			TRnRelease3dSpace &u = *fU;
			TRnRelease3dSpace &v = *fV;
			TRnRelease3dSpace &w = *fW;
			TRnRelease3dSpace &p = *fP;
			//Basic - current time step (End)


			//Previous time step (Begin)
			TRnRelease3dSpace& u_1 = *u_1Ptr;
			TRnRelease3dSpace& v_1 = *v_1Ptr;
			TRnRelease3dSpace& w_1 = *w_1Ptr;
			//Previous time step (End)


			//Intermediate (Begin)
			TRnRelease3dSpace& u_ = *u_Ptr;
			TRnRelease3dSpace& v_ = *v_Ptr;
			TRnRelease3dSpace& w_ = *w_Ptr;
			//Intermediate (End)


			//For control stating (Begin)
			TRnRelease3dSpace& uS = *uSPtr;
			TRnRelease3dSpace& vS = *vSPtr;
			TRnRelease3dSpace& wS = *wSPtr;
			//For control stating (End)

			//Poisson equation right part (Begin)
			TRnRelease3dSpace& rightPart = *rightPartPtr;
			//Poisson equation right part (End)


			//Concentration (Begin)
			T3dNormalGrid& gridC = *gridCPtr;
			TRnRelease3dSpace& cnc_1 = *C_1Ptr;

			TRnRelease3dSpace& cnc = *CPtr;
			//Concentration (End)

			//Viscosity (Begin)
			T3dNormalGrid& gridVis = *gridVisPtr;
			TRnRelease3dSpace& vis = *visPtr;
			//Viscosity (End)

			T3dNormalGrid& gridDensity = *gridDensityPtr;
			TRnRelease3dSpace& density = *densityPtr;

			//Vectors defining (End)


			TRnRelease3dSpace& rpTrX = dynamic_cast<TRnRelease3dSpace &>(u_.CreateInstance());
			TRnRelease3dSpace& rpTrY = dynamic_cast<TRnRelease3dSpace &>(v_.CreateInstance());
			TRnRelease3dSpace& rpTrZ = dynamic_cast<TRnRelease3dSpace &>(w_.CreateInstance());



			//Velocity equations (Begin)
			printf("\n");
			printf("Velocity equations: Starting...\n");
			start = clock();
	#pragma omp parallel num_threads(3)
			{
	#pragma omp for
				for (int eq = 1; eq<4; eq++)
				{
					printf("\nVelocity equation %d: Starting...\n", eq);

					//Current grid (Begin)
					TNodesGrid* gridPtr;
					if (eq == 1) { gridPtr = fGridU; }
					else if (eq == 2) { gridPtr = fGridV; }
					else if (eq == 3) { gridPtr = fGridW; }
					else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }
					T3dNormalGrid& grid = dynamic_cast<T3dNormalGrid&>(*gridPtr);

					int N = grid.GetSeparator1().EndIndex;
					int L = grid.GetSeparator2().EndIndex;
					int M = grid.GetSeparator3().EndIndex;
					const T3dNumberMask& mask = grid.Mask;
					//Current grid (End)

					//Current vectors (Begin)
					TRnRelease3dSpace* vecPrevPtr;
					TRnRelease3dSpace* vecCurPtr;
					if (eq == 1) { vecPrevPtr = &u_1; vecCurPtr = &u_; }
					else if (eq == 2) { vecPrevPtr = &v_1; vecCurPtr = &v_; }
					else if (eq == 3) { vecPrevPtr = &w_1; vecCurPtr = &w_; }
					else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }
					TRnRelease3dSpace& f_1 = *vecPrevPtr;
					TRnRelease3dSpace& f_ = *vecCurPtr;
					//Current vectors (End)

					//Transfer terms (Begin)
					TRnRelease3dSpace& uTr = dynamic_cast<TRnRelease3dSpace &>(f_1.CreateInstance());
					TRnRelease3dSpace& vTr = dynamic_cast<TRnRelease3dSpace &>(f_1.CreateInstance());
					TRnRelease3dSpace& wTr = dynamic_cast<TRnRelease3dSpace &>(f_1.CreateInstance());
										
					TRnRelease3dSpace* rpTrRef;//= dynamic_cast<TRnRelease3dSpace &>(f_1.CreateInstance());

					if (eq == 1) { rpTrRef = &rpTrX; }
					else if (eq == 2) { rpTrRef = &rpTrY; }
					else if (eq == 3) { rpTrRef = &rpTrZ; }
					else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }

					TRnRelease3dSpace& rpTr = *rpTrRef;

					//Transfer terms (End)


					//Viscosity (Begin)
					TRnRelease3dSpace& visTr = dynamic_cast<TRnRelease3dSpace &>(f_1.CreateInstance());
					//Viscosity (End)


					//Current vector - to zero (Begin)
					//Current vector have to be calculated at all points (including border).
					for (int i = 0; i<N + 1; i++)
					{
						for (int j = 0; j<L + 1; j++)
						{
							for (int k = 0; k<M + 1; k++)
							{
								f_[i][j][k] = 0;
							}
						}
					}
					//Current vector - to zero (End)


					//Convective items, right part, border conditions (Begin)
					for (int i = 0; i<N + 1; i++)
					{
						for (int j = 0; j<L + 1; j++)
						{
							for (int k = 0; k<M + 1; k++)
							{
								int mijk = mask[i][j][k];
								if (mijk == TFictivePoint)
								{
									uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0; rpTr[i][j][k] = 0;
									//Viscosity (Begin)
									visTr[i][j][k] = 0;
									//Viscosity (End)
								}
								else if (mijk == TActualPoint)
								{
									if (eq == 1)
									{
										uTr[i][j][k] = u_1[i][j][k];
										vTr[i][j][k] = (v_1[i][j][k] + v_1[i + 1][j][k] + v_1[i][j - 1][k] + v_1[i + 1][j - 1][k]) / 4;
										wTr[i][j][k] = (w_1[i][j][k] + w_1[i + 1][j][k] + w_1[i][j][k - 1] + w_1[i + 1][j][k - 1]) / 4;
										//Viscosity (Begin)
										visTr[i][j][k] = (vis[i][j][k] + vis[i + 1][j][k]) / 2;
										//Viscosity (End)
									}
									else if (eq == 2)
									{
										uTr[i][j][k] = (u_1[i][j][k] + u_1[i][j + 1][k] + u_1[i - 1][j][k] + u_1[i - 1][j + 1][k]) / 4;
										vTr[i][j][k] = v_1[i][j][k];
										wTr[i][j][k] = (w_1[i][j][k] + w_1[i][j + 1][k] + w_1[i][j][k - 1] + w_1[i][j + 1][k - 1]) / 4;
										//Viscosity (Begin)
										visTr[i][j][k] = (vis[i][j][k] + vis[i][j + 1][k]) / 2;
										//Viscosity (End)
									}
									else if (eq == 3)
									{
										uTr[i][j][k] = (u_1[i][j][k] + u_1[i][j][k + 1] + u_1[i - 1][j][k] + u_1[i - 1][j][k + 1]) / 4;
										vTr[i][j][k] = (v_1[i][j][k] + v_1[i][j][k + 1] + v_1[i][j - 1][k] + v_1[i][j - 1][k + 1]) / 4;
										wTr[i][j][k] = w_1[i][j][k];
										//Viscosity (Begin)
										visTr[i][j][k] = (vis[i][j][k] + vis[i][j][k + 1]) / 2;
										//Viscosity (End)
									}
									else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }

									rpTr[i][j][k] = 0; //no external actions
								}
								else if (mijk == TDefinedBorderPoint)
								{
									uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0;
									//Viscosity (Begin)
									visTr[i][j][k] = 0;
									//Viscosity (End)
									double defVal;
									if (eq == 1) { defVal = GetBorderConditionU(i, j, k, timeStepNumber); }
									else if (eq == 2) { defVal = GetBorderConditionV(i, j, k, timeStepNumber); }
									else if (eq == 3) { defVal = GetBorderConditionW(i, j, k, timeStepNumber); }
									else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }
									rpTr[i][j][k] = defVal;
								}
								else if (mijk == TPreDefinedBorderPoint)
								{
									uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0;
									//Viscosity (Begin)
									visTr[i][j][k] = 0;
									//Viscosity (End)
									double defVal;
									if (eq == 1) { defVal = GetBorderConditionU(i, j, k, timeStepNumber); }
									else if (eq == 2) { defVal = GetBorderConditionV(i, j, k, timeStepNumber); }
									else if (eq == 3) { defVal = GetBorderConditionW(i, j, k, timeStepNumber); }
									else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }
									rpTr[i][j][k] = defVal;
								}
								else if (mijk == TNormalBorderPoint)
								{
									//throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity.";

									double nx = (grid.Normals[i][j][k])->X;
									double ny = (grid.Normals[i][j][k])->Y;
									double nz = (grid.Normals[i][j][k])->Z;

									uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0;
									//Viscosity (Begin)
									visTr[i][j][k] = 0;
									//Viscosity (End)


									double defVal;
									if (eq == 1)
									{
										if ((nx == 0) || (ny != 0) || (nz != 0))
										{
											throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity and this equation number and this normal vector.";
										}
										defVal = GetBorderConditionU(i, j, k, timeStepNumber);
									}
									else if (eq == 2)
									{
										//defVal = GetBorderConditionV(i,j,k,timeStepNumber);
										// NOTE: Fixed here for complex conditions
										throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity and this equation number.";
									}
									else if (eq == 3)
									{
										//defVal = GetBorderConditionW(i,j,k,timeStepNumber);
										// NOTE: Fixed here for complex conditions
										throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity and this equation number.";
									}
									else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }
									rpTr[i][j][k] = defVal;
								}
								else if (mijk == TEquationBorderPoint)
								{
									throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity.";

									//Viscosity (Begin)
									visTr[i][j][k] = 0;
									//Viscosity (End)

									double nx = (grid.Normals[i][j][k])->X;
									double ny = (grid.Normals[i][j][k])->Y;
									double nz = (grid.Normals[i][j][k])->Z;

									const double* hx = grid.GetSeparator1().Separation;
									const double* hy = grid.GetSeparator2().Separation;
									const double* hz = grid.GetSeparator2().Separation;

									T3dNormalGrid& gridP = dynamic_cast<T3dNormalGrid&>(*fGridP);
									const double* hxP = gridP.GetSeparator1().Separation;
									const double* hyP = gridP.GetSeparator2().Separation;
									const double* hzP = gridP.GetSeparator3().Separation;


									double h_1X = hx[i - 1];
									double hX = hx[i];
									double  hx2X = hx[i - 1] + hx[i];
									double  hxhxhX = hx[i - 1] * hx[i] * (hx[i - 1] + hx[i]) / 2;

									double h_1Y = hy[j - 1];
									double hY = hy[j];
									double hx2Y = hy[j - 1] + hy[j];
									double hxhxhY = hy[j - 1] * hy[j] * (hy[j - 1] + hy[j]) / 2;

									double h_1Z = hy[k - 1];
									double hZ = hy[k];
									double hx2Z = hy[k - 1] + hy[k];
									double hxhxhZ = hy[k - 1] * hy[k] * (hy[k - 1] + hy[k]) / 2;


									double eqVal = 0;

									if (eq == 1)
									{
										if ((nx == 0) || (ny != 0) || (nz != 0))
										{
											throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity and this equation number and this normal vector.";
										}

										double hxhxhUpX = hx[i] * hx[i + 1] * (hx[i] + hx[i + 1]) / 2;
										double hxhxhDownX = hx[i - 2] * hx[i - 1] * (hx[i - 2] + hx[i - 1]) / 2;

										vTr[i][j][k] = (v_1[i][j][k] + v_1[i + 1][j][k] + v_1[i][j - 1][k] + v_1[i + 1][j - 1][k]) / 4;
										wTr[i][j][k] = (w_1[i][j][k] + w_1[i + 1][j][k] + w_1[i][j][k - 1] + w_1[i + 1][j][k - 1]) / 4;


										eqVal = eqVal - (vTr[i][j][k])*(u_1[i][j + 1][k] - u_1[i][j - 1][k]) / (hx2Y);
										//eqVal = eqVal + (1. / fRe)*(hY*u_1[i][j - 1][k] - hx2Y*u_1[i][j][k] + h_1Y*u_1[i][j + 1][k]) / hxhxhY;
										eqVal = eqVal + vis[i][j][k]*(hY*u_1[i][j - 1][k] - hx2Y*u_1[i][j][k] + h_1Y*u_1[i][j + 1][k]) / hxhxhY;

										eqVal = eqVal - (wTr[i][j][k])*(u_1[i][j][k + 1] - u_1[i][j][k - 1]) / (hx2Z);
										//eqVal = eqVal + (1. / fRe)*(hZ*u_1[i][j][k - 1] - hx2Z*u_1[i][j][k] + h_1Z*u_1[i][j][k + 1]) / hxhxhZ;
										eqVal = eqVal + vis[i][j][k]*(hZ*u_1[i][j][k - 1] - hx2Z*u_1[i][j][k] + h_1Z*u_1[i][j][k + 1]) / hxhxhZ;


										//eqVal = eqVal - (p[i+1][j][k] - p[i][j][k])/hxP[i];


										if (nx<0)
										{
											eqVal = eqVal - (u_1[i][j][k])*(u_1[i + 1][j][k] - u_1[i][j][k]) / (hx[i]);
											//eqVal = eqVal + (1. / fRe)*(hx[i + 1] * u_1[i][j][k] - (hx[i] + hx[i + 1])*u_1[i + 1][j][k] + hx[i] * u_1[i + 2][j][k]) / hxhxhUpX;
											eqVal = eqVal + vis[i][j][k]*(hx[i + 1] * u_1[i][j][k] - (hx[i] + hx[i + 1])*u_1[i + 1][j][k] + hx[i] * u_1[i + 2][j][k]) / hxhxhUpX;
										}
										else
										{
											eqVal = eqVal - (u_1[i][j][k])*(u_1[i][j][k] - u_1[i - 1][j][k]) / (hx[i - 1]);
											//eqVal = eqVal + (1. / fRe)*(hx[i - 1] * u_1[i - 2][j][k] - (hx[i - 2] + hx[i - 1])*u_1[i - 1][j][k] + hx[i - 2] * u_1[i][j][k]) / hxhxhDownX;
											eqVal = eqVal + vis[i][j][k]*(hx[i - 1] * u_1[i - 2][j][k] - (hx[i - 2] + hx[i - 1])*u_1[i - 1][j][k] + hx[i - 2] * u_1[i][j][k]) / hxhxhDownX;
										}

										eqVal = tau*eqVal + u_1[i][j][k];
									}
									else if (eq == 2)
									{
										throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity and this equation number.";
									}
									else if (eq == 3)
									{
										throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity and this equation number.";
									}
									else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }
									uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0;
									rpTr[i][j][k] = eqVal;
								}
								else if ((mijk == TPreNormalBorderPoint) || (mijk == TPreEquationBorderPoint))
								{
									throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity.";
								}
								else
								{
									throw "Error in TVPProblem::SolveBySplitting: Unknown point mask type.";
								}
							} //by k
						} //by j
					} //by i

					TRnRelease3dSpace& rpTrTmp = dynamic_cast<TRnRelease3dSpace &>(f_1.CreateInstance());

					const T3dNumberMask& maskC = (*gridCPtr).Mask;
					for (int i = 1; i<N; i++)
					{
						for (int j = 1; j<L; j++)
						{
							for (int k = 1; k<M; k++)
							{
								maskC[i][j][k] = TActualPoint;
							}
						}
					}

					// eq is equal to TImmersedBoundary::COORD_*
					// +f
					if (eq == 1) {
						GetForce(rpTrTmp, *gridCPtr, grid, uTr, eq, timeStepNumber);
						//InterpolateUpdatePathes(grid, uTr, eq, timeStepNumber);
					}
					else if (eq == 2) {
						GetForce(rpTrTmp, *gridCPtr, grid, vTr, eq, timeStepNumber);
						//InterpolateUpdatePathes(grid, vTr, eq, timeStepNumber);
					}
					else if (eq == 3) {
						GetForce(rpTrTmp, *gridCPtr, grid, wTr, eq, timeStepNumber);
						//InterpolateUpdatePathes(grid, wTr, eq, timeStepNumber);
					}

					for (int i = 0; i<N + 1; i++)
					{
						for (int j = 0; j<L + 1; j++)
						{
							for (int k = 0; k<M + 1; k++)
							{
								int mijk = mask[i][j][k];
								if (mijk == TActualPoint)
								{
									rpTr[i][j][k] = rpTrTmp[i][j][k];
								}
							}
						}
					}


					const double* hx = grid.GetSeparator1().Separation;
					const double* hy = grid.GetSeparator2().Separation;
					const double* hz = grid.GetSeparator3().Separation;

					double maxDVis = 0;
					//-nabla*sigma
					for (int i = 0; i<N + 1; i++)
					{
						for (int j = 0; j<L + 1; j++)
						{
							for (int k = 0; k<M + 1; k++)
							{
								int mijk = mask[i][j][k];
								if (mijk == TActualPoint)
								{
									double dVisX = (visTr[i + 1][j][k] - visTr[i - 1][j][k]) / (hx[i] + hx[i - 1]);
									double dVisY = (visTr[i][j + 1][k] - visTr[i][j - 1][k]) / (hy[j] + hy[j - 1]);
									double dVisZ = (visTr[i][j][k + 1] - visTr[i][j][k - 1]) / (hz[k] + hz[k - 1]);

									if (fabs(dVisX)>maxDVis) { maxDVis = fabs(dVisX); }
									if (fabs(dVisY)>maxDVis) { maxDVis = fabs(dVisY); }
									if (fabs(dVisZ)>maxDVis) { maxDVis = fabs(dVisZ); }

									double dU = 0;
									double dV = 0;
									double dW = 0;

									TRnRelease3dSpace* velTrPtr = NULL;
									if (eq == 1)
									{
										velTrPtr = &uTr;
										dU = (uTr[i + 1][j][k] - uTr[i - 1][j][k]) / (hx[i] + hx[i - 1]);
										dV = (vTr[i + 1][j][k] - vTr[i - 1][j][k]) / (hx[i] + hx[i - 1]);
										dW = (wTr[i + 1][j][k] - wTr[i - 1][j][k]) / (hx[i] + hx[i - 1]);

										//visTr[i][j][k] = (density[i][j][k] + density[i + 1][j][k]) / 2 * visTr[i][j][k];
									}
									else if (eq == 2)
									{
										velTrPtr = &vTr;
										dU = (uTr[i][j + 1][k] - uTr[i][j - 1][k]) / (hy[j] + hy[j - 1]);
										dV = (vTr[i][j + 1][k] - vTr[i][j - 1][k]) / (hy[j] + hy[j - 1]);
										dW = (wTr[i][j + 1][k] - wTr[i][j - 1][k]) / (hy[j] + hy[j - 1]);

										//visTr[i][j][k] = (density[i][j][k] + density[i][j + 1][k]) / 2 * visTr[i][j][k];
									}
									else if (eq == 3)
									{
										velTrPtr = &wTr;
										dU = (uTr[i][j][k + 1] - uTr[i][j][k - 1]) / (hz[k] + hz[k - 1]);
										dV = (vTr[i][j][k + 1] - vTr[i][j][k - 1]) / (hz[k] + hz[k - 1]);
										dW = (wTr[i][j][k + 1] - wTr[i][j][k - 1]) / (hz[k] + hz[k - 1]);

										//visTr[i][j][k] = (density[i][j][k] + density[i][j][k + 1]) / 2 * visTr[i][j][k];
									}
									else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }
									TRnRelease3dSpace& velTr = *velTrPtr;

									double dVelX = (velTr[i + 1][j][k] - velTr[i - 1][j][k]) / (hx[i] + hx[i - 1]);
									double dVelY = (velTr[i][j + 1][k] - velTr[i][j - 1][k]) / (hy[j] + hy[j - 1]);
									double dVelZ = (velTr[i][j][k + 1] - velTr[i][j][k - 1]) / (hz[k] + hz[k - 1]);


									rpTr[i][j][k] = rpTr[i][j][k] + (dVisX*(dVelX + dU) + dVisY*(dVelY + dV) + dVisZ*(dVelZ + dW));
								}
							}
						}
					}
					
					SolveTransferEquation__(f_, f_1, uTr, vTr, wTr, rpTr, grid, visTr, 0);



					delete &uTr;
					delete &vTr;
					delete &wTr;
					delete &rpTrTmp;


					//Viscosity (Begin)
					delete &visTr;
					//Viscosity (End)


					printf("Velocity equation %d: OK.\n", eq);
				} // by equations


			}////#pragma omp parallel

			end = clock();
			diff = (double)(end - start);
			diff = diff / CLOCKS_PER_SEC;
			printf("Velocity equations: OK. (time = %lf s.)\n", diff);
			//Velocity equations (End)





			//Pressure equation (Begin)
			{
				printf("\n");
				printf("Pressure equation: Starting...\n");

				//Preparing variables (Begin)
				T3dNormalGrid& grid = dynamic_cast<T3dNormalGrid&>(*fGridP);
				int N = grid.GetSeparator1().EndIndex;
				int L = grid.GetSeparator2().EndIndex;
				int M = grid.GetSeparator3().EndIndex;
				const T3dNumberMask& mask = grid.Mask;

				T3dNormalGrid& gridU = dynamic_cast<T3dNormalGrid&>(*fGridU);
				const double* hxU = gridU.GetSeparator1().Separation;
				T3dNormalGrid& gridV = dynamic_cast<T3dNormalGrid&>(*fGridV);
				const double* hyV = gridV.GetSeparator2().Separation;
				T3dNormalGrid& gridW = dynamic_cast<T3dNormalGrid&>(*fGridW);
				const double* hzW = gridW.GetSeparator3().Separation;
				//Preparing variables (End)

				//Right part (Begin)
				
				for (int i = 0; i<N + 1; i++)
				{
					for (int j = 0; j<L + 1; j++)
					{
						for (int k = 0; k<M + 1; k++)
						{
							int mijk = mask[i][j][k];
							if (mijk == TFictivePoint)
							{
								rightPart[i][j][k] = 0;
							}
							else if (mijk == TActualPoint)
							{
								double u_i_1 = u_[i - 1][j][k];
								double u_i = u_[i][j][k];
								double v_j_1 = v_[i][j - 1][k];
								double v_j = v_[i][j][k];
								double w_k_1 = w_[i][j][k - 1];
								double w_k = w_[i][j][k];


								double div_ = density[i][j][k]*density[i][j][k]*((u_i - u_i_1) / hxU[i - 1]) + ((v_j - v_j_1) / hyV[j - 1]) + ((w_k - w_k_1) / hzW[k - 1]);
								rightPart[i][j][k] = -(div_ / tau);
							}
							else if (mijk == TPreDefinedBorderPoint)
							{
								rightPart[i][j][k] = GetBorderConditionP(i, j, k, timeStepNumber);
							}
							else if ((mijk == TPreNormalBorderPoint) || (mijk == TDefinedBorderPoint) || (mijk == TNormalBorderPoint) || (mijk == TEquationBorderPoint) || (mijk == TPreEquationBorderPoint))
							{
								throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this border point type of pressure.";
							}
							else
							{
								throw "Error in TVPProblem::SolveBySplitting: Unknown point mask type.";
							}
						} //by k
					} //by j
				} //by i
				//Right part (End)

				//Poisson problem (Begin)
				if (usePackOperator == true)
				{
					T3dVariableDivergentedOperator& A = *APtr;
					A.withDensity(density);

					int RN = p.EndIndex;

					TRnRelease3dSpace& zeroHelp = *(new TRnRelease3dSpace(fSpaceModelNameP));
					TRnRelease3dSpace& rightPartHelp = *(new TRnRelease3dSpace(fSpaceModelNameP));

					zeroHelp.Assign(p);
					zeroHelp |= 0;
					rightPartHelp |= ((TRnLinearOperator&)A)(zeroHelp) - rightPart;
					rightPartHelp |= (-1)*rightPartHelp;

					delete &zeroHelp;


					TRnRelease1dSpace& rightPartRn = *(new TRnRelease1dSpace(spaceModelNameRn));
					TRnRelease1dSpace& approximationRn = *(new TRnRelease1dSpace(spaceModelNameRn));

					for (int rn = 1; rn<RN + 1; rn++)
					{
						rightPartRn[rn - 1] = ((TRnSpace&)rightPartHelp)[rn];
						approximationRn[rn - 1] = ((TRnSpace&)p)[rn];
					}

					delete &rightPartHelp;


					T1dPackOperator& ARn = *ARnPtr;
					TSLAE& slae = *(new TSLAE(ARn, approximationRn, rightPartRn));
					printf("Solving system: Starting...\n");
					start = clock();
					//slae.SolveByBiConjugateGradientsStab();
					slae.SolveByBiConjugateGradientsStab();
					end = clock();
					diff = (double)(end - start);
					diff = diff / CLOCKS_PER_SEC;
					printf("Solving system: OK. (time = %lf s.)\n", diff);
					approximationRn |= slae.Result;
					delete &slae;

					for (int rn = 1; rn<RN + 1; rn++)
					{
						((TRnSpace&)p)[rn] = approximationRn[rn - 1];
					}

					delete &rightPartRn;
					delete &approximationRn;
				}
				else
				{
					//T3dVariableDivergentedOperator& A = NULL;
					if(timeStepNumber == 1)
					{
						//A = *APtr;
						(*APtr).withDensity(density);
						slae = (new TSLAE((*APtr), p, rightPart));
					}else{
						//slae->UpdateA(*APtr);
						(*APtr).withDensity(density);
						slae->UpdateRP(rightPart);
						slae->UpdateP(p);
					}
					printf("Solving system: Starting...\n");
					start = clock();
					slae->SolveByBiConjugateGradientsStab();
					end = clock();
					diff = (double)(end - start);
					diff = diff / CLOCKS_PER_SEC;
					printf("Solving system: OK. (time = %lf s.)\n", diff);
					p |= slae->Result;
					//delete &slae;
				}//if use pack operator
				// Poisson problem (End)
				printf("Pressure equation: OK.\n");
			}
			//Pressure equation (End)



			//Calculating velocity (Begin)
			printf("\n");
			printf("Velocity calculating: Starting...\n");
			start = clock();
	#pragma omp parallel num_threads(3)
			{
				double tau = fTimeSeparator->SeparationValue;
				T3dNormalGrid& gridP = dynamic_cast<T3dNormalGrid&>(*fGridP);
				const double* hxP = gridP.GetSeparator1().Separation;
				const double* hyP = gridP.GetSeparator2().Separation;
				const double* hzP = gridP.GetSeparator3().Separation;

	#pragma omp for
				for (int vel = 1; vel<4; vel++)
				{
					TNodesGrid* gridPtr;
					if (vel == 1) { gridPtr = fGridU; }
					else if (vel == 2) { gridPtr = fGridV; }
					else if (vel == 3) { gridPtr = fGridW; }
					else { throw "Error in TVPProblem::SolveBySplitting: Unknown velocity number."; }

					T3dNormalGrid& grid = dynamic_cast<T3dNormalGrid&>(*gridPtr);
					int N = grid.GetSeparator1().EndIndex;
					int L = grid.GetSeparator2().EndIndex;
					int M = grid.GetSeparator3().EndIndex;
					const T3dNumberMask& mask = grid.Mask;

					TRnRelease3dSpace* vecPrevPtr;
					TRnRelease3dSpace* vecCurPtr;
					if (vel == 1) { vecPrevPtr = &u_; vecCurPtr = &u; }
					else if (vel == 2) { vecPrevPtr = &v_; vecCurPtr = &v; }
					else if (vel == 3) { vecPrevPtr = &w_; vecCurPtr = &w; }
					else { throw "Error in TVPProblem::SolveBySplitting: Unknown velocity number."; }
					TRnRelease3dSpace& vecPrev = *vecPrevPtr;
					TRnRelease3dSpace& vecCur = *vecCurPtr;


					for (int counterPointType = 1; counterPointType<3; counterPointType++)
					{
						for (int i = 0; i<N + 1; i++)
						{
							for (int j = 0; j<L + 1; j++)
							{
								for (int k = 0; k<M + 1; k++)
								{
									int mijk = mask[i][j][k];

									if (counterPointType == 1)
									{
										if (mijk == TPreDefinedBorderPoint) { continue; }
									}
									else if (counterPointType == 2)
									{
										if (mijk != TPreDefinedBorderPoint) { continue; }
									}
									else
									{
										{throw "Error in TVPProblem::SolveBySplitting: Unknown counter of type point."; }
									}


									if (mijk == TFictivePoint)
									{
										vecCur[i][j][k] = 0;
									}
									else if ((mijk == TActualPoint) || (mijk == TNormalBorderPoint) || (mijk == TEquationBorderPoint))
									{
										double gradP = 0;
										if (vel == 1)
										{
											gradP = (p[i + 1][j][k] - p[i][j][k]) / hxP[i];
										}
										else if (vel == 2)
										{
											gradP = (p[i][j + 1][k] - p[i][j][k]) / hyP[j];
										}
										else if (vel == 3)
										{
											gradP = (p[i][j][k + 1] - p[i][j][k]) / hzP[k];
										}
										else { throw "Error in TVPProblem::SolveBySplitting: Unknown velocity number."; }
										vecCur[i][j][k] = vecPrev[i][j][k] - tau*gradP / density[i][j][k];
									}
									else if (mijk == TDefinedBorderPoint)
									{
										double defVal;
										if (vel == 1) { defVal = GetBorderConditionU(i, j, k, timeStepNumber); }
										else if (vel == 2) { defVal = GetBorderConditionV(i, j, k, timeStepNumber); }
										else if (vel == 3) { defVal = GetBorderConditionW(i, j, k, timeStepNumber); }
										else { throw "Error in TVPProblem::SolveBySplitting: Unknown velocity number."; }

										//Checking (Begin)
										if ((vecPrev[i][j][k]) != defVal)
										{
											throw "Error in TVPProblem::SolveBySplitting: Defined border condition checking error.";
										}
										//Checking (End)
										vecCur[i][j][k] = defVal;
									}
									else if (mijk == TPreDefinedBorderPoint)
									{
										double defVal;
										if (vel == 1) { defVal = GetBorderConditionU(i, j, k, timeStepNumber); }
										else if (vel == 2) { defVal = GetBorderConditionV(i, j, k, timeStepNumber); }
										else if (vel == 3) { defVal = GetBorderConditionW(i, j, k, timeStepNumber); }
										else { throw "Error in TVPProblem::SolveBySplitting: Unknown velocity number."; }

										int i_ = -1;
										int j_ = -1;
										int k_ = -1;
										int actualValuesCount = 0;
										//Defining actual value (Begin)
										for (int dir = 1; dir<4; dir++)
										{
											for (int sign = 0; sign<2; sign++)
											{
												int shift = -1 + sign * 2;
												int i__ = i;
												int j__ = j;
												int k__ = k;
												if (dir == 1) { i__ = i + shift; }
												else if (dir == 2) { j__ = j + shift; }
												else if (dir == 3) { k__ = k + shift; }
												else
												{
													throw "Error in TVPProblem::SolveBySplitting: Unknown direction.";
												}

												if ((i__ >= 0) && (i__ <= N) && (j__ >= 0) && (j__ <= L) && (k__ >= 0) && (k__ <= M))
												{

													if (mask[i__][j__][k__] == TActualPoint)
													{
														i_ = i__; j_ = j__; k_ = k__;
														actualValuesCount = actualValuesCount + 1;
													}
												}
											}
										}
										//Defining actual value (End)

										if (actualValuesCount <= 0)
										{
											throw "Error in TVPProblem::SolveBySplitting: PreDefined border condition checking error - actual values count <= 0.";
										}
										if (actualValuesCount>1)
										{
											throw "Error in TVPProblem::SolveBySplitting: PreDefined border condition checking error - actual values count > 1.";
										}


										//Checking (Begin)
										double eps = 1e-10;
										double check = fabs((vecPrev[i][j][k]) - (2 * defVal - vecPrev[i_][j_][k_]));
										if (check>eps)
										{
											throw "Error in TVPProblem::SolveBySplitting: PreDefined border condition checking error - value incorrect.";
										}
										//Checking (End)
										vecCur[i][j][k] = 2 * defVal - vecCur[i_][j_][k_];
									}
									else if ((mijk == TPreNormalBorderPoint) || (mijk == TPreEquationBorderPoint))
									{
										throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity.";
									}
									else
									{
										throw "Error in TVPProblem::SolveBySplitting: Unknown point mask type.";
									}
								} //by k
							} //by j
						} //by i
					} //by counterPointType
				} // by velocity number

			}////#pragma omp parallel
			end = clock();
			diff = (double)(end - start);
			diff = diff / CLOCKS_PER_SEC;
			printf("Velocity calculating: OK. (time = %lf s.)\n", diff);
			//Calculating velocity (End)



			//Control balance (Begin)
			{
				printf("\n");
				printf("Control balance: Starting...\n");

				//Preparing variables (Begin)
				T3dNormalGrid& grid = dynamic_cast<T3dNormalGrid&>(*fGridP);
				int N = grid.GetSeparator1().EndIndex;
				int L = grid.GetSeparator2().EndIndex;
				int M = grid.GetSeparator3().EndIndex;
				const T3dNumberMask& mask = grid.Mask;

				T3dNormalGrid& gridU = dynamic_cast<T3dNormalGrid&>(*fGridU);
				const double* hxU = gridU.GetSeparator1().Separation;
				T3dNormalGrid& gridV = dynamic_cast<T3dNormalGrid&>(*fGridV);
				const double* hyV = gridV.GetSeparator2().Separation;
				T3dNormalGrid& gridW = dynamic_cast<T3dNormalGrid&>(*fGridW);
				const double* hzW = gridW.GetSeparator3().Separation;
				//Preparing variables (End)

				double divMaxValue = 0;
				int iMaxValue = -1;
				int jMaxValue = -1;
				int kMaxValue = -1;

				for (int i = 0; i<N + 1; i++)
				{
					for (int j = 0; j<L + 1; j++)
					{
						for (int k = 0; k<M + 1; k++)
						{
							int mijk = mask[i][j][k];
							if (mijk == TActualPoint)
							{
								double div = ((u[i][j][k] - u[i - 1][j][k]) / hxU[i - 1]) + ((v[i][j][k] - v[i][j - 1][k]) / hyV[j - 1]) + ((w[i][j][k] - w[i][j][k - 1]) / hzW[k - 1]);
								if (fabs(div)>divMaxValue)
								{
									divMaxValue = fabs(div);
									iMaxValue = i;
									jMaxValue = j;
									kMaxValue = k;
								}
							}
						} //by k
					} //by j
				} //by i
				printf("Control balance: OK. (point = %d-%d-%d, value = %.20f)\n", iMaxValue, jMaxValue, kMaxValue, divMaxValue);
			}
			//Control balance (End)


			//Stating (Begin)
			{
				printf("\n");
				printf("Stating calculating: Starting...\n");
				uS |= u - u_1;
				vS |= v - v_1;
				wS |= w - w_1;
				double uSNorma = uS.GetAbsMaxElement();
				double vSNorma = vS.GetAbsMaxElement();
				double wSNorma = wS.GetAbsMaxElement();
				double sNorma = 0;
				if (uSNorma>sNorma) { sNorma = uSNorma; }
				if (vSNorma>sNorma) { sNorma = vSNorma; }
				if (wSNorma>sNorma) { sNorma = wSNorma; }
				printf("Stating calculating: OK. (value = %.20f)\n", sNorma);
			}
			//Stating (End)




			//Concentration (Begin)
			{
				printf("\n");
				printf("Concentration equation: Starting...\n");
				start = clock();

				int N = gridC.GetSeparator1().EndIndex;
				int L = gridC.GetSeparator2().EndIndex;
				int M = gridC.GetSeparator3().EndIndex;
				const T3dNumberMask& mask = gridC.Mask;


				//Transfer terms (Begin)
				TRnRelease3dSpace& uTr = dynamic_cast<TRnRelease3dSpace &>(cnc.CreateInstance());
				TRnRelease3dSpace& vTr = dynamic_cast<TRnRelease3dSpace &>(cnc.CreateInstance());
				TRnRelease3dSpace& wTr = dynamic_cast<TRnRelease3dSpace &>(cnc.CreateInstance());
				TRnRelease3dSpace& rpTr = dynamic_cast<TRnRelease3dSpace &>(cnc.CreateInstance());
				//Transfer terms (End)


				//Current vector - to zero (Begin)
				//Current vector have to be calculated at all points (including border).
				for (int i = 0; i<N + 1; i++)
				{
					for (int j = 0; j<L + 1; j++)
					{
						for (int k = 0; k<M + 1; k++)
						{
							cnc[i][j][k] = 0;
						}
					}
				}
				//Current vector - to zero (End)


				//Convective items, right part, border conditions (Begin)
				for (int i = 0; i<N + 1; i++)
				{
					for (int j = 0; j<L + 1; j++)
					{
						for (int k = 0; k<M + 1; k++)
						{
							int mijk = mask[i][j][k];
							if (mijk == TFictivePoint)
							{
								uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0; rpTr[i][j][k] = 0;
							}
							else if (mijk == TActualPoint)
							{
								double uu = (u[i - 1][j][k] + u[i][j][k]) / 2;
								uTr[i][j][k] = uu;

								double vv = (v[i][j - 1][k] + v[i][j][k]) / 2;
								vTr[i][j][k] = vv;

								double ww = (w[i][j][k - 1] + w[i][j][k]) / 2;
								wTr[i][j][k] = ww;

								rpTr[i][j][k] = 0; //no external actions
							}
							else if (mijk == TPreDefinedBorderPoint)
							{
								uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0;
								if (i == 0)
								{
									rpTr[i][j][k] = ConcentrationInletCondition(i, j, k, N, L, M);
								}
								else
								{
									rpTr[i][j][k] = 0;
								}

								//rpTr[i][j][k] = 0;
							}
							else if (mijk == TPreNormalBorderPoint)
							{
								uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0;
								rpTr[i][j][k] = 0;
							}
							else if ((mijk == TDefinedBorderPoint) || (mijk == TNormalBorderPoint) || (mijk == TEquationBorderPoint) || (mijk == TPreEquationBorderPoint))
							{
								throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of concentration.";
							}
							else
							{
								throw "Error in TVPProblem::SolveBySplitting: Unknown concentration point mask type.";
							}
						} //by k
					} //by j
				} //by i
				  //Convective items, right part, border conditions (End)

				SolveTransferEquation_(cnc, cnc_1, uTr, vTr, wTr, rpTr, gridC, 0.0);
				//SolveTransferEquation__(cnc, cnc_1, uTr, vTr, wTr, rpTr, gridC, 0.0, 1);
				//(TRnRelease3dSpace& lVecCur, TRnRelease3dSpace& lVecPrev,
		        // TRnRelease3dSpace&  lU, TRnRelease3dSpace& lV, TRnRelease3dSpace& lW, TRnRelease3dSpace& lRP,
				// const TNodesGrid& lGrid, TRnRelease3dSpace& lDiffusionCoefficient, int lSchemeType)
				
				//Recalculation (Begin)
				for (int i = 0; i<N + 1; i++)
				{
					for (int j = 0; j<L + 1; j++)
					{
						for (int k = 0; k<M + 1; k++)
						{
							density[i][j][k] = cnc[i][j][k]*(fRho2 - fRho1) + fRho1;
							vis[i][j][k] = cnc[i][j][k]*(fMu2 - fMu1) + fMu1;
						}
					}
				}

				//Recalculation (End)

				delete &uTr;
				delete &vTr;
				delete &wTr;
				delete &rpTr;


				end = clock();
				diff = (double)(end - start);
				diff = diff / CLOCKS_PER_SEC;
				printf("Concentration equation: OK. (time = %lf s.)\n", diff);
			}
			//Concentration (End)


			//Viscosity (Begin)
			/*
			*    if (turbulenceMode==true)
			*    {
			*    printf("\n");
			*    printf("Turbulence model: Starting...\n");
			*    start = clock();
			*
			*    int N = gridVis.GetSeparator1().EndIndex;
			*    int L = gridVis.GetSeparator2().EndIndex;
			*    int M = gridVis.GetSeparator3().EndIndex;
			*    const T3dNumberMask& mask = gridVis.Mask;
			*
			*    const double* hx = gridVis.GetSeparator1().Separation;
			*    const double* hy = gridVis.GetSeparator2().Separation;
			*    const double* hz = gridVis.GetSeparator3().Separation;
			*
			*
			*    //Velocity projecting (Begin)
			*    TRnRelease3dSpace& uTr  = dynamic_cast<TRnRelease3dSpace &>(vis.CreateInstance());
			*    TRnRelease3dSpace& vTr  = dynamic_cast<TRnRelease3dSpace &>(vis.CreateInstance());
			*    TRnRelease3dSpace& wTr  = dynamic_cast<TRnRelease3dSpace &>(vis.CreateInstance());
			*
			*
			*    for (int i=0;i<N+1;i++)
			*    {
			*      for (int j=0;j<L+1;j++)
			*      {
			*        for (int k=0;k<M+1;k++)
			*        {
			*           int mijk = mask[i][j][k];
			*           if (mijk==TFictivePoint)
			*           {
			*                  uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0;
			*           }
			*           else if (mijk==TActualPoint)
			*           {
			*                 double uu = ( u[i-1][j][k] + u[i][j][k])/2;
			*                 uTr[i][j][k] = uu;
			*
			*                 double vv = ( v[i][j-1][k] + v[i][j][k])/2;
			*                 vTr[i][j][k] = vv;
			*
			*                 double ww = ( w[i][j][k-1] + w[i][j][k])/2;
			*                 wTr[i][j][k] = ww;
			*           }
			*           else if ( (mijk==TPreNormalBorderPoint) || (mijk==TPreDefinedBorderPoint) || (mijk==TDefinedBorderPoint) || (mijk==TNormalBorderPoint) || (mijk==TEquationBorderPoint) || (mijk==TPreEquationBorderPoint) )
			*           {
			*                throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of viscosity.";
			*           }
			*           else
			*           {
			*                throw "Error in TVPProblem::SolveBySplitting: Unknown viscosity point mask type.";
			*           }
			*        } //by k
			*      } //by j
			*    } //by i
			*    //Velocity projecting (End)
			*
			*
			*
			*    //Smogarinsky (Begin)
			*    for (int i=0;i<N+1;i++)
			*    {
			*      for (int j=0;j<L+1;j++)
			*      {
			*        for (int k=0;k<M+1;k++)
			*        {
			*           int mijk = mask[i][j][k];
			*           if (mijk==TFictivePoint)
			*           {
			*                 vis[i][j][k] = 0;
			*           }
			*           else if (mijk==TActualPoint)
			*           {
			*                double Cs = TProblem::GlobalTurbulenceParameter;
			*
			*                double delta = powl(hx[i-1]*hy[j-1]*hz[k-1],1.0/3.0);
			*
			*                double Sxx = (uTr[i+1][j][k] - uTr[i-1][j][k])/(hx[i-1]+hx[i]);
			*                double Syy = (vTr[i][j+1][k] - vTr[i][j-1][k])/(hy[j-1]+hy[j]);
			*                double Szz = (wTr[i][j][k+1] - wTr[i][j][k-1])/(hz[k-1]+hz[k]);
			*                double Sxy = ( (vTr[i+1][j][k] - vTr[i-1][j][k])/(hx[i-1]+hx[i]) + (uTr[i][j+1][k] - uTr[i][j-1][k])/(hy[j-1]+hy[j]) )/2;
			*                double Sxz = ( (wTr[i+1][j][k] - wTr[i-1][j][k])/(hx[i-1]+hx[i]) + (uTr[i][j][k+1] - uTr[i][j][k-1])/(hz[k-1]+hz[k]) )/2;
			*                double Syz = ( (vTr[i][j][k+1] - vTr[i][j][k-1])/(hz[k-1]+hz[k]) + (wTr[i][j+1][k] - wTr[i][j-1][k])/(hy[j-1]+hy[j]) )/2;
			*
			*                double Sc = Sxx*Sxx + Syy*Syy + Szz*Szz + 2*Sxy*Sxy + 2*Sxz*Sxz + 2*Syz*Syz;
			*
			*                double InvS = powl(2*Sc,1.0/2.0);
			*
			*                double vSGS = powl(Cs*delta,2)*InvS;
			*
			*                vis[i][j][k] = 1/fRe + vSGS;
			*           }
			*           else if ( (mijk==TPreNormalBorderPoint) || (mijk==TPreDefinedBorderPoint) || (mijk==TDefinedBorderPoint) || (mijk==TNormalBorderPoint) || (mijk==TEquationBorderPoint) || (mijk==TPreEquationBorderPoint) )
			*           {
			*                throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of viscosity.";
			*           }
			*           else
			*           {
			*                throw "Error in TVPProblem::SolveBySplitting: Unknown viscosity point mask type.";
			*           }
			*        } //by k
			*      } //by j
			*    } //by i
			*    //Smogarinsky (End)
			*
			*
			*    delete &uTr;
			*    delete &vTr;
			*    delete &wTr;
			*
			*
			*    end = clock();
			*    diff = (double)(end - start);
			*    diff = diff/CLOCKS_PER_SEC;
			*    printf("Turbulence model: OK. (time = %lf s.)\n",diff);
			*    }
			*/
			//Viscosity (End)





			//Saving results (Begin)
			if (timeStepNumber % TProblem::GlobalSaveStep == 0 || timeStepNumber == 1)
			{
				printf("\n");
				printf("Saving results: Starting...\n");
				start = clock();

				TRnRelease3dSpace& uR = u;
				TRnRelease3dSpace& vR = v;
				TRnRelease3dSpace& wR = w;

				TRnRelease3dSpace& rpTrXR = rpTrX;
				TRnRelease3dSpace& rpTrYR = rpTrY;
				TRnRelease3dSpace& rpTrZR = rpTrZ;
				
				CountBoundaryP(p);
				CountBoundaryConcentration(cnc);
				CountBoundaryViscosity(vis);
				//ComputeStress(u,v,w,timeStepNumber);
				//---------------------------------------------------------------
				// Output to VTK
				//---------------------------------------------------------------


				int len = strlen(TProblem::GlobalCatalog);
				char *vtk_fileName = new char[len + 50];
				string_copy(vtk_fileName, TProblem::GlobalCatalog, len + 50);

				char fileName[40];
				string_print(fileName, 40, "flow%d.vtk", 10000 + timeStepNumber);
				string_concat(vtk_fileName, fileName, len + 50);


				FILE* vtk_ftp = file_open(vtk_fileName, "w");

				fprintf(vtk_ftp, "# vtk DataFile Version 3.0\n");
				fprintf(vtk_ftp, "Flow\nASCII\n\n");
				fprintf(vtk_ftp, "DATASET STRUCTURED_GRID\n");

				T3dNormalGrid& vtk_grid = dynamic_cast<T3dNormalGrid&>(*fGrid);
				int vtk_N = vtk_grid.GetSeparator1().EndIndex;
				int vtk_L = vtk_grid.GetSeparator2().EndIndex;
				int vtk_M = vtk_grid.GetSeparator3().EndIndex;
				int vtk_NML = (vtk_N )*(vtk_L )*(vtk_M );
				fprintf(vtk_ftp, "DIMENSIONS %d %d %d\n", vtk_N , vtk_L , vtk_M );
				fprintf(vtk_ftp, "POINTS %d double\n", vtk_NML);
				const T3dNumberMask& vtk_mask = vtk_grid.Mask;

				const double* vtk_x = vtk_grid.GetSeparator1().Dimension;
				const double* vtk_y = vtk_grid.GetSeparator2().Dimension;
				const double* vtk_z = vtk_grid.GetSeparator3().Dimension;
				//const T3dNumberMask& vtk_mask = grid.Mask;
				for (int k = 1; k < vtk_M + 1; k++)
					for (int j = 1; j < vtk_L + 1; j++)
						for (int i = 1; i < vtk_N + 1; i++)
							fprintf(vtk_ftp, "%lf  %lf  %lf\n", (float)vtk_x[i] * (TProblem::GlobalLength), (float)vtk_y[j] * (TProblem::GlobalLength), (float)vtk_z[k] * (TProblem::GlobalLength));

				T3dNormalGrid& vtk_gridP = dynamic_cast<T3dNormalGrid&>(*fGridP);
				const T3dNumberMask& vtk_maskP = vtk_gridP.Mask;
				fprintf(vtk_ftp, "POINT_DATA %d\n\nVECTORS Velocity double\n", vtk_NML);
				for (int k = 1; k < vtk_M + 1; k++)
				{
					for (int j = 1; j < vtk_L + 1; j++)
					{
						for (int i = 1; i < vtk_N + 1; i++)
						{

							double uu;
							double vv;
							double ww;

							int mijk = vtk_mask[i][j][k];
							if (mijk == TFictivePoint) { uu = 0; vv = 0; ww = 0; }
							else
							{
								uu = (uR[i][j][k] + uR[i][j + 1][k] + uR[i][j][k + 1] + uR[i][j + 1][k + 1]) / 4;
								vv = (vR[i][j][k] + vR[i + 1][j][k] + vR[i][j][k + 1] + vR[i + 1][j][k + 1]) / 4;
								ww = (wR[i][j][k] + wR[i + 1][j][k] + wR[i][j + 1][k] + wR[i + 1][j + 1][k]) / 4;
							}

							//fprintf(f,"%lf,",dnst);
							fprintf(vtk_ftp, "%lf ", (double)uu*(TProblem::GlobalVelocity));
							fprintf(vtk_ftp, "%lf ", (double)vv*(TProblem::GlobalVelocity));
							fprintf(vtk_ftp, "%lf\n ", (double)ww*(TProblem::GlobalVelocity));

						}
					}
				}

				fprintf(vtk_ftp, "\nSCALARS Pressure float\nLOOKUP_TABLE default\n");
				for (int k = 1; k < vtk_M + 1; k++)
				{
					for (int j = 1; j < vtk_L + 1; j++)
					{
						for (int i = 1; i < vtk_N + 1; i++)
						{
							double pp;
							int mijk = vtk_maskP[i][j][k];
							if (mijk == TFictivePoint)
							{
								pp = 0;
							}
							else
							{
								bool avr1 = ((vtk_maskP[i + 1][j][k]) != TFictivePoint);
								bool avr2 = ((vtk_maskP[i + 1][j + 1][k]) != TFictivePoint);
								bool avr3 = ((vtk_maskP[i + 1][j][k + 1]) != TFictivePoint);
								bool avr4 = ((vtk_maskP[i + 1][j + 1][k + 1]) != TFictivePoint);
								bool avr5 = ((vtk_maskP[i][j][k]) != TFictivePoint);
								bool avr6 = ((vtk_maskP[i][j + 1][k]) != TFictivePoint);
								bool avr7 = ((vtk_maskP[i][j][k + 1]) != TFictivePoint);
								bool avr8 = ((vtk_maskP[i][j + 1][k + 1]) != TFictivePoint);
								int avr = (int)avr1 + (int)avr2 + (int)avr3 + (int)avr4 + (int)avr5 + (int)avr6 + (int)avr7 + (int)avr8;

								if (avr == 0)
								{
									throw "Error in TVPProblem::SolveBySplitting: Grid is incorrect while saving pressure results.";
								}

								double p1 = p[i + 1][j][k];
								double p2 = p[i + 1][j + 1][k];
								double p3 = p[i + 1][j][k + 1];
								double p4 = p[i + 1][j + 1][k + 1];
								double p5 = p[i][j][k];
								double p6 = p[i][j + 1][k];
								double p7 = p[i][j][k + 1];
								double p8 = p[i][j + 1][k + 1];

								pp = (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8) / avr;
							}
							fprintf(vtk_ftp, "%lf\n ", (double)pp);
						}
					}
				}

				fprintf(vtk_ftp, "\nSCALARS Density float\nLOOKUP_TABLE default\n");
				const T3dNumberMask& vtk_maskD = gridDensity.Mask;
				for (int k = 1; k < vtk_M + 1; k++)
				{
					for (int j = 1; j < vtk_L + 1; j++)
					{
						for (int i = 1; i < vtk_N + 1; i++)
						{
							double dnst;
							int mijk = vtk_maskD[i][j][k];
							if (mijk == TFictivePoint) { dnst = 0; }
							else
							{
								bool avr1 = ((vtk_maskD[i + 1][j][k]) != TFictivePoint);
								bool avr2 = ((vtk_maskD[i + 1][j + 1][k]) != TFictivePoint);
								bool avr3 = ((vtk_maskD[i + 1][j][k + 1]) != TFictivePoint);
								bool avr4 = ((vtk_maskD[i + 1][j + 1][k + 1]) != TFictivePoint);
								bool avr5 = ((vtk_maskD[i][j][k]) != TFictivePoint);
								bool avr6 = ((vtk_maskD[i][j + 1][k]) != TFictivePoint);
								bool avr7 = ((vtk_maskD[i][j][k + 1]) != TFictivePoint);
								bool avr8 = ((vtk_maskD[i][j + 1][k + 1]) != TFictivePoint);
								int avr = (int)avr1 + (int)avr2 + (int)avr3 + (int)avr4 + (int)avr5 + (int)avr6 + (int)avr7 + (int)avr8;

								if (avr == 0)
								{
									throw "Error in TVPProblem::SolveBySplitting: Grid is incorrect while saving concentration results.";
								}

								double d1 = density[i + 1][j][k];
								double d2 = density[i + 1][j + 1][k];
								double d3 = density[i + 1][j][k + 1];
								double d4 = density[i + 1][j + 1][k + 1];
								double d5 = density[i][j][k];
								double d6 = density[i][j + 1][k];
								double d7 = density[i][j][k + 1];
								double d8 = density[i][j + 1][k + 1];

								dnst = (d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8) / avr;
							}
							fprintf(vtk_ftp, "%lf\n", (double)dnst);

						}
					}
				}

				fprintf(vtk_ftp, "\nSCALARS Concentration float\nLOOKUP_TABLE default\n");
				const T3dNumberMask& vtk_maskC = gridC.Mask;
				//	  fprintf(ftp,"ZONE I=%d, J=%d, K=%d, F=POINT\n", N+1, L+1, M+1);
				for (int k = 1; k < vtk_M + 1; k++)
				{
					for (int j = 1; j < vtk_L + 1; j++)
					{
						for (int i = 1; i < vtk_N + 1; i++)
						{
							double cc;
							int mijk = vtk_maskC[i][j][k];
							if (mijk == TFictivePoint) { cc = 0; }
							else
							{
								bool avr1 = ((vtk_maskC[i + 1][j][k]) != TFictivePoint);
								bool avr2 = ((vtk_maskC[i + 1][j + 1][k]) != TFictivePoint);
								bool avr3 = ((vtk_maskC[i + 1][j][k + 1]) != TFictivePoint);
								bool avr4 = ((vtk_maskC[i + 1][j + 1][k + 1]) != TFictivePoint);
								bool avr5 = ((vtk_maskC[i][j][k]) != TFictivePoint);
								bool avr6 = ((vtk_maskC[i][j + 1][k]) != TFictivePoint);
								bool avr7 = ((vtk_maskC[i][j][k + 1]) != TFictivePoint);
								bool avr8 = ((vtk_maskC[i][j + 1][k + 1]) != TFictivePoint);
								int avr = (int)avr1 + (int)avr2 + (int)avr3 + (int)avr4 + (int)avr5 + (int)avr6 + (int)avr7 + (int)avr8;

								if (avr == 0)
								{
									throw "Error in TVPProblem::SolveBySplitting: Grid is incorrect while saving concentration results.";
								}

								double c1 = cnc[i + 1][j][k];
								double c2 = cnc[i + 1][j + 1][k];
								double c3 = cnc[i + 1][j][k + 1];
								double c4 = cnc[i + 1][j + 1][k + 1];
								double c5 = cnc[i][j][k];
								double c6 = cnc[i][j + 1][k];
								double c7 = cnc[i][j][k + 1];
								double c8 = cnc[i][j + 1][k + 1];

								cc = (c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8) / avr;
							}

							fprintf(vtk_ftp, "%lf\n", (double)cc);
						}
					}
				}
				fprintf(vtk_ftp, "\nSCALARS Viscosity float\nLOOKUP_TABLE default\n");
				//const T3dNumberMask& vtk_maskC = gridC.Mask;
				//	  fprintf(ftp,"ZONE I=%d, J=%d, K=%d, F=POINT\n", N+1, L+1, M+1);
				for (int k = 1; k < vtk_M + 1; k++)
				{
					for (int j = 1; j < vtk_L + 1; j++)
					{
						for (int i = 1; i < vtk_N + 1; i++)
						{
							double cc;
							int mijk = vtk_maskC[i][j][k];
							if (mijk == TFictivePoint) { cc = 0; }
							else
							{
								bool avr1 = ((vtk_maskC[i + 1][j][k]) != TFictivePoint);
								bool avr2 = ((vtk_maskC[i + 1][j + 1][k]) != TFictivePoint);
								bool avr3 = ((vtk_maskC[i + 1][j][k + 1]) != TFictivePoint);
								bool avr4 = ((vtk_maskC[i + 1][j + 1][k + 1]) != TFictivePoint);
								bool avr5 = ((vtk_maskC[i][j][k]) != TFictivePoint);
								bool avr6 = ((vtk_maskC[i][j + 1][k]) != TFictivePoint);
								bool avr7 = ((vtk_maskC[i][j][k + 1]) != TFictivePoint);
								bool avr8 = ((vtk_maskC[i][j + 1][k + 1]) != TFictivePoint);
								int avr = (int)avr1 + (int)avr2 + (int)avr3 + (int)avr4 + (int)avr5 + (int)avr6 + (int)avr7 + (int)avr8;

								if (avr == 0)
								{
									throw "Error in TVPProblem::SolveBySplitting: Grid is incorrect while saving concentration results.";
								}

								double c1 = vis[i + 1][j][k];
								double c2 = vis[i + 1][j + 1][k];
								double c3 = vis[i + 1][j][k + 1];
								double c4 = vis[i + 1][j + 1][k + 1];
								double c5 = vis[i][j][k];
								double c6 = vis[i][j + 1][k];
								double c7 = vis[i][j][k + 1];
								double c8 = vis[i][j + 1][k + 1];

								cc = (c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8) / avr;
							}

							fprintf(vtk_ftp, "%lf\n", (double)cc);
						}
					}
				}
				fclose(vtk_ftp);
				delete[]vtk_fileName;// vtk_fileName = NULL;


				//---------------------------------------------------------------
				// Output to .dat
				//---------------------------------------------------------------

				/*

				len = strlen(TProblem::GlobalCatalog);
				char *zonesName = new char[len + 50];
				string_copy(zonesName, TProblem::GlobalCatalog, len + 50);

				char fileName1[40];
				string_print(fileName1, 40, "full%d.dat", 10000 + timeStepNumber);
				string_concat(zonesName, fileName1, len + 50);

				FILE* f_full = file_open(zonesName, "w");
				T3dNormalGrid& out_grid = dynamic_cast<T3dNormalGrid&>(*fGrid);
				int out_x = out_grid.GetSeparator1().EndIndex;
				int out_y = out_grid.GetSeparator2().EndIndex;
				int out_z = out_grid.GetSeparator3().EndIndex;

				const double* d_x = vtk_grid.GetSeparator1().Dimension;
				const double* d_y = vtk_grid.GetSeparator2().Dimension;
				const double* d_z = vtk_grid.GetSeparator3().Dimension;

				fprintf(f_full, "TITLE = \"try\" \n");
				fprintf(f_full, "VARIABLES = \"X\",\"Y\",\"Z\",\"Pressure\",\"U\",\"V\",\"W\",\"Concentration\",\"Density\",\"Viscosity\"\n");
				fprintf(f_full, "ZONE T=\"Test\", I=%d, J=%d, K=%d, F=POINT\n", out_x, out_y, out_z);

				for (int k = 1; k <= out_z; k++)
				{
					for (int j = 1; j <= out_y; j++)
					{
						for (int i = 1; i <= out_x; i++)
						{
							fprintf(f_full, "%lf %lf %lf ", (float)d_x[i] * (TProblem::GlobalLength), (float)d_y[j] * (TProblem::GlobalLength), (float)d_z[k] * (TProblem::GlobalLength));

							double pp;
							int mijk_p = vtk_maskP[i][j][k];
							if (mijk_p == TFictivePoint)
							{
								pp = 0;
							}
							else
							{
								bool avr1 = ((vtk_maskP[i + 1][j][k]) != TFictivePoint);
								bool avr2 = ((vtk_maskP[i + 1][j + 1][k]) != TFictivePoint);
								bool avr3 = ((vtk_maskP[i + 1][j][k + 1]) != TFictivePoint);
								bool avr4 = ((vtk_maskP[i + 1][j + 1][k + 1]) != TFictivePoint);
								bool avr5 = ((vtk_maskP[i][j][k]) != TFictivePoint);
								bool avr6 = ((vtk_maskP[i][j + 1][k]) != TFictivePoint);
								bool avr7 = ((vtk_maskP[i][j][k + 1]) != TFictivePoint);
								bool avr8 = ((vtk_maskP[i][j + 1][k + 1]) != TFictivePoint);
								int avr = (int)avr1 + (int)avr2 + (int)avr3 + (int)avr4 + (int)avr5 + (int)avr6 + (int)avr7 + (int)avr8;

								if (avr == 0)
								{
									throw "Error in TVPProblem::SolveBySplitting: Grid is incorrect while saving pressure results.";
								}

								double p1 = p[i + 1][j][k];
								double p2 = p[i + 1][j + 1][k];
								double p3 = p[i + 1][j][k + 1];
								double p4 = p[i + 1][j + 1][k + 1];
								double p5 = p[i][j][k];
								double p6 = p[i][j + 1][k];
								double p7 = p[i][j][k + 1];
								double p8 = p[i][j + 1][k + 1];

								pp = (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8) / avr;
							}
							fprintf(f_full, "%lf ", (double)pp);

							double uu;
							double vv;
							double ww;

							int mijk_uvw = vtk_mask[i][j][k];
							if (mijk_uvw == TFictivePoint) { uu = 0; vv = 0; ww = 0; }
							else
							{
								uu = (uR[i][j][k] + uR[i][j + 1][k] + uR[i][j][k + 1] + uR[i][j + 1][k + 1]) / 4;
								vv = (vR[i][j][k] + vR[i + 1][j][k] + vR[i][j][k + 1] + vR[i + 1][j][k + 1]) / 4;
								ww = (wR[i][j][k] + wR[i + 1][j][k] + wR[i][j + 1][k] + wR[i + 1][j + 1][k]) / 4;
							}

							//fprintf(f,"%lf,",dnst);
							fprintf(f_full, "%lf ", (double)uu*(TProblem::GlobalVelocity));
							fprintf(f_full, "%lf ", (double)vv*(TProblem::GlobalVelocity));
							fprintf(f_full, "%lf ", (double)ww*(TProblem::GlobalVelocity));

							double cc;
							double dnst;
							double visc;
							int mijk_conc = vtk_maskC[i][j][k];
							if (mijk_conc == TFictivePoint) { cc = 0; dnst = 0; visc = 0;}
							else
							{
								bool avr1 = ((vtk_maskC[i + 1][j][k]) != TFictivePoint);
								bool avr2 = ((vtk_maskC[i + 1][j + 1][k]) != TFictivePoint);
								bool avr3 = ((vtk_maskC[i + 1][j][k + 1]) != TFictivePoint);
								bool avr4 = ((vtk_maskC[i + 1][j + 1][k + 1]) != TFictivePoint);
								bool avr5 = ((vtk_maskC[i][j][k]) != TFictivePoint);
								bool avr6 = ((vtk_maskC[i][j + 1][k]) != TFictivePoint);
								bool avr7 = ((vtk_maskC[i][j][k + 1]) != TFictivePoint);
								bool avr8 = ((vtk_maskC[i][j + 1][k + 1]) != TFictivePoint);
								int avr = (int)avr1 + (int)avr2 + (int)avr3 + (int)avr4 + (int)avr5 + (int)avr6 + (int)avr7 + (int)avr8;

								if (avr == 0)
								{
									throw "Error in TVPProblem::SolveBySplitting: Grid is incorrect while saving concentration results.";
								}

								double c1 = cnc[i + 1][j][k];
								double c2 = cnc[i + 1][j + 1][k];
								double c3 = cnc[i + 1][j][k + 1];
								double c4 = cnc[i + 1][j + 1][k + 1];
								double c5 = cnc[i][j][k];
								double c6 = cnc[i][j + 1][k];
								double c7 = cnc[i][j][k + 1];
								double c8 = cnc[i][j + 1][k + 1];

								cc = (c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8) / avr;

								double d1 = density[i + 1][j][k];
								double d2 = density[i + 1][j + 1][k];
								double d3 = density[i + 1][j][k + 1];
								double d4 = density[i + 1][j + 1][k + 1];
								double d5 = density[i][j][k];
								double d6 = density[i][j + 1][k];
								double d7 = density[i][j][k + 1];
								double d8 = density[i][j + 1][k + 1];

								dnst = (d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8) / avr;
								
								double v1 = vis[i + 1][j][k];
								double v2 = vis[i + 1][j + 1][k];
								double v3 = vis[i + 1][j][k + 1];
								double v4 = vis[i + 1][j + 1][k + 1];
								double v5 = vis[i][j][k];
								double v6 = vis[i][j + 1][k];
								double v7 = vis[i][j][k + 1];
								double v8 = vis[i][j + 1][k + 1];

								visc = (v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8) / avr;
							}
							fprintf(f_full, "%lf ", (double)cc);
							fprintf(f_full, "%lf ", (double)dnst);
							fprintf(f_full, "%lf\n", (double)visc);
						}
					}
				}

				fclose(f_full);
				delete[]zonesName;// zonesName = NULL;
				//*/
				

				//------------------------------------------------------------------
				// Two format output end
				//------------------------------------------------------------------


				//CountBoundaryP(p);
			#if CNPY_VISUALIZATION
				const unsigned int shape[] = { M + 1,L + 1,N + 1 };

				string_print(fileName, 40, "u_velocity_%03d.npy", timeStepNumber);
				cnpy::npy_save(fileName, U, shape, 3, "w");

				string_print(fileName, 40, "v_velocity_%03d.npy", timeStepNumber);
				cnpy::npy_save(fileName, V, shape, 3, "w");

				string_print(fileName, 40, "w_velocity_%03d.npy", timeStepNumber);
				cnpy::npy_save(fileName, W, shape, 3, "w");

				string_print(fileName, 40, "fx_%03d.npy", timeStepNumber);
				cnpy::npy_save(fileName, FX, shape, 3, "w");

				string_print(fileName, 40, "fy_%03d.npy", timeStepNumber);
				cnpy::npy_save(fileName, FY, shape, 3, "w");

				string_print(fileName, 40, "fz_%03d.npy", timeStepNumber);
				cnpy::npy_save(fileName, FZ, shape, 3, "w");


				string_print(fileName, 40, "pressure_%03d.npy", timeStepNumber);
				cnpy::npy_save(fileName, P, shape, 3, "w");

				string_print(fileName, 40, "concentration_%03d.npy", timeStepNumber);
				cnpy::npy_save(fileName, CNC, shape, 3, "w");

				//snprintf(fileName, 40, "viscosity_%03d.npy", timeStepNumber);
				//cnpy::npy_save(fileName, VIS, shape, 3, "w");

				string_print(fileName, 40, "density_%03d.npy", timeStepNumber);
				cnpy::npy_save(fileName, DNST, shape, 3, "w");

				const unsigned int shape_x[] = { N + 1 };
				const unsigned int shape_y[] = { L + 1 };
				const unsigned int shape_z[] = { M + 1 };

				string_print(fileName, 40, "x_coord.npy", timeStepNumber);
				cnpy::npy_save(fileName, x, shape_x, 1, "w");

				string_print(fileName, 40, "y_coord.npy", timeStepNumber);
				cnpy::npy_save(fileName, y, shape_y, 1, "w");

				string_print(fileName, 40, "z_coord.npy", timeStepNumber);
				cnpy::npy_save(fileName, z, shape_z, 1, "w");

				delete U;
				delete V;
				delete W;
				delete P;

				delete FX;
				delete FY;
				delete FZ;
			#endif


				OutputBoundary(timeStepNumber);
				//ComputeChanging(timeStepNumber);
				//system("./replace.sh");

				delete &rpTrXR;
				delete &rpTrYR;
				delete &rpTrZR;


				end = clock();
				diff = (double)(end - start);
				diff = diff / CLOCKS_PER_SEC;
				printf("Saving rerults: OK. (time = %lf s.)\n", diff);
			}
			//Saving results (End)


			//Reassigning (Begin)
			{
				printf("\n");
				printf("Reassigning: Starting...\n");
				start = clock();

				u_.Assign(u);
				v_.Assign(v);
				w_.Assign(w);

				u_1.Assign(u);
				v_1.Assign(v);
				w_1.Assign(w);

				//Concentration (Begin)
				cnc_1.Assign(cnc);
				//Concentration (End)

				end = clock();
				diff = (double)(end - start);
				diff = diff / CLOCKS_PER_SEC;
				printf("Reassigning: OK. (time = %lf s.)\n", diff);
			}
			//Reassigning (End)


			flagRestruct = RestructStatement(fStatement, fU, fV, fW, fP);
		
		printf("\n");
		printf("Time step = %d: OK.\n", timeStepNumber);
		printf("*****************************************************\n");
	} // by time


	  //Deleting objects (Begin)
	printf("\n\n");
	printf("Deleting temporary objects: Starting...\n");

	delete u_1Ptr;
	delete v_1Ptr;
	delete w_1Ptr;

	delete u_Ptr;
	delete v_Ptr;
	delete w_Ptr;

	delete uSPtr;
	delete vSPtr;
	delete wSPtr;

	delete slae;


	delete rightPartPtr;


	delete APtr;
	if (ARnPtr != NULL) { delete ARnPtr; }


	//Concentration (Begin)
	delete C_1Ptr;
	delete CPtr;
	delete gridCPtr;
	//Concentration (End)

	//Viscosity (Begin)
	delete gridVisPtr;
	delete visPtr;
	//Viscosity (End)


	printf("Deleting temporary objects: OK.\n");
	//Deleting objects (End)

	printf("\n\n");
	printf("Splitting method: OK.\n");
}
//*************************************


//#################################################################



//*************************************
bool TScourProblem::RestructStatement(TMaskGenerator* lStatement)
{
	//Initialization (Begin)
	TMaskGenerator& state = *lStatement;

	TRnRelease3dSpace &U = *fU;
    TRnRelease3dSpace &V = *fV;
    TRnRelease3dSpace &W = *fW;
    TRnRelease3dSpace &P = *fP;
	//Initialization (End)


	//Restructing (Begin)

	// Постоянные
	const double
		d = 2.5e-3,		// диаметр частиц
		mu_d = 0.51,	// коэффициент динамического трения
		mu_s = 0.63,	// коэффициент статического трения
		Y_cr0 = 0.05,	// критическое число Шильдса для горизонтального дна
		g = 9.81,		// ускорение свободного падения
		rho = 998.2,	// плотность жидкости
		rho_s = 2655,	// плотность донного материала
		s = rho_s/rho,	// удельный вес донных частиц
		A = 10,			// параметр зависимости придонной скорости потока от скорости трения
		por = 0.4,		// коэффициент пористости
		k = 0.4,		// константа Кармана
		k_s = 2.5*d,	// высота слоя шероховатости
		z_0 = k_s/30,	// шероховатость дна
		C_D = 0.27;		// Коэффициент пробуксовки ( drag(force) coeff )

	// Массивы
	//glTVector
	//;	**q_b;
	glAllocPlane(q_b,Nx,Ny);

	//int
	//	**gorisont;
	glAllocPlane(gorisont, Nx, Ny);
	glAllocPlane(P_ef_print, Nx, Ny);
	glAllocPlane(U_b_print, Nx, Ny);
	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
		{
			gorisont[i][j] = 0;
			P_ef_print[i][j] = 0;
			U_b_print[i][j].set(0,0,0);
		}

	for(int tt=0; tt<500; ++tt)
	{

	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			q_b[i][j].set(0,0,0);

	// Расчет транспорта через каждый узел
	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
		{
			if(Gp[i][j][ z_indx[i][j]+2 ] != 1)
				continue;

			glTVector
				vel(
					(U[i][j][ z_indx[i][j]+2 ] + U[i-1][j][ z_indx[i][j]+2 ])/2,
					(V[i][j][ z_indx[i][j]+2 ] + V[i][j-1][ z_indx[i][j]+2 ])/2,
					W[i][j][ z_indx[i][j]+2 ]
				);

			double
				u = len( vel ), // скорость в ближайшем к стенке узле
				delta = Hz[ z_indx[i][j]+2 ], // расстояние от стенки до ближайшего узла
				V_fs = k*u/log(delta/z_0), // скорость трения (friction velocity)
				Y = V_fs*V_fs/g/(s-1)/d; // число Шильдса

			// Наибольший уклон (steepest bed slope)
			double
				dx, dy;

			if(Gp[i+1][j][ z_indx[i][j]+2 ] != 1)
				dx = ( h[i][j]-h[i-1][j] )/Hx[i];
			else if(Gp[i-1][j][ z_indx[i][j]+2 ] != 1)
				dx = ( h[i+1][j]-h[i][j] )/Hx[i+1];
			else
				dx = ( h[i+1][j]-h[i-1][j] )/(Hx[i+1]+Hx[i]);

			if(Gp[i][j+1][ z_indx[i][j]+2 ] != 1)
				dy = ( h[i][j]-h[i][j-1] )/Hy[j];
			else if(Gp[i][j-1][ z_indx[i][j]+2 ] != 1)
				dy = ( h[i][j+1]-h[i][j] )/Hy[j+1];
			else
				dy = ( h[i][j+1]-h[i][j-1] )/(Hy[j+1]+Hy[j]);

			glTVector
				n( -dx, -dy, 1 ); // нормаль к поверхности F=0: n = grad F, F = z-h(x,y) <= z=h(x,y)
			normalize(n);

			glTVector
				ortk(0,0,1), // единичный вектор оси Oz
				line, // направляющий вектор прямой пересечения касательной плоскости к поверхности дна с горизонталью
				ksi; // вектор максимального уклона
			cross(n,ortk,line);
			cross(n,line,ksi);
			normalize(ksi);

			// Проекция скорости на касательную плоскость
			/*
vec3 N = cross(OA, OB) - получили нормаль к плоскости.
Знаем, что уравнение плоскости имеет вид A*x + B*y + C*z + D = 0
A, B, C - Это как раз компоненты нормали. Т.е. N = {A, B, C}.
Находим D, подставляя в уравнение координаты точки, принадлежащей плоскости, например точки О:
  D = -(A*O.x + B*O.y + C*O.z)

Теперь ищем проекцию точки Х, назовем Р. Эта точка лежит в плоскости, т.е. должно выполнятся равенство A*Р.x + B*Р.y + C*Р.z + D = 0,
с другой стороны, эта точка лежит на прямой, которая образуется точкой Х и выпущенным из неё вектором нормали N(по определению проекции),
т.е. P = X + K*N, где К - некоторый коэффициент. Подставляем это в уравнение прямой:
  A*(X.x + K*N.x) + B*(X.y + K*N.y) + C*(X.z + K*N.z) + D = 0
  A*X.x + B*X.y + C*X.z + K*(A*N.x + B*N.y + C*N.z) + D = 0

  Не забываем, что  N = {A, B, C}:

  dot(N, X) + K*dot(N, N) + D = 0
  K = -(dot(N, X) + D)/dot(N, N)

  И раз нашли K, то находим точку проекции Р:
  P = X + K*N
			*/
			glTVector
				EndV( vel.x + Cx[i], vel.y + Cy[j], vel.z + h[i][j]); // конец вектора скорости

			double
				D = -dot(n,glTVector(Cx[i],Cy[j],h[i][j])), // коэффициент D касательной плоскости
				K = -dot(n,EndV)-D; // коэффициент масштабирования при проецировании: PV = K*n
			glTVector
				Vprj( vel.x + K*n.x, vel.y + K*n.y, vel.z + K*n.z);

			if( len(Vprj)<1e-10 )
			{
				q_b[i][j].set(0,0,0);
				continue;
			}

			// Доля движущихся частиц и придонная скорость частицы
			double
				P_ef;
			glTVector
				U_b;

			if( len(line) < 1e-10) // дно горизонтальное
			{
				gorisont[i][j] = 0;
				P_ef = pow( 1 + pow(1.0/6*M_PI*mu_d/(Y-Y_cr0),4), -0.25);
				mult(Vprj,0.00001*A*V_fs/len(Vprj),U_b);
			}
			else // дно под наклоном
			{
				gorisont[i][j] = 1;
				// Параметры-углы
				double
					t = dot(Vprj,ksi)/len(Vprj)/len(ksi);
				if(t>1)
					t = 1;
				else if(t<-1)
					t = -1;
				double
					alpha = acos(t); // угол между скоростью потока и наибольшим уклоном

				double
					p = len(n);
				t = dot(n,ortk)/p;
				if(t>1)
					t = 1;
				else if(t<-1)
					t = -1;
				double
					beta = acos( t ); // угол между касательной плоскостью и горизонталью

				// Доля движущихся частиц
				double
					beta1 = beta;
				if(beta1>0.55)
				{
					beta1 = 0.55;
				}
				double
					Y_cr = Y_cr0 * ( cos(beta1)*sqrt(1-sin(alpha)*sin(alpha)*tan(beta1)*tan(beta1)/mu_s/mu_s) - cos(alpha)*sin(beta1)/mu_s ); // критическое число Шильдса
				P_ef = pow( 1 + pow(1.0/6*M_PI*mu_d/(Y-Y_cr),4), -0.25);

				// Силы, действующие на частицу
				double
					W = 1.0/6*M_PI*rho*g*(s-1)*d*d*d, // сила тяжести
					F_D = 0.5*rho*C_D*M_PI_4*d*d; // сила сдвига

				double
					U_r = 1,
					psi1 = 0.1,
					psi = 0.1;

				SolveSys(F_D, W*sin(beta), alpha, mu_d*W*cos(beta), A*V_fs, U_r, psi1, psi);

				double
					U_b_len = A*V_fs*cos(psi) - U_r*cos(psi1), // придонная скорость частицы (длина вектора)
					th = alpha - psi; // угол между вектором придонной скорости и вектором наибольшего уклона в касательной плоскости

				glTVector
					os;
				cross(ksi, Vprj, os);
				if(os.z<0)
					th = -th;

				// матрица поворота вокруг n на угол th
				glTMatrix RotM (
					cos(th) + (1-cos(th))*n.x*n.x,       (1-cos(th))*n.x*n.y - sin(th)*n.z,  (1-cos(th))*n.x*n.z + sin(th)*n.y,
					(1-cos(th))*n.y*n.x + sin(th)*n.z,  cos(th) + (1-cos(th))*n.y*n.y,       (1-cos(th))*n.y*n.z - sin(th)*n.x,
					(1-cos(th))*n.z*n.x - sin(th)*n.y,  (1-cos(th))*n.z*n.y + sin(th)*n.x,  cos(th) + (1-cos(th))*n.z*n.z         );

				// вектор придонной скорости
				glTVector
					U_b_t;
				mult(RotM, ksi, U_b_t);
				mult(U_b_t, U_b_len, U_b);
			}

			// Расчет вектора транспорта q_b
			double
				t = 1.0/6*M_PI*d*P_ef;
			mult(U_b, t, q_b[i][j]);
			P_ef_print[i][j] = P_ef;
			U_b_print[i][j] = U_b;
		}

	// Пересчет h

	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
		{
			if(Gp[i][j][ z_indx[i][j]+2 ] != 1)
				continue;
			double
				dx, dy;

			if(Gp[i+1][j][ z_indx[i][j]+2 ] != 1)
				dx = ( q_b[i][j].x-q_b[i-1][j].x )/Hx[i];
			else if(Gp[i-1][j][ z_indx[i][j]+2 ] != 1)
				dx = ( q_b[i+1][j].x-q_b[i][j].x )/Hx[i+1];
			else
				dx = (q_b[i+1][j].x-q_b[i-1][j].x)/(Hx[i+1]+Hx[i]);

			if(Gp[i][j+1][ z_indx[i][j]+2 ] != 1)
				dy = ( q_b[i][j].y-q_b[i][j-1].y )/Hy[j];
			else if(Gp[i][j-1][ z_indx[i][j]+2 ] != 1)
				dy = ( q_b[i][j+1].y-q_b[i][j].y )/Hy[j+1];
			else
				dy = (q_b[i][j+1].y-q_b[i][j-1].y)/(Hy[j+1]+Hy[j]);

			//h[i][j] = h[i][j] - tau/(1-por)*( dx + dy);
			h[i][j] = h[i][j] - 0.001*( dx + dy);
		}

	if(tt%10 == 1)
	{
		PrintHeight(tt/10+1);
		//printf("                              scour iter #%d\n", tt);
	}
	//PrintHeight();
	}

	//Restructing (End)
	glFreePlane(q_b,Nx,Ny);


	return false;
}
//*************************************
bool TScourProblem::SolveSys(double a, double b, double c, double d, double k, double &x1, double &x2, double &x3)
{
	bool notgoodsolve = true;
	x1 = 1;
	for(double t2=0.1; (t2<M_PI)&&notgoodsolve; t2 += M_PI/6)
		for(double t3=0.1; (t3<M_PI)&&notgoodsolve; t3 += M_PI/6)
		{
			x2 = t2;
			x3 = t3;
			glTVector z;
			int
				ite = 0;
			do
			{
				++ite;
				glTVector F( a*x1*x1*cos(x2) + b*cos(c-x3) - d, a*x1*x1*sin(x2) + b*sin(c-x3), x1*sin(x2) - k*sin(x3) );
				glTMatrix J(
					2*a*cos(x2)*x1, -a*x1*x1*sin(x2),  b*sin(c-x3),
					2*a*sin(x2)*x1,  a*x1*x1*cos(x2), -b*cos(c-x3),
					sin(x2),         x1*cos(x2),      -k*cos(x3)      );
				glTMatrix invJ;
				glTVector s, x(x1,x2,x3), y;
				invers(J,invJ);
				mult(invJ,F,s);
				minus(x,s,y);
				minus(x,y,z);
				x1 = y.x;
				x2 = y.y;
				x3 = y.z;
			} while( (len(z)>1e-10)&&(ite<1000) );
			if(ite<1000)
				notgoodsolve = false;
		}

	return !notgoodsolve;
}
//*************************************
void TScourProblem::PrintHeight(int t)
{
	{
		int
			len = strlen(TProblem::GlobalCatalog);
		char *scourName = new char [len+50];
        string_copy(scourName, TProblem::GlobalCatalog, len+50);

		char fileName[40];
        string_print(fileName, 40, "\\Scour\\ParaView\\Scour%d.vtk", t);
        string_concat(scourName, fileName, len+50);

        FILE *fsc_par = file_open(scourName, "w");
		delete[] scourName;// scourName = NULL;

		fprintf(fsc_par,"# vtk DataFile Version 3.0\n");
		fprintf(fsc_par,"Scour data\n");
		fprintf(fsc_par,"ASCII\n");
		fprintf(fsc_par,"DATASET RECTILINEAR_GRID\n");
		fprintf(fsc_par,"DIMENSIONS %d %d %d\n", Nx, Ny, 1);

		fprintf(fsc_par,"X_COORDINATES %d double\n", Nx);
		for(int i=0; i<Nx; ++i)
			fprintf(fsc_par,"%lf ", (Cx[i])*(TProblem::GlobalLength));
		fprintf(fsc_par,"\n");

		fprintf(fsc_par,"Y_COORDINATES %d double\n", Ny);
		for(int i=0; i<Ny; ++i)
			fprintf(fsc_par,"%lf ", (Cy[i])*(TProblem::GlobalLength));
		fprintf(fsc_par,"\n");

		fprintf(fsc_par,"Z_COORDINATES %d double\n", 1);
		fprintf(fsc_par,"0\n");

		fprintf(fsc_par,"POINT_DATA %d\n", Nx*Ny);
		fprintf(fsc_par,"SCALARS ScourData double\n");
		fprintf(fsc_par,"LOOKUP_TABLE default\n");

		for(int j=0; j<Ny; ++j)
		{
			for(int i=0; i<Nx; ++i)
			{
				fprintf(fsc_par, "%lf ", (h[i][j])*(TProblem::GlobalLength));
			}
			fprintf(fsc_par, "\n");
		}

		fclose(fsc_par);
	}

	// Вывод h
	{
		 int
			 len = strlen(TProblem::GlobalCatalog);
		 char *scourName = new char [len+50];
         string_copy(scourName, TProblem::GlobalCatalog, len+50);

		 char fileName[40];
         string_print(fileName, 40, "\\Scour\\TecPlot\\Scour%d.dat", t);
         string_concat(scourName, fileName, len+50);

         fsc = file_open(scourName, "w");
		 delete[] scourName;// scourName = NULL;
		 fprintf(fsc,"TITLE = Scour\n");
		 fprintf(fsc,"VARIABLES = X,Y,h,q_bx,q_by,uklon,P_ef\n");

		 fprintf(fsc,"ZONE T=lay%d, I=%d, J=%d, F=POINT\n", t, Nx, Ny);
		 for(int j=0; j<Ny; ++j)
		 {
			 for(int i=0; i<Nx; ++i)
			 {

				 fprintf(fsc,"%lf, %lf, %15.12f, %15.12f, %15.12f, %d, %15.12f\n", (Cx[i])*(TProblem::GlobalLength), (Cy[j])*(TProblem::GlobalLength), (h[i][j])*(TProblem::GlobalLength), q_b[i][j].x, q_b[i][j].y, gorisont[i][j], P_ef_print[i][j]);
			 }
		 }

		 fclose(fsc);
	}
}
//*************************************



//#################################################################


//*************************************
TNodesGrid* TBarrelProblem::BuildNormalGrid(TSeparator& lSep1, TSeparator& lSep2, TSeparator& lSep3, const TMask& lMask)
{

   int N = lSep1.EndIndex;
   int L = lSep2.EndIndex;
   int M = lSep3.EndIndex;

   T3dVector**** normals = new T3dVector***[N+1];

   for (int i=0;i<N+1;i++)
   {
     normals[i] = new T3dVector**[L+1];
     for (int j=0;j<L+1;j++)
     {
       normals[i][j] = new T3dVector*[M+1];
       for (int k=0;k<M+1;k++)
       {
         T3dVector* normal = NULL;
		 if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) ) {normal = NULL;}
         else
         {
            if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==0) )
			{
				normal = new T3dVector;
				normal->X = 0;
				normal->Y = 0;
				normal->Z = -1;
			}
            else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==M) )
            {
				normal = new T3dVector;
				normal->X = 0;
				normal->Y = 0;
				normal->Z = 1;
			}

            else if ( (i==0) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
            {
				normal = new T3dVector;
				normal->X = -1;
				normal->Y = 0;
				normal->Z = 0;
			}
            else if ( (i==N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
			{
				normal = new T3dVector;
				normal->X = 1;
				normal->Y = 0;
				normal->Z = 0;
			}

			else if ( (i!=0) && (i!=N) && (j==0) && (k!=0) && (k!=M) )
			{
				normal = new T3dVector;
				normal->X = 0;
				normal->Y = -1;
				normal->Z = 0;
			}
			else if ( (i!=0) && (i!=N) && (j==L) && (k!=0) && (k!=M) )
			{
				normal = new T3dVector;
				normal->X = 0;
				normal->Y = 1;
				normal->Z = 0;
			}
         }
		 normals[i][j][k] = normal;
        }
       }
      }

      TNodesGrid* lGrid = new T3dNormalGrid(&lSep1,&lSep2,&lSep3,lMask,normals);

	  for (int i=0;i<N+1;i++)
      {
       for (int j=0;j<L+1;j++)
	   {
        for (int k=0;k<M+1;k++)
		{
			if ( (normals[i][j][k])!=NULL ) {delete normals[i][j][k];}
		}
		delete[] normals[i][j];
	   }
	   delete[] normals[i];
	  }
	  delete[] normals;// normals = NULL;

	  return lGrid;
}
//*************************************
TNodesGrid* TBarrelProblem::BuildGridU()
{
   double h = fStep;

   TSeparator& sepU1 = *(new TUniformSeparator(0,fRadius*2,h));
   TSeparator& sepU2 = *(new TUniformSeparator(-h/2,fRadius*2+h/2,h));
   TSeparator& sepU3 = *(new TUniformSeparator(-h/2,fHeight+h/2,h));

   TSeparator& sep1 = *(new TRandomSeparator(sepU1.CountNodes,sepU1.Dimension));
   TSeparator& sep2 = *(new TRandomSeparator(sepU2.CountNodes,sepU2.Dimension));
   TSeparator& sep3 = *(new TRandomSeparator(sepU3.CountNodes,sepU3.Dimension));

   delete &sepU1;
   delete &sepU2;
   delete &sepU3;

   int N = sep1.EndIndex;
   int L = sep2.EndIndex;
   int M = sep3.EndIndex;

   T3dNumberMask &mask = *( new T3dNumberMask(N+1,L+1,M+1));

   for (int i=0;i<N+1;i++)
   {
     for (int j=0;j<L+1;j++)
     {
       for (int k=0;k<M+1;k++)
       {
		 mask[i][j][k]=TFictivePoint;
		 if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
         {
             mask[i][j][k]=TActualPoint;
         }
         else
         {
            if ( (i==0) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
            {
				mask[i][j][k]=TDefinedBorderPoint;
			}
            else if ( (i==N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==0) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==0) )
			{
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
            else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==M) )
            {
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
         }
       }
     }
   }

   TNodesGrid* grid = BuildNormalGrid(sep1,sep2,sep3,mask);

   delete &mask;

   return grid;
}
//*************************************
double TBarrelProblem::GetBorderConditionU(int lI, int lJ, int lK, int lN)
{
   return 0;
}
//*************************************
TNodesGrid* TBarrelProblem::BuildGridV()
{
   double h = fStep;

   TSeparator& sepU1 = *(new TUniformSeparator(-h/2,fRadius*2+h/2,h));
   TSeparator& sepU2 = *(new TUniformSeparator(0,fRadius*2,h));
   TSeparator& sepU3 = *(new TUniformSeparator(-h/2,fHeight+h/2,h));

   TSeparator& sep1 = *(new TRandomSeparator(sepU1.CountNodes,sepU1.Dimension));
   TSeparator& sep2 = *(new TRandomSeparator(sepU2.CountNodes,sepU2.Dimension));
   TSeparator& sep3 = *(new TRandomSeparator(sepU3.CountNodes,sepU3.Dimension));

   delete &sepU1;
   delete &sepU2;
   delete &sepU3;

   int N = sep1.EndIndex;
   int L = sep2.EndIndex;
   int M = sep3.EndIndex;

   T3dNumberMask &mask = *( new T3dNumberMask(N+1,L+1,M+1));

   for (int i=0;i<N+1;i++)
   {
     for (int j=0;j<L+1;j++)
     {
       for (int k=0;k<M+1;k++)
       {
		 mask[i][j][k]=TFictivePoint;
		 if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
         {
             mask[i][j][k]=TActualPoint;
         }
         else
         {
            if ( (i==0) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
            {
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
            else if ( (i==N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==0) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==0) )
			{
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
            else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==M) )
            {
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
         }
       }
     }
   }

   TNodesGrid* grid = BuildNormalGrid(sep1,sep2,sep3,mask);

   delete &mask;

   return grid;
}
//*************************************
double TBarrelProblem::GetBorderConditionV(int lI, int lJ, int lK, int lN)
{
   return 0;
}
//*************************************
TNodesGrid* TBarrelProblem::BuildGridW()
{
   double h = fStep;

   TSeparator& sepU1 = *(new TUniformSeparator(-h/2,fRadius*2+h/2,h));
   TSeparator& sepU2 = *(new TUniformSeparator(-h/2,fRadius*2+h/2,h));
   TSeparator& sepU3 = *(new TUniformSeparator(0,fHeight,h));

   TSeparator& sep1 = *(new TRandomSeparator(sepU1.CountNodes,sepU1.Dimension));
   TSeparator& sep2 = *(new TRandomSeparator(sepU2.CountNodes,sepU2.Dimension));
   TSeparator& sep3 = *(new TRandomSeparator(sepU3.CountNodes,sepU3.Dimension));

   delete &sepU1;
   delete &sepU2;
   delete &sepU3;

   int N = sep1.EndIndex;
   int L = sep2.EndIndex;
   int M = sep3.EndIndex;

   T3dNumberMask &mask = *( new T3dNumberMask(N+1,L+1,M+1));

   for (int i=0;i<N+1;i++)
   {
     for (int j=0;j<L+1;j++)
     {
       for (int k=0;k<M+1;k++)
       {
		 mask[i][j][k]=TFictivePoint;
		 if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
         {
             mask[i][j][k]=TActualPoint;
         }
         else
         {
			int tab = 3;

			if ( (i==0) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
            {
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
            else if ( (i==N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==0) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==0) )
			{
				if ( (i<0+tab) || (i>N-tab) || (j<tab) || (j>L-tab) )
				{
                  mask[i][j][k]=TDefinedBorderPoint;
				}
				else
				{
				  mask[i][j][k]=TEquationBorderPoint;
				}
			}
            else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==M) )
            {
				if ( (i<0+tab) || (i>N-tab) || (j<tab) || (j>L-tab) )
				{
                  mask[i][j][k]=TDefinedBorderPoint;
				}
				else
				{
				  mask[i][j][k]=TEquationBorderPoint;
				}
			}
         }
       }
     }
   }

   TNodesGrid* grid = BuildNormalGrid(sep1,sep2,sep3,mask);

   delete &mask;

   return grid;
}
//*************************************
double TBarrelProblem::GetBorderConditionW(int lI, int lJ, int lK, int lN)
{
   return 0;
}
//*************************************
TNodesGrid* TBarrelProblem::BuildGridP()
{
   double h = fStep;
   TSeparator& sepU1 = *(new TUniformSeparator(-h/2,fRadius*2+h/2,h));
   TSeparator& sepU2 = *(new TUniformSeparator(-h/2,fRadius*2+h/2,h));
   TSeparator& sepU3 = *(new TUniformSeparator(-h/2,fHeight+h/2,h));

   TSeparator& sep1 = *(new TRandomSeparator(sepU1.CountNodes,sepU1.Dimension));
   delete &sepU1;
   TSeparator& sep2 = *(new TRandomSeparator(sepU2.CountNodes,sepU2.Dimension));
   delete &sepU2;
   TSeparator& sep3 = *(new TRandomSeparator(sepU3.CountNodes,sepU3.Dimension));
   delete &sepU3;

   int N = sep1.EndIndex;
   int L = sep2.EndIndex;
   int M = sep3.EndIndex;

   T3dNumberMask &mask = *( new T3dNumberMask(N+1,L+1,M+1));

   for (int i=0;i<N+1;i++)
   {
     for (int j=0;j<L+1;j++)
     {
       for (int k=0;k<M+1;k++)
       {
		 mask[i][j][k]=TFictivePoint;
		 if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
         {
             mask[i][j][k]=TActualPoint;
         }
         else
         {
            int tab = 3;

            if ( (i==0) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
            {
				mask[i][j][k]=TPreNormalBorderPoint;
			}
            else if ( (i==N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreNormalBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==0) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreNormalBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreNormalBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==0) )
			{
				if ( (i<0+tab) || (i>N-tab) || (j<tab) || (j>L-tab) )
				{
                  mask[i][j][k]=TPreNormalBorderPoint;
				}
				else
				{
				  mask[i][j][k]=TPreDefinedBorderPoint;
				}
			}
            else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==M) )
            {
				if ( (i<0+tab) || (i>N-tab) || (j<tab) || (j>L-tab) )
				{
                  mask[i][j][k]=TPreNormalBorderPoint;
				}
				else
				{
				  mask[i][j][k]=TPreDefinedBorderPoint;
				}
			}
         }
       }
     }
   }

   TNodesGrid* grid = BuildNormalGrid(sep1,sep2,sep3,mask);

   delete &mask;

   return grid;
}
//*************************************
double TBarrelProblem::GetBorderConditionP(int lI, int lJ, int lK, int lN)
{
   if (lK==0)  {return 1;}
   else {return 0;}
}
//*************************************
TNodesGrid*  TBarrelProblem::BuildGrid()
{
   double h = fStep;
   TSeparator& sepU1 = *(new TUniformSeparator(0,fRadius*2,h));
   TSeparator& sepU2 = *(new TUniformSeparator(0,fRadius*2,h));
   TSeparator& sepU3 = *(new TUniformSeparator(0,fHeight,h));

   TSeparator& sep1 = *(new TRandomSeparator(sepU1.CountNodes,sepU1.Dimension));
   delete &sepU1;
   TSeparator& sep2 = *(new TRandomSeparator(sepU2.CountNodes,sepU2.Dimension));
   delete &sepU2;
   TSeparator& sep3 = *(new TRandomSeparator(sepU3.CountNodes,sepU3.Dimension));
   delete &sepU3;

   int N = sep1.EndIndex;
   int L = sep2.EndIndex;
   int M = sep3.EndIndex;

   T3dNumberMask &mask = *( new T3dNumberMask(N+1,L+1,M+1));

   for (int i=0;i<N+1;i++)
   {
     for (int j=0;j<L+1;j++)
     {
       for (int k=0;k<M+1;k++)
       {
		 mask[i][j][k]=TFictivePoint;
		 if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
         {
             mask[i][j][k]=TActualPoint;
         }
         else
         {
            if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==0) )
			{
				mask[i][j][k]=TBorderPoint;
			}
            else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==M) )
            {
				mask[i][j][k]=TBorderPoint;
			}

            else if ( (i==0) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
            {
				mask[i][j][k]=TBorderPoint;
			}
            else if ( (i==N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TBorderPoint;
			}

			else if ( (i!=0) && (i!=N) && (j==0) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TBorderPoint;
			}
         }
       }
	 }
   }

   TNodesGrid* grid = BuildNormalGrid(sep1,sep2,sep3,mask);

   delete &mask;

   return grid;
}
//*************************************




//#################################################################




//*************************************
TNodesGrid* TCanalProblem::BuildNormalGrid(TSeparator& lSep1, TSeparator& lSep2, TSeparator& lSep3, const TMask& lMask)
{

   int N = lSep1.EndIndex;
   int L = lSep2.EndIndex;
   int M = lSep3.EndIndex;

   T3dVector**** normals = new T3dVector***[N+1];

   for (int i=0;i<N+1;i++)
   {
     normals[i] = new T3dVector**[L+1];
     for (int j=0;j<L+1;j++)
     {
       normals[i][j] = new T3dVector*[M+1];
       for (int k=0;k<M+1;k++)
       {
         T3dVector* normal = NULL;
		 if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) ) {normal = NULL;}
         else
         {
            if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==0) )
			{
				normal = new T3dVector;
				normal->X = 0;
				normal->Y = 0;
				normal->Z = -1;
			}
            else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==M) )
            {
				normal = new T3dVector;
				normal->X = 0;
				normal->Y = 0;
				normal->Z = 1;
			}

            else if ( (i==0) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
            {
				normal = new T3dVector;
				normal->X = -1;
				normal->Y = 0;
				normal->Z = 0;
			}
            else if ( (i==N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
			{
				normal = new T3dVector;
				normal->X = 1;
				normal->Y = 0;
				normal->Z = 0;
			}

			else if ( (i!=0) && (i!=N) && (j==0) && (k!=0) && (k!=M) )
			{
				normal = new T3dVector;
				normal->X = 0;
				normal->Y = -1;
				normal->Z = 0;
			}
			else if ( (i!=0) && (i!=N) && (j==L) && (k!=0) && (k!=M) )
			{
				normal = new T3dVector;
				normal->X = 0;
				normal->Y = 1;
				normal->Z = 0;
			}
         }
		 normals[i][j][k] = normal;
        }
       }
      }

      TNodesGrid* lGrid = new T3dNormalGrid(&lSep1,&lSep2,&lSep3,lMask,normals);

	  for (int i=0;i<N+1;i++)
      {
       for (int j=0;j<L+1;j++)
	   {
        for (int k=0;k<M+1;k++)
		{
			if ( (normals[i][j][k])!=NULL ) {delete normals[i][j][k];}
		}
		delete[] normals[i][j];
	   }
	   delete[] normals[i];
	  }
	  delete[] normals;// normals = NULL;

	  return lGrid;
}
//*************************************
TNodesGrid* TCanalProblem::BuildGridU()
{
   double h = fStep;

   TSeparator& sepU1 = *(new TUniformSeparator(0,fLength,h));
   TSeparator& sepU2 = *(new TUniformSeparator(-h/2,fRadius*2+h/2,h));
   TSeparator& sepU3 = *(new TUniformSeparator(-h/2,fRadius*2+h/2,h));

   TSeparator& sep1 = *(new TRandomSeparator(sepU1.CountNodes,sepU1.Dimension));
   TSeparator& sep2 = *(new TRandomSeparator(sepU2.CountNodes,sepU2.Dimension));
   TSeparator& sep3 = *(new TRandomSeparator(sepU3.CountNodes,sepU3.Dimension));

   delete &sepU1;
   delete &sepU2;
   delete &sepU3;

   int N = sep1.EndIndex;
   int L = sep2.EndIndex;
   int M = sep3.EndIndex;

   T3dNumberMask &mask = *( new T3dNumberMask(N+1,L+1,M+1));

   for (int i=0;i<N+1;i++)
   {
     for (int j=0;j<L+1;j++)
     {
       for (int k=0;k<M+1;k++)
       {
		 mask[i][j][k]=TFictivePoint;
		 if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
         {
             mask[i][j][k]=TActualPoint;
         }
         else
         {
            int tab = 2;
            if ( (i==0) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
            {
				if ( (j<tab) || (j>L-tab) || (k<tab) || (k>M-tab) )
				{
				   mask[i][j][k]=TDefinedBorderPoint;
				}
				else
				{
                   mask[i][j][k]=TEquationBorderPoint;
				}
			}
            else if ( (i==N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
			{
				if ( (j<tab) || (j>L-tab) || (k<tab) || (k>M-tab) )
				{
				   mask[i][j][k]=TDefinedBorderPoint;
				}
				else
				{
                   mask[i][j][k]=TEquationBorderPoint;
				}
			}
			else if ( (i!=0) && (i!=N) && (j==0) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==0) )
			{
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
            else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==M) )
            {
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
         }
       }
     }
   }

   TNodesGrid* grid = BuildNormalGrid(sep1,sep2,sep3,mask);

   delete &mask;

   return grid;
}
//*************************************
double TCanalProblem::GetBorderConditionU(int lI, int lJ, int lK, int lN)
{
   return 0;
}
//*************************************
TNodesGrid* TCanalProblem::BuildGridV()
{
   double h = fStep;

   TSeparator& sepU1 = *(new TUniformSeparator(-h/2,fLength+h/2,h));
   TSeparator& sepU2 = *(new TUniformSeparator(0,fRadius*2,h));
   TSeparator& sepU3 = *(new TUniformSeparator(-h/2,fRadius*2+h/2,h));

   TSeparator& sep1 = *(new TRandomSeparator(sepU1.CountNodes,sepU1.Dimension));
   TSeparator& sep2 = *(new TRandomSeparator(sepU2.CountNodes,sepU2.Dimension));
   TSeparator& sep3 = *(new TRandomSeparator(sepU3.CountNodes,sepU3.Dimension));

   delete &sepU1;
   delete &sepU2;
   delete &sepU3;

   int N = sep1.EndIndex;
   int L = sep2.EndIndex;
   int M = sep3.EndIndex;

   T3dNumberMask &mask = *( new T3dNumberMask(N+1,L+1,M+1));

   for (int i=0;i<N+1;i++)
   {
     for (int j=0;j<L+1;j++)
     {
       for (int k=0;k<M+1;k++)
       {
		 mask[i][j][k]=TFictivePoint;
		 if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
         {
             mask[i][j][k]=TActualPoint;
         }
         else
         {
            if ( (i==0) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
            {
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
            else if ( (i==N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==0) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==0) )
			{
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
            else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==M) )
            {
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
         }
       }
     }
   }

   TNodesGrid* grid = BuildNormalGrid(sep1,sep2,sep3,mask);

   delete &mask;

   return grid;
}
//*************************************
double TCanalProblem::GetBorderConditionV(int lI, int lJ, int lK, int lN)
{
   return 0;
}
//*************************************
TNodesGrid* TCanalProblem::BuildGridW()
{
   double h = fStep;

   TSeparator& sepU1 = *(new TUniformSeparator(-h/2,fLength+h/2,h));
   TSeparator& sepU2 = *(new TUniformSeparator(-h/2,fRadius*2+h/2,h));
   TSeparator& sepU3 = *(new TUniformSeparator(0,fRadius*2,h));

   TSeparator& sep1 = *(new TRandomSeparator(sepU1.CountNodes,sepU1.Dimension));
   TSeparator& sep2 = *(new TRandomSeparator(sepU2.CountNodes,sepU2.Dimension));
   TSeparator& sep3 = *(new TRandomSeparator(sepU3.CountNodes,sepU3.Dimension));

   delete &sepU1;
   delete &sepU2;
   delete &sepU3;

   int N = sep1.EndIndex;
   int L = sep2.EndIndex;
   int M = sep3.EndIndex;

   T3dNumberMask &mask = *( new T3dNumberMask(N+1,L+1,M+1));

   for (int i=0;i<N+1;i++)
   {
     for (int j=0;j<L+1;j++)
     {
       for (int k=0;k<M+1;k++)
       {
		 mask[i][j][k]=TFictivePoint;
		 if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
         {
             mask[i][j][k]=TActualPoint;
         }
         else
         {
			if ( (i==0) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
            {
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
            else if ( (i==N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==0) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreDefinedBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==0) )
			{
                mask[i][j][k]=TDefinedBorderPoint;
			}
            else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==M) )
            {
                mask[i][j][k]=TDefinedBorderPoint;
			}
         }
       }
     }
   }

   TNodesGrid* grid = BuildNormalGrid(sep1,sep2,sep3,mask);

   delete &mask;

   return grid;
}
//*************************************
double TCanalProblem::GetBorderConditionW(int lI, int lJ, int lK, int lN)
{
   return 0;
}
//*************************************
TNodesGrid* TCanalProblem::BuildGridP()
{
   double h = fStep;
   TSeparator& sepU1 = *(new TUniformSeparator(-h/2,fRadius*2+h/2,h));
   TSeparator& sepU2 = *(new TUniformSeparator(-h/2,fRadius*2+h/2,h));
   TSeparator& sepU3 = *(new TUniformSeparator(-h/2,fLength+h/2,h));

   TSeparator& sep1 = *(new TRandomSeparator(sepU1.CountNodes,sepU1.Dimension));
   delete &sepU1;
   TSeparator& sep2 = *(new TRandomSeparator(sepU2.CountNodes,sepU2.Dimension));
   delete &sepU2;
   TSeparator& sep3 = *(new TRandomSeparator(sepU3.CountNodes,sepU3.Dimension));
   delete &sepU3;

   int N = sep1.EndIndex;
   int L = sep2.EndIndex;
   int M = sep3.EndIndex;

   T3dNumberMask &mask = *( new T3dNumberMask(N+1,L+1,M+1));

   for (int i=0;i<N+1;i++)
   {
     for (int j=0;j<L+1;j++)
     {
       for (int k=0;k<M+1;k++)
       {
		 mask[i][j][k]=TFictivePoint;
		 if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
         {
             mask[i][j][k]=TActualPoint;
         }
         else
         {
            int tab = 2;
            if ( (i==0) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
            {
				if ( (j<tab) || (j>L-tab) || (k<tab) || (k>M-tab) )
				{
				   mask[i][j][k]=TPreNormalBorderPoint;
				}
				else
				{
                    mask[i][j][k]=TPreDefinedBorderPoint;
				}
			}
            else if ( (i==N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
			{
				if ( (j<tab) || (j>L-tab) || (k<tab) || (k>M-tab) )
				{
				   mask[i][j][k]=TPreNormalBorderPoint;
				}
				else
				{
                    mask[i][j][k]=TPreDefinedBorderPoint;
				}
			}
			else if ( (i!=0) && (i!=N) && (j==0) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreNormalBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TPreNormalBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==0) )
			{
                mask[i][j][k]=TPreNormalBorderPoint;
			}
            else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==M) )
            {
				mask[i][j][k]=TPreNormalBorderPoint;
			}
         }
       }
     }
   }

   TNodesGrid* grid = BuildNormalGrid(sep1,sep2,sep3,mask);

   delete &mask;

   return grid;
}
//*************************************
double TCanalProblem::GetBorderConditionP(int lI, int lJ, int lK, int lN)
{
   if (lI==0)  {return 1;}
   else {return 0;}
}
//*************************************
TNodesGrid*  TCanalProblem::BuildGrid()
{
   double h = fStep;
   TSeparator& sepU1 = *(new TUniformSeparator(0,fLength,h));
   TSeparator& sepU2 = *(new TUniformSeparator(0,fRadius*2,h));
   TSeparator& sepU3 = *(new TUniformSeparator(0,fRadius*2,h));

   TSeparator& sep1 = *(new TRandomSeparator(sepU1.CountNodes,sepU1.Dimension));
   delete &sepU1;
   TSeparator& sep2 = *(new TRandomSeparator(sepU2.CountNodes,sepU2.Dimension));
   delete &sepU2;
   TSeparator& sep3 = *(new TRandomSeparator(sepU3.CountNodes,sepU3.Dimension));
   delete &sepU3;

   int N = sep1.EndIndex;
   int L = sep2.EndIndex;
   int M = sep3.EndIndex;

   T3dNumberMask &mask = *( new T3dNumberMask(N+1,L+1,M+1));

   for (int i=0;i<N+1;i++)
   {
     for (int j=0;j<L+1;j++)
     {
       for (int k=0;k<M+1;k++)
       {
		 mask[i][j][k]=TFictivePoint;
		 if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
         {
             mask[i][j][k]=TActualPoint;
         }
         else
         {
            if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==0) )
			{
				mask[i][j][k]=TBorderPoint;
			}
            else if ( (i!=0) && (i!=N) && (j!=0) && (j!=L) && (k==M) )
            {
				mask[i][j][k]=TBorderPoint;
			}

            else if ( (i==0) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
            {
				mask[i][j][k]=TBorderPoint;
			}
            else if ( (i==N) && (j!=0) && (j!=L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TBorderPoint;
			}

			else if ( (i!=0) && (i!=N) && (j==0) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TBorderPoint;
			}
			else if ( (i!=0) && (i!=N) && (j==L) && (k!=0) && (k!=M) )
			{
				mask[i][j][k]=TBorderPoint;
			}
         }
       }
	 }
   }

   TNodesGrid* grid = BuildNormalGrid(sep1,sep2,sep3,mask);

   delete &mask;

   return grid;
}
//*************************************




//#################################################################



TMaskGenerator* TCanalUnifyProblem::CreateStatement()
{
    //Preparing (Begin)
	double iPart;
 	int endIndex;

    iPart = 0;
	modf(fLength/fStep,&iPart);
    endIndex = (int)iPart;
	int nx = endIndex + 1;

	iPart = 0;
	modf(fRadius*2/fStep,&iPart);
    endIndex = (int)iPart;
	int ny = endIndex + 1;
	int nz = endIndex + 1;
	//Preparing (End)

	//Creating (Begin)
	return TSimpleWithPressureMaskGenerator::CreateInstance(
		                                                    fLength,fRadius*2,fRadius*2,
		                                                    nx,ny,nz,
															fPressureInput,fPressureOutput
														   );
	//Creating (End)
}



//#################################################################




TMaskGenerator* TFlowPastCubeProblem::CreateStatement()
{
    //Preparing (Begin)
	double iPart;
 	int endIndex;

    iPart = 0;
	modf(fLength/fStep,&iPart);
    endIndex = (int)iPart;
	int nx = endIndex + 1;

	iPart = 0;
	modf(fWidth/fStep,&iPart);
    endIndex = (int)iPart;
	int ny = endIndex + 1;

	iPart = 0;
	modf(fHeight/fStep,&iPart);
    endIndex = (int)iPart;
	int nz = endIndex + 1;

	iPart = 0;
	modf(fCubeSize/fStep,&iPart);
    endIndex = (int)iPart;
	int ns = endIndex + 1;

	int nx_2 = (nx-ns)/2+1;
	int ny_2 = (ny-ns)/2+1;
	int nz_2 = (nz-ns)/2+1;
	//Preparing (End)


	//Creating (Begin)
	return TCubeInsideNewDivWithPressureMaskGenerator::CreateInstance(fLength,fWidth,fHeight,
		                                                              nx_2,ns,nx_2,
																      ny_2,ns,ny_2,
																      nz_2,ns,nz_2,
																      fPressureInput,fPressureOutput);
	//Creating (End)
}




//#################################################################


TMaskGenerator* TBargeProblem::CreateStatement()
{
	double h = 0;

	if (fBargeAngle==0)
	{
	  return TInvertedTrapezeBottomAlongXWithPressureMaskGenerator::CreateInstance
	  (
	   fAreaLength,fAreaWidth,fAreaHeight,fBargeWidth,fBargeLengthBottom,fBargeLengthTop,h,
	   fAnchorX,fAnchorY,fGridSpace,fGridSpace,fGridSpace,fPressureDifference,0
	  );
	}
	else if (fBargeAngle==90)
	{
	  return TInvertedTrapezeBottomAlongYWithPressureMaskGenerator::CreateInstance
	  (
	   fAreaLength,fAreaWidth,fAreaHeight,fBargeWidth,fBargeLengthBottom,fBargeLengthTop,h,
	   fAnchorX,fAnchorY,fGridSpace,fGridSpace,fGridSpace,fPressureDifference,0
	  );
	}
	else
	{
       throw "Error in TBargeProblem::CreateStatement: Unknown barge angle.";
	}
}


//#################################################################


TMaskGenerator* TPlatformProblem::CreateStatement()
{
	 return TDoubleTrapezeBottomWithPressureMaskGenerator::CreateInstance
	  (
	   fAreaLength,fAreaWidth,fAreaHeight,fPlatformSizeBottom,fPlatformSizeTop,
	   fPlatformHeightBottom,fPlatformHeightTop,
	   fAnchorX,fAnchorY,fGridSpace,fGridSpace,fGridSpace,fPressureDifference,0
	  );
}




//#################################################################


TMaskGenerator* TBackwardFacingStepProblem::CreateStatement()
{
    //Preparing (Begin)
	double iPart;
 	int endIndex;

    iPart = 0;
	modf(fLength/fStep,&iPart);
    endIndex = (int)iPart;
	int nx = endIndex + 1;

	iPart = 0;
	modf(fWidth/fStep,&iPart);
    endIndex = (int)iPart;
	int ny = endIndex + 1;

	iPart = 0;
	modf(fHeight/fStep,&iPart);
    endIndex = (int)iPart;
	int nz = endIndex + 1;

	iPart = 0;
	modf(fStepSize/fStep,&iPart);
    endIndex = (int)iPart;
	int ns = endIndex + 1;
	//Preparing (End)


	//Creating (Begin)
	return TBackwardFacingStepMaskGenerator::CreateInstance(fLength,fWidth,fHeight,
		                                                     ns,nx-ns+1,
														     ny,
															 ns,nz-ns+1,
															 fPressureInput,fPressureOutput);

	//Creating (End)
}


TImmersedBoundaryProblem::TImmersedBoundaryProblem(
        double lTimeEndValue,
        double lTimeStepValue,
        double lRe,
        double lStiffness
        ): TVPProblem(lTimeEndValue, lTimeStepValue, lRe),
	fStiffness(lStiffness) {}

TImmersedBoundaryProblem::TImmersedBoundaryProblem(
        double lTimeEndValue,
        double lTimeStepValue,
        double lmu1,
        double lmu2,
        double lrho1,
        double lrho2,
        double lStiffness
        ): TVPProblem(lTimeEndValue, lTimeStepValue, lmu1, lmu2, lrho1, lrho2),
	fStiffness(lStiffness) {}

TImmersedBoundaryProblem::~TImmersedBoundaryProblem()
{
	delete[]slice_x_0;// slice_x_0 = NULL;
	delete[]slice_y_0;// slice_y_0 = NULL;
	delete[]slice_z_0;// slice_z_0 = NULL;
    delete fBoundary;
}

void TImmersedBoundaryProblem::Initialize()
{	
	int lNodesCount = 100;
	double lRadius = 0.1;
	double lCenter_x = 0.5;
	double lCenter_y = 0.25;
	double lCenter_z = 0.25;

    TVPProblem::Initialize();
    fBoundary = (TSphereBoundary *)(new TSphereBoundary())
        ->withStiffness(fStiffness)
        ->withNodesCount(lNodesCount)
        ->withRadius(lRadius)
        ->withCenter(lCenter_x, lCenter_y, lCenter_z)
        ->Initialize();

}

//Here
void TImmersedBoundaryProblem::GetForce(
        TRnRelease3dSpace &force,
        T3dNormalGrid &gridC,
        T3dNormalGrid &grid,
        TRnRelease3dSpace &velocity,
        const int axis,
        const int timeStepNumber
        )
{	
    Interpolate(velocity, &grid, axis, timeStepNumber);
    UpdateBoundaryPosition(axis, timeStepNumber);
    ComputeBoundaryForces(axis, timeStepNumber);
    SpreadForce(force, gridC, &grid, axis, timeStepNumber);
	InterpolateUpdatePathes(velocity, &grid, axis, timeStepNumber);
}

double TImmersedBoundaryProblem::GetDistribution(int _i, int _j, int _k, double hx, double hy, double hz)//, double h)
{
	// int i,j,k
	// h - ?
	// Find point on base mesh xm, ym, zm
	double xm = _i*hx-hx*0.5;
	double ym = _j*hy-hy*0.5;
	double zm = _k*hz-hz*0.5;


	// Find point xs, ys, zs on spline the nearest to xm, ym, zm
	//while (!(xm>=fBoundary->))
	double xs, ys, zs;
	int i = 0;
	double x0, x1, Alpha;
	if (xm < fBoundary->xSpl[0])
	{
		xs = fBoundary->xSpl[0];
		ys = fBoundary->ySpl[0];
		zs = fBoundary->zSpl[0];
	}
	else if (xm > fBoundary->xSpl[fBoundary->nSpl - 1])
	{
		xs = fBoundary->xSpl[fBoundary->nSpl - 1];
		ys = fBoundary->ySpl[fBoundary->nSpl - 1];
		zs = fBoundary->zSpl[fBoundary->nSpl - 1];
	}
	else
	{
	while (!(xm >= fBoundary->xSpl[i] && xm <= fBoundary->xSpl[i+1]) && i < fBoundary->nSpl - 1)
		i++;

		if (xm == fBoundary->xSpl[i])
		{
			xs = fBoundary->xSpl[i];
			ys = fBoundary->ySpl[i];
			zs = fBoundary->zSpl[i];
		}
		else if (xm == fBoundary->xSpl[i+1])
		{
			xs = fBoundary->xSpl[i+1];
			ys = fBoundary->ySpl[i+1];
			zs = fBoundary->zSpl[i+1];
		}
		else
		{
			x0 = fBoundary->xSpl[i];
			x1 = fBoundary->xSpl[i+1];
			Alpha = (xm - x0) / (x1-x0);
			xs = (1.0 - Alpha)*fBoundary->xSpl[i] + Alpha*fBoundary->xSpl[i+1];
			ys = (1.0 - Alpha)*fBoundary->ySpl[i] + Alpha*fBoundary->ySpl[i + 1];
			zs = (1.0 - Alpha)*fBoundary->zSpl[i] + Alpha*fBoundary->zSpl[i + 1];
		}
		//i++;
	}

	//
	double r_sm = sqrt(pow(ym-ys,2)+pow(zm-zs,2));
	//if (r_sm <= fBoundary->radius) printf("Distr: i,j,k=%d,%d,%d = %lf\n", _i, _j, _k, r_sm);

	double r = (fBoundary->radius - r_sm) / fBoundary->radius;
	if (r<0 || r>1) r = 0.00;
	return r;
}


void TImmersedBoundaryProblem::ChangeStifness(const int timeStepNumber)
{
	if (timeStepNumber == 80)
	{
		for (int i = 0; i < fBoundary->nodesCount; i++)
		{
			fBoundary->nodes[i]->stretchingStiffness = fBoundary->stiffness;
		}
	}

}


void TImmersedBoundaryProblem::ComputeBoundaryForces(const int axis, const int timeStepNumber)
{
    std::function<void(TNode *)> compute;
    std::function<void(TNode *)> resetForce;
    switch(axis) {
        case COORD_X:
            compute = [this](TNode *node) {
                double stretchingForce = node->GetStretchingForce(node->neighbors.next, TNode::Ox);
                double bendingForce = node->GetBendingForce(node->neighbors.next, node->neighbors.prev, TNode::Ox);
                double targetForce = node->GetTargetForce(TNode::Ox);

                node->xForce += -(stretchingForce + 2 * bendingForce + targetForce) * fBoundary->GetArea();

                if (node->neighbors.next != NULL)
                {
                    node->neighbors.next->xForce += (stretchingForce + bendingForce) * fBoundary->GetArea();
                }
                else
                {
                    assert(stretchingForce == 0);
                    assert(bendingForce == 0);
                }

                if (node->neighbors.prev != NULL)
                {
                    node->neighbors.prev->xForce += bendingForce * fBoundary->GetArea();
                }
                else
                {
                    assert(bendingForce == 0);
                }
            };
            resetForce = [this](TNode *node) {
                node->xForce = 0;
            };
            break;
        case COORD_Y:
            compute = [this](TNode *node) {
                double stretchingForce = node->GetStretchingForce(node->neighbors.next, TNode::Oy);
                double bendingForce = node->GetBendingForce(node->neighbors.next, node->neighbors.prev, TNode::Oy);
                double targetForce = node->GetTargetForce(TNode::Oy);

                node->yForce += -(stretchingForce + 2 * bendingForce + targetForce) * fBoundary->GetArea();
                if (node->neighbors.next != NULL)
                {
                    node->neighbors.next->yForce += (stretchingForce + bendingForce) * fBoundary->GetArea();
                }
                else
                {
                    assert(stretchingForce == 0);
                    assert(bendingForce == 0);
                }

                if (node->neighbors.prev != NULL)
                {
                    node->neighbors.prev->yForce += bendingForce * fBoundary->GetArea();
                }
                else
                {
                    assert(bendingForce == 0);
                }
            };
            resetForce = [this](TNode *node) {
                node->yForce = 0;
            };
            break;
        case COORD_Z:
            compute = [this](TNode *node) {
                double stretchingForce = node->GetStretchingForce(node->neighbors.next, TNode::Oz);
                double bendingForce = node->GetBendingForce(node->neighbors.next, node->neighbors.prev, TNode::Oz);
                double targetForce = node->GetTargetForce(TNode::Oz);

                node->zForce += -(stretchingForce + 2 * bendingForce + targetForce) * fBoundary->GetArea();
                if (node->neighbors.next != NULL)
                {
                    node->neighbors.next->zForce += (stretchingForce + bendingForce) * fBoundary->GetArea();
                }
                else
                {
                    assert(stretchingForce == 0);
                    assert(bendingForce == 0);
                }

                if (node->neighbors.prev != NULL)
                {
                    node->neighbors.prev->zForce += bendingForce * fBoundary->GetArea();
                }
                else
                {
                    assert(bendingForce == 0);
                }
            };
            resetForce = [this](TNode *node) {
                node->zForce = 0;
            };
            break;
        default:
            throw "Incorrect type of axis";
    }

    #pragma omp parallel for
    for(int n = 0; n < fBoundary->nodesCount; ++n) {
        resetForce(fBoundary->nodes[n]);
    }

    #pragma omp parallel for
    for(int n = 0; n < fBoundary->nodesCount; ++n) {
        compute(fBoundary->nodes[n]);
    }

    return;
}


void TImmersedBoundaryProblem::SpreadForce(TRnRelease3dSpace &force, T3dNormalGrid &gridC, T3dNormalGrid *grid, const int axis, const int timeStepNumber)
{	
    force |= 0.0;
    std::function<double(TNode *)> getNodeForce;
    switch(axis) {
        case COORD_X:
            getNodeForce = [](TNode *node) {
                return node->xForce;
            };
            break;

        case COORD_Y:
            getNodeForce = [](TNode *node) {
                return node->yForce;
            };
            break;

        case COORD_Z:
            getNodeForce = [](TNode *node) {
                return node->zForce;
            };
            break;

        default:
            throw "Incorrect type of axis";
    }

	const T3dNumberMask& mask = gridC.Mask;
    for(int n = 0; n < fBoundary->nodesCount; ++n) {
		int x_int = Index(grid, fBoundary->nodes[n]->x, COORD_X, axis, RANGE, RANGE);
		int y_int = Index(grid, fBoundary->nodes[n]->y, COORD_Y, axis, RANGE, RANGE);
		int z_int = Index(grid, fBoundary->nodes[n]->z, COORD_Z, axis, RANGE, RANGE);

        for(int i = x_int - RANGE; i <= x_int + RANGE; ++i) {
            for(int j = y_int - RANGE; j <= y_int + RANGE; ++j) {
                for(int k = z_int - RANGE; k <= z_int + RANGE; ++k) {


                    const double dist_x = fabs(fBoundary->nodes[n]->x - Coord(grid, i, COORD_X, axis));
                    const double dist_y = fabs(fBoundary->nodes[n]->y - Coord(grid, j, COORD_Y, axis));
                    const double dist_z = fabs(fBoundary->nodes[n]->z - Coord(grid, k, COORD_Z, axis));

                    const double weight_x = Dirichlet(grid, dist_x, COORD_X);
                    const double weight_y = Dirichlet(grid, dist_y, COORD_Y);
                    const double weight_z = Dirichlet(grid, dist_z, COORD_Z);

                    force[i][j][k] += (getNodeForce(fBoundary->nodes[n]) * weight_x * weight_y * weight_z);
                    if (
                        (fBoundary->nodes[n]->type == FIXED_NODE) &&
                        (i == x_int || i == x_int + 1) &&
                        (
                         (j <= 12) && (j == y_int + 1 || j == y_int + 2) ||
                         (j > 12) && (j == y_int || j == y_int + 1)
                        ) &&
                        (
                         (k <= 12) && (k == z_int + 1 || k == z_int + 2) ||
                         (k > 12) && (k == z_int || k == z_int + 1)
                        )
                       )
                    {
                        mask[i][j][k] = TPreNormalBorderPoint;
                    }
                }
            }
        }

    }
}


void TImmersedBoundaryProblem::Interpolate(TRnRelease3dSpace &velocity, T3dNormalGrid *grid, const int axis, const int timeStepNumber)
{
    const TSeparator *separator = NULL;
    std::function<void(TNode *)> resetVelocity;
    std::function<void(TNode *, double)> updateVelocity;

    // update_debug leads to segfault with multiple threads???
    std::function<void(
            TNode *, double,
            double, double, double)> updateDebug;

    switch(axis) {
        case COORD_X:
            separator = &grid->GetSeparator1();
            resetVelocity = [this](TNode *node) {
                node->xVel = 0;
            };
            updateVelocity = [this](TNode *node, double value) {
                node->xVel += value;
		//LOG(INFO) << "x vel = " << node->xVel;
            };
            updateDebug = [this, grid](TNode *node, double vel, int i, int j, int k) {
                tuple<int, int, int> index(i, j, k);

                TNode *neighbor_node = NULL;
                if(node->surround.find(index) == node->surround.end())
                {
                    // there is no TNode with this indexes
                    neighbor_node = new TNode();
                    neighbor_node->xVel = vel;
                    neighbor_node->x = Coord(grid, i, COORD_X);
                    neighbor_node->y = Coord(grid, j, COORD_Y);
                    neighbor_node->z = Coord(grid, k, COORD_Z);

                    node->surround[index] = neighbor_node;
                }
                else
                {
                    neighbor_node = node->surround[index];
                    neighbor_node->xVel = vel;
                }
            };
            break;
        case COORD_Y:
            separator = &grid->GetSeparator2();
            resetVelocity = [this](TNode *node) {
                node->yVel = 0;
            };
            updateVelocity = [this](TNode *node, double value) {
                node->yVel += value;
		//LOG(INFO) << "y vel = " << node->yVel;
            };
            updateDebug = [this, grid](TNode *node, double vel, int i, int j, int k) {
                tuple<int, int, int> index(i, j, k);

                TNode *neighbor_node = NULL;
                if(node->surround.find(index) == node->surround.end())
                {
                    // there is no TNode with this indexes
                    neighbor_node = new TNode();
                    neighbor_node->yVel = vel;
                    neighbor_node->x = Coord(grid, i, COORD_X);
                    neighbor_node->y = Coord(grid, j, COORD_Y);
                    neighbor_node->z = Coord(grid, k, COORD_Z);

                    node->surround[index] = neighbor_node;
                }
                else
                {
                    neighbor_node = node->surround[index];
                    neighbor_node->yVel = vel;
                }
            };
            break;
        case COORD_Z:
            separator = &grid->GetSeparator3();
            resetVelocity = [this](TNode *node) {
                node->zVel = 0;
            };
            updateVelocity = [this](TNode *node, double value) {
                node->zVel += value;
		//LOG(INFO) << "z vel = " << node->zVel;
            };
            updateDebug = [this, grid](TNode *node, double vel, int i, int j, int k) {
                tuple<int, int, int> index(i, j, k);

                TNode *neighbor_node = NULL;
                if(node->surround.find(index) == node->surround.end())
                {
                    // there is no TNode with this indexes
                    neighbor_node = new TNode();
                    neighbor_node->zVel = vel;
                    neighbor_node->x = Coord(grid, i, COORD_X);
                    neighbor_node->y = Coord(grid, j, COORD_Y);
                    neighbor_node->z = Coord(grid, k, COORD_Z);

                    node->surround[index] = neighbor_node;
                }
                else
                {
                    neighbor_node = node->surround[index];
                    neighbor_node->zVel = vel;
                }
            };
            break;
        default:
            throw "Incorrect type of axis";
    }

    for(int n = 0; n < fBoundary->nodesCount; ++n)
    {
        resetVelocity(fBoundary->nodes[n]);
	//int range = 3;
	//if (fBoundary->nodes[n]->type == ELASTIC_NODE)
	//{
	//    range = 2;
	//}

		int x_int = Index(grid, fBoundary->nodes[n]->x, COORD_X, axis, RANGE, RANGE);
		int y_int = Index(grid, fBoundary->nodes[n]->y, COORD_Y, axis, RANGE, RANGE);
		int z_int = Index(grid, fBoundary->nodes[n]->z, COORD_Z, axis, RANGE, RANGE);

        for(int i = x_int - RANGE; i <= x_int + RANGE; ++i)
        {
            for(int j = y_int - RANGE; j <= y_int + RANGE; ++j)
            {
                for(int k = z_int - RANGE; k <= z_int + RANGE; ++k)
                {

                    const double dist_x = fabs(fBoundary->nodes[n]->x - Coord(grid, i, COORD_X, axis));
                    const double dist_y = fabs(fBoundary->nodes[n]->y - Coord(grid, j, COORD_Y, axis));
                    const double dist_z = fabs(fBoundary->nodes[n]->z - Coord(grid, k, COORD_Z, axis));

                    const double weight_x = Dirichlet(grid, dist_x, COORD_X);
                    const double weight_y = Dirichlet(grid, dist_y, COORD_Y);
                    const double weight_z = Dirichlet(grid, dist_z, COORD_Z);

                    // interpolation from staggered grid before computation

                    // I don't know, maybe j-1 is needed (fictive shapes?)
                    //double velocity = velocity[i][j][k];

                    // 1 is a density
                    // this formulas are related to rigid IB method
                    //boundary->nodes[n].x_vel += (
                            //(
                             //velocity_U*CB(Hx[i]) + 0.5 * force_X[i][j][k] / 1.0
                             //) * weight_x * weight_y * weight_z);

                    //boundary->nodes[n].y_vel += (
                            //(
                             //velocity_V*CB(Hy[j]) + 0.5 * force_Y[i][j][k] / 1.0
                            //) * weight_x * weight_y * weight_z);
                    //boundary->nodes[n].z_vel += (
                            //(
                             //velocity_W*CB(Hz[k]) + 0.5 * force_Z[i][j][k] / 1.0
                            //) * weight_x * weight_y * weight_z);

		    //LOG(INFO) << "vel * wx * wy * wz = " << velocity[i][j][k] << " " << weight_x << " " << weight_y << " " << weight_z;
                    updateVelocity(
                            fBoundary->nodes[n],
                            (velocity[i][j][k] * weight_x * weight_y * weight_z)
                            * CB(separator->Separation[0])
                            );
                    //updateDebug(&fBoundary->nodes[n], velocity[i][j][k], i, j, k);
                }
            }
        }
    }

    return;
}

void TImmersedBoundaryProblem::InterpolateUpdatePathes(TRnRelease3dSpace &velocity, T3dNormalGrid *grid, const int axis, const int timeStepNumber)
{
	int div_timestep = 10;
	const TSeparator *separator = NULL;
	double **PCoord;

	switch (axis) {
	case COORD_X:
		separator = &grid->GetSeparator1();
		PCoord = fBoundary->X_pathes;
		break;
	case COORD_Y:
		separator = &grid->GetSeparator2();
		PCoord = fBoundary->Y_pathes;
		break;
	case COORD_Z:
		separator = &grid->GetSeparator3();
		PCoord = fBoundary->Z_pathes;
		break;
	default:
		throw "Incorrect type of axis";
	}

	int tm = timeStepNumber; //% fBoundary->max_path_step;

	//if (tm != 0)
	for (int n = 0; n < fBoundary->n_pathes; ++n)
	{
		double sum_ui = 0.0, sum_wi = 0.0, dist, wi, vel;
		int x_int = Index(grid, fBoundary->X_pathes[n][tm-1], COORD_X, axis, RANGE, RANGE);
		int y_int = Index(grid, fBoundary->Y_pathes[n][tm-1], COORD_Y, axis, RANGE, RANGE);
		int z_int = Index(grid, fBoundary->Z_pathes[n][tm-1], COORD_Z, axis, RANGE, RANGE);

		for (int i = x_int - RANGE; i <= x_int + RANGE; ++i)
		{
			for (int j = y_int - RANGE; j <= y_int + RANGE; ++j)
			{
				for (int k = z_int - RANGE; k <= z_int + RANGE; ++k)
				{

					const double dist_x = fabs(fBoundary->X_pathes[n][tm-1] - Coord(grid, i, COORD_X, axis));
					const double dist_y = fabs(fBoundary->Y_pathes[n][tm-1] - Coord(grid, j, COORD_Y, axis));
					const double dist_z = fabs(fBoundary->Z_pathes[n][tm-1] - Coord(grid, k, COORD_Z, axis));
					dist = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;

					if (dist < 2 * separator->Separation[0])
					{
						wi = 1 / dist;
						sum_wi += wi;
						sum_ui += velocity[i][j][k] * wi;
					}

					vel = sum_ui / sum_wi;
					PCoord[n][tm] = PCoord[n][tm-1] + vel * fTimeSeparator->SeparationValue;


				}
			}
		}
		if (n == 0) printf("tm=%d ax=%d vel=%lf coord=%lf, predcoord=%lf\n", tm, axis, vel, PCoord[n][tm], PCoord[n][tm-1]);
	}
	return;
}


void TImmersedBoundaryProblem::UpdateBoundaryPosition(const int axis, const int timeStepNumber)
{
    std::function<void(int)> update;
    switch(axis) {
        case COORD_X:
            update = [this](int n) {
                fBoundary->nodes[n]->xPrev = fBoundary->nodes[n]->x;
                fBoundary->nodes[n]->x += fBoundary->nodes[n]->xVel * fTimeSeparator->SeparationValue;
		//LOG(INFO) << "x pos = " << fBoundary->nodes[n]->x;
            };
            break;
        case COORD_Y:
            update = [this](int n) {
                fBoundary->nodes[n]->yPrev = fBoundary->nodes[n]->y;
                fBoundary->nodes[n]->y += fBoundary->nodes[n]->yVel * fTimeSeparator->SeparationValue;
		//LOG(INFO) << "y pos = " << fBoundary->nodes[n]->y;
            };
            break;
        case COORD_Z:
            update = [this](int n) {
                fBoundary->nodes[n]->zPrev = fBoundary->nodes[n]->z;
                fBoundary->nodes[n]->z += fBoundary->nodes[n]->zVel * fTimeSeparator->SeparationValue;
		//LOG(INFO) << "z pos = " << fBoundary->nodes[n]->z;
            };
            break;
        default:
            throw "Incorrect type of axis";
    }

    #pragma omp parallel for
    for(int n = 0; n < fBoundary->nodesCount; ++n)
    {
        update(n);
    }

    return;
}

int TImmersedBoundaryProblem::Index(T3dNormalGrid *grid, const double coord, const int type, const int range)
{	
    int index = 0;
    int N = 0;

    // TODO: replace with map and precomputed indexes
    switch(type) {
        case COORD_X:
            index = floor(coord / grid->GetSeparator1().Separation[0]);
            N = grid->GetSeparator1().EndIndex;
            break;
        case COORD_Y:
            index = floor(coord / grid->GetSeparator2().Separation[0]);
            N = grid->GetSeparator2().EndIndex;
            break;
        case COORD_Z:
            index = floor(coord / grid->GetSeparator3().Separation[0]);
            N = grid->GetSeparator3().EndIndex;
            break;
        default:
            throw "Incorrect type of axis";
    }


    if (index > (N - 1) - range) {
        index = (N - 1) - range;
    }
    if (index < range) {
        index = range;
    }


    return index;
}

int TImmersedBoundaryProblem::Index(T3dNormalGrid *grid, const double coord, const int type, const int axis, const int range_down, const int range_up)
{
	int index = 0;
	int N = 0;
	double step = 0;

	switch (type) {
	case COORD_X:
		step = grid->GetSeparator1().Separation[0];
		N = grid->GetSeparator1().EndIndex;
		break;
	case COORD_Y:
		step = grid->GetSeparator2().Separation[0];
		N = grid->GetSeparator2().EndIndex;
		break;
	case COORD_Z:
		step = grid->GetSeparator3().Separation[0];
		N = grid->GetSeparator3().EndIndex;
		break;
	default:
		throw "Incorrect type of axis";
	}

	index = floor((coord + (type != axis ? 0.5 * step : 0)) / step);

	if (index > (N - 1) - range_up)
	{
		index = (N - 1) - range_up;
	}

	if (index < range_down)
	{
		index = range_down;
	}

	return index;
}

double TImmersedBoundaryProblem::Coord(T3dNormalGrid *grid, const int index, const int type)
{
    double step;
    if (type == COORD_X) {
        step = grid->GetSeparator1().Separation[0];
    }
    else if (type == COORD_Y) {
	step = grid->GetSeparator2().Separation[0];
    }
    else if (type == COORD_Z) {
	step = grid->GetSeparator3().Separation[0];
    }
    else {
	throw "Incorrect type of axis";
    }

    return index * step;
}

double TImmersedBoundaryProblem::Coord(T3dNormalGrid *grid, const int index, const int type, const int data_axis)
{
    double step;
    if (type == COORD_X) {
        step = grid->GetSeparator1().Separation[0];
    }
    else if (type == COORD_Y) {
	    step = grid->GetSeparator2().Separation[0];
    }
    else if (type == COORD_Z) {
	    step = grid->GetSeparator3().Separation[0];
    }
    else {
	throw "Incorrect type of axis";
    }

    //return index * step + (type == data_axis ? 0.5 * step : 0);
	return index * step - (type != data_axis ? 0.5 * step : 0);
}

double TImmersedBoundaryProblem::OptimalDirichlet(double distance, const double h)
{
    // NOTE: this algorithm requires to remove * h^3 from force computation
    // and velocity interpolation.
    //
    //distance = distance / h;
    //if(distance >= 1.0 && distance < 2.0)
    //{
    //    return 1.0/8.0 * (5.0 - 2.0 * distance - sqrt(-7.0 + 12.0 * distance - 4.0 * distance * distance));
    //}
    //else if(distance >=0 && distance < 1.0)
    //{
    //    return 1.0/8.0 * (3.0 - 2.0 * distance + sqrt(1.0 + 4.0 * distance - 4.0 * distance * distance));
    //}
    //else
    //{
    //    return 0;
    //}
    if(distance < 2*h) {
        return 1.0/(4 * h) * (1 + cos(M_PI * distance / (2 * h)));
    }
    else {
        return 0;
    }
}

double TImmersedBoundaryProblem::Dirichlet(T3dNormalGrid *grid, const double distance, const int type)
{
    switch(type) {
        case COORD_X:
            return OptimalDirichlet(distance, grid->GetSeparator1().Separation[0]);
        case COORD_Y:
            return OptimalDirichlet(distance, grid->GetSeparator2().Separation[0]);
        case COORD_Z:
            return OptimalDirichlet(distance, grid->GetSeparator3().Separation[0]);
        default:
            throw "Incorrect type of axis";
    }
}

void TImmersedBoundaryProblem::OutputBoundary(int timeStepNumber)
{
	//ComputeStress(timeStepNumber);
	double *boundary_x = new double[fBoundary->nodesCount];
	double *boundary_y = new double[fBoundary->nodesCount];
	double *boundary_z = new double[fBoundary->nodesCount];

	double *boundary_x_next = new double[fBoundary->nodesCount];
	double *boundary_y_next = new double[fBoundary->nodesCount];
	double *boundary_z_next = new double[fBoundary->nodesCount];

	double *boundary_x_prev = new double[fBoundary->nodesCount];
	double *boundary_y_prev = new double[fBoundary->nodesCount];
	double *boundary_z_prev = new double[fBoundary->nodesCount];

	double *boundary_force_x = new double[fBoundary->nodesCount];
	double *boundary_force_y = new double[fBoundary->nodesCount];
	double *boundary_force_z = new double[fBoundary->nodesCount];

	double *boundary_vel_x = new double[fBoundary->nodesCount];
	double *boundary_vel_y = new double[fBoundary->nodesCount];
	double *boundary_vel_z = new double[fBoundary->nodesCount];

	vector<double> surround_x;
	vector<double> surround_y;
	vector<double> surround_z;

	vector<double> surround_x_vel;
	vector<double> surround_y_vel;
	vector<double> surround_z_vel;

	for (int n = 0; n < fBoundary->nodesCount; ++n) {
		boundary_x[n] = fBoundary->nodes[n]->x;
		boundary_y[n] = fBoundary->nodes[n]->y;
		boundary_z[n] = fBoundary->nodes[n]->z;

		if (fBoundary->nodes[n]->neighbors.next != NULL)
		{
			boundary_x_next[n] = fBoundary->nodes[n]->neighbors.next->x;
			boundary_y_next[n] = fBoundary->nodes[n]->neighbors.next->y;
			boundary_z_next[n] = fBoundary->nodes[n]->neighbors.next->z;
		}
		else
		{
			boundary_x_next[n] = 0;
			boundary_y_next[n] = 0;
			boundary_z_next[n] = 0;
		}

		if (fBoundary->nodes[n]->neighbors.prev != NULL)
		{
			boundary_x_prev[n] = fBoundary->nodes[n]->neighbors.prev->x;
			boundary_y_prev[n] = fBoundary->nodes[n]->neighbors.prev->y;
			boundary_z_prev[n] = fBoundary->nodes[n]->neighbors.prev->z;
		}
		else
		{
			boundary_x_prev[n] = 0;
			boundary_y_prev[n] = 0;
			boundary_z_prev[n] = 0;
		}

		boundary_force_x[n] = fBoundary->nodes[n]->xForce;
		boundary_force_y[n] = fBoundary->nodes[n]->yForce;
		boundary_force_z[n] = fBoundary->nodes[n]->zForce;

		boundary_vel_x[n] = fBoundary->nodes[n]->xVel;
		boundary_vel_y[n] = fBoundary->nodes[n]->yVel;
		boundary_vel_z[n] = fBoundary->nodes[n]->zVel;

		for (auto neighbor_node : fBoundary->nodes[n]->surround)
		{
			surround_x.push_back(neighbor_node.second->x);
			surround_y.push_back(neighbor_node.second->y);
			surround_z.push_back(neighbor_node.second->z);

			surround_x_vel.push_back(neighbor_node.second->xVel);
			surround_y_vel.push_back(neighbor_node.second->yVel);
			surround_z_vel.push_back(neighbor_node.second->zVel);
		}
	}

	//-------------------------------------------------------------------
	// Two formats output
	//-------------------------------------------------------------------

	if (timeStepNumber % TProblem::GlobalSaveStep == 0 || timeStepNumber == 1)
	{
		int len = strlen(TProblem::GlobalCatalog);
		int i1, i2, i3, i4;
		int max_cur, max_next;
		int Circles = fBoundary->height_nodes;
		int InCircle = fBoundary->radius_nodes;
		FILE* f = NULL;
		char *zonesName;
		char fileName[40];
		/*
		zonesName = new char[len + 50];
		string_copy(zonesName, TProblem::GlobalCatalog, len + 50);

		
		string_print(fileName, 40, "bound%d.dat", 10000 + timeStepNumber);
		string_concat(zonesName, fileName, len + 50);

		f = file_open(zonesName, "w");

		printf("save Circles = %d inCircle = %d\n", Circles, InCircle);
		fprintf(f, "TITLE = \"CVP\"\n");
		//fprintf(f,"VARIABLES = X,Y,Z,U,V,W,P\n");
		fprintf(f, "VARIABLES = \"X\",\"Y\",\"Z\",\"U\",\"V\",\"W\",\"P\",\"Stress\",\"Mu\",\"Concentration\",\"Stiff\"\n");
		//fprintf(f, "ZONE I=%d, J=%d, K=%d, F=POINT\n", N + 1, L + 1, M + 1);
		fprintf(f, "ZONE N =%d, E =%d, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL\n", Circles*InCircle, InCircle*(Circles - 1));

		for (int cr = 0; cr < Circles; cr++)
			for (int i = 0; i < InCircle; i++)
			{
				fprintf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", (double)boundary_x[cr*InCircle + i], (double)boundary_y[cr*InCircle + i], (double)boundary_z[cr*InCircle + i],
					fBoundary->nodes[cr*InCircle + i]->xVel, fBoundary->nodes[cr*InCircle + i]->yVel, fBoundary->nodes[cr*InCircle + i]->zVel,
					fBoundary->nodes[cr*InCircle + i]->p,fBoundary->nodes[cr*InCircle + i]->ShearStress,
					fBoundary->nodes[cr*InCircle + i]->vis, fBoundary->nodes[cr*InCircle + i]->concentration, fBoundary->nodes[cr*InCircle + i]->stretchingStiffness);
			}


		for (int cr = 0; cr < Circles - 1; ++cr)
			for (int i = 0; i < InCircle; ++i)
			{
				max_cur = (cr + 1)*InCircle;
				max_next = (cr + 2)*InCircle;
				i1 = cr*InCircle + i + 1;
				i2 = i1 + InCircle;
				i3 = i2 + 1;
				if (i3 > max_next)
					i3 -= InCircle;
				i4 = i1 + 1;
				if (i4 > max_cur)
					i4 -= InCircle;
				fprintf(f, "%d %d %d %d\n", i1, i2, i3, i4);
			}
		fclose(f);
		//*/
		//---------------------------------------
		// vtk
		//---------------------------------------
		//len = strlen(TProblem::GlobalCatalog);
		char *vtk_fileName = new char[len + 50];
		string_copy(vtk_fileName, TProblem::GlobalCatalog, len + 50);

		fileName[0] = '\0';

		if (timeStepNumber == 1)
			string_print(fileName, 40, "bound%d.vtk", 10000);
		else
			string_print(fileName, 40, "bound%d.vtk", 10000 + timeStepNumber);
		string_concat(vtk_fileName, fileName, len + 50);

		f = file_open(vtk_fileName, "w");

		fprintf(f, "# vtk DataFile Version 1.0\nBound\nASCII\n\n");
		fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
		fprintf(f, "POINTS %d double\n", Circles*InCircle);

		for (int cr = 0; cr < Circles; ++cr)
			for (int i = 0; i < InCircle; ++i)
			{
				fprintf(f, "%lf %lf %lf\n", (double)boundary_x[cr*InCircle + i], (double)boundary_y[cr*InCircle + i], (double)boundary_z[cr*InCircle + i]);
			}

		fprintf(f, "CELLS %d %d\n", InCircle*(Circles - 1), 5 * InCircle*(Circles - 1)); // 6900 34500
		for (int cr = 0; cr < Circles - 1; ++cr)
			for (int i = 0; i < InCircle; ++i)
			{
				//max_cur = (cr + 1)*InCircle;
				//max_next = (cr + 2)*InCircle;
				i1 = cr*InCircle + i;
				i2 = i1 + 1;
				i3 = i2 + InCircle;
				i4 = i1 + InCircle;
				if (i == InCircle - 1)
				{
					i2 = cr*InCircle;
					i3 = (cr + 1)*InCircle;
				}
				fprintf(f, "%d %d %d %d %d\n", 4, i1, i2, i3, i4);
			}

		fprintf(f, "CELL_TYPES %d\n", InCircle*(Circles - 1)); //  6900
		for (int cr = 0; cr < Circles - 1; ++cr)
			for (int i = 0; i < InCircle; ++i)
				fprintf(f, "%d\n", 9);

		fprintf(f, "POINT_DATA %d\n", InCircle*Circles);//6969
		fprintf(f, "VECTORS Velocity double\n");
		for (int cr = 0; cr < Circles; ++cr)
			for (int i = 0; i < InCircle; ++i)
			{
				fprintf(f, "%lf %lf %lf\n", fBoundary->nodes[cr*InCircle + i]->xVel, fBoundary->nodes[cr*InCircle + i]->yVel,
					fBoundary->nodes[cr*InCircle + i]->zVel);
			}

		fprintf(f, "\n\nSCALARS Pressure float\nLOOKUP_TABLE default\n");
		for (int cr = 0; cr < Circles; ++cr)
			for (int i = 0; i < InCircle; ++i)
			{
				fprintf(f, "%lf\n", fBoundary->nodes[cr*InCircle + i]->p);
			}

		fprintf(f, "\n\nSCALARS Mu float\nLOOKUP_TABLE default\n");
		for (int cr = 0; cr < Circles; ++cr)
			for (int i = 0; i < InCircle; ++i)
			{
				fprintf(f, "%lf\n", fBoundary->nodes[cr*InCircle + i]->vis); // fBoundary->nodes[cr*InCircle + i]->rho);
			}

		fprintf(f, "\n\nSCALARS Concentration float\nLOOKUP_TABLE default\n");
		for (int cr = 0; cr < Circles; ++cr)
			for (int i = 0; i < InCircle; ++i)
			{
				fprintf(f, "%lf\n", fBoundary->nodes[cr*InCircle + i]->concentration); // fBoundary->nodes[cr*InCircle + i]->cc);
			}
		/*fprintf(f, "\n\nSCALARS ShearStress float\nLOOKUP_TABLE default\n");
		for (int cr = 0; cr < Circles; ++cr)
			for (int i = 0; i < InCircle; ++i)
			{
				fprintf(f, "%lf\n", fBoundary->nodes[cr*InCircle + i]->ShearStress); // fBoundary->nodes[cr*InCircle + i]->cc);
			}*/
		//fprintf(f,"VARIABLES = X,Y,Z,U,V,W,P\n");
		//fprintf(f, "VARIABLES = X,Y,Z,U,V,W,P,Rho,CC\n");
		fclose(f);

		delete[]vtk_fileName;// vtk_fileName = NULL;
		delete[]zonesName;// zonesName = NULL;


		
	}

//-------------------------------------------------------------------
// Two formats output ENDs
//-------------------------------------------------------------------
	/*
	if (timeStepNumber % fBoundary->out_path_step == 0)
	{
		int len = strlen(TProblem::GlobalCatalog);
		char *zonsName = new char[len + 50];

		char fileName[40];
		fileName[0] = '\0';

		string_copy(zonsName, TProblem::GlobalCatalog, len + 50);
		string_print(fileName, 40, "bound%d.dat", 10000 + timeStepNumber);
		string_concat(zonsName, fileName, len + 50);

		FILE* f = file_open(zonsName, "w");


		int Circles = fBoundary->height_nodes;
		int InCircle = fBoundary->radius_nodes;
		fprintf(f, "TITLE = \"BOUNDARY AND PATHES\"\n");
		fprintf(f, "VARIABLES = \"X\",\"Y\",\"Z\",\"U\",\"V\",\"W\",\"P\",\"Rho\",\"CC\"\n");
		fprintf(f, "ZONE N =%d, E =%d, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL\n", Circles*InCircle, InCircle*(Circles - 1));

		for (int cr = 0; cr<Circles; ++cr)
		for (int i = 0; i<InCircle; ++i)
		{
		fprintf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", (double)fBoundary->nodes[cr*InCircle + i]->x,
		(double)fBoundary->nodes[cr*InCircle + i]->y, (double)fBoundary->nodes[cr*InCircle + i]->z,
		fBoundary->nodes[cr*InCircle + i]->xVel, fBoundary->nodes[cr*InCircle + i]->yVel, fBoundary->nodes[cr*InCircle + i]->zVel,
		fBoundary->nodes[cr*InCircle + i]->p, 0., 0.);
		}


		int i1, i2, i3, i4;
		int max_cur, max_next;
		for (int cr = 0; cr<Circles - 1; ++cr)
		for (int i = 0; i<InCircle; ++i)
		{
			i1 = cr*InCircle + i;
			i2 = i1 + 1;
			i3 = i2 + InCircle;
			i4 = i1 + InCircle;
			if (i == InCircle - 1)
			{
				i2 = cr*InCircle;
				i3 = (cr + 1)*InCircle;
			}
		    fprintf(f, "%d %d %d %d\n", i1+1, i2+1, i3+1, i4+1);
		}

		for (int ln = 0; ln < fBoundary->n_pathes; ln++)
		{
			fprintf(f, "GEOMETRY X = 0, Y = 0, Z = 0, M = GRID, C = BLACK, T = LINE3D, F = POINT, AAT = BOTH,\n");
			fprintf(f, "AST = PLAIN, ASZ = 0.1, LT=0.3\n");

			fprintf(f, "1\n");
			fprintf(f, "%d\n", timeStepNumber);

			for (int i = 0; i <timeStepNumber; i++)
				fprintf(f, "%lf %lf %lf\n", fBoundary->X_pathes[ln][i], fBoundary->Y_pathes[ln][i], fBoundary->Z_pathes[ln][i]);
		}

		fclose(f);
		delete[]zonsName; zonsName = NULL;
	}
	*/
#if CNPY_VISUALIZATION
    //char fileName[40];

    //double *boundary_x = new double[fBoundary->nodesCount];
    //double *boundary_y = new double[fBoundary->nodesCount];
    //double *boundary_z = new double[fBoundary->nodesCount];

    //double *boundary_x_next = new double[fBoundary->nodesCount];
    //double *boundary_y_next = new double[fBoundary->nodesCount];
    //double *boundary_z_next = new double[fBoundary->nodesCount];

    //double *boundary_x_prev = new double[fBoundary->nodesCount];
    //double *boundary_y_prev = new double[fBoundary->nodesCount];
    //double *boundary_z_prev = new double[fBoundary->nodesCount];

    //double *boundary_force_x = new double[fBoundary->nodesCount];
    //double *boundary_force_y = new double[fBoundary->nodesCount];
    //double *boundary_force_z = new double[fBoundary->nodesCount];

    //double *boundary_vel_x = new double[fBoundary->nodesCount];
    //double *boundary_vel_y = new double[fBoundary->nodesCount];
    //double *boundary_vel_z = new double[fBoundary->nodesCount];

    //vector<double> surround_x;
    //vector<double> surround_y;
    //vector<double> surround_z;

    //vector<double> surround_x_vel;
    //vector<double> surround_y_vel;
    //vector<double> surround_z_vel;

    for(int n = 0; n < fBoundary->nodesCount; ++n) {
        boundary_x[n] = fBoundary->nodes[n]->x;
        boundary_y[n] = fBoundary->nodes[n]->y;
        boundary_z[n] = fBoundary->nodes[n]->z;

        if (fBoundary->nodes[n]->neighbors.next != NULL)
        {
            boundary_x_next[n] = fBoundary->nodes[n]->neighbors.next->x;
            boundary_y_next[n] = fBoundary->nodes[n]->neighbors.next->y;
            boundary_z_next[n] = fBoundary->nodes[n]->neighbors.next->z;
        }
        else
        {
            boundary_x_next[n] = 0;
            boundary_y_next[n] = 0;
            boundary_z_next[n] = 0;
        }

        if (fBoundary->nodes[n]->neighbors.prev != NULL)
        {
            boundary_x_prev[n] = fBoundary->nodes[n]->neighbors.prev->x;
            boundary_y_prev[n] = fBoundary->nodes[n]->neighbors.prev->y;
            boundary_z_prev[n] = fBoundary->nodes[n]->neighbors.prev->z;
        }
        else
        {
            boundary_x_prev[n] = 0;
            boundary_y_prev[n] = 0;
            boundary_z_prev[n] = 0;
        }

        boundary_force_x[n] = fBoundary->nodes[n]->xForce;
        boundary_force_y[n] = fBoundary->nodes[n]->yForce;
        boundary_force_z[n] = fBoundary->nodes[n]->zForce;

        boundary_vel_x[n] = fBoundary->nodes[n]->xVel;
        boundary_vel_y[n] = fBoundary->nodes[n]->yVel;
        boundary_vel_z[n] = fBoundary->nodes[n]->zVel;

        for(auto neighbor_node : fBoundary->nodes[n]->surround)
        {
            surround_x.push_back(neighbor_node.second->x);
            surround_y.push_back(neighbor_node.second->y);
            surround_z.push_back(neighbor_node.second->z);

            surround_x_vel.push_back(neighbor_node.second->xVel);
            surround_y_vel.push_back(neighbor_node.second->yVel);
            surround_z_vel.push_back(neighbor_node.second->zVel);
        }
    }

    const unsigned int shape_boundary[] = {fBoundary->nodesCount};
	string_print(fileName, 40, "boundary_x_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, boundary_x, shape_boundary, 1, "w");

	string_print(fileName, 40, "boundary_y_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, boundary_y, shape_boundary, 1, "w");

	string_print(fileName, 40, "boundary_z_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, boundary_z, shape_boundary, 1, "w");


	string_print(fileName, 40, "boundary_x_next_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, boundary_x_next, shape_boundary, 1, "w");

	string_print(fileName, 40, "boundary_y_next_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, boundary_y_next, shape_boundary, 1, "w");

	string_print(fileName, 40, "boundary_z_next_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, boundary_z_next, shape_boundary, 1, "w");

	string_print(fileName, 40, "boundary_x_prev_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, boundary_x_prev, shape_boundary, 1, "w");

	string_print(fileName, 40, "boundary_y_prev_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, boundary_y_prev, shape_boundary, 1, "w");

	string_print(fileName, 40, "boundary_z_prev_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, boundary_z_prev, shape_boundary, 1, "w");


	string_print(fileName, 40, "boundary_force_x_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, boundary_force_x, shape_boundary, 1, "w");

	string_print(fileName, 40, "boundary_force_y_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, boundary_force_y, shape_boundary, 1, "w");

	string_print(fileName, 40, "boundary_force_z_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, boundary_force_z, shape_boundary, 1, "w");

	string_print(fileName, 40, "boundary_vel_x_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, boundary_vel_x, shape_boundary, 1, "w");

	string_print(fileName, 40, "boundary_vel_y_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, boundary_vel_y, shape_boundary, 1, "w");

	string_print(fileName, 40, "boundary_vel_z_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, boundary_vel_z, shape_boundary, 1, "w");

    const unsigned int surround_shape[] = {(const unsigned int)surround_x.size()};
	string_print(fileName, 40, "surround_x_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, (double *)&surround_x[0], surround_shape, 1, "w");

	string_print(fileName, 40, "surround_y_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, (double *)&surround_y[0], surround_shape, 1, "w");

	string_print(fileName, 40, "surround_z_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, (double *)&surround_z[0], surround_shape, 1, "w");

	string_print(fileName, 40, "surround_x_vel_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, (double *)&surround_x_vel[0], surround_shape, 1, "w");

	string_print(fileName, 40, "surround_y_vel_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, (double *)&surround_y_vel[0], surround_shape, 1, "w");

	string_print(fileName, 40, "surround_z_vel_%03d.npy", timeStepNumber);
    cnpy::npy_save(fileName, (double *)&surround_z_vel[0], surround_shape, 1, "w");

#endif

	delete[]boundary_x;// boundary_x = NULL;
	delete[]boundary_y;// boundary_y = NULL;
	delete[]boundary_z;// boundary_z = NULL;

	delete[]boundary_x_next;// boundary_x_next = NULL;
	delete[]boundary_y_next;// boundary_y_next = NULL;
	delete[]boundary_z_next;// boundary_z_next = NULL;

	delete[]boundary_x_prev;// boundary_x_prev = NULL;
	delete[]boundary_y_prev;// boundary_y_prev = NULL;
	delete[]boundary_z_prev;// boundary_z_prev = NULL;

	delete[]boundary_force_x;// boundary_force_x = NULL;
	delete[]boundary_force_y;// boundary_force_y = NULL;
	delete[]boundary_force_z;// boundary_force_z = NULL;

	delete[]boundary_vel_x;// boundary_vel_x = NULL;
	delete[]boundary_vel_y;// boundary_vel_y = NULL;
	delete[]boundary_vel_z;// boundary_vel_z = NULL;

}

void TImmersedBoundaryProblem::ComputeChanging(int timeStepNumber)
{
	double *boundary_x = new double[fBoundary->nodesCount];
	double *boundary_y = new double[fBoundary->nodesCount];
	double *boundary_z = new double[fBoundary->nodesCount];

	double *boundary_x_next = new double[fBoundary->nodesCount];
	double *boundary_y_next = new double[fBoundary->nodesCount];
	double *boundary_z_next = new double[fBoundary->nodesCount];

	double *boundary_x_prev = new double[fBoundary->nodesCount];
	double *boundary_y_prev = new double[fBoundary->nodesCount];
	double *boundary_z_prev = new double[fBoundary->nodesCount];

	for (int n = 0; n < fBoundary->nodesCount; ++n) {
		boundary_x[n] = fBoundary->nodes[n]->x;
		boundary_y[n] = fBoundary->nodes[n]->y;
		boundary_z[n] = fBoundary->nodes[n]->z;

		if (fBoundary->nodes[n]->neighbors.next != NULL)
		{
			boundary_x_next[n] = fBoundary->nodes[n]->neighbors.next->x;
			boundary_y_next[n] = fBoundary->nodes[n]->neighbors.next->y;
			boundary_z_next[n] = fBoundary->nodes[n]->neighbors.next->z;
		}
		else
		{
			boundary_x_next[n] = 0;
			boundary_y_next[n] = 0;
			boundary_z_next[n] = 0;
		}

		if (fBoundary->nodes[n]->neighbors.prev != NULL)
		{
			boundary_x_prev[n] = fBoundary->nodes[n]->neighbors.prev->x;
			boundary_y_prev[n] = fBoundary->nodes[n]->neighbors.prev->y;
			boundary_z_prev[n] = fBoundary->nodes[n]->neighbors.prev->z;
		}
		else
		{
			boundary_x_prev[n] = 0;
			boundary_y_prev[n] = 0;
			boundary_z_prev[n] = 0;
		}
	}

	//-------------------------------------------------------------------
	// test
	//-------------------------------------------------------------------
	{

		int Circles = fBoundary->height_nodes;
		int InCircle = fBoundary->radius_nodes;

		double *slice_x = new double[Circles+1];
		double *slice_y = new double[Circles+1];
		double *slice_z = new double[Circles+1];

		if (timeStepNumber == 1)
		{
			slice_x_0 = new double[Circles+1];
			slice_y_0 = new double[Circles+1];
			slice_z_0 = new double[Circles+1];
		}

		for (int cr = 0; cr <= Circles; cr++)
		{
			if (timeStepNumber == 1)
			{
				slice_x_0[cr] = boundary_x[cr*InCircle + (InCircle / 2)];
				slice_y_0[cr] = boundary_y[cr*InCircle + (InCircle / 2)];
				slice_z_0[cr] = boundary_z[cr*InCircle + (InCircle / 2)];
			}
			slice_x[cr] = boundary_x[cr*InCircle + (InCircle / 2)];
			slice_y[cr] = boundary_y[cr*InCircle + (InCircle / 2)];
			slice_z[cr] = boundary_z[cr*InCircle + (InCircle / 2)];
		}

		

		double gt_x = 0;
		double gt_y = 0;
		double gt_z = 0;

		double gt_x_1 = 0;
		double gt_y_1 = 0;
		double gt_z_1 = 0;

		int len = strlen(TProblem::GlobalCatalog);
		char *vol_fileName = new char[len + 50];
		char fileName[40];
		string_copy(vol_fileName, TProblem::GlobalCatalog, len + 50);
		string_print(fileName, 40, "Vol_1.txt");
		string_concat(vol_fileName, fileName, len + 50);

		FILE* t1;
		if (timeStepNumber == 1)
			t1 = file_open(vol_fileName, "w");
		else
			t1 = file_open(vol_fileName, "a");

//		printf("%lf\n",slice_x_0[Circles-1]);

		for (int i = 0; i <= Circles - 1; i++)
		{
			gt_x_1 += ( (slice_x[i] - TImmersedBoundaryProblem::slice_x_0[i]) + (slice_x[i + 1] - TImmersedBoundaryProblem::slice_x_0[i + 1]) ) * (fabs(slice_x[i + 1] - slice_x[i])) / 2;
			gt_y_1 += ( (slice_y[i] - TImmersedBoundaryProblem::slice_y_0[i]) + (slice_y[i + 1] - TImmersedBoundaryProblem::slice_y_0[i + 1]) ) * (fabs(slice_x[i + 1] - slice_x[i])) / 2;
			gt_z_1 += ( (slice_z[i] - TImmersedBoundaryProblem::slice_z_0[i]) + (slice_z[i + 1] - TImmersedBoundaryProblem::slice_z_0[i + 1]) ) * (fabs(slice_x[i + 1] - slice_x[i])) / 2;
		}
		//printf("GT_x_1 -> %15.10lf \nGT_y_1 -> %15.10lf \nGT_z_1 -> %15.10lf \n\n", gt_x_1, gt_y_1, gt_z_1);
		fprintf(t1, "%15.10lf %15.10lf %15.10lf\n", gt_x_1, gt_y_1, gt_z_1);
		fclose(t1);

		//*
		char *vol_x_fileName = new char[len + 50];
		char *vol_y_fileName = new char[len + 50];
		char *vol_z_fileName = new char[len + 50];

		char fnx[40];
		char fny[40];
		char fnz[40];

		string_copy(vol_x_fileName, TProblem::GlobalCatalog, len + 50);
		string_copy(vol_y_fileName, TProblem::GlobalCatalog, len + 50);
		string_copy(vol_z_fileName, TProblem::GlobalCatalog, len + 50);

		string_print(fnx, 40, "Vol_x.dat");
		string_print(fny, 40, "Vol_y.dat");
		string_print(fnz, 40, "Vol_z.dat");

		string_concat(vol_x_fileName, fnx, len + 50);
		string_concat(vol_y_fileName, fny, len + 50);
		string_concat(vol_z_fileName, fnz, len + 50);

		FILE* tech_x;
		FILE* tech_y;
		FILE* tech_z;
		if (timeStepNumber == 1)
		{
			tech_x = file_open(vol_x_fileName, "w");
			tech_y = file_open(vol_y_fileName, "w");
			tech_z = file_open(vol_z_fileName, "w");
		}
		else
		{
			tech_x = file_open(vol_x_fileName, "a");
			tech_y = file_open(vol_y_fileName, "a");
			tech_z = file_open(vol_z_fileName, "a");
		}
		if (timeStepNumber == 1)
		{
			fprintf(tech_x, "TITLE=\"BoundX\"\nVARIABLES = \"coord\", \"X\"\nZONE T = \"Test\", I = %d, F = POINT\n",Circles);
			fprintf(tech_y, "TITLE=\"BoundY\"\nVARIABLES = \"coord\", \"Y\"\nZONE T = \"Test\", I = %d, F = POINT\n",Circles);
			fprintf(tech_z, "TITLE=\"BoundZ\"\nVARIABLES = \"coord\", \"Z\"\nZONE T = \"Test\", I = %d, F = POINT\n",Circles);
		}
		else
		{
			fprintf(tech_x, "ZONE T = \"Test\", I = %d, F = POINT\n", Circles+1);
			fprintf(tech_y, "ZONE T = \"Test\", I = %d, F = POINT\n", Circles+1);
			fprintf(tech_z, "ZONE T = \"Test\", I = %d, F = POINT\n",Circles+1);
		}
		for (int i = 0; i < Circles; i++)
		{
			fprintf(tech_x, "%lf %lf\n", slice_x_0[i], slice_x[i]);
			fprintf(tech_y, "%lf %lf\n", slice_x_0[i], slice_y[i]);
			fprintf(tech_z, "%lf %lf\n", slice_x_0[i], slice_z[i]);
		}
		fclose(tech_x);
		fclose(tech_y);
		fclose(tech_z);
		//*/
		delete[]vol_fileName;// vol_fileName = NULL;
		delete[]vol_x_fileName;// vol_x_fileName = NULL;
		delete[]vol_y_fileName;// vol_y_fileName = NULL;
		delete[]vol_z_fileName;// vol_z_fileName = NULL;
	}

	delete[]boundary_x;// boundary_x = NULL;
	delete[]boundary_y;// boundary_y = NULL;
	delete[]boundary_z;// boundary_z = NULL;

	delete[]boundary_x_next;// boundary_x_next = NULL;
	delete[]boundary_y_next;// boundary_y_next = NULL;
	delete[]boundary_z_next;// boundary_z_next = NULL;

	delete[]boundary_x_prev;// boundary_x_prev = NULL;
	delete[]boundary_y_prev;// boundary_y_prev = NULL;
	delete[]boundary_z_prev;// boundary_z_prev = NULL;
}

void TImmersedBoundaryProblem::ComputeStress(TRnRelease3dSpace &U, TRnRelease3dSpace &V, TRnRelease3dSpace &W, int timestepNumber)
{
	printf("Calculating Shear Stress for each imm. boundary node...\n");

	int grid_marg = 2;
	int grid_rad = 3;
	double h = (fStatement->Hx[1] + fStatement->Hy[1] + fStatement->Hz[1]) / 3.;
	double hsqr = h*h;
	double ub = 0, wi = 0, xb = 0, yb = 0, zb = 0,
		xu = 0, yu = 0, zu = 0,
		xv = 0, yv = 0, zv = 0,
		xw = 0, yw = 0, zw = 0,
		xc = 0, yc = 0, zc = 0,
		du = 0, dv = 0, dw = 0, sum_w = 0, dist = 0;
	double Q = 300;
	int Circles = fBoundary->height_nodes;
	int InCircle = fBoundary->radius_nodes;
	//-----------------------------------------------------
	// Counting stress in imm. boundary node
	//-----------------------------------------------------


	for (int il = 0; il < fBoundary->height_nodes; il++)
	{
		for (int ir = 0; ir < fBoundary->radius_nodes; ir++)
		{
			xb = fBoundary->nodes[il*fBoundary->radius_nodes + ir]->x;
			yb = fBoundary->nodes[il*fBoundary->radius_nodes + ir]->y;
			zb = fBoundary->nodes[il*fBoundary->radius_nodes + ir]->z;

			xc = fBoundary->nodes[il*fBoundary->radius_nodes + ir]->x;
			//yc = TCanalIBMWithElasticBoundary->fLengthY/2.;
			//zc = TCanalIBMWithElasticBoundary->fLengthZ/2.;
			yc = fBoundary->GetYCenter();
			zc =  fBoundary->GetZCenter();
			double mu = fBoundary->nodes[il*fBoundary->radius_nodes + ir]->vis;
			double r = sqrt( (xb - xc)*(xb - xc) + (yb - yc)*(yb - yc) + (zb - zc)*(zb - zc) );
			//printf("\n xb > %lf -- yb > %lf -- zb > %lf\n xc > %lf -- yc > %lf -- zc > %lf\nr > %lf\n",xb,yb,zb,xc,yc,zc,r);
			fBoundary->nodes[il*fBoundary->radius_nodes + ir]->ShearStress = (4*mu*Q)/(M_PI*r*r*r);;
		}
	}

	//----------------------------------
	//Debug! Output boundary as 2D plain
	//----------------------------------
	if(0 == 1)
	{
		int len = strlen(TProblem::GlobalCatalog);
		char *zonesName = new char[len + 50];
		string_copy(zonesName, TProblem::GlobalCatalog, len + 50);

		char fileName[40];
		string_print(fileName, 40, "debbound%d.dat", 10000 + timestepNumber);
		string_concat(zonesName, fileName, len + 50);

		FILE* f = file_open(zonesName, "w");

		fprintf(f, "TITLE = \"debbound2d\"\n");
		fprintf(f, "VARIABLES = \"X\",\"Y\",\"U\",\"V\",\"W\",\"P\",\"Stiff\",\"Stress\"\n");
		fprintf(f, "ZONE T=\"Deb%d\", I=%d, J=%d, F=POINT\n", timestepNumber, Circles, InCircle);
		//fprintf(f, "ZONE N =%d, E =%d, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL\n", Circles*InCircle, InCircle*(Circles - 1));

		for (int i = 0; i < InCircle; ++i)
			for (int cr = 0; cr < Circles; ++cr)
			{
				fprintf(f, "%lf %lf %lf %lf %lf %lf %lf %lf\n", cr*fBoundary->step, i*fBoundary->step,
					fBoundary->nodes[cr*InCircle + i]->xVel, fBoundary->nodes[cr*InCircle + i]->yVel, fBoundary->nodes[cr*InCircle + i]->zVel,
					fBoundary->nodes[cr*InCircle + i]->p, fBoundary->nodes[cr*InCircle + i]->stretchingStiffness, fBoundary->nodes[cr*InCircle + i]->ShearStress);
			}

		fclose(f);
		delete[]zonesName;
	}
	printf("Calculating Shear Stress for each imm. boundary node. OK.\n");
}
//********************************************************************************

TMaskGenerator* TCanalUnifyProblemWithImmersedBoundary::CreateStatement()
{
    //Preparing (Begin)
	double iPart;
 	int endIndex;

    iPart = 0;
	modf(fLength/fStep,&iPart);
    endIndex = (int)iPart;
	int nx = endIndex + 1;

	iPart = 0;
	modf(fRadius*2/fStep,&iPart);
    endIndex = (int)iPart;
	int ny = endIndex + 1;
	int nz = endIndex + 1;
	//Preparing (End)

	//Creating (Begin)
	return TSimpleWithPressureMaskGenerator::CreateInstance(
		                                                    fLength,fRadius*2,fRadius*2,
		                                                    nx,ny,nz,
															fPressureInput,fPressureOutput
														   );
	//Creating (End)
}


TSimpleImmersedValveProblem::~TSimpleImmersedValveProblem()
{
    delete fBoundary;
}

void TSimpleImmersedValveProblem::Initialize()
{

    TVPProblem::Initialize();
    fBoundary = (TSimpleValve *)(new TSimpleValve())
        ->withStiffness(fStiffness)
        ->withNodesCount(200)
        ->withRadius(0.1)
        ->withCenter(0.5)
        ->withHeight(0.1)
        ->Initialize();
}


TSimpleImmersedValveProblem::TSimpleImmersedValveProblem(
        double lTimeEndValue,
        double lTimeStepValue,
        double lRe,
        double lStiffness
        ): TImmersedBoundaryProblem(lTimeEndValue, lTimeStepValue, lRe, lStiffness) {}


void TSimpleImmersedValveProblem::GetForce(
        TRnRelease3dSpace &force,
        T3dNormalGrid &gridC,
        T3dNormalGrid &grid,
        TRnRelease3dSpace &velocity,
        const int axis,
        const int timeStepNumber
        )
{
    TImmersedBoundaryProblem::GetForce(force, gridC, grid, velocity, axis, timeStepNumber);
}

double TSimpleImmersedValveProblem::GetDistribution(int i, int j, int k, double hx, double hy, double hz)
{
	return TImmersedBoundaryProblem::GetDistribution(i,j,k,hx,hy,hz);
}


TMaskGenerator* TCanalUnifyProblemWithImmersedSimpleValve::CreateStatement()
{
    //Preparing (Begin)
    double iPart;
     int endIndex;

    iPart = 0;
    modf(fLength/fStep,&iPart);
    endIndex = (int)iPart;
    int nx = endIndex + 1;

    iPart = 0;
    modf(fRadius*2/fStep,&iPart);
    endIndex = (int)iPart;
    int ny = endIndex + 1;
    int nz = endIndex + 1;
    //Preparing (End)

    //Creating (Begin)
    return TSimpleWithPressureMaskGenerator::CreateInstance(
                                                            fLength,fRadius*2,fRadius*2,
                                                            nx,ny,nz,
                                                            fPressureInput,fPressureOutput
                                                           );
    //Creating (End)
}


TCylinderBoundaryProblem::~TCylinderBoundaryProblem()
{
    //delete fBoundary;
}

void TCylinderBoundaryProblem::Initialize()
{

    TVPProblem::Initialize();
    int count = 14400;
    fBoundary = (new TCylinderBoundary())
        ->withStiffness(fStiffness)
        ->withNodesCount(count)
        ->withRadius(0.11)
        ->withCenter(0.18, 0.25)
        ->Initialize();
}


TCylinderBoundaryProblem::TCylinderBoundaryProblem(
        double lTimeEndValue,
        double lTimeStepValue,
        double lRe,
        double lStiffness
        ): TImmersedBoundaryProblem(lTimeEndValue, lTimeStepValue, lRe, lStiffness) {}

TCylinderBoundaryProblem::TCylinderBoundaryProblem(
        double lTimeEndValue,
        double lTimeStepValue,
        double lmu1,
        double lmu2,
        double lrho1,
        double lrho2,
        double lStiffness
        ): TImmersedBoundaryProblem(lTimeEndValue, lTimeStepValue, lmu1, lmu2, lrho1, lrho2, lStiffness) {}
		

void TCylinderBoundaryProblem::GetForce(
        TRnRelease3dSpace &force,
        T3dNormalGrid &gridC,
        T3dNormalGrid &grid,
        TRnRelease3dSpace &velocity,
        const int axis,
        const int timeStepNumber
        )
{
    TImmersedBoundaryProblem::GetForce(force, gridC, grid, velocity, axis, timeStepNumber);
}

double TCylinderBoundaryProblem::GetDistribution(int i, int j, int k, double hx, double hy, double hz)
{
	return TImmersedBoundaryProblem::GetDistribution(i,j,k,hx,hy,hz);
}

void TCylinderBoundaryProblem::ChangeStifness(
	const int timeStepNumber
	)
{
	TImmersedBoundaryProblem::ChangeStifness(timeStepNumber);
}


TMaskGenerator* TCanalUnifyProblemWithImmersedCylinder::CreateStatement()
{
    //Preparing (Begin)
    double iPart;
     int endIndex;

    iPart = 0;
    modf(fLength/fStep,&iPart);
    endIndex = (int)iPart;
    int nx = endIndex + 1;

    iPart = 0;
    modf(fRadius*2/fStep,&iPart);
    endIndex = (int)iPart;
    int ny = endIndex + 1;
    int nz = endIndex + 1;
    //Preparing (End)

    std::function<bool(int, int, int, int, int, int)> condition = [](int i, int j, int k, int N, int L, int M) {
        return (j - L/4)*(j - L/4) + (k - M/2)*(k - M/2) > (M/6)*(M/6);
    };


    //Creating (Begin)
    return TSimpleWithCylinderPressureMaskGenerator::CreateInstance(
                                                            fLength,fRadius*2,fRadius*2,
                                                            nx,ny,nz,
                                                            fPressureInput,fPressureOutput,
                                                            condition
                                                           );
    //Creating (End)
}

TMaskGenerator* TCanalIBMWithSourceProblem::CreateStatement()
{
	
    //Preparing (Begin)
    double iPart;
     int endIndex;

    iPart = 0;
    modf(fLength/fStep,&iPart);
    endIndex = (int)iPart;
    int nx = endIndex + 1;

    iPart = 0;
    modf(fRadius*2/fStep,&iPart);
    endIndex = (int)iPart;
    int ny = endIndex + 1;
    int nz = endIndex + 1;
    //Preparing (End)

    std::function<bool(int, int, int, int, int, int)> condition = [](int i, int j, int k, int N, int L, int M) {
        return (j - L/4 - 3)*(j - L/4 - 3) + (k - M/2)*(k - M/2) > (M/6)*(M/6);
    };

    //Creating (Begin)
    return TSimpleWithCylinderPressureMaskGenerator::CreateInstance(
                                                            fLength,fRadius*2,fRadius*2,
                                                            nx,ny,nz,
                                                            fPressureInput,fPressureOutput,
                                                            condition
                                                           );
    //Creating (End)
}

bool TCanalIBMWithSourceProblem::CylinderMask(int i, int j, int k, int L, int M, int K)
{
      return i == 0 && (j - L/4 + 3)*(j - L/4 + 3) + (k - M/2)*(k - M/2) < (L/6)*(M/6);
}

bool TCanalIBMWithSourceProblem::InitialDensityDistribution(int i, int j, int k, int L, int M, int K)
{
    return this->CylinderMask(i, j, k, L, M, K);
}

bool TCanalIBMWithSourceProblem::InitialConcentrationDistribution(int i, int j, int k, int L, int M, int K)
{
    return this->CylinderMask(i, j, k, L, M, K);
}

bool TCanalIBMWithSourceProblem::InitialViscosityDistribution(int i, int j, int k, int L, int M, int K)
{
    return this->CylinderMask(i, j, k, L, M, K);
}

bool TCanalIBMWithSourceProblem::ConcentrationInletMask(int i, int j, int k, int L, int M, int K)
{
    return this->CylinderMask(i, j, k, L, M, K);
}

const double TCanalIBMWithSourceProblem::ConcentrationInletCondition(int i, int j, int k, int L, int M, int K)
{
    return this->CylinderMask(i, j, k, L, M, K) ? 0.5 : 0.0;
}

TMaskGenerator* TCanalIBMWithThrombus::CreateStatement()
{
    //Preparing (Begin)
    double iPart;
     int endIndex;

    iPart = 0;
    modf(fLength/fStep,&iPart);
    endIndex = (int)iPart;
    int nx = endIndex + 1;

    iPart = 0;
    modf(fRadius*2/fStep,&iPart);
    endIndex = (int)iPart;
    int ny = endIndex + 1;
    int nz = endIndex + 1;
    //Preparing (End)

    std::function<bool(int, int, int, int, int, int)> condition = [](int i, int j, int k, int N, int L, int M) {
        return (j - L/4)*(j - L/4) + (k - M/2)*(k - M/2) > (M/6)*(M/6);
    };

    //Creating (Begin)
    return TSimpleWithCylinderPressureMaskGenerator::CreateInstance(
                                                            fLength,fRadius*2,fRadius*2,
                                                            nx,ny,nz,
                                                            fPressureInput,fPressureOutput,
                                                            condition
                                                           );
    //Creating (End)
}

bool TCanalIBMWithThrombus::CylinderMask(int i, int j, int k, int L, int M, int K)
{
      return (i-L/2)*(i-L/2) + (j-M/2)*(j-M/2) + (k-K/2)*(k-K/2) < 20;
}

bool TCanalIBMWithThrombus::InitialDensityDistribution(int i, int j, int k, int L, int M, int K)
{
    return this->CylinderMask(i, j, k, L, M, K);
}

bool TCanalIBMWithThrombus::InitialConcentrationDistribution(int i, int j, int k, int L, int M, int K)
{
    return this->CylinderMask(i, j, k, L, M, K);
}

bool TCanalIBMWithThrombus::InitialViscosityDistribution(int i, int j, int k, int L, int M, int K)
{
    return this->CylinderMask(i, j, k, L, M, K);
}

bool TCanalIBMWithThrombus::ConcentrationInletMask(int i, int j, int k, int L, int M, int K)
{
    return false;
}

const double TCanalIBMWithThrombus::ConcentrationInletCondition(int i, int j, int k, int L, int M, int K)
{
    return 0.0;
}

void TCanalIBMWithThrombus::Initialize()
{
    TVPProblem::Initialize();
    fBoundary = (TDeformedCylinderBoundary *)(new TDeformedCylinderBoundary())
        ->withStiffness(fStiffness)
        ->withNodesCount(14400)
        ->withRadius(0.11)
        ->withCenter(0.25, 0.25)
        ->withDeformation(0.02)
        ->Initialize();
}


TMaskGenerator* TCanalIBMWithElasticBoundary::CreateStatement()
{
    //Preparing (Begin)
    double iPart;
    int endIndex;

    iPart = 0;
    modf(fLengthX/fStep,&iPart);
    endIndex = (int)iPart;
    int nx = endIndex + 1;

    iPart = 0;
	modf(fLengthY/ fStep, &iPart);
    endIndex = (int)iPart;
    int ny = endIndex + 1;

	iPart = 0;
	modf(fLengthZ/ fStep, &iPart);
	endIndex = (int)iPart;
    int nz = endIndex + 1;
	
	
    //Preparing (End)
	
    std::function<bool(int, int, int, int, int, int)> condition = [](int i, int j, int k, int L, int M, int K) {
		
		return	(double(j) - double(M) / 2.)*(double(j) - double(M) / 2.) + (double(k) - double(K) / 2.)*(double(k) - double(K) / 2.) >
			(double(M) / 10.1)*(double(K) / 10.1);

    };

    //Creating (Begin)
    return TSimpleWithCylinderPressureMaskGenerator::CreateInstance(
		                                                    fLengthX, fLengthY, fLengthZ,
                                                            nx,ny,nz,
                                                            fPressureInput,fPressureOutput,
                                                            condition
                                                           );
    //Creating (End)
}

bool TCanalIBMWithElasticBoundary::CylinderMask(int i, int j, int k, int L, int M, int K)  //dira
{

	return ( (i == 0) && ((double(j) - double(M) / 2. - 1)*(double(j) - double(M) / 2. -1) + (double(k) - double(K) / 2. -1)*(double(k) -1 - double(K) / 2.)
	< (double(M) / 10.1)*(double(K) / 10.1)) );

}

bool TCanalIBMWithElasticBoundary::InitialDensityDistribution(int i, int j, int k, int L, int M, int K)
{
   // i = 0;
    return this->CylinderMask(i, j, k, L, M, K);
}

bool TCanalIBMWithElasticBoundary::InitialConcentrationDistribution(int i, int j, int k, int L, int M, int K)
{
   // i = 0;
    return this->CylinderMask(i, j, k, L, M, K);
}

bool TCanalIBMWithElasticBoundary::InitialViscosityDistribution(int i, int j, int k, int L, int M, int K)
{
   // i = 0;
    return this->CylinderMask(i, j, k, L, M, K);
}

bool TCanalIBMWithElasticBoundary::ConcentrationInletMask(int i, int j, int k, int L, int M, int K)
{
	//i = 0;
    return this->CylinderMask(i, j, k, L, M, K);
}

const double TCanalIBMWithElasticBoundary::ConcentrationInletCondition(int i, int j, int k, int L, int M, int K)
{
    //return this->CylinderMask(i, j, k, L, M, K) ? 0.5 : 0.0;
	return this->CylinderMask(i, j, k, L, M, K) ? 0.45 : 0.0;
}

void TCanalIBMWithElasticBoundary::Initialize()
{
    TVPProblem::Initialize();
    //int count = 14400 + 3282;
   // int count = 14400;
    fBoundary = (new TCylinderElasticBoundary())
        ->withStiffness(fStiffness)
        //->withNodesCount(count)
        ->withRadius(fRadius)
		->withHeight(fLengthX )//dlina pogr gran
		->withStep(fStep)
        ->withCenter(fLengthY/2, fLengthZ/2)
        ->Initialize();
}

void TImmersedBoundaryProblem::CountBoundaryP(TRnRelease3dSpace &pressure)
{
	int grid_marg = 2;
	int grid_rad = 3;
	double h = (fStatement->Hx[1] + fStatement->Hy[1] + fStatement->Hz[1]) / 3.;
	double hsqr = h*h;
	double ub = 0, pb = 0, wi = 0, xb = 0, yb = 0, zb = 0, xp = 0, yp = 0, zp = 0, sum_u = 0, sum_w = 0, dist = 0;

	for (int il = 0; il < fBoundary->height_nodes; il++)
	{
		xb = fBoundary->nodes[il*fBoundary->radius_nodes]->x;
		// find layer near point xc
		int x_near = 0;
		for (int i = 0; i < fStatement->NxU; i++)
		if (fabs(fStatement->XU[i] - xb) <= h)
		{
			x_near = i;
			break;
		}
		//printf("\nil=%d x_near=%d", il, x_near);
		for (int ir = 0; ir < fBoundary->radius_nodes; ir++)
		{
			xb = fBoundary->nodes[il*fBoundary->radius_nodes + ir]->x;
			yb = fBoundary->nodes[il*fBoundary->radius_nodes + ir]->y;
			zb = fBoundary->nodes[il*fBoundary->radius_nodes + ir]->z;
			//printf("\nir=%d", ir);

			//-----------------------------------------------------
			// Counting p in imm. boundary node
			//-----------------------------------------------------
			pb = 0.;
			sum_w = 0.;
			sum_u = 0.;
			for (int ix = x_near - grid_marg>0 ? (x_near - grid_marg) : 0; ix < x_near + grid_marg + grid_rad && ix<fStatement->NxP; ix++)
			{
				//printf("\nix >> %d",ix);
				bool no_near_point = true;
				bool first_near_point = false;
				for (int iy = 1; iy < fStatement->NyU - 1; iy++)
					for (int iz = 1; iz < fStatement->NzU - 1; iz++)
					{
						xp = fStatement->XP[ix];
						yp = fStatement->YP[iy];
						zp = fStatement->ZP[iz];
						dist = pow(xb - xp, 2) + pow(yb - yp, 2) + pow(zb - zp, 2);
						if (dist < hsqr)
						{
							if (!first_near_point) first_near_point = true;
							no_near_point = false;
							wi = 1. / dist;
							//printf("\npoint %d %d %d p=%lf d=%lf wi=%lf", ix, iy, iz, (double) pressure[ix][iy][iz],dist,wi);

							sum_w += wi*10.e-2;
							sum_u += pressure[ix][iy][iz] * wi*10.e-2;
							//printf("\nres sum_u=%lf  sum_w=%lf", sum_u, sum_w);
						}

				}
				if (first_near_point && no_near_point) break;
			}
			pb += sum_u / sum_w;
			fBoundary->nodes[il*fBoundary->radius_nodes + ir]->p = pb;
		//	printf("\nOK:res sum_u=%lf  sum_w=%lf p=%lf\n",sum_u,sum_w,pb);


		}
	}

	//approxP(pressure);

	//approxUP(TRnRelease3dSpace&, TRnRelease3dSpace&);

}

void TImmersedBoundaryProblem::CountBoundaryConcentration(TRnRelease3dSpace &c)
{
	int grid_marg = 2;
	int grid_rad = 3;
	double h = (fStatement->Hx[1] + fStatement->Hy[1] + fStatement->Hz[1]) / 3.;
	double hsqr = h*h;
	double ub = 0, pb = 0, wi = 0, xb = 0, yb = 0, zb = 0, xp = 0, yp = 0, zp = 0, sum_u = 0, sum_w = 0, dist = 0;

	for (int il = 0; il < fBoundary->height_nodes; il++)
	{
		xb = fBoundary->nodes[il*fBoundary->radius_nodes]->x;
		// find layer near point xc
		int x_near = 0;
		for (int i = 0; i < fStatement->NxU; i++)
		if (fabs(fStatement->XU[i] - xb) <= h)
		{
			x_near = i;
			break;
		}
		//printf("\nil=%d x_near=%d", il, x_near);
		for (int ir = 0; ir < fBoundary->radius_nodes; ir++)
		{
			xb = fBoundary->nodes[il*fBoundary->radius_nodes + ir]->x;
			yb = fBoundary->nodes[il*fBoundary->radius_nodes + ir]->y;
			zb = fBoundary->nodes[il*fBoundary->radius_nodes + ir]->z;
			//printf("\nir=%d", ir);

			//-----------------------------------------------------
			// Counting p in imm. boundary node
			//-----------------------------------------------------
			pb = 0.;
			sum_w = 0.;
			sum_u = 0.;
			for (int ix = x_near - grid_marg>0 ? (x_near - grid_marg) : 0; ix < x_near + grid_marg + grid_rad && ix<fStatement->NxP; ix++)
			{
				//printf("\nix >> %d",ix);
				bool no_near_point = true;
				bool first_near_point = false;
				for (int iy = 1; iy < fStatement->NyU - 1; iy++)
					for (int iz = 1; iz < fStatement->NzU - 1; iz++)
					{
						xp = fStatement->XP[ix];
						yp = fStatement->YP[iy];
						zp = fStatement->ZP[iz];
						dist = pow(xb - xp, 2) + pow(yb - yp, 2) + pow(zb - zp, 2);
						if (dist < hsqr)
						{
							if (!first_near_point) first_near_point = true;
							no_near_point = false;
							wi = 1. / dist;
							//printf("\npoint %d %d %d p=%lf d=%lf wi=%lf", ix, iy, iz, (double) c[ix][iy][iz],dist,wi);

							sum_w += wi*10.e-2;
							sum_u += c[ix][iy][iz] * wi*10.e-2;
							//printf("\nres sum_u=%lf  sum_w=%lf", sum_u, sum_w);
						}

				}
				if (first_near_point && no_near_point) break;
			}
			pb += sum_u / sum_w;
			fBoundary->nodes[il*fBoundary->radius_nodes + ir]->concentration = pb;
		//	printf("\nOK:res sum_u=%lf  sum_w=%lf p=%lf\n",sum_u,sum_w,pb);


		}
	}

	//approxP(c);

	//approxUP(TRnRelease3dSpace&, TRnRelease3dSpace&);

}

void TImmersedBoundaryProblem::CountBoundaryViscosity(TRnRelease3dSpace &vis)
{
	int grid_marg = 2;
	int grid_rad = 3;
	double h = (fStatement->Hx[1] + fStatement->Hy[1] + fStatement->Hz[1]) / 3.;
	double hsqr = h*h;
	double ub = 0, pb = 0, wi = 0, xb = 0, yb = 0, zb = 0, xp = 0, yp = 0, zp = 0, sum_u = 0, sum_w = 0, dist = 0;

	for (int il = 0; il < fBoundary->height_nodes; il++)
	{
		xb = fBoundary->nodes[il*fBoundary->radius_nodes]->x;
		// find layer near point xc
		int x_near = 0;
		for (int i = 0; i < fStatement->NxU; i++)
		if (fabs(fStatement->XU[i] - xb) <= h)
		{
			x_near = i;
			break;
		}
		//printf("\nil=%d x_near=%d", il, x_near);
		for (int ir = 0; ir < fBoundary->radius_nodes; ir++)
		{
			xb = fBoundary->nodes[il*fBoundary->radius_nodes + ir]->x;
			yb = fBoundary->nodes[il*fBoundary->radius_nodes + ir]->y;
			zb = fBoundary->nodes[il*fBoundary->radius_nodes + ir]->z;
			//printf("\nir=%d", ir);

			//-----------------------------------------------------
			// Counting p in imm. boundary node
			//-----------------------------------------------------
			pb = 0.;
			sum_w = 0.;
			sum_u = 0.;
			for (int ix = x_near - grid_marg>0 ? (x_near - grid_marg) : 0; ix < x_near + grid_marg + grid_rad && ix<fStatement->NxP; ix++)
			{
				//printf("\nix >> %d",ix);
				bool no_near_point = true;
				bool first_near_point = false;
				for (int iy = 1; iy < fStatement->NyU - 1; iy++)
					for (int iz = 1; iz < fStatement->NzU - 1; iz++)
					{
						xp = fStatement->XP[ix];
						yp = fStatement->YP[iy];
						zp = fStatement->ZP[iz];
						dist = pow(xb - xp, 2) + pow(yb - yp, 2) + pow(zb - zp, 2);
						if (dist < hsqr)
						{
							if (!first_near_point) first_near_point = true;
							no_near_point = false;
							wi = 1. / dist;
							//printf("\npoint %d %d %d p=%lf d=%lf wi=%lf", ix, iy, iz, (double) vis[ix][iy][iz],dist,wi);

							sum_w += wi*10.e-2;
							sum_u += vis[ix][iy][iz] * wi*10.e-2;
							//printf("\nres sum_u=%lf  sum_w=%lf", sum_u, sum_w);
						}

				}
				if (first_near_point && no_near_point) break;
			}
			pb += sum_u / sum_w;
			fBoundary->nodes[il*fBoundary->radius_nodes + ir]->vis = pb;
		//	printf("\nOK:res sum_u=%lf  sum_w=%lf p=%lf\n",sum_u,sum_w,pb);


		}
	}

	//approxP(pressure);

	//approxUP(TRnRelease3dSpace&, TRnRelease3dSpace&);

}