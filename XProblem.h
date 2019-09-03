#ifndef ProblemsH
#define ProblemsH

#include "XGrid.h"
#include "XSpace.h"
#include "XStatement.h"
#include "XImmersedBoundary.h"
#include <csignal>

#include <math.h>


//#################################################################

class TProblem
{
    public:
	 TProblem() {};
	 virtual ~TProblem() {};
     virtual void Solve()=0;
     
	 //Temporary step (Begin)
	 static bool GlobalTurbulenceMode;
	 static double GlobalTurbulenceParameter;
	 static bool GlobalUsePackOperator;
	 static char* GlobalCatalog;
	 static double GlobalLength;
	 static double GlobalVelocity;
	 static double GlobalPressure;
	 static int GlobalSaveStep;
	 
	 //Temporary step (End)
    protected:      
    private:
      TProblem(const TProblem &);
      TProblem &operator=(const TProblem &);
};

//#################################################################

class T1dPoissonProblem: public TProblem
{
    public:
	 T1dPoissonProblem(double lH);
	 virtual ~T1dPoissonProblem();
     virtual void Solve();
    protected:      
    private:
      double fH;
      T1dNormalGrid* fGrid;	  

	  double SolutionFunction(double x);
	  double SolutionFunctionDX(double x);
	  double RightPartFunction(double x);

      T1dPoissonProblem(const T1dPoissonProblem &);
      T1dPoissonProblem &operator=(const T1dPoissonProblem &);
};

//#################################################################

class T3dPoissonProblem: public TProblem
{
    public:
	 T3dPoissonProblem(double lH);
	 virtual ~T3dPoissonProblem();
     virtual void Solve();
    protected:      
    private:
      double fH;
      T3dNormalGrid* fGrid;	  
	  T3dNormalGrid* fGridGeometry;	  

	  double SolutionFunction(double x, double y, double z);
	  double SolutionFunctionDX(double x, double y, double z);
	  double SolutionFunctionDY(double x, double y, double z);
	  double SolutionFunctionDZ(double x, double y, double z);
	  double RightPartFunction(double x, double y, double z);

      T3dPoissonProblem(const T3dPoissonProblem &);
      T3dPoissonProblem &operator=(const T3dPoissonProblem &);
};

//#################################################################

class THydrodynamicsProblem: public TProblem
{
    public:	
     THydrodynamicsProblem(double lTimeEndValue, double lTimeStepValue, double lRe):TProblem(),
		                                                                                      fTimeEndValue(lTimeEndValue),fTimeStepValue(lTimeStepValue),
																							  fRe(lRe),
																							  fGrid(NULL)                                                                                         
	 {
		 if (fTimeEndValue<=0) {throw "Error in THydrodynamicsProblem::Constructor: Time end value is incorrect.";}
		 if  (fTimeStepValue<=0) {throw "Error in THydrodynamicsProblem::Constructor: Time step value is incorrect.";}
		 if (fRe<=0) {throw "Error in THydrodynamicsProblem::Constructor: Reinolds value is incorrect.";}
	 };

	THydrodynamicsProblem(double lTimeEndValue, double lTimeStepValue, double lmu1, double lmu2, double lrho1, double lrho2):TProblem(),
		                                                                                      fTimeEndValue(lTimeEndValue),fTimeStepValue(lTimeStepValue),
																							  fMu1(lmu1),
																							  fMu2(lmu2),
																							  fRho1(lrho1),
																							  fRho2(lrho2),
																							  fGrid(NULL)                                                                                         
	{
		if (fTimeEndValue<=0) {throw "Error in THydrodynamicsProblem::Constructor: Time end value is incorrect.";}
		if  (fTimeStepValue<=0) {throw "Error in THydrodynamicsProblem::Constructor: Time step value is incorrect.";}
		if (fMu1<=0) {throw "Error in THydrodynamicsProblem::Constructor: First mu value is incorrect.";}
		if (fMu2<=0) {throw "Error in THydrodynamicsProblem::Constructor: First mu value is incorrect.";}
		if (fRho1<=0) {throw "Error in THydrodynamicsProblem::Constructor: First mu value is incorrect.";}
		if (fRho2<=0) {throw "Error in THydrodynamicsProblem::Constructor: First mu value is incorrect.";}
	};
	 
	virtual ~THydrodynamicsProblem() {delete fGrid;}     
    protected:      
     
	virtual TNodesGrid* BuildGrid() = 0;

	 void ReBuildGrid()
	 {
		 delete fGrid;
		 fGrid = BuildGrid();
	 }

	 void Initialize() 
	 {
		 fTimeSeparator = new TUniformSeparator(0,fTimeEndValue,fTimeStepValue);
		 fGrid = BuildGrid();
	 }
	 	 
	 TUniformSeparator* fTimeSeparator;
     double fRe;
     double fMu1;
     double fMu2;
     double fRho1;
     double fRho2;
	 TNodesGrid* fGrid; 

	 //Note: may be any dimension, not only 3d - to be continue...
	 void SolveTransferEquation(TRnRelease3dSpace& lVecCur, TRnRelease3dSpace& lVecPrev, 
		                         TRnRelease3dSpace&  lU, TRnRelease3dSpace& lV, TRnRelease3dSpace& lW, TRnRelease3dSpace& lRP, 
								 const TNodesGrid& lGrid);

	 void SolveTransferEquation(TRnRelease3dSpace& lVecCur, TRnRelease3dSpace& lVecPrev, 
		                         TRnRelease3dSpace&  lU, TRnRelease3dSpace& lV, TRnRelease3dSpace& lW, TRnRelease3dSpace& lRP, 
								 const TNodesGrid& lGrid, double lDiffusionCoefficient);

	 void SolveTransferEquation_(TRnRelease3dSpace& lVecCur, TRnRelease3dSpace& lVecPrev, 
		                         TRnRelease3dSpace&  lU, TRnRelease3dSpace& lV, TRnRelease3dSpace& lW, TRnRelease3dSpace& lRP, 
								 const TNodesGrid& lGrid, double lDiffusionCoefficient);

	 void SolveTransferEquation__(TRnRelease3dSpace& lVecCur, TRnRelease3dSpace& lVecPrev, 
		                         TRnRelease3dSpace&  lU, TRnRelease3dSpace& lV, TRnRelease3dSpace& lW, TRnRelease3dSpace& lRP, 
								 const TNodesGrid& lGrid, TRnRelease3dSpace& lDiffusionCoefficient, int lSchemeType);

    private:            
      double fTimeEndValue; 
	  double fTimeStepValue;

      THydrodynamicsProblem(const THydrodynamicsProblem &);
      THydrodynamicsProblem &operator=(const THydrodynamicsProblem &);
};

//#################################################################

const double PI = 3.14159265358979323846;
const double alfa = 1;
class T3dTranferProblem: public THydrodynamicsProblem
{
    public:
     static T3dTranferProblem* CreateInstance(double lTimeEndValue, double lTimeStepValue, double lRe, double lSpaceStepValue) 
	 {
       T3dTranferProblem* tp = new T3dTranferProblem(lTimeEndValue,lTimeStepValue,lRe,lSpaceStepValue);
	   tp->Initialize();
	   return tp;
	 }
     
	 virtual ~T3dTranferProblem(); 

	 virtual void Solve();
    protected:   
      T3dTranferProblem(double lTimeEndValue, double lTimeStepValue, double lRe, double lSpaceStepValue):
			            THydrodynamicsProblem(lTimeEndValue,lTimeStepValue,lRe), 
		                fStep(lSpaceStepValue) {};

	  //From THydrodynamicsProblem (Begin)
	  virtual TNodesGrid* BuildGrid();
	  //From THydrodynamicsProblem (End)
      
      virtual TNodesGrid* BuildGridSolution();

	  virtual TNodesGrid* BuildNormalGrid(TSeparator& lSep1, TSeparator& lSep2, TSeparator& lSep3, const TMask& lMask);

      void PrepareVariables();

	  void Initialize() 
	  {
		 THydrodynamicsProblem::Initialize();

		 fGridSolution = BuildGridSolution();		 
		 PrepareVariables();
	  } 

      TNodesGrid* fGridSolution; 
	  
	  TRnRelease3dSpace* fSolution;	  
    private:                        
      double fStep;	 

	  static double U(double x, double y, double z, double t) {return 1;}
	  static double V(double x, double y, double z, double t) {return 1;}
	  static double W(double x, double y, double z, double t) {return 1;}

	  static double SolutionFunction(double x, double y, double z, double t) {return sin(t)*sin(alfa*x)*sin(alfa*y)*sin(alfa*z);}

	  static double DT(double x, double y, double z, double t) {return cos(t)*sin(alfa*x)*sin(alfa*y)*sin(alfa*z);}

	  static double DX(double x, double y, double z, double t) {return sin(t)*alfa*cos(alfa*x)*sin(alfa*y)*sin(alfa*z);}
	  static double DY(double x, double y, double z, double t) {return sin(t)*sin(alfa*x)*alfa*cos(alfa*y)*sin(alfa*z);}
	  static double DZ(double x, double y, double z, double t) {return sin(t)*sin(alfa*x)*sin(alfa*y)*alfa*cos(alfa*z);}

	  static double D2X(double x, double y, double z, double t) {return sin(t)*alfa*alfa*(-sin(alfa*x))*sin(alfa*y)*sin(alfa*z);}
	  static double D2Y(double x, double y, double z, double t) {return sin(t)*sin(alfa*x)*alfa*alfa*(-sin(alfa*y))*sin(alfa*z);}
	  static double D2Z(double x, double y, double z, double t) {return sin(t)*sin(alfa*x)*sin(alfa*y)*alfa*alfa*(-sin(alfa*z));}
	  
	  double RightPartFunction(double x, double y, double z, double t) 
	  {
         double Re = fRe; 
		 
		 double u = U(x,y,z,t);
		 double v = V(x,y,z,t);
		 double w = W(x,y,z,t);

		 double dt = DT(x,y,z,t);
		 
		 double dx = DX(x,y,z,t);
		 double dy = DY(x,y,z,t);
		 double dz = DZ(x,y,z,t);

		 double d2x = D2X(x,y,z,t);
		 double d2y = D2Y(x,y,z,t);
		 double d2z = D2Z(x,y,z,t);


	     return  dt + u*dx + v*dy + w*dz - (1/Re)*(d2x + d2y + d2z);
	  }

      T3dTranferProblem(const T3dTranferProblem &);
      T3dTranferProblem &operator=(const T3dTranferProblem &);
};


//#################################################################

class TVPProblem: public THydrodynamicsProblem
{
    public:	
	TVPProblem(double lTimeEndValue, double lTimeStepValue, double lRe):
		           THydrodynamicsProblem(lTimeEndValue,lTimeStepValue,lRe) 
				   {fStatement = NULL; fMethod = 1;}

	TVPProblem(double lTimeEndValue, double lTimeStepValue, double lmu1, double lmu2, double lrho1, double lrho2):
		           THydrodynamicsProblem(lTimeEndValue,lTimeStepValue,lmu1,lmu2,lrho1,lrho2) 
				   {fStatement = NULL; fMethod = 1;}
	
	virtual ~TVPProblem();

	 void Solve() 
	 {
		 if (fMethod==1) {SolveBySplitting();}
		 else {throw "Error in TVPProblem::Solve: Unknown method.";}
	 }
    protected:
       
	  virtual TMaskGenerator* CreateStatement() = 0;
	  virtual bool RestructStatement(TMaskGenerator* lStatement) {return false;}
	  virtual bool RestructStatement(TMaskGenerator* lStatement, TRnRelease3dSpace* lU, TRnRelease3dSpace* lV, TRnRelease3dSpace* lW, TRnRelease3dSpace* lP) {return false;}

	  virtual TNodesGrid* BuildLocalGrid(
		                                 double *lX, double *lY, double *lZ,		                                 
										 int lNx, int lNy, int lNz,
										 glTMaskWithNormals ***lMask
										);
		

	  //From THydrodynamicsProblem (Begin)
		virtual TNodesGrid* BuildGrid();
	  //From THydrodynamicsProblem (End)

	  virtual TNodesGrid* BuildGridU();
	  virtual double GetBorderConditionU(int lI, int lJ, int lK, int lN);
	  virtual TNodesGrid* BuildGridV();
	  virtual double GetBorderConditionV(int lI, int lJ, int lK, int lN);
	  virtual TNodesGrid* BuildGridW();
	  virtual double GetBorderConditionW(int lI, int lJ, int lK, int lN);
	  virtual TNodesGrid* BuildGridP();
	  virtual double GetBorderConditionP(int lI, int lJ, int lK, int lN);

	  void LocalPrepareVariables();

	  virtual void PrepareVariables();


	  void Initialize() 
	  {		 
		 fStatement = CreateStatement();
		  
		 THydrodynamicsProblem::Initialize();

		 fGridU = BuildGridU();
		 fGridV = BuildGridV();
		 fGridW = BuildGridW();
		 fGridP = BuildGridP();

		 PrepareVariables();
	  }
      

	  void ReBuildGrids()
	  {
		 THydrodynamicsProblem::ReBuildGrid();
		 delete fGridU;
		 delete fGridV;
		 delete fGridW;
		 delete fGridP;
		 fGridU = BuildGridU();
		 fGridV = BuildGridV();
		 fGridW = BuildGridW();
		 fGridP = BuildGridP();
	  }





	  TMaskGenerator* fStatement;

	  
	  int fSpaceModelNameU;
	  TRnRelease3dSpace* fU;
	  TNodesGrid* fGridU;

	  int fSpaceModelNameV;
	  TRnRelease3dSpace* fV;
	  TNodesGrid* fGridV;

	  int fSpaceModelNameW;
	  TRnRelease3dSpace* fW;
	  TNodesGrid* fGridW; 

	  int fSpaceModelNameP;
	  TRnRelease3dSpace* fP;
	  TNodesGrid* fGridP; 

	  int fMethod;
	  
	  void SolveBySplitting();
	  virtual void GetForce(
                TRnRelease3dSpace &force,
                T3dNormalGrid &gridC,
                T3dNormalGrid &grid,
                TRnRelease3dSpace &velocity,
                const int axis,
                const int timeStepNumber
                )
      {
          force |= 0;
      }

	  virtual double GetDistribution(int, int, int,double, double, double)
	  { 
		  return 0.0; 
	  }

	  virtual void ChangeStifness(const int timeStepNumber)
	  {

	  }

	 
	  virtual void InterpolateUpdatePathes(TRnRelease3dSpace &velocity, T3dNormalGrid *grid, const int axis, const int timeStepNumber){}
	  virtual void OutputBoundary(int timeStepNumber) {}
	  virtual void CountBoundaryP(TRnRelease3dSpace &pressure){}
	  virtual void CountBoundaryConcentration(TRnRelease3dSpace &c){}
	  virtual void CountBoundaryViscosity(TRnRelease3dSpace &vis){}
	  virtual bool InitialDensityDistribution(int i, int j, int k, int N, int L, int M) {return false;};
	  virtual bool InitialConcentrationDistribution(int i, int j, int k, int N, int L, int M) {return false;};
	  virtual bool InitialViscosityDistribution(int i, int j, int k, int N, int L, int M) {return false;};
	  virtual bool ConcentrationInletMask(int i, int j, int k, int N, int L, int K) {return false;};
	  virtual const double ConcentrationInletCondition(int i, int j, int k, int N, int L, int M) {return false;};

	  //============
	  virtual void ComputeChanging(int timeStepNumber) {}
	  virtual void ComputeStress() {}
	  virtual void ComputeStress(TRnRelease3dSpace &U, TRnRelease3dSpace &V, TRnRelease3dSpace &W, int timestepNumber) {}
	  //============

    private:            
      TVPProblem(const TVPProblem &);
      TVPProblem &operator=(const TVPProblem &);
};


//#################################################################


class TScour
{
    public:
	  virtual char* GetName() = 0;
};

class TMorphologicScour: public TScour
{
   public:
	    virtual char* GetName() {return NULL;};

        double fDiameter;
        double fDynamicFriction;
        double fStaticFriction;
        double fShieldsNumber;
        double fDensity;
        double fPorosity;
        double fSlip;
};

class TScourProblem: public TVPProblem
{
    public:	
	  TScourProblem(double lTimeEndValue, double lTimeStepValue, double lRe, double lFluidDensity, const TMorphologicScour& lScour):
		           TVPProblem(lTimeEndValue,lTimeStepValue,lRe),
				   fFluidDensity(lFluidDensity),
				   fDiameter(lScour.fDiameter),fDynamicFriction(lScour.fDynamicFriction),fStaticFriction(lScour.fStaticFriction),fShieldsNumber(lScour.fShieldsNumber),
				   fPorosity(lScour.fPorosity),fSlip(lScour.fSlip)
				   {}
		

      void PrintHeight(int t);

	  virtual ~TScourProblem() {}


    protected:
      
	  //From TVPProblem (Begin)
      virtual void PrepareVariables()
	  {
         TVPProblem::PrepareVariables();
		 
		 Nx = fStatement->NxP;
		 Ny = fStatement->NyP;
		 Nz = fStatement->NzP;
		 Hx = fStatement->HxP;
		 Hy = fStatement->HyP;
		 Hz = fStatement->HzP;
		 Cx = fStatement->XP;
		 Cy = fStatement->YP;
		 Cz = fStatement->ZP;

		 glAllocPlane(z_indx, Nx, Ny);
		 glAllocPlane(h, Nx, Ny);

		 glAllocMatrix(Gp, Nx, Ny, Nz);

		 for(int i=0; i<Nx; ++i)
			for(int j=0; j<Ny; ++j)
				for(int k=0; k<Nz; ++k)
					Gp[i][j][k] = fStatement->MaskP[i][j][k].mask;

		 tau = *fTimeSeparator->Separation;

         
		 PrepareScourFields();		 
	  }
	        
	  virtual bool RestructStatement(TMaskGenerator* lStatement);
	  //From TVPProblem (Begin)

	  virtual void PrepareScourFields() = 0;
	  //Scour fields (Begin)
	  int
		Nx, Ny, Nz,
		**z_indx;

	  double
		**h;

	   glTVector
		**q_b;

	  int
	   **gorisont;

	  double
		**P_ef_print;
	  glTVector
		**U_b_print;
	  //Scour fields (End)
    private: 
      //Help fields (Begin)
	  double
	  tau;

	  int
		***Gp;

	  double
		*Hx, *Hy, *Hz,
		*Cx, *Cy, *Cz;

	  FILE *fsc;
      //Help fields (End)

	  //Foundation parameters (Begin)
      double fFluidDensity;	  
	  double fDiameter;
      double fDynamicFriction;
      double fStaticFriction;
      double fShieldsNumber;
      double fDensity;
      double fPorosity;
      double fSlip;
      //Foundation parameters (End)

	  bool SolveSys(double a, double b, double c, double d, double k, double &x1, double &x2, double &x3);

      TScourProblem(const TScourProblem &);
      TScourProblem &operator=(const TScourProblem &);
};


//#################################################################


class TBarrelProblem: public TVPProblem
{
    public:	
	 static TBarrelProblem* CreateInstance(double lTimeEndValue, double lTimeStepValue, double lRe, double lRadius, double lHeight, double lSpaceStepValue) 
	 {
       TBarrelProblem* bp = new TBarrelProblem(lTimeEndValue,lTimeStepValue,lRe,lRadius,lHeight,lSpaceStepValue);
	   bp->Initialize();
	   return bp;
	 }
	 virtual ~TBarrelProblem() {};	 
    protected:  		
		TBarrelProblem(double lTimeEndValue, double lTimeStepValue, double lRe, double lRadius, double lHeight, double lSpaceStepValue):
		               TVPProblem(lTimeEndValue,lTimeStepValue,lRe),fRadius(lRadius),fHeight(lHeight),fStep(lSpaceStepValue) {}
					
		//From THydrodynamicsProblem (Begin)
		virtual TNodesGrid* BuildGrid();
		//From THydrodynamicsProblem (End)

		//From TVPProblem (Begin)
		virtual TMaskGenerator* CreateStatement() {return NULL;} // own grids - no statement		

		virtual TNodesGrid* BuildGridU();
	    virtual double GetBorderConditionU(int lI, int lJ, int lK, int lN);
	    virtual TNodesGrid* BuildGridV();
	    virtual double GetBorderConditionV(int lI, int lJ, int lK, int lN);
	    virtual TNodesGrid* BuildGridW();
	    virtual double GetBorderConditionW(int lI, int lJ, int lK, int lN);
	    virtual TNodesGrid* BuildGridP();
	    virtual double GetBorderConditionP(int lI, int lJ, int lK, int lN);
		//From TVPProblem (End)

		virtual TNodesGrid* BuildNormalGrid(TSeparator& lSep1, TSeparator& lSep2, TSeparator& lSep3, const TMask& lMask);

    private:
      double fRadius;
	  double fHeight;
	  double fStep;
      
      TBarrelProblem(const TBarrelProblem &);
      TBarrelProblem &operator=(const TBarrelProblem &);
};


//#################################################################


class TCanalProblem: public TVPProblem
{
    public:	
	 static TCanalProblem* CreateInstance(double lTimeEndValue, double lTimeStepValue, double lRe, double lRadius, double lLength, double lSpaceStepValue) 
	 {
       TCanalProblem* cp = new TCanalProblem(lTimeEndValue,lTimeStepValue,lRe,lRadius,lLength,lSpaceStepValue);
	   cp->Initialize();
	   return cp;
	 }
	 virtual ~TCanalProblem() {};	 
    protected:  		
		TCanalProblem(double lTimeEndValue, double lTimeStepValue, double lRe, double lRadius, double lLength, double lSpaceStepValue):
		               TVPProblem(lTimeEndValue,lTimeStepValue,lRe),fRadius(lRadius),fLength(lLength),fStep(lSpaceStepValue) {}
					
		//From THydrodynamicsProblem (Begin)
		virtual TNodesGrid* BuildGrid();
		//From THydrodynamicsProblem (End)

		//From TVPProblem (Begin)
		virtual TMaskGenerator* CreateStatement() {return NULL;} // own grids - no statement		

		virtual TNodesGrid* BuildGridU();
	    virtual double GetBorderConditionU(int lI, int lJ, int lK, int lN);
	    virtual TNodesGrid* BuildGridV();
	    virtual double GetBorderConditionV(int lI, int lJ, int lK, int lN);
	    virtual TNodesGrid* BuildGridW();
	    virtual double GetBorderConditionW(int lI, int lJ, int lK, int lN);
	    virtual TNodesGrid* BuildGridP();
	    virtual double GetBorderConditionP(int lI, int lJ, int lK, int lN);
		//From TVPProblem (End)

		virtual TNodesGrid* BuildNormalGrid(TSeparator& lSep1, TSeparator& lSep2, TSeparator& lSep3, const TMask& lMask);

    private:
      double fRadius;
	  double fLength;
	  double fStep;
      
      TCanalProblem(const TCanalProblem &);
      TCanalProblem &operator=(const TCanalProblem &);
};


//#################################################################


class TCanalUnifyProblem: public TVPProblem
{
    public:	
	 static TCanalUnifyProblem* CreateInstance(double lTimeEndValue, double lTimeStepValue, double lRe, double lRadius, double lLength, double lPressureInput, double lPressureOutput, double lSpaceStepValue) 
	 {
       TCanalUnifyProblem* cup = new TCanalUnifyProblem(lTimeEndValue,lTimeStepValue,lRe,lRadius,lLength,lPressureInput,lPressureOutput,lSpaceStepValue);
	   cup->Initialize();
	   return cup;
	 }
	 virtual ~TCanalUnifyProblem() {};	 
    protected:  		
		TCanalUnifyProblem(double lTimeEndValue, double lTimeStepValue, double lRe, double lRadius, double lLength, double lPressureInput, double lPressureOutput, double lSpaceStepValue):
		               TVPProblem(lTimeEndValue,lTimeStepValue,lRe),
					   fRadius(lRadius),fLength(lLength),
					   fPressureInput(lPressureInput),fPressureOutput(lPressureOutput),
					   fStep(lSpaceStepValue) {}
							
		//From TVPProblem (Begin)
		virtual TMaskGenerator* CreateStatement();		
		//From TVPProblem (End)
		
    private:
      double fRadius;
	  double fLength;

	  double fPressureInput;
	  double fPressureOutput;
      
	  double fStep;
      
      TCanalUnifyProblem(const TCanalUnifyProblem &);
      TCanalUnifyProblem &operator=(const TCanalUnifyProblem &);
};


//#################################################################


class TFlowPastCubeProblem: public TVPProblem
{
    public:	
	 static TFlowPastCubeProblem* CreateInstance(double lTimeEndValue, double lTimeStepValue, double lRe, double lLength, double lWidth, double lHeight, double lCubeSize, double lPressureInput, double lPressureOutput, double lSpaceStepValue) 
	 {
       TFlowPastCubeProblem* fpcp = new TFlowPastCubeProblem(lTimeEndValue,lTimeStepValue,lRe,lLength,lWidth,lHeight,lCubeSize,lPressureInput,lPressureOutput,lSpaceStepValue);
	   fpcp->Initialize();
	   return fpcp;
	 }
	 virtual ~TFlowPastCubeProblem() {};	 
    protected:  		
		TFlowPastCubeProblem(double lTimeEndValue, double lTimeStepValue, double lRe, double lLength, double lWidth, double lHeight, double lCubeSize, double lPressureInput, double lPressureOutput, double lSpaceStepValue):
		               TVPProblem(lTimeEndValue,lTimeStepValue,lRe),fLength(lLength),fWidth(lWidth),fHeight(lHeight),
						                                            fCubeSize(lCubeSize), 
																	fPressureInput(lPressureInput),fPressureOutput(lPressureOutput),
																	fStep(lSpaceStepValue) {}					                                               
			

		//From TVPProblem (Begin)
		virtual TMaskGenerator* CreateStatement();		
		//From TVPProblem (End)
		
    private:
      double fLength;
	  double fWidth;
	  double fHeight;

	  double fCubeSize;
	  
	  double fStep;

	  double fPressureInput;
	  double fPressureOutput;
      
      TFlowPastCubeProblem(const TFlowPastCubeProblem &);
      TFlowPastCubeProblem &operator=(const TFlowPastCubeProblem &);
};


//#################################################################



class TObstacle
{
    public:
	  virtual char* GetName() = 0;
};

class TIniParameters
{
   public:
	    double fAreaLength;
	    double fAreaWidth;
	    double fAreaHeight;

	    double fAnchorX;
	    double fAnchorY;
	    double fAnchorZ;
	  
		TObstacle* fObstacle;
	    
	    double fTimeFinish;

	    double fFluidDensity; 
        double fFluidViscosity;

	    double fPressureDifference;
	            
        TScour* fScour;
 	 	  
	    double fGridSpace;
	    double fGridTime;
};


//#################################################################


class TBargeObstacle: public TObstacle
{
    public:
	   virtual char* GetName() {return NULL;}  

	   double fWidth;
	   double fLengthBottom;
	   double fLengthTop;
	   double fAngle;
};

class TBargeProblem: public TScourProblem
{
    public:	
	 static TBargeProblem* CreateInstance(const TIniParameters& fIniParameters) 
	 {       
       const TIniParameters& params = fIniParameters;

	   const TBargeObstacle* bargeP = dynamic_cast<const TBargeObstacle*>(params.fObstacle);
	   const TMorphologicScour* scourP = dynamic_cast<const TMorphologicScour*>(params.fScour);
	   const TBargeObstacle& barge = *bargeP;
	   const TMorphologicScour& scour = *scourP;


	   //UnSizing (Begin)
	   double L_C = params.fAreaWidth;

	   double areaLength = params.fAreaLength/L_C;
	   double areaWidth = params.fAreaWidth/L_C;
	   double areaHeight = params.fAreaHeight/L_C;

	   double anchorX = params.fAnchorX/L_C;
	   double anchorY = params.fAnchorY/L_C;
	   double anchorZ = params.fAnchorZ/L_C;
	  
	   double bargeWidth = barge.fWidth/L_C;	   
	   double bargeLengthBottom = barge.fLengthBottom/L_C;
	   double bargeLengthTop = barge.fLengthTop/L_C;
	   double bargeAngle = barge.fAngle; 

	   double gridSpace = params.fGridSpace/L_C;


	   double U_C = powl((params.fPressureDifference)/(params.fFluidDensity),0.5);
	   double T_C = L_C/U_C;

	   double timeFinish = params.fTimeFinish/T_C;
	   double gridTime = params.fGridTime/T_C; 

	   double pressureDifference = 1;

	   double re = (params.fFluidDensity)*(L_C*U_C)/(params.fFluidViscosity);

	   printf("Reinolds number: %lf\n",re);
	   //if (re>1000) {re=1000;}
       //UnSizing (End)

	   TProblem::GlobalLength=L_C;
       TProblem::GlobalVelocity=U_C;
       TProblem::GlobalPressure=params.fPressureDifference;

	   TBargeProblem* bargeProblem = new TBargeProblem(
		                                               timeFinish,gridTime,re,
		                                               areaLength,areaWidth,areaHeight,
													   anchorX,anchorY,anchorZ,
													   bargeWidth,bargeLengthBottom,bargeLengthTop,bargeAngle,
													   pressureDifference,gridSpace,
													   params.fFluidDensity,scour
													   );
	   bargeProblem->Initialize();
	   return bargeProblem;
	 }
	 virtual ~TBargeProblem() {};	 
    protected:  		
		TBargeProblem(
			          double lTimeEndValue, double lTimeStepValue, double lRe, 
			          double lAreaLength, double lAreaWidth, double lAreaHeight, 
					  double lAnchorX, double lAnchorY, double lAnchorZ,
                      double lBargeWidth, double lBargeLengthBottom, double lBargeLengthTop, double lBargeAngle, 
					  double lPressureDifference,double lGridSpace,
					  double lFluidDensity, const TMorphologicScour& lScour
					  ):
		              TScourProblem(lTimeEndValue,lTimeStepValue,lRe,lFluidDensity,lScour),
					  fAreaLength(lAreaLength),fAreaWidth(lAreaWidth),fAreaHeight(lAreaHeight),						                                            																	
					  fAnchorX(lAnchorX),fAnchorY(lAnchorY),fAnchorZ(lAnchorZ),
					  fBargeWidth(lBargeWidth),fBargeLengthBottom(lBargeLengthBottom),fBargeLengthTop(lBargeLengthTop),fBargeAngle(lBargeAngle),
                      fPressureDifference(lPressureDifference),fGridSpace(lGridSpace)
					  {}					                                               
			

		//From TVPProblem (Begin)
		virtual TMaskGenerator* CreateStatement();		
		//From TVPProblem (End)
		
		//From TScourProblem (Begin)
		virtual void PrepareScourFields() 
		{
			for(int i=0; i<Nx; ++i)
			{
			  for(int j=0; j<Ny; ++j)
			  {
					 z_indx[i][j] = 0;
					 h[i][j] = 0;
			  }
			}
		}
		//From TScourProblem (End)


    private:
        double fAreaLength;
	    double fAreaWidth;
	    double fAreaHeight;

	    double fAnchorX;
	    double fAnchorY;
	    double fAnchorZ;
	  
	    double fBargeWidth;
	    double fBargeLengthBottom;
	    double fBargeLengthTop;
	    double fBargeAngle;

		double fPressureDifference;
	    
		double fGridSpace;
	      
        TBargeProblem(const TBargeProblem &);
        TBargeProblem &operator=(const TBargeProblem &);
};


//#################################################################


class TPlatformObstacle: public TObstacle
{
    public:
	   virtual char* GetName() {return NULL;}  

	   double fSizeBottom;
	   double fSizeTop;
	   double fHeightBottom;
	   double fHeightTop;
};

class TPlatformProblem: public TScourProblem
{
    public:	
	 static TPlatformProblem* CreateInstance(const TIniParameters& fIniParameters) 
	 {       
       const TIniParameters& params = fIniParameters;

	   const TPlatformObstacle* platformP = dynamic_cast<const TPlatformObstacle*>(params.fObstacle);
	   const TMorphologicScour* scourP = dynamic_cast<const TMorphologicScour*>(params.fScour);
	   const TPlatformObstacle& platform = *platformP;
	   const TMorphologicScour& scour = *scourP;


	   //UnSizing (Begin)
	   double L_C = params.fAreaLength;

	   double areaLength = params.fAreaLength/L_C;
	   double areaWidth = params.fAreaWidth/L_C;
	   double areaHeight = params.fAreaHeight/L_C;

	   double anchorX = params.fAnchorX/L_C;
	   double anchorY = params.fAnchorY/L_C;
	   double anchorZ = params.fAnchorZ/L_C;
	  
	   double platformSizeBottom = platform.fSizeBottom/L_C;	   
	   double platformSizeTop = platform.fSizeTop/L_C;
	   double platformHeightBottom = platform.fHeightBottom/L_C;
	   double platformHeightTop = platform.fHeightTop/L_C; 

	   double gridSpace = params.fGridSpace/L_C;


	   double U_C = powl(2*(params.fPressureDifference)/(params.fFluidDensity),0.5);
	   double T_C = L_C/U_C;

	   double timeFinish = params.fTimeFinish/T_C;
	   double gridTime = params.fGridTime/T_C; 

	   double pressureDifference = 1;

	   double re = (params.fFluidDensity)*(L_C*U_C)/(params.fFluidViscosity);

	   printf("Reinolds number: %lf\n",re);
	   //if (re>1000) {re=1000;}
       //UnSizing (End)


	   TProblem::GlobalLength=L_C;
       TProblem::GlobalVelocity=U_C;
       TProblem::GlobalPressure=(params.fPressureDifference);

	   TPlatformProblem* platformProblem = new TPlatformProblem(
		                                               timeFinish,gridTime,re,
		                                               areaLength,areaWidth,areaHeight,
													   anchorX,anchorY,anchorZ,
													   platformSizeBottom,platformSizeTop,platformHeightBottom,platformHeightTop,
													   pressureDifference,gridSpace,
													   params.fFluidDensity,scour
													   );
	   platformProblem->Initialize();
	   return platformProblem;
	 }
	 virtual ~TPlatformProblem() {};	 
    protected:  		
		TPlatformProblem(
			          double lTimeEndValue, double lTimeStepValue, double lRe, 
			          double lAreaLength, double lAreaWidth, double lAreaHeight, 
					  double lAnchorX, double lAnchorY, double lAnchorZ,
                      double lPlatformSizeBottom, double lPlatformSizeTop, double lPlatformHeightBottom, double lPlatformHeightTop, 
					  double lPressureDifference,double lGridSpace,
					  double lFluidDensity, const TMorphologicScour& lScour
					  ):
		              TScourProblem(lTimeEndValue,lTimeStepValue,lRe,lFluidDensity,lScour),
					  fAreaLength(lAreaLength),fAreaWidth(lAreaWidth),fAreaHeight(lAreaHeight),						                                            																	
					  fAnchorX(lAnchorX),fAnchorY(lAnchorY),fAnchorZ(lAnchorZ),
					  fPlatformSizeBottom(lPlatformSizeBottom),fPlatformSizeTop(lPlatformSizeTop),fPlatformHeightBottom(lPlatformHeightBottom),fPlatformHeightTop(lPlatformHeightTop),
                      fPressureDifference(lPressureDifference),fGridSpace(lGridSpace)
					  {}					                                               
			

		//From TVPProblem (Begin)
		virtual TMaskGenerator* CreateStatement();		
		//From TVPProblem (End)
		
		//From TScourProblem (Begin)
		virtual void PrepareScourFields() 
		{
			for(int i=0; i<Nx; ++i)
			{
			  for(int j=0; j<Ny; ++j)
			  {
					 z_indx[i][j] = 0;
					 h[i][j] = 0;
			  }
			}
		}
		//From TScourProblem (End)


    private:
        double fAreaLength;
	    double fAreaWidth;
	    double fAreaHeight;

	    double fAnchorX;
	    double fAnchorY;
	    double fAnchorZ;
	  
	    double fPlatformSizeBottom;
	    double fPlatformSizeTop;
	    double fPlatformHeightBottom;
	    double fPlatformHeightTop;

		double fPressureDifference;
	    
		double fGridSpace;
	      
        TPlatformProblem(const TPlatformProblem &);
        TPlatformProblem &operator=(const TPlatformProblem &);
};


//#################################################################


class TBackwardFacingStepProblem: public TVPProblem
{
    public:	
	 static TBackwardFacingStepProblem* CreateInstance(double lTimeEndValue, double lTimeStepValue, double lRe, double lLength, double lWidth, double lHeight, double lStepSize, double lPressureInput, double lPressureOutput, double lSpaceStepValue) 
	 {
       TBackwardFacingStepProblem* bfsp = new TBackwardFacingStepProblem(lTimeEndValue,lTimeStepValue,lRe,lLength,lWidth,lHeight,lStepSize,lPressureInput,lPressureOutput,lSpaceStepValue);
	   bfsp->Initialize();
	   return bfsp;
	 }
	 virtual ~TBackwardFacingStepProblem() {};	 
    protected:  		
		TBackwardFacingStepProblem(double lTimeEndValue, double lTimeStepValue, double lRe, double lLength, double lWidth, double lHeight, double lStepSize, double lPressureInput, double lPressureOutput, double lSpaceStepValue):
		               TVPProblem(lTimeEndValue,lTimeStepValue,lRe),fLength(lLength),fWidth(lWidth),fHeight(lHeight),
						                                            fStepSize(lStepSize), 
																	fPressureInput(lPressureInput),fPressureOutput(lPressureOutput),
																	fStep(lSpaceStepValue) {}					                                               
			

		//From TVPProblem (Begin)
		virtual TMaskGenerator* CreateStatement();		
		//From TVPProblem (End)
		
    private:
      double fLength;
	  double fWidth;
	  double fHeight;

	  double fStepSize;
	  
	  double fStep;

	  double fPressureInput;
	  double fPressureOutput;
      
      TBackwardFacingStepProblem(const TBackwardFacingStepProblem &);
      TBackwardFacingStepProblem &operator=(const TBackwardFacingStepProblem &);
};


//#################################################################

#define CB(x) ((x) * (x) * (x))

class TImmersedBoundaryProblem: public TVPProblem
{
    public:
        TImmersedBoundaryProblem(
                double lTimeEndValue,
                double lTimeStepValue,
                double lRe,
                double lStiffness
                );
				
		TImmersedBoundaryProblem(
                double lTimeEndValue,
                double lTimeStepValue,
                double lmu1,
                double lmu2,
                double lrho1,
                double lrho2,
                double lStiffness
                );
				
        virtual ~TImmersedBoundaryProblem();
        void Initialize();
		void CountBoundaryP(TRnRelease3dSpace &pressure);
		void CountBoundaryConcentration(TRnRelease3dSpace &c);
		void CountBoundaryViscosity(TRnRelease3dSpace &vis);
    protected:
        void GetForce(
                TRnRelease3dSpace &force,
                T3dNormalGrid &gridC,
                T3dNormalGrid &grid,
                TRnRelease3dSpace &velocity,
                const int axis,
                const int timeStepNumber
                );
		double GetDistribution(int, int, int, double, double, double);
		void ChangeStifness(const int timeStepNumber);

		void InterpolateUpdatePathes(TRnRelease3dSpace &velocity, T3dNormalGrid *grid, const int axis, const int timeStepNumber);
        void OutputBoundary(int timeStepNumber);

        double fStiffness;
        TImmersedBoundary *fBoundary;

        static const int COORD_X = 1;
        static const int COORD_Y = 2;
        static const int COORD_Z = 3;

        void ComputeBoundaryForces(const int axis, const int timeStepNumber);
        void SpreadForce(TRnRelease3dSpace &force, T3dNormalGrid &gridC, T3dNormalGrid *grid, const int axis, const int timeStepNumber);
        void Interpolate(TRnRelease3dSpace &velocity, T3dNormalGrid *grid, const int axis, const int timeStepNumber);
        void UpdateBoundaryPosition(const int axis, const int timeStepNumber);

		//==========
		/*
		double *slice_x_0 = new double [101];
		double *slice_y_0 = new double [101];
		double *slice_z_0 = new double [101];


		//*/
		//*
		double *slice_x_0;
		double *slice_y_0;
		double *slice_z_0;
		//*/
		void ComputeChanging(int timeStepNumber);
		void ComputeStress(int timestepNumber);
		void ComputeStress(TRnRelease3dSpace &U, TRnRelease3dSpace &V, TRnRelease3dSpace &W, int timestepNumber);
        //==========
		
		//TNode* UpdateNodeType(TNode *node);

        /*
         * @brief
         *      Return integer index of node by the coordinate.
         *      Assumed temporary that spatial steps is equals
         *      by Ox, Oy, Oz and grid is uniform.
         * @param coord
         *      Coordinate of node by Ox/Oy/Oz
         * @param type
         *      One of COORD_X, COORD_Y, COORD_Z
         * @return 
         *      Index of fluid node, that coordinate is lower than coord by selected axis
        */
        int Index(T3dNormalGrid *grid, const double coord, const int type, const int range);
		int Index(T3dNormalGrid *grid, const double coord, const int type, const int axis, const int range_down, const int range_up);
        /*
         * @brief
         *      Return coordinate of node by the index.
         *      Assumed temporary that spatial steps is equals
         *      by Ox, Oy, Oz and grid is uniform.
         * @param index
         *      Index of node by Ox/Oy/Oz
         * @param type
         *      One of COORD_X, COORD_Y, COORD_Z
	 * @param data_axis
         *      specifies particulat axis, along which data(velocity) was shifted
         *      because of staggered grid.
         * @return 
         *      Coordinate of fluid node
        */
        double Coord(T3dNormalGrid *grid, const int index, const int type);
        double Coord(T3dNormalGrid *grid, const int index, const int type, const int data_axis);
		
        //double ModuleDirichlet(T3dNormalGrid *grid, const double distance, const int type);
        double OptimalDirichlet(const double distance, const double h);
        double Dirichlet(T3dNormalGrid *grid, const double distance, const int type);
};


class TSimpleImmersedValveProblem: public TImmersedBoundaryProblem
{
    public:
        TSimpleImmersedValveProblem(
                double lTimeEndValue,
                double lTimeStepValue,
                double lRe,
                double lStiffness
                );
        virtual ~TSimpleImmersedValveProblem();
        void Initialize();

        void GetForce(
                TRnRelease3dSpace &force,
                T3dNormalGrid &gridC,
                T3dNormalGrid &grid,
                TRnRelease3dSpace &velocity,
                const int axis,
                const int timeStepNumber
                );
		double GetDistribution(int, int, int, double, double, double);
};


class TCylinderBoundaryProblem: public TImmersedBoundaryProblem
{
    public:
        TCylinderBoundaryProblem(
                double lTimeEndValue,
                double lTimeStepValue,
                double lRe,
                double lStiffness
                );
		TCylinderBoundaryProblem(
                double lTimeEndValue,
                double lTimeStepValue,
                double lmu1,
                double lmu2,
                double lrho1,
                double lrho2,
                double lStiffness
                );
        virtual ~TCylinderBoundaryProblem();
        void Initialize();

        void GetForce(
                TRnRelease3dSpace &force,
                T3dNormalGrid &gridC,
                T3dNormalGrid &grid,
                TRnRelease3dSpace &velocity,
                const int axis,
                const int timeStepNumber
                );
		double GetDistribution(int, int, int, double, double, double);
		void ChangeStifness(const int timeStepNumber);
};


class TCanalUnifyProblemWithImmersedBoundary: public TImmersedBoundaryProblem
{
    public:	
        static TCanalUnifyProblemWithImmersedBoundary* CreateInstance(
                double lTimeEndValue,
                double lTimeStepValue,
                double lRe,
                double lStiffness,
                double lRadius,
                double lLength,
                double lPressureInput,
                double lPressureOutput,
                double lSpaceStepValue
                ) 
        {
            TCanalUnifyProblemWithImmersedBoundary* cup = new TCanalUnifyProblemWithImmersedBoundary(
                    lTimeEndValue,
                    lTimeStepValue,
                    lRe,
                    lStiffness,
                    lRadius,
                    lLength,
                    lPressureInput,
                    lPressureOutput,
                    lSpaceStepValue
                    );
            cup->Initialize();
            return cup;
        }
        virtual ~TCanalUnifyProblemWithImmersedBoundary() {};	 
    protected:  		
        TCanalUnifyProblemWithImmersedBoundary(
                double lTimeEndValue,
                double lTimeStepValue,
                double lRe,
                double lStiffness,
                double lRadius,
                double lLength,
                double lPressureInput,
                double lPressureOutput,
                double lSpaceStepValue
                ):
            TImmersedBoundaryProblem(
                    lTimeEndValue,
                    lTimeStepValue,
                    lRe,
                    lStiffness
                    ),
            fRadius(lRadius),fLength(lLength),
            fPressureInput(lPressureInput),fPressureOutput(lPressureOutput),
            fStep(lSpaceStepValue) {}

        //From TVPProblem (Begin)
        virtual TMaskGenerator* CreateStatement();
        //From TVPProblem (End)

    private:
        double fRadius;
        double fLength;

        double fPressureInput;
        double fPressureOutput;

        double fStep;

        TCanalUnifyProblemWithImmersedBoundary(const TCanalUnifyProblemWithImmersedBoundary &);
        TCanalUnifyProblemWithImmersedBoundary &operator=(const TCanalUnifyProblemWithImmersedBoundary &);
};


class TCanalUnifyProblemWithImmersedSimpleValve: public TSimpleImmersedValveProblem
{
    public:	
        static TCanalUnifyProblemWithImmersedSimpleValve* CreateInstance(
                double lTimeEndValue,
                double lTimeStepValue,
                double lRe,
                double lStiffness,
                double lRadius,
                double lLength,
                double lPressureInput,
                double lPressureOutput,
                double lSpaceStepValue
                ) 
        {
            TCanalUnifyProblemWithImmersedSimpleValve* cup = new TCanalUnifyProblemWithImmersedSimpleValve(
                    lTimeEndValue,
                    lTimeStepValue,
                    lRe,
                    lStiffness,
                    lRadius,
                    lLength,
                    lPressureInput,
                    lPressureOutput,
                    lSpaceStepValue
                    );
            cup->Initialize();
            return cup;
        }
        virtual ~TCanalUnifyProblemWithImmersedSimpleValve() {};	 
    protected:  		
        TCanalUnifyProblemWithImmersedSimpleValve(
                double lTimeEndValue,
                double lTimeStepValue,
                double lRe,
                double lStiffness,
                double lRadius,
                double lLength,
                double lPressureInput,
                double lPressureOutput,
                double lSpaceStepValue
                ):
            TSimpleImmersedValveProblem(
                    lTimeEndValue,
                    lTimeStepValue,
                    lRe,
                    lStiffness
                    ),
            fRadius(lRadius),fLength(lLength),
            fPressureInput(lPressureInput),fPressureOutput(lPressureOutput),
            fStep(lSpaceStepValue) {}

        //From TVPProblem (Begin)
        virtual TMaskGenerator* CreateStatement();
        //From TVPProblem (End)

    private:
        double fRadius;
        double fLength;

        double fPressureInput;
        double fPressureOutput;

        double fStep;

        TCanalUnifyProblemWithImmersedSimpleValve(const TCanalUnifyProblemWithImmersedSimpleValve &);
        TCanalUnifyProblemWithImmersedSimpleValve &operator=(const TCanalUnifyProblemWithImmersedSimpleValve &);
};


/*
 * @brief Required space step 0.02, 0.25x1 area, 40000 base stiffness, 900 boundary nodes,
 * 0.15 radius, (0.25, 0.25) center, 10 stiffness coeff by Ox, 0.01 stiffness coeff by Oy
*/
class TCanalUnifyProblemWithImmersedCylinder: public TCylinderBoundaryProblem
{
    public:	
        static TCanalUnifyProblemWithImmersedCylinder* CreateInstance(
                double lTimeEndValue,
                double lTimeStepValue,
                double lRe,
                double lStiffness,
                double lRadius,
                double lLength,
                double lPressureInput,
                double lPressureOutput,
                double lSpaceStepValue
                ) 
        {
            TCanalUnifyProblemWithImmersedCylinder* cup = new TCanalUnifyProblemWithImmersedCylinder(
                    lTimeEndValue,
                    lTimeStepValue,
                    lRe,
                    lStiffness,
                    lRadius,
                    lLength,
                    lPressureInput,
                    lPressureOutput,
                    lSpaceStepValue
                    );
            cup->Initialize();
            return cup;
        }
        virtual ~TCanalUnifyProblemWithImmersedCylinder() {};
    protected:  		
        TCanalUnifyProblemWithImmersedCylinder(
                double lTimeEndValue,
                double lTimeStepValue,
                double lRe,
                double lStiffness,
                double lRadius,
                double lLength,
                double lPressureInput,
                double lPressureOutput,
                double lSpaceStepValue
                ):
            TCylinderBoundaryProblem(
                    lTimeEndValue,
                    lTimeStepValue,
                    lRe,
                    lStiffness
                    ),
            fRadius(lRadius),fLength(lLength),
            fPressureInput(lPressureInput),fPressureOutput(lPressureOutput),
            fStep(lSpaceStepValue) {}

        //From TVPProblem (Begin)
        virtual TMaskGenerator* CreateStatement();
        //From TVPProblem (End)

    private:
        double fRadius;
        double fLength;

        double fPressureInput;
        double fPressureOutput;

        double fStep;

        TCanalUnifyProblemWithImmersedCylinder(const TCanalUnifyProblemWithImmersedCylinder &);
        TCanalUnifyProblemWithImmersedCylinder &operator=(const TCanalUnifyProblemWithImmersedCylinder &);
};

class TCanalIBMWithSourceProblem: public TCylinderBoundaryProblem
{
    public:	
        static TCanalIBMWithSourceProblem* CreateInstance(
                double lTimeEndValue,
                double lTimeStepValue,
                double lRe,
                double lStiffness,
                double lRadius,
                double lLength,
                double lPressureInput,
                double lPressureOutput,
                double lSpaceStepValue
                ) 
        {
            TCanalIBMWithSourceProblem* cup = new TCanalIBMWithSourceProblem(
                    lTimeEndValue,
                    lTimeStepValue,
                    lRe,
                    lStiffness,
                    lRadius,
                    lLength,
                    lPressureInput,
                    lPressureOutput,
                    lSpaceStepValue
                    );
            cup->Initialize();
            return cup;
        }
        virtual ~TCanalIBMWithSourceProblem() {};
    protected:  		
        TCanalIBMWithSourceProblem(
                double lTimeEndValue,
                double lTimeStepValue,
                double lRe,
                double lStiffness,
                double lRadius,
                double lLength,
                double lPressureInput,
                double lPressureOutput,
                double lSpaceStepValue
                ):
            TCylinderBoundaryProblem(
                    lTimeEndValue,
                    lTimeStepValue,
                    lRe,
                    lStiffness
                    ),
            fRadius(lRadius),fLength(lLength),
            fPressureInput(lPressureInput),fPressureOutput(lPressureOutput),
            fStep(lSpaceStepValue) {}

        //From TVPProblem (Begin)
        virtual TMaskGenerator* CreateStatement();
        //From TVPProblem (End)

        double fRadius;
        double fLength;

        double fPressureInput;
        double fPressureOutput;

        double fStep;

        TCanalIBMWithSourceProblem(const TCanalIBMWithSourceProblem &);
        TCanalIBMWithSourceProblem &operator=(const TCanalIBMWithSourceProblem &);

        bool CylinderMask(int i, int j, int k, int N, int L, int M);
        bool InitialDensityDistribution(int i, int j, int k, int N, int L, int M);
        bool InitialConcentrationDistribution(int i, int j, int k, int N, int L, int M);
        bool InitialViscosityDistribution(int i, int j, int k, int N, int L, int M);
        bool ConcentrationInletMask(int i, int j, int k, int N, int L, int M);
        const double ConcentrationInletCondition(int i, int j, int k, int N, int L, int M);
};

class TCanalIBMWithThrombus: public TCylinderBoundaryProblem
{
    public:	
        static TCanalIBMWithThrombus* CreateInstance(
                double lTimeEndValue,
                double lTimeStepValue,
                double lRe,
                double lStiffness,
                double lRadius,
                double lLength,
                double lPressureInput,
                double lPressureOutput,
                double lSpaceStepValue
                ) 
        {
            TCanalIBMWithThrombus* cup = new TCanalIBMWithThrombus(
                    lTimeEndValue,
                    lTimeStepValue,
                    lRe,
                    lStiffness,
                    lRadius,
                    lLength,
                    lPressureInput,
                    lPressureOutput,
                    lSpaceStepValue
                    );
            cup->Initialize();
            return cup;
        }
        virtual ~TCanalIBMWithThrombus() {};
    protected:  		
        TCanalIBMWithThrombus(
                double lTimeEndValue,
                double lTimeStepValue,
                double lRe,
                double lStiffness,
                double lRadius,
                double lLength,
                double lPressureInput,
                double lPressureOutput,
                double lSpaceStepValue
                ):
            TCylinderBoundaryProblem(
                    lTimeEndValue,
                    lTimeStepValue,
                    lRe,
                    lStiffness
                    ),
            fRadius(lRadius),fLength(lLength),
            fPressureInput(lPressureInput),fPressureOutput(lPressureOutput),
            fStep(lSpaceStepValue) {}

        //From TVPProblem (Begin)
        virtual TMaskGenerator* CreateStatement();
        //From TVPProblem (End)

        void Initialize();

    private:
        double fRadius;
        double fLength;

        double fPressureInput;
        double fPressureOutput;

        double fStep;

        TCanalIBMWithThrombus(const TCanalIBMWithThrombus &);
        TCanalIBMWithThrombus &operator=(const TCanalIBMWithThrombus &);

        virtual bool CylinderMask(int i, int j, int k, int N, int L, int M);
        virtual bool InitialDensityDistribution(int i, int j, int k, int N, int L, int M);
        virtual bool InitialConcentrationDistribution(int i, int j, int k, int N, int L, int M);
        virtual bool InitialViscosityDistribution(int i, int j, int k, int N, int L, int M);
        virtual bool ConcentrationInletMask(int i, int j, int k, int N, int L, int M);
        virtual const double ConcentrationInletCondition(int i, int j, int k, int N, int L, int M);
};

class TCanalIBMWithElasticBoundary: public TCylinderBoundaryProblem
{
    public:	
        static TCanalIBMWithElasticBoundary* CreateInstance(
                double lTimeEndValue,
                double lTimeStepValue,
                double lRe,
                double lStiffness,
                double lRadius,
                double lLengthX,
				double lLengthY,
				double lLengthZ,
                double lPressureInput,
                double lPressureOutput,
                double lSpaceStepValue
                ) 
        {
            TCanalIBMWithElasticBoundary* cup = new TCanalIBMWithElasticBoundary(
                    lTimeEndValue,
                    lTimeStepValue,
                    lRe,
                    lStiffness,
                    lRadius,
                    lLengthX,
					lLengthY,
					lLengthZ,
                    lPressureInput,
                    lPressureOutput,
                    lSpaceStepValue
                    );
            cup->Initialize();
            return cup;
        }
		static TCanalIBMWithElasticBoundary* CreateInstance(
                double lTimeEndValue,
                double lTimeStepValue,
				double lmu1,
				double lmu2,
				double lrho1,
				double lrho2,
                double lStiffness,
                double lRadius,
                double lLengthX,
				double lLengthY,
				double lLengthZ,
                double lPressureInput,
                double lPressureOutput,
                double lSpaceStepValue
                ) 
        {
            TCanalIBMWithElasticBoundary* cup = new TCanalIBMWithElasticBoundary(
                    lTimeEndValue,
                    lTimeStepValue,
                    lmu1,
                    lmu2,
                    lrho1,
                    lrho2,
                    lStiffness,
                    lRadius,
                    lLengthX,
					lLengthY,
					lLengthZ,
                    lPressureInput,
                    lPressureOutput,
                    lSpaceStepValue
                    );
            cup->Initialize();
            return cup;
        }
        virtual ~TCanalIBMWithElasticBoundary() {};
    protected:  		
        TCanalIBMWithElasticBoundary(						//base contructor with Re
                double lTimeEndValue,
                double lTimeStepValue,
                double lRe,
                double lStiffness,
                double lRadius,
                double lLengthX,
				double lLengthY,
				double lLengthZ,
                double lPressureInput,
                double lPressureOutput,
                double lSpaceStepValue
                ):
            TCylinderBoundaryProblem(
                    lTimeEndValue,
                    lTimeStepValue,
                    lRe,
                    lStiffness
                    ),
            fRadius(lRadius),
			fLengthX(lLengthX), fLengthY(lLengthY), fLengthZ(lLengthZ), 
            fPressureInput(lPressureInput),fPressureOutput(lPressureOutput),
            fStep(lSpaceStepValue) {};
		
		TCanalIBMWithElasticBoundary(						//contructor with changing mu and rho
				double lTimeEndValue,
                double lTimeStepValue,
                double lmu1,					//mu1
                double lmu2,					//mu2
                double lrho1,					//rho1
                double lrho2,					//rho2
                double lStiffness,
                double lRadius,
                double lLengthX,
				double lLengthY,
				double lLengthZ,
                double lPressureInput,
                double lPressureOutput,
                double lSpaceStepValue
                ):
            TCylinderBoundaryProblem(
                    lTimeEndValue,
                    lTimeStepValue,
                    lmu1,
                    lmu2,
                    lrho1,
                    lrho2,
                    lStiffness
                    ),
            fRadius(lRadius),
			fLengthX(lLengthX), fLengthY(lLengthY), fLengthZ(lLengthZ), 
            fPressureInput(lPressureInput),fPressureOutput(lPressureOutput),
            fStep(lSpaceStepValue) {};
		
        //From TVPProblem (Begin)
        virtual TMaskGenerator* CreateStatement();
        //From TVPProblem (End)

        void Initialize();

        //double fRadius;
        double fLengthX;
		double fLengthY;
		double fLengthZ;
		double fRadius;
        double fPressureInput;
        double fPressureOutput;

        double fStep;

        TCanalIBMWithElasticBoundary(const TCanalIBMWithElasticBoundary &);
        TCanalIBMWithElasticBoundary &operator=(const TCanalIBMWithElasticBoundary &);

        bool CylinderMask(int i, int j, int k, int N, int L, int M);
        bool InitialDensityDistribution(int i, int j, int k, int N, int L, int M);
        bool InitialConcentrationDistribution(int i, int j, int k, int N, int L, int M);
        bool InitialViscosityDistribution(int i, int j, int k, int N, int L, int M);
        bool ConcentrationInletMask(int i, int j, int k, int N, int L, int M);
        const double ConcentrationInletCondition(int i, int j, int k, int N, int L, int M);
		double GetLengthX() {return this->fLengthX;}
		double GetLengthY() {return this->fLengthY;}
		double GetLengthZ() {return this->fLengthZ;}
};


#endif
