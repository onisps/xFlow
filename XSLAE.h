#ifndef SLAEH
#define SLAEH

#include "XOperator.h"
#include "property.h"

//#################################################################

class TSLAE
{
  public:
     TSLAE(const TRnLinearOperator &lA,TRnSpace &lx,TRnSpace &lb);
	 TSLAE(const TRnLinearOperator &lA,TRnSpace &lx, TRnSpace &lu, TRnSpace &lv, TRnSpace &lw, TRnSpace &lb);
     virtual ~TSLAE();
   
     virtual const TRnSpace &GetResult() const;       
     void SetResult(TRnSpace const& Result) {
         throw "Not implemented";
     }
     Property<TSLAE, const TRnSpace &, &TSLAE::GetResult, &TSLAE::SetResult> Result;
	 //__declspec(property(get = GetResult)) const TRnSpace &Result;

	 static void SolveByScalarRunning(int lM, double* lU, const double* lA, const double* lB, const double* lC, const double* lF);
	 static void SolveByGauss(const double * const * const lA,const double * const lf,double * const lres,const int lN);

	 void SolveByGauss();

	 void UpdateP(TRnSpace &lp);
	 void UpdateRP(TRnSpace &lrp);
	 void UpdateA(TRnLinearOperator &lA);
	 
	 void SolveByMinimalResiduals();
	 void SolveByConjugateGradients();
	 void SolveByBiConjugateGradients();
	 void SolveByBiConjugateGradientsStab();

  protected:     
      const TRnLinearOperator *fA = NULL;
      TRnSpace *fx = NULL;
	  TRnSpace *fu = NULL;
	  TRnSpace *fv = NULL;
	  TRnSpace *fw = NULL;
      TRnSpace *fb = NULL;
      TRnSpace *fResult = NULL;      
  private:
      TSLAE(const TSLAE &);
      TSLAE& operator=(const TSLAE &);
};

//#################################################################

#endif
