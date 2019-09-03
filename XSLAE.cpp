#pragma hdrstop

#include <stdio.h>
#include <math.h>
#include "XSLAE.h"

//#################################################################

//*************************************
TSLAE::TSLAE(const TRnLinearOperator &lA,TRnSpace &lx,TRnSpace &lb):fA(&lA), fx(&lx), fb(&lb), Result(this)
{
    fResult = &lx.Copy();
	*fResult|= 0;
}
TSLAE::TSLAE(const TRnLinearOperator &lA,TRnSpace &lx, TRnSpace &lu, TRnSpace &lv, TRnSpace &lw, TRnSpace &lb):fA(&lA), fx(&lx), fu(&lu), fv(&lv), fw(&lw), fb(&lb), Result(this)
{
    fResult = &lx.Copy();
	*fResult|= 0;
}
//*************************************
TSLAE::~TSLAE()
{
	if(fA != NULL)	{printf("deleting fA\n");	delete this->fA;}
	if(fx != NULL)  {printf("deleting fx\n"); 	delete this->fx;}
	if(fu != NULL)	{printf("deleting fu\n");	delete this->fu;}
	if(fv != NULL)	{printf("deleting fv\n");	delete this->fv;}
	if(fw != NULL)	{printf("deleting fw\n");	delete this->fw;}
	if(fb != NULL)  {printf("deleting fb\n");	delete this->fb;}
	if(fResult != NULL)   {printf("deleting fResult\n");delete this->fResult;}
}
//*************************************
const TRnSpace &TSLAE::GetResult() const
{
    return *fResult;
}
//*************************************
void TSLAE::UpdateP(TRnSpace &lp)
{
	fx->Assign(lp);
}
	 
//*************************************
void TSLAE::UpdateRP(TRnSpace &lrp)
{
	fb->Assign(lrp);
}
	 
//*************************************
void TSLAE::UpdateA(TRnLinearOperator &lA)
{
	throw "TSLAE::UpdateA too lazy to make this!";
}
//*************************************
//*************************************

void TSLAE::SolveByConjugateGradients()
{
          const TRnLinearOperator &A = *fA;
          TRnSpace &x = *fx;
          TRnSpace &b = *fb;
          TRnSpace &result = *fResult;


          double eps = 1e-5;

          TRnSpace &xk_1 = x.Copy();		

          TRnSpace &rk_1 = x.CreateInstance();
          rk_1|= A(xk_1) - b;
          double rk_1Xrk_1 = (rk_1,rk_1);
          result|= xk_1;


		  printf("I:");
		  printf("%d",0);
	      printf("   ");
	      printf("R:");
          printf("%20.8lf ", (double) sqrtl(rk_1Xrk_1));
		  printf("\n");
		  //_getch();


          if (rk_1Xrk_1!=0)
          {

          double tauk = rk_1Xrk_1/(A(rk_1),rk_1);


          TRnSpace &xk = x.Copy();
          xk|= xk_1 - tauk*rk_1;


          TRnSpace &rk = x.CreateInstance();
          rk|= A(xk) - b;
          double rkXrk = (rk,rk);
          result|= xk;


          if (rkXrk!=0)
          {


          double alfak = 1;

          int iteration = 0;



          double epsCur = sqrtl(rkXrk);

          double tauk1;
          double alfak1;

          while (epsCur>eps)
          {
             
             iteration = iteration+1;

             tauk1 = rkXrk/(A(rk),rk);
             alfak1 = 1-(1/alfak)*(tauk1/tauk)*(rkXrk/rk_1Xrk_1);
             alfak1 = 1/alfak1;


             result|= alfak1*xk+(1-alfak1)*xk_1;
             result|= result - (tauk1*alfak1)*rk;

             xk_1|= xk;
             rk_1|= rk;
             rk_1Xrk_1 = rkXrk;
             tauk = tauk1;
             alfak = alfak1;


             xk|= result;
             rk|= A(xk) - b;
             rkXrk = (rk,rk);


             epsCur = sqrtl(rkXrk);

			 if (iteration%1000==0)
			 {
			  printf("I:");
		      printf("%d",iteration);
			  printf("\n");
			  printf("R:");
              printf("%20.8lf ", (double) epsCur);
			  printf("\n");
			 }
          }
          }
          delete &xk;
          delete &rk;
          }     
          delete &xk_1;
          delete &rk_1;
}
//*************************************
void TSLAE::SolveByBiConjugateGradients()
{
	const TRnLinearOperator &A = *fA;
	TRnSpace &x = *fx;
	TRnSpace &b = *fb;
	TRnSpace &result = *fResult;

	double eps = 1e-8;

	
	TRnSpace &xk = x.Copy();		
	TRnSpace &xk_new = x.Copy();		
	xk|= 0;

	
	TRnSpace &rk = x.CreateInstance();
	TRnSpace &rt = x.CreateInstance();
	TRnSpace &rk_new = x.CreateInstance();
	TRnSpace &rt_new = x.CreateInstance();
	rk|= A(xk) - b;
	rt|=(-1)*rk;
	rk|= rt;

	
	TRnSpace &pk = x.CreateInstance();
	TRnSpace &pt = x.CreateInstance();
	TRnSpace &pk_new = x.CreateInstance();
	TRnSpace &pt_new = x.CreateInstance();
	pk|= rk;
	pt|= rt;


	double rkXrk = (rk,rk);
	result|= xk;

	
	double alpha, beta;
	int iteration = 0;
	double epsCur = sqrtl(rkXrk);


	printf("BCG metod: Iteration = %d,  Resudial = %20.8lf\n", iteration, (double) epsCur);	

	double alpha1, beta1;
	double alpha2, beta2;

	const TRnLinearOperator &A_ = A.Copy();
	~A_;

	#pragma omp parallel num_threads(5)
	{	    
	while (epsCur>eps)
	{
		#pragma omp master
        {
		  ++iteration;
		}  

		#pragma omp sections
		{
		   #pragma omp section
           {
                alpha1 = (rk, rt);    
           }
           #pragma omp section
           {
                alpha2 = (A(pk), pt);     
           }   
		}
		
		#pragma omp barrier

		#pragma omp master
		{
		  alpha =  alpha1/alpha2; 		  
		}
		
		#pragma omp barrier
		
		#pragma omp sections
		{
		   #pragma omp section
           {
			 xk_new|= xk + alpha*pk;
		   }
		   #pragma omp section
           {         
		     rk_new|= A(pk);
		     rk_new|=(-alpha)*rk_new + rk;
		   }
		   #pragma omp section
           {         		   
		     rt_new|= A_(pt);
		     rt_new|=(-alpha)*rt_new + rt;
		   }
		}
	
		#pragma omp barrier

		#pragma omp sections
		{
		   #pragma omp section
           {
                beta1 = (rk_new, rt_new);    
           }
           #pragma omp section
           {
                beta2 = (rk, rt);     
           }   
		}

		#pragma omp barrier

		#pragma omp master
		{		   
		   if(beta1==0) {throw "Error in TSLAE::SolveByBiConjugateGradients: Beta1 is zero.";}
		   else if(beta2==0) {throw "Error in TSLAE::SolveByBiConjugateGradients: Beta2 is zero.";}
		   beta =  beta1/beta2;		   
		} 
        
		#pragma omp barrier

		#pragma omp sections
		{
		   #pragma omp section
           {
		     pk_new|= rk_new + beta*pk;
		   }	 
           #pragma omp section
           {
		     pt_new|= rt_new + beta*pt;
		   }
		} 
		
		#pragma omp barrier
		
		#pragma omp sections
		{
		   #pragma omp section
           {
		     xk|= xk_new;
		   }
		   #pragma omp section
           {
		     rk|= rk_new;
		   }
		   #pragma omp section
           {
		     rt|= rt_new;
		   }
		   #pragma omp section
           {
		     pk|= pk_new;
		   }
		   #pragma omp section
           {
		     pt|= pt_new;
		   }
		}

		#pragma omp barrier

		#pragma omp master
        {            
		  rkXrk = (rk,rk);
		  epsCur = sqrtl(rkXrk);       
		  if (iteration%1000==0) {printf("BCG metod: Iteration = %d,  Resudial = %20.8lf\n", iteration, (double) epsCur);}
		}

		#pragma omp barrier
	} //while (epsCur>eps)	
	}//#pragma omp parallel

	printf("BCG metod: Iteration = %d,  Resudial = %20.8lf\n", iteration, (double) epsCur);		

	result|= xk;

	delete &xk;
	delete &rk;
	delete &rt;
	delete &pk;
	delete &pt;
	delete &xk_new;
	delete &rk_new;
	delete &rt_new;
	delete &pk_new;
	delete &pt_new;

	delete &A_;
}
//*************************************
void TSLAE::SolveByMinimalResiduals()
{        
	const TRnLinearOperator &A = *fA;
    TRnSpace &x = *fx;
    TRnSpace &b = *fb;
    TRnSpace &result = *fResult;

          
    double eps = 0.0000001;

	TRnSpace &xk = x.Copy();

    TRnSpace &rk = x.CreateInstance();
    rk|= A(xk) - b;
    double rkXrk = (rk,rk);
    result|= xk;

    TRnSpace &Ark = x.CreateInstance();

    double alfa;


    if (rkXrk!=0)
    {

		double epsCur = sqrtl(rkXrk);
		int iteration = 0;          

		printf("Minimal Residual method: Iteration = %d Residual = %8.20lf \n", iteration, epsCur);

		while (epsCur>eps)
		{
			iteration = iteration+1;

			Ark|=A(rk);
			alfa = (-1)*(Ark,rk)/(Ark,Ark);

			xk|=xk+alfa*rk;
			
			rk|=A(xk) - b;

			rkXrk = (rk,rk);
			epsCur = sqrtl(rkXrk);             

			if (iteration%1000==0)
			{
				printf("Minimal Residual method: Iteration = %d Residual = %8.20lf \n", iteration, epsCur);
				/*printf("%d",iteration);
				printf("\t");
				printf("R:");
				printf("%8.20f",epsCur);
				printf("\n");
				getchar();*/
			}
		}
		printf("Minimal Residual method: Iteration = %d Residual = %8.20lf \n", iteration, epsCur);
	}


    result|=xk;

    delete &xk;
    delete &rk;
    delete &Ark;        
}
//*************************************
void TSLAE::SolveByScalarRunning(int lM, double* lU, const double* lA, const double* lB, const double* lC, const double* lF)
{
   int M = lM; 
   double* u = lU;
   const double* a = lA;
   const double* b= lB;
   const double* c = lC;
   const double* f = lF;
   
   double *A = new double[M+1];
   double *B = new double[M+1];   
   
   A[0] = 0;
   B[0] = u[0];
   A[M] = B[M] = 0;

   for (int m=1; m<M; m++)
   {
      A[m] = (-c[m])/(a[m]*A[m-1] + b[m]) ;
      B[m] = (f[m] - a[m]*B[m-1])/(a[m]*A[m-1] + b[m]);
   }

   for (int m=M-1; m>0; m--)
   {
      u[m] = A[m]*u[m+1] + B[m];
   }

   delete[] A;
   delete[] B;         
}
//*************************************
void TSLAE::SolveByGauss(const double * const * const lA,const double * const lf,double * const lres,const int lN)
{
       int i; int j;

       double ** A = new double *[lN];
       double *f = new double [lN];
       for (i=0;i<lN;i++)
        {
           f[i] = lf[i];
           A[i] = new double [lN];
           for (j=0;j<lN;j++)
            {
              A[i][j] = lA[i][j];
            }
        }

       double *res = lres;

       int N = lN-1;
       int step;
       double tempValue; double tempValue1;
       double maxElement; int numOfmaxElement;

       double resi;


	   //Direct process (Begin)
       for(step=0;step<N;step++) 
       {
                   //Find max. element and rows exchange (Begin)
                   maxElement = fabs(A[step][step]);
                   numOfmaxElement = step;
                   for(i=step;i<N+1;i++)
                   {
                       if (fabs(A[i][step])>maxElement)
                       {
                          maxElement = fabs(A[i][step]);
                          numOfmaxElement = i;
                       }
                   }
                   if (maxElement==0)
                   {
                      throw "Error in TSLAE::SolveByGauss: The max. element in current column = 0";
                   }
                   if (numOfmaxElement!=step)
                   {
                        for(j=step;j<N+1;j++)
                        {
                            tempValue = A[step][j];
                            A[step][j] = A[numOfmaxElement][j];
                            A[numOfmaxElement][j] = tempValue;
                        }
                        tempValue = f[step];
                        f[step] = f[numOfmaxElement];
                        f[numOfmaxElement] = tempValue;
                   } 
                   //Find max. element and rows exchange (End)

                   //Combination - to zero (Begin) 
                   tempValue = A[step][step];
                   for(i=step+1;i<N+1;i++)
                   {
                       tempValue1 = A[i][step];
                       if (tempValue1!=0)
                       {
                           for(j=step;j<N+1;j++)
                           {
                               A[i][j] = A[i][j] - (A[step][j]*tempValue1/tempValue);
                           }
                           f[i] = f[i]-(f[step]*tempValue1/tempValue);
                       }
                   } 
				   //Combination - to zero (End) 
       } 
       //Direct process (End)

       //Diagonal analysis (Begin)
	   double max=fabs(A[0][0]);
       double min=fabs(A[0][0]);
       int maxI = 0;
       int minI = 0;
       for (i=0;i<lN;i++)
       {
          if ( fabs(A[i][i]) > max ) {max = fabs(A[i][i]); maxI = i;}
          if ( fabs(A[i][i]) < min ) {min = fabs(A[i][i]); minI = i;}
       }                
       //Diagonal analysis (End)

       //Compatible analysis (Begin)
       double *constVec = new double [lN];
       double *resConstVec = new double [lN];
       
	   for (i=0;i<lN;i++)
       {
          constVec[i] = 5; //for example
       }

       for (i=0;i<lN;i++)
       {
          resConstVec[i] = 0;
          for (j=0;j<lN;j++)
          {
             resConstVec[i] = resConstVec[i] + lA[i][j]*constVec[j];
          }
       }

       double ker = 0;
       for (i=0;i<lN;i++)
       {
           ker = ker + fabs(resConstVec[i]);
       }       
                
	   double comp = 0;
       for (i=0;i<lN;i++)
       {
          comp = comp + lf[i]*constVec[i];
       }
       //Compatible analysis (Begin)

	   delete[] constVec;
	   delete[] resConstVec;

       //Back process (Begin)
       for(i=N;i>-1;i--)
       {
          resi = f[i];
          for(j=i+1;j<N+1;j++)
          {
             resi = resi-(A[i][j]*res[j]);
          }
          res[i] = resi/A[i][i];
        }
        //Back process (End)


       delete[] f;
       for (i=0;i<lN;i++)
        {
           delete[] A[i];
        }
       delete[] A;
}
//*************************************
void TSLAE::SolveByGauss()
{       
       const TRnLinearOperator &A = *fA;
       if (A.Matrix==NULL)
        {
          throw "Error in TSLAE::SolveByGauss: The operator matrix is not ready.";
        }
       TRnSpace &x = *fx;
       TRnSpace &b = *fb;
       TRnSpace &result = *fResult;

       int N = x.EndIndex;

       const double * const * const GaussA = A.Matrix;

       double *GaussRes = new double [N];
       double *GaussCheck = new double [N];
      
       for (int i=1;i<N+1;i++)
        {
              GaussRes[i-1] = 0;
              GaussCheck[i-1] = 0;
        }

       double *GaussF = new double [N];
       
	   //Constructing right part (Begin)
       TRnSpace &helpB = x.CreateInstance();
       TRnSpace &helpZero = x.Copy();

       helpZero|=0;
       helpB|=A(helpZero) - b;
       for (int i=1;i<N+1;i++)
        {
          GaussF[i-1] = (-1)*helpB[i];
        }

       delete &helpB;
       delete &helpZero;
       //Constructing right part (End)

       const double * const GaussFConst = GaussF;
       double * const GaussResConst = GaussRes;
       SolveByGauss(GaussA,GaussFConst,GaussResConst,N); 


       
	   //Checking solution (Begin)
       /*
       double ch=0;

       for (i=0;i<N;i++)
       {
         for (j=0;j<N;j++)
         {
           GaussCheck[i] = GaussCheck[i] + GaussA[i][j]*GaussRes[j];
         }

         GaussCheck[i] = GaussCheck[i] - GaussF[i];
         ch = ch + fabs(GaussCheck[i]);
       }
       */
       //Checking solution (End)


       for (int i=1;i<N+1;i++)
        {
           result[i] = GaussRes[i-1];
        }

  
       delete[] GaussF;
       delete[] GaussRes;
       delete[] GaussCheck;       
}
//*************************************

void TSLAE::SolveByBiConjugateGradientsStab()
{
	double eps = 1e-5;

	const TRnLinearOperator &A = *fA;
	TRnSpace &x = *fx;
	TRnSpace &b = *fb;
	TRnSpace &result = *fResult;

	TRnSpace &Un = x.Copy();

	TRnSpace &tmp1 = x.CreateInstance();
	TRnSpace &tmp2 = x.CreateInstance();
	TRnSpace &Rn = x.CreateInstance();
	Rn |= A(Un);
	Rn |= b - Rn;

	// для векторов реализованы операции сложения/вычитания и умножения/деления на скаляр

	TRnSpace &sRn = x.CreateInstance();
	sRn |= Rn;

	
	TRnSpace &Pn = x.CreateInstance();
	TRnSpace &Pn1 = x.CreateInstance();
	TRnSpace &Vn = x.CreateInstance();
	TRnSpace &Sn = x.CreateInstance();
	TRnSpace &Tn = x.CreateInstance();

	// заполнение вектора нулями
	Pn |= 0;
	Vn |= 0;

	double pn = 1.;
	double an = 1.;
	double bn = 1.;
	double wn = 1.;
	
	double pn1;

	

	int iteration = 0;
	//double norm = sqrt((Rn, Rn));
	//double norm_0 = 1.0 / (Rn,Rn);
	double norm = Rn.GetAbsMaxElement();

    //eps = pow(eps,2);
	printf("BCGstab method: Iteration = %d,  Resudial = %22.10lf\n", iteration, norm);
	while ( norm > eps)
	{
		// дальше алгоритм с википедии (модифицированный)
		iteration++;
		pn1 = (sRn, Rn);
		bn = (pn1 / pn) * (an / wn);

		Pn1 |= wn*Vn;
		Pn1 |= Pn - Pn1;

		Pn1 |= bn*Pn1;
		Pn1 |= Rn + Pn1;

		Vn |= A(Pn1);
		an = pn1 / (sRn, Vn);

		Sn |= an*Vn;
		Sn |= Rn - Sn;

		Tn |= A(Sn);
		wn = (Tn, Sn) / (Tn, Tn);

		tmp1 |= an*Pn1;
		tmp2 |= tmp1 + wn*Sn;

		Un |= Un + tmp2;

		Rn |= wn*Tn;
		Rn |= Sn - Rn;

		pn = pn1;
		Pn |= Pn1;
		
		norm = Rn.GetAbsMaxElement();

		if(iteration%1000 == 0)
			printf("BCGstab method: Iteration = %d Residual = %22.10lf \n", iteration, norm);
	}
	result |= Un;
	printf("BCGstab method: Iteration = %d,  Resudial = %22.10lf\n", iteration, norm);

	delete &Un;
	delete &Rn;
	delete &sRn;
	delete &Pn;
	delete &Pn1;
	delete &Vn;
	delete &Sn;
	delete &Tn;
	delete &tmp1;
	delete &tmp2;
}
//#################################################################
