#pragma hdrstop

#include "XMask.h"

//#################################################################


//*************************************
  T1dNumberMask::T1dNumberMask(const int lN)
                              : TMask(), fN(lN), N(this), Mask(this)
    {
       fMask = new int[fN];
       for (int i=0;i<fN;i++)
        {
          fMask[i] = 0;
        }

    }
//*************************************
  T1dNumberMask::~T1dNumberMask()
   {
	  delete[] fMask;// fMask = NULL;
   }
//*************************************
  int &T1dNumberMask::operator[](const int i) const
   {
     return fMask[i];
   }
//*************************************
  T1dNumberMask &T1dNumberMask::Copy() const
   {
      T1dNumberMask &res = *(new T1dNumberMask(fN));
      for (int i=0;i<fN;i++)
       {
          res.fMask[i] = fMask[i];
       }
      return res;
   }
//*************************************


//#################################################################


//*************************************
  T2dNumberMask::T2dNumberMask(const int lN, const int lL)
                              : TMask(), fN(lN), fL(lL), N(this), L(this), Mask(this)
    {
       fMask = new int*[fN];
       for (int i=0;i<fN;i++)
        {
          fMask[i] = new int[fL];
          for (int j=0;j<fL;j++)
           {
            fMask[i][j] = 0;
           } 
        }

    }
//*************************************
  T2dNumberMask::~T2dNumberMask()
   {
      for (int i=0;i<fN;i++)
        {
          delete[] fMask[i];
        }
	  delete[] fMask;// fMask = NULL;
   }
//*************************************
  int * T2dNumberMask::operator[](const int i) const
   {
     return fMask[i];
   }
//*************************************
  T2dNumberMask &T2dNumberMask::Copy() const
   {
      T2dNumberMask &res = *(new T2dNumberMask(fN,fL));
      for (int i=0;i<fN;i++)
       {
          for (int j=0;j<fL;j++)
           {
              res.fMask[i][j] = fMask[i][j];
           }
       }
      return res;
   }
//*************************************


//#################################################################


//*************************************
  T3dNumberMask::T3dNumberMask(const int lN, const int lL, const int lM)
                              : TMask(), fN(lN), fL(lL), fM(lM), N(this), L(this), M(this), Mask(this)
    {
       fMask = new int**[fN];
       for (int i=0;i<fN;i++)
        {
          fMask[i] = new int*[fL];
          for (int j=0;j<fL;j++)
           {
            fMask[i][j] = new int[fM];
            for (int k=0;k<fM;k++)
             {
               fMask[i][j][k] = 0;
             }
           } 
        }

    }
//*************************************
  T3dNumberMask::~T3dNumberMask()
   {
      for (int i=0;i<fN;i++)
        {
          for (int j=0;j<fL;j++)
           {
             delete[] fMask[i][j];
           }
           delete[] fMask[i];
        }
	  delete[] fMask;// fMask = NULL;
   }
//*************************************
  int * const * T3dNumberMask::operator[](const int i) const
   {
     return fMask[i];
   }
//*************************************
  T3dNumberMask &T3dNumberMask::Copy() const
   {
      T3dNumberMask &res = *(new T3dNumberMask(fN,fL,fM));
      for (int i=0;i<fN;i++)
       {
          for (int j=0;j<fL;j++)
           {
              for (int k=0;k<fM;k++)
                {
                  res.fMask[i][j][k] = fMask[i][j][k];
                }
           }
       }
      return res;
   }
//*************************************


//#################################################################
