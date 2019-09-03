#pragma hdrstop

#include <stdio.h>
#include <math.h>
#include "XSeparator.h"

//#################################################################

//*************************************
TSeparator::TSeparator(const double lStartValue, const double lEndValue)
                           : fStartValue(lStartValue), fEndValue(lEndValue),
                           StartValue(this), EndValue(this), EndIndex(this),
                           CountNodes(this), Dimension(this), Separation(this)
       {
         if (fStartValue>=fEndValue) {throw "Error in TSeparator constructor: Start value is equal or more then end value.";}
         fEndIndex = 0;
         fDimension = NULL;
         fSeparation = NULL;
       }
//*************************************

//#################################################################

//*************************************
TUniformSeparator::TUniformSeparator(const double lStartValue, const double lEndValue, const double lSeparationValue)
                                         : TSeparator(lStartValue, lEndValue), fSeparationValue(lSeparationValue), SeparationValue(this)
       {
          /*
          Не всегда возможно разделить отрезок на равные целые части с заданным шагом.
          Поэтому скорректируем шаг в сторону уменьшения (если необходимо) для
          достижения этой цели.
          */

          if (fSeparationValue<=0) {throw "Error in TUniformSeparator constructor: Separation value is equal or less then zero.";}


          double resDev = (fEndValue-fStartValue)/fSeparationValue;
          double iPart = 0;

          modf(resDev,&iPart);

          fEndIndex = (int)iPart;
		  double eps = 0.00000000001;
          if ( fabs(resDev-iPart)>eps ) {fEndIndex = fEndIndex + 1;}

          fSeparationValue = (fEndValue-fStartValue)/((double)fEndIndex);

          int i;

          fDimension = new double[fEndIndex+1];
          fSeparation = new double[fEndIndex];

          fDimension[0] = fStartValue;
          fDimension[fEndIndex] = fEndValue;

          fSeparation[0] = SeparationValue;

          for (i=1;i<fEndIndex;i++)
           {
            fDimension[i] = fDimension[i-1] + fSeparationValue;
            fSeparation[i] = SeparationValue;
           }

       }
//*************************************
TUniformSeparator::~TUniformSeparator()
      {
        delete[] fDimension;
        delete[] fSeparation;
      }
//*************************************
TUniformSeparator &TUniformSeparator::Copy() const
      {
         return *(new TUniformSeparator(fStartValue, fEndValue, fSeparationValue));
      }
//*************************************

//#################################################################

//*************************************
TRandomSeparator::TRandomSeparator(int lCountNodes, const double *lDimension)
                                         : TSeparator(lDimension[0], lDimension[lCountNodes-1])
       {          
          fEndIndex = lCountNodes -1;
          
          fDimension = new double[fEndIndex+1];
          fSeparation = new double[fEndIndex];

          fDimension[0] = fStartValue;
          fDimension[fEndIndex] = fEndValue;

          fSeparation[0] = lDimension[1]-lDimension[0];

          for (int i=1;i<fEndIndex;i++)
           {
            fDimension[i] = lDimension[i];
            fSeparation[i] = lDimension[i+1]-lDimension[i];
           }
       }
//*************************************
TRandomSeparator::~TRandomSeparator()
      {
        delete[] fDimension;
        delete[] fSeparation;
      }
//*************************************
TRandomSeparator &TRandomSeparator::Copy() const
      {
         return *(new TRandomSeparator(GetCountNodes(),fDimension));
      }
//*************************************

//#################################################################
