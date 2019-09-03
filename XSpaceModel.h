#ifndef SpacesModelsH
#define SpacesModelsH

#include "XMask.h"
#include "property.h"
#include <cstddef>
#include <stdio.h>
#include <math.h>

//------------------Abstract (Begin)------------------------------------------


//#################################################################

class TCorrespondence
{
       public:
        TCorrespondence() {}
        virtual ~TCorrespondence() {}
       protected:
        virtual TCorrespondence &Copy() const = 0;
       private:
        TCorrespondence(const TCorrespondence &) {}
        TCorrespondence &operator=(const TCorrespondence &) {return *this;}
};

//#################################################################

class TInverseCorrespondence
     {
       public:
        TInverseCorrespondence() {}
        virtual ~TInverseCorrespondence() {}
       protected:
        virtual TInverseCorrespondence &Copy() const = 0;
       private:
        TInverseCorrespondence(const TInverseCorrespondence &) {}
        TInverseCorrespondence &operator=(const TInverseCorrespondence &) {return *this;}
     };

//#################################################################     

class TRnSpaceModel
{
      public:
       TRnSpaceModel(): EndIndex(this), InverseCorrespondence(this) {}
       virtual ~TRnSpaceModel() {}
       virtual TRnSpaceModel &Copy() = 0;
	   
	   int GetEndIndex() const {return fEndIndex;}
	   void SetEndIndex(const int& N) {
           throw "Not implemented";
       }
       Property<TRnSpaceModel, int, &TRnSpaceModel::GetEndIndex, &TRnSpaceModel::SetEndIndex> EndIndex;

       const TCorrespondence & GetCorrespondence(const int lIndex) const
        {
          if ( (lIndex<1) || (lIndex>fEndIndex) )
           {
             throw "Error in TRnSpaceModel::GetCorrespondence: The index is out of range.";
           }
          return  *(fCorrespondence[lIndex-1]);
        }
       const TInverseCorrespondence & GetInverseCorrespondence() const
        {
          return *fInverseCorrespondence;
        }

       void SetInverseCorrespondence(const TInverseCorrespondence& InverseCorrespondence){
           throw "Not implemented";
       }
       Property<TRnSpaceModel,
           const TInverseCorrespondence &,
           &TRnSpaceModel::GetInverseCorrespondence,
           &TRnSpaceModel::SetInverseCorrespondence> InverseCorrespondence;
      protected:       
       int fEndIndex;
       TCorrespondence **fCorrespondence;
       TInverseCorrespondence *fInverseCorrespondence;
      private:
       TRnSpaceModel(TRnSpaceModel &): EndIndex(this), InverseCorrespondence(this) {}
       TRnSpaceModel &operator=(TRnSpaceModel &) {return *this;}
};

//#################################################################     


//---------------------Abstract (End)-----------------------------------------







//---------------------MultiData (Begin)-----------------------------------------


//#################################################################     

class TMultiDataCorrespondence: public TCorrespondence
{
       public:
        TMultiDataCorrespondence(const int lInd): TCorrespondence(), fInd(lInd), Ind(this)
         {
             if (fInd<0)
              {
                throw "MultiDataCorrespondence constructor error: The index value is incorrect.";
              }
         }

		int GetInd() const {return fInd;}
        void SetInd(const int& Ind) {
            throw "Not implemented";
        }
        Property<TMultiDataCorrespondence,
            int,
            &TMultiDataCorrespondence::GetInd,
            &TMultiDataCorrespondence::SetInd> Ind;
       protected:
        const int fInd;
       private:
        TMultiDataCorrespondence(const TMultiDataCorrespondence &):fInd(0), Ind(this) {}
        TMultiDataCorrespondence &operator=(const TMultiDataCorrespondence &) {return *this;}
};

//#################################################################     

class TMultiDataInverseCorrespondence: public TInverseCorrespondence
{
       public:
        TMultiDataInverseCorrespondence(const int lDataCount):TInverseCorrespondence(), fDataCount(lDataCount)
          {
            if (fDataCount<=0)
              {
                throw "TMultiDataInverseCorrespondence constructor error: The data count value is incorrect.";
              }
          }
        virtual ~TMultiDataInverseCorrespondence() {};
       protected:
        const int fDataCount;
       private:
        TMultiDataInverseCorrespondence(const TMultiDataInverseCorrespondence &): fDataCount(0) {}
        TMultiDataInverseCorrespondence &operator=(const TMultiDataInverseCorrespondence &) {return *this;}
};

//#################################################################     

class TRnReleaseMultiDataSpaceModel: public TRnSpaceModel
{
      public:
        TRnReleaseMultiDataSpaceModel(const int lDataCount);
        TRnReleaseMultiDataSpaceModel(const int lDataCount, const TMask * const * lMasks);
        virtual ~TRnReleaseMultiDataSpaceModel();

		int GetDataCount() const {return fDataCount;}
        void SetDataCount(const int& N) {
            throw "Not implemented";
        }
        Property<TRnReleaseMultiDataSpaceModel,
            int,
            &TRnReleaseMultiDataSpaceModel::GetDataCount,
            &TRnReleaseMultiDataSpaceModel::SetDataCount> DataCount;

        TMask ** GetMasks() const {return fMasks;}
        void SetMasks(TMask ** const& Masks) {
            throw "Not implemented";
        }
        Property<TRnReleaseMultiDataSpaceModel,
            TMask **,
            &TRnReleaseMultiDataSpaceModel::GetMasks,
            &TRnReleaseMultiDataSpaceModel::SetMasks> Masks;
      protected:
        const int fDataCount;
        TMask * * fMasks;
      private:
        TRnReleaseMultiDataSpaceModel(const TRnReleaseMultiDataSpaceModel &): fDataCount(0), fMasks(NULL),
            DataCount(this), Masks(this) {}
        TRnReleaseMultiDataSpaceModel &operator=(const TRnReleaseMultiDataSpaceModel &) {return *this;}
};

//#################################################################     


//---------------------MultiData (End)-----------------------------------------






//---------------------1D (Begin)-----------------------------------------------


//#################################################################     

class T1dCorrespondence: public TMultiDataCorrespondence
{
       public:
        T1dCorrespondence(const int lInd, const int lI): TMultiDataCorrespondence(lInd),fI(lI), I(this)
          {
             if (fI<0)
              {
                throw "T1dCorrespondence constructor error: The element ID value is incorrect.";
              }
          }

		int GetI() const {return fI;}
        void SetI(const int& I) {
            throw "Not implemented";
        }
        Property<T1dCorrespondence, int, &T1dCorrespondence::GetI, &T1dCorrespondence::SetI> I;
       protected:
        virtual T1dCorrespondence &Copy() const
         {
           throw "Error in T1dCorrespondence::Copy: The function is not implemented.";
         }
        const int fI;
       private:
        T1dCorrespondence(const T1dCorrespondence &):TMultiDataCorrespondence(0), fI(0), I(this) {}
        T1dCorrespondence &operator=(const T1dCorrespondence &) {return *this;}
};

//#################################################################

class T1dInverseCorrespondence: public TMultiDataInverseCorrespondence
{
       friend class TRnRelease1dSpaceModel;
       public:        
        const int *operator[](const int lInd) const
         {
           if ( (lInd<0) || (lInd>fDataCount-1) )
            {
              throw "Error in T1dInverseCorrespondence::operator[]: The index value is out of range.";
            }
           return fInvCor[lInd];
         }
        T1dInverseCorrespondence(const int lDataCount,const int lN): TMultiDataInverseCorrespondence(lDataCount),fN(lN)
          {

              if (fN<=0)
              {
                throw "T1dInverseCorrespondence constructor error: The elements count value is incorrect.";
              }

              fInvCor = new int * [fDataCount];
              int ind;
              int i;
              for (ind=0; ind<fDataCount; ind++)
                {
                  fInvCor[ind] = new int [fN];
                }

              for (ind=0; ind<fDataCount; ind++)
                {
                    for (i=0; i<fN; i++)
                       {
                         fInvCor[ind][i] = 0;
                       }
                }


          }
       ~T1dInverseCorrespondence()
          {
              int lInd;
              for (lInd=0; lInd<fDataCount; lInd++)
                {
                  delete[] fInvCor[lInd];
                }
			  delete[] fInvCor;// fInvCor = NULL;
          }
       protected:
        virtual T1dInverseCorrespondence &Copy() const
         {
           throw "Error in T1dInverseCorrespondence::Copy: The function is not implemented.";
         }
        const int fN;
       private:
        T1dInverseCorrespondence(const T1dInverseCorrespondence &):TMultiDataInverseCorrespondence(0),fN(0) {}
        T1dInverseCorrespondence &operator=(const T1dInverseCorrespondence &) {return *this;}

        int **fInvCor;
     };

//#################################################################

class TRnRelease1dSpaceModel: public TRnReleaseMultiDataSpaceModel
    {
      public:
       TRnRelease1dSpaceModel(const int lDataCount, const TMask * const * lMasks);
       TRnRelease1dSpaceModel(const int lDataCount,const int lN,const unsigned char *lMask);
       virtual ~TRnRelease1dSpaceModel();
       virtual TRnRelease1dSpaceModel &Copy();

       const T1dNumberMask & GetMask(const int lInd);
      protected:
      private:
       TRnRelease1dSpaceModel(const TRnRelease1dSpaceModel &):TRnReleaseMultiDataSpaceModel(0,NULL) {}
       TRnRelease1dSpaceModel &operator=(const TRnRelease1dSpaceModel &) {return *this;}
       
	   //Help constructor functions (Begin)
       void ReleaseConstructor();
       static const TMask * const * BuildMasks(const int lDataCount, const int lN,const unsigned char *lMask);       
       static TMask * * builtMasks;
       static TMask *  builtMask;
	   //Help constructor functions (End)
    };

//#################################################################


//---------------------1D (End)-----------------------------------------






//---------------------2D (Begin)-----------------------------------------


//---------------------2D (End)-----------------------------------------







//---------------------3D (Begin)-----------------------------------------

//#################################################################

class T3dCorrespondence: public TMultiDataCorrespondence
{
       public:
        T3dCorrespondence(const int lInd, const int lI, const int lJ, const int lK):
            TMultiDataCorrespondence(lInd),fI(lI), fJ(lJ), fK(lK),
            I(this), J(this), K(this)
          {
            if (fI<0)
              {
                throw "T2dCorrespondence constructor error: The first dimension element ID value is incorrect.";
              }
             if (fJ<0)
              {
                throw "T2dCorrespondence constructor error: The second dimension element ID value is incorrect.";
              }
             if (fK<0)
              {
                throw "T2dCorrespondence constructor error: The third dimension element ID value is incorrect.";
              }
          }

        int GetI() const {return fI;}
        void SetI(const int& I) {
            throw "Not implemented";
        }
        Property<T3dCorrespondence, int, &T3dCorrespondence::GetI, &T3dCorrespondence::SetI> I;

        int GetJ() const {return fJ;}
        void SetJ(const int& J) {
            throw "Not implemented";
        }
        Property<T3dCorrespondence, int, &T3dCorrespondence::GetJ, &T3dCorrespondence::SetJ> J;

        int GetK() const {return fK;}
        void SetK(const int& K) {
            throw "Not implemented";
        }
        Property<T3dCorrespondence, int, &T3dCorrespondence::GetK, &T3dCorrespondence::SetK> K;
       protected:
        virtual T3dCorrespondence &Copy()const
         {
           throw "Error in T3dCorrespondence::Copy: The function is not implemented.";
         }
        int fI;  
        int fJ;
        int fK;
       private:
        T3dCorrespondence(const T3dCorrespondence &):TMultiDataCorrespondence(0), fI(0), fJ(0), fK(0),
               I(this), J(this), K(this) {}
        T3dCorrespondence &operator=(const T3dCorrespondence &) {return *this;}
};

//#################################################################

class T3dInverseCorrespondence: public TMultiDataInverseCorrespondence
     {
       /*
       Этот класс должен использовать закрытую информацию,
       поскольку он то эту информацию и формирует.
       */
       friend class TRnRelease3dSpaceModel;

       public:
        /*
          Внимание! Операция чтения данных не проверяет корректность
          переданных индексов.
        */
        const int * const * const * operator[](int lInd) const
          {
            if ( (lInd<0) || (lInd>fDataCount-1) )
            {
              throw "Error in T3dInverseCorrespondence::operator[]: The index value is out of range.";
            }
            return fInvCor[lInd];
          }
        T3dInverseCorrespondence(const int lDataCount, const int lN, const int lL, const int lM): TMultiDataInverseCorrespondence(lDataCount),fN(lN),fL(lL),fM(lM)
          {
              if (fN<=0)
              {
                throw "T3dInverseCorrespondence constructor error: The first dimension elements count value is incorrect.";
              }

              if (fL<=0)
              {
                throw "T3dInverseCorrespondence constructor error: The second dimension elements count value is incorrect.";
              }

              if (fM<=0)
              {
                throw "T3dInverseCorrespondence constructor error: The third dimension elements count value is incorrect.";
              }


              fInvCor = new int *** [fDataCount];
              int ind;
              int i;
              int j;
              int k;
              for (ind=0; ind<fDataCount; ind++)
                {
                  fInvCor[ind] = new int **[fN];
                  for (i=0; i<fN; i++)
                    {
                       fInvCor[ind][i] = new int *[fL];
                       for (j=0; j<fL; j++)
                         {
                            fInvCor[ind][i][j] = new int [fM];
                         }

                    }
                }

              for (ind=0; ind<fDataCount; ind++)
                {
                  for (i=0; i<fN; i++)
                    {
                       for (j=0; j<fL; j++)
                         {
                            for (k=0; k<fM; k++)
                               {
                                  fInvCor[ind][i][j][k] = 0;
                               }
                         }

                    }
                }



          }
       ~T3dInverseCorrespondence()
          {
              int ind;
              int i;
              int j;
              for (ind=0; ind<fDataCount; ind++)
                {
                  for (i=0; i<fN; i++)
                   {
                     for (j=0; j<fL; j++)
                      {
                        delete[] fInvCor[ind][i][j];
                      }
                     delete[] fInvCor[ind][i];
                   }
                  delete[] fInvCor[ind];
                }
              delete[] fInvCor;// fInvCor = NULL;
          }
       protected:
        virtual T3dInverseCorrespondence &Copy() const
         {
           throw "Error in T3dInverseCorrespondence::Copy: The function is not implemented.";
         }
        int fN; int fL; int fM; //количества данных
       private:
        T3dInverseCorrespondence(const T3dInverseCorrespondence &):TMultiDataInverseCorrespondence(0),fN(0),fL(0),fM(0) {}
        T3dInverseCorrespondence &operator=(const T3dInverseCorrespondence &) {return *this;}

        int ****fInvCor;
     };

//#################################################################

class TRnRelease3dSpaceModel: public TRnReleaseMultiDataSpaceModel
    {
      public:
       TRnRelease3dSpaceModel(const int lDataCount, const TMask * const * lMasks);
       /*
         Следующий конструктор означает, что нужно создать маску по переданным
         параметрам и сопоставить ее копии всем блокам мультиданных.
       */
       TRnRelease3dSpaceModel(const int lDataCount, const int lN, const int lL, const int lM, unsigned char ***lMask);
       virtual ~TRnRelease3dSpaceModel();
       virtual TRnRelease3dSpaceModel &Copy();

	   const T3dNumberMask & GetMask(const int lInd);
      protected:
      private:
       TRnRelease3dSpaceModel(const TRnRelease3dSpaceModel &):TRnReleaseMultiDataSpaceModel(0,NULL) {}
       TRnRelease3dSpaceModel &operator=(const TRnRelease3dSpaceModel &) {return *this;}
       
       /*
         Вспомогательная функция для реализации
         основных действий конструкторов.
       */
       void ReleaseConstructor();
       /*
         Вспомогательная статическая функция для корректного
         функционирования второго конструктора.
       */
       static const TMask * const * BuildMasks(const int lDataCount, const int lN, const int lL, const int lM, const unsigned char * const * const * lMask);
       /*
         Вспомогательные статические поля для корректного
         функционирования второго конструктора.
       */
       static TMask * * builtMasks;
       static TMask *  builtMask;
    };

//#################################################################


//---------------------3D (End)-----------------------------------------
 

#endif
