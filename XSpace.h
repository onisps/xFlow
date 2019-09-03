#ifndef SpacesH
#define SpacesH

#include "XSpaceModel.h"
#include "property.h"

//---------------------Abstract (Begin)------------------------------


//#################################################################

class TSpace
{
    public:
     TSpace();
     virtual ~TSpace();
     TSpace &operator|=(const TSpace &s);
     TSpace &CreateInstance();
     TSpace &Copy();
	 TSpace &Assign(const TSpace &s);
     void SetAsPermanent();
     void SetAsIntermediate();
     bool IsIntermediate();

     static void DeleteIntermediate(const TSpace &);

     //Временное - вспомогательное (Начало)
       static int CountInstances;  //текущее количество экземпляров класса
       static bool ShowRemainCountInstances; //показывать ли в дестр-ре ост. кол. экз. класса
     //Временное - вспомогательное (Конец)

    protected:
     static bool deleteIntermediateItems;
     bool intermediate;
	 virtual TSpace &NewItem() = 0;
     virtual void ReleaseOperatorEqual(const TSpace &,bool isCopying)=0;
    private:
      TSpace(const TSpace &);
      TSpace &operator=(const TSpace &s);
};

//#################################################################

class TMetricSpace : public virtual TSpace
{
     public:
      TMetricSpace();
      virtual ~TMetricSpace();
      double Metric(const TSpace &s, int type);
     protected:
      virtual double ReleaseMetric(const TSpace &s,int type)=0;
};

//#################################################################

class TFullSpace : public TMetricSpace
{
     public:
      TFullSpace();
      virtual ~TFullSpace();
};

//#################################################################


//---------------------Abstract (End)------------------------------






//---------------------Linear (Begin)------------------------------

class TLinearCombinationWith2Operands;
class TLinearCombinationWith2LinearOperatorOperands;

class TRnRelease1dSpace;
class TRnRelease2dSpace;
class TRnRelease3dSpace;
  
#include "XOperator.h"

//#################################################################

class TLinearSpace : public virtual TSpace
{
     public:
       TLinearSpace();
       virtual ~TLinearSpace();

       TSpace &operator|=(const TLinearCombinationWith2Operands &);
       TSpace &operator|=(const TLinearCombinationWith2LinearOperatorOperands &);
       TSpace &operator|=(const int);

       TLinearCombinationWith2Operands &operator+(const TSpace &);
       TLinearCombinationWith2Operands &operator-(const TSpace &);

       TLinearCombinationWith2Operands &operator+(TLinearCombinationWith2Operands &);
       TLinearCombinationWith2Operands &operator-(TLinearCombinationWith2Operands &);

       TLinearCombinationWith2Operands &operator*(double);


     protected:
       virtual void ToZero()=0;
       virtual void ReleaseAssignSpeedLinearCombination(double p1,const TSpace *s1,
                                      double p2,const TSpace *s2,bool operation=true)=0;
       virtual void ReleaseAssignSpeedLinearCombination(const TLinearOperator *A1,
                                      const TSpace *s2,bool operation=true)=0;
};

//#################################################################   

class TLinearCombinationWith2Operands
{
             friend class TLinearSpace;

             public:
                TLinearCombinationWith2Operands(double lp1, const TSpace &ls1);
                TLinearCombinationWith2Operands(const TSpace &ls1, const TSpace &ls2,
                                                bool loperation=true);
                virtual ~TLinearCombinationWith2Operands();

                const TSpace *GetFirstOperand() const;
                const TSpace *GetSecondOperand() const;

                void AddSecondOperand(double lp, const TSpace &ls,bool loperation=true);
                void AddSecondOperand(const TSpace &ls,bool loperation=true);

                void InsertSecondOperand(double lp, const TSpace &ls,bool loperation=true);
                void InsertSecondOperand(const TSpace &ls,bool loperation=true);



                TLinearCombinationWith2Operands &operator+(const TLinearCombinationWith2Operands &);
                TLinearCombinationWith2Operands &operator-(const TLinearCombinationWith2Operands &);

                TLinearCombinationWith2Operands &operator+(const TSpace &);
                TLinearCombinationWith2Operands &operator-(const TSpace &);

             private:
                double p1,p2;
                const TSpace *s1,*s2;
                bool operation;
};

//#################################################################

TLinearCombinationWith2Operands &operator*(double, TLinearSpace &);

class TEuclidSpace;

//#################################################################

class TLinearCombinationWith2LinearOperatorOperands
{
             friend class TLinearSpace;
             friend class TEuclidSpace;

             public:
                TLinearCombinationWith2LinearOperatorOperands(const TLinearOperator &lA1);
                virtual ~TLinearCombinationWith2LinearOperatorOperands();

                const TSpace *GetSecondOperand() const;

                void AddSecondOperand(const TSpace &ls,bool loperation=true);


                TLinearCombinationWith2LinearOperatorOperands &operator+(const TSpace &);
                TLinearCombinationWith2LinearOperatorOperands &operator-(const TSpace &);

             private:
                const TLinearOperator *A1;
                const TSpace *s2;
                bool operation;
};

//#################################################################


//---------------------Linear (End)------------------------------





//---------------------Abstract (Begin)------------------------------


//#################################################################

class TNormalizedSpace : public TLinearSpace
{
     public:
       TNormalizedSpace();
       virtual ~TNormalizedSpace();
       double Norma(int type);
       double operator!();
     protected:
       virtual double ReleaseNorma(int type)=0;
};

//#################################################################

class TEuclidSpace : public TNormalizedSpace
{
     public:
      TEuclidSpace();
      virtual ~TEuclidSpace();
      double Innerproduct(const TSpace &s, int type);
      double Innerproduct(const TLinearCombinationWith2LinearOperatorOperands &lc, int type);
      double operator,(const TSpace &s);
      double operator,(const TLinearCombinationWith2LinearOperatorOperands &lc);
      double NormaAccordingInnerproduct(int type);
     protected:
      virtual double ReleaseInnerproduct(const TSpace &s, int type)=0;
      virtual double ReleaseInnerproduct(const TLinearOperator *lA1,const TSpace *s2, int type)=0;
};

double operator,(const TLinearCombinationWith2LinearOperatorOperands &lc,const TEuclidSpace &s);

//#################################################################

class TBanahSpace : public TEuclidSpace, public TFullSpace
{
     public:
      TBanahSpace();
      virtual ~TBanahSpace();
      double MetricAccordingNorma(const TSpace &s, int type);
};

//#################################################################


//---------------------Abstract (End)------------------------------








//---------------------Rn (Begin)------------------------------


//#################################################################

class TRnSpaceModelAndItsName
{
     friend class TRnSpace;
     public:
      TRnSpaceModelAndItsName(int lName,TRnSpaceModel *lModel);
      virtual ~TRnSpaceModelAndItsName();
     private:
      int fName;
      TRnSpaceModel *fModel;
      int ReferencesCount;

      TRnSpaceModelAndItsName(TRnSpaceModelAndItsName &);
      void operator=(TRnSpaceModelAndItsName &);
};

//#################################################################

class TRnSpace : public TBanahSpace
{
     friend class TRnLinearOperator;

     public:
       static void AddModel(int lName,TRnSpaceModel *lModel,bool deleteModelAfterAdding);
       static void RemoveModel(int lName);
       static void RemoveAllModels();

       TRnSpace(const int lModelName);
       virtual ~TRnSpace();

       //Сервис (Начало)
       double GetAbsMaxElement();
       double GetAbsMaxElement(int &elementIndex);
       double GetAbsMinElement();
       double GetAbsMinElement(int &elementIndex);
       //Сервис(Конец)


       //Индексация (Начало)
	   int GetEndIndex() const;
	   void SetEndIndex(const int& EndIndex) {
           throw "Not implemented";
       }
       Property<TRnSpace, int, &TRnSpace::GetEndIndex, &TRnSpace::SetEndIndex> EndIndex;

       double &operator[](int i) const;

	   int GetActualIndexesCount() const {return fActualIndexesCount;}
	   void SetActualIndexesCount(const int& lCount);
       Property<TRnSpace, int, &TRnSpace::GetActualIndexesCount, &TRnSpace::SetActualIndexesCount> ActualIndexesCount;

	   int * GetActualIndexes() const {return fActualIndexes;}
	   void SetActualIndexes(int * const& ActualIndexes) {
           throw "Not implemented";
       }
       Property<TRnSpace, int *, &TRnSpace::GetActualIndexes, &TRnSpace::SetActualIndexes> ActualIndexes;
    
	   bool GetAutomaticActualIndexesAssigning() const {return fAutomaticActualIndexesAssigning;}
	   void SetAutomaticActualIndexesAssigning(bool const& lAutomaticActualIndexesAssigning) {fAutomaticActualIndexesAssigning = lAutomaticActualIndexesAssigning;}
       Property<TRnSpace,
           bool,
           &TRnSpace::GetAutomaticActualIndexesAssigning,
           &TRnSpace::SetAutomaticActualIndexesAssigning> AutomaticActualIndexesAssigning;
	   //Индексация (Конец)



       //Переопределяемые от TSpace (Начало)
       TRnSpace &CreateInstance();
       TRnSpace &Copy();
       TRnSpace &operator|=(const TRnSpace &s);
	   TRnSpace &Assign(const TRnSpace &s);
       //Переопределяемые от TSpace (Конец)


       //Переопределяемые от TMetricSpace (Начало)
         double Metric(const TRnSpace &s, int type);
       //Переопределяемые от TMetricSpace (Конец)


       //Переопределяемые от TLinearSpace (Начало)
         TRnSpace &operator|=(const TLinearCombinationWith2Operands &);
         TRnSpace &operator|=(const TLinearCombinationWith2LinearOperatorOperands &);
         TRnSpace &operator|=(const int);

         TLinearCombinationWith2Operands &operator+(const TRnSpace &);
         TLinearCombinationWith2Operands &operator-(const TRnSpace &);

         TLinearCombinationWith2Operands &operator+(TLinearCombinationWith2Operands &);
         TLinearCombinationWith2Operands &operator-(TLinearCombinationWith2Operands &);

         TLinearCombinationWith2Operands &operator*(double);
       //Переопределяемые от TLinearSpace (Конец)


       //Переопределяемые от TEuclidSpace (Начало)
         double Innerproduct(const TRnSpace &s, int type);
       //Переопределяемые от TEuclidSpace (Конец)



     protected:
       static TRnSpaceModelAndItsName **fModels;
       static unsigned char fModelsCount;
       static TRnSpaceModel *GetModel(int lName);
       static void IncreaseReferencesCountToModel(int lName);
       static void DecreaseReferencesCountToModel(int lName);	   

       //Сервис (Начало)
       virtual double GetSpecificElement(int type,int &elementIndex) const = 0;
       //Сервис (Конец)

       //Индексация (Начало)
       virtual double &GetByRnIndex(int i) const =0;
       int fActualIndexesCount;
       int *fActualIndexes;
       bool fAutomaticActualIndexesAssigning;
       //Индексация (Конец)

       int fModelName;
       TRnSpaceModel *fModel;

     private:
};

TLinearCombinationWith2Operands &operator*(double,TRnSpace &);

//#################################################################

class TRnData
{
       public:
         TRnData();
         virtual ~TRnData();
       protected:
         virtual void HelpVirtualFunction() const = 0;
};

//#################################################################

class TRnReleaseMultiDataSpace: public TRnSpace
{
      public:
        TRnReleaseMultiDataSpace(const int lModelName);
        virtual ~TRnReleaseMultiDataSpace();

		int GetActiveDataIndex() const {return fActiveDataIndex;}
		void SetActiveDataIndex(const int&);
        Property<TRnReleaseMultiDataSpace,
            int,
            &TRnReleaseMultiDataSpace::GetActiveDataIndex,
            &TRnReleaseMultiDataSpace::SetActiveDataIndex> ActiveDataIndex;

		TRnData** GetMultiData() const {return fMultiData;}
		void SetMultiData(TRnData ** const &);
        Property<TRnReleaseMultiDataSpace,
            TRnData **,
            &TRnReleaseMultiDataSpace::GetMultiData,
            &TRnReleaseMultiDataSpace::SetMultiData> MultiData;

      protected:

       //Переопределяемые от TSpace (Начало)
       virtual void ReleaseOperatorEqual(const TSpace &s,bool isCopying);
       //Переопределяемые от TSpace (Конец)
       virtual void LocalReleaseOperatorEqual(const int id,const TRnData &s,bool isCopying)=0;

       //Переопределяемые от TMetricSpace (Начало)
       virtual double ReleaseMetric(const TSpace &s,int type);
       //Переопределяемые от TMetricSpace (Конец)
       virtual double LocalReleaseMetric(const int id,const TRnData &s,int type)=0;

       //Переопределяемые от TLinearSpace (Начало)
       virtual void ToZero();
       virtual void LocalToZero(const int id)=0;

       virtual void ReleaseAssignSpeedLinearCombination(double p1,const TSpace *s1,
                                      double p2,const TSpace *s2,bool operation=true);

       virtual void LocalReleaseAssignSpeedLinearCombination(const int id,double p1,const TRnData *s1,
                                      double p2,const TRnData *s2,bool operation=true)=0;


       virtual void ReleaseAssignSpeedLinearCombination(const TLinearOperator *A1,
                                      const TSpace *s2,bool operation=true);

       virtual void LocalReleaseAssignSpeedLinearCombination(const int id,const TLinearOperator *A1,
                                      const TRnData *s2,bool operation=true)=0;

       //Переопределяемые от TLinearSpace (Конец)

       //Переопределяемые от TNormalizedSpace (Начало)
       virtual double ReleaseNorma(int type);
       virtual double LocalReleaseNorma(const int id,int type)=0;
       //Переопределяемые от TNormalizedSpace (Конец)

       //Переопределяемые от TEuclidSpace (Начало)
       virtual double ReleaseInnerproduct(const TSpace &s, int type);
       virtual double LocalReleaseInnerproduct(const int id,const TRnData &s, int type)=0;

       virtual double ReleaseInnerproduct(const TLinearOperator *lA1,const TSpace *s2, int type);
       virtual double LocalReleaseInnerproduct(const int id,const TLinearOperator *lA1,const TRnData *s2, int type)=0;
       //Переопределяемые от TEuclidSpace (Конец)


        //Перекрываем Rn...
        virtual double GetSpecificElement(int type,int &elementIndex) const;
        virtual double LocalGetSpecificElement(const int id,int type,int &elementIndex) const = 0;
        int fActiveDataIndex;
        TRnData **fMultiData;
      private:
        TRnReleaseMultiDataSpaceModel *fCastModel;
};

//#################################################################


//---------------------Rn (End)------------------------------





//---------------------Rn1 (Begin)------------------------------



//#################################################################

class TRnRelease1dData: public TRnData
{
      friend class TRnRelease1dSpace;

      public:
       TRnRelease1dData(const int lN);
       virtual ~TRnRelease1dData();

	   double *GetData() const {return fData;}
       void SetData(double * const & Data) {
           throw "Not implemented";
       }
       Property<TRnRelease1dData,
           double *,
           &TRnRelease1dData::GetData,
           &TRnRelease1dData::SetData> Data;

      protected:
       virtual void HelpVirtualFunction() const;
      private:
       int fN;
       double *fData;
};

//#################################################################

class TRnRelease1dSpace : public TRnReleaseMultiDataSpace
{

     public:
       TRnRelease1dSpace(const int lModelName);
       virtual ~TRnRelease1dSpace();

       double &operator[](int);



       
     protected:

       virtual TRnSpace &NewItem();

       virtual void LocalReleaseOperatorEqual(const int id,const TRnData &s,bool isCopying);
       virtual double LocalReleaseMetric(const int id,const TRnData &s,int type);
       virtual void LocalToZero(const int id);
       virtual void LocalReleaseAssignSpeedLinearCombination(const int id,double p1,const TRnData *s1,
                                      double p2,const TRnData *s2,bool operation=true);
       virtual void LocalReleaseAssignSpeedLinearCombination(const int id,const TLinearOperator *A1,
                                      const TRnData *s2,bool operation=true);

       virtual double LocalReleaseNorma(const int id,int type);
       virtual double LocalReleaseInnerproduct(const int id,const TRnData &s, int type);
       virtual double LocalReleaseInnerproduct(const int id,const TLinearOperator *lA1,const TRnData *s2, int type);


       //Перекрываем из Rn...

       virtual double &GetByRnIndex(int i) const;

       //Сервис (Начало)
       virtual double LocalGetSpecificElement(const int id,int type,int &elementIndex) const;
       //Сервис (Конец)


     private:
       TRnRelease1dSpaceModel *fCastModel;
};

//#################################################################



//---------------------Rn1 (End)------------------------------





//---------------------Rn2 (Begin)------------------------------

//---------------------Rn2 (End)--------------------------------



//---------------------Rn3 (Begin)------------------------------


//#################################################################

class TRnRelease3dData: public TRnData
{
      friend class TRnRelease3dSpace;

      public:
       TRnRelease3dData(const int lN,const int lL,const int lM);
       virtual ~TRnRelease3dData();

	   double*** GetData() const {return fData;}
       void SetData(double *** const & Data) {
           throw "Not implemented";
       }
       Property<TRnRelease3dData,
           double ***,
           &TRnRelease3dData::GetData,
           &TRnRelease3dData::SetData> Data;

      protected:
       virtual void HelpVirtualFunction() const;
      private:
       int fN;
       int fL;
       int fM;
       double ***fData;
};

//#################################################################

class TRnRelease3dSpace : public TRnReleaseMultiDataSpace
{

     public:
       TRnRelease3dSpace(const int lModelName);
       virtual ~TRnRelease3dSpace();

       double **operator[](int);

       
       
       

     protected:

	   virtual TRnSpace &NewItem();	  

       virtual void LocalReleaseOperatorEqual(const int id,const TRnData &s,bool isCopying);
       virtual double LocalReleaseMetric(const int id,const TRnData &s,int type);
       virtual void LocalToZero(const int id);
       virtual void LocalReleaseAssignSpeedLinearCombination(const int id,double p1,const TRnData *s1,
                                      double p2,const TRnData *s2,bool operation=true);
       virtual void LocalReleaseAssignSpeedLinearCombination(const int id,const TLinearOperator *A1,
                                      const TRnData *s2,bool operation=true);

       virtual double LocalReleaseNorma(const int id,int type);
       virtual double LocalReleaseInnerproduct(const int id,const TRnData &s, int type);
       virtual double LocalReleaseInnerproduct(const int id,const TLinearOperator *lA1,const TRnData *s2, int type);


       //Перекрываем из Rn...

       virtual double &GetByRnIndex(int i) const; //Перекрываем из Rn

       virtual double LocalGetSpecificElement(const int id,int type,int &elementIndex) const;

     private:
       TRnRelease3dSpaceModel *fCastModel;
};

//#################################################################


//---------------------Rn3 (End)--------------------------------
#endif
