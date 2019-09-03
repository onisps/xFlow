#ifndef OperatorsH
#define OperatorsH

class TLinearOperator;

#include "XSpace.h"
#include "XGrid.h"
#include "property.h"
#include <cstddef>
#include <stdio.h>
#include <math.h>

const int TDefinedBorderPoint = 3;
const int TPreDefinedBorderPoint = 4;
const int TNormalBorderPoint = 5;
const int TPreNormalBorderPoint = 6;
const int TEquationBorderPoint = 7;
const int TPreEquationBorderPoint = 8;

//#################################################################

class TLinearOperator
{
     /*
       Данные классы в перегруженных операциях обязаны иметь дступ
       к важному полю - fArgument, хранящее ссылку на вектор, на котором
       вызывается оператор. Это не очень хорошо.
     */
     friend class TLinearSpace;
     friend class TEuclidSpace;
     public:
       TLinearOperator():fArgument(NULL), Argument(this) {};
       virtual ~TLinearOperator() {}

       virtual int GetClassName() const = 0;

       /*
         Очень выжный метод - операция вызова. Ее действие заключается в
         формировании fArgument и конструировании нового объекта класса
         линейной комбинации.
       */
       TLinearCombinationWith2LinearOperatorOperands &operator()(const TSpace &) const;

       //Не очень хорошо, что возвращается не константная ссылка.
	   TSpace &GetArgument() const;
	   void SetArgument(TSpace & Argument) {
           throw "Not implemented";
       }
       Property<TLinearOperator, TSpace &, &TLinearOperator::GetArgument, &TLinearOperator::SetArgument> Argument;

    protected:
         //Не очень хорошо, что mutable.
         mutable TSpace *fArgument;
         /*
         Метод вызываемый операцией вызова и необходимая для проверки
         корректности вызова.
         */
         virtual void CheckCorrectCall() const = 0;
     private:              
       TLinearOperator(const TLinearOperator &);
       TLinearOperator &operator=(const TLinearOperator &) {return *this;}
};

//#################################################################

class TRnLinearOperator: public TLinearOperator
{
     public:
       TRnLinearOperator(const int lModelName);
       virtual ~TRnLinearOperator();
       	   
	   virtual void operator~() const = 0; //использовать сопряженный
	   
	   
	   virtual void ConstructActivityTables(); //Автозаполнение таблиц активности
	   bool IsActivityTableConstructed() {return (fActivityTable!=NULL);}
	   bool IsActivityTableConjugateConstructed() {return (fActivityTableConjugate!=NULL);}

	   int *GetRelatedIndexes(int li) const;
       
       
	   void ConstructIndexesGroupsTable();//Автозаполнение таблицы групп индексов
       void ConstructMatrix(); //Построение матрицы оператора.

       //Временно-потом сделать нормально (Начало)
       TSpace **helpService;
       int helpServiceCount;

       double *helpService1;
       int helpService1Count;
       //Временно-потом сделать нормально (Конец)

	   
	   int** GetiActivityTable() {return fActivityTable;}
	   int** GetiActivityTableConjugate() {return fActivityTableConjugate;}

	   void SetActivityTable(int** lActivityTable);
	   void SetActivityTableConjugate(int** lActivityTableConjugate);

	   int GetActualIndexesCount() const {return fActualIndexesCount;}
	   void SetActualIndexesCount(const int& N) {
           throw "Not implemented";
       }
       Property<TRnLinearOperator,
           int,
           &TRnLinearOperator::GetActualIndexesCount,
           &TRnLinearOperator::SetActualIndexesCount> ActualIndexesCount;

	   int* GetActualIndexes() const {return fActualIndexes;}
	   void SetActualIndexes(int * const & N) {
           throw "Not implemented";
       }
       Property<TRnLinearOperator,
           int *,
           &TRnLinearOperator::GetActualIndexes,
           &TRnLinearOperator::SetActualIndexes> ActualIndexes;

	   int GetMaxRelatedIndexesCount() const;
	   void SetMaxRelatedIndexesCount(const int& N) {
           throw "Not implemented";
       }
       Property<TRnLinearOperator,
           int,
           &TRnLinearOperator::GetMaxRelatedIndexesCount,
           &TRnLinearOperator::SetMaxRelatedIndexesCount> MaxRelatedIndexesCount;

	   int GetIndexesGroupsCount() const {return fIndexesGroupsCount;}
	   void SetIndexesGroupsCount(const int& N) {
           throw "Not implemented";
       }
       Property<TRnLinearOperator,
           int,
           &TRnLinearOperator::GetIndexesGroupsCount,
           &TRnLinearOperator::SetIndexesGroupsCount> IndexesGroupsCount;

	   int** GetIndexesGroupsTable() const {return fIndexesGroupsTable;}
	   void SetIndexesGroupsTable(int ** const & N) {
           throw "Not implemented";
       }
       Property<TRnLinearOperator,
           int **,
           &TRnLinearOperator::GetIndexesGroupsTable,
           &TRnLinearOperator::SetIndexesGroupsTable> IndexesGroupsTable;

       /*
        * TODO: This field is needed?
		*__declspec(property(get = GetMatrixElement)) double MatrixElement[][];
        */
	   double GetMatrixElement(const int li,const int lj) const; //Получить элемент матрицы.

	   const double * const * const GetMatrix() const; //Получить матрицу.
	   void SetMatrix(const double * const * const & N) {
           throw "Not implemented";
       }
       Property<TRnLinearOperator,
           const double * const * const,
           &TRnLinearOperator::GetMatrix,
           &TRnLinearOperator::SetMatrix> Matrix;

	   const int GetModelName() const {return fModelName;}       
	   void SetModelName(const int& N) {
           throw "Not implemented";
       }
       Property<TRnLinearOperator,
           const int,
           &TRnLinearOperator::GetModelName,
           &TRnLinearOperator::SetModelName> ModelName;

	   virtual TSpace &NewArgument() const = 0;
	   virtual TRnLinearOperator &Copy() const = 0;
    protected:
       int fModelName; //Имя модели решения для корректности мат. операций.
       TRnSpaceModel *fModel; //Модель решения для корректности мат. операций.
       
	   int **fActivityTable; //Таблица активности оператора.
	   bool fExternalActivityTable; //Используется внешняя таблица активности оператора.
	   int ** fActivityTableConjugate; //Таблица активности сопряженного оператора.
       bool fExternalActivityTableConjugate; //Используется внешняя таблица активности сопряженного оператор.
	   void DeleteTable(bool lExternalTable, int ** lTable);

       int **fIndexesGroupsTable; //Табл. непересек. групп индексов.
       int fIndexesGroupsCount; //Количество групп.
       //Нехорошо, тчо mutable
       mutable int fActualIndexesCount; //Количество актуальных индесов.
       mutable int *fActualIndexes; //Актуальные индексы.
       double **fMatrix; //Матрица оператора.

       const TCorrespondence *GetModelCorrespondence(int li) const;
       const TInverseCorrespondence *GetInverseModelCorrespondence() const;       
       virtual void CheckCorrectCall() const;
       virtual void LocalCheckCorrectCall() const = 0;
     private:
       TRnLinearOperator(const TRnLinearOperator &):TLinearOperator(),
            ActualIndexesCount(this), ActualIndexes(this), 
            MaxRelatedIndexesCount(this), IndexesGroupsCount(this),
            IndexesGroupsTable(this), Matrix(this), ModelName(this) {}
       TRnLinearOperator &operator=(const TRnLinearOperator &) {return *this;}       
};

//#################################################################

class TRnRelease1dLinearOperator: public virtual TRnLinearOperator
{
     public:
       TRnRelease1dLinearOperator(const int lModelName);
       virtual ~TRnRelease1dLinearOperator() {}
       virtual double operator()(int idData,int i) const = 0;

	   virtual TSpace &NewArgument() const;
     protected:
         mutable TRnRelease1dSpace *fCastArgument;
         TRnRelease1dSpaceModel *fCastModel;

         virtual void LocalCheckCorrectCall() const;
         virtual void ComfortableCall() const = 0;		 
     private:
       TRnRelease1dLinearOperator(const TRnRelease1dLinearOperator &):TRnLinearOperator(0) {}
       TRnRelease1dLinearOperator &operator=(const TRnRelease1dLinearOperator &) {return *this;}       
   };

//#################################################################

class TRnRelease3dLinearOperator: public virtual TRnLinearOperator
{
     public:
       TRnRelease3dLinearOperator(const int lModelName);
	   virtual ~TRnRelease3dLinearOperator() {};
       virtual double operator()(int idData,int i,int j, int k) const = 0;

	   virtual TSpace &NewArgument() const;
     protected:
         mutable TRnRelease3dSpace *fCastArgument;
         TRnRelease3dSpaceModel *fCastModel;

         virtual void LocalCheckCorrectCall() const;
         virtual void ComfortableCall() const = 0;		 
     private:
       TRnRelease3dLinearOperator(const TRnRelease3dLinearOperator &);
       TRnRelease3dLinearOperator &operator=(const TRnRelease3dLinearOperator &);
};

//#################################################################

class TFiniteDifferenceOperator: public virtual TRnLinearOperator
{
     public:
       TFiniteDifferenceOperator(const int lModelName, const int lGridsCount, const TNodesGrid* const* lGrids);
       virtual ~TFiniteDifferenceOperator();
    protected:
        const int fGridsCount;
		TNodesGrid **fGrids;
     private:
       TFiniteDifferenceOperator(const TFiniteDifferenceOperator &): TRnLinearOperator(0),fGridsCount(0),fGrids(NULL) {}
       TFiniteDifferenceOperator &operator=(const TFiniteDifferenceOperator &) {return *this;}
};

//#################################################################

class T1dLaplaceOperator: public TRnRelease1dLinearOperator,TFiniteDifferenceOperator
{
     public:       
       T1dLaplaceOperator(const int lModelName, const int lGridsCount, const TNodesGrid* const* lGrids);
	   virtual ~T1dLaplaceOperator() {};
       
       virtual int GetClassName() const {return 1;}

       virtual double operator()(int idData,int i) const;

	   virtual void operator~() const {throw "Error in T1dLaplaceOperator::operator~: Operation is not implemented";}

	   virtual TSpace &NewArgument() const {return TRnRelease1dLinearOperator::NewArgument();}
	   virtual TRnLinearOperator &Copy() const {return *(new T1dLaplaceOperator(fModelName,fGridsCount,fGrids));}
     protected:
          const T1dNumberMask *mask;
          const T1dVector* const* fNormals;

		  bool fIsUniformGrid;
          double hx;
          double hxhx;
         
		  virtual void LocalCheckCorrectCall() const {TRnRelease1dLinearOperator::LocalCheckCorrectCall();}

          mutable double *f;
          virtual void ComfortableCall() const;		  
     private:      
       T1dLaplaceOperator(const T1dLaplaceOperator &);
       T1dLaplaceOperator &operator=(const T1dLaplaceOperator &);     
   };

//#################################################################

class T3dLaplaceOperator: public TRnRelease3dLinearOperator,TFiniteDifferenceOperator
{
     public:       
       T3dLaplaceOperator(const int lModelName, const int lGridsCount, const TNodesGrid* const* lGrids);
	   virtual ~T3dLaplaceOperator() {};
       
       virtual int GetClassName() const {return 2;}

       virtual double operator()(int idData,int i,int j,int k) const;

	   virtual void operator~() const {throw "Error in T3dLaplaceOperator::operator~: Operation is not implemented";}

	   virtual TSpace &NewArgument() const {return TRnRelease3dLinearOperator::NewArgument();}
	   virtual TRnLinearOperator &Copy() const {return *(new T3dLaplaceOperator(fModelName,fGridsCount,fGrids));}
	   virtual void ConstructActivityTables(); //Ручное заполнение таблицы активности
     protected:
          const T3dNormalGrid* fGrid;
		  const T3dNumberMask *fMask;
          const T3dVector* const* const* const* fNormals;
		  
		  virtual void LocalCheckCorrectCall() const {TRnRelease3dLinearOperator::LocalCheckCorrectCall();}

          mutable double ***f;
          virtual void ComfortableCall() const;		  
     private:      
       T3dLaplaceOperator(const T3dLaplaceOperator &);
       T3dLaplaceOperator &operator=(const T3dLaplaceOperator &);     
   };

//#################################################################

class T1dPackOperator: public TRnRelease1dLinearOperator
{
     public:       
       T1dPackOperator(const int lModelName, const TRnLinearOperator &lA, double** lVector);
	   virtual ~T1dPackOperator();
       
       virtual int GetClassName() const {return 4;}

       virtual double operator()(int idData,int i) const;

	   virtual void operator~() const {fUseTransposeMatrix = !fUseTransposeMatrix;} 

	   virtual TSpace &NewArgument() const {return TRnRelease1dLinearOperator::NewArgument();}
	   virtual TRnLinearOperator &Copy() const 
	   {
	     //throw "Error in T1dPackOperator::Copy: Method is not ready for cipying.";
         return *(new T1dPackOperator(
		                            fModelName,
		                            fRowsCount, fRowsElementsCount,  fMatrixPack,  fColumns,
					                fRowsCountT,fRowsElementsCountT, fMatrixPackT, fColumnsT
					               ));	   
	   }
     protected:
          		  
		  //Matrix (Begin)
		  bool fExternalMatrixPack;
		  int fRowsCount;
		  int* fRowsElementsCount;
	      double** fMatrixPack;
	      int** fColumns;     
		  //Matrix (End)


		  //Transpose matrix (Begin)
		  bool fExternalMatrixPackT;
		  int fRowsCountT;
		  int* fRowsElementsCountT;
	      double** fMatrixPackT;
	      int** fColumnsT;     
		  //Transpose matrix (End)


		  virtual void LocalCheckCorrectCall() const {TRnRelease1dLinearOperator::LocalCheckCorrectCall();}


		  mutable bool fUseTransposeMatrix;

		  //Comfortable сall variables (Begin)
          mutable double *f;

		  //Actual matrix (Begin)
		  mutable int fRowsCountAct;
		  mutable int* fRowsElementsCountAct;
	      mutable double** fMatrixPackAct;
	      mutable int** fColumnsAct;     
		  //Actual matrix (End)

		  //Comfortable сall variables (End)

          virtual void ComfortableCall() const;		  
     private:
       T1dPackOperator(
		               const int lModelName,
		               int lRowsCount,  int* lRowsElementsCount,  double** lMatrixPack,  int** lColumns,
					   int lRowsCountT, int* lRowsElementsCountT, double** lMatrixPackT, int** lColumnsT
					   );	   
       T1dPackOperator(const T1dPackOperator &);
       T1dPackOperator &operator=(const T1dPackOperator &);     
   };


class T3dVariableDensityOperator: public T3dLaplaceOperator
{
public:
    T3dVariableDensityOperator(
        const int lModelName,
        const int lGridsCount,
        const TNodesGrid* const* lGrids
    ): T3dLaplaceOperator(lModelName, lGridsCount, lGrids), TRnLinearOperator(lModelName) {};
    virtual double operator()(int idData, int i, int j, int k) const;
    void withDensity(TRnRelease3dSpace &density);
private:
    double ***density;
};

//#################################################################

class T3dVariableDivergentedOperator : public T3dLaplaceOperator
{
public:
	T3dVariableDivergentedOperator(
		const int lModelName,
		const int lGridsCount,
		const TNodesGrid* const* lGrids
	) : T3dLaplaceOperator(lModelName, lGridsCount, lGrids), TRnLinearOperator(lModelName) {};
	virtual double operator()(int idData, int i, int j, int k) const;
	void withDensity(TRnRelease3dSpace &density);
	 ~T3dVariableDivergentedOperator();
private:
	double ***density = NULL;
};

//#################################################################

class T3dVariableVsStreamOperator : public T3dLaplaceOperator
{
public:
	T3dVariableVsStreamOperator(
		const int lModelName,
		const int lGridsCount,
		const TNodesGrid* const* lGrids
	) : T3dLaplaceOperator(lModelName, lGridsCount, lGrids), TRnLinearOperator(lModelName) {};
	virtual double operator()(int idData, TRnRelease3dSpace* U, TRnRelease3dSpace* V, TRnRelease3dSpace* W, int i, int j, int k) const;
	void withDensity(TRnRelease3dSpace &density);
private:
	double ***density;
};

//#################################################################
#endif
