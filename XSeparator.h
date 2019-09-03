#ifndef SeparatorsH
#define SeparatorsH

#include "property.h"
//#################################################################

class TSeparator
{
    public:
     TSeparator(const double lStartValue, const double lEndValue);
     virtual ~TSeparator() {}
     virtual TSeparator &Copy() const = 0;

	 double GetStartValue() const {return fStartValue;}
     void SetStartValue(const double& N) {
         throw "Not implemented";
     }
     Property<TSeparator,
         double,
         &TSeparator::GetStartValue,
         &TSeparator::SetStartValue> StartValue;

	 double GetEndValue() const {return fEndValue;}
     void SetEndValue(const double& N) {
         throw "Not implemented";
     }
     Property<TSeparator,
         double,
         &TSeparator::GetEndValue,
         &TSeparator::SetEndValue> EndValue;

	 int GetEndIndex() const {return fEndIndex;}
     void SetEndIndex(const int& EndIndex) {
         this->fEndIndex = EndIndex;
     }
     Property<TSeparator,
         int,
         &TSeparator::GetEndIndex,
         &TSeparator::SetEndIndex> EndIndex;

	 int GetCountNodes() const {return fEndIndex+1;}
     void SetCountNodex(const int& CountNodes) {
         throw "Not implemented";
     }
     Property<TSeparator,
         int,
         &TSeparator::GetCountNodes,
         &TSeparator::SetCountNodex> CountNodes;
     
     const double *GetDimension() const {return fDimension;}
     void SetDimension(const double * const & Dimension) {
         throw "Not implemented";
     }
     Property<TSeparator,
         const double *,
         &TSeparator::GetDimension,
         &TSeparator::SetDimension> Dimension;

     const double *GetSeparation() const {return fSeparation;}
     void SetSeparation(const double * const & Separation) {
         throw "Not implemented";
     }
     Property<TSeparator,
         const double *,
         &TSeparator::GetSeparation,
         &TSeparator::SetSeparation> Separation;

    protected:
       const double fStartValue; //Начало отрезка. Совпадает с первым узлом.
       const double fEndValue; //Конец отрезка. Совпадает с последним узлом.
       int fEndIndex; //Номер последнего узла. Номер первого узла - 0.
       double *fDimension; //Узлы разбиения.
       double *fSeparation; //Части разбиения.

    private:
      TSeparator(const TSeparator &):fStartValue(0),fEndValue(0),
           StartValue(this), EndValue(this), EndIndex(this),
           CountNodes(this), Dimension(this), Separation(this) {}
      TSeparator &operator=(const TSeparator &) {return *this;}      
};

//#################################################################

class TUniformSeparator: public TSeparator
{
    public:
     TUniformSeparator(const double lStartValue, const double lEndValue, const double lSeparationValue);
     virtual ~TUniformSeparator();
     virtual TUniformSeparator &Copy() const;

     double GetSeparationValue() const {return fSeparationValue;}
     void SetSeparationValue(const double& SeparationValue) {
         this->fSeparationValue = SeparationValue;
     }
     Property<TUniformSeparator,
         double,
         &TUniformSeparator::GetSeparationValue,
         &TUniformSeparator::SetSeparationValue> SeparationValue;

    protected:
     double fSeparationValue; //Шаг разбиения.
    private:
      TUniformSeparator(const TUniformSeparator &): TSeparator(0,0), fSeparationValue(0), SeparationValue(this) {}
      TUniformSeparator &operator=(const TUniformSeparator &) {return *this;}
};

//#################################################################

class TRandomSeparator: public TSeparator
{
    public:
     TRandomSeparator(int lCountNodes, const double *lDimension);
     virtual ~TRandomSeparator();
     virtual TRandomSeparator &Copy() const;
    private:
      TRandomSeparator(const TRandomSeparator &): TSeparator(0,0) {}
      TRandomSeparator &operator=(const TRandomSeparator &) {return *this;}
};

//#################################################################

#endif
