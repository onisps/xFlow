#ifndef MasksH
#define MasksH

#include "property.h"

const int TFictivePoint = 0;
const int TActualPoint = 1;


//#################################################################

class TMask
{
      public:
       TMask() {}
       virtual ~TMask() {}
       virtual TMask &Copy() const = 0;
      protected:
      private:
       TMask(const TMask &) {}
       TMask &operator=(const TMask &) {return *this;}
};

//#################################################################

class T1dNumberMask: public TMask
{
      public:
       T1dNumberMask(const int lN);
       virtual ~T1dNumberMask();
       virtual T1dNumberMask &Copy() const;
       int &operator[](const int i) const;

	   int GetN() const {return fN;}
	   void SetN(const int& N) {this->fN = N;}
       Property<T1dNumberMask, int, &T1dNumberMask::GetN, &T1dNumberMask::SetN> N;

       int * GetMask() const {return fMask;}
       void SetMask(int * const& Mask) {
           throw "Not implemented";
       }
       Property<T1dNumberMask, int *, &T1dNumberMask::GetMask, &T1dNumberMask::SetMask> Mask;
      protected:
       int fN;	   
       int *fMask;
      private:
       T1dNumberMask(const T1dNumberMask &): TMask(), fN(0), N(this), Mask(this) {}
       T1dNumberMask &operator=(const T1dNumberMask &) {return *this;}
};

//#################################################################

class T2dNumberMask: public TMask
{
      public:
       T2dNumberMask(const int lN, const int lL);
       virtual ~T2dNumberMask();
       virtual T2dNumberMask &Copy() const;
       int * operator[](const int i) const;

	   int GetN() const {return fN;}
	   void SetN(const int& N) {
           throw "Not implemented";
       }
       Property<T2dNumberMask, int, &T2dNumberMask::GetN, &T2dNumberMask::SetN> N;

	   int GetL() const {return fL;}
	   void SetL(const int& L) {
           throw "Not implemented";
       }
       Property<T2dNumberMask, int, &T2dNumberMask::GetL, &T2dNumberMask::SetL> L;

       int ** GetMask() const {return fMask;}
       void SetMask(int ** const& Mask) {
           throw "Not implemented";
       }
       Property<T2dNumberMask, int **, &T2dNumberMask::GetMask, &T2dNumberMask::SetMask> Mask;
      protected:
       const int fN;
       const int fL;
       int **fMask;
      private:
       T2dNumberMask(const T2dNumberMask &): TMask(), fN(0), fL(0),
           N(this), L(this), Mask(this) {}
       T2dNumberMask &operator=(const T2dNumberMask &) {return *this;}
};

//#################################################################

class T3dNumberMask: public TMask
{
      public:
       T3dNumberMask(const int lN, const int lL, const int lM);
       virtual ~T3dNumberMask();
       virtual T3dNumberMask &Copy() const;
       int * const * operator[](const int i) const;

	   int GetN() const {return fN;}
	   void SetN(const int& N) {
           throw "Not implemented";
       }
       Property<T3dNumberMask, int, &T3dNumberMask::GetN, &T3dNumberMask::SetN> N;

	   int GetL() const {return fL;}
	   void SetL(const int& L) {
           throw "Not implemented";
       }
       Property<T3dNumberMask, int, &T3dNumberMask::GetL, &T3dNumberMask::SetL> L;

	   int GetM() const {return fM;}
	   void SetM(const int& M) {
           throw "Not implemented";
       }
       Property<T3dNumberMask, int, &T3dNumberMask::GetM, &T3dNumberMask::SetM> M;

       int *** GetMask() const {return fMask;}
       void SetMask(int *** const& Mask) {
           throw "Not implemented";
       }
       Property<T3dNumberMask, int ***, &T3dNumberMask::GetMask, &T3dNumberMask::SetMask> Mask;
      protected:
       const int fN;
       const int fL;
       const int fM;
       int ***fMask;
      private:
       T3dNumberMask(const T3dNumberMask &): TMask(),
           fN(0), fL(0), fM(0),
           N(this), L(this), M(this), Mask(this) {}
       T3dNumberMask &operator=(const T3dNumberMask &) {return *this;}
};

//#################################################################

#endif
