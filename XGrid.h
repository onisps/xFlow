#ifndef NodesGridsH
#define NodesGridsH

#include "XSeparator.h"
#include "XMask.h"
#include "property.h"
#include <cstddef>
#include <stdio.h>
#include <math.h>

const int TBorderPoint = 2;

//#################################################################

class TNodesGrid
{
    public:
     TNodesGrid(): CountNodes(this) {};
     virtual ~TNodesGrid() {};
     virtual TNodesGrid &Copy() const = 0;

     virtual int GetCountNodes() const = 0;
     void SetCountNodes(const int & CountNodes) {
        throw "Not implemented";
     };
     Property<TNodesGrid, int, &TNodesGrid::GetCountNodes, &TNodesGrid::SetCountNodes> CountNodes;
    protected:
    private:
      TNodesGrid(const TNodesGrid &): CountNodes(this) {}
      TNodesGrid &operator=(const TNodesGrid &) {return *this;}
};

//#################################################################

class T1dGrid: public virtual TNodesGrid
{
    public:
     T1dGrid(TSeparator * lSeparator1, bool lDeleteSeparators=true);
     virtual ~T1dGrid();

     const TSeparator &GetSeparator1() const;
     void SetSeparator1(const TSeparator& Separator) {
         throw "Not implemented";
     }
     Property<T1dGrid,
         const TSeparator&,
         &T1dGrid::GetSeparator1,
         &T1dGrid::SetSeparator1> Separator1;

	 virtual int GetCountNodes() const;
    protected:
     const bool fDeleteSeparators;
     TSeparator *fSeparator1;       
    private:
     T1dGrid(const T1dGrid &):TNodesGrid(), fDeleteSeparators(false), Separator1(this) {}
     T1dGrid &operator=(const T1dGrid &){return *this;}
};

//#################################################################

class T2dGrid: public T1dGrid
{
    public:
     T2dGrid(TSeparator * lSeparator1, TSeparator * lSeparator2, bool lDeleteSeparators=true);
     virtual ~T2dGrid();

	 const TSeparator &GetSeparator2() const;
     void SetSeparator2(const TSeparator& Separator) {
         throw "Not implemented";
     }
     Property<T2dGrid,
         const TSeparator&,
         &T2dGrid::GetSeparator2,
         &T2dGrid::SetSeparator2> Separator2;

	 virtual int GetCountNodes() const;       
    protected:
       TSeparator *fSeparator2;       
    private:
      T2dGrid(const T2dGrid &):T1dGrid(NULL, false), Separator2(this) {}
      T2dGrid &operator=(const T2dGrid &){return *this;}
};

//#################################################################

class T3dGrid: public T2dGrid
{
    public:
     T3dGrid(TSeparator * lSeparator1, TSeparator * lSeparator2, TSeparator * lSeparator3, bool lDeleteSeparators=true);
     virtual ~T3dGrid();

	 const TSeparator &GetSeparator3() const;
     void SetSeparator3(const TSeparator& Separator) {
         throw "Not implemented";
     }
     Property<T3dGrid,
         const TSeparator&,
         &T3dGrid::GetSeparator3,
         &T3dGrid::SetSeparator3> Separator3;
     virtual int GetCountNodes() const;
    protected:
     TSeparator *fSeparator3;     
    private:
     T3dGrid(const T3dGrid &):T2dGrid(NULL, NULL, false), Separator3(this) {}
     T3dGrid &operator=(const T3dGrid &){return *this;}
};

//#################################################################

class TMaskGrid: public virtual TNodesGrid
{
    public:
      TMaskGrid(const TMask &lMask);
      virtual ~TMaskGrid();

      TMask &CopyMask() const;

      const TMask &GetMask() const { return *fMask;}
      void SetMask(const TMask& Mask) {
          throw "Not implemented";
      }
      Property<TMaskGrid, const TMask &, &TMaskGrid::GetMask, &TMaskGrid::SetMask> Mask;
    protected:
      TMaskGrid(const TMask *lMask): Mask(this) {}
      TMask * fMask;
    private:
      TMaskGrid(const TMaskGrid &):TNodesGrid(), Mask(this) {}
      TMaskGrid &operator=(const TMaskGrid &) {return *this;}     
};

//#################################################################

class T1dMaskGrid: public T1dGrid, public TMaskGrid
{
    public:
     T1dMaskGrid(TSeparator *lSeparator1, const TMask &lMask, bool lDeleteSeparators=true);
	 T1dMaskGrid(TSeparator *lSeparator1, const int *lMask, bool lDeleteSeparators=true);
     virtual ~T1dMaskGrid() {};


     virtual T1dMaskGrid &Copy() const;
     T1dNumberMask &CopyMask() const {return dynamic_cast<T1dNumberMask &>(TMaskGrid::CopyMask());}
	 virtual int GetCountNodes() const {return T1dGrid::GetCountNodes();}
     int *CopyMaskAsArray() const;

	 const T1dNumberMask &GetMask() const {return *( dynamic_cast<T1dNumberMask *>(fMask) );}
     void SetMask(const T1dNumberMask& Mask) {
         throw "Not implemented";
     }
     Property<T1dMaskGrid, const T1dNumberMask &, &T1dMaskGrid::GetMask, &T1dMaskGrid::SetMask> Mask;
    protected:
    private:
      T1dMaskGrid(const T1dMaskGrid &):T1dGrid(NULL, false),TMaskGrid(NULL),TNodesGrid(), Mask(this) {}
      T1dMaskGrid &operator=(const T1dMaskGrid &) {return *this;}
      
      static const TMask & BuildMasks(const int lN, const int *lMask);      
      static TMask *  builtMask;
};

//#################################################################

class T2dMaskGrid: public T2dGrid, public TMaskGrid
{
    public:
     T2dMaskGrid(TSeparator *lSeparator1, TSeparator *lSeparator2, const TMask &lMask, bool lDeleteSeparators=true);
     T2dMaskGrid(TSeparator *lSeparator1, TSeparator *lSeparator2, const int * const * lMask, bool lDeleteSeparators=true);
     virtual ~T2dMaskGrid() {};


     virtual T2dMaskGrid &Copy() const;
     T2dNumberMask &CopyMask() const {return dynamic_cast<T2dNumberMask &>(TMaskGrid::CopyMask());}
	 virtual int GetCountNodes() const {return T2dGrid::GetCountNodes();}
     int * * CopyMaskAsArray() const;

	 const T2dNumberMask &GetMask() const {return *( dynamic_cast<T2dNumberMask *>(fMask) );}
    protected:
    private:
      T2dMaskGrid(const T2dMaskGrid &):T2dGrid(NULL, NULL, false),TMaskGrid(NULL),TNodesGrid() {}
      T2dMaskGrid &operator=(const T2dMaskGrid &) {return *this;}

      static const TMask & BuildMasks(const int lN, const int lL, const int * const * lMask);
      static TMask *  builtMask;
};

//#################################################################

 class T3dMaskGrid: public T3dGrid, public TMaskGrid
{
    public:
     T3dMaskGrid(TSeparator *lSeparator1, TSeparator *lSeparator2, TSeparator *lSeparator3, const TMask &lMask, bool lDeleteSeparators=true);
     T3dMaskGrid(TSeparator *lSeparator1, TSeparator *lSeparator2, TSeparator *lSeparator3, const int * const * const * lMask, bool lDeleteSeparators=true);
     virtual ~T3dMaskGrid() {};


     virtual T3dMaskGrid &Copy() const;
     T3dNumberMask &CopyMask() const {return dynamic_cast<T3dNumberMask &>(TMaskGrid::CopyMask());}
	 virtual int GetCountNodes() const {return T3dGrid::GetCountNodes();}
     int * * * CopyMaskAsArray() const;

	 const T3dNumberMask &GetMask() const {return *( dynamic_cast<T3dNumberMask *>(fMask) );}
     void SetMask(const T3dNumberMask& Mask) {
         throw "Not implemented";
     }
     Property<T3dMaskGrid, const T3dNumberMask &, &T3dMaskGrid::GetMask, &T3dMaskGrid::SetMask> Mask;
    protected:
    private:
      T3dMaskGrid(const T3dMaskGrid &):T3dGrid(NULL,NULL,NULL,false),TMaskGrid(NULL),TNodesGrid(), Mask(this) {}
      T3dMaskGrid &operator=(const T3dMaskGrid &) {return *this;}
     
	  static const TMask & BuildMasks(const int lN, const int lL, const int lM, const int * const * const * lMask);
      static TMask *  builtMask;
};

//#################################################################

struct T1dVector
{
	double X;
};

class T1dNormalGrid: public T1dMaskGrid
{
    public:
     T1dNormalGrid(TSeparator *lSeparator1, const TMask &lMask, const T1dVector* const* lNormals, bool lDeleteSeparators=true);
     virtual ~T1dNormalGrid();
     virtual T1dNormalGrid &Copy() const;

	 const T1dVector* const* GetNormals() const {return fNormals;}
     void SetNormals(const T1dVector * const * const& Normals) {
         throw "Not implemented";
     }
     Property<T1dNormalGrid,
         const T1dVector * const *,
         &T1dNormalGrid::GetNormals,
         &T1dNormalGrid::SetNormals> Normals;
    protected:
     T1dVector** fNormals;
    private:
      T1dNormalGrid(const T1dNormalGrid &):T1dMaskGrid(NULL,NULL,false), Normals(this) {}
      T1dNormalGrid &operator=(const T1dNormalGrid &) {return *this;}
};

//#################################################################

struct T2dVector: public T1dVector
{
	double Y;
};

class T2dNormalGrid: public T2dMaskGrid
{
    public:
     T2dNormalGrid(TSeparator *lSeparator1, TSeparator *lSeparator2, const TMask &lMask, const T2dVector* const* const* lNormals, bool lDeleteSeparators=true);
     virtual ~T2dNormalGrid();
     virtual T2dNormalGrid &Copy() const;

	 const T2dVector* const* const* GetNormals() const {return fNormals;}
     void SetNormals(const T2dVector * const * const * const& Normals) {
         throw "Not implemented";
     }
     Property<T2dNormalGrid,
         const T2dVector * const * const *,
         &T2dNormalGrid::GetNormals,
         &T2dNormalGrid::SetNormals> Normals;
    protected:
     T2dVector*** fNormals;
    private:
      T2dNormalGrid(const T2dNormalGrid &):T2dMaskGrid(NULL,NULL,NULL,false), Normals(this) {}
      T2dNormalGrid &operator=(const T2dNormalGrid &) {return *this;}
};

//#################################################################

struct T3dVector: public T2dVector
{
	double Z;
};

class T3dNormalGrid: public T3dMaskGrid
{
    public:
     T3dNormalGrid(TSeparator *lSeparator1, TSeparator *lSeparator2, TSeparator *lSeparator3, const TMask &lMask, const T3dVector* const* const* const* lNormals, bool lDeleteSeparators=true);
     virtual ~T3dNormalGrid();
     virtual T3dNormalGrid &Copy() const;

	 const T3dVector* const* const* const* GetNormals() const {return fNormals;}
     void SetNormals(const T3dVector * const * const * const * const& Normals) {
         throw "Not implemented";
     }
     Property<T3dNormalGrid,
         const T3dVector * const * const * const *,
         &T3dNormalGrid::GetNormals,
         &T3dNormalGrid::SetNormals> Normals;
    protected:
     T3dVector**** fNormals;
    private:
      T3dNormalGrid(const T3dNormalGrid &):T3dMaskGrid(NULL,NULL,NULL,NULL,false), Normals(this) {}
      T3dNormalGrid &operator=(const T3dNormalGrid &) {return *this;}
};

//#################################################################

#endif
