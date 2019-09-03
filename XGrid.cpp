#pragma hdrstop

#include "XGrid.h"

//#################################################################

//*************************************
     T1dGrid::T1dGrid(TSeparator * lSeparator1, bool lDeleteSeparators):
                       TNodesGrid(), fDeleteSeparators(lDeleteSeparators), Separator1(this)
      {
         fSeparator1 = &(lSeparator1->Copy());
         if (fDeleteSeparators==true)
          {
            delete lSeparator1;
          }
      }
//*************************************
     T1dGrid::~T1dGrid()
      {
         delete fSeparator1;
      }
//*************************************
     int T1dGrid::GetCountNodes() const
      {
        return fSeparator1->CountNodes;
      }
//*************************************
    const TSeparator &T1dGrid::GetSeparator1() const
     {
        return *fSeparator1;
     }
//*************************************

//#################################################################

//*************************************
     T2dGrid::T2dGrid(TSeparator *lSeparator1,TSeparator *lSeparator2,bool lDeleteSeparators):
                      T1dGrid(lSeparator1,lDeleteSeparators), Separator2(this)
     {
        fSeparator2 = &(lSeparator2->Copy());
        if (fDeleteSeparators==true)
          {
            delete lSeparator2;
          }
     }
//*************************************
     T2dGrid::~T2dGrid()
      {
            delete fSeparator2;
      }
//*************************************
     int T2dGrid::GetCountNodes() const
      {
        return (fSeparator1->CountNodes)*(fSeparator2->CountNodes);
      }
//*************************************
    const TSeparator &T2dGrid::GetSeparator2() const
     {
        return *fSeparator2;
     }
//*************************************

//#################################################################

//*************************************
     T3dGrid::T3dGrid(TSeparator * lSeparator1, TSeparator * lSeparator2, TSeparator * lSeparator3, bool lDeleteSeparators):
                      T2dGrid(lSeparator1,lSeparator2,lDeleteSeparators), Separator3(this)
      {
          fSeparator3 = &(lSeparator3->Copy());
          if (fDeleteSeparators==true)
           {
             delete lSeparator3;
           }
      }
//*************************************
     T3dGrid::~T3dGrid()
      {
            delete fSeparator3;
      }
//*************************************
     int T3dGrid::GetCountNodes() const
      {
        return (fSeparator1->CountNodes)*(fSeparator2->CountNodes)*(fSeparator3->CountNodes);
      }
//*************************************
    const TSeparator &T3dGrid::GetSeparator3() const
     {
        return *fSeparator3;
     }
//*************************************

//#################################################################

//*************************************
  TMaskGrid::TMaskGrid(const TMask &lMask): Mask(this)
   {
      fMask = &( lMask.Copy() );
   }
//*************************************
  TMaskGrid::~TMaskGrid()
   {
       delete fMask;
   }
//*************************************
  TMask &TMaskGrid::CopyMask() const
   {
      return fMask->Copy();
   }
//*************************************

//#################################################################

//*************************************
     T1dMaskGrid::T1dMaskGrid(TSeparator *lSeparator1, const TMask &lMask, bool lDeleteSeparators):
                 T1dGrid(lSeparator1,lDeleteSeparators), TMaskGrid(lMask), TNodesGrid(), Mask(this)
       {
          T1dNumberMask * p = dynamic_cast<T1dNumberMask *>(fMask);
          if (p==NULL)
           {
             throw "T1dMaskGrid constructor error: Could not cast the mask to T1dNumberMask * type.";
           }
       }
//*************************************
     TMask * T1dMaskGrid::builtMask = NULL;
//*************************************
     const TMask & T1dMaskGrid::BuildMasks(const int lN, const int * lMask)
      {
           T1dNumberMask & mask = *(new T1dNumberMask(lN));
           for (int i=0;i<lN;i++)
            {
               mask[i] = lMask[i];
            }

           builtMask = &mask;
           
           return mask;
      }
//*************************************
	 T1dMaskGrid::T1dMaskGrid(TSeparator *lSeparator1, const int *lMask,bool lDeleteSeparators):
							 T1dGrid(lSeparator1,lDeleteSeparators), TMaskGrid(BuildMasks(lSeparator1->EndIndex+1,lMask)),
                             TNodesGrid(), Mask(this)
	  {
		 delete builtMask;
		 //builtMask = NULL;

		 T1dNumberMask * p = dynamic_cast<T1dNumberMask *>(fMask);
		 if (p==NULL)
		   {
			 throw "T1dMaskGrid constructor error: Could not cast the mask to T1dNumberMask * type.";
		   }
	  }
//*************************************
      int *T1dMaskGrid::CopyMaskAsArray() const
       {
          int fN = fSeparator1->EndIndex;
          int * lMask = new int [fN+1];

          for (int i=0;i<fN+1;i++)
            {
               lMask[i]=(dynamic_cast<T1dNumberMask &>(*fMask))[i];
            }

          return lMask;
       }

//*************************************
      T1dMaskGrid &T1dMaskGrid::Copy() const
       {
          return *(new T1dMaskGrid(fSeparator1, *fMask, false));
       }
//*************************************

//#################################################################

//*************************************
     T2dMaskGrid::T2dMaskGrid(TSeparator *lSeparator1, TSeparator *lSeparator2, const TMask &lMask, bool lDeleteSeparators):
                 T2dGrid(lSeparator1,lSeparator2,lDeleteSeparators), TMaskGrid(lMask), TNodesGrid()
       {
          T2dNumberMask * p = dynamic_cast<T2dNumberMask *>(fMask);
          if (p==NULL)
           {
             throw "T2dMaskGrid constructor error: Could not cast the mask to T2dNumberMask * type.";
           }
       }
//*************************************
     TMask * T2dMaskGrid::builtMask = NULL;
//*************************************
	 const TMask & T2dMaskGrid::BuildMasks(const int lN, const int lL, const int * const * lMask)
	  {
		   T2dNumberMask & mask = *(new T2dNumberMask(lN,lL));
		   for (int i=0;i<lN;i++)
			{
			  for (int j=0;j<lL;j++)
			   {
				  mask[i][j] = lMask[i][j];
			   }
			}

		   builtMask = &mask;

		   return mask;
	  }
//*************************************
	 T2dMaskGrid::T2dMaskGrid(TSeparator *lSeparator1, TSeparator *lSeparator2, const int * const * lMask,bool lDeleteSeparators):
							 T2dGrid(lSeparator1,lSeparator2,lDeleteSeparators),
                             TMaskGrid(BuildMasks(lSeparator1->EndIndex+1,lSeparator2->EndIndex+1,lMask)),
                             TNodesGrid()
	  {
		 delete builtMask;
		 //builtMask = NULL;

		 T2dNumberMask * p = dynamic_cast<T2dNumberMask *>(fMask);
		 if (p==NULL)
		   {
			 throw "T2dMaskGrid constructor error: Could not cast the mask to T2dNumberMask * type.";
		   }
	  }
//*************************************
      int * * T2dMaskGrid::CopyMaskAsArray() const
       {
          int fN = fSeparator1->EndIndex;
          int fL = fSeparator1->EndIndex;
          int * * lMask = new int *[fN+1];

          for (int i=0;i<fN+1;i++)
            {
               lMask[i] = new int [fL+1];
               for (int j=0;j<fL+1;j++)
                {
                  lMask[i][j]=(dynamic_cast<T2dNumberMask &>(*fMask))[i][j];
                }  
            }

          return lMask;
       }

//*************************************
      T2dMaskGrid &T2dMaskGrid::Copy() const
       {
          return *(new T2dMaskGrid(fSeparator1, fSeparator2, *fMask, false));
       }
//*************************************

//#################################################################

//*************************************
     T3dMaskGrid::T3dMaskGrid(TSeparator *lSeparator1, TSeparator *lSeparator2, TSeparator *lSeparator3, const TMask &lMask, bool lDeleteSeparators):
                 T3dGrid(lSeparator1,lSeparator2,lSeparator3,lDeleteSeparators), TMaskGrid(lMask), TNodesGrid(), Mask(this)
       {
          T3dNumberMask * p = dynamic_cast<T3dNumberMask *>(fMask);
          if (p==NULL)
           {
             throw "T3dMaskGrid constructor error: Could not cast the mask to T3dNumberMask * type.";
           }
       }
//*************************************
     TMask * T3dMaskGrid::builtMask = NULL;
//*************************************
     const TMask & T3dMaskGrid::BuildMasks(const int lN, const int lL, const int lM, const int * const * const * lMask)
      {
           T3dNumberMask & mask = *(new T3dNumberMask(lN,lL,lM));
           for (int i=0;i<lN;i++)
            {
              for (int j=0;j<lL;j++)
               {
                 for (int k=0;k<lM;k++)
                  {
                    mask[i][j][k] = lMask[i][j][k];
                  }  
               }
            }

           builtMask = &mask;
           
           return mask;
      }
//*************************************
	 T3dMaskGrid::T3dMaskGrid(TSeparator *lSeparator1, TSeparator *lSeparator2, TSeparator *lSeparator3, const int * const * const * lMask,bool lDeleteSeparators):
							 T3dGrid(lSeparator1,lSeparator2,lSeparator3,lDeleteSeparators),
                             TMaskGrid(BuildMasks(lSeparator1->EndIndex+1,lSeparator2->EndIndex+1,lSeparator3->EndIndex+1,lMask)),
                             TNodesGrid(), Mask(this)
	  {
		 delete builtMask;
		 //builtMask = NULL;

		 T3dNumberMask * p = dynamic_cast<T3dNumberMask *>(fMask);
		 if (p==NULL)
		   {
			 throw "T3dMaskGrid constructor error: Could not cast the mask to T3dNumberMask * type.";
		   }
	  }
//*************************************
	  int * * * T3dMaskGrid::CopyMaskAsArray() const
	   {
		  int fN = fSeparator1->EndIndex;
		  int fL = fSeparator2->EndIndex;
		  int fM = fSeparator3->EndIndex;
		  int * * * lMask = new int **[fN+1];

		  for (int i=0;i<fN+1;i++)
			{
			   lMask[i] = new int *[fL+1];
               for (int j=0;j<fL+1;j++)
                {
                  lMask[i][j] = new int [fM+1];
                  for (int k=0;k<fM+1;k++)
                   {
                     lMask[i][j][k]=(dynamic_cast<T3dNumberMask &>(*fMask))[i][j][k];
                   }
                }  
            }
          return lMask;
       }

//*************************************
      T3dMaskGrid &T3dMaskGrid::Copy() const
       {
          return *(new T3dMaskGrid(fSeparator1, fSeparator2, fSeparator3, *fMask, false));
       }
//*************************************

//#################################################################

//*************************************
     T1dNormalGrid::T1dNormalGrid(TSeparator *lSeparator1,const TMask &lMask, const T1dVector* const* lNormals, bool lDeleteSeparators):
                    T1dMaskGrid(lSeparator1,lMask,lDeleteSeparators), Normals(this)
       {
		   int lN = fSeparator1->EndIndex;
		   fNormals = new T1dVector*[lN+1];
		   for (int i=0;i<lN+1;i++)
		   {
			   fNormals[i] = NULL;
			   if ( (lNormals[i])!=NULL )
			   {
			    fNormals[i] = new T1dVector;
			    *(fNormals[i]) = *(lNormals[i]);
			   }
		   }
       }
//*************************************
      T1dNormalGrid::~T1dNormalGrid()
	   {
           int lN = fSeparator1->EndIndex;
		   for (int i=0;i<lN+1;i++)
		   {			   
			   if ( (fNormals[i])!=NULL )
			   {			   
			     delete fNormals[i];			   
			   }
		   }
		   delete[] fNormals;// fNormals = NULL;
	   }
//*************************************
      T1dNormalGrid &T1dNormalGrid::Copy() const
       {
          return *(new T1dNormalGrid(fSeparator1, *fMask, fNormals, false));
       }
//*************************************

//#################################################################

//*************************************
     T2dNormalGrid::T2dNormalGrid(TSeparator *lSeparator1,
             TSeparator *lSeparator2,
             const TMask &lMask,
             const T2dVector* const* const* lNormals,
             bool lDeleteSeparators): T2dMaskGrid(lSeparator1,lSeparator2,lMask,lDeleteSeparators), Normals(this)
       {
		   int lN = fSeparator1->EndIndex;
		   int lL = fSeparator2->EndIndex;
		   fNormals = new T2dVector**[lN+1];
		   for (int i=0;i<lN+1;i++)
		   {
			   fNormals[i] = new T2dVector*[lL+1];
			   for (int j=0;j<lL+1;j++)
		       {
                 fNormals[i][j] = NULL;
			     if ( (lNormals[i][j])!=NULL )
			     {			   
			       fNormals[i][j] = new T2dVector;   
			       *(fNormals[i][j]) = *(lNormals[i][j]);
				 }
			   }
		   }
       }
//*************************************
      T2dNormalGrid::~T2dNormalGrid()
	   {
           int lN = fSeparator1->EndIndex;
           int lL = fSeparator2->EndIndex;
		   for (int i=0;i<lN+1;i++)
		   {
			   for (int j=0;j<lL+1;j++)
		       {				  
			      if ( (fNormals[i][j])!=NULL )
			      {			   
				    delete fNormals[i][j];			   
				  }  
			   }
			   delete[] fNormals[i];			   
		   }
		   delete[] fNormals;// fNormals = NULL;
	   }
//*************************************
      T2dNormalGrid &T2dNormalGrid::Copy() const
       {
          return *(new T2dNormalGrid(fSeparator1, fSeparator2, *fMask, fNormals, false));
       }
//*************************************

//#################################################################

//*************************************
     T3dNormalGrid::T3dNormalGrid(
             TSeparator *lSeparator1,
             TSeparator *lSeparator2,
             TSeparator *lSeparator3,
             const TMask &lMask,
             const T3dVector* const* const* const* lNormals,
             bool lDeleteSeparators
             ): T3dMaskGrid(lSeparator1,lSeparator2,lSeparator3,lMask,lDeleteSeparators), Normals(this)
       {
		   int lN = fSeparator1->EndIndex;
		   int lL = fSeparator2->EndIndex;
		   int lM = fSeparator3->EndIndex;
		   fNormals = new T3dVector***[lN+1];
		   for (int i=0;i<lN+1;i++)
		   {
			   fNormals[i] = new T3dVector**[lL+1];
			   for (int j=0;j<lL+1;j++)
		       {
			     fNormals[i][j] = new T3dVector*[lM+1];
				 for (int k=0;k<lM+1;k++)
		         {
					 fNormals[i][j][k] = NULL;
			         if ( (lNormals[i][j][k])!=NULL )
			         {
					  fNormals[i][j][k] = new T3dVector;
					  *(fNormals[i][j][k]) = *(lNormals[i][j][k]);
					 }
				 } 			     
			   }
		   }
       }
//*************************************
      T3dNormalGrid::~T3dNormalGrid()
	   {
           int lN = fSeparator1->EndIndex;
           int lL = fSeparator2->EndIndex;
		   int lM = fSeparator3->EndIndex;
		   for (int i=0;i<lN+1;i++)
		   {
			   for (int j=0;j<lL+1;j++)
		       {
                  for (int k=0;k<lM+1;k++)
				  {
                     if ( (fNormals[i][j][k])!=NULL )   
					 {  
                        delete fNormals[i][j][k];			   
					 }
				  }
				  delete[] fNormals[i][j];			   
			   }
			   delete[] fNormals[i];			   
		   }
		   delete[] fNormals;// fNormals = NULL;
	   }
//*************************************
      T3dNormalGrid &T3dNormalGrid::Copy() const
       {
          return *(new T3dNormalGrid(fSeparator1, fSeparator2, fSeparator3, *fMask, fNormals, false));
       }
//*************************************

//#################################################################
