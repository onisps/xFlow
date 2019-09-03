#pragma hdrstop
#include <cstddef>
#include "XSpaceModel.h"

//---------------------MultiData (Begin)-----------------------------------------


//#################################################################

//*************************************
  TRnReleaseMultiDataSpaceModel::TRnReleaseMultiDataSpaceModel(const int lDataCount):
                                 TRnSpaceModel(), fDataCount(lDataCount), DataCount(this), Masks(this)
   {
     if (fDataCount<=0)
      {
        throw "TRnReleaseMultiDataSpaceModel constructor error: The data count value is incorrect.";
      }

      fMasks = NULL;
   }
//*************************************
  TRnReleaseMultiDataSpaceModel::TRnReleaseMultiDataSpaceModel(const int lDataCount, const TMask * const * lMasks):
                                 TRnSpaceModel(), fDataCount(lDataCount), DataCount(this), Masks(this)
   {
       if (fDataCount<=0)
        {
          throw "TRnReleaseMultiDataSpaceModel constructor error: The data count value is incorrect.";
        }

       fMasks = new TMask * [fDataCount];
       for (int i=0;i<fDataCount;i++)
        {
           fMasks[i] = &( (lMasks[i])->Copy() );
        }
   }
//*************************************
  TRnReleaseMultiDataSpaceModel::~TRnReleaseMultiDataSpaceModel()
   {
      if (fMasks!=NULL)
       {
         for (int i=0;i<fDataCount;i++)
          {
            delete fMasks[i];
          }
         delete[] fMasks;
       }
   }

//#################################################################
  

//---------------------MultiData (End)-----------------------------------------





//---------------------1D (Begin)-----------------------------------------------


//#################################################################


//*************************************
   void TRnRelease1dSpaceModel::ReleaseConstructor()
    {
       //�������� �� ������������ ����� � �������� ����� (������)
          /*
            ��� ����� ������ ����������� � ���� ���������� ��������
            �����.
          */

		//������-���������: ����� ������������ ���� ��� ������� ������ ��� ����������� ������� �����������, �.�. ������������ - �� ���������
		//� �� ���������� = 1        
        if (fDataCount>1)
        {
          throw "TRnReleaseMultiDataSpaceModel constructor error: The data count value is incorrect.";
        }

        for (int ind=0;ind<fDataCount;ind++)
         {
                T1dNumberMask * pMask = dynamic_cast<T1dNumberMask * >(fMasks[ind]);
                if (pMask==NULL)
                 {
                    throw "TRnRelease1dSpaceModel constructor error: Could not cast one of the masks to T1dNumberMask * type.";
                 }

                 //�������� �� ������������ �������� ����� (������)
                   /*
                    ��� �������� ��������� ����� ������ ���� ���� 0 ���� 1.
                    ���� �������� �������� ����� ���� 1, ������ ���������������
                    ������� ������� ������������, �.�. ���������� � ���. � ��.
                    ���������. ���� �� �������� �������� ����� ���� 0, ��
                    ��������������� ������� ������� ���������.
                   */

                   T1dNumberMask & mask = *pMask;
                   int n = mask.N;
                   for (int i=0;i<n;i++)
                   {
                     
                     //������-���������: ����� ������������ ���� ��� ������� ������ ��� ����������� ������� �����������, �.�. ����� - �� ��������� � ��
                     //������ ��������� 0 
                     if (mask[i]!=1)
                     {
                         throw "TRnRelease1dSpaceModel constructor error:One of the mask element value is incorrect. Must be 1.";
                     }
					 /*
                     if ( (mask[i]!=0) && (mask[i]!=1) )
                     {
                         throw "TRnRelease1dSpaceModel constructor error:One of the mask element value is incorrect. Must be 0 or 1.";
                     }
					 */
                   }
                 //�������� �� ������������ �������� ����� (�����)

         }  // ���� �� ������ ������������.
         
        //�������� �� ������������ ����� � �������� ����� (�����)




        //������ �� ���������� ���������� ����� Rn (������)
         /*
           ���������� ��������� ����� ���� ��� "����� ��������� ������� Rn",
           � ����� �������� ������ � ��������� ������ � �������� ������������. 
         */
         int ind;
         fEndIndex = 0;
         int maxN = 0;
         for (ind=0;ind<fDataCount;ind++)
          {
              T1dNumberMask * pMask = dynamic_cast<T1dNumberMask * >(fMasks[ind]);
              T1dNumberMask & mask = *pMask;
              int n = mask.N;
              if (n>maxN) {maxN = n;}
              for (int i=0;i<n;i++)
               {
                 if (mask[i]>0)
                  {
                    fEndIndex = fEndIndex +1;
                  }
               }
          }

          fCorrespondence = new TCorrespondence *[fEndIndex];


          /*
           ��������� �� ���� ��������� �������������� ������� ����� ������������
           ��� ���������� ��������� ������������. �.�. ������ ���������� ������ ����
           ������ �� ����� �����, �������������� ����������� ������� �����. ������, �.�.
           ����� ������������ � ����������� ������� �� ������ ���������� ���� �� ����� ��
           ����� ����������, �� �� �� �������� ��������� � ������, ���� ��������� ������ ����������
           ������������ �����������.
          */


          T1dInverseCorrespondence &lInvCor = *new T1dInverseCorrespondence(fDataCount,maxN);

          /*
            ��� ���������� ��������� ������������ ����� ��� ����, ��� ��� ����������
            ������� ��� ��������. �.�. ����� ��������� ������ ������ ��������.
          */

          int num = 0;

          for (ind=0; ind<fDataCount;ind++)
          {
              T1dNumberMask * pMask = dynamic_cast<T1dNumberMask * >(fMasks[ind]);
              T1dNumberMask & mask = *pMask;
              int n = mask.N;

              for (int i=0;i<n;i++)
               {
                 if (mask[i]>0)
                  {
                     num = num +1;
                     fCorrespondence[num-1] = new T1dCorrespondence(ind,i);
                     lInvCor.fInvCor[ind][i] = num;
                  }
               }
          }

       fInverseCorrespondence = &lInvCor;

      //������ �� ���������� ���������� ����� Rn (�����)
    }
//*************************************
   TRnRelease1dSpaceModel::TRnRelease1dSpaceModel(const int lDataCount, const TMask * const * lMasks):
                           TRnReleaseMultiDataSpaceModel(lDataCount,lMasks)
     {
        ReleaseConstructor();
     }
//*************************************
   TMask * * TRnRelease1dSpaceModel::builtMasks = NULL;
   TMask *  TRnRelease1dSpaceModel::builtMask = NULL;
//*************************************
   const TMask * const * TRnRelease1dSpaceModel::BuildMasks(const int lDataCount, const int lN, const unsigned char *lMask)
     {
           T1dNumberMask & mask = *(new T1dNumberMask(lN));
           for (int i=0;i<lN;i++)
            {
               mask[i] = lMask[i];
            }

           builtMask = &mask;

           TMask * * masks = new TMask * [lDataCount];
           for (int ind=0; ind<lDataCount;ind++)
            {
               masks[ind] = &mask;
            }

           builtMasks = masks;

           return masks;
     }
//*************************************
   TRnRelease1dSpaceModel::TRnRelease1dSpaceModel(const int lDataCount,const int lN,const unsigned char *lMask):
                           TRnReleaseMultiDataSpaceModel(lDataCount,BuildMasks(lDataCount, lN,lMask))
    {

         delete builtMask;
         delete[] builtMasks;
		// builtMask = NULL;
		// builtMasks = NULL;
         ReleaseConstructor();
    }
//*************************************
   const T1dNumberMask & TRnRelease1dSpaceModel::GetMask(const int lInd)
    {
       if ( (lInd<0) || (lInd>fDataCount-1) )
            {
              throw "Error in TRnRelease1dSpaceModel::GetMask: The index value is out of range.";
            }
       return *(dynamic_cast<T1dNumberMask *>(fMasks[lInd]));
    }
//*************************************
   TRnRelease1dSpaceModel::~TRnRelease1dSpaceModel()
   {
     for (int i=0;i<fEndIndex;i++)
      {
         delete fCorrespondence[i];
      }
     delete[] fCorrespondence;

     delete fInverseCorrespondence;
   }
//*************************************
   TRnRelease1dSpaceModel &TRnRelease1dSpaceModel::Copy()
   {
        const TMask * const * lMasks = fMasks;
        return *( new TRnRelease1dSpaceModel(fDataCount,lMasks) );
   }
//*************************************

//#################################################################



//---------------------1D (End)-----------------------------------------------






//---------------------3D (Begin)-----------------------------------------------


//#################################################################

//*************************************
   void TRnRelease3dSpaceModel::ReleaseConstructor()
    {
      //�������� �� ������������ ����� � �������� ����� (������)
          /*
            ��� ����� ������ ����������� � ���� ���������� ��������
            �����.
          */
        for (int ind=0;ind<fDataCount;ind++)
         {
                T3dNumberMask * pMask = dynamic_cast<T3dNumberMask * >(fMasks[ind]);
                if (pMask==NULL)
                 {
                    throw "TRnRelease3dSpaceModel constructor error: Could not cast one of the masks to T3dNumberMask * type.";
                 }

                 //�������� �� ������������ �������� ����� (������)
                   /*
                    ��� �������� ��������� ����� ������ ���� ���� 0 ���� 1.
                    ���� �������� �������� ����� ���� 1, ������ ���������������
                    ������� ������� ������������, �.�. ���������� � ���. � ��.
                    ���������. ���� �� �������� �������� ����� ���� 0, ��
                    ��������������� ������� ������� ���������.
                   */

                   T3dNumberMask & mask = *pMask;
                   int n = mask.N;
                   int l = mask.L;
                   int m = mask.M;
                   for (int i=0;i<n;i++)
                    {
                     for (int j=0;j<l;j++)
                       {
                         for (int k=0;k<m;k++)
                           {
                             if ( (mask[i][j][k]!=0) && (mask[i][j][k]!=1) )
                               {
                                 throw "TRnRelease3dSpaceModel constructor error:One of the mask element value is incorrect. Must be 0 or 1.";
                               }
                           }
                       }
                    }
                 //�������� �� ������������ �������� ����� (�����)

         }  // ���� �� ������ ������������.

        //�������� �� ������������ ����� � �������� ����� (�����)


        //������ �� ���������� ���������� ����� Rn (������)
         /*
           ���������� ��������� ����� ���� ��� "����� ��������� ������� Rn",
           � ����� �������� ������ � ��������� ������ � �������� ������������.
         */
         int ind;
         fEndIndex = 0;
         int maxN = 0;
         int maxL = 0;
         int maxM = 0;
         for (ind=0; ind<fDataCount;ind++)
          {
              T3dNumberMask * pMask = dynamic_cast<T3dNumberMask * >(fMasks[ind]);
              T3dNumberMask & mask = *pMask;
              int n = mask.N;
              int l = mask.L;
              int m = mask.M;
              if (n>maxN) {maxN = n;}
              if (l>maxL) {maxL = l;}
              if (m>maxM) {maxM = m;}
              for (int i=0;i<n;i++)
               {
                 for (int j=0;j<l;j++)
                  {
                    for (int k=0;k<m;k++)
                     {
                        if (mask[i][j][k]>0)
                          {
                            fEndIndex = fEndIndex +1;
                          }
                     }
                  }

               }
          }


          fCorrespondence = new TCorrespondence *[fEndIndex];

          /*
           ��������� �� ���� ��������� �������������� ������� ����� ������������
           ��� ���������� ��������� ������������. �.�. ������, ������� � ��������� ����������� ������ ����
           ������� �� ����� �����, �������������� ����������� ������� �����. ������, �.�.
           ����� ������������ � ����������� ������� �� ������ ���������� ���� �� ����� ��
           ����� ����������, �� �� �� �������� ��������� � ������, ���� ��������� ������,
           ������� � ��������� ����������� ������������ �����������.
          */


          T3dInverseCorrespondence &lInvCor = *new T3dInverseCorrespondence(DataCount,maxN,maxL,maxM);

          /*
            ��� ���������� ��������� ������������ ����� ��� ����, ��� ��� ����������
            ������� ��� ��������. �.�. ����� ��������� ������ ������ ��������.
          */


          int num = 0;

          for (ind=0; ind<fDataCount;ind++)
          {
              T3dNumberMask * pMask = dynamic_cast<T3dNumberMask * >(fMasks[ind]);
              T3dNumberMask & mask = *pMask;
              int n = mask.N;
              int l = mask.L;
              int m = mask.M;

              for (int i=0;i<n;i++)
               {
                 for (int j=0;j<l;j++)
                  {
                     for (int k=0;k<m;k++)
                       {
                         if (mask[i][j][k]>0)
                           {
                             num = num +1;
                             fCorrespondence[num-1] = new T3dCorrespondence(ind,i,j,k);
                             lInvCor.fInvCor[ind][i][j][k] = num;
                           }
                       }
                  }
               }
          }

       fInverseCorrespondence = &lInvCor;

      //������ �� ���������� ���������� ����� Rn (������)
    }
//*************************************
   TRnRelease3dSpaceModel::TRnRelease3dSpaceModel(const int lDataCount, const TMask * const * lMasks):
                           TRnReleaseMultiDataSpaceModel(lDataCount,lMasks)
     {
        ReleaseConstructor();
     }
//*************************************
    TMask * * TRnRelease3dSpaceModel::builtMasks = NULL;
    TMask * TRnRelease3dSpaceModel::builtMask = NULL;
//*************************************
    const TMask * const * TRnRelease3dSpaceModel::BuildMasks(const int lDataCount, const int lN, const int lL, const int lM, const unsigned char * const * const * lMask)
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

           TMask * * masks = new TMask * [lDataCount];
           for (int ind=0; ind<lDataCount;ind++)
            {
               masks[ind] = &mask;
            }

           builtMasks = masks;

           return masks;
     }
//*************************************
   TRnRelease3dSpaceModel::TRnRelease3dSpaceModel(const int lDataCount, const int lN, const int lL, const int lM, unsigned char *** lMask)
                                                  : TRnReleaseMultiDataSpaceModel(lDataCount,BuildMasks(lDataCount,lN,lL,lM,lMask)) {
         delete builtMask;
         delete[] builtMasks;
         builtMask = NULL;
         builtMasks = NULL;

         ReleaseConstructor();
    }
//*************************************
   const T3dNumberMask & TRnRelease3dSpaceModel::GetMask(const int lInd)
    {
       if ( (lInd<0) || (lInd>fDataCount-1) )
            {
              throw "Error in TRnRelease3dSpaceModel::GetMask: The index value is out of range.";
            }
       return *(dynamic_cast<T3dNumberMask *>(fMasks[lInd]));
    }
//*************************************

   TRnRelease3dSpaceModel::~TRnRelease3dSpaceModel()
   {
     for (int i=0;i<fEndIndex;i++)
      {
         delete fCorrespondence[i];
      }
     delete[] fCorrespondence;

     delete fInverseCorrespondence;
   }
//*************************************
   TRnRelease3dSpaceModel &TRnRelease3dSpaceModel::Copy()
   {
      const TMask * const * lMasks = fMasks;
      return *( new TRnRelease3dSpaceModel(fDataCount,lMasks) );
   }
//*************************************

//#################################################################


//---------------------3D (End)-----------------------------------------------
