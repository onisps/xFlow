#pragma hdrstop

#include "XOperator.h"

#include <cstddef>
#include <stdio.h>
#include <math.h>
#include <omp.h>

//#################################################################

//*************************************
 /*
  *TLinearOperator::TLinearOperator(), Argument(this)
  * {
  *   fArgument = NULL;
  * }
  */
//*************************************
 TLinearCombinationWith2LinearOperatorOperands &TLinearOperator::operator()(const TSpace &s) const
   {
     if (&s==NULL)
      {
       throw "Error in TLinearOperator::operator(): The parameter pased is NULL.";
      }

     fArgument = const_cast<TSpace *>(&s);
     CheckCorrectCall();
     return * new TLinearCombinationWith2LinearOperatorOperands(*this);
   }
//*************************************
 TSpace &TLinearOperator::GetArgument() const
	{
     return *fArgument;
   }
//*************************************   

//#################################################################

//*************************************
 TRnLinearOperator::TRnLinearOperator(const int lModelName):TLinearOperator(),
        ActualIndexesCount(this), ActualIndexes(this), MaxRelatedIndexesCount(this),
        IndexesGroupsCount(this), IndexesGroupsTable(this), Matrix(this), ModelName(this)
   {


    //Временно-потом сделать нормально (Начало)
       helpService = NULL;
       helpServiceCount = 0;

       helpService1 = NULL;
       helpService1Count = 0;
    //Временно-потом сделать нормально (Конец)




     /*
      Необходимо будет создать механизм конструирования объектов класса
      через статическую функцию, с целью вызова при конструировании
      объектов безопасных функций инициализации, таких как заполнение таблицы
      активности, таблицы групп индексов и т.д.
     */

     fModelName = lModelName;
     TRnSpaceModel *model;
     model=TRnSpace::GetModel(lModelName);
     if (model!=NULL)
        {
          fModel=model;
        }
     else
        {
          throw "Error in TRnLinearOperator constructor: Model not found.";
        }

     int N = fModel->EndIndex;
     

	 fActivityTable = NULL;
     fExternalActivityTable = false;
	 fActivityTableConjugate = NULL;
     fExternalActivityTableConjugate = false;

     fIndexesGroupsCount = 0;
     fIndexesGroupsTable = NULL;

     fActualIndexes = NULL;

     fMatrix = NULL;
   }
//*************************************
 void TRnLinearOperator::DeleteTable(bool lExternalTable, int ** lTable)
 {     
	 if ( (lExternalTable==false) && (lTable!=NULL)	)
	 {
        int N = fModel->EndIndex;
        for(int i=0;i<N;i++)
        {
          if (lTable[i]!=NULL)
          {
              delete[] lTable[i];
          }
        }	  
		delete[] lTable;// lTable = NULL;
	 }
 }
 //*************************************
 TRnLinearOperator::~TRnLinearOperator()
   { 
	  if (fActivityTableConjugate == fActivityTable) {fActivityTableConjugate=NULL;}
	  DeleteTable(fExternalActivityTable,fActivityTable);
	  DeleteTable(fExternalActivityTableConjugate,fActivityTableConjugate);
	
      int count = fIndexesGroupsCount;
      if (count>0)
       {
         for(int i=0;i<count;i++)
          {
            delete[] fIndexesGroupsTable[i];
          }

		 delete[] fIndexesGroupsTable; fIndexesGroupsTable = NULL;
       }


      if (fActualIndexes!=NULL)
        {
         delete[] fActualIndexes;//  fActualIndexes = NULL;
        }

	  int N = fModel->EndIndex;

      if (fMatrix!=NULL)
       {
          for(int i=0;i<N;i++)
           {
              if (fMatrix[i]!=NULL)
                {
                   delete[] fMatrix[i];
                }
           }
		  delete[] fMatrix;// fMatrix = NULL;
       }


       //Временно-потом сделать нормально (Начало)
       if (helpServiceCount>0)
        {
          for (int i=0;i<helpServiceCount;i++)
           {
             delete helpService[i];
           }
		  delete[] helpService;// helpService = NULL;
        }

       if (helpServiceCount>0)
        {
		   delete[] helpService1;// helpService1 = NULL;
        }  
       //Временно-потом сделать нормально (Начало)

   }
//*************************************
 void TRnLinearOperator::ConstructActivityTables()
  {
	  throw "Error in TRnLinearOperator::ConstructActivityTables: The method is deprecated.";

      TRnLinearOperator &A = *this;
      TSpace &s = NewArgument();
      /*
       Далее работаем только с полученным вектором и не используем NewArgument...
       Конечно надо подумать как вообще избавиться от этого метода - некое дублирование
       наблюдается...
      */
      TRnSpace *spRn = dynamic_cast<TRnSpace *>(&s);
      if (spRn==NULL)
       {
         throw "Error in TRnLinearOperator::ConstructActivityTables: Could not cast new argument to TRnSpace * type.";
       }
      TRnSpace &sRn = *spRn;

      //Теперь имеем линейный оператор, действующий над Rn и вектор Rn для всех действий...
      //Можно приступать к созданию таблиц...

       int N = sRn.EndIndex;
	   int i;


	   if (fActivityTableConjugate == fActivityTable) {fActivityTableConjugate=NULL;}
	   DeleteTable(fExternalActivityTable,fActivityTable);
	   DeleteTable(fExternalActivityTableConjugate,fActivityTableConjugate);
	

	   fActivityTable = new int *[N];       
       for(i=0;i<N;i++)
       {
         fActivityTable[i] = NULL;
       }



	   TRnSpace &z = sRn.CreateInstance();
       TRnSpace &Az = sRn.CreateInstance();

       for (int i=1;i<N+1;i++) {z[i]=0;Az[i]=0;}

       for (int i=1;i<N+1;i++)
        {
              z[i]  = 1;
              if (i!=1) {z[i-1] =0; }
              Az|=A(z);

              int countAct = 0;
              for (int j=1;j<N+1;j++)
               {
                  if (Az[j]!=0) {countAct=countAct+1;}
               }

               
              fActivityTable[i-1] = new int[countAct+1];
              fActivityTable[i-1][0]  = countAct;

              int num = 0;
              for (int j=1;j<N+1;j++)
               {
                  if (Az[j]!=0)
                   {
                     num = num + 1;
                     fActivityTable[i-1][num] = j;
                   }
               }

             //FormMain->StatusBar1->SimpleText = IntToStr(i);

        }

       delete &z;
       delete &Az;

       delete &s;

	   fExternalActivityTable = false;

	   fActivityTableConjugate = fActivityTable;
	   fExternalActivityTableConjugate = false;
  }
//*************************************
void TRnLinearOperator::SetActivityTable(int** lActivityTable)
{
	 DeleteTable(fExternalActivityTable,fActivityTable);	
	 fActivityTable = lActivityTable; 
	 fExternalActivityTable = true;
}  
//*************************************
void TRnLinearOperator::SetActivityTableConjugate(int** lActivityTableConjugate)
{
	 DeleteTable(fExternalActivityTableConjugate,fActivityTableConjugate);
	 fActivityTableConjugate = lActivityTableConjugate; 
	 fExternalActivityTableConjugate = true;
}  
//*************************************
 void TRnLinearOperator::ConstructIndexesGroupsTable()
  {
      throw "Error in TRnLinearOperator::ConstructIndexesGroupsTable: The method is deprecated.";
	  /*
        Реализован механизм автоматического заолнения таблицы
        групп индексов с помощью алгоритма уплотнения.
      */

       //Подсчет количества групп...
       int N = fModel->EndIndex;

       int * baseVector = new int[N+1]; //базовый вектор для построения групп.
       baseVector[0]=-1; //не используется.

       int * groupVector = new int[N+1]; //вектор для каждой группы
       groupVector[0]=-1; //не используется...

       for (int i=1;i<N+1;i++)
        {
           baseVector[i]=0;
           groupVector[i]=0;
        }

       int ** AT = fActivityTable;
       int groupsCount = 0;
       int i; int j; int num;

       for (i=1;i<N+1;i++)
        {
          if (baseVector[i]==0)
           {
              groupsCount = groupsCount + 1;
              for (j=1;j<N+1;j++) {groupVector[j]=0;}
              int countIndexes = 0;
              int firstIndex;
              for (j=i;j<N+1;j++)
               {
                 if (baseVector[j]==0)
                   {
                     int countRI = AT[j-1][0];
                     bool flagAcceptIntoGroup=true;
                     for (num=1;num<countRI+1;num++)
                      {
                        if (groupVector[AT[j-1][num]]!=0) {flagAcceptIntoGroup=false;}
                      }//Цикл по связанным индексам из ТА оператора
                     if (flagAcceptIntoGroup==true)
                      {
                        countIndexes = countIndexes + 1;
                        if (countIndexes==1) {firstIndex=j;}
                        baseVector[j]=groupsCount;
                        for (num=1;num<countRI+1;num++)
                         {
                           groupVector[AT[j-1][num]]=1;
                         }//Цикл по связанным индексам из ТА оператора
                      }//Если можно включить в тек. группу
                   }//Если еще не принадлежит никакой другой группе
               }//Цикл по групповому вектору
              /*
               Теперь надо где-то зафиксировать количество индексов в текущей группе...
              */
              baseVector[firstIndex] = N*groupsCount+countIndexes;
           } //Если начинается новая группа
        } //Цикл по базовому вектору

       /*
         Теперь на основе базового вектора заполним
         таблицу групп индексов...
       */
       fIndexesGroupsTable = new int *[groupsCount];

       int currentCroup=0;
       for (i=1;i<N+1;i++)
        {
           if (baseVector[i]>N) //Начинается новая группа
            {
              currentCroup = currentCroup +1;
              int M = baseVector[i];
              int gr = M/N;
              int countIndexes = M-N*gr;

              if (gr!=currentCroup)
               {
                 throw "Error in TRnLinearOperator::ConstructIndexesGroupsTable: Groups numbers are not equal.";
               }//Перестрахуемся - номера групп должны совпадать

              fIndexesGroupsTable[gr-1] = new int[countIndexes+1];
              fIndexesGroupsTable[gr-1][0] = countIndexes;
              fIndexesGroupsTable[gr-1][1] = i;

              num = 1;
              for (j=i+1;j<N+1;j++)
               {
                 if (baseVector[j]==gr)
                  {
                    num=num+1;
                    fIndexesGroupsTable[gr-1][num] = j;
                  }
               }

              if (num!=countIndexes)
               {
                  throw "Error in TRnLinearOperator::ConstructIndexesGroupsTable: Indexes counts are not equal.";
               }//Перестрахуемся - количества индексов должны совпадать

            }//Если новая группа
        }//Цикл по базовому вектору

        if (currentCroup!=groupsCount)
          {
            throw "Error in TRnLinearOperator::ConstructIndexesGroupsTable: Groups counts are not equal.";
          }//Перестрахуемся - количества групп должны совпадать

       fIndexesGroupsCount = groupsCount;

       //ShowMessage("Auto: "+IntToStr(groupsCount));

	   delete[] baseVector;// baseVector = NULL;
	   delete[] groupVector;// groupVector = NULL;
  }
//*************************************
 void TRnLinearOperator::ConstructMatrix() //Построение матрицы оператора.
  {      	  
	  TRnLinearOperator &A = *this;
      TSpace &s = NewArgument();
      /*
       Далее работаем только с полученным вектором и не используем NewArgument...
       Конечно надо подумать как вообще избавиться от этого метода - некое дублирование
       наблюдается...
      */
      TRnSpace *spRn = dynamic_cast<TRnSpace *>(&s);
      if (spRn==NULL)
       {
         throw "Error in TRnLinearOperator::ConstructMatrix: Could not cast new argument to TRnSpace * type.";
       }
      TRnSpace &sRn = *spRn;

      //Теперь имеем линейный оператор, действующий над Rn и вектор Rn для всех действий...
      //Можно приступать к созданию таблиц...

       int N = sRn.EndIndex;

       TRnSpace &z = sRn.CreateInstance();
       TRnSpace &Az = sRn.CreateInstance();


       int i; int j;

	   if (fMatrix!=NULL)
       {
          for(i=0;i<N;i++)
           {
              if (fMatrix[i]!=NULL)
                {
                   delete[] fMatrix[i];
                }
           }
		  delete[] fMatrix;// fMatrix = NULL;
       }
	  


       fMatrix = new double *[N];
       for (i=0;i<N;i++)
        {
          fMatrix[i] = new double [N];
        }

       for (i=1;i<N+1;i++)
        {
              z[i] = 1;
              if (i!=1) {z[i-1] =0;}
              Az|=A(z);
              for (j=1;j<N+1;j++)
               {
                 fMatrix[j-1][i-1] = Az[j];
               }
        }

        //Временно вывод матрицы (Начало)
        /*
        TStringGrid &sg = *(FormMatrix->StringGridMatrix);
        sg.RowCount = N+1;
        sg.ColCount = N+1;
        for (i=1;i<N+1;i++)
        {
              for (j=1;j<N+1;j++)
               {
                 double value = fMatrix[i-1][j-1];
                 sg.Cells[j][i] = FloatToStr(value);
               }
        }
        ShowMessage("Matrix dimension: "+IntToStr(N));

        for (i=1;i<N+1;i++)
        {
           const TCorrespondence &cor = *(A.GetModelCorrespondence(i));

           try
           {
           //2d (Начало)
           const T2dCorrespondence &cor2d = dynamic_cast<const T2dCorrespondence &>(cor);
           sg.Cells[0][i]=FloatToStr(cor2d.I)+":"+FloatToStr(cor2d.J);
           sg.Cells[i][0]=FloatToStr(cor2d.I)+":"+FloatToStr(cor2d.J);
           //2d (Конец)
           }
           catch (...)
           {
           //3d (Начало)
           const T3dCorrespondence &cor3d = dynamic_cast<const T3dCorrespondence &>(cor);
           sg.Cells[0][i]=FloatToStr(cor3d.I)+":"+FloatToStr(cor3d.J)+":"+FloatToStr(cor3d.K);
           sg.Cells[i][0]=FloatToStr(cor3d.I)+":"+FloatToStr(cor3d.J)+":"+FloatToStr(cor3d.K);
           //3d (Конец)
           }


        }

        FormMatrix->Show();
        */
        //Временно вывод матрицы (Конец)

        delete &z;
        delete &Az;

        delete &s;
  }
//*************************************
 const double * const * const TRnLinearOperator::GetMatrix() const //Получить матрицу.
  {
     return fMatrix;
  }
//*************************************
 double TRnLinearOperator::GetMatrixElement(const int li,const int lj) const //Получить элемент матрицы.
  {
     int N = fModel->EndIndex;
     if ( (li<1) || (li>N) || (lj<1) || (lj>N) )
      {
        throw "TRnLinearOperator::GetMatrixElement: The index is out of range.";
      }

     if (fMatrix==NULL)
      {
        throw "TRnLinearOperator::GetMatrixElement: The matrix is not constructed.";
      }

     return fMatrix[li-1][lj-1];
  }
//*************************************
 void TRnLinearOperator::CheckCorrectCall() const
  {
     const TRnSpace *s = NULL;
     //Попытка приведения аргумента к TRnSpace *
         try {s=dynamic_cast<const TRnSpace *>(fArgument);}
         catch (...) {throw "Error in TRnLinearOperator::CheckCorrectCall: Could not cast call argument to TRnSpace & type.";}
     //Конец Попытка приведения аргумента к TRnSpace *

     if (fModel!=s->fModel)
     {throw "TRnLinearOperator::CheckCorrectCall: models mismatch.";}


     /*
       Ниже производится заполнение актуальных индексов - это надо конечно
       же перенести в другой метод...
     */

     fActualIndexesCount = 0;
     if (fActualIndexes!=NULL)
      {
        delete[] fActualIndexes; fActualIndexes = NULL;
      }

     int i; int j;
     int aic = s->ActualIndexesCount;
     int *ai = s->ActualIndexes;
     if (aic>0)
     {
		 if (fActivityTable==NULL)
		 {
           throw "TRnLinearOperator::CheckCorrectCall: Actual indexes can't be unsed. Activity table is not ready.";
		 }

         for (i=0;i<aic;i++)
          {
            fActualIndexesCount = fActualIndexesCount + fActivityTable[ai[i]-1][0];
          }
     }

     if (fActualIndexesCount>0)
      {
        fActualIndexes = new int[fActualIndexesCount];
        int num=-1;
        for (i=0;i<aic;i++)
          {
             int tc = fActivityTable[ai[i]-1][0];
             int *ti = fActivityTable[ai[i]-1];
             for (j=1;j<=tc;j++)
              {
                num = num + 1; 
                fActualIndexes[num] = ti[j]; // ? Вообще желательно проверять на
                                             //повторения... - надо сделать !!!
              }
          }
      }

     LocalCheckCorrectCall();
  }
//*************************************
 const TCorrespondence *TRnLinearOperator::GetModelCorrespondence(int li) const
    {
      if ( (li<1) || (li>fModel->EndIndex) )
       {
        throw "TRnLinearOperator::GetModelCorrespondence: index is out of range.";
       }

      return &(fModel->GetCorrespondence(li));
    }
//*************************************
 const TInverseCorrespondence *TRnLinearOperator::GetInverseModelCorrespondence() const
    {
      return &(fModel->GetInverseCorrespondence());
    }
//*************************************
 int *TRnLinearOperator::GetRelatedIndexes(int li) const
   {
     if ( (li<1) || (li>fModel->EndIndex) )
       {
        throw "TRnLinearOperator::GetRelatedIndexes: index is out of range.";
       }
     return fActivityTable[li-1];
   }
//*************************************
 int TRnLinearOperator::GetMaxRelatedIndexesCount() const
   {
      int res = 0;
      int i;
      int N = fModel->EndIndex;
      int localRes;
      for (i=1;i<N+1;i++)
       {
          localRes = fActivityTable[i-1][0];
          if (localRes>res) {res = localRes;}
       }
      return res;
   }
//*************************************

//#################################################################

//*************************************
 TRnRelease1dLinearOperator::TRnRelease1dLinearOperator(const int lModelName):TRnLinearOperator(lModelName)
  {
     fCastArgument = NULL;

     //Попытка приведения модели к TRnRelease1dSpaceModel *
         try {fCastModel=dynamic_cast<TRnRelease1dSpaceModel *>(fModel);}
         catch (...) {throw "Error in TRnRelease1dLinearOperator constructor: Could not cast model to TRnRelease1dSpaceModel * type.";}
     //Конец Попытка приведения модели к TRnRelease1dSpaceModel *
  }
//*************************************
 void TRnRelease1dLinearOperator::LocalCheckCorrectCall()  const
   {
     //Попытка приведения аргумента к TRnRelease1dSpace *
         try {fCastArgument=dynamic_cast<TRnRelease1dSpace *>(fArgument);}
         catch (...) {throw "Error in TRnRelease1dLinearOperator::LocalCheckCorrectCall: Could not cast call argument to TRnRelease1dSpace * type.";}
     //Конец Попытка приведения аргумента к TRnRelease1dSpace *

     ComfortableCall();
   }
//*************************************
 TSpace &TRnRelease1dLinearOperator::NewArgument() const
  {
     return *(new TRnRelease1dSpace(fModelName));
  }
//*************************************

//#################################################################

//*************************************
 TRnRelease3dLinearOperator::TRnRelease3dLinearOperator(const int lModelName):TRnLinearOperator(lModelName)
  {
     fCastArgument = NULL;

     //Попытка приведения модели к TRnRelease3dSpaceModel *
         try {fCastModel=dynamic_cast<TRnRelease3dSpaceModel *>(fModel);}
         catch (...) {throw "Error in TRnRelease3dLinearOperator constructor: Could not cast model to TRnRelease3dSpaceModel * type.";}
     //Конец Попытка приведения модели к TRnRelease3dSpaceModel *
  }
//*************************************
 void TRnRelease3dLinearOperator::LocalCheckCorrectCall()  const
   {
     //Попытка приведения аргумента к TRnRelease3dSpace *
         try {fCastArgument=dynamic_cast<TRnRelease3dSpace *>(fArgument);}
         catch (...) {throw "Error in TRnRelease3dLinearOperator::LocalCheckCorrectCall: Could not cast call argument to TRnRelease3dSpace * type.";}
     //Конец Попытка приведения аргумента к TRnRelease3dSpace *

     ComfortableCall();
   }
//*************************************
 TSpace &TRnRelease3dLinearOperator::NewArgument() const
  {
     return *(new TRnRelease3dSpace(fModelName));
  }
//*************************************

//#################################################################

//*************************************
TFiniteDifferenceOperator::TFiniteDifferenceOperator(const int lModelName, const int lGridsCount, const TNodesGrid* const* lGrids):
                           TRnLinearOperator(lModelName),fGridsCount(lGridsCount)
{
      if (fGridsCount<=0)
       {
         throw "Error in TFiniteDifferenceOperator constructor: GridsCount<=0.";
       }

	  fGrids = new TNodesGrid*[fGridsCount];
	  for (int i=0;i<fGridsCount;i++)
	  {
		  try
		  {
			  TNodesGrid& gr = lGrids[i]->Copy();
			  fGrids[i] = &gr;
		  }
		  catch(...)
		  {
             throw "Error in TFiniteDifferenceOperator constructor: Error while processing grids.";
		  }
	  }
}
//*************************************
TFiniteDifferenceOperator::~TFiniteDifferenceOperator()
{
	  for (int i=0;i<fGridsCount;i++)
	  {
		  delete fGrids[i];
	  }
	  delete fGrids;
}
//*************************************

//#################################################################

//*************************************
T1dLaplaceOperator::T1dLaplaceOperator(const int lModelName, const int lGridsCount, const TNodesGrid* const* lGrids):
                        TRnRelease1dLinearOperator(lModelName),
                        TFiniteDifferenceOperator(lModelName,lGridsCount,lGrids),
                        TRnLinearOperator(lModelName)
  {
    if (fGridsCount!=1) 
	{
		throw "Error in T1dLaplaceXOperator constructor: GridCount<>1.";
	}

	const T1dNormalGrid* grP = NULL;
    try 
	{ 
		const T1dNormalGrid& gr = dynamic_cast<const T1dNormalGrid&>(*(fGrids[0])); 
		mask = &(gr.GetMask());
	    fNormals = gr.Normals;
		grP = &gr;
	}
    catch (...)
    {
      throw "Error in T1dLaplaceXOperator constructor: Could not cast grid to T1dNormalGrid* type.";
    }


	const T1dNormalGrid& gr = *grP;


	fIsUniformGrid = true;    
   
	const TSeparator &sep1 = gr.Separator1;
    try 
	{ 
		const TUniformSeparator& su = dynamic_cast<const TUniformSeparator& >(sep1); 
		hx = su.SeparationValue;
        hxhx = hx*hx;
	}
    catch (...)
    {
      fIsUniformGrid = false;
    }
        	
	if (fIsUniformGrid == false)
	{
		hx = 0;
        hxhx = 0;        
	}
  }
//*************************************
 void T1dLaplaceOperator::ComfortableCall() const
  {
    TRnRelease1dData *Rn1dData = dynamic_cast<TRnRelease1dData *>(fCastArgument->MultiData[0]);
    f = Rn1dData->Data;
  }
//*************************************
 double T1dLaplaceOperator::operator()(int idData,int i) const
  {
   //operatr is: -L(f) for inner points
   //            df/dn for normal border points, where n - is unitary external normal vector	             
   if (fIsUniformGrid == false)
   {
	 throw "Error in T1dLaplaceOperator::operator(): Operator is not implemented for non uniform grid.";
   }

   if (idData==0)
   {
    const int mijk = (*mask)[i];
    if (mijk == TActualPoint)
    {
     return
     -(

       ( (f[i+1]-2*f[i]+f[i-1])/(hxhx) ) 	 

     );	   
    }
    else if  (mijk == TNormalBorderPoint)
    {
	   const T1dVector* normal = fNormals[i];
	   if (normal==NULL) {throw "Error in T1dLaplaceOperator::operator(): Normal vector is undefined for this point.";}
	   const double nx = normal->X;	   
	   
	   double dfx=0;
	   if (nx>0) {dfx = (f[i]-f[i-1])/hxhx;} 
	   else if (nx<0) {dfx = (f[i+1]-f[i])/hxhx;}

	   return dfx*nx;
    }
    else
    {
       throw "Error in T1dLaplaceOperator::operator(): Operator is not implemented for this mask value.";
    }
   }
   else
   {
       throw "Error in T1dLaplaceOperator::operator(): Operator is not implemented for multidata index <> 0.";
   }

  }
//*************************************

//#################################################################

//*************************************
T3dLaplaceOperator::T3dLaplaceOperator(const int lModelName, const int lGridsCount, const TNodesGrid* const* lGrids):
                        TRnRelease3dLinearOperator(lModelName),
                        TFiniteDifferenceOperator(lModelName,lGridsCount,lGrids),
                        TRnLinearOperator(lModelName)
{ 
    if (fGridsCount!=1) 
	{
		throw "Error in T3dLaplaceOperator constructor: GridCount<>1.";
	}	
    try 
	{ 
		const T3dNormalGrid& gr = dynamic_cast<const T3dNormalGrid&>(*(fGrids[0])); 
		fGrid = &gr;
		fMask = &(gr.GetMask());
	    fNormals = gr.Normals;		
	}
    catch (...)
    {
      throw "Error in T3dLaplaceXOperator constructor: Could not cast grid to T3dNormalGrid* type.";
    }
}
//*************************************
 void T3dLaplaceOperator::ComfortableCall() const
  {
    TRnRelease3dData *Rn3dData = dynamic_cast<TRnRelease3dData *>(fCastArgument->MultiData[0]);
    f = Rn3dData->Data;
  }
//*************************************
void T3dLaplaceOperator::ConstructActivityTables()
{      
	 int RN = fModel->EndIndex;
	 
	 //Deleting old (Begin)
	 if (fActivityTableConjugate == fActivityTable) {fActivityTableConjugate=NULL;}
	 DeleteTable(fExternalActivityTable,fActivityTable);
	 DeleteTable(fExternalActivityTableConjugate,fActivityTableConjugate);	
     //Deleting old (End)

	 //Preparing new (Begin)

	 int** table = new int *[RN];       
     for(int i=0;i<RN;i++)
     {
         table[i] = NULL;
     }     
	 //Preparing new (End)


	 //Constructing (Begin)
	 TRnRelease3dSpaceModel& fCastModel = dynamic_cast<TRnRelease3dSpaceModel &>(*fModel);
	 
	 int elementsCount = 7; //7-point pattern of operator

	 for(int rn=1;rn<=RN;rn++)
     {
        
        const TCorrespondence* cor = GetModelCorrespondence(rn);
		const TInverseCorrespondence* invcor = GetInverseModelCorrespondence();
		
		//Do not check correction of cast types
		const T3dCorrespondence &cor3d = dynamic_cast<const T3dCorrespondence &>(*cor);
		const T3dInverseCorrespondence &invcor3d = dynamic_cast<const T3dInverseCorrespondence &>(*invcor);
        
		int ind = cor3d.Ind;                                       
        if (ind!=0) {throw "Error in T3dLaplaceOperator::ConstructActivityTables(): Unknown index od multidata.";}
		int i = cor3d.I;
        int j = cor3d.J;
        int k = cor3d.K;

		const T3dNumberMask & mask = fCastModel.GetMask(ind);
		
		int N = mask.N;
		int L = mask.L;
		int M = mask.M;

        table[rn-1] = new int[elementsCount+1];

	    int count = 1; 			
		table[rn-1][count] = rn;        
		for (int loopVar=1;loopVar<=elementsCount-1;loopVar++)
		{
           int i_1=-1; int j_1=-1; int k_1=-1;
		   
		   if (loopVar==1) {i_1=i-1; j_1=j; k_1=k;}
		   else if (loopVar==2) {i_1=i+1; j_1=j; k_1=k;}

		   else if (loopVar==3) {i_1=i; j_1=j-1; k_1=k;}
		   else if (loopVar==4) {i_1=i; j_1=j+1; k_1=k;}

		   else if (loopVar==5) {i_1=i; j_1=j; k_1=k-1;}
		   else if (loopVar==6) {i_1=i; j_1=j; k_1=k+1;}
		   else {throw "Error in T3dLaplaceOperator::ConstructActivityTables(): Unknown loop variable.";}

		   if ( (i_1>=0) && (i_1<N) && (j_1>=0) && (j_1<L) && (k_1>=0) && (k_1<M) )
		   {
			   if (mask[i_1][j_1][k_1]>0)
			   {
                    count++;
					table[rn-1][count] = invcor3d[ind][i_1][j_1][k_1];
			   }
		   }		   		   
		}

		table[rn-1][0] = count;
     }//by Rn index	 

	 //Constructing (End)

	 fActivityTableConjugate = table;
	 fExternalActivityTableConjugate = false;
	 

	 bool zeroSymmetry = true; //temporary forced

	 if (zeroSymmetry==true) {fActivityTable = fActivityTableConjugate;}
	 else
	 {
		  //Notes:
          // 1. Build from fActivityTableConjugate.
		  // 2. Parallel.
		  // 3. Check correct both tables.
		  // 4. Compare tables - it is interesting:
          //	 What operator is generating different tables?
		  //     How many elements we need to allocate for each row in this case?
	 }//if zero symmetry

	 fExternalActivityTable = false;
}
//*************************************
double T3dLaplaceOperator::operator()(int idData,int i,int j,int k) const
{
   if (idData==0)
   {
    const double* hx = fGrid->GetSeparator1().Separation;
	const double* hy = fGrid->GetSeparator2().Separation;
	const double* hz = fGrid->GetSeparator3().Separation;

	const T3dNumberMask& mask = *fMask;
	const int mijk = (*fMask)[i][j][k];	
	

	int stepI_;
	int stepI_1;
	double hxl;
	   
	int stepJ_;
	int stepJ_1;	   
    double hyl;
	   
	int stepK_;
	int stepK_1;	   	   
	double hzl;

	double fl;
	double fl_1;	   

	double dfx;
	double dfy;
	double dfz;

	double res;

    if (mijk == TActualPoint)
    {              
	   double h_1X; double hX; double hx2X; double hx_2X; double hxhxhX;  
	   h_1X = hx[i-1];
	   hX = hx[i];
	   hx2X = hx[i-1]+hx[i];
	   hx_2X = (hx[i-1]+hx[i])/2;
	   hxhxhX = hx[i-1]*hx[i]*(hx[i-1]+hx[i])/2;

	   double h_1Y; double hY; double hx2Y; double hx_2Y; double hxhxhY;  
	   h_1Y = hy[j-1];
	   hY = hy[j];
	   hx2Y = hy[j-1]+hy[j];
	   hx_2Y = (hy[j-1]+hy[j])/2;
	   hxhxhY = hy[j-1]*hy[j]*(hy[j-1]+hy[j])/2;

       double h_1Z; double hZ; double hx2Z; double hx_2Z; double hxhxhZ;  
	   h_1Z = hz[k-1];
       hZ = hz[k];
	   hx2Z = hz[k-1]+hz[k];
	   hx_2Z = (hz[k-1]+hz[k])/2;
	   hxhxhZ = hz[k-1]*hz[k]*(hz[k-1]+hz[k])/2;

       double f_i_1 = (f[i][j][k] - f[i-1][j][k])/(h_1X);
	   double f_i = (f[i+1][j][k] - f[i][j][k])/(hX);;
	   double f_j_1 = (f[i][j][k] - f[i][j-1][k])/(h_1Y);
	   double f_j = (f[i][j+1][k] - f[i][j][k])/(hY);				
       double f_k_1 = (f[i][j][k] - f[i][j][k-1])/(h_1Z);
	   double f_k = (f[i][j][k+1] - f[i][j][k])/(hZ); 
				

	   //Pre-border points (Begin)
	   if (mask[i-1][j][k]==TFictivePoint) {f_i_1 = 0;}
	   if (mask[i+1][j][k]==TFictivePoint) {f_i = 0;}
	   if (mask[i][j-1][k]==TFictivePoint) {f_j_1 = 0;}
	   if (mask[i][j+1][k]==TFictivePoint) {f_j = 0;}
	   if (mask[i][j][k-1]==TFictivePoint) {f_k_1 = 0;}
	   if (mask[i][j][k+1]==TFictivePoint) {f_k = 0;}
	   //Pre-border points (End)

	   //*
	   return
       -(
	       (f_i - f_i_1)/hx_2X
		   +
		   (f_j - f_j_1)/hx_2Y
		   +
		   (f_k - f_k_1)/hx_2Z
	    );
	   //*/


	   /*	
	   return
       -(
	       (hX*f[i-1][j][k] - hx2X*f[i][j][k] + h_1X*f[i+1][j][k])/hxhxhX 
		   +
		   (hY*f[i][j-1][k] - hx2Y*f[i][j][k] + h_1Y*f[i][j+1][k])/hxhxhY
		   +
		   (hZ*f[i][j][k-1] - hx2Z*f[i][j][k] + h_1Z*f[i][j][k+1])/hxhxhZ
	    );
	   //*/
    }
	else if  (mijk == TPreDefinedBorderPoint)
	{
	   const T3dVector* normal = fNormals[i][j][k];
	   if (normal==NULL) {throw "Error in T3dLaplaceOperator::operator(): Normal vector is undefined for this point.";}
	   const double nx = normal->X;
	   const double ny = normal->Y;
	   const double nz = normal->Z;


	   if ( (nx==0) && (ny==0) && (nz==0) ) {throw "Error in T3dLaplaceOperator::operator(): Normal vector is zero for this point.";}
	   else if (nx!=0)
	   {
         if ( (ny!=0) || (nz!=0) || (fabs(nx)!=1) )  {throw "Error in T3dLaplaceOperator::operator(): Normal vector is incorrect for this point.";}
	   }
	   else if (ny!=0)
	   {
         if ( (nx!=0) || (nz!=0) || (fabs(ny)!=1) ) {throw "Error in T3dLaplaceOperator::operator(): Normal vector is incorrect for this point.";}
	   }
	   else if (nz!=0)
	   {
         if ( (nx!=0) || (ny!=0) || (fabs(nz)!=1) )  {throw "Error in T3dLaplaceOperator::operator(): Normal vector is incorrect for this point.";}
	   }

	   double semiSum;

	   if (nx>0) {semiSum = (f[i][j][k]+f[i-1][j][k])/2;} 
	   else if (nx<0) {semiSum = (f[i+1][j][k]+f[i][j][k])/2;}

	   if (ny>0) {semiSum = (f[i][j][k]+f[i][j-1][k])/2;}
	   else if (ny<0) {semiSum = (f[i][j+1][k]+f[i][j][k])/2;}
	   
	   if (nz>0) {semiSum = (f[i][j][k]+f[i][j][k-1])/2;} 
	   else if (nz<0) {semiSum = (f[i][j][k+1]+f[i][j][k])/2;}

	   return semiSum;
	}
    else if  (mijk == TPreNormalBorderPoint)
    {
	   //throw "Error in T3dLaplaceOperator::operator(): Operator is not implemented for TPreNormalBorderPoint mask value.";

	   const T3dVector* normal = fNormals[i][j][k];
	   if (normal==NULL) {throw "Error in T3dLaplaceOperator::operator(): Normal vector is undefined for this point.";}
	   const double nx = normal->X;
	   const double ny = normal->Y;
	   const double nz = normal->Z;

	   if ( (nx==0) && (ny==0) && (nz==0) ) {throw "Error in T3dLaplaceOperator::operator(): Normal vector is zero for this point.";}

	   
	   stepI_=0;
	   stepI_1=0;
	   hxl=1;
	   if (nx>0) {stepI_=0; stepI_1=-1; hxl=hx[i-1];}
	   else if (nx<0) {stepI_=1; stepI_1=0; hxl=hx[i];}

	   
	   stepJ_=0;
	   stepJ_1=0;	   
	   hyl=1;
	   if (ny>0) {stepJ_=0; stepJ_1=-1; hyl=hy[j-1];}
	   else if (ny<0) {stepJ_=1; stepJ_1=0; hyl=hy[j];}

	   stepK_=0;
	   stepK_1=0;	   	   
	   hzl=1;
	   if (nz>0) {stepK_=0; stepK_1=-1; hzl=hz[k-1];}
	   else if (nz<0) {stepK_=1; stepK_1=0; ;hzl=hz[k];}

	   
	   
	   fl   = (f[i+stepI_][j][k] + f[i+stepI_][j+stepJ_][k] + f[i+stepI_][j][k+stepK_] + f[i+stepI_][j+stepJ_][k+stepK_])/4;
	   fl_1 = (f[i+stepI_1][j][k] + f[i+stepI_1][j+stepJ_][k] + f[i+stepI_1][j][k+stepK_] + f[i+stepI_1][j+stepJ_][k+stepK_])/4;
	   dfx = (fl - fl_1)/hxl;
	   
	   fl   = (f[i][j+stepJ_][k] + f[i+stepI_][j+stepJ_][k] + f[i][j+stepJ_][k+stepK_] + f[i+stepI_][j+stepJ_][k+stepK_])/4;
	   fl_1 = (f[i][j+stepJ_1][k] + f[i+stepI_][j+stepJ_1][k] + f[i][j+stepJ_1][k+stepK_] + f[i+stepI_][j+stepJ_1][k+stepK_])/4;
	   dfy = (fl - fl_1)/hyl;
	   
	   fl   = (f[i][j][k+stepK_] + f[i+stepI_][j][k+stepK_] + f[i][j+stepJ_][k+stepK_] + f[i+stepI_][j+stepJ_][k+stepK_])/4;
	   fl_1 = (f[i][j][k+stepK_1] + f[i+stepI_][j][k+stepK_1] + f[i][j+stepJ_][k+stepK_1] + f[i+stepI_][j+stepJ_][k+stepK_1])/4;
	   dfz = (fl - fl_1)/hzl;
	   
	   res = dfx*nx + dfy*ny + dfz*nz;

	   return res;
    }
    else
    {
       throw "Error in T3dLaplaceOperator::operator(): Operator is not implemented for this mask value.";
    }
   }
   else
   {
       throw "Error in T3dLaplaceOperator::operator(): Operator is not implemented for multidata index <> 0.";
   }
}
//*************************************

//#################################################################

//*************************************
T1dPackOperator::T1dPackOperator(const int lModelName, const TRnLinearOperator &lA, double** lVector):
                        TRnRelease1dLinearOperator(lModelName),
                        TRnLinearOperator(lModelName)						
  {	
	 TRnLinearOperator &op = const_cast<TRnLinearOperator &>(lA);	 
     TSpace &s = op.NewArgument();      
     
	 TRnSpace *spRn = dynamic_cast<TRnSpace *>(&s);
     if (spRn==NULL)
     {
         throw "Error in T1dPackOperator constructor: Could not cast new argument to TRnSpace* type.";
     }
     TRnSpace &sRn = *spRn;

	 int N = sRn.EndIndex;
     	 
	 if (N!=(fModel->EndIndex))
	 {
       throw "Error in T1dPackOperator constructor: Contradiction between size of model and size of constructing matrix.";
	 }


     //Now we have linear operator under Rn and vector from Rn. It is all we need for actions.
	 //Turn on time control.
     //clock_t start = clock();
	 //clock_t end = clock();
	 double diff = 0;



	 printf("Constructing operator activity tables: Starting...\n");
	 op.ConstructActivityTables();	 
	 printf("Constructing operator activity tables: OK.\n");
	 


     //Transpose matrix (Begin)          	 
	 int* columnsElementsCount = new int [N];
	 double** matrixPackC = new double*[N];
	 int** rows = new int*[N];     
     
	 printf("Generating transpose matrix (dim = %d): Starting...\n",N);
	 //start = clock();
	 #pragma omp parallel
	 {
     
     TRnLinearOperator &A = const_cast<TRnLinearOperator &>(op.Copy());
	 A.SetActivityTable(op.GetiActivityTable());	 
	 
	 
     TRnSpace &z = sRn.CreateInstance();
     TRnSpace &Az = sRn.CreateInstance();
	 z|=0;
	 Az|=0; 

	 double* colVector = new double[N];   
	 int* colIndex = new int[N];   


	 bool useOptimization = A.IsActivityTableConstructed();	 

	 if (useOptimization) 
	 {
		 z.ActualIndexesCount = 1;
		 z.AutomaticActualIndexesAssigning = false;
		 Az.AutomaticActualIndexesAssigning = false;
	 }

	 int j_1 = 0;
     	 
	 #pragma omp for 
     for (int j=1;j<N+1;j++)
     {
		 if ((j)%100000==0) {printf("Generating transpose matrix: %d of %d\n",j,N);}
         
		 z[j] = 1;
         if (useOptimization) {z.ActualIndexes[0]=j;}
		 
		 if (j_1!=0) {z[j_1]=0;}
		 j_1=j;
         		 

		 Az|=A(z); //matrix column j
		 
		 //Attention! If using actual indexes, then after |= only theese indexes can be used.
		 //Other indexes are form previous calculations and corresponding values are incorrect.

		 int elementsCount = 0;
		 int aic = N;
		 if (useOptimization) {aic = A.ActualIndexesCount;}
				 
         for (int loopVar=1;loopVar<aic+1;loopVar++)
         {
			 int i = loopVar;
			 if (useOptimization) {i = A.ActualIndexes[loopVar-1];} 
			 if (Az[i]!=0) 
			 {
				 colVector[elementsCount] = Az[i];
				 colIndex[elementsCount] = i-1;
				 elementsCount++;				 
			 }
         }
		

		 if (elementsCount==0) {throw "Error in T1dPackOperator constructor: Zero row was found while transpose matrix constructing.";}

		 columnsElementsCount[j-1] = elementsCount;
         matrixPackC[j-1] = new double[elementsCount];
         rows[j-1] = new int[elementsCount];
		 		
		 for (int i=0;i<elementsCount;i++)
         {			 
		     matrixPackC[j-1][i] = colVector[i];
		     rows[j-1][i] = colIndex[i]; 			
		 }		 
     }


	 delete[] colVector;// colVector = NULL;
	 delete[] colIndex;// colIndex = NULL;

	 delete &A;
	 delete &z;
     delete &Az;
	 
	 }//#pragma omp parallel

	 //end = clock();
	 //diff = (double)(end - start);
	 //diff = diff/CLOCKS_PER_SEC;
	 //printf("Generating transponse matrix: OK. (time = %lf s.)\n",diff);	 


	 fRowsCountT = N;
	 fRowsElementsCountT = columnsElementsCount;
	 fMatrixPackT = matrixPackC;
	 fColumnsT = rows; 
	 //Transpose matrix (End)          	 



	 delete &s;


	 //Original matrix (Begin)          	 
	 fRowsCount = N;
	 fRowsElementsCount = new int[N];
	 fMatrixPack = new double*[N];
	 fColumns = new int*[N];     	
	  
	 		 
	        
     
     
	 printf("Generating original matrix (dim = %d): Starting...\n",N);

	 //start = clock();	   	 
	 #pragma omp parallel	
	 {         
       bool useOptimization = op.IsActivityTableConjugateConstructed();	  	   
	   int** activityTable = op.GetiActivityTableConjugate();

       double* colVector = new double[N];   
	   int* colIndex = new int[N];   	   
	   
	   #pragma omp for 
       for (int i=0;i<N;i++)
	   {
         if ((i+1)%100000==0) {printf("Generating original matrix: %d of %d\n",i+1,N);}

		 //work with i+1 unit vector

		 int elementsCount = 0;
		 
		 int aic = N;
		 if (useOptimization) {aic = activityTable[i][0];}
				 
         for (int loopVar=0;loopVar<aic;loopVar++)
         {
			 int j = loopVar;
			 if (useOptimization) {j = (activityTable[i][loopVar+1])-1;} 

			 double value = 0; 
			 int count = fRowsElementsCountT[j];
			 for (int k=0;k<count;k++)
			 {
				 if (fColumnsT[j][k]==i) {value = fMatrixPackT[j][k]; break;}
			 }

			 if (value!=0) 
			 {
				 colVector[elementsCount] = value;
				 colIndex[elementsCount] = j;
				 elementsCount++;				 
			 }
         }		


		 if (elementsCount==0)
		 {
			 throw "Error in T1dPackOperator constructor: Zero row was found while original matrix constructing.";
		 }

         fRowsElementsCount[i] = elementsCount;
		 fMatrixPack[i] = new double[elementsCount];
         fColumns[i] = new int[elementsCount];
		 
		 for (int j=0;j<elementsCount;j++) 
		 {
            fMatrixPack[i][j] = colVector[j];
			fColumns[i][j] = colIndex[j];; 
		 }		 
	   } //by rows

	   delete[] colVector;// colVector = NULL;
	   delete[] colIndex;// colIndex = NULL;

	 }//#pragma omp parallel
	   
	 //end = clock();
	 //diff = (double)(end - start);
	 //diff = diff/CLOCKS_PER_SEC;
	 //printf("Generating original matrix: OK. (time = %lf s.)\n",diff);	 

	 fUseTransposeMatrix = false;

	 fRowsCountAct = fRowsCount;
	 fRowsElementsCountAct = fRowsElementsCount;
	 fMatrixPackAct = fMatrixPack;
	 fColumnsAct = fColumns;     		 	   
	 //Original matrix (End)          	 

	 fExternalMatrixPack = false;
	 fExternalMatrixPackT = false;


     //********************************************************************************************
     //follow - deprecated 

	 if (*lVector!=NULL)	 
	 {	
       throw "Error in T1dPackOperator constructor: Method is not ready for Gausse transformation.";
	   
	   double* colVector = new double[N];   
	   
	   //Transforming vetor (Begin)	 
	   double* newVector = new double[N];
	   for (int i=0;i<N;i++)
	   {
		 double value = 0;		 
		 for (int j=0;j<columnsElementsCount[i];j++)
		 {
			 value = value + (matrixPackC[i][j])*((*lVector)[rows[i][j]]);
		 }
		 newVector[i] = value;
	   }	 
	   delete *lVector;
	   *lVector = newVector;
	   //Transforming vetor (End)


	   //Generating ultimate matrix (Begin)	
	   double* tempVector = new double[N];

	   for (int i=0;i<N;i++)
	   {
         if (i%1000==0)
		 {		 
          printf("Generating ultimate matrix: ");
		  printf("%d",i);
		  printf("\n");
		 }

		 int elementsCount = 0;
		 for (int j=0;j<N;j++) {colVector[j]=0;}

		 for (int j=0;j<N;j++)
		 {
		    double value = 0;
			for (int k=0;k<N;k++) {tempVector[k]=0;}
			for (int k=0;k<columnsElementsCount[j];k++)
			{
				tempVector[rows[j][k]]=matrixPackC[j][k];
			}	
			for (int k=0;k<columnsElementsCount[i];k++)
			{
				int row = rows[i][k];
				value = value + (matrixPackC[i][k])*tempVector[row];
			}

			//We have (i,j) element
			if (value!=0)
			{
				elementsCount++;
				colVector[j] = value;
			}
		 } //by columns		


		 if (elementsCount==0)
		 {
			 throw "Error in T1dPackOperator constructor: Zero row was found.";
		 }

         fRowsElementsCount[i] = elementsCount;
		 fMatrixPack[i] = new double[elementsCount];
         fColumns[i] = new int[elementsCount];

		 int elementNumber = -1;
		 for (int j=0;j<N;j++) 
		 {
			 double value = colVector[j];
			 if (value!=0)
			 {
                elementNumber++;
			    fMatrixPack[i][elementNumber] = value;
				fColumns[i][elementNumber] = j;
			 }
		 }
	   } // by rows

	   delete[] tempVector;// tempVector = NULL;
	   //Generating ultimate matrix (End)

	   delete[] colVector;// colVector = NULL;

     }// if transforming is need 
  }
//*************************************
T1dPackOperator::T1dPackOperator(
		                          const int lModelName,
		                          int lRowsCount,  int* lRowsElementsCount,  double** lMatrixPack,  int** lColumns,
				                  int lRowsCountT, int* lRowsElementsCountT, double** lMatrixPackT, int** lColumnsT
			                    ):TRnRelease1dLinearOperator(lModelName), TRnLinearOperator(lModelName)						
                        
{
	fRowsCount = lRowsCount;
	fRowsElementsCount = lRowsElementsCount;
	fMatrixPack = lMatrixPack;
	fColumns = lColumns;

	fRowsCountT = lRowsCountT;
	fRowsElementsCountT = lRowsElementsCountT;
	fMatrixPackT = lMatrixPackT;
	fColumnsT = lColumnsT;

	fUseTransposeMatrix = false;

	fRowsCountAct = fRowsCount;
	fRowsElementsCountAct = fRowsElementsCount;
	fMatrixPackAct = fMatrixPack;
	fColumnsAct = fColumns;     		 	   
	 

	fExternalMatrixPack = true;
	fExternalMatrixPackT = true;

}
//*************************************
T1dPackOperator::~T1dPackOperator()
{	

	if (fExternalMatrixPack == false)
	{

	for (int i=0;i<fRowsCount;i++)
	{		 
	    delete[] fMatrixPack[i];
		delete[] fColumns[i];
	}

	delete[] fMatrixPack;// fMatrixPack = NULL;
	delete[] fColumns;// fColumns = NULL;
	delete[] fRowsElementsCount;// fRowsElementsCount = NULL;

	}//if (fExternalMatrixPack == false)


	if (fExternalMatrixPackT == false)
	{

	for (int i=0;i<fRowsCountT;i++)
	{		 
	    delete[] fMatrixPackT[i];
		delete[] fColumnsT[i];
	}

	delete[] fMatrixPack;// fMatrixPack = NULL;
	delete[] fColumns;// fColumns = NULL;
	delete[] fRowsElementsCount;// fRowsElementsCount = NULL;

	}//if (fExternalMatrixPackT == false)
}
//*************************************
 void T1dPackOperator::ComfortableCall() const
  {
    TRnRelease1dData *Rn1dData = dynamic_cast<TRnRelease1dData *>(fCastArgument->MultiData[0]);
    f = Rn1dData->Data;

	if (fUseTransposeMatrix == false)
	{
	   fRowsCountAct = fRowsCount;
	   fRowsElementsCountAct = fRowsElementsCount;
	   fMatrixPackAct = fMatrixPack;
	   fColumnsAct = fColumns;     		 
	}
	else
	{
	   fRowsCountAct = fRowsCountT;
	   fRowsElementsCountAct = fRowsElementsCountT;
	   fMatrixPackAct = fMatrixPackT;
	   fColumnsAct = fColumnsT;     		 
	}
  }
//*************************************
 double T1dPackOperator::operator()(int idData,int i) const
 {
   //Жертва-ускорение	   
   //if (idData==0)
   //{    
     double res = 0;
	 int count = fRowsElementsCountAct[i];
	 double* matr = fMatrixPackAct[i];
	 int* colms = fColumnsAct[i];
	 for (int j=0;j<count;j++)
	 {
        double mij = matr[j];
		int col = colms[j];
        res = res + (mij)*(f[col]);
	 }
     return res;     	      
   //}
   //else
   //{
       //throw "Error in T1dPackOperator::operator(): Operator is not implemented for multidata index <> 0.";
   //}
 }
//*************************************

double T3dVariableDensityOperator::operator()(int idData, int i, int j, int k) const
{
    if (idData == 0)
    {
        const double* hx = fGrid->GetSeparator1().Separation;
        const double* hy = fGrid->GetSeparator2().Separation;
        const double* hz = fGrid->GetSeparator3().Separation;

        const T3dNumberMask& mask = *fMask;
        const int mijk = (*fMask)[i][j][k];


        int stepI_;
        int stepI_1;
        double hxl;

        int stepJ_;
        int stepJ_1;
        double hyl;

        int stepK_;
        int stepK_1;
        double hzl;

        double fl;
        double fl_1;

        double dfx;
        double dfy;
        double dfz;

        if (mijk == TActualPoint)
        {
            double h_1X;
            double hX;
            double hx2X;
            double hx_2X;
            double hxhxhX;
            h_1X = hx[i - 1];
            hX = hx[i];
            hx2X = hx[i - 1] + hx[i];
            hx_2X = (hx[i - 1] + hx[i]) / 2;
            hxhxhX = hx[i - 1] * hx[i] * (hx[i - 1] + hx[i]) / 2;

            double h_1Y;
            double hY;
            double hx2Y;
            double hx_2Y;
            double hxhxhY;
            h_1Y = hy[j - 1];
            hY = hy[j];
            hx2Y = hy[j - 1] + hy[j];
            hx_2Y = (hy[j - 1] + hy[j]) / 2;
            hxhxhY = hy[j - 1] * hy[j] * (hy[j - 1] + hy[j]) / 2;

            double h_1Z;
            double hZ;
            double hx2Z;
            double hx_2Z;
            double hxhxhZ;
            h_1Z = hz[k - 1];
            hZ = hz[k];
            hx2Z = hz[k - 1] + hz[k];
            hx_2Z = (hz[k - 1] + hz[k]) / 2;
            hxhxhZ = hz[k - 1] * hz[k] * (hz[k - 1] + hz[k]) / 2;

            double f_i_1 = (f[i][j][k] - f[i - 1][j][k]) / (h_1X);
            double f_i = (f[i + 1][j][k] - f[i][j][k]) / (hX);
            double f_j_1 = (f[i][j][k] - f[i][j - 1][k]) / (h_1Y);
            double f_j = (f[i][j + 1][k] - f[i][j][k]) / (hY);
            double f_k_1 = (f[i][j][k] - f[i][j][k - 1]) / (h_1Z);
            double f_k = (f[i][j][k + 1] - f[i][j][k]) / (hZ);

            double d_i = (density[i+1][j][k] - density[i][j][k]) / (h_1X);
            double d_j = (density[i][j+1][k] - density[i][j][k]) / (h_1Y);
            double d_k = (density[i][j][k+1] - density[i][j][k]) / (h_1Z);
            //printf("di = density[i+1][j][k] - density[i][j][k]  (%LF %LF) / %LF \n", density[i+1][j][k], density[i][j][k], h_1X);

            //Pre-border points (Begin)
            if (mask[i - 1][j][k] == TFictivePoint)
            {
                f_i_1 = 0;
            }
            if (mask[i + 1][j][k] == TFictivePoint)
            {
                f_i = 0;
                d_i = 0;
            }
            if (mask[i][j - 1][k] == TFictivePoint)
            {
                f_j_1 = 0;
            }
            if (mask[i][j + 1][k] == TFictivePoint)
            {
                f_j = 0;
                d_j = 0;
            }
            if (mask[i][j][k - 1] == TFictivePoint)
            {
                f_k_1 = 0;
            }
            if (mask[i][j][k + 1] == TFictivePoint)
            {
                f_k = 0;
                d_k = 0;
            }
            //Pre-border points (End)

            double laplace_p = (
                (f_i - f_i_1) / hx_2X
                +
                (f_j - f_j_1) / hx_2Y
                +
                (f_k - f_k_1) / hx_2Z
            );

            double divergence_p = (
                f_i + f_j + f_k
            );

            double divergence_density = (
                d_i + d_j + d_k
            );

            //printf("density[%d][%d][%d] = %LF\t divergence_density = %LF\t divergence_p = %LF\n",
                    //i, j, k, density[i][j][k], divergence_density, divergence_p);
            return -(density[i][j][k] * laplace_p - divergence_density * divergence_p);
        }
        else if  (mijk == TPreDefinedBorderPoint)
        {
            const T3dVector* normal = fNormals[i][j][k];
            if (normal == NULL)
            {
                throw "Error in T3dLaplaceOperator::operator(): Normal vector is undefined for this point.";
            }
            const double nx = normal->X;
            const double ny = normal->Y;
            const double nz = normal->Z;


            if ( (nx == 0) && (ny == 0) && (nz == 0) )
            {
                throw "Error in T3dLaplaceOperator::operator(): Normal vector is zero for this point.";
            }
            else if (nx != 0)
            {
                if ( (ny != 0) || (nz != 0) || (fabs(nx) != 1) )
                {
                    throw "Error in T3dLaplaceOperator::operator(): Normal vector is incorrect for this point.";
                }
            }
            else if (ny != 0)
            {
                if ( (nx != 0) || (nz != 0) || (fabs(ny) != 1) )
                {
                    throw "Error in T3dLaplaceOperator::operator(): Normal vector is incorrect for this point.";
                }
            }
            else if (nz != 0)
            {
                if ( (nx != 0) || (ny != 0) || (fabs(nz) != 1) )
                {
                    throw "Error in T3dLaplaceOperator::operator(): Normal vector is incorrect for this point.";
                }
            }

            double semiSum;

            if (nx > 0)
            {
                semiSum = (f[i][j][k] + f[i - 1][j][k]) / 2;
            }
            else if (nx < 0)
            {
                semiSum = (f[i + 1][j][k] + f[i][j][k]) / 2;
            }

            if (ny > 0)
            {
                semiSum = (f[i][j][k] + f[i][j - 1][k]) / 2;
            }
            else if (ny < 0)
            {
                semiSum = (f[i][j + 1][k] + f[i][j][k]) / 2;
            }

            if (nz > 0)
            {
                semiSum = (f[i][j][k] + f[i][j][k - 1]) / 2;
            }
            else if (nz < 0)
            {
                semiSum = (f[i][j][k + 1] + f[i][j][k]) / 2;
            }

            return semiSum;
        }
        else if  (mijk == TPreNormalBorderPoint)
        {
            //throw "Error in T3dLaplaceOperator::operator(): Operator is not implemented for TPreNormalBorderPoint mask value.";

            const T3dVector* normal = fNormals[i][j][k];
            if (normal == NULL)
            {
                throw "Error in T3dLaplaceOperator::operator(): Normal vector is undefined for this point.";
            }
            const double nx = normal->X;
            const double ny = normal->Y;
            const double nz = normal->Z;

            if ( (nx == 0) && (ny == 0) && (nz == 0) )
            {
                throw "Error in T3dLaplaceOperator::operator(): Normal vector is zero for this point.";
            }


            stepI_ = 0;
            stepI_1 = 0;
            hxl = 1;
            if (nx > 0)
            {
                stepI_ = 0;
                stepI_1 = -1;
                hxl = hx[i - 1];
            }
            else if (nx < 0)
            {
                stepI_ = 1;
                stepI_1 = 0;
                hxl = hx[i];
            }


            stepJ_ = 0;
            stepJ_1 = 0;
            hyl = 1;
            if (ny > 0)
            {
                stepJ_ = 0;
                stepJ_1 = -1;
                hyl = hy[j - 1];
            }
            else if (ny < 0)
            {
                stepJ_ = 1;
                stepJ_1 = 0;
                hyl = hy[j];
            }

            stepK_ = 0;
            stepK_1 = 0;
            hzl = 1;
            if (nz > 0)
            {
                stepK_ = 0;
                stepK_1 = -1;
                hzl = hz[k - 1];
            }
            else if (nz < 0)
            {
                stepK_ = 1;
                stepK_1 = 0; ;
                hzl = hz[k];
            }



            fl   = (f[i + stepI_][j][k] + f[i + stepI_][j + stepJ_][k] + f[i + stepI_][j][k + stepK_] + f[i + stepI_][j + stepJ_][k + stepK_]) / 4;
            fl_1 = (f[i + stepI_1][j][k] + f[i + stepI_1][j + stepJ_][k] + f[i + stepI_1][j][k + stepK_] + f[i + stepI_1][j + stepJ_][k + stepK_]) / 4;
            dfx = (fl - fl_1) / hxl;

            fl   = (f[i][j + stepJ_][k] + f[i + stepI_][j + stepJ_][k] + f[i][j + stepJ_][k + stepK_] + f[i + stepI_][j + stepJ_][k + stepK_]) / 4;
            fl_1 = (f[i][j + stepJ_1][k] + f[i + stepI_][j + stepJ_1][k] + f[i][j + stepJ_1][k + stepK_] + f[i + stepI_][j + stepJ_1][k + stepK_]) / 4;
            dfy = (fl - fl_1) / hyl;

            fl   = (f[i][j][k + stepK_] + f[i + stepI_][j][k + stepK_] + f[i][j + stepJ_][k + stepK_] + f[i + stepI_][j + stepJ_][k + stepK_]) / 4;
            fl_1 = (f[i][j][k + stepK_1] + f[i + stepI_][j][k + stepK_1] + f[i][j + stepJ_][k + stepK_1] + f[i + stepI_][j + stepJ_][k + stepK_1]) / 4;
            dfz = (fl - fl_1) / hzl;

            // TODO: not used?
            //double d_x = density[i+1][j][k] - density[i][j][k] / (h_1X);
            //double d_y = density[i][j+1][k] - density[i][j][k] / (h_1Y);
            //double d_z = density[i][j][k+1] - density[i][j][k] / (h_1Z);

            double laplace_p = dfx * nx + dfy * ny + dfz * nz;
            // TODO: not used?
            //double divergence_p = dfx + dfy + dfz;
            //double divergence_density = d_x + d_y + d_z;

            // why without minus?
            //return density * laplace_p - divergence_density * divergence_p;
            return laplace_p;
        }
        else
        {
            throw "Error in T3dLaplaceOperator::operator(): Operator is not implemented for this mask value.";
        }
    }
    else
    {
        throw "Error in T3dLaplaceOperator::operator(): Operator is not implemented for multidata index <> 0.";
    }
}

void T3dVariableDensityOperator::withDensity(TRnRelease3dSpace &density)
{
    int N = fGrid->GetSeparator1().EndIndex;
    int L = fGrid->GetSeparator2().EndIndex;
    int M = fGrid->GetSeparator3().EndIndex;

    this->density = new double **[N + 1];
    for (int i = 0; i < N + 1; i++)
    {
        this->density[i] = new double *[L + 1];
        for (int j = 0; j < L + 1; j++)
        {
            this->density[i][j] = new double [M + 1];
            for (int k = 0; k < M + 1; k++)
            {
                this->density[i][j][k] = density[i][j][k];
            }
        }
    }
}

//***********************************************

double T3dVariableDivergentedOperator::operator()(int idData, int i, int j, int k) const
{
	if (idData == 0)
	{
		//printf("%d %d %d \n",i,j,k);
		//const double* hx = fGrid->GetSeparator1().Separation;
		//const double* hy = fGrid->GetSeparator2().Separation;
		//const double* hz = fGrid->GetSeparator3().Separation;

		//const T3dNumberMask& (*fMask) = *fMask;
		const int mijk = (*fMask)[i][j][k];


		int stepI_;
		int stepI_1;
		double hxl;

		int stepJ_;
		int stepJ_1;
		double hyl;

		int stepK_;
		int stepK_1;
		double hzl;

		double fl;
		double fl_1;

		double dfx;
		double dfy;
		double dfz;

		if (mijk == TActualPoint)
		{
			double h_1X;
			double hX;
			double hx2X;
			double hx_2X;
			double hxhxhX;
			h_1X = fGrid->GetSeparator1().Separation[i - 1];
			hX = fGrid->GetSeparator1().Separation[i];
			hx2X = fGrid->GetSeparator1().Separation[i - 1] + fGrid->GetSeparator1().Separation[i];
			hx_2X = (fGrid->GetSeparator1().Separation[i - 1] + fGrid->GetSeparator1().Separation[i]) / 2;
			hxhxhX = fGrid->GetSeparator1().Separation[i - 1] * fGrid->GetSeparator1().Separation[i] * (fGrid->GetSeparator1().Separation[i - 1] + fGrid->GetSeparator1().Separation[i]) / 2;

			double h_1Y;
			double hY;
			double hx2Y;
			double hx_2Y;
			double hxhxhY;
			h_1Y = fGrid->GetSeparator2().Separation[j - 1];
			hY = fGrid->GetSeparator2().Separation[j];
			hx2Y = fGrid->GetSeparator2().Separation[j - 1] + fGrid->GetSeparator2().Separation[j];
			hx_2Y = (fGrid->GetSeparator2().Separation[j - 1] + fGrid->GetSeparator2().Separation[j]) / 2;
			hxhxhY = fGrid->GetSeparator2().Separation[j - 1] * fGrid->GetSeparator2().Separation[j] * (fGrid->GetSeparator2().Separation[j - 1] + fGrid->GetSeparator2().Separation[j]) / 2;

			double h_1Z;
			double hZ;
			double hx2Z;
			double hx_2Z;
			double hxhxhZ;
			h_1Z = fGrid->GetSeparator3().Separation[k - 1];
			hZ = fGrid->GetSeparator3().Separation[k];
			hx2Z = fGrid->GetSeparator3().Separation[k - 1] + fGrid->GetSeparator3().Separation[k];
			hx_2Z = (fGrid->GetSeparator3().Separation[k - 1] + fGrid->GetSeparator3().Separation[k]) / 2;
			hxhxhZ = fGrid->GetSeparator3().Separation[k - 1] * fGrid->GetSeparator3().Separation[k] * (fGrid->GetSeparator3().Separation[k - 1] + fGrid->GetSeparator3().Separation[k]) / 2;

			double f_i_1 = density[i - 1][j][k] * (f[i][j][k] - f[i - 1][j][k]) / (h_1X);
			double f_i   = density[i + 1][j][k] * (f[i + 1][j][k] - f[i][j][k]) / (hX);
			double f_j_1 = density[i][j - 1][k] * (f[i][j][k] - f[i][j - 1][k]) / (h_1Y);
			double f_j   = density[i][j + 1][k] * (f[i][j + 1][k] - f[i][j][k]) / (hY);
			double f_k_1 = density[i][j][k - 1] * (f[i][j][k] - f[i][j][k - 1]) / (h_1Z);
			double f_k   = density[i][j][k + 1] * (f[i][j][k + 1] - f[i][j][k]) / (hZ);

			double d_i = (density[i + 1][j][k] - density[i][j][k]) / (h_1X);
			double d_j = (density[i][j + 1][k] - density[i][j][k]) / (h_1Y);
			double d_k = (density[i][j][k + 1] - density[i][j][k]) / (h_1Z);
			//printf("di = density[i+1][j][k] - density[i][j][k]  (%LF %LF) / %LF \n", density[i+1][j][k], density[i][j][k], h_1X);

			//Pre-border points (Begin)
			if ((*fMask)[i - 1][j][k] == TFictivePoint)
			{
				f_i_1 = 0;
			}
			if ((*fMask)[i + 1][j][k] == TFictivePoint)
			{
				f_i = 0;
				d_i = 0;
			}
			if ((*fMask)[i][j - 1][k] == TFictivePoint)
			{
				f_j_1 = 0;
			}
			if ((*fMask)[i][j + 1][k] == TFictivePoint)
			{
				f_j = 0;
				d_j = 0;
			}
			if ((*fMask)[i][j][k - 1] == TFictivePoint)
			{
				f_k_1 = 0;
			}
			if ((*fMask)[i][j][k + 1] == TFictivePoint)
			{
				f_k = 0;
				d_k = 0;
			}
			//Pre-border points (End)

            double laplace_p = (
                (f_i - f_i_1) / hx_2X
                +
                (f_j - f_j_1) / hx_2Y
                +
                (f_k - f_k_1) / hx_2Z
            );

			double divergence_p = f_i + f_j + f_k;

			double divergence_density = d_i + d_j + d_k;

			//printf("density[%d][%d][%d] = %LF\t divergence_density = %LF\t divergence_p = %LF\n",
			//i, j, k, density[i][j][k], divergence_density, divergence_p);
			//return -(laplace_p - divergence_density * divergence_p);
			return -(density[i][j][k]*laplace_p - divergence_density/**divergence_p);//*/*f[i][j][k]);
		}
		else if (mijk == TPreDefinedBorderPoint)
		{
			//const T3dVector* normal = fNormals[i][j][k];
			if (fNormals[i][j][k] == NULL)
			{
				throw "Error in T3dLaplaceOperator::operator(): Normal vector is undefined for this point.";
			}
			const double nx = fNormals[i][j][k]->X;
			const double ny = fNormals[i][j][k]->Y;
			const double nz = fNormals[i][j][k]->Z;


			if ((nx == 0) && (ny == 0) && (nz == 0))
			{
				throw "Error in T3dLaplaceOperator::operator(): Normal vector is zero for this point.";
			}
			else if (nx != 0)
			{
				if ((ny != 0) || (nz != 0) || (fabs(nx) != 1))
				{
					throw "Error in T3dLaplaceOperator::operator(): Normal vector is incorrect for this point.";
				}
			}
			else if (ny != 0)
			{
				if ((nx != 0) || (nz != 0) || (fabs(ny) != 1))
				{
					throw "Error in T3dLaplaceOperator::operator(): Normal vector is incorrect for this point.";
				}
			}
			else if (nz != 0)
			{
				if ((nx != 0) || (ny != 0) || (fabs(nz) != 1))
				{
					throw "Error in T3dLaplaceOperator::operator(): Normal vector is incorrect for this point.";
				}
			}

			double semiSum;

			if (nx > 0)
			{
				semiSum = (f[i][j][k] + f[i - 1][j][k]) / 2;
			}
			else if (nx < 0)
			{
				semiSum = (f[i + 1][j][k] + f[i][j][k]) / 2;
			}

			if (ny > 0)
			{
				semiSum = (f[i][j][k] + f[i][j - 1][k]) / 2;
			}
			else if (ny < 0)
			{
				semiSum = (f[i][j + 1][k] + f[i][j][k]) / 2;
			}

			if (nz > 0)
			{
				semiSum = (f[i][j][k] + f[i][j][k - 1]) / 2;
			}
			else if (nz < 0)
			{
				semiSum = (f[i][j][k + 1] + f[i][j][k]) / 2;
			}
			
			return semiSum;
		}
		else if (mijk == TPreNormalBorderPoint)
		{
			//throw "Error in T3dLaplaceOperator::operator(): Operator is not implemented for TPreNormalBorderPoint mask value.";
			double h_1X;
			double hX;
			double hx2X;
			double hx_2X;
			double hxhxhX;
			h_1X = fGrid->GetSeparator1().Separation[i - 1];
			hX = fGrid->GetSeparator1().Separation[i];
			hx2X = fGrid->GetSeparator1().Separation[i - 1] + fGrid->GetSeparator1().Separation[i];
			hx_2X = (fGrid->GetSeparator1().Separation[i - 1] + fGrid->GetSeparator1().Separation[i]) / 2;
			hxhxhX = fGrid->GetSeparator1().Separation[i - 1] * fGrid->GetSeparator1().Separation[i] * (fGrid->GetSeparator1().Separation[i - 1] + fGrid->GetSeparator1().Separation[i]) / 2;

			double h_1Y;
			double hY;
			double hx2Y;
			double hx_2Y;
			double hxhxhY;
			h_1Y = fGrid->GetSeparator2().Separation[j - 1];
			hY = fGrid->GetSeparator2().Separation[j];
			hx2Y = fGrid->GetSeparator2().Separation[j - 1] + fGrid->GetSeparator2().Separation[j];
			hx_2Y = (fGrid->GetSeparator2().Separation[j - 1] + fGrid->GetSeparator2().Separation[j]) / 2;
			hxhxhY = fGrid->GetSeparator2().Separation[j - 1] * fGrid->GetSeparator2().Separation[j] * (fGrid->GetSeparator2().Separation[j - 1] + fGrid->GetSeparator2().Separation[j]) / 2;

			double h_1Z;
			double hZ;
			double hx2Z;
			double hx_2Z;
			double hxhxhZ;
			h_1Z = fGrid->GetSeparator3().Separation[k - 1];
			hZ = fGrid->GetSeparator3().Separation[k];
			hx2Z = fGrid->GetSeparator3().Separation[k - 1] + fGrid->GetSeparator3().Separation[k];
			hx_2Z = (fGrid->GetSeparator3().Separation[k - 1] + fGrid->GetSeparator3().Separation[k]) / 2;
			hxhxhZ = fGrid->GetSeparator3().Separation[k - 1] * fGrid->GetSeparator3().Separation[k] * (fGrid->GetSeparator3().Separation[k - 1] + fGrid->GetSeparator3().Separation[k]) / 2;
			
			//const T3dVector* normal = fNormals[i][j][k];
			if (fNormals[i][j][k] == NULL)
			{
				throw "Error in T3dLaplaceOperator::operator(): Normal vector is undefined for this point.";
			}
			const double nx = fNormals[i][j][k]->X;
			const double ny = fNormals[i][j][k]->Y;
			const double nz = fNormals[i][j][k]->Z;

			if ((nx == 0) && (ny == 0) && (nz == 0))
			{
				throw "Error in T3dLaplaceOperator::operator(): Normal vector is zero for this point.";
			}


			stepI_ = 0;
			stepI_1 = 0;
			hxl = 1;
			if (nx > 0)
			{
				stepI_ = 0;
				stepI_1 = -1;
				hxl = fGrid->GetSeparator1().Separation[i - 1];
			}
			else if (nx < 0)
			{
				stepI_ = 1;
				stepI_1 = 0;
				hxl = fGrid->GetSeparator1().Separation[i];
			}


			stepJ_ = 0;
			stepJ_1 = 0;
			hyl = 1;
			if (ny > 0)
			{
				stepJ_ = 0;
				stepJ_1 = -1;
				hyl = fGrid->GetSeparator2().Separation[j - 1];
			}
			else if (ny < 0)
			{
				stepJ_ = 1;
				stepJ_1 = 0;
				hyl = fGrid->GetSeparator2().Separation[j];
			}

			stepK_ = 0;
			stepK_1 = 0;
			hzl = 1;
			if (nz > 0)
			{
				stepK_ = 0;
				stepK_1 = -1;
				hzl = fGrid->GetSeparator3().Separation[k - 1];
			}
			else if (nz < 0)
			{
				stepK_ = 1;
				stepK_1 = 0; ;
				hzl = fGrid->GetSeparator3().Separation[k];
			}



			fl = (f[i + stepI_][j][k] + f[i + stepI_][j + stepJ_][k] + f[i + stepI_][j][k + stepK_] + f[i + stepI_][j + stepJ_][k + stepK_]) / 4;
			fl_1 = (f[i + stepI_1][j][k] + f[i + stepI_1][j + stepJ_][k] + f[i + stepI_1][j][k + stepK_] + f[i + stepI_1][j + stepJ_][k + stepK_]) / 4;
			dfx = (fl - fl_1) / hxl;

			fl = (f[i][j + stepJ_][k] + f[i + stepI_][j + stepJ_][k] + f[i][j + stepJ_][k + stepK_] + f[i + stepI_][j + stepJ_][k + stepK_]) / 4;
			fl_1 = (f[i][j + stepJ_1][k] + f[i + stepI_][j + stepJ_1][k] + f[i][j + stepJ_1][k + stepK_] + f[i + stepI_][j + stepJ_1][k + stepK_]) / 4;
			dfy = (fl - fl_1) / hyl;

			fl = (f[i][j][k + stepK_] + f[i + stepI_][j][k + stepK_] + f[i][j + stepJ_][k + stepK_] + f[i + stepI_][j + stepJ_][k + stepK_]) / 4;
			fl_1 = (f[i][j][k + stepK_1] + f[i + stepI_][j][k + stepK_1] + f[i][j + stepJ_][k + stepK_1] + f[i + stepI_][j + stepJ_][k + stepK_1]) / 4;
			dfz = (fl - fl_1) / hzl;

			// TODO: not used?
			double d_x = density[i+1][j][k] - density[i][j][k] / (h_1X);
			double d_y = density[i][j+1][k] - density[i][j][k] / (h_1Y);
			double d_z = density[i][j][k+1] - density[i][j][k] / (h_1Z);

			double laplace_p = dfx * nx + dfy * ny + dfz * nz;
			// TODO: not used?
			//double divergence_p = dfx + dfy + dfz;
			double divergence_density = d_x + d_y + d_z;

			// why without minus?
			return density[i][j][k] * laplace_p - divergence_density * f[i][j][k];
			//return laplace_p;
		}
		else
		{
			throw "Error in T3dLaplaceOperator::operator(): Operator is not implemented for this mask value.";
		}
	}
	else
	{
		throw "Error in T3dLaplaceOperator::operator(): Operator is not implemented for multidata index <> 0.";
	}
}

void T3dVariableDivergentedOperator::withDensity(TRnRelease3dSpace &density)
{
	int N = fGrid->GetSeparator1().EndIndex;
	int L = fGrid->GetSeparator2().EndIndex;
	int M = fGrid->GetSeparator3().EndIndex;

	//if(this->density == NULL)
	//{
		this->density = new double **[N + 1];
		for (int i = 0; i < N + 1; i++)
		{
			this->density[i] = new double *[L + 1];
			for (int j = 0; j < L + 1; j++)
			{
				this->density[i][j] = new double[M + 1];
				for (int k = 0; k < M + 1; k++)
				{
					this->density[i][j][k] = density[i][j][k];
				}
			}
		}
	/*}
	else
		for (int i = 0; i < N + 1; i++)
		{
			for (int j = 0; j < L + 1; j++)
			{
				for (int k = 0; k < M + 1; k++)
				{
					this->density[i][j][k] = density[i][j][k];
				}
			}
		}*/
}

T3dVariableDivergentedOperator::~T3dVariableDivergentedOperator()
{
	int N = fGrid->GetSeparator1().EndIndex;
	int L = fGrid->GetSeparator2().EndIndex;
	//int M = fGrid->GetSeparator3().EndIndex;
//printf("Deleting density in operator\n");
	//if(density != NULL)
	//{

		for (int i = 0; i < N + 1; i++)
		{
			for (int j = 0; j < L + 1; j++)
			{
				delete[] this->density[i][j];			
			}
			delete[] this->density[i];
		}

		delete[] this->density;
		//delete density;
	//	density = NULL;
		//printf("Done\n");
	//}
}

//***********************************************

double T3dVariableVsStreamOperator::operator()(int idData, TRnRelease3dSpace* U, TRnRelease3dSpace* V, TRnRelease3dSpace* W, int i, int j, int k) const
{
	if (idData == 0)
	{
		const double* hx = fGrid->GetSeparator1().Separation;
		const double* hy = fGrid->GetSeparator2().Separation;
		const double* hz = fGrid->GetSeparator3().Separation;

		const T3dNumberMask& mask = *fMask;
		const int mijk = (*fMask)[i][j][k];


		int stepI_;
		int stepI_1;
		double hxl;

		int stepJ_;
		int stepJ_1;
		double hyl;

		int stepK_;
		int stepK_1;
		double hzl;

		double fl;
		double fl_1;

		double dfx;
		double dfy;
		double dfz;

		if (mijk == TActualPoint)
		{
			double h_1X;
			double hX;
			double hx2X;
			double hx_2X;
			double hxhxhX;
			h_1X = hx[i - 1];
			hX = hx[i];
			hx2X = hx[i - 1] + hx[i];
			hx_2X = (hx[i - 1] + hx[i]) / 2;
			hxhxhX = hx[i - 1] * hx[i] * (hx[i - 1] + hx[i]) / 2;

			double h_1Y;
			double hY;
			double hx2Y;
			double hx_2Y;
			double hxhxhY;
			h_1Y = hy[j - 1];
			hY = hy[j];
			hx2Y = hy[j - 1] + hy[j];
			hx_2Y = (hy[j - 1] + hy[j]) / 2;
			hxhxhY = hy[j - 1] * hy[j] * (hy[j - 1] + hy[j]) / 2;

			double h_1Z;
			double hZ;
			double hx2Z;
			double hx_2Z;
			double hxhxhZ;
			double f_i_1;
			double f_i;
			double f_j_1;
			double f_j;
			double f_k_1;
			double f_k; 
			h_1Z = hz[k - 1];
			hZ = hz[k];
			hx2Z = hz[k - 1] + hz[k];
			hx_2Z = (hz[k - 1] + hz[k]) / 2;
			hxhxhZ = hz[k - 1] * hz[k] * (hz[k - 1] + hz[k]) / 2;
			printf("Debug!\n");
			if(U[i][j][k] > 0)
				f_i_1 = density[i - 1][j][k] * (f[i][j][k] - f[i - 1][j][k]) / (h_1X);
			else
				f_i_1 = density[i - 1][j][k] * (f[i - 1][j][k] - f[i][j][k]) / (h_1X);
			
			if(V[i][j][k] > 0)
				f_j_1 = density[i][j - 1][k] * (f[i][j][k] - f[i][j - 1][k]) / (h_1Y);
			else
				f_i   = density[i + 1][j][k] * (f[i + 1][j][k] - f[i][j][k]) / (hX);
			
			if(W[i][j][k] > 0)
				f_k_1 = density[i][j][k - 1] * (f[i][j][k] - f[i][j][k - 1]) / (h_1Z);
			else
				f_j   = density[i][j + 1][k] * (f[i][j + 1][k] - f[i][j][k]) / (hY);
			
			f_k   = density[i][j][k + 1] * (f[i][j][k + 1] - f[i][j][k]) / (hZ);

			double d_i = (density[i + 1][j][k] - density[i][j][k]) / (h_1X);
			double d_j = (density[i][j + 1][k] - density[i][j][k]) / (h_1Y);
			double d_k = (density[i][j][k + 1] - density[i][j][k]) / (h_1Z);
			//printf("di = density[i+1][j][k] - density[i][j][k]  (%LF %LF) / %LF \n", density[i+1][j][k], density[i][j][k], h_1X);

			//Pre-border points (Begin)
			if (mask[i - 1][j][k] == TFictivePoint)
			{
				f_i_1 = 0;
			}
			if (mask[i + 1][j][k] == TFictivePoint)
			{
				f_i = 0;
				d_i = 0;
			}
			if (mask[i][j - 1][k] == TFictivePoint)
			{
				f_j_1 = 0;
			}
			if (mask[i][j + 1][k] == TFictivePoint)
			{
				f_j = 0;
				d_j = 0;
			}
			if (mask[i][j][k - 1] == TFictivePoint)
			{
				f_k_1 = 0;
			}
			if (mask[i][j][k + 1] == TFictivePoint)
			{
				f_k = 0;
				d_k = 0;
			}
			//Pre-border points (End)
			double dd_dx, dd_dy, dd_dz;
			
			if(U[i][j][k] > 0)
				dd_dx = (f_i - f_i_1) / hx_2X;
			else
				dd_dx = (f_i_1 - f_i) / hx_2X;
			
			if(V[i][j][k] > 0)
				dd_dy = (f_j - f_j_1) / hx_2Y;
			else
				dd_dy = (f_j_1 - f_j) / hx_2Y;
			
			if(W[i][j][k] > 0)
				dd_dz = (f_k - f_k_1) / hx_2Z;
			else
				dd_dz = (f_k_1 - f_k) / hx_2Z;
			
			double laplace_p = ( dd_dx + dd_dy + dd_dy );

			double divergence_p = (
				f_i + f_j + f_k
				);

			double divergence_density = (
				d_i + d_j + d_k
				);

			//printf("density[%d][%d][%d] = %LF\t divergence_density = %LF\t divergence_p = %LF\n",
			//i, j, k, density[i][j][k], divergence_density, divergence_p);
			return -(laplace_p - divergence_density * divergence_p);
		}
		else if (mijk == TPreDefinedBorderPoint)
		{
			const T3dVector* normal = fNormals[i][j][k];
			if (normal == NULL)
			{
				throw "Error in T3dLaplaceOperator::operator(): Normal vector is undefined for this point.";
			}
			const double nx = normal->X;
			const double ny = normal->Y;
			const double nz = normal->Z;


			if ((nx == 0) && (ny == 0) && (nz == 0))
			{
				throw "Error in T3dLaplaceOperator::operator(): Normal vector is zero for this point.";
			}
			else if (nx != 0)
			{
				if ((ny != 0) || (nz != 0) || (fabs(nx) != 1))
				{
					throw "Error in T3dLaplaceOperator::operator(): Normal vector is incorrect for this point.";
				}
			}
			else if (ny != 0)
			{
				if ((nx != 0) || (nz != 0) || (fabs(ny) != 1))
				{
					throw "Error in T3dLaplaceOperator::operator(): Normal vector is incorrect for this point.";
				}
			}
			else if (nz != 0)
			{
				if ((nx != 0) || (ny != 0) || (fabs(nz) != 1))
				{
					throw "Error in T3dLaplaceOperator::operator(): Normal vector is incorrect for this point.";
				}
			}

			double semiSum;

			if (nx > 0)
			{
				semiSum = (f[i][j][k] + f[i - 1][j][k]) / 2;
			}
			else if (nx < 0)
			{
				semiSum = (f[i + 1][j][k] + f[i][j][k]) / 2;
			}

			if (ny > 0)
			{
				semiSum = (f[i][j][k] + f[i][j - 1][k]) / 2;
			}
			else if (ny < 0)
			{
				semiSum = (f[i][j + 1][k] + f[i][j][k]) / 2;
			}

			if (nz > 0)
			{
				semiSum = (f[i][j][k] + f[i][j][k - 1]) / 2;
			}
			else if (nz < 0)
			{
				semiSum = (f[i][j][k + 1] + f[i][j][k]) / 2;
			}

			return semiSum;
		}
		else if (mijk == TPreNormalBorderPoint)
		{
			//throw "Error in T3dLaplaceOperator::operator(): Operator is not implemented for TPreNormalBorderPoint mask value.";

			const T3dVector* normal = fNormals[i][j][k];
			if (normal == NULL)
			{
				throw "Error in T3dLaplaceOperator::operator(): Normal vector is undefined for this point.";
			}
			const double nx = normal->X;
			const double ny = normal->Y;
			const double nz = normal->Z;

			if ((nx == 0) && (ny == 0) && (nz == 0))
			{
				throw "Error in T3dLaplaceOperator::operator(): Normal vector is zero for this point.";
			}


			stepI_ = 0;
			stepI_1 = 0;
			hxl = 1;
			if (nx > 0)
			{
				stepI_ = 0;
				stepI_1 = -1;
				hxl = hx[i - 1];
			}
			else if (nx < 0)
			{
				stepI_ = 1;
				stepI_1 = 0;
				hxl = hx[i];
			}


			stepJ_ = 0;
			stepJ_1 = 0;
			hyl = 1;
			if (ny > 0)
			{
				stepJ_ = 0;
				stepJ_1 = -1;
				hyl = hy[j - 1];
			}
			else if (ny < 0)
			{
				stepJ_ = 1;
				stepJ_1 = 0;
				hyl = hy[j];
			}

			stepK_ = 0;
			stepK_1 = 0;
			hzl = 1;
			if (nz > 0)
			{
				stepK_ = 0;
				stepK_1 = -1;
				hzl = hz[k - 1];
			}
			else if (nz < 0)
			{
				stepK_ = 1;
				stepK_1 = 0; ;
				hzl = hz[k];
			}



			fl = (f[i + stepI_][j][k] + f[i + stepI_][j + stepJ_][k] + f[i + stepI_][j][k + stepK_] + f[i + stepI_][j + stepJ_][k + stepK_]) / 4;
			fl_1 = (f[i + stepI_1][j][k] + f[i + stepI_1][j + stepJ_][k] + f[i + stepI_1][j][k + stepK_] + f[i + stepI_1][j + stepJ_][k + stepK_]) / 4;
			dfx = (fl - fl_1) / hxl;

			fl = (f[i][j + stepJ_][k] + f[i + stepI_][j + stepJ_][k] + f[i][j + stepJ_][k + stepK_] + f[i + stepI_][j + stepJ_][k + stepK_]) / 4;
			fl_1 = (f[i][j + stepJ_1][k] + f[i + stepI_][j + stepJ_1][k] + f[i][j + stepJ_1][k + stepK_] + f[i + stepI_][j + stepJ_1][k + stepK_]) / 4;
			dfy = (fl - fl_1) / hyl;

			fl = (f[i][j][k + stepK_] + f[i + stepI_][j][k + stepK_] + f[i][j + stepJ_][k + stepK_] + f[i + stepI_][j + stepJ_][k + stepK_]) / 4;
			fl_1 = (f[i][j][k + stepK_1] + f[i + stepI_][j][k + stepK_1] + f[i][j + stepJ_][k + stepK_1] + f[i + stepI_][j + stepJ_][k + stepK_1]) / 4;
			dfz = (fl - fl_1) / hzl;

			// TODO: not used?
			//double d_x = density[i+1][j][k] - density[i][j][k] / (h_1X);
			//double d_y = density[i][j+1][k] - density[i][j][k] / (h_1Y);
			//double d_z = density[i][j][k+1] - density[i][j][k] / (h_1Z);

			double laplace_p = dfx * nx + dfy * ny + dfz * nz;
			// TODO: not used?
			//double divergence_p = dfx + dfy + dfz;
			//double divergence_density = d_x + d_y + d_z;

			// why without minus?
			//return density * laplace_p - divergence_density * divergence_p;
			return laplace_p;
		}
		else
		{
			throw "Error in T3dLaplaceOperator::operator(): Operator is not implemented for this mask value.";
		}
	}
	else
	{
		throw "Error in T3dLaplaceOperator::operator(): Operator is not implemented for multidata index <> 0.";
	}
}

void T3dVariableVsStreamOperator::withDensity(TRnRelease3dSpace &density)
{
	int N = fGrid->GetSeparator1().EndIndex;
	int L = fGrid->GetSeparator2().EndIndex;
	int M = fGrid->GetSeparator3().EndIndex;

	this->density = new double **[N + 1];
	for (int i = 0; i < N + 1; i++)
	{
		this->density[i] = new double *[L + 1];
		for (int j = 0; j < L + 1; j++)
		{
			this->density[i][j] = new double[M + 1];
			for (int k = 0; k < M + 1; k++)
			{
				this->density[i][j][k] = density[i][j][k];
			}
		}
	}
}

//#################################################################
