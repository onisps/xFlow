#pragma hdrstop

#include "XSpace.h"


//---------------------Abstract (Begin)------------------------------


//#################################################################

//*************************************
   bool TSpace::deleteIntermediateItems=true; //По умолчанию Удалять промежуточные объекты
   int  TSpace::CountInstances=0;
   bool  TSpace::ShowRemainCountInstances=false;
//*************************************
   TSpace::TSpace() //Конструктор
    {
      intermediate=false;
      CountInstances=CountInstances+1;
      //ShowMessage("Count objects is "+IntToStr(CountInstances));
    }
//*************************************
   TSpace &TSpace::operator|=(const TSpace &s)
    {
        if (this!=&s)
         {
           ReleaseOperatorEqual(s,false);
           DeleteIntermediate(s);
         }
        return *this;
    }
//*************************************
   TSpace &TSpace::Assign(const TSpace &s)
    {
        if (this!=&s)
         {
           ReleaseOperatorEqual(s,true);
           DeleteIntermediate(s);
         }
        return *this;
    }
//*************************************
   TSpace::TSpace(const TSpace &) //Запрещенный конструктор копии
    {

    }
//*************************************
   TSpace &TSpace::operator=(const TSpace &s) //Запрещенная операция =
    {
           return *this;
    }
//*************************************
   TSpace::~TSpace() //Деструктор
    {
      CountInstances=CountInstances-1;
      if (ShowRemainCountInstances==true)
        {
          //ShowMessage("Remain count objects is "+IntToStr(CountInstances));
        }
    }
//*************************************
   void TSpace::DeleteIntermediate(const TSpace &s) //Удаление промежуточных объектов
    {
      if ((s.intermediate==true) && (deleteIntermediateItems==true)) delete &s;
    }
//*************************************
   bool TSpace::IsIntermediate()
    {
     return intermediate;
    }
   void TSpace::SetAsPermanent() //Закрепление объекта
    {
       intermediate=false;
    }
//*************************************
   void TSpace::SetAsIntermediate() //Пометка объекта как промежуточного
    {
      intermediate=true;
    }
//*************************************
TSpace &TSpace::CreateInstance() 
{
   TSpace &res=NewItem();
   return res;
}
//*************************************
   TSpace &TSpace::Copy()   //Генерации копии
    {
     TSpace &res=NewItem();
     
     //res|=*this;
     res.ReleaseOperatorEqual(*this,true);
     DeleteIntermediate(*this);

     return res;
    }
//*************************************

//#################################################################

//*************************************
  TMetricSpace::TMetricSpace():TSpace() //Конструктор
   {
   }
//*************************************
  TMetricSpace::~TMetricSpace() //Деструктор
   {
   }
//*************************************
  double TMetricSpace::Metric(const TSpace &s, int type)
   {
     double m = ReleaseMetric(s,type);
     if (&s!=this)
      {
       DeleteIntermediate(s);
      }
     DeleteIntermediate(*this);
     return m;
   }
//*************************************   

//#################################################################


  TFullSpace::TFullSpace():TMetricSpace() //Конструктор
   {
   }
//*************************************
  TFullSpace::~TFullSpace() //Деструктор
   {
   }
//*************************************

//#################################################################


//---------------------Abstract (End)------------------------------





//---------------------Linear (Begin)------------------------------


//#################################################################

//*************************************
  TLinearSpace::TLinearSpace():TSpace() //Конструктор
   {
   }
//*************************************
  TLinearSpace::~TLinearSpace() //Деструктор
   {
   }
//*************************************
  TSpace &TLinearSpace::operator|=(const int zero)
    {
         if (zero!=0)
           {
            throw "Error in TLinearSpace::operator|=: The parametr must be equal zero.";
           }
         ToZero();
         return *this;
    }
//*************************************
  TSpace &TLinearSpace::operator|=(const TLinearCombinationWith2Operands &lc)
    {
        ReleaseAssignSpeedLinearCombination(lc.p1,lc.s1,lc.p2,lc.s2,lc.operation);
        if ( (this!=lc.s1) && (lc.s1!=NULL)  )
         {
           DeleteIntermediate(*lc.s1);
         }
         if ( (this!=lc.s2) && (lc.s2!=NULL) )
         {
           DeleteIntermediate(*lc.s2);
         }
         delete &lc;
         return *this;
    }
//*************************************
  TSpace &TLinearSpace::operator|=(const TLinearCombinationWith2LinearOperatorOperands &lc)
    {
        ReleaseAssignSpeedLinearCombination(lc.A1,lc.s2,lc.operation);
        TSpace& arg = lc.A1->Argument;
        TSpace *s1 = &arg;
        if ( (this!=s1) && (s1!=NULL)  )
         {
           DeleteIntermediate(*s1);
           ( const_cast<TLinearOperator *>(lc.A1) )->fArgument = NULL;
         }
         if ( (this!=lc.s2) && (lc.s2!=NULL) )
         {
           DeleteIntermediate(*lc.s2);
         }
         delete &lc;
         return *this;
    }
//*************************************

//#################################################################

  TLinearCombinationWith2Operands &TLinearSpace::operator+(TLinearCombinationWith2Operands &lc)
    {
       lc.AddSecondOperand(*this,true);
       return lc;
    }
//*************************************
  TLinearCombinationWith2Operands &TLinearSpace::operator-(TLinearCombinationWith2Operands &lc)
    {
        lc.InsertSecondOperand(*this,false);
        return lc;
    }
//*************************************
  TLinearCombinationWith2Operands &TLinearSpace::operator+(const TSpace &s)
    {
       return *(new TLinearCombinationWith2Operands(*this,s,true));
    }
//*************************************
  TLinearCombinationWith2Operands &TLinearSpace::operator-(const TSpace &s)
    {
       return *(new TLinearCombinationWith2Operands(*this,s,false));
    }
//*************************************
   TLinearCombinationWith2Operands &TLinearSpace::operator*(double p)
    {
       return *(new TLinearCombinationWith2Operands(p,*this));
    }
//*************************************


  //Вспомогательные функции и классы

      //Класс линейная комбинация с двумя операндами
        //*************************************
          TLinearCombinationWith2Operands::TLinearCombinationWith2Operands(double lp1,
                                                                             const TSpace &ls1)
              {
                   p1=lp1;
                   s1=&ls1;

                   p2=0;
                   s2=NULL;

                   operation=true;
              }
        //*************************************
           TLinearCombinationWith2Operands::TLinearCombinationWith2Operands(const TSpace &ls1,
                                                                            const TSpace &ls2,
                                                                            bool loperation)
              {

                  p1=1;
                  s1=&ls1;

                  p2=1;
                  s2=&ls2;

                  operation=loperation;

              }
        //*************************************
           TLinearCombinationWith2Operands::~TLinearCombinationWith2Operands()
              {
                 //Ничего дополнительного делать не надо....
              }
        //*************************************
           const TSpace *TLinearCombinationWith2Operands::GetFirstOperand() const
            {
               return s1;
            }
        //*************************************
           const TSpace *TLinearCombinationWith2Operands::GetSecondOperand() const
            {
               return s2;
            }
        //*************************************
           void TLinearCombinationWith2Operands::AddSecondOperand(double lp,
                                                                  const TSpace &ls,
                                                                  bool loperation)
              {
                 if( (p2!=0) || (s2!=NULL) )
                  {
                     throw "Error in TLinearCombinationWith2Operands::AddSecondOperand: The linear combination is already full.";
                  }

                  p2=lp;
                  s2=&ls;

                  operation=loperation;
              }
        //*************************************
           void TLinearCombinationWith2Operands::AddSecondOperand(const TSpace &ls,
                                                                  bool loperation)
              {
                  if( (p2!=0) || (s2!=NULL) )
                  {
                     throw "Error in TLinearCombinationWith2Operands::AddSecondOperand: The linear combination is already full.";
                  }

                  p2=1;
                  s2=&ls;

                  operation=loperation;
              }
        //*************************************
           void TLinearCombinationWith2Operands::InsertSecondOperand(double lp,
                                                                  const TSpace &ls,
                                                                  bool loperation)
              {
                 if( (p2!=0) || (s2!=NULL) )
                  {
                     throw "Error in TLinearCombinationWith2Operands::AddSecondOperand: The linear combination is already full.";
                  }


                  double helpLD = p1;
                  const TSpace * helpSpace = s1;


                  p1 = lp;
                  s1=&ls;

                  p2=helpLD;
                  s2=helpSpace;

                  operation=loperation;
              }
        //*************************************
           void TLinearCombinationWith2Operands::InsertSecondOperand(const TSpace &ls,
                                                                  bool loperation)
              {
                  if( (p2!=0) || (s2!=NULL) )
                  {
                     throw "Error in TLinearCombinationWith2Operands::AddSecondOperand: The linear combination is already full.";
                  }



                  double helpLD = p1;
                  const TSpace * helpSpace = s1;


                  p1 = 1;
                  s1=&ls;

                  p2=helpLD;
                  s2=helpSpace;

                  operation=loperation;
              }
        //*************************************
           TLinearCombinationWith2Operands &TLinearCombinationWith2Operands::operator+(
                                                          const TLinearCombinationWith2Operands & lc)
              {
                  if( (lc.p2!=0) || (lc.s2!=NULL) )
                  {
                     throw "Error in TLinearCombinationWith2Operands::AddSecondOperand: The linear combination is already full.";
                  }
                  AddSecondOperand(lc.p1,*lc.s1,true);
                  delete &lc;
                  return *this;
              }
        //*************************************
           TLinearCombinationWith2Operands &TLinearCombinationWith2Operands::operator-(
                                                        const TLinearCombinationWith2Operands &lc)
              {
                 if( (lc.p2!=0) || (lc.s2!=NULL) )
                  {
                     throw "Error in TLinearCombinationWith2Operands::AddSecondOperand: The linear combination is already full.";
                  }
                 AddSecondOperand(lc.p1,*lc.s1,false);
                 delete &lc;
                 return *this;
              }
        //*************************************
           TLinearCombinationWith2Operands &TLinearCombinationWith2Operands::operator+(
                                                                             const TSpace &s)
              {
                   AddSecondOperand(s,true);
                   return *this;
              }
        //*************************************
            TLinearCombinationWith2Operands &TLinearCombinationWith2Operands::operator-(
                                                                               const TSpace &s)
              {
                   AddSecondOperand(s,false);
                   return *this;
              }
        //*************************************
      //Конец Класс линейная комбинация с двумя операндами


   //*************************************

     TLinearCombinationWith2Operands &operator*(double p, TLinearSpace &s)
       {
         return s*p;
       }
   //*************************************

//#################################################################

        //*************************************
          TLinearCombinationWith2LinearOperatorOperands::TLinearCombinationWith2LinearOperatorOperands(const TLinearOperator &lA1)

              {
                   A1=&lA1;
                   s2 = NULL;
                   operation=true;
              }
        //*************************************
           TLinearCombinationWith2LinearOperatorOperands::~TLinearCombinationWith2LinearOperatorOperands()
              {
                 //Ничего дополнительного делать не надо....
              }
        //*************************************
           const TSpace *TLinearCombinationWith2LinearOperatorOperands::GetSecondOperand() const
            {
               return s2;
            }
        //*************************************
           void TLinearCombinationWith2LinearOperatorOperands::AddSecondOperand(const TSpace &ls,
                                                                  bool loperation)
              {
                  if (s2!=NULL)
                  {
                     throw "Error in TLinearCombinationWith2LinearOperatorOperands::AddSecondOperand: The linear combination is already full.";
                  }

                  s2=&ls;

                  operation=loperation;
              }
        //*************************************
           TLinearCombinationWith2LinearOperatorOperands &TLinearCombinationWith2LinearOperatorOperands::operator+(
                                                                             const TSpace &s)
              {
                   AddSecondOperand(s,true);
                   return *this;
              }
        //*************************************
            TLinearCombinationWith2LinearOperatorOperands &TLinearCombinationWith2LinearOperatorOperands::operator-(
                                                                               const TSpace &s)
              {
                   AddSecondOperand(s,false);
                   return *this;
              }
        //*************************************

//#################################################################



//---------------------Linear (End)------------------------------






//---------------------Abstract (End)------------------------------


//#################################################################

//*************************************
  TNormalizedSpace::TNormalizedSpace():TLinearSpace() //Конструктор
   {
   }
//*************************************
  TNormalizedSpace::~TNormalizedSpace() //Деструктор
   {
   }
//*************************************
   double TNormalizedSpace::Norma(int type)
   {
     double p=ReleaseNorma(type);
     DeleteIntermediate(*this);
     return p;
   }
//*************************************
   double TNormalizedSpace::operator!()
    {
      return Norma(0);
    }
//*************************************

//#################################################################

//*************************************
 TEuclidSpace::TEuclidSpace():TNormalizedSpace() //Конструктор
  {
  }
//*************************************
 TEuclidSpace::~TEuclidSpace() //Деструктор
  {
  }
//*************************************
 double TEuclidSpace::Innerproduct(const TSpace &s, int type)
  {
     double m = ReleaseInnerproduct(s,type);
     if (&s!=this)
      {
       DeleteIntermediate(s);
      }
     DeleteIntermediate(*this);
     return m;
  }
//*************************************
 double TEuclidSpace::operator,(const TSpace &s)
  {
    return Innerproduct(s,0);
  }
//*************************************
 double TEuclidSpace::Innerproduct(const TLinearCombinationWith2LinearOperatorOperands &lc, int type)
  {
     double m = ReleaseInnerproduct(lc.A1,lc.s2,type);

     TSpace& arg = lc.A1->Argument;
     TSpace *s1 = &arg;
     if ( (this!=s1) && (s1!=NULL)  )
         {
           DeleteIntermediate(*s1);
           ( const_cast<TLinearOperator *>(lc.A1) )->fArgument = NULL;
         }


     if ( (this!=lc.s2) && (lc.s2!=NULL)  )
         {
           DeleteIntermediate(*lc.s2);
         }

     DeleteIntermediate(*this);

     delete &lc;

     return m;
  }
//*************************************
 double TEuclidSpace::operator,(const TLinearCombinationWith2LinearOperatorOperands &lc)
  {
    return Innerproduct(lc,0);
  }
//*************************************
 double operator,(const TLinearCombinationWith2LinearOperatorOperands &lc,const TEuclidSpace &s)
  {
    TEuclidSpace &sNonConst = const_cast<TEuclidSpace &>(s);
    return (sNonConst,lc);
  }
//*************************************

//*************************************
 double TEuclidSpace::NormaAccordingInnerproduct(int type) //Норма, согласованная со скалярным произведением
  {
    return sqrtl(Innerproduct(*this,type));
  }
//*************************************

//#################################################################

//*************************************
  TBanahSpace::TBanahSpace():TEuclidSpace(),TFullSpace() //Конструктор
   {
   }
//*************************************
  TBanahSpace::~TBanahSpace() //Деструктор
   {
   }
//*************************************
  double TBanahSpace::MetricAccordingNorma(const TSpace &s,int type)
   {
      TBanahSpace &temp=dynamic_cast<TBanahSpace &>(this->NewItem());
      temp|=(*this-s);
      double n=temp.Norma(type);
      delete &temp;
      return n;
   }
//*************************************

//#################################################################


//---------------------Abstract (End)------------------------------






//---------------------Rn (Begin)------------------------------

//#################################################################

//*************************************
   TRnSpaceModelAndItsName::TRnSpaceModelAndItsName(int lName,TRnSpaceModel *lModel)
    {
      fName=lName;
      fModel=lModel;
      ReferencesCount=0;
    }
//*************************************    
   TRnSpaceModelAndItsName::~TRnSpaceModelAndItsName() {};
//*************************************
   TRnSpaceModelAndItsName::TRnSpaceModelAndItsName(TRnSpaceModelAndItsName &) {}; //Пустышка
   void TRnSpaceModelAndItsName::operator=(TRnSpaceModelAndItsName &) {}; //Пустышка
//*************************************

//#################################################################

//*************************************
  TRnSpaceModelAndItsName **TRnSpace::fModels=NULL;
  unsigned char TRnSpace::fModelsCount=0;
//*************************************
  TRnSpaceModel *TRnSpace::GetModel(int lName)
   {
      TRnSpaceModel *res=NULL;
      bool flag=false;
      int i;
      for (i=0; (i<fModelsCount) && (flag==false); i++)
       {
          if (fModels[i]->fName==lName)
           {
             flag=true;
             res=fModels[i]->fModel;
           }
       }
      return res;
   }
//*************************************
  void TRnSpace::IncreaseReferencesCountToModel(int lName)
   {
      TRnSpaceModelAndItsName *res=NULL;
      int i;
      bool flag=false;
      for (i=0; (i<fModelsCount) && (flag==false); i++)
       {
          if (fModels[i]->fName==lName)
           {
             flag=true;
             res=fModels[i];
           }
       }
      if (res==NULL)
       {
         throw "TRnSpace::IncreaseReferencesCountToModel: Model not found.";
       }
      else {res->ReferencesCount=res->ReferencesCount+1;}
   }
//*************************************
  void TRnSpace::DecreaseReferencesCountToModel(int lName)
   {
      TRnSpaceModelAndItsName *res=NULL;
      int i;
      bool flag=false;
      for (i=0; (i<fModelsCount) && (flag==false); i++)
       {
          if (fModels[i]->fName==lName)
           {
             flag=true;
             res=fModels[i];
           }
       }
      if (res==NULL)
       {
         throw "TRnSpace::DecreaseReferencesCountToModel: Model not found.";
       }
      else {res->ReferencesCount=res->ReferencesCount-1;}
   }

//*************************************

  void TRnSpace::AddModel(int lName,TRnSpaceModel *lModel,bool deleteModelAfterAdding)
   {
      if (GetModel(lName)==NULL)
       {
          unsigned char helpModelsCount=fModelsCount;
          fModelsCount=fModelsCount+1;
          TRnSpaceModelAndItsName **helpModels=fModels;
          fModels=new TRnSpaceModelAndItsName *[fModelsCount];
          if (helpModels==fModels) {/*ShowMessage("1  !==!");*/}
          int i;
          for (i=0;i<helpModelsCount;i++)
           {
             fModels[i]=helpModels[i];
           }
           if (helpModelsCount>0) {delete[] helpModels;}
           fModels[fModelsCount-1]=new TRnSpaceModelAndItsName(lName,&(lModel->Copy()));
           if (deleteModelAfterAdding==true)
            {
               delete lModel;
            }

       }
      else {throw "Error in TRnSpace::AddModel: model is already exists.";}
   }
//*************************************
  void TRnSpace::RemoveAllModels()
   {
     int i;
     int help;
     //Проверка возможности удаления моделей
       for (i=0;i<fModelsCount;i++)
           {
             help=fModels[i]->ReferencesCount;
             if (help>0)
              {
                throw "Error in TRnSpace::RemoveAllModels:Can't remove model. There are references.";
              }
           }
     //Конец Проверка возможности удаления моделей

     for (i=0;i<fModelsCount;i++)
           {
             delete fModels[i]->fModel;
             delete fModels[i];
           }
     delete[] fModels;
     fModels=NULL;
     fModelsCount=0;
   }
//*************************************
   void TRnSpace::RemoveModel(int lName)
    {
       TRnSpaceModel *res = GetModel(lName);
       if (res!=NULL)
        {
          unsigned char helpModelsCount=fModelsCount;
          fModelsCount=fModelsCount-1;
          if (fModelsCount>0)
           {
             TRnSpaceModelAndItsName **helpModels=fModels;
             fModels=new TRnSpaceModelAndItsName *[fModelsCount];
             if (helpModels==fModels) {/*ShowMessage("1  !==!");*/}
             int i;
             int j = -1;
             for (i=0;i<helpModelsCount;i++)
              {
                if (helpModels[i]->fModel!=res)
                 {
                   j = j+1;
                   fModels[j]=helpModels[i];
                 }
                else
                 {
                   delete helpModels[i]->fModel;
                   delete helpModels[i];
                 }
              }
             if (helpModelsCount>0) {delete[] helpModels;}
           else
            {
              delete fModels[0]->fModel;
              delete fModels[0];
              delete fModels;// fModels = NULL;
            }
        }
       else
       {
         throw "Error in TRnSpace::RemoveModel: model not found.";
       }

     }

    }
//*************************************

  TRnSpace::TRnSpace(const int lModelName):TBanahSpace(),
        EndIndex(this), ActualIndexesCount(this), ActualIndexes(this), AutomaticActualIndexesAssigning(this) //Конструктор
   {
       TRnSpaceModel *model;
       model=GetModel(lModelName);
       if (model!=NULL)
        {
          fModelName=lModelName;
          fModel=model;
          IncreaseReferencesCountToModel(fModelName);
        }
        else
        {
          throw "Error in TRnSpace constructor: Model not found.";
        }

        fActualIndexesCount = 0;
        fActualIndexes = NULL;
        fAutomaticActualIndexesAssigning = true;

   }
//*************************************
  TRnSpace::~TRnSpace() //Деструктор
   {
     DecreaseReferencesCountToModel(fModelName);
   }
//*************************************
  TRnSpace &TRnSpace::CreateInstance()
   {
     return dynamic_cast<TRnSpace &>(TSpace::CreateInstance());
   }
//*************************************
  TRnSpace &TRnSpace::Copy()
   {
     return dynamic_cast<TRnSpace &>(TSpace::Copy());
   }
//*************************************
  TRnSpace &TRnSpace::operator|=(const TRnSpace &s)
   {
      const TRnSpace &sCast=s;
      if (fModel!=sCast.fModel)
      {throw "Error in TRnSpace::operator=: models mismatch.";}
      return dynamic_cast<TRnSpace &>(TSpace::operator|=(s));
   }
//*************************************
  TRnSpace &TRnSpace::Assign(const TRnSpace &s)
   {
      const TRnSpace &sCast=s;
      if (fModel!=sCast.fModel)
      {throw "Error in TRnSpace::Assign: models mismatch.";}
      return dynamic_cast<TRnSpace &>(TSpace::Assign(s));
   }
//*************************************
  TRnSpace &TRnSpace::operator|=(const TLinearCombinationWith2Operands &lc)
   {
      const TRnSpace &s1Cast=dynamic_cast<const TRnSpace &>(*lc.GetFirstOperand());
      if (fModel!=s1Cast.fModel)
      {throw "Error in TRnSpace::operator=: models mismatch.";}

      if (lc.GetSecondOperand()!=NULL)
       {
         const TRnSpace &s2Cast=dynamic_cast<const TRnSpace &>(*lc.GetSecondOperand());
         if (fModel!=s2Cast.fModel)
         {throw "Error in TRnSpace::operator=: models mismatch.";}
       }

      return dynamic_cast<TRnSpace &>(TLinearSpace::operator|=(lc));
   }
//*************************************
  TRnSpace &TRnSpace::operator|=(const TLinearCombinationWith2LinearOperatorOperands &lc)
   {
      if (lc.GetSecondOperand()!=NULL)
       {
         const TRnSpace &s2Cast=dynamic_cast<const TRnSpace &>(*lc.GetSecondOperand());
         if (fModel!=s2Cast.fModel)
         {throw "Error in TRnSpace::operator=: models mismatch.";}
       }

       //!!! И надо для скалярного произведения модели проверять...

      //Попытка приведения оператора к оператору над Rn
        //Это надо сделать, как и проверку приведеения пространств во всех функциях...
      //Конец Попытка приведения оператора к оператору над Rn

      return dynamic_cast<TRnSpace &>(TLinearSpace::operator|=(lc));
   }
//*************************************
  TRnSpace &TRnSpace::operator|=(const int zero)
   {
      return dynamic_cast<TRnSpace &>(TLinearSpace::operator|=(zero));
   }
//*************************************
  double TRnSpace::Metric(const TRnSpace &s, int type)
   {
       const TRnSpace &sCast=s;
       if (fModel!=sCast.fModel)
       {throw "Error in TRnSpace::Metric: models mismatch.";}
       return TMetricSpace::Metric(s,type);
   }
//*************************************
  TLinearCombinationWith2Operands &TRnSpace::operator+(const TRnSpace &s)
   {
      return TLinearSpace::operator+(s);
   }
//*************************************
  TLinearCombinationWith2Operands &TRnSpace::operator-(const TRnSpace &s)
   {
      return TLinearSpace::operator-(s);
   }
//*************************************
  TLinearCombinationWith2Operands &TRnSpace::operator+(TLinearCombinationWith2Operands &lc)
   {
      return TLinearSpace::operator+(lc);
   }
//*************************************
  TLinearCombinationWith2Operands &TRnSpace::operator-(TLinearCombinationWith2Operands &lc)
   {
      return TLinearSpace::operator-(lc);
   }
//*************************************
  TLinearCombinationWith2Operands &TRnSpace::operator*(double p)
   {
      return TLinearSpace::operator*(p);
   }
//*************************************
  double TRnSpace::Innerproduct(const TRnSpace &s, int type)
   {
       const TRnSpace &sCast=s;
       if (fModel!=sCast.fModel)
       {throw "Error in TRnSpace::Innerproduct: models mismatch.";}
       return TEuclidSpace::Innerproduct(s,type);
   }
//*************************************

  //Вспомогательные функции и классы
     TLinearCombinationWith2Operands &operator*(double p, TRnSpace &s)
      {
         return s*p;
      }
  //Конец Вспомогательные функции и классы

//*************************************

  double &TRnSpace::operator[](int i) const
    {
       if ( (i<1) || (i>EndIndex) )
        {
          throw "Error in TRnSpace::operator[]: The index is out of range.";
        }

        return GetByRnIndex(i);

    }
//*************************************

  int TRnSpace::GetEndIndex() const
    {
        return fModel->EndIndex;

    }
//*************************************
  void TRnSpace::SetActualIndexesCount(const int& lCount)
    {
      if (lCount<0)
        {
          throw "TRnSpace::SetActualIndexesCount: The count is less then zero.";
        }
      delete[] fActualIndexes;// fActualIndexes = NULL;
      fActualIndexesCount = lCount;
      if (lCount>0)
       {
         fActualIndexes = new int[fActualIndexesCount];
         int i;
         for (i=0;i<fActualIndexesCount;i++)
           {
             fActualIndexes[i]=0;
           }
       }
    }
//*************************************
 double TRnSpace::GetAbsMaxElement()
  {
    int index = 1;
    double res;
    res = GetSpecificElement(1,index);
    if ( (index<1) || (index>EndIndex) || (res<0) )
     {
       throw "Error in TRnSpace::GetAbsMaxElement: Abnormal situation.";
     }
    
    return res;
  }
//*************************************
 double TRnSpace::GetAbsMaxElement(int &elementIndex)
  {
    int index = 1;
    double res;
    res = GetSpecificElement(1,index);
    if ( (index<1) || (index>EndIndex) || (res<0) )
     {
       throw "Error in TRnSpace::GetAbsMaxElement: Abnormal situation.";
     }

    elementIndex = index;
    return res;
  }
//*************************************
 double TRnSpace::GetAbsMinElement()
  {
    int index = 1;
    double res;
    res = GetSpecificElement(2,index);
    if ( (index<1) || (index>EndIndex) || (res<0) )
     {
       throw "Error in TRnSpace::GetAbsMinElement: Abnormal situation.";
     }

    return res;
  }
//*************************************
 double TRnSpace::GetAbsMinElement(int &elementIndex)
  {
    int index = 1;
    double res;
    res = GetSpecificElement(2,index);
    if ( (index<1) || (index>EndIndex) || (res<0) )
     {
       throw "Error in TRnSpace::GetAbsMinElement: Abnormal situation.";
     }

    elementIndex = index;
    return res;
  }
//*************************************

//#################################################################

TRnData::TRnData() {}
TRnData::~TRnData() {}

//#################################################################

//*************************************
  TRnReleaseMultiDataSpace::TRnReleaseMultiDataSpace(const int lModelName): TRnSpace(lModelName),
        ActiveDataIndex(this), MultiData(this) //Конструктор
   {
      //Приведение модели к TRnRelease2dSpaceModel
       try
        {
          fCastModel=dynamic_cast<TRnReleaseMultiDataSpaceModel *>(fModel);
          if (fCastModel==NULL) throw "Error in TRnReleaseMultiDataSpace constructor: Bad cast.";
        }
       catch (...) {throw "Error in TRnReleaseMultiDataSpace constructor: Could not cast model to TRnReleaseMultiDataSpaceModel type.";}
      //Конец Приведение модели к TRnRelease2dSpaceModel

      fMultiData = new TRnData * [fCastModel->DataCount];
      fActiveDataIndex = 0;
   }
//*************************************
  TRnReleaseMultiDataSpace::~TRnReleaseMultiDataSpace()
   {
      delete fMultiData;// fMultiData = NULL;
   }
//*************************************
  void TRnReleaseMultiDataSpace::SetActiveDataIndex(const int& lIndex)
   {
     if ( (lIndex<0) || (lIndex>fCastModel->DataCount-1) )
      {
        throw "Error in TRnReleaseMultiDataSpace SetActiveDataIndex: Index is out of range.";
      }
     fActiveDataIndex = lIndex;
   }
//*************************************
  void TRnReleaseMultiDataSpace::ReleaseOperatorEqual(const TSpace &s,bool isCopying)
   {
      //Попытка приведения параметра к TRnReleaseMultiDataSpace
         try {dynamic_cast<const TRnReleaseMultiDataSpace &>(s);}
         catch (...) {throw "Error in TRnReleaseMultiDataSpace::ReleaseOperatorEqual: Could not cast parametr pased to TRnReleaseMultiDataSpace & type.";}
      //Конец Попытка приведения параметра к TRnReleaseMultiDataSpace

      int dataCount = fCastModel->DataCount;
      int i;
      for (i=0;i<dataCount;i++)
       {
          LocalReleaseOperatorEqual(i,*dynamic_cast<const TRnReleaseMultiDataSpace &>(s).fMultiData[i],isCopying);
       }

   }
//*************************************
  double TRnReleaseMultiDataSpace::ReleaseMetric(const TSpace &s,int type)
   {
      //Попытка приведения параметра к TRnReleaseMultiDataSpace
         try {dynamic_cast<const TRnReleaseMultiDataSpace &>(s);}
         catch (...) {throw "Error in TRnReleaseMultiDataSpace::ReleaseMetric: Could not cast parametr pased to TRnReleaseMultiDataSpace & type.";}
      //Конец Попытка приведения параметра к TRnReleaseMultiDataSpace

      double res = 0;
      int dataCount = fCastModel->DataCount;
      int i;
      for (i=0;i<dataCount;i++)
       {
          res = res + LocalReleaseMetric(i,*dynamic_cast<const TRnReleaseMultiDataSpace &>(s).fMultiData[i],type);
       }

       res = sqrtl(res);

       return res;
   }
//*************************************
  void TRnReleaseMultiDataSpace::ToZero()
   {
      int dataCount = fCastModel->DataCount;
      int i;
      for (i=0;i<dataCount;i++)
       {
          LocalToZero(i);
       }
   }
//*************************************
  void TRnReleaseMultiDataSpace::ReleaseAssignSpeedLinearCombination(double p1,const TSpace *s1,
                                      double p2,const TSpace *s2,bool operation)
   {
     //Попытка приведения параметра к TRnReleaseMultiDataSpace
      if (s1!=NULL)
       {
         try {dynamic_cast<const TRnReleaseMultiDataSpace *>(s1);}
         catch (...) {throw "Error in TRnReleaseMultiDataSpace::ReleaseAssignSpeedLinearCombination: Could not cast parametr pased to TRnReleaseMultiDataSpace * type.";}
       }
     //Конец Попытка приведения параметра к TRnReleaseMultiDataSpace

     //Попытка приведения параметра к TRnReleaseMultiDataSpace
      if (s2!=NULL)
       {
         try {dynamic_cast<const TRnReleaseMultiDataSpace *>(s2);}
         catch (...) {throw "Error in TRnReleaseMultiDataSpace::ReleaseAssignSpeedLinearCombination: Could not cast parametr pased to TRnReleaseMultiDataSpace * type.";}
       }
     //Конец Попытка приведения параметра к TRnReleaseMultiDataSpace


      int dataCount = fCastModel->DataCount;
      int i;
      for (i=0;i<dataCount;i++)
       {
         if ( (s1!=NULL) && (s2!=NULL) )
          {
            LocalReleaseAssignSpeedLinearCombination(i,p1,dynamic_cast<const TRnReleaseMultiDataSpace *>(s1)->fMultiData[i],p2,dynamic_cast<const TRnReleaseMultiDataSpace *>(s2)->fMultiData[i],operation);
          }
         else if  ( (s1!=NULL) && (s2==NULL) )
          {
            LocalReleaseAssignSpeedLinearCombination(i,p1,dynamic_cast<const TRnReleaseMultiDataSpace *>(s1)->fMultiData[i],0,NULL,operation);
          }
         else if  ( (s1==NULL) && (s2!=NULL) )
          {
            LocalReleaseAssignSpeedLinearCombination(i,0,NULL,p2,dynamic_cast<const TRnReleaseMultiDataSpace *>(s2)->fMultiData[i],operation);
          }
         else
          {
            LocalReleaseAssignSpeedLinearCombination(i,0,NULL,0,NULL,operation);
          }
       }
   }
//*************************************
   void TRnReleaseMultiDataSpace::ReleaseAssignSpeedLinearCombination(const TLinearOperator *A1,
                                            const TSpace *s2,bool operation)
    {
      //Попытка приведения параметра к TRnReleaseMultiDataSpace
      if (s2!=NULL)
       {
         try {dynamic_cast<const TRnReleaseMultiDataSpace *>(s2);}
         catch (...) {throw "Error in TRnReleaseMultiDataSpace::ReleaseAssignSpeedLinearCombination: Could not cast parametr pased to TRnReleaseMultiDataSpace * type.";}
       }
     //Конец Попытка приведения параметра к TRnReleaseMultiDataSpace



      //Работа с актуальными индексами (Начало)
        const TRnLinearOperator *A1Cast = NULL;
        if (A1 == NULL)
       {
          throw "Error in TRnReleaseMultiDataSpace::ReleaseAssignSpeedLinearCombination: First operator is NULL.";
       }

       //Попытка приведения параметра к TRnLinearOperator
       try {A1Cast = dynamic_cast<const TRnLinearOperator *>(A1);}
       catch (...) {throw "Error in TRnReleaseMultiDataSpace::ReleaseAssignSpeedLinearCombination: Could not cast parametr pased to TRnLinearOperator * type.";}
       if (A1Cast == NULL)
        {
          throw "Error in TRnReleaseMultiDataSpace::ReleaseAssignSpeedLinearCombination: Could not cast parametr pased to TRnLinearOperator * type.";
        }
       //Конец Попытка приведения параметра к TRnLinearOperator


       if (this->AutomaticActualIndexesAssigning==true)
        {
          if (this->ActualIndexesCount!=A1Cast->ActualIndexesCount)
           {
             this->SetActualIndexesCount(A1Cast->ActualIndexesCount);

           }

          int i;
          int count = this->ActualIndexesCount;
             for (i=0;i<count;i++)
               {
                 this->ActualIndexes[i] = A1Cast->ActualIndexes[i];
               }
        }
      //Работа с актуальными индексами (Конец)



      int dataCount = fCastModel->DataCount;
      int i;
      for (i=0;i<dataCount;i++)
       {

         if ( (A1!=NULL) && (s2!=NULL) )
          {
            LocalReleaseAssignSpeedLinearCombination(i,A1,dynamic_cast<const TRnReleaseMultiDataSpace *>(s2)->fMultiData[i],operation);
          }
         else if  ( (A1!=NULL) && (s2==NULL) )
          {
            LocalReleaseAssignSpeedLinearCombination(i,A1,NULL,operation);
          }
         else if  ( (A1==NULL) && (s2!=NULL) )
          {
            LocalReleaseAssignSpeedLinearCombination(i,NULL,dynamic_cast<const TRnReleaseMultiDataSpace *>(s2)->fMultiData[i],operation);
          }
         else
          {
            LocalReleaseAssignSpeedLinearCombination(i,NULL,NULL,operation);
          }
       }

    }
//*************************************

  double TRnReleaseMultiDataSpace::ReleaseNorma(int type)
   {
      double res = 0;
      int dataCount = fCastModel->DataCount;
      int i;
      for (i=0;i<dataCount;i++)
       {
          res = res + LocalReleaseNorma(i,type);
       }

       res = sqrtl(res);

       return res;
   }
//*************************************
  double TRnReleaseMultiDataSpace::ReleaseInnerproduct(const TSpace &s, int type)
   {
      //Попытка приведения параметра к TRnReleaseMultiDataSpace
         try {dynamic_cast<const TRnReleaseMultiDataSpace &>(s);}
         catch (...) {throw "Error in TRnReleaseMultiDataSpace::ReleaseInnerproduct: Could not cast parametr pased to TRnReleaseMultiDataSpace & type.";}
      //Конец Попытка приведения параметра к TRnReleaseMultiDataSpace

      double res = 0;
      int dataCount = fCastModel->DataCount;
      int i;
      for (i=0;i<dataCount;i++)
       {
          res = res + LocalReleaseInnerproduct(i,*dynamic_cast<const TRnReleaseMultiDataSpace &>(s).fMultiData[i],type);
       }

      return res;

   }
//*************************************
 double TRnReleaseMultiDataSpace::ReleaseInnerproduct(const TLinearOperator *lA1,const TSpace *s2, int type)
   {
      const TLinearOperator *A1 = lA1;
      //Попытка приведения параметра к TRnReleaseMultiDataSpace
      if (s2!=NULL)
       {
         try {dynamic_cast<const TRnReleaseMultiDataSpace *>(s2);}
         catch (...) {throw "Error in RnReleaseMultiDataSpace::ReleaseInnerproduct: Could not cast parametr pased to TRnReleaseMultiDataSpace * type.";}
       }
     //Конец Попытка приведения параметра к TRnReleaseMultiDataSpace

      double res = 0;
      int dataCount = fCastModel->DataCount;
      int i;
      for (i=0;i<dataCount;i++)
       {

         if ( (A1!=NULL) && (s2!=NULL) )
          {
            res = res + LocalReleaseInnerproduct(i,A1,dynamic_cast<const TRnReleaseMultiDataSpace *>(s2)->fMultiData[i],type);
          }
         else if  ( (A1!=NULL) && (s2==NULL) )
          {
            res = res + LocalReleaseInnerproduct(i,A1,NULL,type);
          }
         else if  ( (A1==NULL) && (s2!=NULL) )
          {
            res = res + LocalReleaseInnerproduct(i,NULL,dynamic_cast<const TRnReleaseMultiDataSpace *>(s2)->fMultiData[i],type);
          }
         else
          {
            res = res + LocalReleaseInnerproduct(i,NULL,NULL,type);
          }
       }

       return res;
   }
//*************************************
 double TRnReleaseMultiDataSpace::GetSpecificElement(int type,int &elementIndex) const
  {

      double res = 0;

      int dataCount = fCastModel->DataCount;
      int i;

      double localRes;
      int localIndex = 1;
      for (i=0;i<dataCount;i++)
       {
          localRes = fabs(LocalGetSpecificElement(i,type,localIndex));

          if (type==1)
           {
            if (localRes>res) {res = localRes; elementIndex = localIndex;}
           }
          else if (type==2)
           {
            if (localRes<res) {res = localRes; elementIndex = localIndex;}
           }
          else
           {
             throw "Error in TRnReleaseMultiDataSpace::GetSpecificElement: The fuction is not implemented for this type.";
           }
       }

      return res;

  }
 //*************************************

 //#################################################################


 //---------------------Rn (End)------------------------------





//---------------------Rn1 (Begin)------------------------------


 //#################################################################

   //*************************************
       TRnRelease1dData::TRnRelease1dData(const int lN): TRnData(),fN(lN), Data(this)
        {
           //Выделение памяти
             int n=fN;
             int i;
             fData=new double [n];
           //Конец Выделение памяти

           //Обнуление всех элементов (Начало)

             for (i=0;i<n;i++)
               {
                    fData[i]=0;
               }
           //Обнуление всех элементов (Конец)
        }
   //*************************************
       TRnRelease1dData::~TRnRelease1dData()
        {
          delete[] fData;
          fData=NULL;
        }
   //*************************************
       void TRnRelease1dData::HelpVirtualFunction() const
        {

        }
   //*************************************

   //#################################################################

//*************************************
  TRnRelease1dSpace::TRnRelease1dSpace(const int lModelName): TRnReleaseMultiDataSpace(lModelName) //Конструктор
   {
      //Приведение модели к TRnRelease1dSpaceModel
      try {fCastModel=dynamic_cast<TRnRelease1dSpaceModel *>(fModel);}
      catch (...) {throw "Error in TRnRelease1dSpace constructor: Could not cast model to TRnRelease1dSpaceModel type.";}
      //Конец Приведение модели к TRnRelease1dSpaceModel

      int dataCount = fCastModel->DataCount;
      for (int i=0;i<dataCount;i++)
       {
          int n = (fCastModel->GetMask(i)).N;
          fMultiData[i] = new TRnRelease1dData(n);
       }
   }
//*************************************
  TRnRelease1dSpace::~TRnRelease1dSpace() //Деструктор
   {
      int dataCount = fCastModel->DataCount;
      for (int i=0;i<dataCount;i++)
       {
          delete fMultiData[i];// fMultiData[i] = NULL;
       }
   }
//*************************************
  double &TRnRelease1dSpace::operator[](int i)
   {
     return dynamic_cast<TRnRelease1dData *>(fMultiData[fActiveDataIndex])->Data[i];
   }
//*************************************
  TRnSpace &TRnRelease1dSpace::NewItem()
   {
     return *(new TRnRelease1dSpace(fModelName));
   }
//*************************************
   void TRnRelease1dSpace::LocalReleaseOperatorEqual(const int id,const TRnData &s,bool isCopying)
   {

        //Попытка приведения параметра к TRnRelease1dData
         try {dynamic_cast<const TRnRelease1dData &>(s);}
         catch (...) {throw "Error in TRnRelease1dData::LocalReleaseOperatorEqual: Could not cast parametr pased to TRnRelease1dData & type.";}
        //Конец Попытка приведения параметра к TRnRelease1dData

        const TRnRelease1dData &sCast=dynamic_cast<const TRnRelease1dData &>(s);
        //Выполнение присваивания
        const T1dNumberMask & sMask = fCastModel->GetMask(id);
        const int *mask = sMask.Mask;
        int n = sMask.N;
        int i;

        double *sCastData=sCast.fData;
        double *fCastData=dynamic_cast<const TRnRelease1dData &>(*fMultiData[id]).fData;

        if (isCopying==false)
         {
           for (i=0;i<n;i++)
            {
				   //Жертва-ускорение: будем использовать этот тип вектора только для упакованных плотных пространств, т.е. маска - не актуальна и не
                   //должна содержать 0 
                   fCastData[i]=sCastData[i];
                   /*if (mask[i]>0) {fCastData[i]=sCastData[i];}*/
            }
         }
        else
         {
           for (i=0;i<n;i++)
            {
                   fCastData[i]=sCastData[i];
            }
         }
        //Конец Выполнение присваивания
   }
//*************************************
   double TRnRelease1dSpace::LocalReleaseMetric(const int id,const TRnData &s,int type)
   {
        //Попытка приведения параметра к TRnRelease1dData
         try {dynamic_cast<const TRnRelease1dData &>(s);}
         catch (...) {throw "Error in TRnRelease1dData::LocalReleaseMetric: Could not cast parametr pased to TRnRelease1dData & type.";}
        //Конец Попытка приведения параметра к TRnRelease1dData

        const TRnRelease1dData &sCast=dynamic_cast<const TRnRelease1dData &>(s);

        double res=0;
        //Выполнение вычисления метрики
        const T1dNumberMask & sMask = fCastModel->GetMask(id);
        const int *mask = sMask.Mask;
        int n = sMask.N;
        int i;

        double *sCastData=sCast.fData;
        double *fCastData=dynamic_cast<const TRnRelease1dData &>(*fMultiData[id]).fData;

        double help;

        for (i=0;i<n;i++)
          {
			    //Жертва-ускорение: будем использовать этот тип вектора только для упакованных плотных пространств, т.е. маска - не актуальна и не
                //должна содержать 0
			    help=fCastData[i]-sCastData[i];
                res=res+help*help;
			    /*
                if (mask[i]>0)
                 {
                  help=fCastData[i]-sCastData[i];
                  res=res+help*help;
                 }
                */
          }
        //Конец Выполнение вычисления метрики
        return res;
   }
//*************************************
   void TRnRelease1dSpace::LocalToZero(const int id)
   {
        const T1dNumberMask & sMask = fCastModel->GetMask(id);
        const int *mask = sMask.Mask;
        int n = sMask.N;
        int i;

        double *fCastData=dynamic_cast<const TRnRelease1dData &>(*fMultiData[id]).fData;

        int actIndsCount = ActualIndexesCount;
        int *actInds = ActualIndexes;
        int loopVar;
        int indexRn;

        if (actIndsCount>0)
         {
            for (loopVar=0;loopVar<actIndsCount;loopVar++)
                        {
                          indexRn = actInds[loopVar];
                          const TCorrespondence &cor = fCastModel->GetCorrespondence(indexRn);
                          const T1dCorrespondence &cor1d = dynamic_cast<const T1dCorrespondence &>(cor);
                          int ind = cor1d.Ind;
                          if (ind==id)
                           {
                             i = cor1d.I;
							 //Жертва-ускорение: будем использовать этот тип вектора только для упакованных плотных пространств, т.е. маска - не актуальна и не
                             //должна содержать 0
							 fCastData[i]=0;
							 /*
                             if (mask[i]>0)
                               {
                                 fCastData[i]=0;
                               }
                             */ 
                           } //Если индексы равны
                        } //Цикл по актуальным индексам
                    } //Если есть акт. индексы
        else
        {

           for (i=0;i<n;i++)
            {
                //Жертва-ускорение: будем использовать этот тип вектора только для упакованных плотных пространств, т.е. маска - не актуальна и не
                //должна содержать 0 
				fCastData[i]=0;
				/*if (mask[i]>0) {fCastData[i]=0;}*/
            }

        }
   }
//*************************************
   void TRnRelease1dSpace::LocalReleaseAssignSpeedLinearCombination(const int id,double p1,const TRnData *s1,
                                      double p2,const TRnData *s2,bool operation)
   {
       //Может быть несколько лишне организована проверка на то, что
       //один из сомножителей =1 - может это делает компилятор (хотя как ?)


        int variant; //-1 не использовать операнды
                     //0 - оба операнда
                     //1 - использовать первый
                     //2 - использовать второй



        if  (s1!=NULL)
        {
        //Попытка приведения параметра к TRnRelease1dData
         try {dynamic_cast<const TRnRelease1dData *>(s1);}
         catch (...) {throw "Error in TRnRelease1dData::LocalReleaseAssignSpeedLinearCombination: Could not cast parametr pased to TRnRelease1dData * type.";}
        //Конец Попытка приведения параметра к TRnRelease1dData
        }

        if  (s2!=NULL)
        {
        //Попытка приведения параметра к TRnRelease1dData
         try {dynamic_cast<const TRnRelease1dData *>(s2);}
         catch (...) {throw "Error in TRnRelease1dData::LocalReleaseAssignSpeedLinearCombination: Could not cast parametr pased to TRnRelease1dData * type.";}
        //Конец Попытка приведения параметра к TRnRelease1dData
        }


       //Подстраховка ...

         if ( ( (p1!=0)&&(s1==NULL) ) || ( (p2!=0)&&(s2==NULL) ) )
          {
            throw "Error in TRnRelease1dData::LocalReleaseAssignSpeedLinearCombination: Incorrect combination: p<>0 but s=NULL.";
          }

       //Конец Подстраховка ...



       //Определение варианта
         if (p1!=0)
          {
            if (p2!=0) {variant=0;}
            else {variant=1;}
          }
         else
          {
             if (p2!=0) {variant=2;}
             else {variant=-1;}

          }
       //Конец Определение варианта



       const T1dNumberMask & sMask = fCastModel->GetMask(id);
       const int *mask = sMask.Mask;
       int n = sMask.N;
       int i;

       double *resCastData=dynamic_cast<const TRnRelease1dData &>(*fMultiData[id]).fData;

       //Отсюда начать...

       if (variant==-1) {LocalToZero(id); return;}
       else if (variant==1)
         {
           if (p1!=1)
             {
                 //Выполнение быстрой линейной комбинации
                const TRnRelease1dData &s1Cast=dynamic_cast<const TRnRelease1dData &>(*s1);
                double *s1CastData=s1Cast.fData;
                for (i=0;i<n;i++)
                       {
						       //Жертва-ускорение
                               //if (mask[i]>0)
                                  //{
                                    resCastData[i]=p1*s1CastData[i];
                                  //}

                        }
                //Конец Выполнение быстрой линейной комбинации
             }
           else
             {
                 //Выполнение быстрой линейной комбинации
                const TRnRelease1dData &s1Cast=dynamic_cast<const TRnRelease1dData &>(*s1);
                double *s1CastData=s1Cast.fData;
                for (i=0;i<n;i++)
                       {
						       //Жертва-ускорение
                               //if (mask[i]>0)
                                  //{
                                    resCastData[i]=s1CastData[i];
                                  //}

                        }
                //Конец Выполнение быстрой линейной комбинации
             }
         }
       else if (variant==2)
         {
            if (p2!=1)
             {
                //Выполнение быстрой линейной комбинации
                const TRnRelease1dData &s2Cast=dynamic_cast<const TRnRelease1dData &>(*s2);
                double *s2CastData=s2Cast.fData;
                if (operation==true)
                   {
                     for (i=0;i<n;i++)
                       {
						       //Жертва-ускорение
                               //if (mask[i]>0)
                                  //{
                                    resCastData[i]=p2*s2CastData[i];
                                  //}

                        }
                   }
                else
                   {
                     for (i=0;i<n;i++)
                       {
						       //Жертва-ускорение
                               //if (mask[i]>0)
                                   //{
                                     resCastData[i]=-p2*s2CastData[i];
                                   //}
                       }
                    }
               //Конец Выполнение быстрой линейной комбинации

             }
           else
             {
                 //Выполнение быстрой линейной комбинации
                const TRnRelease1dData &s2Cast=dynamic_cast<const TRnRelease1dData &>(*s2);
                double *s2CastData=s2Cast.fData;
                if (operation==true)
                   {
                     for (i=0;i<n;i++)
                       {
						       //Жертва-ускорение
                               //if (mask[i]>0)
                                  //{
                                    resCastData[i]=s2CastData[i];
                                  //}
                        }
                   }
                else
                   {
                     for (i=0;i<n;i++)
                       {
						       //Жертва-ускорение
                               //if (mask[i]>0)
                                   //{
                                     resCastData[i]=-s2CastData[i];
                                   //}
                       }
                    }
               //Конец Выполнение быстрой линейной комбинации
             }
         }
       else
         {

           if ( (p1!=1) && (p2!=1) )
             {

                //Выполнение быстрой линейной комбинации
                const TRnRelease1dData &s1Cast=dynamic_cast<const TRnRelease1dData &>(*s1);
                const TRnRelease1dData &s2Cast=dynamic_cast<const TRnRelease1dData &>(*s2);
                double *s1CastData=s1Cast.fData;
                double *s2CastData=s2Cast.fData;
                if (operation==true)
                   {
                     for (i=0;i<n;i++)
                       {
						       //Жертва-ускорение
                               //if (mask[i]>0)
                                  //{
                                    resCastData[i]=p1*s1CastData[i]+p2*s2CastData[i];
                                  //}
                        }
                   }
                else
                   {
                     for (i=0;i<n;i++)
                       {
						        //Жертва-ускорение
                                //if (mask[i]>0)
                                  //{
                                     resCastData[i]=p1*s1CastData[i]-p2*s2CastData[i];
                                  //}
                       }
                    }
               //Конец Выполнение быстрой линейной комбинации

           }
           else if ( (p1!=1) && (p2==1) )
             {
                //Выполнение быстрой линейной комбинации
                const TRnRelease1dData &s1Cast=dynamic_cast<const TRnRelease1dData &>(*s1);
                const TRnRelease1dData &s2Cast=dynamic_cast<const TRnRelease1dData &>(*s2);
                double *s1CastData=s1Cast.fData;
                double *s2CastData=s2Cast.fData;
                if (operation==true)
                   {
                     for (i=0;i<n;i++)
                       {
						       //Жертва-ускорение
                               //if (mask[i]>0)
                                  //{
                                    resCastData[i]=p1*s1CastData[i]+s2CastData[i];
                                  //}


                        }
                   }
                else
                   {
                     for (i=0;i<n;i++)
                       {
						         //Жертва-ускорение
                                 //if (mask[i]>0)
                                   //{
                                     resCastData[i]=p1*s1CastData[i]-s2CastData[i];
                                   //}
                       }
                    }
               //Конец Выполнение быстрой линейной комбинации
             }
           else if ( (p1==1) && (p2!=1) )
             {
                 //Выполнение быстрой линейной комбинации
                const TRnRelease1dData &s1Cast=dynamic_cast<const TRnRelease1dData &>(*s1);
                const TRnRelease1dData &s2Cast=dynamic_cast<const TRnRelease1dData &>(*s2);
                double *s1CastData=s1Cast.fData;
                double *s2CastData=s2Cast.fData;
                if (operation==true)
                   {
                     for (i=0;i<n;i++)
                       {
						       //Жертва-ускорение
                               //if (mask[i]>0)
                                  //{
                                    resCastData[i]=s1CastData[i]+p2*s2CastData[i];
                                  //}
                        }
                   }
                else
                   {
                     for (i=0;i<n;i++)
                       {
						       //Жертва-ускорение
                               //if (mask[i]>0)
                                   //{
                                     resCastData[i]=s1CastData[i]-p2*s2CastData[i];
                                   //}
                       }
                    }
               //Конец Выполнение быстрой линейной комбинации
             }
           else
             {
                 //Выполнение быстрой линейной комбинации
                const TRnRelease1dData &s1Cast=dynamic_cast<const TRnRelease1dData &>(*s1);
                const TRnRelease1dData &s2Cast=dynamic_cast<const TRnRelease1dData &>(*s2);
                double *s1CastData=s1Cast.fData;
                double *s2CastData=s2Cast.fData;
                if (operation==true)
                   {
                     for (i=0;i<n;i++)
                       {
						       //Жертва-ускорение
                               //if (mask[i]>0)
                                  //{
                                    resCastData[i]=s1CastData[i]+s2CastData[i];
                                  //}
                        }
                   }
                else
                   {
                     for (i=0;i<n;i++)
                       {
						       //Жертва-ускорение
                               //if (mask[i]>0)
                                   //{
                                     resCastData[i]=s1CastData[i]-s2CastData[i];
                                   //}
                       }
                    }
               //Конец Выполнение быстрой линейной комбинации
             }

         }  // Конец если variant==


   }

//*************************************
  void TRnRelease1dSpace::LocalReleaseAssignSpeedLinearCombination(const int id,const TLinearOperator *A1,
                                      const TRnData *s2,bool operation)
   {
      if (A1 == NULL)
       {
          throw "TRnRelease1dSpace::LocalReleaseAssignSpeedLinearCombination: First operator is NULL.";
       }

       //Попытка приведения параметра к TRnRelease1dLinearOperator
       try {dynamic_cast<const TRnRelease1dLinearOperator *>(A1);}
       catch (...) {throw "Error in TRnRelease1dData::LocalReleaseAssignSpeedLinearCombination: Could not cast parametr pased to TRnRelease1dLinearOperator * type.";}
       //Конец Попытка приведения параметра к TRnRelease1dLinearOperator


        if  (s2!=NULL)
        {
         //Попытка приведения параметра к TRnRelease1dData
          try {dynamic_cast<const TRnRelease1dData *>(s2);}
          catch (...) {throw "Error in TRnRelease1dData::LocalReleaseAssignSpeedLinearCombination: Could not cast parametr pased to TRnRelease1dData * type.";}
         //Конец Попытка приведения параметра к TRnRelease1dData
        }

       const T1dNumberMask & sMask = fCastModel->GetMask(id);
       const int *mask = sMask.Mask;
       int n = sMask.N;
       int i;

       double *resCastData=dynamic_cast<const TRnRelease1dData &>(*fMultiData[id]).fData;


       const TRnRelease1dLinearOperator &A1Cast=dynamic_cast<const TRnRelease1dLinearOperator &>(*A1);

       int actIndsCount = A1Cast.ActualIndexesCount;
       int *actInds = A1Cast.ActualIndexes;
       int loopVar;
       int indexRn;




       if (s2==NULL)
        {
           //Выполнение быстрой линейной комбинации
            if (actIndsCount>0)
                    {
                      for (loopVar=0;loopVar<actIndsCount;loopVar++)
                        {
                          indexRn = actInds[loopVar];
                          const TCorrespondence &cor = fCastModel->GetCorrespondence(indexRn);
                          const T1dCorrespondence &cor1d = dynamic_cast<const T1dCorrespondence &>(cor);
                          int ind = cor1d.Ind;
                          if (ind==id)
                           {
                             i = cor1d.I;
							 //Жертва-ускорение
                             //if (mask[i]>0)
                               //{
                                 resCastData[i]=A1Cast(id,i);
                               //}
                           } //Если индексы равны
                        } //Цикл по актуальным индексам
                    } //Если есть акт. индексы
            else
                 {
                     for (i=0;i<n;i++)
                       {
						       //Жертва-ускорение
                               //if (mask[i]>0)
                                   //{
                                     resCastData[i]=A1Cast(id,i);
                                   //}
                       }

                 }
               //Конец Выполнение быстрой линейной комбинации

        } //Конец нет второго операнда
      else
        {
           //Выполнение быстрой линейной комбинации
             const TRnRelease1dLinearOperator &A1Cast=dynamic_cast<const TRnRelease1dLinearOperator &>(*A1);
             const TRnRelease1dData &s2Cast=dynamic_cast<const TRnRelease1dData &>(*s2);
             double *s2CastData=s2Cast.fData;
             if (operation==true)
                 {
                   if (actIndsCount>0)
                    {
                      for (loopVar=0;loopVar<actIndsCount;loopVar++)
                        {
                          indexRn = actInds[loopVar];
                          const TCorrespondence &cor = fCastModel->GetCorrespondence(indexRn);
                          const T1dCorrespondence &cor1d = dynamic_cast<const T1dCorrespondence &>(cor);
                          int ind = cor1d.Ind;
                          if (ind==id)
                           {
                             i = cor1d.I;
							 //Жертва-ускорение
                             //if (mask[i]>0)
                               //{
                                 resCastData[i]=A1Cast(id,i)+s2CastData[i];
                               //}
                           } //Если индексы равны
                        } //Цикл по актуальным индексам
                    } //Если есть акт. индексы
                   else
                    {
                     for (i=0;i<n;i++)
                       {
						       //Жертва-ускорение
                               //if (mask[i]>0)
                                  //{
                                    resCastData[i]=A1Cast(id,i)+s2CastData[i];
                                  //}
                        }

                    }
                  }
                else
                   {
                    if (actIndsCount>0)
                    {
                      for (loopVar=0;loopVar<actIndsCount;loopVar++)
                        {
                          indexRn = actInds[loopVar];
                          const TCorrespondence &cor = fCastModel->GetCorrespondence(indexRn);
                          const T1dCorrespondence &cor1d = dynamic_cast<const T1dCorrespondence &>(cor);
                          int ind = cor1d.Ind;
                          if (ind==id)
                           {
                             i = cor1d.I;
							 //Жертва-ускорение
                             //if (mask[i]>0)
                               //{
                                 resCastData[i]=A1Cast(id,i)-s2CastData[i];
                               //}
                           } //Если индексы равны
                        } //Цикл по актуальным индексам
                    } //Если есть акт. индексы
                   else
                    {
                     for (i=0;i<n;i++)
                       {
						       //Жертва-ускорение
                               //if (mask[i]>0)
                                   //{
                                     resCastData[i]=A1Cast(id,i)-s2CastData[i];
                                   //}
                       }

                    }
                   }
               //Конец Выполнение быстрой линейной комбинации
        } //Конец есть второй операнд

   }
//*************************************
   double TRnRelease1dSpace::LocalReleaseNorma(const int id,int type)
   {
        double res=0;
        double *fCastData=dynamic_cast<const TRnRelease1dData &>(*fMultiData[id]).fData;
        //Выполнение вычисления нормы
        const T1dNumberMask & sMask = fCastModel->GetMask(id);
        const int *mask = sMask.Mask;
        int n = sMask.N;
        int i;

        double help;
        for (i=0;i<n;i++)
          {
			    //Жертва-ускорение
                //if (mask[i]>0)
                 //{
                  help=fCastData[i];
                  res=res+help*help;
                 //}
          }
        //Конец Выполнение вычисления нормы
        return res;
   }
//*************************************
   double TRnRelease1dSpace::LocalReleaseInnerproduct(const int id,const TRnData &s, int type)
   {

        //Попытка приведения параметра к TRnRelease1dData
         try {dynamic_cast<const TRnRelease1dData &>(s);}
         catch (...) {throw "Error in TRnRelease1dData::LocalReleaseInnerproduct: Could not cast parametr pased to TRnRelease1dData & type.";}
        //Конец Попытка приведения параметра к TRnRelease1dData

        const TRnRelease1dData &sCast=dynamic_cast<const TRnRelease1dData &>(s);

        double res=0;
        //Выполнение вычисления скалярного произведения
        const T1dNumberMask & sMask = fCastModel->GetMask(id);
        const int *mask = sMask.Mask;
        int n = sMask.N;
        int i;

        double *sCastData=sCast.fData;
        double *fCastData=dynamic_cast<const TRnRelease1dData &>(*fMultiData[id]).fData;

        int actIndsCount = ActualIndexesCount;
        int *actInds = ActualIndexes;
        int loopVar;
        int indexRn;


        if (actIndsCount>0)
         {
            for (loopVar=0;loopVar<actIndsCount;loopVar++)
                        {
                          indexRn = actInds[loopVar];
                          const TCorrespondence &cor = fCastModel->GetCorrespondence(indexRn);
                          const T1dCorrespondence &cor1d = dynamic_cast<const T1dCorrespondence &>(cor);
                          int ind = cor1d.Ind;
                          if (ind==id)
                           {
                             i = cor1d.I;
							 //Жертва-ускорение
                             //if (mask[i]>0)
                               //{
                                 res=res+fCastData[i]*sCastData[i];
                               //}
                           } //Если индексы равны
                        } //Цикл по актуальным индексам
         } //Если есть акт. индексы
        else
        {


        for (i=0;i<n;i++)
          {
			    //Жертва-ускорение
                //if (mask[i]>0)
                 //{
                  res=res+fCastData[i]*sCastData[i];
                 //}
          }


         }
        //Конец Выполнение вычисления скалярного произведения
        return res;
   }
//*************************************
 double TRnRelease1dSpace::LocalReleaseInnerproduct(const int id,const TLinearOperator *lA1,const TRnData *s2, int type)
  {
      const TLinearOperator *A1 = lA1;

      if (A1 == NULL)
       {
          throw "TRnRelease1dSpace::LocalReleaseInnerproduct: First operator is NULL.";
       }

       //Попытка приведения параметра к TRnRelease1dLinearOperator
       try {dynamic_cast<const TRnRelease1dLinearOperator *>(A1);}
       catch (...) {throw "Error in TRnRelease1dSpace::LocalReleaseInnerproduct: Could not cast parametr pased to TRnRelease1dLinearOperator * type.";}
       //Конец Попытка приведения параметра к TRnRelease1dLinearOperator


        if  (s2!=NULL)
        {
         //Попытка приведения параметра к TRnRelease1dData
          try {dynamic_cast<const TRnRelease1dData *>(s2);}
          catch (...) {throw "TRnRelease1dSpace::LocalReleaseInnerproduct: Could not cast parametr pased to TRnRelease1dData * type.";}
         //Конец Попытка приведения параметра к TRnRelease1dData
        }

       double res=0;

       const T1dNumberMask & sMask = fCastModel->GetMask(id);
       const int *mask = sMask.Mask;
       int n = sMask.N;
       int i;

       double *fCastData=dynamic_cast<const TRnRelease1dData &>(*fMultiData[id]).fData;

       if (s2==NULL)
        {
           //Выполнение быстрой линейной комбинации
             const TRnRelease1dLinearOperator &A1Cast=dynamic_cast<const TRnRelease1dLinearOperator &>(*A1);
             for (i=0;i<n;i++)
              {
				      //Жертва-ускорение
                      //if (mask[i]>0)
                       //{
                          res=res+fCastData[i]*A1Cast(id,i);
                       //}

              }
           //Конец Выполнение быстрой линейной комбинации

        } //Конец нет второго операнда

      else
        {
           throw "TRnRelease1dSpace::LocalReleaseInnerproduct: The function is not implemented for non-NULL second operand.";
        } //Конец есть второй операнд

      return res;
  }
//*************************************

   double &TRnRelease1dSpace::GetByRnIndex(int i) const
    {
       const T1dCorrespondence &c = dynamic_cast<const T1dCorrespondence &>(fCastModel->GetCorrespondence(i));
       double *fCastData=dynamic_cast<const TRnRelease1dData &>(*fMultiData[c.Ind]).fData;
       return fCastData[c.I];
    }
//*************************************
   double TRnRelease1dSpace::LocalGetSpecificElement(const int id,int type,int &elementIndex) const
    {
       double res=0;

       double *fCastData=dynamic_cast<const TRnRelease1dData &>(*fMultiData[id]).fData;
        //Выполнение вычисления
        const T1dNumberMask & sMask = fCastModel->GetMask(id);
        const int *mask = sMask.Mask;
        int n = sMask.N;
        int i;

        double localRes;
        int localIndex = 1;
        bool assignLocalIndex;
        const T1dInverseCorrespondence &invCor = dynamic_cast<const T1dInverseCorrespondence &>(fCastModel->GetInverseCorrespondence());
        for (i=0;i<n;i++)
          {
			    //Жертва-ускорение
                //if (mask[i]>0)
                 //{
                   assignLocalIndex = false;
                   localRes=fabs(fCastData[i]);
                   if (type == 1)
                    {
                      if (localRes>res) {res = localRes; }
                    }
                   else if (type == 2)
                    {
                      if (localRes<res) {res = localRes; }
                    }
                   else
                    {
                       throw "Error in TRnRelease1dSpace::GetSpecificElement: The fuction is not implemented for this type.";
                    }
                   if (assignLocalIndex==true)
                    {
                      localIndex = invCor[id][i];
                    }
                 //}
          }
        //Конец Выполнение вычисления
        elementIndex = localIndex;
        return res;
    }
//*************************************

//#################################################################



//---------------------Rn1 (End)------------------------------





//---------------------Rn3 (Begin)------------------------------



//#################################################################

   //*************************************
       TRnRelease3dData::TRnRelease3dData(const int lN,const int lL,const int lM): TRnData(),fN(lN), fL(lL),fM(lM), Data(this)
        {
           //Выделение памяти
             int n=fN;
             int l=fL;
             int m=fM;
             int i,j,k;
             fData=new double **[n];
             for (i=0;i<n;i++)
               {
                 fData[i]=new double *[l];
                 for (j=0;j<l;j++)
                  {
                    fData[i][j]=new double [m];
                  }
               }
           //Конец Выделение памяти

           //Обнуление всех элементов (Начало)

             for (i=0;i<n;i++)
               {
                 for (j=0; j<l; j++)
                  {
                    for (k=0;k<m;k++)
                     {
                       fData[i][j][k]=0;
                     }  
                  }
               }
           //Обнуление всех элементов (Конец)
        }
   //*************************************
       TRnRelease3dData::~TRnRelease3dData()
        {
          int i; int n=fN;
          int j; int l=fL;
          for (i=0;i<n;i++)
            {
              for (j=0;j<l;j++)
                {
                  delete[] fData[i][j];
                }
              delete[] fData[i];
            }
          delete[] fData;
          fData=NULL;
        }
   //*************************************

       void TRnRelease3dData::HelpVirtualFunction() const
        {

        }
   //*************************************

//#################################################################

//*************************************
  TRnRelease3dSpace::TRnRelease3dSpace(const int lModelName): TRnReleaseMultiDataSpace(lModelName) //Конструктор
   {
      //Приведение модели к TRnRelease3dSpaceModel
      try {fCastModel=dynamic_cast<TRnRelease3dSpaceModel *>(fModel);}
      catch (...) {throw "Error in TRnRelease3dSpace constructor: Could not cast model to TRnRelease3dSpaceModel type.";}
      //Конец Приведение модели к TRnRelease3dSpaceModel

      int dataCount = fCastModel->DataCount;
      int i;
      for (i=0;i<dataCount;i++)
       {
          int n = (fCastModel->GetMask(i)).N;
          int l = (fCastModel->GetMask(i)).L;
          int m = (fCastModel->GetMask(i)).M;

          fMultiData[i] = new TRnRelease3dData(n,l,m);
       }
   }
//*************************************
  TRnRelease3dSpace::~TRnRelease3dSpace() //Деструктор
   {
      int dataCount = fCastModel->DataCount;
      int i;
      for (i=0;i<dataCount;i++)
       {
          delete fMultiData[i];// fMultiData[i] = NULL;
       }
   }
//*************************************
  double **TRnRelease3dSpace::operator[](int i)
   {
     return dynamic_cast<TRnRelease3dData *>(fMultiData[fActiveDataIndex])->Data[i];
   }
//*************************************
  TRnSpace &TRnRelease3dSpace::NewItem()
   {
     return *(new TRnRelease3dSpace(fModelName));
   }
//*************************************
   void TRnRelease3dSpace::LocalReleaseOperatorEqual(const int id,const TRnData &s,bool isCopying)
   {

        //Попытка приведения параметра к TRnRelease3dData
         try {dynamic_cast<const TRnRelease3dData &>(s);}
         catch (...) {throw "Error in TRnRelease3dData::LocalReleaseOperatorEqual: Could not cast parametr pased to TRnRelease3dData & type.";}
        //Конец Попытка приведения параметра к TRnRelease3dData

        const TRnRelease3dData &sCast=dynamic_cast<const TRnRelease3dData &>(s);

        //Выполнение присваивания
        const T3dNumberMask & sMask = fCastModel->GetMask(id);
        const int * const * const * mask = sMask.Mask;
        int n = sMask.N;
        int l = sMask.L;
        int m = sMask.M;
        int i; int j; int k;

        double ***sCastData=sCast.fData;
        double ***fCastData=dynamic_cast<const TRnRelease3dData &>(*fMultiData[id]).fData;

        
        
        double **helpData;
        double *helpDataF;
        double **helpCastData;
        double *helpCastDataF;

        if (isCopying==false)
         {
           for (i=0;i<n;i++)
            {
              const int * const * helpMask = mask[i];
              helpData=fCastData[i];
              helpCastData=sCastData[i];

              for (j=0;j<l;j++)
                {
                   const int * helpMaskF = helpMask[j];
                   helpDataF=helpData[j];
                   helpCastDataF=helpCastData[j];

                   for (k=0;k<m;k++)
                    {
                     if (helpMaskF[k]>0) {helpDataF[k]=helpCastDataF[k];}
                    }
                }

            }
         }
        else
         {
           for (i=0;i<n;i++)
            {
              helpData=fCastData[i];
              helpCastData=sCastData[i];

              for (j=0;j<l;j++)
                {
                   helpDataF=helpData[j];
                   helpCastDataF=helpCastData[j];
                   for (k=0;k<m;k++)
                    {
                      helpDataF[k]=helpCastDataF[k];
                    }
                }

            }
         }
        //Конец Выполнение присваивания
   }
//*************************************
   double TRnRelease3dSpace::LocalReleaseMetric(const int id,const TRnData &s,int type)
   {
        //Попытка приведения параметра к TRnRelease3dData
         try {dynamic_cast<const TRnRelease3dData &>(s);}
         catch (...) {throw "Error in TRnRelease3dData::LocalReleaseMetric: Could not cast parametr pased to TRnRelease3dData & type.";}
        //Конец Попытка приведения параметра к TRnRelease3dData

        const TRnRelease3dData &sCast=dynamic_cast<const TRnRelease3dData &>(s);

        double res=0;
        //Выполнение вычисления метрики
        const T3dNumberMask & sMask = fCastModel->GetMask(id);
        const int * const * const * mask = sMask.Mask;
        int n = sMask.N;
        int l = sMask.L;
        int m = sMask.M;
        int i; int j; int k;

        double ***sCastData=sCast.fData;
        double ***fCastData=dynamic_cast<const TRnRelease3dData &>(*fMultiData[id]).fData;


        
        
        double **helpData;
        double *helpDataF;
        double **helpCastData;
        double *helpCastDataF;



        double help;

        for (i=0;i<n;i++)
          {
            const int * const * helpMask = mask[i];
            helpData=fCastData[i];
            helpCastData=sCastData[i];

            for (j=0;j<l;j++)
              {
                   const int * helpMaskF = helpMask[j];
                   helpDataF=helpData[j];
                   helpCastDataF=helpCastData[j];

                   for (k=0;k<m;k++)
                    {
                       if (helpMaskF[k]>0)
                          {
                            help=helpDataF[k]-helpCastDataF[k];
                            res=res+help*help;
                          }
                 }
              }
          }
        //Конец Выполнение вычисления метрики
        return res;
   }
//*************************************
   void TRnRelease3dSpace::LocalToZero(const int id)
   {
        const T3dNumberMask & sMask = fCastModel->GetMask(id);
        const int * const * const * mask = sMask.Mask;
        int n = sMask.N;
        int l = sMask.L;
        int m = sMask.M;
        int i; int j; int k;

        double ***fCastData=dynamic_cast<const TRnRelease3dData &>(*fMultiData[id]).fData;

        

        double **helpData;
        double *helpDataF;


        int actIndsCount = ActualIndexesCount;
        int *actInds = ActualIndexes;
        int loopVar;
        int indexRn;

        if (actIndsCount>0)
         {
            for (loopVar=0;loopVar<actIndsCount;loopVar++)
                        {
                          indexRn = actInds[loopVar];
                          const TCorrespondence &cor = fCastModel->GetCorrespondence(indexRn);
                          const T3dCorrespondence &cor3d = dynamic_cast<const T3dCorrespondence &>(cor);
                          int ind = cor3d.Ind;
                          if (ind==id)
                           {
                             i = cor3d.I;
                             j = cor3d.J;
                             k = cor3d.K;
                             if (mask[i][j][k]>0)
                               {
                                 fCastData[i][j][k]=0;
                               }
                           } //Если индексы равны
                        } //Цикл по актуальным индексам
                    } //Если есть акт. индексы
        else
        {


           for (i=0;i<n;i++)
            {
              const int * const * helpMask = mask[i];
              helpData=fCastData[i];
              for (j=0;j<l;j++)
               {
                 const int * helpMaskF = helpMask[j];
                 helpDataF=helpData[j];

                 for (k=0;k<m;k++)
                    {
                       if (helpMaskF[k]>0) {helpDataF[k]=0;}
                    }
               }
            }
        } //Если есть акт. индексы

   }
//*************************************
   void TRnRelease3dSpace::LocalReleaseAssignSpeedLinearCombination(const int id,double p1,const TRnData *s1,
                                      double p2,const TRnData *s2,bool operation)
   {
       //Может быть несколько лишне организована проверка на то, что
       //один из сомножителей =1 - может это делает компилятор (хотя как ?)


        int variant; //-1 не использовать операнды
                     //0 - оба операнда
                     //1 - использовать первый
                     //2 - использовать второй



        if  (s1!=NULL)
        {
        //Попытка приведения параметра к TRnRelease3dData
         try {dynamic_cast<const TRnRelease3dData *>(s1);}
         catch (...) {throw "Error in TRnRelease3dData::LocalReleaseAssignSpeedLinearCombination: Could not cast parametr pased to TRnRelease3dData * type.";}
        //Конец Попытка приведения параметра к TRnRelease3dData
        }

        if  (s2!=NULL)
        {
        //Попытка приведения параметра к TRnRelease3dData
         try {dynamic_cast<const TRnRelease3dData *>(s2);}
         catch (...) {throw "Error in TRnRelease3dData::LocalReleaseAssignSpeedLinearCombination: Could not cast parametr pased to TRnRelease3dData * type.";}
        //Конец Попытка приведения параметра к TRnRelease3dData
        }


       //Подстраховка ...

         if ( ( (p1!=0)&&(s1==NULL) ) || ( (p2!=0)&&(s2==NULL) ) )
          {
            throw "Error in TRnRelease3dData::LocalReleaseAssignSpeedLinearCombination: Incorrect combination: p<>0 but s=NULL.";
          }

       //Конец Подстраховка ...



       //Определение варианта
         if (p1!=0)
          {
            if (p2!=0) {variant=0;}
            else {variant=1;}
          }
         else
          {
             if (p2!=0) {variant=2;}
             else {variant=-1;}

          }
       //Конец Определение варианта



       const T3dNumberMask & sMask = fCastModel->GetMask(id);
       const int * const * const * mask = sMask.Mask;
       int n = sMask.N;
       int l = sMask.L;
       int m = sMask.M;
       int i; int j; int k;


       double ***resCastData=dynamic_cast<const TRnRelease3dData &>(*fMultiData[id]).fData;

       
       
       double **help1CastData;
       double *help1CastDataF;
       double **help2CastData;
       double *help2CastDataF;
       double **helpResCastData;
       double *helpResCastDataF;



       if (variant==-1) {LocalToZero(id); return;}
       else if (variant==1)
         {
           if (p1!=1)
             {
               //Выполнение быстрой линейной комбинации
                const TRnRelease3dData &s1Cast=dynamic_cast<const TRnRelease3dData &>(*s1);
                double ***s1CastData=s1Cast.fData;
                for (i=0;i<n;i++)
                       {
                          const int * const * helpMask = mask[i];
                          help1CastData=s1CastData[i];
                          helpResCastData=resCastData[i];

                          for (j=0;j<l;j++)
                            {
                               const int * helpMaskF = helpMask[j];
                               help1CastDataF=help1CastData[j];
                               helpResCastDataF=helpResCastData[j];

                               for (k=0;k<m;k++)
                                  {
                                      if (helpMaskF[k]>0)
                                         {
                                           helpResCastDataF[k]=p1*help1CastDataF[k];
                                         }
                                  }
                            }

                        }
                //Конец Выполнение быстрой линейной комбинации
             }
           else
             {
                 //Выполнение быстрой линейной комбинации
                const TRnRelease3dData &s1Cast=dynamic_cast<const TRnRelease3dData &>(*s1);
                double ***s1CastData=s1Cast.fData;
                for (i=0;i<n;i++)
                       {
                          const int * const * helpMask = mask[i];
                          help1CastData=s1CastData[i];
                          helpResCastData=resCastData[i];

                          for (j=0;j<l;j++)
                            {
                              const int * helpMaskF = helpMask[j];
                              help1CastDataF=help1CastData[j];
                              helpResCastDataF=helpResCastData[j];

                               for (k=0;k<m;k++)
                                  {
                                      if (helpMaskF[k]>0)
                                         {
                                            helpResCastDataF[k]=help1CastDataF[k];
                                          }
                                  }
                            }

                        }
                //Конец Выполнение быстрой линейной комбинации
             }
         }
       else if (variant==2)
         {
            if (p2!=1)
             {
                //Выполнение быстрой линейной комбинации
                const TRnRelease3dData &s2Cast=dynamic_cast<const TRnRelease3dData &>(*s2);
                double ***s2CastData=s2Cast.fData;
                if (operation==true)
                   {
                     for (i=0;i<n;i++)
                       {
                          const int * const * helpMask = mask[i];
                          help2CastData=s2CastData[i];
                          helpResCastData=resCastData[i];

                          for (j=0;j<l;j++)
                            {
                               const int * helpMaskF = helpMask[j];
                               help2CastDataF=help2CastData[j];
                               helpResCastDataF=helpResCastData[j];

                               for (k=0;k<m;k++)
                                  {
                                    if (helpMaskF[k]>0)
                                     {
                                       helpResCastDataF[k]=p2*help2CastDataF[k];
                                     }
                                  }
                            }

                        }
                   }
                else
                   {
                     for (i=0;i<n;i++)
                       {
                         const int * const * helpMask = mask[i];
                         help2CastData=s2CastData[i];
                         helpResCastData=resCastData[i];

                         for (j=0;j<l;j++)
                            {
                               const int * helpMaskF = helpMask[j];
                               help2CastDataF=help2CastData[j];
                               helpResCastDataF=helpResCastData[j];

                               for (k=0;k<m;k++)
                                  {
                                      if (helpMaskF[k]>0)
                                        {
                                           helpResCastDataF[k]=-p2*help2CastDataF[k];
                                        }
                                  }
                             }

                       }
                    }
               //Конец Выполнение быстрой линейной комбинации

             }
           else
             {
                 //Выполнение быстрой линейной комбинации
                const TRnRelease3dData &s2Cast=dynamic_cast<const TRnRelease3dData &>(*s2);
                double ***s2CastData=s2Cast.fData;
                if (operation==true)
                   {
                     for (i=0;i<n;i++)
                       {
                          const int * const * helpMask = mask[i];
                          help2CastData=s2CastData[i];
                          helpResCastData=resCastData[i];

                          for (j=0;j<l;j++)
                            {
                               const int * helpMaskF = helpMask[j];
                               help2CastDataF=help2CastData[j];
                               helpResCastDataF=helpResCastData[j];

                               for (k=0;k<m;k++)
                                  {
                                     if (helpMaskF[k]>0)
                                       {
                                         helpResCastDataF[k]=help2CastDataF[k];
                                       }
                                  }
                            }

                        }
                   }
                else
                   {
                     for (i=0;i<n;i++)
                       {
                         const int * const * helpMask = mask[i];
                         help2CastData=s2CastData[i];
                         helpResCastData=resCastData[i];

                         for (j=0;j<l;j++)
                            {
                               const int * helpMaskF = helpMask[j];
                               help2CastDataF=help2CastData[j];
                               helpResCastDataF=helpResCastData[j];

                               for (k=0;k<m;k++)
                                  {
                                     if (helpMaskF[k]>0)
                                       {
                                          helpResCastDataF[k]=-help2CastDataF[k];
                                       }
                                  }
                             }

                       }
                    }
               //Конец Выполнение быстрой линейной комбинации
             }
         }
       else
         {

           if ( (p1!=1) && (p2!=1) )
             {

                //Выполнение быстрой линейной комбинации
                const TRnRelease3dData &s1Cast=dynamic_cast<const TRnRelease3dData &>(*s1);
                const TRnRelease3dData &s2Cast=dynamic_cast<const TRnRelease3dData &>(*s2);
                double ***s1CastData=s1Cast.fData;
                double ***s2CastData=s2Cast.fData;
                if (operation==true)
                   {
                     for (i=0;i<n;i++)
                       {
                          const int * const * helpMask = mask[i];
                          help1CastData=s1CastData[i];
                          help2CastData=s2CastData[i];
                          helpResCastData=resCastData[i];

                          for (j=0;j<l;j++)
                            {
                               const int * helpMaskF = helpMask[j];
                               help1CastDataF=help1CastData[j];
                               help2CastDataF=help2CastData[j];
                               helpResCastDataF=helpResCastData[j];

                               for (k=0;k<m;k++)
                                  {
                                    if (helpMaskF[k]>0)
                                      {
                                        helpResCastDataF[k]=p1*help1CastDataF[k]+p2*help2CastDataF[k];
                                      }
                                  }
                            }

                        }
                   }
                else
                   {
                     for (i=0;i<n;i++)
                       {

                         const int * const * helpMask = mask[i];
                         help1CastData=s1CastData[i];
                         help2CastData=s2CastData[i];
                         helpResCastData=resCastData[i];

                         for (j=0;j<l;j++)
                            {
                               const int * helpMaskF = helpMask[j];
                               help1CastDataF=help1CastData[j];
                               help2CastDataF=help2CastData[j];
                               helpResCastDataF=helpResCastData[j];

                               for (k=0;k<m;k++)
                                  {
                                    if (helpMaskF[k]>0)
                                      {
                                        helpResCastDataF[k]=p1*help1CastDataF[k]-p2*help2CastDataF[k];
                                      }
                                  }
                             }

                       }
                    }
               //Конец Выполнение быстрой линейной комбинации

             }
           else if ( (p1!=1) && (p2==1) )
             {
                //Выполнение быстрой линейной комбинации
                const TRnRelease3dData &s1Cast=dynamic_cast<const TRnRelease3dData &>(*s1);
                const TRnRelease3dData &s2Cast=dynamic_cast<const TRnRelease3dData &>(*s2);
                double ***s1CastData=s1Cast.fData;
                double ***s2CastData=s2Cast.fData;
                if (operation==true)
                   {
                     for (i=0;i<n;i++)
                       {
                          const int * const * helpMask = mask[i];
                          help1CastData=s1CastData[i];
                          help2CastData=s2CastData[i];
                          helpResCastData=resCastData[i];

                          for (j=0;j<l;j++)
                            {
                               const int * helpMaskF = helpMask[j];
                               help1CastDataF=help1CastData[j];
                               help2CastDataF=help2CastData[j];
                               helpResCastDataF=helpResCastData[j];

                               for (k=0;k<m;k++)
                                  {
                                    if (helpMaskF[k]>0)
                                      {
                                         helpResCastDataF[k]=p1*help1CastDataF[k]+help2CastDataF[k];
                                      }
                                  }
                            }

                        }
                   }
                else
                   {
                     for (i=0;i<n;i++)
                       {
                         const int * const * helpMask = mask[i];
                         help1CastData=s1CastData[i];
                         help2CastData=s2CastData[i];
                         helpResCastData=resCastData[i];

                         for (j=0;j<l;j++)
                            {
                               const int * helpMaskF = helpMask[j];
                               help1CastDataF=help1CastData[j];
                               help2CastDataF=help2CastData[j];
                               helpResCastDataF=helpResCastData[j];

                               for (k=0;k<m;k++)
                                  {
                                      if (helpMaskF[k]>0)
                                         {
                                           helpResCastDataF[k]=p1*help1CastDataF[k]-help2CastDataF[k];
                                         }
                                  }
                             }

                       }
                    }
               //Конец Выполнение быстрой линейной комбинации
             }
           else if ( (p1==1) && (p2!=1) )
             {
                 //Выполнение быстрой линейной комбинации
                const TRnRelease3dData &s1Cast=dynamic_cast<const TRnRelease3dData &>(*s1);
                const TRnRelease3dData &s2Cast=dynamic_cast<const TRnRelease3dData &>(*s2);
                double ***s1CastData=s1Cast.fData;
                double ***s2CastData=s2Cast.fData;
                if (operation==true)
                   {
                     for (i=0;i<n;i++)
                       {
                          const int * const * helpMask = mask[i];
                          help1CastData=s1CastData[i];
                          help2CastData=s2CastData[i];
                          helpResCastData=resCastData[i];

                          for (j=0;j<l;j++)
                            {
                               const int * helpMaskF = helpMask[j];
                               help1CastDataF=help1CastData[j];
                               help2CastDataF=help2CastData[j];
                               helpResCastDataF=helpResCastData[j];

                               for (k=0;k<m;k++)
                                  {
                                     if (helpMaskF[k]>0)
                                       {
                                          helpResCastDataF[k]=help1CastDataF[k]+p2*help2CastDataF[k];
                                       }
                                  }
                            }

                        }
                   }
                else
                   {
                     for (i=0;i<n;i++)
                       {
                         const int * const * helpMask = mask[i];
                         help1CastData=s1CastData[i];
                         help2CastData=s2CastData[i];
                         helpResCastData=resCastData[i];

                         for (j=0;j<l;j++)
                            {
                               const int * helpMaskF = helpMask[j];
                               help1CastDataF=help1CastData[j];
                               help2CastDataF=help2CastData[j];
                               helpResCastDataF=helpResCastData[j];

                               for (k=0;k<m;k++)
                                  {
                                     if (helpMaskF[k]>0)
                                      {
                                        helpResCastDataF[k]=help1CastDataF[k]-p2*help2CastDataF[k];
                                      }
                                  }
                             }

                       }
                    }
               //Конец Выполнение быстрой линейной комбинации
             }
           else
             {
                 //Выполнение быстрой линейной комбинации
                const TRnRelease3dData &s1Cast=dynamic_cast<const TRnRelease3dData &>(*s1);
                const TRnRelease3dData &s2Cast=dynamic_cast<const TRnRelease3dData &>(*s2);
                double ***s1CastData=s1Cast.fData;
                double ***s2CastData=s2Cast.fData;
                if (operation==true)
                   {
                     for (i=0;i<n;i++)
                       {
                          const int * const * helpMask = mask[i];
                          help1CastData=s1CastData[i];
                          help2CastData=s2CastData[i];
                          helpResCastData=resCastData[i];

                          for (j=0;j<l;j++)
                            {
                               const int * helpMaskF = helpMask[j];
                               help1CastDataF=help1CastData[j];
                               help2CastDataF=help2CastData[j];
                               helpResCastDataF=helpResCastData[j];

                               for (k=0;k<m;k++)
                                  {
                                      if (helpMaskF[k]>0)
                                        {
                                            helpResCastDataF[k]=help1CastDataF[k]+help2CastDataF[k];
                                        }
                                  }
                            }

                        }
                   }
                else
                   {
                     for (i=0;i<n;i++)
                       {
                         const int * const * helpMask = mask[i];
                         help1CastData=s1CastData[i];
                         help2CastData=s2CastData[i];
                         helpResCastData=resCastData[i];

                         for (j=0;j<l;j++)
                            {
                               const int * helpMaskF = helpMask[j];
                               help1CastDataF=help1CastData[j];
                               help2CastDataF=help2CastData[j];
                               helpResCastDataF=helpResCastData[j];

                               for (k=0;k<m;k++)
                                  {
                                     if (helpMaskF[k]>0)
                                       {
                                          helpResCastDataF[k]=help1CastDataF[k]-help2CastDataF[k];
                                       }
                                  }
                             }

                       }
                    }
               //Конец Выполнение быстрой линейной комбинации
             }

         }  // Конец если variant==


   }

//*************************************
  void TRnRelease3dSpace::LocalReleaseAssignSpeedLinearCombination(const int id,const TLinearOperator *A1,
                                      const TRnData *s2,bool operation)
   {
      if (A1 == NULL)
       {
          throw "TRnRelease3dSpace::LocalReleaseAssignSpeedLinearCombination: First operator is NULL.";
       }

       //Попытка приведения параметра к TRnRelease3dLinearOperator
       try {dynamic_cast<const TRnRelease3dLinearOperator *>(A1);}
       catch (...) {throw "Error in TRnRelease3dData::LocalReleaseAssignSpeedLinearCombination: Could not cast parametr pased to TRnRelease3dLinearOperator * type.";}
       //Конец Попытка приведения параметра к TRnRelease3dLinearOperator


        if  (s2!=NULL)
        {
         //Попытка приведения параметра к TRnRelease3dData
          try {dynamic_cast<const TRnRelease3dData *>(s2);}
          catch (...) {throw "Error in TRnRelease3dData::LocalReleaseAssignSpeedLinearCombination: Could not cast parametr pased to TRnRelease3dData * type.";}
         //Конец Попытка приведения параметра к TRnRelease3dData
        }

       const T3dNumberMask & sMask = fCastModel->GetMask(id);
       const int * const * const * mask = sMask.Mask;
       int n = sMask.N;
       int l = sMask.L;
       int m = sMask.M;
       int i; int j; int k;


       double ***resCastData=dynamic_cast<const TRnRelease3dData &>(*fMultiData[id]).fData;

       
       
       //double **help1CastData;
       //double *help1CastDataF;
       double **help2CastData;
       double *help2CastDataF;
       double **helpResCastData;
       double *helpResCastDataF;


      const TRnRelease3dLinearOperator &A1Cast=dynamic_cast<const TRnRelease3dLinearOperator &>(*A1);
      int actIndsCount = A1Cast.ActualIndexesCount;
      int *actInds =A1Cast.ActualIndexes;
      int loopVar;
      int indexRn;


      if (s2==NULL)
        {
           //Выполнение быстрой линейной комбинации
                   if (actIndsCount>0)
                    {
                      for (loopVar=0;loopVar<actIndsCount;loopVar++)
                        {
                          indexRn = actInds[loopVar];
                          const TCorrespondence &cor = fCastModel->GetCorrespondence(indexRn);
                          const T3dCorrespondence &cor3d = dynamic_cast<const T3dCorrespondence &>(cor);
                          int ind = cor3d.Ind;
                          if (ind==id)
                           {
                             i = cor3d.I;
                             j = cor3d.J;
                             k = cor3d.K;
                             if (mask[i][j][k]>0)
                               {
                                 resCastData[i][j][k]=A1Cast(id,i,j,k);
                               }
                           } //Если индексы равны
                        } //Цикл по актуальным индексам
                    } //Если есть акт. индексы
                   else
                    {
                     for (i=0;i<n;i++)
                       {
                          const int * const * helpMask = mask[i];
                          helpResCastData=resCastData[i];

                          for (j=0;j<l;j++)
                            {
                               const int * helpMaskF = helpMask[j];
                               helpResCastDataF=helpResCastData[j];

                               for (k=0;k<m;k++)
                                  {
                                     if (helpMaskF[k]>0)
                                       {
                                         helpResCastDataF[k]=A1Cast(id,i,j,k);
                                       }
                                  }
                            }

                        } //Цикл по первому измерению

                    } //Если есть акт. индексы

               //Конец Выполнение быстрой линейной комбинации

        } //Конец нет второго операнда
      else
        {
           //Выполнение быстрой линейной комбинации
             const TRnRelease3dData &s2Cast=dynamic_cast<const TRnRelease3dData &>(*s2);
             double ***s2CastData=s2Cast.fData;
             if (operation==true)
                 {
                   if (actIndsCount>0)
                    {
                      for (loopVar=0;loopVar<actIndsCount;loopVar++)
                        {
                          indexRn = actInds[loopVar];
                          const TCorrespondence &cor = fCastModel->GetCorrespondence(indexRn);
                          const T3dCorrespondence &cor3d = dynamic_cast<const T3dCorrespondence &>(cor);
                          int ind = cor3d.Ind;
                          if (ind==id)
                           {
                             i = cor3d.I;
                             j = cor3d.J;
                             k = cor3d.K;
                             if (mask[i][j][k]>0)
                               {
                                 resCastData[i][j][k]=A1Cast(id,i,j,k)+s2CastData[i][j][k];
                               }
                           } //Если индексы равны
                        } //Цикл по актуальным индексам
                    } //Если есть акт. индексы
                   else
                    {
                     for (i=0;i<n;i++)
                       {
                          const int * const * helpMask = mask[i];
                          help2CastData=s2CastData[i];
                          helpResCastData=resCastData[i];

                          for (j=0;j<l;j++)
                            {
                               const int * helpMaskF = helpMask[j];
                               help2CastDataF=help2CastData[j];
                               helpResCastDataF=helpResCastData[j];

                               for (k=0;k<m;k++)
                                  {
                                     if (helpMaskF[k]>0)
                                      {
                                        helpResCastDataF[k]=A1Cast(id,i,j,k)+help2CastDataF[k];
                                      }
                                  }
                            }
                        } //Цикл по первому измерению
                    } //Если есть акт. индексы
                  }//если суммировать
                else
                   {
                    if (actIndsCount>0)
                    {
                      for (loopVar=0;loopVar<actIndsCount;loopVar++)
                        {
                          indexRn = actInds[loopVar];
                          const TCorrespondence &cor = fCastModel->GetCorrespondence(indexRn);
                          const T3dCorrespondence &cor3d = dynamic_cast<const T3dCorrespondence &>(cor);
                          int ind = cor3d.Ind;
                          if (ind==id)
                           {
                             i = cor3d.I;
                             j = cor3d.J;
                             k = cor3d.K;
                             if (mask[i][j][k]>0)
                               {
                                 resCastData[i][j][k]=A1Cast(id,i,j,k) - s2CastData[i][j][k];
                               }
                           } //Если индексы равны
                        } //Цикл по актуальным индексам
                    } //Если есть акт. индексы
                   else
                    {
                     for (i=0;i<n;i++)
                       {
                         const int * const * helpMask = mask[i];
                         help2CastData=s2CastData[i];
                         helpResCastData=resCastData[i];

                         for (j=0;j<l;j++)
                            {
                               const int * helpMaskF = helpMask[j];
                               help2CastDataF=help2CastData[j];
                               helpResCastDataF=helpResCastData[j];

                               for (k=0;k<m;k++)
                                  {
                                     if (helpMaskF[k]>0)
                                      {
                                         helpResCastDataF[k]=A1Cast(id,i,j,k)-help2CastDataF[k];
                                      }
                                  }
                             }

                       } //Цикл по первому измерению
                    } //Если есть акт. индексы

                   }//если ссумировать
               //Конец Выполнение быстрой линейной комбинации
        } //Конец есть второй операнд

   }
//*************************************
   double TRnRelease3dSpace::LocalReleaseNorma(const int id,int type)
   {
        double res=0;
        double ***fCastData=dynamic_cast<const TRnRelease3dData &>(*fMultiData[id]).fData;
        //Выполнение вычисления нормы
        const T3dNumberMask & sMask = fCastModel->GetMask(id);
        const int * const * const * mask = sMask.Mask;
        int n = sMask.N;
        int l = sMask.L;
        int m = sMask.M;
        int i; int j; int k;

        
        
        double **helpData;
        double *helpDataF;

        double help;
        for (i=0;i<n;i++)
          {
            const int * const * helpMask = mask[i];
            helpData=fCastData[i];

            for (j=0;j<l;j++)
             {
                const int * helpMaskF = helpMask[j];
                helpDataF=helpData[j];

                for (k=0;k<m;k++)
                  {
                      if (helpMaskF[k]>0)
                        {
                           help=helpDataF[k];
                           res=res+help*help;
                        }
                  }
             }
          }
        //Конец Выполнение вычисления нормы
        return res;
   }
//*************************************
   double TRnRelease3dSpace::LocalReleaseInnerproduct(const int id,const TRnData &s, int type)
   {

        //Попытка приведения параметра к TRnRelease3dData
         try {dynamic_cast<const TRnRelease3dData &>(s);}
         catch (...) {throw "Error in TRnRelease3dData::LocalReleaseInnerproduct: Could not cast parametr pased to TRnRelease3dData & type.";}
        //Конец Попытка приведения параметра к TRnRelease3dData

        const TRnRelease3dData &sCast=dynamic_cast<const TRnRelease3dData &>(s);

        double res=0;
        //Выполнение вычисления скалярного произведения
        const T3dNumberMask & sMask = fCastModel->GetMask(id);
        const int * const * const * mask = sMask.Mask;
        int n = sMask.N;
        int l = sMask.L;
        int m = sMask.M;
        int i; int j; int k;


        double ***sCastData=sCast.fData;
        double ***fCastData=dynamic_cast<const TRnRelease3dData &>(*fMultiData[id]).fData;

        
        
        double **helpCastData;
        double *helpCastDataF;
        double **helpData;
        double *helpDataF;

        int actIndsCount = ActualIndexesCount;
        int *actInds = ActualIndexes;
        int loopVar;
        int indexRn;


        if (actIndsCount>0)
         {
            for (loopVar=0;loopVar<actIndsCount;loopVar++)
                        {
                          indexRn = actInds[loopVar];
                          const TCorrespondence &cor = fCastModel->GetCorrespondence(indexRn);
                          const T3dCorrespondence &cor3d = dynamic_cast<const T3dCorrespondence &>(cor);
                          int ind = cor3d.Ind;
                          if (ind==id)
                           {
                             i = cor3d.I;
                             j = cor3d.J;
                             k = cor3d.K;
                             if (mask[i][j][k]>0)
                               {
                                 res=res+fCastData[i][j][k]*sCastData[i][j][k];
                               }
                           } //Если индексы равны
                        } //Цикл по актуальным индексам
                    } //Если есть акт. индексы
        else
        {
         for (i=0;i<n;i++)
          {
            const int * const * helpMask = mask[i];
            helpData=fCastData[i];
            helpCastData=sCastData[i];

            for (j=0;j<l;j++)
             {
                const int * helpMaskF = helpMask[j];
                helpDataF=helpData[j];
                helpCastDataF=helpCastData[j];
                for (k=0;k<m;k++)
                  {
                    if (helpMaskF[k]>0)
                      {
                         res=res+helpDataF[k]*helpCastDataF[k];
                      }
                  }
             }
          } //цикл по первому измерению
        } //если есть акт. индексы
        //Конец Выполнение вычисления скалярного произведения
        return res;
   }
//*************************************
 double TRnRelease3dSpace::LocalReleaseInnerproduct(const int id,const TLinearOperator *lA1,const TRnData *s2, int type)
  {
      const TLinearOperator *A1 = lA1;

      if (A1 == NULL)
       {
          throw "TRnRelease3dSpace::LocalReleaseInnerproduct: First operator is NULL.";
       }

       //Попытка приведения параметра к TRnRelease3dLinearOperator
       try {dynamic_cast<const TRnRelease3dLinearOperator *>(A1);}
       catch (...) {throw "Error in TRnRelease3dSpace::LocalReleaseInnerproduct: Could not cast parametr pased to TRnRelease3dLinearOperator * type.";}
       //Конец Попытка приведения параметра к TRnRelease3dLinearOperator


        if  (s2!=NULL)
        {
         //Попытка приведения параметра к TRnRelease3dData
          try {dynamic_cast<const TRnRelease3dData *>(s2);}
          catch (...) {throw "TRnRelease3dSpace::LocalReleaseInnerproduct: Could not cast parametr pased to TRnRelease3dData * type.";}
         //Конец Попытка приведения параметра к TRnRelease3dData
        }

       double res=0;

       const T3dNumberMask & sMask = fCastModel->GetMask(id);
       const int * const * const * mask = sMask.Mask;
       int n = sMask.N;
       int l = sMask.L;
       int m = sMask.M;
       int i; int j; int k;


       double ***fCastData=dynamic_cast<const TRnRelease3dData &>(*fMultiData[id]).fData;
       

       
       double **helpfCastData;
       double *helpfCastDataF;


      if (s2==NULL)
        {
           //Выполнение быстрой линейной комбинации
             const TRnRelease3dLinearOperator &A1Cast=dynamic_cast<const TRnRelease3dLinearOperator &>(*A1);
             for (i=0;i<n;i++)
              {
                 const int * const * helpMask = mask[i];
                 helpfCastData=fCastData[i];

                 for (j=0;j<l;j++)
                   {
                     const int * helpMaskF = helpMask[j];
                     helpfCastDataF=helpfCastData[j];

                     for (k=0;k<m;k++)
                       {
                          if (helpMaskF[k]>0)
                            {
                               res=res+helpfCastDataF[k]*A1Cast(id,i,j,k);
                            }
                       }
                   }

              }
           //Конец Выполнение быстрой линейной комбинации

        } //Конец нет второго операнда

      else
        {
           throw "TRnRelease3dSpace::LocalReleaseInnerproduct: The function is not implemented for non-NULL second operand.";
        } //Конец есть второй операнд

      return res;
  }
//*************************************
   double &TRnRelease3dSpace::GetByRnIndex(int i) const
    {
       const T3dCorrespondence &c = dynamic_cast<const T3dCorrespondence &>(fCastModel->GetCorrespondence(i));
       double ***fCastData=dynamic_cast<const TRnRelease3dData &>(*fMultiData[c.Ind]).fData;
       return fCastData[c.I][c.J][c.K];
    }
//*************************************
  double TRnRelease3dSpace::LocalGetSpecificElement(const int id,int type,int &elementIndex) const
    {
       double res=0;

        double ***fCastData=dynamic_cast<const TRnRelease3dData &>(*fMultiData[id]).fData;
        //Выполнение вычисления
        const T3dNumberMask & sMask = fCastModel->GetMask(id);
        const int * const * const * mask = sMask.Mask;
        int n = sMask.N;
        int l = sMask.L;
        int m = sMask.M;
        int i; int j; int k;

        double localRes;
        int localIndex = 1;
        bool assignLocalIndex;
        const T3dInverseCorrespondence &invCor = dynamic_cast<const T3dInverseCorrespondence &>(fCastModel->GetInverseCorrespondence());
        for (i=0;i<n;i++)
          {
            for (j=0;j<l;j++)
             {
               for (k=0;k<m;k++)
                {

                  if (mask[i][j][k]>0)
                    {
                      assignLocalIndex = false;
                      localRes=fabs(fCastData[i][j][k]);
                      if (type == 1)
                       {
                         if (localRes>res)
                          {
                            res = localRes;
                            assignLocalIndex = true;
                          }
                       }
                      else if (type == 2)
                       {
                         if (localRes<res)
                          {
                            res = localRes;
                            assignLocalIndex = true;
                          }
                       }
                      else
                       {
                         throw "Error in TRnRelease3dSpace::GetSpecificElement: The fuction is not implemented for this type.";
                       }
                      if (assignLocalIndex==true)
                       {
                         localIndex = invCor[id][i][j][k];
                       }
                    }
                }
             }
          }
        //Конец Выполнение вычисления
        elementIndex = localIndex;
        return res;
    }
//*************************************

//#################################################################



//---------------------Rn3 (End)--------------------------------
