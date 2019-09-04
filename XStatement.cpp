#pragma hdrstop

#include "XStatement.h"
#include "XProblem.h"
#include "XOperator.h"
#include <csignal>
#include "math.h"

//#define M_PI 3.14159265358979323846
//#define M_PI_4 M_PI/4.

//#################################################################


//*************************************
void TMaskGenerator::InitMemory()
{
	// Узлы сетки и шаги
	glAllocVector(X, Nx);
	glAllocVector(Y, Ny);
	glAllocVector(Z, Nz);
	glAllocVector(Hx, Nx);
	glAllocVector(Hy, Ny);
	glAllocVector(Hz, Nz);

	// Маска
	glAllocMatrix(Mask, Nx, Ny, Nz);

	// смещенные сетки
	glAllocVector(XU, NxU);
	glAllocVector(YU, NyU);
	glAllocVector(ZU, NzU);
	glAllocVector(HxU, NxU);
	glAllocVector(HyU, NyU);
	glAllocVector(HzU, NzU);

	glAllocVector(XV, NxV);
	glAllocVector(YV, NyV);
	glAllocVector(ZV, NzV);
	glAllocVector(HxV, NxV);
	glAllocVector(HyV, NyV);
	glAllocVector(HzV, NzV);

	glAllocVector(XW, NxW);
	glAllocVector(YW, NyW);
	glAllocVector(ZW, NzW);
	glAllocVector(HxW, NxW);
	glAllocVector(HyW, NyW);
	glAllocVector(HzW, NzW);

	glAllocVector(XP, NxP);
	glAllocVector(YP, NyP);
	glAllocVector(ZP, NzP);
	glAllocVector(HxP, NxP);
	glAllocVector(HyP, NyP);
	glAllocVector(HzP, NzP);

	glAllocMatrix(MaskU, NxU, NyU, NzU);
	glAllocMatrix(MaskV, NxV, NyV, NzV);
	glAllocMatrix(MaskW, NxW, NyW, NzW);
	glAllocMatrix(MaskP, NxP, NyP, NzP);
}
//*************************************



//#################################################################



//*************************************
void TSimpleWithPressureMaskGenerator::CreateGeometryNodes()
{
	// Узлы, шаги сетки
	for(int i=0; i<Nx; ++i)
		X[i] = LengthX*i/(Nx-1);
	for(int i=1; i<Nx; ++i)
		Hx[i] = X[i]-X[i-1];
	Hx[0] = Hx[1];

	for(int j=0; j<Ny; ++j)
		Y[j] = LengthY*j/(Ny-1);
	for(int j=1; j<Ny; ++j)
		Hy[j] = Y[j]-Y[j-1];
	Hy[0] = Hy[1];

	for(int k=0; k<Nz; ++k)
		Z[k] = LengthZ*k/(Nz-1);
	for(int k=1; k<Nz; ++k)
		Hz[k] = Z[k]-Z[k-1];
	Hy[0] = Hy[1];
}
//*************************************
void TSimpleWithPressureMaskGenerator::CreateGeometryGrid()
{
	// Сначала все точки - расчетные
	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz; ++k)
			{
				Mask[i][j][k].mask = TActualPoint;
				Mask[i][j][k].normal = 0;
			}

	// Левая и правая границы
	for(int j=1; j<Ny-1; ++j)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[Nx-1][j][k].mask = Mask[0][j][k].mask = TBorderPoint;
			Mask[0][j][k].normal = new glTVector(-1,0,0);
			Mask[Nx-1][j][k].normal = new glTVector(1,0,0);
		}

	// Ближняя и дальняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[i][Ny-1][k].mask = Mask[i][0][k].mask = TBorderPoint;
			Mask[i][0][k].normal = new glTVector(0,-1,0);
			Mask[i][Ny-1][k].normal = new glTVector(0,1,0);
		}

	// Нижняя и верхняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int j=1; j<Ny-1; ++j)
		{
			Mask[i][j][Nz-1].mask = Mask[i][j][0].mask = TBorderPoint;
			Mask[i][j][0].normal = new glTVector(0,0,-1);
			Mask[i][j][Nz-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<Nx; ++i)
	{
		Mask[i][0][0].mask = TFictivePoint;
		Mask[i][0][Nz-1].mask = TFictivePoint;
		Mask[i][Ny-1][0].mask = TFictivePoint;
		Mask[i][Ny-1][Nz-1].mask = TFictivePoint;
	}

	for(int j=0; j<Ny; ++j)
	{
		Mask[0][j][0].mask = TFictivePoint;
		Mask[0][j][Nz-1].mask = TFictivePoint;
		Mask[Nx-1][j][0].mask = TFictivePoint;
		Mask[Nx-1][j][Nz-1].mask = TFictivePoint;
	}

	for(int k=0; k<Nz; ++k)
	{
		Mask[0][0][k].mask = TFictivePoint;
		Mask[0][Ny-1][k].mask = TFictivePoint;
		Mask[Nx-1][0][k].mask = TFictivePoint;
		Mask[Nx-1][Ny-1][k].mask = TFictivePoint;
	}
}
//*************************************
void TSimpleWithPressureMaskGenerator::CreateMultyNodes()
{
	// Узлы для компоненты U скорости
	for(int i=0; i<NxU; ++i)
		XU[i] = X[i];

	for(int j=1; j<NyU-1; ++j)
		YU[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YU[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YU[NyU-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzU-1; ++k)
		ZU[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZU[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZU[NzU-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxU; ++i)
		HxU[i] = XU[i] - XU[i-1];
	HxU[0] = HxU[1];

	for(int j=1; j<NyU; ++j)
		HyU[j] = YU[j] - YU[j-1];
	HyU[0] = HyU[1];

	for(int k=1; k<NzU; ++k)
		HzU[k] = ZU[k] - ZU[k-1];
	HzU[0] = HzU[1];

	// Узлы для компоненты V скорости
	for(int i=1; i<NxV-1; ++i)
		XV[i] = ( X[i-1] + X[i] ) / 2.0;
	XV[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XV[NxV-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=0; j<NyV; ++j)
		YV[j] = Y[j];

	for(int k=1; k<NzV-1; ++k)
		ZV[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZV[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZV[NzV-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxV; ++i)
		HxV[i] = XV[i] - XV[i-1];
	HxV[0] = HxV[1];

	for(int j=1; j<NyV; ++j)
		HyV[j] = YV[j] - YV[j-1];
	HyV[0] = HyV[1];

	for(int k=1; k<NzV; ++k)
		HzV[k] = ZV[k] - ZV[k-1];
	HzV[0] = HzV[1];

	// Узлы для компоненты W скорости
	for(int i=1; i<NxW-1; ++i)
		XW[i] = ( X[i-1] + X[i] ) / 2.0;
	XW[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XW[NxW-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyW-1; ++j)
		YW[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YW[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YW[NyW-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=0; k<NzW; ++k)
		ZW[k] = Z[k];

	for(int i=1; i<NxW; ++i)
		HxW[i] = XW[i] - XW[i-1];
	HxW[0] = HxW[1];

	for(int j=1; j<NyW; ++j)
		HyW[j] = YW[j] - YW[j-1];
	HyW[0] = HyW[1];

	for(int k=1; k<NzW; ++k)
		HzW[k] = ZW[k] - ZW[k-1];
	HzW[0] = HzW[1];

	// Узлы для давления
	for(int i=1; i<NxP-1; ++i)
		XP[i] = ( X[i-1] + X[i] ) / 2.0;
	XP[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XP[NxP-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyP-1; ++j)
		YP[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YP[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YP[NyP-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzP-1; ++k)
		ZP[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZP[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZP[NzP-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxP; ++i)
		HxP[i] = XP[i] - XP[i-1];
	HxP[0] = HxP[1];

	for(int j=1; j<NyP; ++j)
		HyP[j] = YP[j] - YP[j-1];
	HyP[0] = HyP[1];

	for(int k=1; k<NzP; ++k)
		HzP[k] = ZP[k] - ZP[k-1];
	HzP[0] = HzP[1];
}
//*************************************
void TSimpleWithPressureMaskGenerator::CreateMultyGrids()
{
		// Заполнение маски для U //

	for(int i=0; i<NxU; ++i)
		for(int j=0; j<NyU; ++j)
			for(int k=0; k<NzU; ++k)
			{
				MaskU[i][j][k].mask = TActualPoint;
				MaskU[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyU-1; ++j)
		for(int k=1; k<NzU-1; ++k)
		{			
			MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TNormalBorderPoint;
			//MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TEquationBorderPoint;
			MaskU[0][j][k].normal = new glTVector(-1,0,0);
			MaskU[NxU-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[i][0][k].mask = MaskU[i][NyU-1][k].mask = TPreDefinedBorderPoint;
			MaskU[i][0][k].normal =  new glTVector(0,-1,0);
			MaskU[i][NyU-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int j=1; j<NyU-1; ++j)
		{
			MaskU[i][j][0].mask = MaskU[i][j][NzU-1].mask = TPreDefinedBorderPoint;
			MaskU[i][j][0].normal = new glTVector(0,0,-1);
			MaskU[i][j][NzU-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxU; ++i)
	{
		MaskU[i][0][0].mask = TFictivePoint;
		MaskU[i][0][NzU-1].mask = TFictivePoint;
		MaskU[i][NyU-1][0].mask = TFictivePoint;
		MaskU[i][NyU-1][NzU-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyU; ++j)
	{
		MaskU[0][j][0].mask = TFictivePoint;
		MaskU[0][j][NzU-1].mask = TFictivePoint;
		MaskU[NxU-1][j][0].mask = TFictivePoint;
		MaskU[NxU-1][j][NzU-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzU; ++k)
	{
		MaskU[0][0][k].mask = TFictivePoint;
		MaskU[0][NyU-1][k].mask = TFictivePoint;
		MaskU[NxU-1][0][k].mask = TFictivePoint;
		MaskU[NxU-1][NyU-1][k].mask = TFictivePoint;
	}


	// Заполнение маски для V //

	for(int i=0; i<NxV; ++i)
		for(int j=0; j<NyV; ++j)
			for(int k=0; k<NzV; ++k)
			{
				MaskV[i][j][k].mask = TActualPoint;
				MaskV[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyV-1; ++j)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[0][j][k].mask = MaskV[NxV-1][j][k].mask = TPreDefinedBorderPoint;
			MaskV[0][j][k].normal = new glTVector(-1,0,0);
			MaskV[NxV-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[i][0][k].mask = MaskV[i][NyV-1][k].mask = TDefinedBorderPoint;
			MaskV[i][0][k].normal =  new glTVector(0,-1,0);
			MaskV[i][NyV-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int j=1; j<NyV-1; ++j)
		{
			MaskV[i][j][0].mask = MaskV[i][j][NzV-1].mask = TPreDefinedBorderPoint;
			MaskV[i][j][0].normal = new glTVector(0,0,-1);
			MaskV[i][j][NzV-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxV; ++i)
	{
		MaskV[i][0][0].mask = TFictivePoint;
		MaskV[i][0][NzV-1].mask = TFictivePoint;
		MaskV[i][NyV-1][0].mask = TFictivePoint;
		MaskV[i][NyV-1][NzV-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyV; ++j)
	{
		MaskV[0][j][0].mask = TFictivePoint;
		MaskV[0][j][NzV-1].mask = TFictivePoint;
		MaskV[NxV-1][j][0].mask = TFictivePoint;
		MaskV[NxV-1][j][NzV-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzV; ++k)
	{
		MaskV[0][0][k].mask = TFictivePoint;
		MaskV[0][NyV-1][k].mask = TFictivePoint;
		MaskV[NxV-1][0][k].mask = TFictivePoint;
		MaskV[NxV-1][NyV-1][k].mask = TFictivePoint;
	}


	// Заполнение маски для W //

	for(int i=0; i<NxW; ++i)
		for(int j=0; j<NyW; ++j)
			for(int k=0; k<NzW; ++k)
			{
				MaskW[i][j][k].mask = TActualPoint;
				MaskW[i][j][k].normal = 0;
			}
	// Границы области
	for(int j=1; j<NyW-1; ++j)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[0][j][k].mask = MaskW[NxW-1][j][k].mask = TPreDefinedBorderPoint; 
			MaskW[0][j][k].normal = new glTVector(-1,0,0);
			MaskW[NxW-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[i][0][k].mask = MaskW[i][NyW-1][k].mask = TPreDefinedBorderPoint;
			MaskW[i][0][k].normal =  new glTVector(0,-1,0);
			MaskW[i][NyW-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int j=1; j<NyW-1; ++j)
		{
			MaskW[i][j][0].mask = MaskW[i][j][NzW-1].mask = TDefinedBorderPoint;
			MaskW[i][j][0].normal = new glTVector(0,0,-1);
			MaskW[i][j][NzW-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxW; ++i)
	{
		MaskW[i][0][0].mask = TFictivePoint;
		MaskW[i][0][NzW-1].mask = TFictivePoint;
		MaskW[i][NyW-1][0].mask = TFictivePoint;
		MaskW[i][NyW-1][NzW-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyW; ++j)
	{
		MaskW[0][j][0].mask = TFictivePoint;
		MaskW[0][j][NzW-1].mask = TFictivePoint;
		MaskW[NxW-1][j][0].mask = TFictivePoint;
		MaskW[NxW-1][j][NzW-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzW; ++k)
	{
		MaskW[0][0][k].mask = TFictivePoint;
		MaskW[0][NyW-1][k].mask = TFictivePoint;
		MaskW[NxW-1][0][k].mask = TFictivePoint;
		MaskW[NxW-1][NyW-1][k].mask = TFictivePoint;
	}


	// Заполнение маски для P
	for(int i=0; i<NxP; ++i)
		for(int j=0; j<NyP; ++j)
			for(int k=0; k<NzP; ++k)
			{
				MaskP[i][j][k].mask = TActualPoint;
				MaskP[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyP-1; ++j)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreDefinedBorderPoint;
			MaskP[0][j][k].normal = new glTVector(-1,0,0);
			MaskP[NxP-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[i][0][k].mask = MaskP[i][NyP-1][k].mask = TFictivePoint;
			MaskP[i][0][k].normal =  new glTVector(0,-1,0);
			MaskP[i][NyP-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int j=1; j<NyP-1; ++j)
		{
			MaskP[i][j][0].mask = MaskP[i][j][NzP-1].mask = TFictivePoint;
			MaskP[i][j][0].normal = new glTVector(0,0,-1);
			MaskP[i][j][NzP-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxP; ++i)
	{
		MaskP[i][0][0].mask = TFictivePoint;
		MaskP[i][0][NzP-1].mask = TFictivePoint;
		MaskP[i][NyP-1][0].mask = TFictivePoint;
		MaskP[i][NyP-1][NzP-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyP; ++j)
	{
		MaskP[0][j][0].mask = TFictivePoint;
		MaskP[0][j][NzP-1].mask = TFictivePoint;
		MaskP[NxP-1][j][0].mask = TFictivePoint;
		MaskP[NxP-1][j][NzP-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzP; ++k)
	{
		MaskP[0][0][k].mask = TFictivePoint;
		MaskP[0][NyP-1][k].mask = TFictivePoint;
		MaskP[NxP-1][0][k].mask = TFictivePoint;
		MaskP[NxP-1][NyP-1][k].mask = TFictivePoint;
	}
}
//*************************************
double TSimpleWithPressureMaskGenerator::getUCondition(int i, int j, int k, int n)
{	
	return 0;
}
//*************************************
double TSimpleWithPressureMaskGenerator::getVCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TSimpleWithPressureMaskGenerator::getWCondition(int i, int j, int k, int n)
{	
	return 0;
}
//*************************************
double TSimpleWithPressureMaskGenerator::getPCondition(int i, int j, int k, int n)
{
	if(i==0)
		return PLeft;
	else if(i==NxP-1)
		return PRight;
	else
		return 0;
}
//*************************************



//#################################################################



//*************************************
void TCubeInsideWithPressureMaskGenerator::CreateGeometryNodes()
{
	for(int i=0; i<Nx; ++i)
		X[i] = LengthX*i/(Nx-1);
	for(int i=1; i<Nx; ++i)
		Hx[i] = X[i]-X[i-1];
	Hx[0] = Hx[1];

	for(int j=0; j<Ny; ++j)
		Y[j] = LengthY*j/(Ny-1);
	for(int j=1; j<Ny; ++j)
		Hy[j] = Y[j]-Y[j-1];
	Hy[0] = Hy[1];

	for(int k=0; k<Nz; ++k)
		Z[k] = LengthZ*k/(Nz-1);
	for(int k=1; k<Nz; ++k)
		Hz[k] = Z[k]-Z[k-1];
	Hy[0] = Hy[1];
}
//*************************************
void TCubeInsideWithPressureMaskGenerator::CreateGeometryGrid()
{
	// Сначала все точки - расчетные
	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz; ++k)
			{
				Mask[i][j][k].mask = TActualPoint;
				Mask[i][j][k].normal = 0;
			}

	// Левая и правая границы
	for(int j=1; j<Ny-1; ++j)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[Nx-1][j][k].mask = Mask[0][j][k].mask = TBorderPoint;
			Mask[0][j][k].normal = new glTVector(-1,0,0);
			Mask[Nx-1][j][k].normal = new glTVector(1,0,0);
		}

	// Ближняя и дальняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[i][Ny-1][k].mask = Mask[i][0][k].mask = TBorderPoint;
			Mask[i][0][k].normal = new glTVector(0,-1,0);
			Mask[i][Ny-1][k].normal = new glTVector(0,1,0);
		}

	// Нижняя и верхняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int j=1; j<Ny-1; ++j)
		{
			Mask[i][j][Nz-1].mask = Mask[i][j][0].mask = TBorderPoint;
			Mask[i][j][0].normal = new glTVector(0,0,-1);
			Mask[i][j][Nz-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<Nx; ++i)
	{
		Mask[i][0][0].mask = TFictivePoint;
		Mask[i][0][Nz-1].mask = TFictivePoint;
		Mask[i][Ny-1][0].mask = TFictivePoint;
		Mask[i][Ny-1][Nz-1].mask = TFictivePoint;
	}

	for(int j=0; j<Ny; ++j)
	{
		Mask[0][j][0].mask = TFictivePoint;
		Mask[0][j][Nz-1].mask = TFictivePoint;
		Mask[Nx-1][j][0].mask = TFictivePoint;
		Mask[Nx-1][j][Nz-1].mask = TFictivePoint;
	}

	for(int k=0; k<Nz; ++k)
	{
		Mask[0][0][k].mask = TFictivePoint;
		Mask[0][Ny-1][k].mask = TFictivePoint;
		Mask[Nx-1][0][k].mask = TFictivePoint;
		Mask[Nx-1][Ny-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1; i<Nx1+Nx2-2; ++i)
		for(int j=Ny1; j<Ny1+Ny2-2; ++j)
			for(int k=Nz1; k<Nz1+Nz2-2; ++k)
				Mask[i][j][k].mask = TFictivePoint;
	
	// Грани препятствия
	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		for(int k=Nz1-1; k<Nz1+Nz2-1; ++k)
		{
			Mask[Nx1-1][j][k].mask = Mask[Nx1+Nx2-2][j][k].mask = TBorderPoint;
			Mask[Nx1-1][j][k].normal = new glTVector(1,0,0);
			Mask[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int k=Nz1-1; k<Nz1+Nz2-1; ++k)
		{
			Mask[i][Ny1-1][k].mask = Mask[i][Ny1+Ny2-2][k].mask = TBorderPoint;
			Mask[i][Ny1-1][k].normal = new glTVector(0,1,0);
			Mask[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		{
			Mask[i][j][Nz1-1].mask = Mask[i][j][Nz1+Nz2-2].mask = TBorderPoint;
			Mask[i][j][Nz1-1].normal = new glTVector(0,0,1);
			Mask[i][j][Nz1+Nz2-2].normal = new glTVector(0,0,-1);
		}

	// Ребра препятствия
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
	{
		Mask[i][Ny1-1][Nz1-1].normal->set(1,1,0);
		Mask[i][Ny1+Ny2-2][Nz1-1].normal->set(-1,1,0);
		Mask[i][Ny1-1][Nz1+Nz2-2].normal->set(1,-1,0);
		Mask[i][Ny1+Ny2-2][Nz1+Nz2-2].normal->set(-1,-1,0);
	}

	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
	{
		Mask[Nx1-1][j][Nz1-1].normal->set(1,0,1);
		Mask[Nx1+Nx2-2][j][Nz1-1].normal->set(-1,0,1);
		Mask[Nx1-1][j][Nz1+Nz2-2].normal->set(1,0,-1);
		Mask[Nx1+Nx2-2][j][Nz1+Nz2-2].normal->set(-1,0,-1);
	}

	for(int k=Nz1-1; k<Nz1+Nz2-1; ++k)
	{
		Mask[Nx1-1][Ny1-1][k].normal->set(0,1,1);
		Mask[Nx1+Nx2-2][Ny1-1][k].normal->set(0,-1,1);
		Mask[Nx1-1][Ny1+Ny2-2][k].normal->set(0,1,-1);
		Mask[Nx1+Nx2-2][Ny1+Ny2-2][k].normal->set(0,-1,-1);
	}

	// Угловые точки препятствия
	Mask[Nx1-1][Ny1-1][Nz1-1].normal->set(1,1,1);
	Mask[Nx1+Nx2-2][Ny1-1][Nz1-1].normal->set(-1,1,1);
	Mask[Nx1-1][Ny1-1][Nz1+Nz2-2].normal->set(1,1,-1);
	Mask[Nx1+Nx2-2][Ny1-1][Nz1+Nz2-2].normal->set(-1,1,-1);
	Mask[Nx1-1][Ny1+Ny2-2][Nz1-1].normal->set(1,-1,1);
	Mask[Nx1+Nx2-2][Ny1+Ny2-2][Nz1-1].normal->set(-1,-1,1);
	Mask[Nx1-1][Ny1+Ny2-2][Nz1+Nz2-2].normal->set(1,-1,-1);
	Mask[Nx1+Nx2-2][Ny1+Ny2-2][Nz1+Nz2-2].normal->set(-1,-1,-1);
}
//*************************************
void TCubeInsideWithPressureMaskGenerator::CreateMultyNodes()
{
	// Узлы для компоненты U скорости
	for(int i=0; i<NxU; ++i)
		XU[i] = X[i];

	for(int j=1; j<NyU-1; ++j)
		YU[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YU[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YU[NyU-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzU-1; ++k)
		ZU[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZU[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZU[NzU-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxU; ++i)
		HxU[i] = XU[i] - XU[i-1];
	HxU[0] = HxU[1];

	for(int j=1; j<NyU; ++j)
		HyU[j] = YU[j] - YU[j-1];
	HyU[0] = HyU[1];

	for(int k=1; k<NzU; ++k)
		HzU[k] = ZU[k] - ZU[k-1];
	HzU[0] = HzU[1];

	// Узлы для компоненты V скорости
	for(int i=1; i<NxV-1; ++i)
		XV[i] = ( X[i-1] + X[i] ) / 2.0;
	XV[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XV[NxV-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=0; j<NyV; ++j)
		YV[j] = Y[j];

	for(int k=1; k<NzV-1; ++k)
		ZV[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZV[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZV[NzV-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxV; ++i)
		HxV[i] = XV[i] - XV[i-1];
	HxV[0] = HxV[1];

	for(int j=1; j<NyV; ++j)
		HyV[j] = YV[j] - YV[j-1];
	HyV[0] = HyV[1];

	for(int k=1; k<NzV; ++k)
		HzV[k] = ZV[k] - ZV[k-1];
	HzV[0] = HzV[1];

	// Узлы для компоненты W скорости
	for(int i=1; i<NxW-1; ++i)
		XW[i] = ( X[i-1] + X[i] ) / 2.0;
	XW[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XW[NxW-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyW-1; ++j)
		YW[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YW[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YW[NyW-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=0; k<NzW; ++k)
		ZW[k] = Z[k];

	for(int i=1; i<NxW; ++i)
		HxW[i] = XW[i] - XW[i-1];
	HxW[0] = HxW[1];

	for(int j=1; j<NyW; ++j)
		HyW[j] = YW[j] - YW[j-1];
	HyW[0] = HyW[1];

	for(int k=1; k<NzW; ++k)
		HzW[k] = ZW[k] - ZW[k-1];
	HzW[0] = HzW[1];

	// Узлы для давления
	for(int i=1; i<NxP-1; ++i)
		XP[i] = ( X[i-1] + X[i] ) / 2.0;
	XP[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XP[NxP-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyP-1; ++j)
		YP[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YP[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YP[NyP-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzP-1; ++k)
		ZP[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZP[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZP[NzP-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxP; ++i)
		HxP[i] = XP[i] - XP[i-1];
	HxP[0] = HxP[1];

	for(int j=1; j<NyP; ++j)
		HyP[j] = YP[j] - YP[j-1];
	HyP[0] = HyP[1];

	for(int k=1; k<NzP; ++k)
		HzP[k] = ZP[k] - ZP[k-1];
	HzP[0] = HzP[1];
}
//*************************************
void TCubeInsideWithPressureMaskGenerator::CreateMultyGrids()
{
	// Заполнение маски для U //

	for(int i=0; i<NxU; ++i)
		for(int j=0; j<NyU; ++j)
			for(int k=0; k<NzU; ++k)
			{
				MaskU[i][j][k].mask = TActualPoint;
				MaskU[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyU-1; ++j)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TDefinedBorderPoint;
			MaskU[0][j][k].normal = new glTVector(-1,0,0);
			MaskU[NxU-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int j=2; j<NyU-2; ++j)
		for(int k=2; k<NzU-2; ++k)
		{
			MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TEquationBorderPoint;
			MaskU[0][j][k].normal = new glTVector(-1,0,0);
			MaskU[NxU-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[i][0][k].mask = MaskU[i][NyU-1][k].mask = TPreDefinedBorderPoint;
			MaskU[i][0][k].normal =  new glTVector(0,-1,0);
			MaskU[i][NyU-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int j=1; j<NyU-1; ++j)
		{
			MaskU[i][j][0].mask = MaskU[i][j][NzU-1].mask = TPreDefinedBorderPoint;
			MaskU[i][j][0].normal = new glTVector(0,0,-1);
			MaskU[i][j][NzU-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxU; ++i)
	{
		MaskU[i][0][0].mask = TFictivePoint;
		MaskU[i][0][NzU-1].mask = TFictivePoint;
		MaskU[i][NyU-1][0].mask = TFictivePoint;
		MaskU[i][NyU-1][NzU-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyU; ++j)
	{
		MaskU[0][j][0].mask = TFictivePoint;
		MaskU[0][j][NzU-1].mask = TFictivePoint;
		MaskU[NxU-1][j][0].mask = TFictivePoint;
		MaskU[NxU-1][j][NzU-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzU; ++k)
	{
		MaskU[0][0][k].mask = TFictivePoint;
		MaskU[0][NyU-1][k].mask = TFictivePoint;
		MaskU[NxU-1][0][k].mask = TFictivePoint;
		MaskU[NxU-1][NyU-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1; i<Nx1+Nx2-2; ++i)
		for(int j=Ny1+1; j<Ny1+Ny2-2; ++j)
			for(int k=Nz1+1; k<Nz1+Nz2-2; ++k)
				MaskU[i][j][k].mask = TFictivePoint;

	// Грани препятствия
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskU[i][Ny1][k].mask = MaskU[i][Ny1+Ny2-2][k].mask = TPreDefinedBorderPoint;
			MaskU[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskU[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			MaskU[i][j][Nz1].mask = MaskU[i][j][Nz1+Nz2-2].mask = TPreDefinedBorderPoint;
			MaskU[i][j][Nz1].normal = new glTVector(0,0,1);
			MaskU[i][j][Nz1+Nz2-2].normal = new glTVector(0,0,-1);
		}

	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskU[Nx1-1][j][k].mask = MaskU[Nx1+Nx2-2][j][k].mask = TDefinedBorderPoint;
			MaskU[Nx1-1][j][k].normal = new glTVector(1,0,0);
			MaskU[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}

	// Заполнение маски для V //

	for(int i=0; i<NxV; ++i)
		for(int j=0; j<NyV; ++j)
			for(int k=0; k<NzV; ++k)
			{
				MaskV[i][j][k].mask = TActualPoint;
				MaskV[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyV-1; ++j)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[0][j][k].mask = MaskV[NxV-1][j][k].mask = TPreDefinedBorderPoint;
			MaskV[0][j][k].normal = new glTVector(-1,0,0);
			MaskV[NxV-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[i][0][k].mask = MaskV[i][NyV-1][k].mask = TDefinedBorderPoint;
			MaskV[i][0][k].normal =  new glTVector(0,-1,0);
			MaskV[i][NyV-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int j=1; j<NyV-1; ++j)
		{
			MaskV[i][j][0].mask = MaskV[i][j][NzV-1].mask = TPreDefinedBorderPoint;
			MaskV[i][j][0].normal = new glTVector(0,0,-1);
			MaskV[i][j][NzV-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxV; ++i)
	{
		MaskV[i][0][0].mask = TFictivePoint;
		MaskV[i][0][NzV-1].mask = TFictivePoint;
		MaskV[i][NyV-1][0].mask = TFictivePoint;
		MaskV[i][NyV-1][NzV-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyV; ++j)
	{
		MaskV[0][j][0].mask = TFictivePoint;
		MaskV[0][j][NzV-1].mask = TFictivePoint;
		MaskV[NxV-1][j][0].mask = TFictivePoint;
		MaskV[NxV-1][j][NzV-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzV; ++k)
	{
		MaskV[0][0][k].mask = TFictivePoint;
		MaskV[0][NyV-1][k].mask = TFictivePoint;
		MaskV[NxV-1][0][k].mask = TFictivePoint;
		MaskV[NxV-1][NyV-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1+1; i<Nx1+Nx2-2; ++i)
		for(int j=Ny1; j<Ny1+Ny2-2; ++j)
			for(int k=Nz1+1; k<Nz1+Nz2-2; ++k)
				MaskV[i][j][k].mask = TFictivePoint;
	
	// Грани препятствия
	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskV[Nx1][j][k].mask = MaskV[Nx1+Nx2-2][j][k].mask = TPreDefinedBorderPoint;
			MaskV[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskV[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		{
			MaskV[i][j][Nz1].mask = MaskV[i][j][Nz1+Nz2-2].mask = TPreDefinedBorderPoint;
			MaskV[i][j][Nz1].normal = new glTVector(0,0,1);
			MaskV[i][j][Nz1+Nz2-2].normal = new glTVector(0,0,-1);
		}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskV[i][Ny1-1][k].mask = MaskV[i][Ny1+Ny2-2][k].mask = TDefinedBorderPoint;
			MaskV[i][Ny1-1][k].normal = new glTVector(0,1,0);
			MaskV[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	// Заполнение маски для W //

	for(int i=0; i<NxW; ++i)
		for(int j=0; j<NyW; ++j)
			for(int k=0; k<NzW; ++k)
			{
				MaskW[i][j][k].mask = TActualPoint;
				MaskW[i][j][k].normal = 0;
			}
	// Границы области
	for(int j=1; j<NyW-1; ++j)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[0][j][k].mask = MaskW[NxW-1][j][k].mask = TPreDefinedBorderPoint; 
			MaskW[0][j][k].normal = new glTVector(-1,0,0);
			MaskW[NxW-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[i][0][k].mask = MaskW[i][NyW-1][k].mask = TPreDefinedBorderPoint;
			MaskW[i][0][k].normal =  new glTVector(0,-1,0);
			MaskW[i][NyW-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int j=1; j<NyW-1; ++j)
		{
			MaskW[i][j][0].mask = MaskW[i][j][NzW-1].mask = TDefinedBorderPoint;
			MaskW[i][j][0].normal = new glTVector(0,0,-1);
			MaskW[i][j][NzW-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxW; ++i)
	{
		MaskW[i][0][0].mask = TFictivePoint;
		MaskW[i][0][NzW-1].mask = TFictivePoint;
		MaskW[i][NyW-1][0].mask = TFictivePoint;
		MaskW[i][NyW-1][NzW-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyW; ++j)
	{
		MaskW[0][j][0].mask = TFictivePoint;
		MaskW[0][j][NzW-1].mask = TFictivePoint;
		MaskW[NxW-1][j][0].mask = TFictivePoint;
		MaskW[NxW-1][j][NzW-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzW; ++k)
	{
		MaskW[0][0][k].mask = TFictivePoint;
		MaskW[0][NyW-1][k].mask = TFictivePoint;
		MaskW[NxW-1][0][k].mask = TFictivePoint;
		MaskW[NxW-1][NyW-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1+1; i<Nx1+Nx2-2; ++i)
		for(int j=Ny1+1; j<Ny1+Ny2-2; ++j)
			for(int k=Nz1; k<Nz1+Nz2-2; ++k)
				MaskW[i][j][k].mask = TFictivePoint;

	// Грани препятствия
	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		for(int k=Nz1-1; k<Nz1+Nz2-1; ++k)
		{
			MaskW[Nx1][j][k].mask = MaskW[Nx1+Nx2-2][j][k].mask = TPreDefinedBorderPoint;
			MaskW[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskW[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int k=Nz1-1; k<Nz1+Nz2-1; ++k)
		{
			MaskW[i][Ny1][k].mask = MaskW[i][Ny1+Ny2-2][k].mask = TPreDefinedBorderPoint;
			MaskW[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskW[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			MaskW[i][j][Nz1-1].mask = MaskW[i][j][Nz1+Nz2-2].mask = TDefinedBorderPoint;
			MaskW[i][j][Nz1-1].normal = new glTVector(0,0,1);
			MaskW[i][j][Nz1+Nz2-2].normal = new glTVector(0,0,-1);
		}

	// Заполнение маски для P
	for(int i=0; i<NxP; ++i)
		for(int j=0; j<NyP; ++j)
			for(int k=0; k<NzP; ++k)
			{
				MaskP[i][j][k].mask = TActualPoint;
				MaskP[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyP-1; ++j)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreNormalBorderPoint;
			MaskP[0][j][k].normal = new glTVector(-1,0,0);
			MaskP[NxP-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int j=2; j<NyP-2; ++j)
		for(int k=2; k<NzP-2; ++k)
		{
			MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreDefinedBorderPoint; 
			MaskP[0][j][k].normal = new glTVector(-1,0,0);
			MaskP[NxP-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[i][0][k].mask = MaskP[i][NyP-1][k].mask = TPreNormalBorderPoint;
			MaskP[i][0][k].normal =  new glTVector(0,-1,0);
			MaskP[i][NyP-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int j=1; j<NyP-1; ++j)
		{
			MaskP[i][j][0].mask = MaskP[i][j][NzP-1].mask = TPreNormalBorderPoint;
			MaskP[i][j][0].normal = new glTVector(0,0,-1);
			MaskP[i][j][NzP-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxP; ++i)
	{
		MaskP[i][0][0].mask = TFictivePoint;
		MaskP[i][0][NzP-1].mask = TFictivePoint;
		MaskP[i][NyP-1][0].mask = TFictivePoint;
		MaskP[i][NyP-1][NzP-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyP; ++j)
	{
		MaskP[0][j][0].mask = TFictivePoint;
		MaskP[0][j][NzP-1].mask = TFictivePoint;
		MaskP[NxP-1][j][0].mask = TFictivePoint;
		MaskP[NxP-1][j][NzP-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzP; ++k)
	{
		MaskP[0][0][k].mask = TFictivePoint;
		MaskP[0][NyP-1][k].mask = TFictivePoint;
		MaskP[NxP-1][0][k].mask = TFictivePoint;
		MaskP[NxP-1][NyP-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1+1; i<Nx1+Nx2-2; ++i)
		for(int j=Ny1+1; j<Ny1+Ny2-2; ++j)
			for(int k=Nz1+1; k<Nz1+Nz2-2; ++k)
				MaskP[i][j][k].mask = TFictivePoint;
	
	// Грани препятствия
	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskP[Nx1][j][k].mask = MaskP[Nx1+Nx2-2][j][k].mask = TPreNormalBorderPoint;
			MaskP[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskP[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskP[i][Ny1][k].mask = MaskP[i][Ny1+Ny2-2][k].mask = TPreNormalBorderPoint;
			MaskP[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskP[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			MaskP[i][j][Nz1].mask = MaskP[i][j][Nz1+Nz2-2].mask = TPreNormalBorderPoint;
			MaskP[i][j][Nz1].normal = new glTVector(0,0,1);
			MaskP[i][j][Nz1+Nz2-2].normal = new glTVector(0,0,-1);
		}

	// Ребра препятствия
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	{
		MaskP[i][Ny1][Nz1].normal->set(1,1,0);
		MaskP[i][Ny1+Ny2-2][Nz1].normal->set(-1,1,0);
		MaskP[i][Ny1][Nz1+Nz2-2].normal->set(1,-1,0);
		MaskP[i][Ny1+Ny2-2][Nz1+Nz2-2].normal->set(-1,-1,0);
	}

	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	{
		MaskP[Nx1][j][Nz1].normal->set(1,0,1);
		MaskP[Nx1+Nx2-2][j][Nz1].normal->set(-1,0,1);
		MaskP[Nx1][j][Nz1+Nz2-2].normal->set(1,0,-1);
		MaskP[Nx1+Nx2-2][j][Nz1+Nz2-2].normal->set(-1,0,-1);
	}

	for(int k=Nz1; k<Nz1+Nz2-1; ++k)
	{
		MaskP[Nx1][Ny1][k].normal->set(0,1,1);
		MaskP[Nx1+Nx2-2][Ny1][k].normal->set(0,-1,1);
		MaskP[Nx1][Ny1+Ny2-2][k].normal->set(0,1,-1);
		MaskP[Nx1+Nx2-2][Ny1+Ny2-2][k].normal->set(0,-1,-1);
	}

	// Угловые точки препятствия
	MaskP[Nx1][Ny1][Nz1].normal->set(1,1,1);
	MaskP[Nx1+Nx2-2][Ny1][Nz1].normal->set(-1,1,1);
	MaskP[Nx1][Ny1][Nz1+Nz2-2].normal->set(1,1,-1);
	MaskP[Nx1+Nx2-2][Ny1][Nz1+Nz2-2].normal->set(-1,1,-1);
	MaskP[Nx1][Ny1+Ny2-2][Nz1].normal->set(1,-1,1);
	MaskP[Nx1+Nx2-2][Ny1+Ny2-2][Nz1].normal->set(-1,-1,1);
	MaskP[Nx1][Ny1+Ny2-2][Nz1+Nz2-2].normal->set(1,-1,-1);
	MaskP[Nx1+Nx2-2][Ny1+Ny2-2][Nz1+Nz2-2].normal->set(-1,-1,-1);
}
//*************************************
double TCubeInsideWithPressureMaskGenerator::getUCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TCubeInsideWithPressureMaskGenerator::getVCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TCubeInsideWithPressureMaskGenerator::getWCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TCubeInsideWithPressureMaskGenerator::getPCondition(int i, int j, int k, int n)
{
	if(i==0)
		return PInner;
	else if(i==NxP-1)
		return POuter;
	else
		return 0;
}
//*************************************




//#################################################################




//*************************************
void TTwoParallelepipedsWithPressureMaskGenerator::CreateGeometryNodes()
{
	// Узлы, шаги сетки
	for(int i=0; i<Nx; ++i)
		X[i] = LengthX*i/(Nx-1);
	for(int i=1; i<Nx; ++i)
		Hx[i] = X[i]-X[i-1];
	Hx[0] = Hx[1];

	for(int j=0; j<Ny; ++j)
		Y[j] = LengthY*j/(Ny-1);
	for(int j=1; j<Ny; ++j)
		Hy[j] = Y[j]-Y[j-1];
	Hy[0] = Hy[1];

	for(int k=0; k<Nz; ++k)
		Z[k] = LengthZ*k/(Nz-1);
	for(int k=1; k<Nz; ++k)
		Hz[k] = Z[k]-Z[k-1];
	Hy[0] = Hy[1];
}
//*************************************
void TTwoParallelepipedsWithPressureMaskGenerator::CreateGeometryGrid()
{
	// Сначала все точки - расчетные
	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz; ++k)
			{
				Mask[i][j][k].mask = TActualPoint;
				Mask[i][j][k].normal = 0;
			}

	// Левая и правая границы
	for(int j=1; j<Ny-1; ++j)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[Nx-1][j][k].mask = Mask[0][j][k].mask = TBorderPoint;
			Mask[0][j][k].normal = new glTVector(-1,0,0);
			Mask[Nx-1][j][k].normal = new glTVector(1,0,0);
		}

	// Ближняя и дальняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[i][Ny-1][k].mask = Mask[i][0][k].mask = TBorderPoint;
			Mask[i][0][k].normal = new glTVector(0,-1,0);
			Mask[i][Ny-1][k].normal = new glTVector(0,1,0);
		}

	// Нижняя и верхняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int j=1; j<Ny-1; ++j)
		{
			Mask[i][j][Nz-1].mask = Mask[i][j][0].mask = TBorderPoint;
			Mask[i][j][0].normal = new glTVector(0,0,-1);
			Mask[i][j][Nz-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<Nx; ++i)
	{
		Mask[i][0][0].mask = TFictivePoint;
		Mask[i][0][Nz-1].mask = TFictivePoint;
		Mask[i][Ny-1][0].mask = TFictivePoint;
		Mask[i][Ny-1][Nz-1].mask = TFictivePoint;
	}

	for(int j=0; j<Ny; ++j)
	{
		Mask[0][j][0].mask = TFictivePoint;
		Mask[0][j][Nz-1].mask = TFictivePoint;
		Mask[Nx-1][j][0].mask = TFictivePoint;
		Mask[Nx-1][j][Nz-1].mask = TFictivePoint;
	}

	for(int k=0; k<Nz; ++k)
	{
		Mask[0][0][k].mask = TFictivePoint;
		Mask[0][Ny-1][k].mask = TFictivePoint;
		Mask[Nx-1][0][k].mask = TFictivePoint;
		Mask[Nx-1][Ny-1][k].mask = TFictivePoint;
	}

	// Внутренний параллелепипед
	for(int i=Nx1+1; i<Nx1+Nx2-2; ++i)
		for(int j=Ny1+1; j<Ny1+Ny2-2; ++j)
			for(int k=0; k<NzP; ++k)
				Mask[i][j][k].mask = TFictivePoint; // внутренность

	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		for(int k=0; k<NzP; ++k)
		{
			Mask[Nx1][j][k].mask = Mask[Nx1+Nx2-2][j][k].mask = TBorderPoint;
			Mask[Nx1][j][k].normal = new glTVector(1,0,0);
			Mask[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int k=0; k<NzP; ++k)
		{
			Mask[i][Ny1][k].mask = Mask[i][Ny1+Ny2-2][k].mask = TBorderPoint;
			Mask[i][Ny1][k].normal = new glTVector(0,1,0);
			Mask[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	// Ребра препятствия
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	{
		Mask[i][Ny1][0].mask = TFictivePoint;
		Mask[i][Ny1+Ny2-2][0].mask = TFictivePoint;
		Mask[i][Ny1][Nz-1].mask = TFictivePoint;
		Mask[i][Ny1+Ny2-2][Nz-1].mask = TFictivePoint;
	}

	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	{
		Mask[Nx1][j][0].mask = TFictivePoint;
		Mask[Nx1+Nx2-2][j][0].mask = TFictivePoint;
		Mask[Nx1][j][Nz-1].mask = TFictivePoint;
		Mask[Nx1+Nx2-2][j][Nz-1].mask = TFictivePoint;
	}

	for(int k=1; k<Nz-1; ++k)
	{
		Mask[Nx1][Ny1][k].normal = new glTVector(0,1,1);
		Mask[Nx1+Nx2-2][Ny1][k].normal = new glTVector(0,-1,1);
		Mask[Nx1][Ny1+Ny2-2][k].normal = new glTVector(0,1,-1);
		Mask[Nx1+Nx2-2][Ny1+Ny2-2][k].normal = new glTVector(0,-1,-1);
	}

	// Отверстия на внутреннем цилиндре
	for(int j=Ny4; j<Ny4+Ny5-1; ++j)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			Mask[Nx1][j][k].mask = Mask[Nx1+Nx2-2][j][k].mask = TBorderPoint;
			Mask[Nx1][j][k].normal = new glTVector(1,0,0);
			Mask[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx4; i<Nx4+Nx5-1; ++i)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			Mask[i][Ny1][k].mask = Mask[i][Ny1+Ny2-2][k].mask = TBorderPoint;
			Mask[i][Ny1][k].normal = new glTVector(0,1,0);
			Mask[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	// Отверстия на внешнем цилиндре
	for(int j=Ny6; j<Ny6+Ny7-1; ++j)
		for(int k=Nz3; k<Nz3+Nz4-1; ++k)
		{
			Mask[0][j][k].mask = Mask[Nx-1][j][k].mask = TBorderPoint;
			Mask[0][j][k].normal = new glTVector(1,0,0);
			Mask[Nx-1][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx6; i<Nx6+Nx7-1; ++i)
		for(int k=Nz3; k<Nz3+Nz4-1; ++k)
		{
			Mask[i][0][k].mask = Mask[i][Ny-1][k].mask = TBorderPoint;
			Mask[i][0][k].normal = new glTVector(0,1,0);
			Mask[i][Ny-1][k].normal = new glTVector(0,-1,0);
		}
}
//*************************************
void TTwoParallelepipedsWithPressureMaskGenerator::CreateMultyNodes()
{
	// Узлы для U
	for(int i=1; i<NxU-1; ++i)
		XU[i] = (X[i]+X[i-1])/2;
	XU[0] = X[0] - (X[1]-X[0])/2;
	XU[NxU-1] = X[NxU-2] + (X[NxU-2]-X[NxU-3])/2;
	for(int i=1; i<NxU; ++i)
		HxU[i] = XU[i]-XU[i-1];
	HxU[0] = HxU[1];
	for(int j=0; j<NyU; ++j)
	{
		YU[j] = Y[j];
		HyU[j] = Hy[j];
	}
	for(int k=0; k<NzU; ++k)
	{
		ZU[k] = Z[k];
		HzU[k] = Hz[k];
	}

	// Узлы для V
	for(int i=0; i<NxV; ++i)
	{
		XV[i] = X[i];
		HxV[i] = Hx[i];
	}
	for(int j=1; j<NyV-1; ++j)
		YV[j] = (Y[j]+Y[j-1])/2;
	YV[0] = Y[0] - (Y[1]-Y[0])/2;
	YV[NyV-1] = Y[NyV-2] + (Y[NyV-2]-Y[NyV-3])/2;
	for(int j=1; j<NyV; ++j)
		HyV[j] = YV[j]-YV[j-1];
	HyV[0] = HyV[1];
	for(int k=0; k<NzV; ++k)
	{
		ZV[k] = Z[k];
		HzV[k] = Hz[k];
	}

	// Узлы для W
	for(int i=0; i<NxW; ++i)
	{
		XW[i] = X[i];
		HxW[i] = Hx[i];
	}
	for(int j=0; j<NyW; ++j)
	{
		YW[j] = Y[j];
		HyW[j] = Hy[j];
	}
	for(int k=1; k<NzW-1; ++k)
		ZW[k] = (Z[k]+Z[k-1])/2;
	ZW[0] = Z[0] - (Z[1]-Z[0])/2;
	ZW[NzW-1] = Z[NzW-2] + (Z[NzW-2]-Z[NzW-3])/2;
	for(int k=1; k<NzW; ++k)
		HzW[k] = ZW[k]-ZW[k-1];
	HzW[0] = HzW[1];
}
//*************************************
void TTwoParallelepipedsWithPressureMaskGenerator::CreateMultyGrids()
{
	// Заполнение маски для U //

	for(int i=0; i<NxU; ++i)
		for(int j=0; j<NyU; ++j)
			for(int k=0; k<NzU; ++k)
			{
				MaskU[i][j][k].mask = TActualPoint;
				MaskU[i][j][k].normal = 0;
			}

	// Границы большого параллелепипеда
	for(int j=1; j<NyU-1; ++j)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TDefinedBorderPoint;
			MaskU[0][j][k].normal = new glTVector(-1,0,0);
			MaskU[NxU-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[i][0][k].mask = MaskU[i][NyU-1][k].mask = TPreDefinedBorderPoint;
			MaskU[i][0][k].normal =  new glTVector(0,-1,0);
			MaskU[i][NyU-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int j=1; j<NyU-1; ++j)
		{
			MaskU[i][j][0].mask = MaskU[i][j][NzU-1].mask = TPreDefinedBorderPoint;
			MaskU[i][j][0].normal = new glTVector(0,0,-1);
			MaskU[i][j][NzU-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxU; ++i)
	{
		MaskU[i][0][0].mask = TFictivePoint;
		MaskU[i][0][NzU-1].mask = TFictivePoint;
		MaskU[i][NyU-1][0].mask = TFictivePoint;
		MaskU[i][NyU-1][NzU-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyU; ++j)
	{
		MaskU[0][j][0].mask = TFictivePoint;
		MaskU[0][j][NzU-1].mask = TFictivePoint;
		MaskU[NxU-1][j][0].mask = TFictivePoint;
		MaskU[NxU-1][j][NzU-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzU; ++k)
	{
		MaskU[0][0][k].mask = TFictivePoint;
		MaskU[0][NyU-1][k].mask = TFictivePoint;
		MaskU[NxU-1][0][k].mask = TFictivePoint;
		MaskU[NxU-1][NyU-1][k].mask = TFictivePoint;
	}

	// Внутренний параллелепипед
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
			for(int k=0; k<NzU; ++k)
				MaskU[i][j][k].mask = TFictivePoint; // внутренность

	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[i][Ny1][k].mask = MaskU[i][Ny1+Ny2-2][k].mask = TPreDefinedBorderPoint;
			MaskU[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskU[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[Nx1-1][j][k].mask = MaskU[Nx1+Nx2-2][j][k].mask = TDefinedBorderPoint;
			MaskU[Nx1-1][j][k].normal = new glTVector(1,0,0);
			MaskU[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	// Отверстия на внутреннем цилиндре
	for(int j=Ny4; j<Ny4+Ny5-1; ++j)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskU[Nx1-1][j][k].mask = MaskU[Nx1+Nx2-2][j][k].mask = TNormalBorderPoint;
			MaskU[Nx1-1][j][k].normal = new glTVector(1,0,0);
			MaskU[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx4-1; i<Nx4+Nx5-1; ++i)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskU[i][Ny1][k].mask = MaskU[i][Ny1+Ny2-2][k].mask = TPreDefinedBorderPoint;
			MaskU[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskU[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	// Отверстия на внешнем цилиндре
	for(int j=Ny6; j<Ny6+Ny7-1; ++j)
		for(int k=Nz3; k<Nz3+Nz4-1; ++k)
		{
			MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TNormalBorderPoint;
			MaskU[0][j][k].normal = new glTVector(1,0,0);
			MaskU[NxU-1][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx6-1; i<Nx6+Nx7-1; ++i)
		for(int k=Nz3; k<Nz3+Nz4-1; ++k)
		{
			MaskU[i][0][k].mask = MaskU[i][NyU-1][k].mask = TPreDefinedBorderPoint;
			MaskU[i][0][k].normal = new glTVector(0,1,0);
			MaskU[i][NyU-1][k].normal = new glTVector(0,-1,0);
		}

	// Заполнение маски для V //

	for(int i=0; i<NxV; ++i)
		for(int j=0; j<NyV; ++j)
			for(int k=0; k<NzV; ++k)
			{
				MaskV[i][j][k].mask = TActualPoint;
				MaskV[i][j][k].normal = 0;
			}

	// Границы большого параллелепипеда
	for(int j=1; j<NyV-1; ++j)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[0][j][k].mask = MaskV[NxV-1][j][k].mask = TPreDefinedBorderPoint;
			MaskV[0][j][k].normal = new glTVector(-1,0,0);
			MaskV[NxV-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[i][0][k].mask = MaskV[i][NyV-1][k].mask = TDefinedBorderPoint;
			MaskV[i][0][k].normal =  new glTVector(0,-1,0);
			MaskV[i][NyV-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int j=1; j<NyV-1; ++j)
		{
			MaskV[i][j][0].mask = MaskV[i][j][NzV-1].mask = TPreDefinedBorderPoint;
			MaskV[i][j][0].normal = new glTVector(0,0,-1);
			MaskV[i][j][NzV-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxV; ++i)
	{
		MaskV[i][0][0].mask = TFictivePoint;
		MaskV[i][0][NzV-1].mask = TFictivePoint;
		MaskV[i][NyV-1][0].mask = TFictivePoint;
		MaskV[i][NyV-1][NzV-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyV; ++j)
	{
		MaskV[0][j][0].mask = TFictivePoint;
		MaskV[0][j][NzV-1].mask = TFictivePoint;
		MaskV[NxV-1][j][0].mask = TFictivePoint;
		MaskV[NxV-1][j][NzV-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzV; ++k)
	{
		MaskV[0][0][k].mask = TFictivePoint;
		MaskV[0][NyV-1][k].mask = TFictivePoint;
		MaskV[NxV-1][0][k].mask = TFictivePoint;
		MaskV[NxV-1][NyV-1][k].mask = TFictivePoint;
	}

	// Внутренний параллелепипед
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
			for(int k=0; k<NzV; ++k)
				MaskV[i][j][k].mask = TFictivePoint; // внутренность

	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[Nx1][j][k].mask = MaskV[Nx1+Nx2-2][j][k].mask = TPreDefinedBorderPoint;
			MaskV[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskV[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[i][Ny1-1][k].mask = MaskV[i][Ny1+Ny2-2][k].mask = TDefinedBorderPoint;
			MaskV[i][Ny1-1][k].normal = new glTVector(0,1,0);
			MaskV[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	// Отверстия на внутреннем цилиндре
	for(int j=Ny4-1; j<Ny4+Ny5-1; ++j)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskV[Nx1][j][k].mask = MaskV[Nx1+Nx2-2][j][k].mask = TPreDefinedBorderPoint;
			MaskV[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskV[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx4; i<Nx4+Nx5-1; ++i)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskV[i][Ny1-1][k].mask = MaskV[i][Ny1+Ny2-2][k].mask = TNormalBorderPoint;
			MaskV[i][Ny1-1][k].normal = new glTVector(0,1,0);
			MaskV[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	// Отверстия на внешнем цилиндре
	for(int j=Ny6-1; j<Ny6+Ny7-1; ++j)
		for(int k=Nz3; k<Nz3+Nz4-1; ++k)
		{
			MaskV[0][j][k].mask = MaskV[NxV-1][j][k].mask = TPreDefinedBorderPoint;
			MaskV[0][j][k].normal = new glTVector(1,0,0);
			MaskV[NxV-1][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx6; i<Nx6+Nx7-1; ++i)
		for(int k=Nz3; k<Nz3+Nz4-1; ++k)
		{
			MaskV[i][0][k].mask = MaskV[i][NyV-1][k].mask = TNormalBorderPoint;
			MaskV[i][0][k].normal = new glTVector(0,1,0);
			MaskV[i][NyV-1][k].normal = new glTVector(0,-1,0);
		}

	// Заполнение маски для W

	for(int i=0; i<NxW; ++i)
		for(int j=0; j<NyW; ++j)
			for(int k=0; k<NzW; ++k)
			{
				MaskW[i][j][k].mask = TActualPoint;
				MaskW[i][j][k].normal = 0;
			}

	// Границы большого параллелепипеда
	for(int j=1; j<NyW-1; ++j)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[0][j][k].mask = MaskW[NxW-1][j][k].mask = TPreDefinedBorderPoint; 
			MaskW[0][j][k].normal = new glTVector(-1,0,0);
			MaskW[NxW-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[i][0][k].mask = MaskW[i][NyW-1][k].mask = TPreDefinedBorderPoint;
			MaskW[i][0][k].normal =  new glTVector(0,-1,0);
			MaskW[i][NyW-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int j=1; j<NyW-1; ++j)
		{
			MaskW[i][j][0].mask = MaskW[i][j][NzW-1].mask = TDefinedBorderPoint;
			MaskW[i][j][0].normal = new glTVector(0,0,-1);
			MaskW[i][j][NzW-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxW; ++i)
	{
		MaskW[i][0][0].mask = TFictivePoint;
		MaskW[i][0][NzW-1].mask = TFictivePoint;
		MaskW[i][NyW-1][0].mask = TFictivePoint;
		MaskW[i][NyW-1][NzW-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyW; ++j)
	{
		MaskW[0][j][0].mask = TFictivePoint;
		MaskW[0][j][NzW-1].mask = TFictivePoint;
		MaskW[NxW-1][j][0].mask = TFictivePoint;
		MaskW[NxW-1][j][NzW-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzW; ++k)
	{
		MaskW[0][0][k].mask = TFictivePoint;
		MaskW[0][NyW-1][k].mask = TFictivePoint;
		MaskW[NxW-1][0][k].mask = TFictivePoint;
		MaskW[NxW-1][NyW-1][k].mask = TFictivePoint;
	}

	// Внутренний параллелепипед
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
			for(int k=0; k<NzW; ++k)
				MaskW[i][j][k].mask = TFictivePoint; // внутренность

	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[Nx1][j][k].mask = MaskW[Nx1+Nx2-2][j][k].mask = TPreDefinedBorderPoint;
			MaskW[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskW[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[i][Ny1][k].mask = MaskW[i][Ny1+Ny2-2][k].mask = TPreDefinedBorderPoint;
			MaskW[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskW[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	// Отверстия на внутреннем цилиндре
	for(int j=Ny4; j<Ny4+Ny5-1; ++j)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskW[Nx1][j][k].mask = MaskW[Nx1+Nx2-2][j][k].mask = TPreDefinedBorderPoint;
			MaskW[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskW[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx4; i<Nx4+Nx5-1; ++i)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskW[i][Ny1][k].mask = MaskW[i][Ny1+Ny2-2][k].mask = TPreDefinedBorderPoint;
			MaskW[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskW[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	// Отверстия на внешнем цилиндре
	for(int j=Ny6; j<Ny6+Ny7-1; ++j)
		for(int k=Nz3; k<Nz3+Nz4-1; ++k)
		{
			MaskW[0][j][k].mask = MaskW[NxW-1][j][k].mask = TPreDefinedBorderPoint;
			MaskW[0][j][k].normal = new glTVector(1,0,0);
			MaskW[NxW-1][j][k].normal = new glTVector(-1,0,0);
	}
		
	for(int i=Nx6; i<Nx6+Nx7-1; ++i)
		for(int k=Nz3; k<Nz3+Nz4-1; ++k)
		{
			MaskW[i][0][k].mask = MaskW[i][NyW-1][k].mask = TPreDefinedBorderPoint;
			MaskW[i][0][k].normal = new glTVector(0,1,0);
			MaskW[i][NyW-1][k].normal = new glTVector(0,-1,0);
		}

	// Заполнение маски для P
	for(int i=0; i<NxP; ++i)
		for(int j=0; j<NyP; ++j)
			for(int k=0; k<NzP; ++k)
			{
				MaskP[i][j][k].mask = TActualPoint;
				MaskP[i][j][k].normal = 0;
			}

	// Границы большого параллелепипеда
	for(int j=1; j<NyP-1; ++j)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreNormalBorderPoint; 
			MaskP[0][j][k].normal = new glTVector(-1,0,0);
			MaskP[NxP-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[i][0][k].mask = MaskP[i][NyP-1][k].mask = TPreNormalBorderPoint;
			MaskP[i][0][k].normal =  new glTVector(0,-1,0);
			MaskP[i][NyP-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int j=1; j<NyP-1; ++j)
		{
			MaskP[i][j][0].mask = MaskP[i][j][NzP-1].mask = TPreNormalBorderPoint;
			MaskP[i][j][0].normal = new glTVector(0,0,-1);
			MaskP[i][j][NzP-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxP; ++i)
	{
		MaskP[i][0][0].mask = TFictivePoint;
		MaskP[i][0][NzP-1].mask = TFictivePoint;
		MaskP[i][NyP-1][0].mask = TFictivePoint;
		MaskP[i][NyP-1][NzP-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyP; ++j)
	{
		MaskP[0][j][0].mask = TFictivePoint;
		MaskP[0][j][NzP-1].mask = TFictivePoint;
		MaskP[NxP-1][j][0].mask = TFictivePoint;
		MaskP[NxP-1][j][NzP-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzP; ++k)
	{
		MaskP[0][0][k].mask = TFictivePoint;
		MaskP[0][NyP-1][k].mask = TFictivePoint;
		MaskP[NxP-1][0][k].mask = TFictivePoint;
		MaskP[NxP-1][NyP-1][k].mask = TFictivePoint;
	}

	// Внутренний параллелепипед
	for(int i=Nx1+1; i<Nx1+Nx2-2; ++i)
		for(int j=Ny1+1; j<Ny1+Ny2-2; ++j)
			for(int k=0; k<NzP; ++k)
				MaskP[i][j][k].mask = TFictivePoint; // внутренность

	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		for(int k=0; k<NzP; ++k)
		{
			MaskP[Nx1][j][k].mask = MaskP[Nx1+Nx2-2][j][k].mask = TPreNormalBorderPoint;
			MaskP[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskP[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int k=0; k<NzP; ++k)
		{
			MaskP[i][Ny1][k].mask = MaskP[i][Ny1+Ny2-2][k].mask = TPreNormalBorderPoint;
			MaskP[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskP[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	// Ребра препятствия
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	{
		MaskP[i][Ny1][0].mask = TFictivePoint;
		MaskP[i][Ny1+Ny2-2][0].mask = TFictivePoint;
		MaskP[i][Ny1][NzP-1].mask = TFictivePoint;
		MaskP[i][Ny1+Ny2-2][NzP-1].mask = TFictivePoint;
	}

	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	{
		MaskP[Nx1][j][Nz1].mask = TFictivePoint;
		MaskP[Nx1+Nx2-2][j][Nz1].mask = TFictivePoint;
		MaskP[Nx1][j][NzP-1].mask = TFictivePoint;
		MaskP[Nx1+Nx2-2][j][NzP-1].mask = TFictivePoint;
	}

	for(int k=1; k<NzP-1; ++k)
	{
		MaskP[Nx1][Ny1][k].normal = new glTVector(0,1,1);
		MaskP[Nx1+Nx2-2][Ny1][k].normal = new glTVector(0,-1,1);
		MaskP[Nx1][Ny1+Ny2-2][k].normal = new glTVector(0,1,-1);
		MaskP[Nx1+Nx2-2][Ny1+Ny2-2][k].normal = new glTVector(0,-1,-1);
	}

	// Отверстия на внутреннем цилиндре
	for(int j=Ny4; j<Ny4+Ny5-1; ++j)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskP[Nx1][j][k].mask = MaskP[Nx1+Nx2-2][j][k].mask = TPreDefinedBorderPoint;
			MaskP[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskP[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx4; i<Nx4+Nx5-1; ++i)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskP[i][Ny1][k].mask = MaskP[i][Ny1+Ny2-2][k].mask = TPreDefinedBorderPoint;
			MaskP[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskP[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	// Отверстия на внешнем цилиндре
	for(int j=Ny6; j<Ny6+Ny7-1; ++j)
		for(int k=Nz3; k<Nz3+Nz4-1; ++k)
		{
			MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreDefinedBorderPoint;
			MaskP[0][j][k].normal = new glTVector(1,0,0);
			MaskP[NxP-1][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx6; i<Nx6+Nx7-1; ++i)
		for(int k=Nz3; k<Nz3+Nz4-1; ++k)
		{
			MaskP[i][0][k].mask = MaskP[i][NyP-1][k].mask = TPreDefinedBorderPoint;
			MaskP[i][0][k].normal = new glTVector(0,1,0);
			MaskP[i][NyP-1][k].normal = new glTVector(0,-1,0);
		}
}
//*************************************
double TTwoParallelepipedsWithPressureMaskGenerator::getUCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TTwoParallelepipedsWithPressureMaskGenerator::getVCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TTwoParallelepipedsWithPressureMaskGenerator::getWCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TTwoParallelepipedsWithPressureMaskGenerator::getPCondition(int i, int j, int k, int n)
{
	if((k>=Nz1)&&(k<Nz1+Nz2-1))
		return PInner;
	else if((k>=Nz3)&&(k<Nz3+Nz4-1))
		return POuter;
	else
		return 0;
}
//*************************************



//#################################################################



//*************************************
void TCubeBottomWithPressureMaskGenerator::CreateGeometryNodes()
{
	// Узлы, шаги сетки
	for(int i=0; i<Nx; ++i)
		X[i] = LengthX*i/(Nx-1);
	for(int i=1; i<Nx; ++i)
		Hx[i] = X[i]-X[i-1];
	Hx[0] = Hx[1];

	for(int j=0; j<Ny; ++j)
		Y[j] = LengthY*j/(Ny-1);
	for(int j=1; j<Ny; ++j)
		Hy[j] = Y[j]-Y[j-1];
	Hy[0] = Hy[1];

	for(int k=0; k<Nz; ++k)
		Z[k] = LengthZ*k/(Nz-1);
	for(int k=1; k<Nz; ++k)
		Hz[k] = Z[k]-Z[k-1];
	Hy[0] = Hy[1];
}
//*************************************
void TCubeBottomWithPressureMaskGenerator::CreateGeometryGrid()
{
	// Сначала все точки - расчетные
	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz; ++k)
			{
				Mask[i][j][k].mask = TActualPoint;
				Mask[i][j][k].normal = 0;
			}

	// Левая и правая границы
	for(int j=1; j<Ny-1; ++j)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[Nx-1][j][k].mask = Mask[0][j][k].mask = TBorderPoint;
			Mask[0][j][k].normal = new glTVector(-1,0,0);
			Mask[Nx-1][j][k].normal = new glTVector(1,0,0);
		}

	// Ближняя и дальняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[i][Ny-1][k].mask = Mask[i][0][k].mask = TBorderPoint;
			Mask[i][0][k].normal = new glTVector(0,-1,0);
			Mask[i][Ny-1][k].normal = new glTVector(0,1,0);
		}

	// Нижняя и верхняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int j=1; j<Ny-1; ++j)
		{
			Mask[i][j][Nz-1].mask = Mask[i][j][0].mask = TBorderPoint;
			Mask[i][j][0].normal = new glTVector(0,0,-1);
			Mask[i][j][Nz-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<Nx; ++i)
	{
		Mask[i][0][0].mask = TFictivePoint;
		Mask[i][0][Nz-1].mask = TFictivePoint;
		Mask[i][Ny-1][0].mask = TFictivePoint;
		Mask[i][Ny-1][Nz-1].mask = TFictivePoint;
	}

	for(int j=0; j<Ny; ++j)
	{
		Mask[0][j][0].mask = TFictivePoint;
		Mask[0][j][Nz-1].mask = TFictivePoint;
		Mask[Nx-1][j][0].mask = TFictivePoint;
		Mask[Nx-1][j][Nz-1].mask = TFictivePoint;
	}

	for(int k=0; k<Nz; ++k)
	{
		Mask[0][0][k].mask = TFictivePoint;
		Mask[0][Ny-1][k].mask = TFictivePoint;
		Mask[Nx-1][0][k].mask = TFictivePoint;
		Mask[Nx-1][Ny-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1; i<Nx1+Nx2-2; ++i)
		for(int j=Ny1; j<Ny1+Ny2-2; ++j)
			for(int k=0; k<Nz1-1; ++k)
				Mask[i][j][k].mask = TFictivePoint;
	
	// Грани препятствия
	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		for(int k=0; k<Nz1; ++k)
		{
			Mask[Nx1-1][j][k].mask = Mask[Nx1+Nx2-2][j][k].mask = TBorderPoint;
			Mask[Nx1-1][j][k].normal = new glTVector(1,0,0);
			Mask[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int k=0; k<Nz1; ++k)
		{
			Mask[i][Ny1-1][k].mask = Mask[i][Ny1+Ny2-2][k].mask = TBorderPoint;
			Mask[i][Ny1-1][k].normal = new glTVector(0,1,0);
			Mask[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		{
			Mask[i][j][Nz1-1].mask = TBorderPoint;
			Mask[i][j][Nz1-1].normal = new glTVector(0,0,-1);
		}

	// Ребра препятствия
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
	{
		Mask[i][Ny1-1][0].mask = TFictivePoint;
		Mask[i][Ny1+Ny2-2][0].mask = TFictivePoint;
		Mask[i][Ny1-1][Nz1-1].normal->set(1,-1,0);
		Mask[i][Ny1+Ny2-2][Nz1-1].normal->set(-1,-1,0);
	}

	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
	{
		Mask[Nx1-1][j][0].mask = TFictivePoint;
		Mask[Nx1+Nx2-2][j][0].mask = TFictivePoint;
		Mask[Nx1-1][j][Nz1-1].normal->set(1,0,-1);
		Mask[Nx1+Nx2-2][j][Nz1-1].normal->set(-1,0,-1);
	}

	for(int k=0; k<Nz1; ++k)
	{
		Mask[Nx1-1][Ny1-1][k].normal->set(0,1,1);
		Mask[Nx1+Nx2-2][Ny1-1][k].normal->set(0,-1,1);
		Mask[Nx1-1][Ny1+Ny2-2][k].normal->set(0,1,-1);
		Mask[Nx1+Nx2-2][Ny1+Ny2-2][k].normal->set(0,-1,-1);
	}

	// Угловые точки препятствия
	Mask[Nx1-1][Ny1-1][Nz1-1].normal->set(1,1,-1);
	Mask[Nx1+Nx2-2][Ny1-1][Nz1-1].normal->set(-1,1,-1);
	Mask[Nx1-1][Ny1+Ny2-2][Nz1-1].normal->set(1,-1,-1);
	Mask[Nx1+Nx2-2][Ny1+Ny2-2][Nz1-1].normal->set(-1,-1,-1);
}
//*************************************
void TCubeBottomWithPressureMaskGenerator::CreateMultyNodes()
{
	// Узлы для компоненты U скорости
	for(int i=0; i<NxU; ++i)
		XU[i] = X[i];

	for(int j=1; j<NyU-1; ++j)
		YU[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YU[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YU[NyU-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzU-1; ++k)
		ZU[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZU[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZU[NzU-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxU; ++i)
		HxU[i] = XU[i] - XU[i-1];
	HxU[0] = HxU[1];

	for(int j=1; j<NyU; ++j)
		HyU[j] = YU[j] - YU[j-1];
	HyU[0] = HyU[1];

	for(int k=1; k<NzU; ++k)
		HzU[k] = ZU[k] - ZU[k-1];
	HzU[0] = HzU[1];

	// Узлы для компоненты V скорости
	for(int i=1; i<NxV-1; ++i)
		XV[i] = ( X[i-1] + X[i] ) / 2.0;
	XV[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XV[NxV-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=0; j<NyV; ++j)
		YV[j] = Y[j];

	for(int k=1; k<NzV-1; ++k)
		ZV[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZV[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZV[NzV-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxV; ++i)
		HxV[i] = XV[i] - XV[i-1];
	HxV[0] = HxV[1];

	for(int j=1; j<NyV; ++j)
		HyV[j] = YV[j] - YV[j-1];
	HyV[0] = HyV[1];

	for(int k=1; k<NzV; ++k)
		HzV[k] = ZV[k] - ZV[k-1];
	HzV[0] = HzV[1];

	// Узлы для компоненты W скорости
	for(int i=1; i<NxW-1; ++i)
		XW[i] = ( X[i-1] + X[i] ) / 2.0;
	XW[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XW[NxW-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyW-1; ++j)
		YW[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YW[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YW[NyW-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=0; k<NzW; ++k)
		ZW[k] = Z[k];

	for(int i=1; i<NxW; ++i)
		HxW[i] = XW[i] - XW[i-1];
	HxW[0] = HxW[1];

	for(int j=1; j<NyW; ++j)
		HyW[j] = YW[j] - YW[j-1];
	HyW[0] = HyW[1];

	for(int k=1; k<NzW; ++k)
		HzW[k] = ZW[k] - ZW[k-1];
	HzW[0] = HzW[1];

	// Узлы для давления
	for(int i=1; i<NxP-1; ++i)
		XP[i] = ( X[i-1] + X[i] ) / 2.0;
	XP[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XP[NxP-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyP-1; ++j)
		YP[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YP[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YP[NyP-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzP-1; ++k)
		ZP[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZP[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZP[NzP-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxP; ++i)
		HxP[i] = XP[i] - XP[i-1];
	HxP[0] = HxP[1];

	for(int j=1; j<NyP; ++j)
		HyP[j] = YP[j] - YP[j-1];
	HyP[0] = HyP[1];

	for(int k=1; k<NzP; ++k)
		HzP[k] = ZP[k] - ZP[k-1];
	HzP[0] = HzP[1];
}
//*************************************
void TCubeBottomWithPressureMaskGenerator::CreateMultyGrids()
{
	// Заполнение маски для U //

	for(int i=0; i<NxU; ++i)
		for(int j=0; j<NyU; ++j)
			for(int k=0; k<NzU; ++k)
			{
				MaskU[i][j][k].mask = TActualPoint;
				MaskU[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyU-1; ++j)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TDefinedBorderPoint;
			MaskU[0][j][k].normal = new glTVector(-1,0,0);
			MaskU[NxU-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int j=2; j<NyU-2; ++j)
		for(int k=2; k<NzU-2; ++k)
		{
			MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TNormalBorderPoint;
			MaskU[0][j][k].normal = new glTVector(-1,0,0);
			MaskU[NxU-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[i][0][k].mask = MaskU[i][NyU-1][k].mask = TPreDefinedBorderPoint;
			MaskU[i][0][k].normal =  new glTVector(0,-1,0);
			MaskU[i][NyU-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int j=1; j<NyU-1; ++j)
		{
			MaskU[i][j][0].mask = MaskU[i][j][NzU-1].mask = TPreDefinedBorderPoint;
			MaskU[i][j][0].normal = new glTVector(0,0,-1);
			MaskU[i][j][NzU-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxU; ++i)
	{
		MaskU[i][0][0].mask = TFictivePoint;
		MaskU[i][0][NzU-1].mask = TFictivePoint;
		MaskU[i][NyU-1][0].mask = TFictivePoint;
		MaskU[i][NyU-1][NzU-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyU; ++j)
	{
		MaskU[0][j][0].mask = TFictivePoint;
		MaskU[0][j][NzU-1].mask = TFictivePoint;
		MaskU[NxU-1][j][0].mask = TFictivePoint;
		MaskU[NxU-1][j][NzU-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzU; ++k)
	{
		MaskU[0][0][k].mask = TFictivePoint;
		MaskU[0][NyU-1][k].mask = TFictivePoint;
		MaskU[NxU-1][0][k].mask = TFictivePoint;
		MaskU[NxU-1][NyU-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1; i<Nx1+Nx2-2; ++i)
		for(int j=Ny1+1; j<Ny1+Ny2-2; ++j)
			for(int k=0; k<Nz1-1; ++k)
				MaskU[i][j][k].mask = TFictivePoint;

	// Грани препятствия
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int k=0; k<Nz1; ++k)
		{
			MaskU[i][Ny1][k].mask = MaskU[i][Ny1+Ny2-2][k].mask = TPreDefinedBorderPoint;
			MaskU[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskU[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			MaskU[i][j][Nz1-1].mask = TPreDefinedBorderPoint;
			MaskU[i][j][Nz1-1].normal = new glTVector(0,0,-1);
		}

	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		for(int k=0; k<Nz1; ++k)
		{
			MaskU[Nx1-1][j][k].mask = MaskU[Nx1+Nx2-2][j][k].mask = TDefinedBorderPoint;
			MaskU[Nx1-1][j][k].normal = new glTVector(1,0,0);
			MaskU[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}

	// Заполнение маски для V //

	for(int i=0; i<NxV; ++i)
		for(int j=0; j<NyV; ++j)
			for(int k=0; k<NzV; ++k)
			{
				MaskV[i][j][k].mask = TActualPoint;
				MaskV[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyV-1; ++j)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[0][j][k].mask = MaskV[NxV-1][j][k].mask = TPreDefinedBorderPoint;
			MaskV[0][j][k].normal = new glTVector(-1,0,0);
			MaskV[NxV-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[i][0][k].mask = MaskV[i][NyV-1][k].mask = TDefinedBorderPoint;
			MaskV[i][0][k].normal =  new glTVector(0,-1,0);
			MaskV[i][NyV-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int j=1; j<NyV-1; ++j)
		{
			MaskV[i][j][0].mask = MaskV[i][j][NzV-1].mask = TPreDefinedBorderPoint;
			MaskV[i][j][0].normal = new glTVector(0,0,-1);
			MaskV[i][j][NzV-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxV; ++i)
	{
		MaskV[i][0][0].mask = TFictivePoint;
		MaskV[i][0][NzV-1].mask = TFictivePoint;
		MaskV[i][NyV-1][0].mask = TFictivePoint;
		MaskV[i][NyV-1][NzV-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyV; ++j)
	{
		MaskV[0][j][0].mask = TFictivePoint;
		MaskV[0][j][NzV-1].mask = TFictivePoint;
		MaskV[NxV-1][j][0].mask = TFictivePoint;
		MaskV[NxV-1][j][NzV-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzV; ++k)
	{
		MaskV[0][0][k].mask = TFictivePoint;
		MaskV[0][NyV-1][k].mask = TFictivePoint;
		MaskV[NxV-1][0][k].mask = TFictivePoint;
		MaskV[NxV-1][NyV-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1+1; i<Nx1+Nx2-2; ++i)
		for(int j=Ny1; j<Ny1+Ny2-2; ++j)
			for(int k=0; k<Nz1-1; ++k)
				MaskV[i][j][k].mask = TFictivePoint;
	
	// Грани препятствия
	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		for(int k=0; k<Nz1; ++k)
		{
			MaskV[Nx1][j][k].mask = MaskV[Nx1+Nx2-2][j][k].mask = TPreDefinedBorderPoint;
			MaskV[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskV[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		{
			MaskV[i][j][Nz1-1].mask = TPreDefinedBorderPoint;
			MaskV[i][j][Nz1-1].normal = new glTVector(0,0,-1);
		}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int k=0; k<Nz1; ++k)
		{
			MaskV[i][Ny1-1][k].mask = MaskV[i][Ny1+Ny2-2][k].mask = TDefinedBorderPoint;
			MaskV[i][Ny1-1][k].normal = new glTVector(0,1,0);
			MaskV[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	// Заполнение маски для W //

	for(int i=0; i<NxW; ++i)
		for(int j=0; j<NyW; ++j)
			for(int k=0; k<NzW; ++k)
			{
				MaskW[i][j][k].mask = TActualPoint;
				MaskW[i][j][k].normal = 0;
			}
	// Границы области
	for(int j=1; j<NyW-1; ++j)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[0][j][k].mask = MaskW[NxW-1][j][k].mask = TPreDefinedBorderPoint; 
			MaskW[0][j][k].normal = new glTVector(-1,0,0);
			MaskW[NxW-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[i][0][k].mask = MaskW[i][NyW-1][k].mask = TPreDefinedBorderPoint;
			MaskW[i][0][k].normal =  new glTVector(0,-1,0);
			MaskW[i][NyW-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int j=1; j<NyW-1; ++j)
		{
			MaskW[i][j][0].mask = MaskW[i][j][NzW-1].mask = TDefinedBorderPoint;
			MaskW[i][j][0].normal = new glTVector(0,0,-1);
			MaskW[i][j][NzW-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxW; ++i)
	{
		MaskW[i][0][0].mask = TFictivePoint;
		MaskW[i][0][NzW-1].mask = TFictivePoint;
		MaskW[i][NyW-1][0].mask = TFictivePoint;
		MaskW[i][NyW-1][NzW-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyW; ++j)
	{
		MaskW[0][j][0].mask = TFictivePoint;
		MaskW[0][j][NzW-1].mask = TFictivePoint;
		MaskW[NxW-1][j][0].mask = TFictivePoint;
		MaskW[NxW-1][j][NzW-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzW; ++k)
	{
		MaskW[0][0][k].mask = TFictivePoint;
		MaskW[0][NyW-1][k].mask = TFictivePoint;
		MaskW[NxW-1][0][k].mask = TFictivePoint;
		MaskW[NxW-1][NyW-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1+1; i<Nx1+Nx2-2; ++i)
		for(int j=Ny1+1; j<Ny1+Ny2-2; ++j)
			for(int k=0; k<Nz1-1; ++k)
				MaskW[i][j][k].mask = TFictivePoint;

	// Грани препятствия
	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		for(int k=0; k<Nz1; ++k)
		{
			MaskW[Nx1][j][k].mask = MaskW[Nx1+Nx2-2][j][k].mask = TPreDefinedBorderPoint;
			MaskW[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskW[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int k=0; k<Nz1; ++k)
		{
			MaskW[i][Ny1][k].mask = MaskW[i][Ny1+Ny2-2][k].mask = TPreDefinedBorderPoint;
			MaskW[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskW[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			MaskW[i][j][Nz1-1].mask = TDefinedBorderPoint;
			MaskW[i][j][Nz1-1].normal = new glTVector(0,0,-1);
		}

	// Заполнение маски для P
	for(int i=0; i<NxP; ++i)
		for(int j=0; j<NyP; ++j)
			for(int k=0; k<NzP; ++k)
			{
				MaskP[i][j][k].mask = TActualPoint;
				MaskP[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyP-1; ++j)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreNormalBorderPoint;
			MaskP[0][j][k].normal = new glTVector(-1,0,0);
			MaskP[NxP-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int j=2; j<NyP-2; ++j)
		for(int k=2; k<NzP-2; ++k)
		{
			MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreDefinedBorderPoint; 
			//MaskP[0][j][k].normal = new glTVector(-1,0,0);
			//MaskP[NxP-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[i][0][k].mask = MaskP[i][NyP-1][k].mask = TPreNormalBorderPoint;
			MaskP[i][0][k].normal =  new glTVector(0,-1,0);
			MaskP[i][NyP-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int j=1; j<NyP-1; ++j)
		{
			MaskP[i][j][0].mask = MaskP[i][j][NzP-1].mask = TPreNormalBorderPoint;
			MaskP[i][j][0].normal = new glTVector(0,0,-1);
			MaskP[i][j][NzP-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxP; ++i)
	{
		MaskP[i][0][0].mask = TFictivePoint;
		MaskP[i][0][NzP-1].mask = TFictivePoint;
		MaskP[i][NyP-1][0].mask = TFictivePoint;
		MaskP[i][NyP-1][NzP-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyP; ++j)
	{
		MaskP[0][j][0].mask = TFictivePoint;
		MaskP[0][j][NzP-1].mask = TFictivePoint;
		MaskP[NxP-1][j][0].mask = TFictivePoint;
		MaskP[NxP-1][j][NzP-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzP; ++k)
	{
		MaskP[0][0][k].mask = TFictivePoint;
		MaskP[0][NyP-1][k].mask = TFictivePoint;
		MaskP[NxP-1][0][k].mask = TFictivePoint;
		MaskP[NxP-1][NyP-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1+1; i<Nx1+Nx2-2; ++i)
		for(int j=Ny1+1; j<Ny1+Ny2-2; ++j)
			for(int k=0; k<Nz1-1; ++k)
				MaskP[i][j][k].mask = TFictivePoint;
	
	// Грани препятствия
	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		for(int k=0; k<Nz1; ++k)
		{
			MaskP[Nx1][j][k].mask = MaskP[Nx1+Nx2-2][j][k].mask = TPreNormalBorderPoint;
			MaskP[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskP[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int k=0; k<Nz1; ++k)
		{
			MaskP[i][Ny1][k].mask = MaskP[i][Ny1+Ny2-2][k].mask = TPreNormalBorderPoint;
			MaskP[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskP[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			MaskP[i][j][Nz1-1].mask = TPreNormalBorderPoint;
			MaskP[i][j][Nz1-1].normal = new glTVector(0,0,-1);
		}

	// Ребра препятствия
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	{
		MaskP[i][Ny1][0].mask = TFictivePoint;
		MaskP[i][Ny1+Ny2-2][0].mask = TFictivePoint;
		MaskP[i][Ny1][Nz1-1].normal->set(1,-1,0);
		MaskP[i][Ny1+Ny2-2][Nz1-1].normal->set(-1,-1,0);
	}

	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	{
		MaskP[Nx1][j][0].mask = TFictivePoint;
		MaskP[Nx1+Nx2-2][j][0].mask = TFictivePoint;
		MaskP[Nx1][j][Nz1-1].normal->set(1,0,-1);
		MaskP[Nx1+Nx2-2][j][Nz1-1].normal->set(-1,0,-1);
	}

	for(int k=0; k<Nz1; ++k)
	{
		MaskP[Nx1][Ny1][k].normal->set(0,1,1);
		MaskP[Nx1+Nx2-2][Ny1][k].normal->set(0,-1,1);
		MaskP[Nx1][Ny1+Ny2-2][k].normal->set(0,1,-1);
		MaskP[Nx1+Nx2-2][Ny1+Ny2-2][k].normal->set(0,-1,-1);
	}

	// Угловые точки препятствия
	MaskP[Nx1][Ny1][Nz1-1].normal->set(1,1,-1);
	MaskP[Nx1+Nx2-2][Ny1][Nz1-1].normal->set(-1,1,-1);
	MaskP[Nx1][Ny1+Ny2-2][Nz1-1].normal->set(1,-1,-1);
	MaskP[Nx1+Nx2-2][Ny1+Ny2-2][Nz1-1].normal->set(-1,-1,-1);
}
//*************************************
double TCubeBottomWithPressureMaskGenerator::getUCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TCubeBottomWithPressureMaskGenerator::getVCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TCubeBottomWithPressureMaskGenerator::getWCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TCubeBottomWithPressureMaskGenerator::getPCondition(int i, int j, int k, int n)
{
	if(i==0)
		return PInner;
	else if(i==NxP-1)
		return POuter;
	else
		return 0;
}
//*************************************




//#################################################################




//*************************************
void TTrapezeBottomWithPressureMaskGenerator::CreateGeometryNodes()
{
	// Узлы, шаги сетки
	for(int i=0; i<Nx; ++i)
		X[i] = LengthX*i/(Nx-1);
	for(int i=1; i<Nx; ++i)
		Hx[i] = X[i]-X[i-1];
	Hx[0] = Hx[1];

	for(int j=0; j<Ny; ++j)
		Y[j] = LengthY*j/(Ny-1);
	for(int j=1; j<Ny; ++j)
		Hy[j] = Y[j]-Y[j-1];
	Hy[0] = Hy[1];

	for(int k=0; k<Nz; ++k)
		Z[k] = LengthZ*k/(Nz-1);
	for(int k=1; k<Nz; ++k)
		Hz[k] = Z[k]-Z[k-1];
	Hy[0] = Hy[1];
}
//*************************************
void TTrapezeBottomWithPressureMaskGenerator::CreateGeometryGrid()
{
	// Сначала все точки - расчетные
	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz; ++k)
			{
				Mask[i][j][k].mask = TActualPoint;
				Mask[i][j][k].normal = 0;
			}

	// Левая и правая границы
	for(int j=1; j<Ny-1; ++j)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[Nx-1][j][k].mask = Mask[0][j][k].mask = TBorderPoint;
			Mask[0][j][k].normal = new glTVector(-1,0,0);
			Mask[Nx-1][j][k].normal = new glTVector(1,0,0);
		}

	// Ближняя и дальняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[i][Ny-1][k].mask = Mask[i][0][k].mask = TBorderPoint;
			Mask[i][0][k].normal = new glTVector(0,-1,0);
			Mask[i][Ny-1][k].normal = new glTVector(0,1,0);
		}

	// Нижняя и верхняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int j=1; j<Ny-1; ++j)
		{
			Mask[i][j][Nz-1].mask = Mask[i][j][0].mask = TBorderPoint;
			Mask[i][j][0].normal = new glTVector(0,0,-1);
			Mask[i][j][Nz-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<Nx; ++i)
	{
		Mask[i][0][0].mask = TFictivePoint;
		Mask[i][0][Nz-1].mask = TFictivePoint;
		Mask[i][Ny-1][0].mask = TFictivePoint;
		Mask[i][Ny-1][Nz-1].mask = TFictivePoint;
	}

	for(int j=0; j<Ny; ++j)
	{
		Mask[0][j][0].mask = TFictivePoint;
		Mask[0][j][Nz-1].mask = TFictivePoint;
		Mask[Nx-1][j][0].mask = TFictivePoint;
		Mask[Nx-1][j][Nz-1].mask = TFictivePoint;
	}

	for(int k=0; k<Nz; ++k)
	{
		Mask[0][0][k].mask = TFictivePoint;
		Mask[0][Ny-1][k].mask = TFictivePoint;
		Mask[Nx-1][0][k].mask = TFictivePoint;
		Mask[Nx-1][Ny-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		{
			int p = i - (Nx1-1), q;
			if(p < Nz2)
				q = p;
			else if(p >= Nx2-Nz2)
				q = Nx2-1-p;
			else
				q = Nz2;

			for(int k=0; k<Nz1-1+q; ++k)
				MaskP[i][j][k].mask = TFictivePoint;
		}

	// Грани препятствия
	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
	{
		for(int k=0; k<Nz2-1; ++k)
		{
			MaskP[(Nx1-1) + k][j][Nz1-1+k].mask = MaskP[(Nx1-1) + Nx2-1-k][j][Nz1-1+k].mask = TBorderPoint;
			MaskP[(Nx1-1) + k][j][Nz1-1+k].normal = new glTVector(1,0,1);
			MaskP[(Nx1-1) + Nx2-1-k][j][Nz1-1+k].normal = new glTVector(-1,0,1);
		}
		for(int k=0; k<Nz1-1; ++k)
		{
			MaskP[Nx1-1][j][k].mask = MaskP[Nx1+Nx2-2][j][k].mask = TBorderPoint;
			MaskP[Nx1-1][j][k].normal = new glTVector(1,0,1);
			MaskP[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,1);
		}
	}

	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
	{
		int p = i - (Nx1-1), q;
		if(p < Nz2)
			q = p;
		else if(p >= Nx2-Nz2)
			q = Nx2-1-p;
		else
			q = Nz2;

		for(int k=0; k<Nz1-1+q; ++k)
		{
			MaskP[i][Ny1-1][k].mask = MaskP[i][Ny1+Ny2-2][k].mask = TBorderPoint;
			MaskP[i][Ny1-1][k].normal = new glTVector(0,1,0);
			MaskP[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}
	}

	//// Ребра препятствия
	//for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
	//{
	//	MaskP[i][Ny1-1][Nz+Nz2-2].normal->set(0,1,-1);
	//	MaskP[i][Ny1+Ny2-2][Nz+Nz2-2].normal->set(0,-1,-1);
	//}
	//for(int i=(Nx1-1)+Nz-1; i<(Nx1-1)+Nx2-Nz+1; ++i)
	//{
	//	MaskP[i][Ny1-1][0].mask = TFictivePoint;
	//	MaskP[i][Ny1+Ny2-2][0].mask = TFictivePoint;
	//}

	//for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
	//{
	//	MaskP[Nx1-1][j][Nz+Nz2-2].normal->set(1,0,-1);
	//	MaskP[Nx1+Nx2-2][j][Nz+Nz2-2].normal->set(-1,0,-1);
	//}
	//for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
	//{
	//	MaskP[(Nx1-1)+Nz-1][j][0].mask = TFictivePoint;
	//	MaskP[(Nx1-1)+Nx2-Nz][j][0].mask = TFictivePoint;
	//}

	//for(int k=0; k<Nz; ++k)
	//{
	//	Mask[Nx1-1][Ny1-1][k].normal->set(1,1,0);
	//	Mask[Nx1+Nx2-2][Ny1-1][k].normal->set(-1,1,0);
	//	Mask[Nx1-1][Ny1+Ny2-2][k].normal->set(1,-1,0);
	//	Mask[Nx1+Nx2-2][Ny1+Ny2-2][k].normal->set(-1,-1,0);
	//}

	//// Угловые точки препятствия
	//Mask[Nx1-1][Ny1-1][Nz-1].normal->set(1,1,-1);
	//Mask[Nx1+Nx2-2][Ny1-1][Nz-1].normal->set(-1,1,-1);
	//Mask[Nx1-1][Ny1+Ny2-2][Nz-1].normal->set(1,-1,-1);
	//Mask[Nx1+Nx2-2][Ny1+Ny2-2][Nz-1].normal->set(-1,-1,-1);
}
//*************************************
void TTrapezeBottomWithPressureMaskGenerator::CreateMultyNodes()
{
	// Узлы для компоненты U скорости
	for(int i=0; i<NxU; ++i)
		XU[i] = X[i];

	for(int j=1; j<NyU-1; ++j)
		YU[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YU[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YU[NyU-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzU-1; ++k)
		ZU[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZU[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZU[NzU-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxU; ++i)
		HxU[i] = XU[i] - XU[i-1];
	HxU[0] = HxU[1];

	for(int j=1; j<NyU; ++j)
		HyU[j] = YU[j] - YU[j-1];
	HyU[0] = HyU[1];

	for(int k=1; k<NzU; ++k)
		HzU[k] = ZU[k] - ZU[k-1];
	HzU[0] = HzU[1];

	// Узлы для компоненты V скорости
	for(int i=1; i<NxV-1; ++i)
		XV[i] = ( X[i-1] + X[i] ) / 2.0;
	XV[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XV[NxV-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=0; j<NyV; ++j)
		YV[j] = Y[j];

	for(int k=1; k<NzV-1; ++k)
		ZV[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZV[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZV[NzV-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxV; ++i)
		HxV[i] = XV[i] - XV[i-1];
	HxV[0] = HxV[1];

	for(int j=1; j<NyV; ++j)
		HyV[j] = YV[j] - YV[j-1];
	HyV[0] = HyV[1];

	for(int k=1; k<NzV; ++k)
		HzV[k] = ZV[k] - ZV[k-1];
	HzV[0] = HzV[1];

	// Узлы для компоненты W скорости
	for(int i=1; i<NxW-1; ++i)
		XW[i] = ( X[i-1] + X[i] ) / 2.0;
	XW[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XW[NxW-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyW-1; ++j)
		YW[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YW[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YW[NyW-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=0; k<NzW; ++k)
		ZW[k] = Z[k];

	for(int i=1; i<NxW; ++i)
		HxW[i] = XW[i] - XW[i-1];
	HxW[0] = HxW[1];

	for(int j=1; j<NyW; ++j)
		HyW[j] = YW[j] - YW[j-1];
	HyW[0] = HyW[1];

	for(int k=1; k<NzW; ++k)
		HzW[k] = ZW[k] - ZW[k-1];
	HzW[0] = HzW[1];

	// Узлы для давления
	for(int i=1; i<NxP-1; ++i)
		XP[i] = ( X[i-1] + X[i] ) / 2.0;
	XP[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XP[NxP-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyP-1; ++j)
		YP[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YP[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YP[NyP-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzP-1; ++k)
		ZP[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZP[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZP[NzP-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxP; ++i)
		HxP[i] = XP[i] - XP[i-1];
	HxP[0] = HxP[1];

	for(int j=1; j<NyP; ++j)
		HyP[j] = YP[j] - YP[j-1];
	HyP[0] = HyP[1];

	for(int k=1; k<NzP; ++k)
		HzP[k] = ZP[k] - ZP[k-1];
	HzP[0] = HzP[1];
}
//*************************************
void TTrapezeBottomWithPressureMaskGenerator::CreateMultyGrids()
{
	// Заполнение маски для U //

	for(int i=0; i<NxU; ++i)
		for(int j=0; j<NyU; ++j)
			for(int k=0; k<NzU; ++k)
			{
				MaskU[i][j][k].mask = TActualPoint;
				MaskU[i][j][k].normal = 0;
			}

	// Границы области
	//for(int j=1; j<NyU-1; ++j)
	//	for(int k=1; k<NzU-1; ++k)
	//	{
	//		MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TDefinedBorderPoint;
	//		MaskU[0][j][k].normal = new glTVector(-1,0,0);
	//		MaskU[NxU-1][j][k].normal = new glTVector(1,0,0);
	//	}

	for(int j=1; j<NyU-1; ++j)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TNormalBorderPoint;
			MaskU[0][j][k].normal = new glTVector(-1,0,0);
			MaskU[NxU-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[i][0][k].mask = MaskU[i][NyU-1][k].mask = TPreDefinedBorderPoint;
			MaskU[i][0][k].normal =  new glTVector(0,-1,0);
			MaskU[i][NyU-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int j=1; j<NyU-1; ++j)
		{
			MaskU[i][j][0].mask = MaskU[i][j][NzU-1].mask = TPreDefinedBorderPoint;
			MaskU[i][j][0].normal = new glTVector(0,0,-1);
			MaskU[i][j][NzU-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxU; ++i)
	{
		MaskU[i][0][0].mask = TFictivePoint;
		MaskU[i][0][NzU-1].mask = TFictivePoint;
		MaskU[i][NyU-1][0].mask = TFictivePoint;
		MaskU[i][NyU-1][NzU-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyU; ++j)
	{
		MaskU[0][j][0].mask = TFictivePoint;
		MaskU[0][j][NzU-1].mask = TFictivePoint;
		MaskU[NxU-1][j][0].mask = TFictivePoint;
		MaskU[NxU-1][j][NzU-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzU; ++k)
	{
		MaskU[0][0][k].mask = TFictivePoint;
		MaskU[0][NyU-1][k].mask = TFictivePoint;
		MaskU[NxU-1][0][k].mask = TFictivePoint;
		MaskU[NxU-1][NyU-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			int p = i-(Nx1-1), q;
			if(p < Nz2-2)
				q = p+1;
			else if(p > Nx2-Nz2+1)
				q = Nx2-p;
			else
				q = Nz2;
			for(int k=0; k<Nz1+q; ++k)
				MaskU[i][j][k].mask = TFictivePoint;
		}

	// Грани препятствия
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
	{
		int p = i-(Nx1-1), q;
		if(p < Nz2-2)
			q = p+1;
		else if(p > Nx2-Nz2+1)
			q = Nx2-p;
		else
			q = Nz2;
		for(int k=0; k<Nz1+q; ++k)
		{
			MaskU[i][Ny1][k].mask = MaskU[i][Ny1+Ny2-2][k].mask = TPreDefinedBorderPoint;
			MaskU[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskU[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}
	}

	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	{
		for(int k=1; k<Nz1; ++k)
		{
			MaskU[Nx1-1][j][k].mask = MaskU[Nx1+Nx2-2][j][k].mask = TDefinedBorderPoint;
			MaskU[Nx1-1][j][k].normal = new glTVector(1,0,0);
			MaskU[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		for(int k=1; k<Nz2; ++k)
		{
			MaskU[(Nx1-1) + k-1][j][Nz1-1+k].mask = MaskU[(Nx1-1) + (Nx2+1)-k-1][j][Nz1-1+k].mask = TDefinedBorderPoint;
			MaskU[(Nx1-1) + k-1][j][Nz1-1+k].normal = new glTVector(1,0,1);
			MaskU[(Nx1-1) + (Nx2+1)-k-1][j][Nz1-1+k].normal = new glTVector(-1,0,1);
		}
	}

	// Заполнение маски для V //

	for(int i=0; i<NxV; ++i)
		for(int j=0; j<NyV; ++j)
			for(int k=0; k<NzV; ++k)
			{
				MaskV[i][j][k].mask = TActualPoint;
				MaskV[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyV-1; ++j)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[0][j][k].mask = MaskV[NxV-1][j][k].mask = TPreDefinedBorderPoint;
			MaskV[0][j][k].normal = new glTVector(-1,0,0);
			MaskV[NxV-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[i][0][k].mask = MaskV[i][NyV-1][k].mask = TDefinedBorderPoint;
			MaskV[i][0][k].normal =  new glTVector(0,-1,0);
			MaskV[i][NyV-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int j=1; j<NyV-1; ++j)
		{
			MaskV[i][j][0].mask = MaskV[i][j][NzV-1].mask = TPreDefinedBorderPoint;
			MaskV[i][j][0].normal = new glTVector(0,0,-1);
			MaskV[i][j][NzV-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxV; ++i)
	{
		MaskV[i][0][0].mask = TFictivePoint;
		MaskV[i][0][NzV-1].mask = TFictivePoint;
		MaskV[i][NyV-1][0].mask = TFictivePoint;
		MaskV[i][NyV-1][NzV-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyV; ++j)
	{
		MaskV[0][j][0].mask = TFictivePoint;
		MaskV[0][j][NzV-1].mask = TFictivePoint;
		MaskV[NxV-1][j][0].mask = TFictivePoint;
		MaskV[NxV-1][j][NzV-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzV; ++k)
	{
		MaskV[0][0][k].mask = TFictivePoint;
		MaskV[0][NyV-1][k].mask = TFictivePoint;
		MaskV[NxV-1][0][k].mask = TFictivePoint;
		MaskV[NxV-1][NyV-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		{
			int p = i-(Nx1-1), q;
			if(p < Nz2-1)
				q = p;
			else if(p > Nx2-Nz2+1)
				q = Nx2-p;
			else
				q = Nz2;
			for(int k=0; k<Nz1+q; ++k)
				MaskV[i][j][k].mask = TFictivePoint;
		}

	// Грани препятствия
	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
	{
		for(int k=1; k<Nz1; ++k)
		{
			MaskV[Nx1][j][k].mask = MaskV[Nx1+Nx2-2][j][k].mask = TPreDefinedBorderPoint;
			MaskV[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskV[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		for(int k=1; k<Nz2; ++k)
		{
			MaskV[(Nx1-1) + k][j][Nz1-1+k].mask = MaskV[(Nx1-1) + Nx2-k][j][Nz1-1+k].mask = TPreDefinedBorderPoint;
			MaskV[(Nx1-1) + k][j][Nz1-1+k].normal = new glTVector(1,0,1);
			MaskV[(Nx1-1) + Nx2-k][j][Nz1-1+k].normal = new glTVector(-1,0,1);
		}
	}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	{
		int p = i-(Nx1-1), q;
		if(p < Nz2-1)
			q = p;
		else if(p > Nx2-Nz2+1)
			q = Nx2-p;
		else
			q = Nz2;
		for(int k=0; k<Nz1+q; ++k)
		{
			MaskV[i][Ny1-1][k].mask = MaskV[i][Ny1+Ny2-2][k].mask = TDefinedBorderPoint;
			MaskV[i][Ny1-1][k].normal = new glTVector(0,1,0);
			MaskV[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}
	}

	// Заполнение маски для W //

	for(int i=0; i<NxW; ++i)
		for(int j=0; j<NyW; ++j)
			for(int k=0; k<NzW; ++k)
			{
				MaskW[i][j][k].mask = TActualPoint;
				MaskW[i][j][k].normal = 0;
			}
	// Границы области
	for(int j=1; j<NyW-1; ++j)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[0][j][k].mask = MaskW[NxW-1][j][k].mask = TPreDefinedBorderPoint; 
			MaskW[0][j][k].normal = new glTVector(-1,0,0);
			MaskW[NxW-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[i][0][k].mask = MaskW[i][NyW-1][k].mask = TPreDefinedBorderPoint;
			MaskW[i][0][k].normal =  new glTVector(0,-1,0);
			MaskW[i][NyW-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int j=1; j<NyW-1; ++j)
		{
			MaskW[i][j][0].mask = MaskW[i][j][NzW-1].mask = TDefinedBorderPoint;
			MaskW[i][j][0].normal = new glTVector(0,0,-1);
			MaskW[i][j][NzW-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxW; ++i)
	{
		MaskW[i][0][0].mask = TFictivePoint;
		MaskW[i][0][NzW-1].mask = TFictivePoint;
		MaskW[i][NyW-1][0].mask = TFictivePoint;
		MaskW[i][NyW-1][NzW-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyW; ++j)
	{
		MaskW[0][j][0].mask = TFictivePoint;
		MaskW[0][j][NzW-1].mask = TFictivePoint;
		MaskW[NxW-1][j][0].mask = TFictivePoint;
		MaskW[NxW-1][j][NzW-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzW; ++k)
	{
		MaskW[0][0][k].mask = TFictivePoint;
		MaskW[0][NyW-1][k].mask = TFictivePoint;
		MaskW[NxW-1][0][k].mask = TFictivePoint;
		MaskW[NxW-1][NyW-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			int p = i-(Nx1-1), q;
			if(p < Nz2-1)
				q = p;
			else if(p > Nx2-Nz2+1)
				q = Nx2-p;
			else
				q = Nz2-1;
			for(int k=0; k<Nz1+q; ++k)
				MaskW[i][j][k].mask = TFictivePoint;
		}

	// Грани препятствия
	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	{
		for(int k=1; k<Nz1; ++k)
		{
			MaskW[Nx1][j][k].mask = MaskW[Nx1+Nx2-2][j][k].mask = TDefinedBorderPoint;
			MaskW[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskW[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		for(int k=1; k<Nz2-1; ++k)
		{
			MaskW[(Nx1-1) + k][j][Nz1-1+k].mask = MaskW[(Nx1-1) + Nx2-k][j][Nz1-1+k].mask = TDefinedBorderPoint;
			MaskW[(Nx1-1) + k][j][Nz1-1+k].normal = new glTVector(1,0,1);
			MaskW[(Nx1-1) + Nx2-k][j][Nz1-1+k].normal = new glTVector(-1,0,1);
		}
	}	

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	{
		int p = i-(Nx1-1), q;
		if(p < Nz2-1)
			q = p;
		else if(p > Nx2-Nz2+1)
			q = Nx2-p;
		else
			q = Nz2-1;
		for(int k=0; k<Nz1+q; ++k)
		{
			MaskW[i][Ny1][k].mask = MaskW[i][Ny1+Ny2-2][k].mask = TPreDefinedBorderPoint;
			MaskW[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskW[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}
	}


	// Заполнение маски для P
	for(int i=0; i<NxP; ++i)
		for(int j=0; j<NyP; ++j)
			for(int k=0; k<NzP; ++k)
			{
				MaskP[i][j][k].mask = TActualPoint;
				MaskP[i][j][k].normal = 0;
			}

	// Границы области
	//for(int j=1; j<NyP-1; ++j)
	//	for(int k=1; k<NzP-1; ++k)
	//	{
	//		MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreNormalBorderPoint;
	//		MaskP[0][j][k].normal = new glTVector(-1,0,0);
	//		MaskP[NxP-1][j][k].normal = new glTVector(1,0,0);
	//	}

	for(int j=1; j<NyP-1; ++j)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreDefinedBorderPoint; 
			MaskP[0][j][k].normal = new glTVector(-1,0,0);
			MaskP[NxP-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[i][0][k].mask = MaskP[i][NyP-1][k].mask = TPreNormalBorderPoint;
			MaskP[i][0][k].normal =  new glTVector(0,-1,0);
			MaskP[i][NyP-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int j=1; j<NyP-1; ++j)
		{
			MaskP[i][j][0].mask = MaskP[i][j][NzP-1].mask = TPreNormalBorderPoint;
			MaskP[i][j][0].normal = new glTVector(0,0,-1);
			MaskP[i][j][NzP-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxP; ++i)
	{
		MaskP[i][0][0].mask = TFictivePoint;
		MaskP[i][0][NzP-1].mask = TFictivePoint;
		MaskP[i][NyP-1][0].mask = TFictivePoint;
		MaskP[i][NyP-1][NzP-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyP; ++j)
	{
		MaskP[0][j][0].mask = TFictivePoint;
		MaskP[0][j][NzP-1].mask = TFictivePoint;
		MaskP[NxP-1][j][0].mask = TFictivePoint;
		MaskP[NxP-1][j][NzP-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzP; ++k)
	{
		MaskP[0][0][k].mask = TFictivePoint;
		MaskP[0][NyP-1][k].mask = TFictivePoint;
		MaskP[NxP-1][0][k].mask = TFictivePoint;
		MaskP[NxP-1][NyP-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			int p = i - (Nx1-1), q;
			if(p < Nz2-1)
				q = p;
			else if(p > Nx2-Nz2+1)
				q = Nx2-p;
			else
				q = Nz2;
			for(int k=0; k<Nz1+q; ++k)
				MaskP[i][j][k].mask = TFictivePoint;
		}

	// Грани препятствия
	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	{
		for(int k=1; k<Nz1; ++k)
		{
			MaskP[Nx1][j][k].mask = MaskP[Nx1+Nx2-2][j][k].mask = TPreNormalBorderPoint;
			MaskP[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskP[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		for(int k=1; k<Nz2; ++k)
		{
			MaskP[(Nx1-1) + k][j][Nz1-1+k].mask = MaskP[(Nx1-1) + Nx2-k][j][Nz1-1+k].mask = TPreNormalBorderPoint;
			MaskP[(Nx1-1) + k][j][Nz1-1+k].normal = new glTVector(1,0,1);
			MaskP[(Nx1-1) + Nx2-k][j][Nz1-1+k].normal = new glTVector(-1,0,1);
		}
	}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	{
		int p = i - (Nx1-1), q;
		if(p < Nz2-1)
			q = p;
		else if(p > Nx2-Nz2+1)
			q = Nx2-p;
		else
			q = Nz2;
		for(int k=0; k<Nz1+q; ++k)
		{
			MaskP[i][Ny1][k].mask = MaskP[i][Ny1+Ny2-2][k].mask = TPreNormalBorderPoint;
			MaskP[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskP[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}
	}

	//// Ребра препятствия
	//for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	//{
	//	MaskP[i][Ny1][Nz+Nz2-2].normal->set(0,1,-1);
	//	MaskP[i][Ny1+Ny2-2][Nz+Nz2-2].normal->set(0,-1,-1);
	//}
	//for(int i=(Nx1-1)+Nz; i<(Nx1-1)+Nx2-Nz+1; ++i)
	//{
	//	MaskP[i][Ny1][0].mask = TFictivePoint;
	//	MaskP[i][Ny1+Ny2-2][0].mask = TFictivePoint;
	//}

	//for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	//{
	//	MaskP[Nx1][j][Nz+Nz2-2].normal->set(1,0,-1);
	//	MaskP[Nx1+Nx2-2][j][Nz+Nz2-2].normal->set(-1,0,-1);
	//}
	//for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	//{
	//	MaskP[(Nx1-1)+Nz][j][0].mask = TFictivePoint;
	//	MaskP[(Nx1-1)+Nx2-Nz][j][0].mask = TFictivePoint;
	//}

	//for(int k=0; k<Nz; ++k)
	//{
	//	MaskP[Nx1][Ny1][k].normal->set(1,1,0);
	//	MaskP[Nx1+Nx2-2][Ny1][k].normal->set(-1,1,0);
	//	MaskP[Nx1][Ny1+Ny2-2][k].normal->set(1,-1,0);
	//	MaskP[Nx1+Nx2-2][Ny1+Ny2-2][k].normal->set(-1,-1,0);
	//}

	// Угловые точки препятствия
	//MaskP[Nx1][Ny1][Nz-1].normal->set(1,1,-1);
	//MaskP[Nx1+Nx2-2][Ny1][Nz-1].normal->set(-1,1,-1);
	//MaskP[Nx1][Ny1+Ny2-2][Nz-1].normal->set(1,-1,-1);
	//MaskP[Nx1+Nx2-2][Ny1+Ny2-2][Nz-1].normal->set(-1,-1,-1);

	{

     /*

		FILE *f = fopen("u.txt", "w");
		for(int k=0; k<NzU; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyU-1; j>=0; --j)
			{
				for(int i=0; i<NxU; ++i)
					fprintf(f, "%4d ", MaskU[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);


		f = fopen("v.txt", "w");
		for(int k=0; k<NzV; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyV-1; j>=0; --j)
			{
				for(int i=0; i<NxV; ++i)
					fprintf(f, "%4d ", MaskV[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);

		f = fopen("w.txt", "w");
		for(int k=0; k<NzW; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyW-1; j>=0; --j)
			{
				for(int i=0; i<NxW; ++i)
					fprintf(f, "%4d ", MaskW[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);

		f = fopen("p.txt", "w");
		for(int k=0; k<NzP; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyP-1; j>=0; --j)
			{
				for(int i=0; i<NxP; ++i)
					fprintf(f, "%4d ", MaskP[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);
    */

	}

}
//*************************************
double TTrapezeBottomWithPressureMaskGenerator::getUCondition(int i, int j, int k, int n)
{

	//if(k==NzU-1)
	//{
	//	return 0.25;
	//}

	return 0;
}
//*************************************
double TTrapezeBottomWithPressureMaskGenerator::getVCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TTrapezeBottomWithPressureMaskGenerator::getWCondition(int i, int j, int k, int n)
{
	
	//if(  (k>=Nz1)&&(k<Nz1+Nz2-1)&&(i==(Nx1-1) + k - (Nz1-1)) ||
	//	 (k>=2)&&(k<Nz1)&&(i==Nx1)
	//	  )
	//{
	//	return
	//		-0.25;
	//}


		//for(int k=1; k<Nz1; ++k)
		//{
		//	MaskW[Nx1][j][k].mask = MaskW[Nx1+Nx2-2][j][k].mask = TDefinedBorderPoint;
		//	MaskW[Nx1][j][k].normal = new glTVector(1,0,0);
		//	MaskW[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		//}

		//for(int k=1; k<Nz2-1; ++k)
		//{
		//	MaskW[(Nx1-1) + k][j][Nz1-1+k].mask = MaskW[(Nx1-1) + Nx2-k][j][Nz1-1+k].mask = TDefinedBorderPoint;
		//	MaskW[(Nx1-1) + k][j][Nz1-1+k].normal = new glTVector(1,0,1);
		//	MaskW[(Nx1-1) + Nx2-k][j][Nz1-1+k].normal = new glTVector(-1,0,1);
		//}



	//if ( ((i==Nx1-1)||(i==Nx1))&&(j>Ny1-1)&&(j<Ny1+Ny2-1) )
	//{
	//	return
	//		-0.25;
	//}

	return 0;
}
//*************************************
double TTrapezeBottomWithPressureMaskGenerator::getPCondition(int i, int j, int k, int n)
{
	if(i==0)
		return PInner;
	else if(i==NxP-1)
		return POuter;
	else
		return 0;
}
//*************************************




//#################################################################




//*************************************
void TInvertedTrapezeBottomAlongXWithPressureMaskGenerator::CreateGeometryNodes()
{
	// Узлы, шаги сетки
	for(int i=0; i<Nx; ++i)
		X[i] = LengthX*i/(Nx-1);
	for(int i=1; i<Nx; ++i)
		Hx[i] = X[i]-X[i-1];
	Hx[0] = Hx[1];

	for(int j=0; j<Ny; ++j)
		Y[j] = LengthY*j/(Ny-1);
	for(int j=1; j<Ny; ++j)
		Hy[j] = Y[j]-Y[j-1];
	Hy[0] = Hy[1];

	for(int k=0; k<Nz; ++k)
		Z[k] = LengthZ*k/(Nz-1);
	for(int k=1; k<Nz; ++k)
		Hz[k] = Z[k]-Z[k-1];
	Hy[0] = Hy[1];
}
//*************************************
void TInvertedTrapezeBottomAlongXWithPressureMaskGenerator::CreateGeometryGrid()
{
	// Сначала все точки - расчетные
	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz; ++k)
			{
				Mask[i][j][k].mask = TActualPoint;
				Mask[i][j][k].normal = 0;
			}

	// Левая и правая границы
	for(int j=1; j<Ny-1; ++j)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[Nx-1][j][k].mask = Mask[0][j][k].mask = TBorderPoint;
			Mask[0][j][k].normal = new glTVector(-1,0,0);
			Mask[Nx-1][j][k].normal = new glTVector(1,0,0);
		}

	// Ближняя и дальняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[i][Ny-1][k].mask = Mask[i][0][k].mask = TBorderPoint;
			Mask[i][0][k].normal = new glTVector(0,-1,0);
			Mask[i][Ny-1][k].normal = new glTVector(0,1,0);
		}

	// Нижняя и верхняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int j=1; j<Ny-1; ++j)
		{
			Mask[i][j][Nz-1].mask = Mask[i][j][0].mask = TBorderPoint;
			Mask[i][j][0].normal = new glTVector(0,0,-1);
			Mask[i][j][Nz-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<Nx; ++i)
	{
		Mask[i][0][0].mask = TFictivePoint;
		Mask[i][0][Nz-1].mask = TFictivePoint;
		Mask[i][Ny-1][0].mask = TFictivePoint;
		Mask[i][Ny-1][Nz-1].mask = TFictivePoint;
	}

	for(int j=0; j<Ny; ++j)
	{
		Mask[0][j][0].mask = TFictivePoint;
		Mask[0][j][Nz-1].mask = TFictivePoint;
		Mask[Nx-1][j][0].mask = TFictivePoint;
		Mask[Nx-1][j][Nz-1].mask = TFictivePoint;
	}

	for(int k=0; k<Nz; ++k)
	{
		Mask[0][0][k].mask = TFictivePoint;
		Mask[0][Ny-1][k].mask = TFictivePoint;
		Mask[Nx-1][0][k].mask = TFictivePoint;
		Mask[Nx-1][Ny-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		{
			int p = i - (Nx1-1), q;
			if(p < Nz)
				q = Nz-1 - p;
			else if(p >= Nx2-Nz)
				q = p - (Nx2-Nz);
			else
				q = 0;

			for(int k=q; k<Nz; ++k)
				MaskP[i][j][k].mask = TFictivePoint;
		}

	// Грани препятствия
	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
	{
		for(int k=0; k<Nz-1; ++k)
		{
			MaskP[(Nx1-1) + Nz-1-k][j][k].mask = MaskP[(Nx1-1) + (Nx2-Nz)+k][j][k].mask = TBorderPoint;
			MaskP[(Nx1-1) + Nz-1-k][j][k].normal = new glTVector(1,0,1);
			MaskP[(Nx1-1) + (Nx2-Nz)+k][j][k].normal = new glTVector(-1,0,1);
		}
	}

	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
	{
		int p = i - (Nx1-1), q;
		if(p < Nz)
			q = Nz-1 - p;
		else if(p >= Nx2-Nz)
			q = p - (Nx2-Nz);
		else
			q = 0;

		for(int k=q; k<Nz-1; ++k)
		{
			MaskP[i][Ny1-1][k].mask = MaskP[i][Ny1+Ny2-2][k].mask = TBorderPoint;
			MaskP[i][Ny1-1][k].normal = new glTVector(0,1,0);
			MaskP[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}
	}

	//// Ребра препятствия
	//for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
	//{
	//	MaskP[i][Ny1-1][Nz+Nz2-2].normal->set(0,1,-1);
	//	MaskP[i][Ny1+Ny2-2][Nz+Nz2-2].normal->set(0,-1,-1);
	//}
	//for(int i=(Nx1-1)+Nz-1; i<(Nx1-1)+Nx2-Nz+1; ++i)
	//{
	//	MaskP[i][Ny1-1][0].mask = TFictivePoint;
	//	MaskP[i][Ny1+Ny2-2][0].mask = TFictivePoint;
	//}

	//for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
	//{
	//	MaskP[Nx1-1][j][Nz+Nz2-2].normal->set(1,0,-1);
	//	MaskP[Nx1+Nx2-2][j][Nz+Nz2-2].normal->set(-1,0,-1);
	//}
	//for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
	//{
	//	MaskP[(Nx1-1)+Nz-1][j][0].mask = TFictivePoint;
	//	MaskP[(Nx1-1)+Nx2-Nz][j][0].mask = TFictivePoint;
	//}

	//for(int k=0; k<Nz; ++k)
	//{
	//	Mask[Nx1-1][Ny1-1][k].normal->set(1,1,0);
	//	Mask[Nx1+Nx2-2][Ny1-1][k].normal->set(-1,1,0);
	//	Mask[Nx1-1][Ny1+Ny2-2][k].normal->set(1,-1,0);
	//	Mask[Nx1+Nx2-2][Ny1+Ny2-2][k].normal->set(-1,-1,0);
	//}

	//// Угловые точки препятствия
	//Mask[Nx1-1][Ny1-1][Nz-1].normal->set(1,1,-1);
	//Mask[Nx1+Nx2-2][Ny1-1][Nz-1].normal->set(-1,1,-1);
	//Mask[Nx1-1][Ny1+Ny2-2][Nz-1].normal->set(1,-1,-1);
	//Mask[Nx1+Nx2-2][Ny1+Ny2-2][Nz-1].normal->set(-1,-1,-1);
}
//*************************************
void TInvertedTrapezeBottomAlongXWithPressureMaskGenerator::CreateMultyNodes()
{
	// Узлы для компоненты U скорости
	for(int i=0; i<NxU; ++i)
		XU[i] = X[i];

	for(int j=1; j<NyU-1; ++j)
		YU[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YU[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YU[NyU-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzU-1; ++k)
		ZU[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZU[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZU[NzU-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxU; ++i)
		HxU[i] = XU[i] - XU[i-1];
	HxU[0] = HxU[1];

	for(int j=1; j<NyU; ++j)
		HyU[j] = YU[j] - YU[j-1];
	HyU[0] = HyU[1];

	for(int k=1; k<NzU; ++k)
		HzU[k] = ZU[k] - ZU[k-1];
	HzU[0] = HzU[1];

	// Узлы для компоненты V скорости
	for(int i=1; i<NxV-1; ++i)
		XV[i] = ( X[i-1] + X[i] ) / 2.0;
	XV[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XV[NxV-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=0; j<NyV; ++j)
		YV[j] = Y[j];

	for(int k=1; k<NzV-1; ++k)
		ZV[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZV[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZV[NzV-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxV; ++i)
		HxV[i] = XV[i] - XV[i-1];
	HxV[0] = HxV[1];

	for(int j=1; j<NyV; ++j)
		HyV[j] = YV[j] - YV[j-1];
	HyV[0] = HyV[1];

	for(int k=1; k<NzV; ++k)
		HzV[k] = ZV[k] - ZV[k-1];
	HzV[0] = HzV[1];

	// Узлы для компоненты W скорости
	for(int i=1; i<NxW-1; ++i)
		XW[i] = ( X[i-1] + X[i] ) / 2.0;
	XW[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XW[NxW-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyW-1; ++j)
		YW[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YW[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YW[NyW-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=0; k<NzW; ++k)
		ZW[k] = Z[k];

	for(int i=1; i<NxW; ++i)
		HxW[i] = XW[i] - XW[i-1];
	HxW[0] = HxW[1];

	for(int j=1; j<NyW; ++j)
		HyW[j] = YW[j] - YW[j-1];
	HyW[0] = HyW[1];

	for(int k=1; k<NzW; ++k)
		HzW[k] = ZW[k] - ZW[k-1];
	HzW[0] = HzW[1];

	// Узлы для давления
	for(int i=1; i<NxP-1; ++i)
		XP[i] = ( X[i-1] + X[i] ) / 2.0;
	XP[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XP[NxP-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyP-1; ++j)
		YP[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YP[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YP[NyP-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzP-1; ++k)
		ZP[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZP[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZP[NzP-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxP; ++i)
		HxP[i] = XP[i] - XP[i-1];
	HxP[0] = HxP[1];

	for(int j=1; j<NyP; ++j)
		HyP[j] = YP[j] - YP[j-1];
	HyP[0] = HyP[1];

	for(int k=1; k<NzP; ++k)
		HzP[k] = ZP[k] - ZP[k-1];
	HzP[0] = HzP[1];
}
//*************************************
void TInvertedTrapezeBottomAlongXWithPressureMaskGenerator::CreateMultyGrids()
{
	// Заполнение маски для U //

	for(int i=0; i<NxU; ++i)
		for(int j=0; j<NyU; ++j)
			for(int k=0; k<NzU; ++k)
			{
				MaskU[i][j][k].mask = TActualPoint;
				MaskU[i][j][k].normal = 0;
			}

	// Границы области
	//for(int j=1; j<NyU-1; ++j)
	//	for(int k=1; k<NzU-1; ++k)
	//	{
	//		MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TDefinedBorderPoint;
	//		MaskU[0][j][k].normal = new glTVector(-1,0,0);
	//		MaskU[NxU-1][j][k].normal = new glTVector(1,0,0);
	//	}

	for(int j=1; j<NyU-1; ++j)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TNormalBorderPoint;
			MaskU[0][j][k].normal = new glTVector(-1,0,0);
			MaskU[NxU-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[i][0][k].mask = MaskU[i][NyU-1][k].mask = TPreDefinedBorderPoint;
			MaskU[i][0][k].normal =  new glTVector(0,-1,0);
			MaskU[i][NyU-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int j=1; j<NyU-1; ++j)
		{
			MaskU[i][j][0].mask = MaskU[i][j][NzU-1].mask = TPreDefinedBorderPoint;
			MaskU[i][j][0].normal = new glTVector(0,0,-1);
			MaskU[i][j][NzU-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxU; ++i)
	{
		MaskU[i][0][0].mask = TFictivePoint;
		MaskU[i][0][NzU-1].mask = TFictivePoint;
		MaskU[i][NyU-1][0].mask = TFictivePoint;
		MaskU[i][NyU-1][NzU-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyU; ++j)
	{
		MaskU[0][j][0].mask = TFictivePoint;
		MaskU[0][j][NzU-1].mask = TFictivePoint;
		MaskU[NxU-1][j][0].mask = TFictivePoint;
		MaskU[NxU-1][j][NzU-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzU; ++k)
	{
		MaskU[0][0][k].mask = TFictivePoint;
		MaskU[0][NyU-1][k].mask = TFictivePoint;
		MaskU[NxU-1][0][k].mask = TFictivePoint;
		MaskU[NxU-1][NyU-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			int p = i-(Nx1-1), q;
			if(p < Nz-2)
				q = (Nz-1) - p;
			else if(p > Nx2-Nz+1)
				q = p - (Nx2-Nz);
			else
				q = 0;
			for(int k= q; k<Nz+1; ++k)
				MaskU[i][j][k].mask = TFictivePoint;
		}

	// Грани препятствия
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
	{
		int p = i-(Nx1-1), q;
		if(p < Nz-2)
			q = (Nz-1) - p;
		else if(p > Nx2-Nz+1)
			q = p - (Nx2-Nz);
		else
			q = 0;
		for(int k=q; k<Nz; ++k)
		{
			MaskU[i][Ny1][k].mask = MaskU[i][Ny1+Ny2-2][k].mask = TPreDefinedBorderPoint;
			MaskU[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskU[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}
	}

	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	{
		for(int k=1; k<Nz; ++k)
		{
			MaskU[(Nx1-1) + (Nz-1)-k][j][k].mask = MaskU[(Nx1-1) + (Nx2-Nz)+k][j][k].mask = TDefinedBorderPoint;
			MaskU[(Nx1-1) + (Nz-1)-k][j][k].normal = new glTVector(1,0,1);
			MaskU[(Nx1-1) + (Nx2-Nz)+k][j][k].normal = new glTVector(-1,0,1);
		}
	}

	// Заполнение маски для V //

	for(int i=0; i<NxV; ++i)
		for(int j=0; j<NyV; ++j)
			for(int k=0; k<NzV; ++k)
			{
				MaskV[i][j][k].mask = TActualPoint;
				MaskV[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyV-1; ++j)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[0][j][k].mask = MaskV[NxV-1][j][k].mask = TPreDefinedBorderPoint;
			MaskV[0][j][k].normal = new glTVector(-1,0,0);
			MaskV[NxV-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[i][0][k].mask = MaskV[i][NyV-1][k].mask = TDefinedBorderPoint;
			MaskV[i][0][k].normal =  new glTVector(0,-1,0);
			MaskV[i][NyV-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int j=1; j<NyV-1; ++j)
		{
			MaskV[i][j][0].mask = MaskV[i][j][NzV-1].mask = TPreDefinedBorderPoint;
			MaskV[i][j][0].normal = new glTVector(0,0,-1);
			MaskV[i][j][NzV-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxV; ++i)
	{
		MaskV[i][0][0].mask = TFictivePoint;
		MaskV[i][0][NzV-1].mask = TFictivePoint;
		MaskV[i][NyV-1][0].mask = TFictivePoint;
		MaskV[i][NyV-1][NzV-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyV; ++j)
	{
		MaskV[0][j][0].mask = TFictivePoint;
		MaskV[0][j][NzV-1].mask = TFictivePoint;
		MaskV[NxV-1][j][0].mask = TFictivePoint;
		MaskV[NxV-1][j][NzV-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzV; ++k)
	{
		MaskV[0][0][k].mask = TFictivePoint;
		MaskV[0][NyV-1][k].mask = TFictivePoint;
		MaskV[NxV-1][0][k].mask = TFictivePoint;
		MaskV[NxV-1][NyV-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		{
			int p = i-(Nx1-1), q;
			if(p < Nz-1)
				q = Nz - p;
			else if(p > Nx2-Nz+1)
				q = p - (Nx2-Nz);
			else
				q = 0;
			for(int k=q; k<Nz+1; ++k)
				MaskV[i][j][k].mask = TFictivePoint;
		}

	// Грани препятствия
	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
	{
		for(int k=1; k<Nz; ++k)
		{
			MaskV[(Nx1-1) + Nz-k][j][k].mask = MaskV[(Nx1-1) + (Nx2-Nz)+k][j][k].mask = TPreDefinedBorderPoint;
			MaskV[(Nx1-1) + Nz-k][j][k].normal = new glTVector(1,0,1);
			MaskV[(Nx1-1) + (Nx2-Nz)+k][j][k].normal = new glTVector(-1,0,1);
		}
	}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	{
		int p = i-(Nx1-1), q;
		if(p < Nz-1)
			q = Nz - p;
		else if(p > Nx2-Nz+1)
			q = p - (Nx2-Nz);
		else
			q = 0;
		for(int k=q; k<Nz; ++k)
		{
			MaskV[i][Ny1-1][k].mask = MaskV[i][Ny1+Ny2-2][k].mask = TDefinedBorderPoint;
			MaskV[i][Ny1-1][k].normal = new glTVector(0,1,0);
			MaskV[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}
	}

	// Заполнение маски для W //

	for(int i=0; i<NxW; ++i)
		for(int j=0; j<NyW; ++j)
			for(int k=0; k<NzW; ++k)
			{
				MaskW[i][j][k].mask = TActualPoint;
				MaskW[i][j][k].normal = 0;
			}
	// Границы области
	for(int j=1; j<NyW-1; ++j)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[0][j][k].mask = MaskW[NxW-1][j][k].mask = TPreDefinedBorderPoint; 
			MaskW[0][j][k].normal = new glTVector(-1,0,0);
			MaskW[NxW-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[i][0][k].mask = MaskW[i][NyW-1][k].mask = TPreDefinedBorderPoint;
			MaskW[i][0][k].normal =  new glTVector(0,-1,0);
			MaskW[i][NyW-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int j=1; j<NyW-1; ++j)
		{
			MaskW[i][j][0].mask = MaskW[i][j][NzW-1].mask = TDefinedBorderPoint;
			MaskW[i][j][0].normal = new glTVector(0,0,-1);
			MaskW[i][j][NzW-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxW; ++i)
	{
		MaskW[i][0][0].mask = TFictivePoint;
		MaskW[i][0][NzW-1].mask = TFictivePoint;
		MaskW[i][NyW-1][0].mask = TFictivePoint;
		MaskW[i][NyW-1][NzW-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyW; ++j)
	{
		MaskW[0][j][0].mask = TFictivePoint;
		MaskW[0][j][NzW-1].mask = TFictivePoint;
		MaskW[NxW-1][j][0].mask = TFictivePoint;
		MaskW[NxW-1][j][NzW-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzW; ++k)
	{
		MaskW[0][0][k].mask = TFictivePoint;
		MaskW[0][NyW-1][k].mask = TFictivePoint;
		MaskW[NxW-1][0][k].mask = TFictivePoint;
		MaskW[NxW-1][NyW-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			int p = i-(Nx1-1), q;
			if(p < Nz-1)
				q = (Nz-1) - p;
			else if(p > Nx2-Nz+1)
				q = p - (Nx2-Nz+1);
			else
				q = 0;
			for(int k=q; k<Nz; ++k)
				MaskW[i][j][k].mask = TFictivePoint;
		}

	// Грани препятствия
	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	{
		for(int k=1; k<Nz-1; ++k)
		{
			MaskW[(Nx1-1) + (Nz-1)-k][j][k].mask = MaskW[(Nx1-1) + (Nx2-Nz+1)+k][j][k].mask = TDefinedBorderPoint;
			MaskW[(Nx1-1) + (Nz-1)-k][j][k].normal = new glTVector(1,0,1);
			MaskW[(Nx1-1) + (Nx2-Nz+1)+k][j][k].normal = new glTVector(-1,0,1);
		}
	}	

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	{
		int p = i-(Nx1-1), q;
		if(p < Nz-1)
			q = (Nz-1) - p;
		else if(p > Nx2-Nz+1)
			q = p - (Nx2-Nz+1);
		else
			q = 0;
		for(int k=q; k<Nz-1; ++k)
		{
			MaskW[i][Ny1][k].mask = MaskW[i][Ny1+Ny2-2][k].mask = TPreDefinedBorderPoint;
			MaskW[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskW[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}
	}


	// Заполнение маски для P
	for(int i=0; i<NxP; ++i)
		for(int j=0; j<NyP; ++j)
			for(int k=0; k<NzP; ++k)
			{
				MaskP[i][j][k].mask = TActualPoint;
				MaskP[i][j][k].normal = 0;
			}

	// Границы области
	//for(int j=1; j<NyP-1; ++j)
	//	for(int k=1; k<NzP-1; ++k)
	//	{
	//		MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreNormalBorderPoint;
	//		MaskP[0][j][k].normal = new glTVector(-1,0,0);
	//		MaskP[NxP-1][j][k].normal = new glTVector(1,0,0);
	//	}

	for(int j=1; j<NyP-1; ++j)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreDefinedBorderPoint; 
			MaskP[0][j][k].normal = new glTVector(-1,0,0);
			MaskP[NxP-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[i][0][k].mask = MaskP[i][NyP-1][k].mask = TPreNormalBorderPoint;
			MaskP[i][0][k].normal =  new glTVector(0,-1,0);
			MaskP[i][NyP-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int j=1; j<NyP-1; ++j)
		{
			MaskP[i][j][0].mask = MaskP[i][j][NzP-1].mask = TPreNormalBorderPoint;
			MaskP[i][j][0].normal = new glTVector(0,0,-1);
			MaskP[i][j][NzP-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxP; ++i)
	{
		MaskP[i][0][0].mask = TFictivePoint;
		MaskP[i][0][NzP-1].mask = TFictivePoint;
		MaskP[i][NyP-1][0].mask = TFictivePoint;
		MaskP[i][NyP-1][NzP-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyP; ++j)
	{
		MaskP[0][j][0].mask = TFictivePoint;
		MaskP[0][j][NzP-1].mask = TFictivePoint;
		MaskP[NxP-1][j][0].mask = TFictivePoint;
		MaskP[NxP-1][j][NzP-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzP; ++k)
	{
		MaskP[0][0][k].mask = TFictivePoint;
		MaskP[0][NyP-1][k].mask = TFictivePoint;
		MaskP[NxP-1][0][k].mask = TFictivePoint;
		MaskP[NxP-1][NyP-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			int p = i - (Nx1-1), q;
			if(p < Nz-1)
				q = Nz - p;
			else if(p > Nx2-Nz+1)
				q = p - (Nx2-Nz);
			else
				q = 0;
			for(int k=q; k<Nz+1; ++k)
				MaskP[i][j][k].mask = TFictivePoint;
		}

	// Грани препятствия
	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	{
		for(int k=1; k<Nz; ++k)
		{
			MaskP[(Nx1-1) + Nz-k][j][k].mask = MaskP[(Nx1-1) + (Nx2-Nz)+k][j][k].mask = TPreNormalBorderPoint;
			MaskP[(Nx1-1) + Nz-k][j][k].normal = new glTVector(1,0,1);
			MaskP[(Nx1-1) + (Nx2-Nz)+k][j][k].normal = new glTVector(-1,0,1);
		}
	}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	{
		int p = i - (Nx1-1), q;
		if(p < Nz-1)
			q = Nz - p;
		else if(p > Nx2-Nz+1)
			q = p - (Nx2-Nz);
		else
			q = 0;
		for(int k=q; k<Nz; ++k)
		{
			MaskP[i][Ny1][k].mask = MaskP[i][Ny1+Ny2-2][k].mask = TPreNormalBorderPoint;
			MaskP[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskP[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}
	}

	//// Ребра препятствия
	//for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	//{
	//	MaskP[i][Ny1][Nz+Nz2-2].normal->set(0,1,-1);
	//	MaskP[i][Ny1+Ny2-2][Nz+Nz2-2].normal->set(0,-1,-1);
	//}
	//for(int i=(Nx1-1)+Nz; i<(Nx1-1)+Nx2-Nz+1; ++i)
	//{
	//	MaskP[i][Ny1][0].mask = TFictivePoint;
	//	MaskP[i][Ny1+Ny2-2][0].mask = TFictivePoint;
	//}

	//for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	//{
	//	MaskP[Nx1][j][Nz+Nz2-2].normal->set(1,0,-1);
	//	MaskP[Nx1+Nx2-2][j][Nz+Nz2-2].normal->set(-1,0,-1);
	//}
	//for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	//{
	//	MaskP[(Nx1-1)+Nz][j][0].mask = TFictivePoint;
	//	MaskP[(Nx1-1)+Nx2-Nz][j][0].mask = TFictivePoint;
	//}

	//for(int k=0; k<Nz; ++k)
	//{
	//	MaskP[Nx1][Ny1][k].normal->set(1,1,0);
	//	MaskP[Nx1+Nx2-2][Ny1][k].normal->set(-1,1,0);
	//	MaskP[Nx1][Ny1+Ny2-2][k].normal->set(1,-1,0);
	//	MaskP[Nx1+Nx2-2][Ny1+Ny2-2][k].normal->set(-1,-1,0);
	//}

	// Угловые точки препятствия
	//MaskP[Nx1][Ny1][Nz-1].normal->set(1,1,-1);
	//MaskP[Nx1+Nx2-2][Ny1][Nz-1].normal->set(-1,1,-1);
	//MaskP[Nx1][Ny1+Ny2-2][Nz-1].normal->set(1,-1,-1);
	//MaskP[Nx1+Nx2-2][Ny1+Ny2-2][Nz-1].normal->set(-1,-1,-1);

	{

     /*

		FILE *f = fopen("u.txt", "w");
		for(int k=0; k<NzU; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyU-1; j>=0; --j)
			{
				for(int i=0; i<NxU; ++i)
					fprintf(f, "%4d ", MaskU[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);


		f = fopen("v.txt", "w");
		for(int k=0; k<NzV; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyV-1; j>=0; --j)
			{
				for(int i=0; i<NxV; ++i)
					fprintf(f, "%4d ", MaskV[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);

		f = fopen("w.txt", "w");
		for(int k=0; k<NzW; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyW-1; j>=0; --j)
			{
				for(int i=0; i<NxW; ++i)
					fprintf(f, "%4d ", MaskW[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);

		f = fopen("p.txt", "w");
		for(int k=0; k<NzP; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyP-1; j>=0; --j)
			{
				for(int i=0; i<NxP; ++i)
					fprintf(f, "%4d ", MaskP[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);
    */

	}

}
//*************************************
double TInvertedTrapezeBottomAlongXWithPressureMaskGenerator::getUCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TInvertedTrapezeBottomAlongXWithPressureMaskGenerator::getVCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TInvertedTrapezeBottomAlongXWithPressureMaskGenerator::getWCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TInvertedTrapezeBottomAlongXWithPressureMaskGenerator::getPCondition(int i, int j, int k, int n)
{
	if(i==0)
		return PInner;
	else if(i==NxP-1)
		return POuter;
	else
		return 0;
}
//*************************************



//#################################################################




//*************************************
void TInvertedTrapezeBottomAlongYWithPressureMaskGenerator::CreateGeometryNodes()
{
	// Узлы, шаги сетки
	for(int i=0; i<Nx; ++i)
		X[i] = LengthX*i/(Nx-1);
	for(int i=1; i<Nx; ++i)
		Hx[i] = X[i]-X[i-1];
	Hx[0] = Hx[1];

	for(int j=0; j<Ny; ++j)
		Y[j] = LengthY*j/(Ny-1);
	for(int j=1; j<Ny; ++j)
		Hy[j] = Y[j]-Y[j-1];
	Hy[0] = Hy[1];

	for(int k=0; k<Nz; ++k)
		Z[k] = LengthZ*k/(Nz-1);
	for(int k=1; k<Nz; ++k)
		Hz[k] = Z[k]-Z[k-1];
	Hy[0] = Hy[1];
}
//*************************************
void TInvertedTrapezeBottomAlongYWithPressureMaskGenerator::CreateGeometryGrid()
{
	// Сначала все точки - расчетные
	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz; ++k)
			{
				Mask[i][j][k].mask = TActualPoint;
				Mask[i][j][k].normal = 0;
			}

	// Левая и правая границы
	for(int j=1; j<Ny-1; ++j)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[Nx-1][j][k].mask = Mask[0][j][k].mask = TBorderPoint;
			Mask[0][j][k].normal = new glTVector(-1,0,0);
			Mask[Nx-1][j][k].normal = new glTVector(1,0,0);
		}

	// Ближняя и дальняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[i][Ny-1][k].mask = Mask[i][0][k].mask = TBorderPoint;
			Mask[i][0][k].normal = new glTVector(0,-1,0);
			Mask[i][Ny-1][k].normal = new glTVector(0,1,0);
		}

	// Нижняя и верхняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int j=1; j<Ny-1; ++j)
		{
			Mask[i][j][Nz-1].mask = Mask[i][j][0].mask = TBorderPoint;
			Mask[i][j][0].normal = new glTVector(0,0,-1);
			Mask[i][j][Nz-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<Nx; ++i)
	{
		Mask[i][0][0].mask = TFictivePoint;
		Mask[i][0][Nz-1].mask = TFictivePoint;
		Mask[i][Ny-1][0].mask = TFictivePoint;
		Mask[i][Ny-1][Nz-1].mask = TFictivePoint;
	}

	for(int j=0; j<Ny; ++j)
	{
		Mask[0][j][0].mask = TFictivePoint;
		Mask[0][j][Nz-1].mask = TFictivePoint;
		Mask[Nx-1][j][0].mask = TFictivePoint;
		Mask[Nx-1][j][Nz-1].mask = TFictivePoint;
	}

	for(int k=0; k<Nz; ++k)
	{
		Mask[0][0][k].mask = TFictivePoint;
		Mask[0][Ny-1][k].mask = TFictivePoint;
		Mask[Nx-1][0][k].mask = TFictivePoint;
		Mask[Nx-1][Ny-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		{
			int p = j - (Ny1-1), q;
			if(p < Nz)
				q = Nz-1 - p;
			else if(p >= Ny2-Nz)
				q = p - (Ny2-Nz);
			else
				q = 0;

			for(int k=q; k<Nz; ++k)
				MaskP[i][j][k].mask = TFictivePoint;
		}

	// Грани препятствия
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
	{
		for(int k=0; k<Nz; ++k)
		{
			MaskP[i][(Ny1-1) + Nz-1-k][k].mask = MaskP[i][(Ny1-1) + (Ny2-Nz)+k][k].mask = TBorderPoint;
			MaskP[i][(Ny1-1) + Nz-1-k][k].normal = new glTVector(0,1,1);
			MaskP[i][(Ny1-1) + (Ny2-Nz)+k][k].normal = new glTVector(0,-1,1);
		}
	}

	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
	{
		int p = j - (Ny1-1), q;
		if(p < Nz)
			q = Nz-1 - p;
		else if(p >= Ny2-Nz)
			q = p - (Ny2-Nz);
		else
			q = 0;

		for(int k=q; k<Nz-1; ++k)
		{
			MaskP[Nx1-1][j][k].mask = MaskP[Nx1+Nx2-2][j][k].mask = TBorderPoint;
			MaskP[Nx1-1][j][k].normal = new glTVector(1,0,0);
			MaskP[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
	}

	//// Ребра препятствия
	//for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
	//{
	//	MaskP[i][Ny1-1][Nz1+Nz2-2].normal->set(0,1,-1);
	//	MaskP[i][Ny1+Ny2-2][Nz1+Nz2-2].normal->set(0,-1,-1);
	//}
	//for(int i=(Nx1-1)+Nz1-1; i<(Nx1-1)+Nx2-Nz1+1; ++i)
	//{
	//	MaskP[i][Ny1-1][0].mask = TFictivePoint;
	//	MaskP[i][Ny1+Ny2-2][0].mask = TFictivePoint;
	//}

	//for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
	//{
	//	MaskP[Nx1-1][j][Nz1+Nz2-2].normal->set(1,0,-1);
	//	MaskP[Nx1+Nx2-2][j][Nz1+Nz2-2].normal->set(-1,0,-1);
	//}
	//for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
	//{
	//	MaskP[(Nx1-1)+Nz1-1][j][0].mask = TFictivePoint;
	//	MaskP[(Nx1-1)+Nx2-Nz1][j][0].mask = TFictivePoint;
	//}

	//for(int k=0; k<Nz1; ++k)
	//{
	//	Mask[Nx1-1][Ny1-1][k].normal->set(1,1,0);
	//	Mask[Nx1+Nx2-2][Ny1-1][k].normal->set(-1,1,0);
	//	Mask[Nx1-1][Ny1+Ny2-2][k].normal->set(1,-1,0);
	//	Mask[Nx1+Nx2-2][Ny1+Ny2-2][k].normal->set(-1,-1,0);
	//}

	//// Угловые точки препятствия
	//Mask[Nx1-1][Ny1-1][Nz1-1].normal->set(1,1,-1);
	//Mask[Nx1+Nx2-2][Ny1-1][Nz1-1].normal->set(-1,1,-1);
	//Mask[Nx1-1][Ny1+Ny2-2][Nz1-1].normal->set(1,-1,-1);
	//Mask[Nx1+Nx2-2][Ny1+Ny2-2][Nz1-1].normal->set(-1,-1,-1);
}
//*************************************
void TInvertedTrapezeBottomAlongYWithPressureMaskGenerator::CreateMultyNodes()
{
	// Узлы для компоненты U скорости
	for(int i=0; i<NxU; ++i)
		XU[i] = X[i];

	for(int j=1; j<NyU-1; ++j)
		YU[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YU[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YU[NyU-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzU-1; ++k)
		ZU[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZU[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZU[NzU-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxU; ++i)
		HxU[i] = XU[i] - XU[i-1];
	HxU[0] = HxU[1];

	for(int j=1; j<NyU; ++j)
		HyU[j] = YU[j] - YU[j-1];
	HyU[0] = HyU[1];

	for(int k=1; k<NzU; ++k)
		HzU[k] = ZU[k] - ZU[k-1];
	HzU[0] = HzU[1];

	// Узлы для компоненты V скорости
	for(int i=1; i<NxV-1; ++i)
		XV[i] = ( X[i-1] + X[i] ) / 2.0;
	XV[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XV[NxV-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=0; j<NyV; ++j)
		YV[j] = Y[j];

	for(int k=1; k<NzV-1; ++k)
		ZV[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZV[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZV[NzV-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxV; ++i)
		HxV[i] = XV[i] - XV[i-1];
	HxV[0] = HxV[1];

	for(int j=1; j<NyV; ++j)
		HyV[j] = YV[j] - YV[j-1];
	HyV[0] = HyV[1];

	for(int k=1; k<NzV; ++k)
		HzV[k] = ZV[k] - ZV[k-1];
	HzV[0] = HzV[1];

	// Узлы для компоненты W скорости
	for(int i=1; i<NxW-1; ++i)
		XW[i] = ( X[i-1] + X[i] ) / 2.0;
	XW[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XW[NxW-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyW-1; ++j)
		YW[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YW[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YW[NyW-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=0; k<NzW; ++k)
		ZW[k] = Z[k];

	for(int i=1; i<NxW; ++i)
		HxW[i] = XW[i] - XW[i-1];
	HxW[0] = HxW[1];

	for(int j=1; j<NyW; ++j)
		HyW[j] = YW[j] - YW[j-1];
	HyW[0] = HyW[1];

	for(int k=1; k<NzW; ++k)
		HzW[k] = ZW[k] - ZW[k-1];
	HzW[0] = HzW[1];

	// Узлы для давления
	for(int i=1; i<NxP-1; ++i)
		XP[i] = ( X[i-1] + X[i] ) / 2.0;
	XP[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XP[NxP-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyP-1; ++j)
		YP[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YP[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YP[NyP-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzP-1; ++k)
		ZP[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZP[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZP[NzP-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxP; ++i)
		HxP[i] = XP[i] - XP[i-1];
	HxP[0] = HxP[1];

	for(int j=1; j<NyP; ++j)
		HyP[j] = YP[j] - YP[j-1];
	HyP[0] = HyP[1];

	for(int k=1; k<NzP; ++k)
		HzP[k] = ZP[k] - ZP[k-1];
	HzP[0] = HzP[1];
}
//*************************************
void TInvertedTrapezeBottomAlongYWithPressureMaskGenerator::CreateMultyGrids()
{
	// Заполнение маски для U //

	for(int i=0; i<NxU; ++i)
		for(int j=0; j<NyU; ++j)
			for(int k=0; k<NzU; ++k)
			{
				MaskU[i][j][k].mask = TActualPoint;
				MaskU[i][j][k].normal = 0;
			}

	// Границы области
	//for(int j=1; j<NyU-1; ++j)
	//	for(int k=1; k<NzU-1; ++k)
	//	{
	//		MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TDefinedBorderPoint;
	//		MaskU[0][j][k].normal = new glTVector(-1,0,0);
	//		MaskU[NxU-1][j][k].normal = new glTVector(1,0,0);
	//	}

	for(int j=1; j<NyU-1; ++j)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TNormalBorderPoint;
			MaskU[0][j][k].normal = new glTVector(-1,0,0);
			MaskU[NxU-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[i][0][k].mask = MaskU[i][NyU-1][k].mask = TPreDefinedBorderPoint;
			MaskU[i][0][k].normal =  new glTVector(0,-1,0);
			MaskU[i][NyU-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int j=1; j<NyU-1; ++j)
		{
			MaskU[i][j][0].mask = MaskU[i][j][NzU-1].mask = TPreDefinedBorderPoint;
			MaskU[i][j][0].normal = new glTVector(0,0,-1);
			MaskU[i][j][NzU-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxU; ++i)
	{
		MaskU[i][0][0].mask = TFictivePoint;
		MaskU[i][0][NzU-1].mask = TFictivePoint;
		MaskU[i][NyU-1][0].mask = TFictivePoint;
		MaskU[i][NyU-1][NzU-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyU; ++j)
	{
		MaskU[0][j][0].mask = TFictivePoint;
		MaskU[0][j][NzU-1].mask = TFictivePoint;
		MaskU[NxU-1][j][0].mask = TFictivePoint;
		MaskU[NxU-1][j][NzU-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzU; ++k)
	{
		MaskU[0][0][k].mask = TFictivePoint;
		MaskU[0][NyU-1][k].mask = TFictivePoint;
		MaskU[NxU-1][0][k].mask = TFictivePoint;
		MaskU[NxU-1][NyU-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		{
			int p = j-(Ny1-1), q;
			if(p < Nz-1)
				q = Nz - p;
			else if(p > Ny2-Nz+1)
				q = p - (Ny2-Nz);
			else
				q = 0;
			for(int k=q; k<Nz+1; ++k)
				MaskU[i][j][k].mask = TFictivePoint;
		}

	// Грани препятствия
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
	{
		for(int k=1; k<Nz; ++k)
		{
			MaskU[i][(Ny1-1) + Nz-k][k].mask = MaskU[i][(Ny1-1) + (Ny2-Nz)+k][k].mask = TPreDefinedBorderPoint;
			MaskU[i][(Ny1-1) + Nz-k][k].normal = new glTVector(0,1,1);
			MaskU[i][(Ny1-1) + (Ny2-Nz)+k][k].normal = new glTVector(0,-1,1);
		}
	}

	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	{
		int p = j-(Ny1-1), q;
		if(p < Nz-1)
			q = Nz - p;
		else if(p > Ny2-Nz+1)
			q = p - (Ny2-Nz);
		else
			q = 0;
		for(int k=q; k<Nz; ++k)
		{
			MaskU[Nx1-1][j][k].mask = MaskU[Nx1+Nx2-2][j][k].mask = TDefinedBorderPoint;
			MaskU[Nx1-1][j][k].normal = new glTVector(1,0,0);
			MaskU[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
	}

	// Заполнение маски для V //

	for(int i=0; i<NxV; ++i)
		for(int j=0; j<NyV; ++j)
			for(int k=0; k<NzV; ++k)
			{
				MaskV[i][j][k].mask = TActualPoint;
				MaskV[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyV-1; ++j)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[0][j][k].mask = MaskV[NxV-1][j][k].mask = TPreDefinedBorderPoint;
			MaskV[0][j][k].normal = new glTVector(-1,0,0);
			MaskV[NxV-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[i][0][k].mask = MaskV[i][NyV-1][k].mask = TDefinedBorderPoint;
			MaskV[i][0][k].normal =  new glTVector(0,-1,0);
			MaskV[i][NyV-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int j=1; j<NyV-1; ++j)
		{
			MaskV[i][j][0].mask = MaskV[i][j][NzV-1].mask = TPreDefinedBorderPoint;
			MaskV[i][j][0].normal = new glTVector(0,0,-1);
			MaskV[i][j][NzV-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxV; ++i)
	{
		MaskV[i][0][0].mask = TFictivePoint;
		MaskV[i][0][NzV-1].mask = TFictivePoint;
		MaskV[i][NyV-1][0].mask = TFictivePoint;
		MaskV[i][NyV-1][NzV-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyV; ++j)
	{
		MaskV[0][j][0].mask = TFictivePoint;
		MaskV[0][j][NzV-1].mask = TFictivePoint;
		MaskV[NxV-1][j][0].mask = TFictivePoint;
		MaskV[NxV-1][j][NzV-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzV; ++k)
	{
		MaskV[0][0][k].mask = TFictivePoint;
		MaskV[0][NyV-1][k].mask = TFictivePoint;
		MaskV[NxV-1][0][k].mask = TFictivePoint;
		MaskV[NxV-1][NyV-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		{
			int p = j-(Ny1-1), q;
			if(p < Nz-2)
				q = (Nz-1) - p;
			else if(p > Ny2-Nz+1)
				q = p - (Ny2-Nz);
			else
				q = 0;
			for(int k= q; k<Nz+1; ++k)
				MaskV[i][j][k].mask = TFictivePoint;
		}

	// Грани препятствия
	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
	{
		int p = j-(Ny1-1), q;
		if(p < Nz-2)
			q = (Nz-1) - p;
		else if(p > Ny2-Nz+1)
			q = p - (Ny2-Nz);
		else
			q = 0;
		for(int k=q; k<Nz; ++k)
		{
			MaskV[Nx1][j][k].mask = MaskV[Nx1+Nx2-2][j][k].mask = TPreDefinedBorderPoint;
			MaskV[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskV[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
	}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	{
		for(int k=1; k<Nz; ++k)
		{
			MaskV[i][(Ny1-1) + (Nz-1)-k][k].mask = MaskV[i][(Ny1-1) + (Ny2-Nz)+k][k].mask = TDefinedBorderPoint;
			MaskV[i][(Ny1-1) + (Nz-1)-k][k].normal = new glTVector(0,1,1);
			MaskV[i][(Ny1-1) + (Ny2-Nz)+k][k].normal = new glTVector(0,-1,1);
		}
	}



	// Заполнение маски для W //

	for(int i=0; i<NxW; ++i)
		for(int j=0; j<NyW; ++j)
			for(int k=0; k<NzW; ++k)
			{
				MaskW[i][j][k].mask = TActualPoint;
				MaskW[i][j][k].normal = 0;
			}
	// Границы области
	for(int j=1; j<NyW-1; ++j)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[0][j][k].mask = MaskW[NxW-1][j][k].mask = TPreDefinedBorderPoint; 
			MaskW[0][j][k].normal = new glTVector(-1,0,0);
			MaskW[NxW-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[i][0][k].mask = MaskW[i][NyW-1][k].mask = TPreDefinedBorderPoint;
			MaskW[i][0][k].normal =  new glTVector(0,-1,0);
			MaskW[i][NyW-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int j=1; j<NyW-1; ++j)
		{
			MaskW[i][j][0].mask = MaskW[i][j][NzW-1].mask = TDefinedBorderPoint;
			MaskW[i][j][0].normal = new glTVector(0,0,-1);
			MaskW[i][j][NzW-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxW; ++i)
	{
		MaskW[i][0][0].mask = TFictivePoint;
		MaskW[i][0][NzW-1].mask = TFictivePoint;
		MaskW[i][NyW-1][0].mask = TFictivePoint;
		MaskW[i][NyW-1][NzW-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyW; ++j)
	{
		MaskW[0][j][0].mask = TFictivePoint;
		MaskW[0][j][NzW-1].mask = TFictivePoint;
		MaskW[NxW-1][j][0].mask = TFictivePoint;
		MaskW[NxW-1][j][NzW-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzW; ++k)
	{
		MaskW[0][0][k].mask = TFictivePoint;
		MaskW[0][NyW-1][k].mask = TFictivePoint;
		MaskW[NxW-1][0][k].mask = TFictivePoint;
		MaskW[NxW-1][NyW-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		{
			int p = j-(Ny1-1), q;
			if(p < Nz-1)
				q = (Nz-1) - p;
			else if(p > Ny2-Nz+1)
				q = p - (Ny2-Nz+1);
			else
				q = 0;
			for(int k=q; k<Nz; ++k)
				MaskW[i][j][k].mask = TFictivePoint;
		}

	// Грани препятствия
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	{
		for(int k=1; k<Nz-1; ++k)
		{
			MaskW[i][(Ny1-1) + (Nz-1)-k][k].mask = MaskW[i][(Ny1-1) + (Ny2-Nz+1)+k][k].mask = TDefinedBorderPoint;
			MaskW[i][(Ny1-1) + (Nz-1)-k][k].normal = new glTVector(0,1,1);
			MaskW[i][(Ny1-1) + (Ny2-Nz+1)+k][k].normal = new glTVector(0,-1,1);
		}
	}	

	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	{
		int p = j-(Ny1-1), q;
		if(p < Nz-1)
			q = (Nz-1) - p;
		else if(p > Ny2-Nz+1)
			q = p - (Ny2-Nz+1);
		else
			q = 0;
		for(int k=q; k<Nz; ++k)
		{
			MaskW[Nx1][j][k].mask = MaskW[Nx1+Nx2-2][j][k].mask = TPreDefinedBorderPoint;
			MaskW[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskW[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
	}

	// Заполнение маски для P
	for(int i=0; i<NxP; ++i)
		for(int j=0; j<NyP; ++j)
			for(int k=0; k<NzP; ++k)
			{
				MaskP[i][j][k].mask = TActualPoint;
				MaskP[i][j][k].normal = 0;
			}

	// Границы области
	//for(int j=1; j<NyP-1; ++j)
	//	for(int k=1; k<NzP-1; ++k)
	//	{
	//		MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreNormalBorderPoint;
	//		MaskP[0][j][k].normal = new glTVector(-1,0,0);
	//		MaskP[NxP-1][j][k].normal = new glTVector(1,0,0);
	//	}

	for(int j=1; j<NyP-1; ++j)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreDefinedBorderPoint; 
			MaskP[0][j][k].normal = new glTVector(-1,0,0);
			MaskP[NxP-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[i][0][k].mask = MaskP[i][NyP-1][k].mask = TPreNormalBorderPoint;
			MaskP[i][0][k].normal =  new glTVector(0,-1,0);
			MaskP[i][NyP-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int j=1; j<NyP-1; ++j)
		{
			MaskP[i][j][0].mask = MaskP[i][j][NzP-1].mask = TPreNormalBorderPoint;
			MaskP[i][j][0].normal = new glTVector(0,0,-1);
			MaskP[i][j][NzP-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxP; ++i)
	{
		MaskP[i][0][0].mask = TFictivePoint;
		MaskP[i][0][NzP-1].mask = TFictivePoint;
		MaskP[i][NyP-1][0].mask = TFictivePoint;
		MaskP[i][NyP-1][NzP-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyP; ++j)
	{
		MaskP[0][j][0].mask = TFictivePoint;
		MaskP[0][j][NzP-1].mask = TFictivePoint;
		MaskP[NxP-1][j][0].mask = TFictivePoint;
		MaskP[NxP-1][j][NzP-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzP; ++k)
	{
		MaskP[0][0][k].mask = TFictivePoint;
		MaskP[0][NyP-1][k].mask = TFictivePoint;
		MaskP[NxP-1][0][k].mask = TFictivePoint;
		MaskP[NxP-1][NyP-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		{
			int p = j - (Ny1-1), q;
			if(p < Nz-1)
				q = Nz - p;
			else if(p > Ny2-Nz+1)
				q = p - (Ny2-Nz);
			else
				q = 0;
			for(int k=q; k<Nz; ++k)
				MaskP[i][j][k].mask = TFictivePoint;
		}

	// Грани препятствия
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	{
		for(int k=1; k<Nz; ++k)
		{
			MaskP[i][(Ny1-1) + Nz-k][k].mask = MaskP[i][(Ny1-1) + (Ny2-Nz)+k][k].mask = TPreNormalBorderPoint;
			MaskP[i][(Ny1-1) + Nz-k][k].normal = new glTVector(0,1,1);
			MaskP[i][(Ny1-1) + (Ny2-Nz)+k][k].normal = new glTVector(0,-1,1);
		}
	}

	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	{
		int p = j - (Ny1-1), q;
		if(p < Nz-1)
			q = Nz - p;
		else if(p > Ny2-Nz+1)
			q = p - (Ny2-Nz);
		else
			q = 0;
		for(int k=q; k<Nz; ++k)
		{
			MaskP[Nx1][j][k].mask = MaskP[Nx1+Nx2-2][j][k].mask = TPreNormalBorderPoint;
			MaskP[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskP[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
	}

	//// Ребра препятствия
	//for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	//{
	//	MaskP[i][Ny1][Nz1+Nz2-2].normal->set(0,1,-1);
	//	MaskP[i][Ny1+Ny2-2][Nz1+Nz2-2].normal->set(0,-1,-1);
	//}
	//for(int i=(Nx1-1)+Nz1; i<(Nx1-1)+Nx2-Nz1+1; ++i)
	//{
	//	MaskP[i][Ny1][0].mask = TFictivePoint;
	//	MaskP[i][Ny1+Ny2-2][0].mask = TFictivePoint;
	//}

	//for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	//{
	//	MaskP[Nx1][j][Nz1+Nz2-2].normal->set(1,0,-1);
	//	MaskP[Nx1+Nx2-2][j][Nz1+Nz2-2].normal->set(-1,0,-1);
	//}
	//for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	//{
	//	MaskP[(Nx1-1)+Nz1][j][0].mask = TFictivePoint;
	//	MaskP[(Nx1-1)+Nx2-Nz1][j][0].mask = TFictivePoint;
	//}

	//for(int k=0; k<Nz1; ++k)
	//{
	//	MaskP[Nx1][Ny1][k].normal->set(1,1,0);
	//	MaskP[Nx1+Nx2-2][Ny1][k].normal->set(-1,1,0);
	//	MaskP[Nx1][Ny1+Ny2-2][k].normal->set(1,-1,0);
	//	MaskP[Nx1+Nx2-2][Ny1+Ny2-2][k].normal->set(-1,-1,0);
	//}

	// Угловые точки препятствия
	//MaskP[Nx1][Ny1][Nz1-1].normal->set(1,1,-1);
	//MaskP[Nx1+Nx2-2][Ny1][Nz1-1].normal->set(-1,1,-1);
	//MaskP[Nx1][Ny1+Ny2-2][Nz1-1].normal->set(1,-1,-1);
	//MaskP[Nx1+Nx2-2][Ny1+Ny2-2][Nz1-1].normal->set(-1,-1,-1);

	{

      /*

		FILE *f = fopen("u.txt", "w");
		for(int k=0; k<NzU; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyU-1; j>=0; --j)
			{
				for(int i=0; i<NxU; ++i)
					fprintf(f, "%4d ", MaskU[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);


		f = fopen("v.txt", "w");
		for(int k=0; k<NzV; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyV-1; j>=0; --j)
			{
				for(int i=0; i<NxV; ++i)
					fprintf(f, "%4d ", MaskV[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);

		f = fopen("w.txt", "w");
		for(int k=0; k<NzW; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyW-1; j>=0; --j)
			{
				for(int i=0; i<NxW; ++i)
					fprintf(f, "%4d ", MaskW[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);

		f = fopen("p.txt", "w");
		for(int k=0; k<NzP; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyP-1; j>=0; --j)
			{
				for(int i=0; i<NxP; ++i)
					fprintf(f, "%4d ", MaskP[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);

       */
	}

}
//*************************************
double TInvertedTrapezeBottomAlongYWithPressureMaskGenerator::getUCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TInvertedTrapezeBottomAlongYWithPressureMaskGenerator::getVCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TInvertedTrapezeBottomAlongYWithPressureMaskGenerator::getWCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TInvertedTrapezeBottomAlongYWithPressureMaskGenerator::getPCondition(int i, int j, int k, int n)
{
	if(i==0)
		return PInner;
	else if(i==NxP-1)
		return POuter;
	else
		return 0;
}
//*************************************




//#################################################################


//*************************************
void TCubeInsideNewDivWithPressureMaskGenerator::CreateGeometryNodes()
{
	for(int i=0; i<Nx; ++i)
		X[i] = LengthX*i/(Nx-1);
	for(int i=1; i<Nx; ++i)
		Hx[i] = X[i]-X[i-1];
	Hx[0] = Hx[1];

	for(int j=0; j<Ny; ++j)
		Y[j] = LengthY*j/(Ny-1);
	for(int j=1; j<Ny; ++j)
		Hy[j] = Y[j]-Y[j-1];
	Hy[0] = Hy[1];

	for(int k=0; k<Nz; ++k)
		Z[k] = LengthZ*k/(Nz-1);
	for(int k=1; k<Nz; ++k)
		Hz[k] = Z[k]-Z[k-1];
	Hy[0] = Hy[1];
}
//*************************************
void TCubeInsideNewDivWithPressureMaskGenerator::CreateGeometryGrid()
{
	// Сначала все точки - расчетные
	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz; ++k)
			{
				Mask[i][j][k].mask = TActualPoint;
				Mask[i][j][k].normal = 0;
			}

	// Левая и правая границы
	for(int j=1; j<Ny-1; ++j)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[Nx-1][j][k].mask = Mask[0][j][k].mask = TBorderPoint;
			Mask[0][j][k].normal = new glTVector(-1,0,0);
			Mask[Nx-1][j][k].normal = new glTVector(1,0,0);
		}

	// Ближняя и дальняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[i][Ny-1][k].mask = Mask[i][0][k].mask = TBorderPoint;
			Mask[i][0][k].normal = new glTVector(0,-1,0);
			Mask[i][Ny-1][k].normal = new glTVector(0,1,0);
		}

	// Нижняя и верхняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int j=1; j<Ny-1; ++j)
		{
			Mask[i][j][Nz-1].mask = Mask[i][j][0].mask = TBorderPoint;
			Mask[i][j][0].normal = new glTVector(0,0,-1);
			Mask[i][j][Nz-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<Nx; ++i)
	{
		Mask[i][0][0].mask = TFictivePoint;
		Mask[i][0][Nz-1].mask = TFictivePoint;
		Mask[i][Ny-1][0].mask = TFictivePoint;
		Mask[i][Ny-1][Nz-1].mask = TFictivePoint;
	}

	for(int j=0; j<Ny; ++j)
	{
		Mask[0][j][0].mask = TFictivePoint;
		Mask[0][j][Nz-1].mask = TFictivePoint;
		Mask[Nx-1][j][0].mask = TFictivePoint;
		Mask[Nx-1][j][Nz-1].mask = TFictivePoint;
	}

	for(int k=0; k<Nz; ++k)
	{
		Mask[0][0][k].mask = TFictivePoint;
		Mask[0][Ny-1][k].mask = TFictivePoint;
		Mask[Nx-1][0][k].mask = TFictivePoint;
		Mask[Nx-1][Ny-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1; i<Nx1+Nx2-2; ++i)
		for(int j=Ny1; j<Ny1+Ny2-2; ++j)
			for(int k=Nz1; k<Nz1+Nz2-2; ++k)
				Mask[i][j][k].mask = TFictivePoint;
	
	// Грани препятствия
	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		for(int k=Nz1-1; k<Nz1+Nz2-1; ++k)
		{
			Mask[Nx1-1][j][k].mask = Mask[Nx1+Nx2-2][j][k].mask = TBorderPoint;
			Mask[Nx1-1][j][k].normal = new glTVector(1,0,0);
			Mask[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int k=Nz1-1; k<Nz1+Nz2-1; ++k)
		{
			Mask[i][Ny1-1][k].mask = Mask[i][Ny1+Ny2-2][k].mask = TBorderPoint;
			Mask[i][Ny1-1][k].normal = new glTVector(0,1,0);
			Mask[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		{
			Mask[i][j][Nz1-1].mask = Mask[i][j][Nz1+Nz2-2].mask = TBorderPoint;
			Mask[i][j][Nz1-1].normal = new glTVector(0,0,1);
			Mask[i][j][Nz1+Nz2-2].normal = new glTVector(0,0,-1);
		}

	// Ребра препятствия
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
	{
		Mask[i][Ny1-1][Nz1-1].normal->set(1,1,0);
		Mask[i][Ny1+Ny2-2][Nz1-1].normal->set(-1,1,0);
		Mask[i][Ny1-1][Nz1+Nz2-2].normal->set(1,-1,0);
		Mask[i][Ny1+Ny2-2][Nz1+Nz2-2].normal->set(-1,-1,0);
	}

	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
	{
		Mask[Nx1-1][j][Nz1-1].normal->set(1,0,1);
		Mask[Nx1+Nx2-2][j][Nz1-1].normal->set(-1,0,1);
		Mask[Nx1-1][j][Nz1+Nz2-2].normal->set(1,0,-1);
		Mask[Nx1+Nx2-2][j][Nz1+Nz2-2].normal->set(-1,0,-1);
	}

	for(int k=Nz1-1; k<Nz1+Nz2-1; ++k)
	{
		Mask[Nx1-1][Ny1-1][k].normal->set(0,1,1);
		Mask[Nx1+Nx2-2][Ny1-1][k].normal->set(0,-1,1);
		Mask[Nx1-1][Ny1+Ny2-2][k].normal->set(0,1,-1);
		Mask[Nx1+Nx2-2][Ny1+Ny2-2][k].normal->set(0,-1,-1);
	}

	// Угловые точки препятствия
	Mask[Nx1-1][Ny1-1][Nz1-1].normal->set(1,1,1);
	Mask[Nx1+Nx2-2][Ny1-1][Nz1-1].normal->set(-1,1,1);
	Mask[Nx1-1][Ny1-1][Nz1+Nz2-2].normal->set(1,1,-1);
	Mask[Nx1+Nx2-2][Ny1-1][Nz1+Nz2-2].normal->set(-1,1,-1);
	Mask[Nx1-1][Ny1+Ny2-2][Nz1-1].normal->set(1,-1,1);
	Mask[Nx1+Nx2-2][Ny1+Ny2-2][Nz1-1].normal->set(-1,-1,1);
	Mask[Nx1-1][Ny1+Ny2-2][Nz1+Nz2-2].normal->set(1,-1,-1);
	Mask[Nx1+Nx2-2][Ny1+Ny2-2][Nz1+Nz2-2].normal->set(-1,-1,-1);
}
//*************************************
void TCubeInsideNewDivWithPressureMaskGenerator::CreateMultyNodes()
{
	// Узлы для компоненты U скорости
	for(int i=0; i<NxU; ++i)
		XU[i] = X[i];

	for(int j=1; j<NyU-1; ++j)
		YU[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YU[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YU[NyU-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzU-1; ++k)
		ZU[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZU[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZU[NzU-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxU; ++i)
		HxU[i] = XU[i] - XU[i-1];
	HxU[0] = HxU[1];

	for(int j=1; j<NyU; ++j)
		HyU[j] = YU[j] - YU[j-1];
	HyU[0] = HyU[1];

	for(int k=1; k<NzU; ++k)
		HzU[k] = ZU[k] - ZU[k-1];
	HzU[0] = HzU[1];

	// Узлы для компоненты V скорости
	for(int i=1; i<NxV-1; ++i)
		XV[i] = ( X[i-1] + X[i] ) / 2.0;
	XV[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XV[NxV-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=0; j<NyV; ++j)
		YV[j] = Y[j];

	for(int k=1; k<NzV-1; ++k)
		ZV[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZV[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZV[NzV-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxV; ++i)
		HxV[i] = XV[i] - XV[i-1];
	HxV[0] = HxV[1];

	for(int j=1; j<NyV; ++j)
		HyV[j] = YV[j] - YV[j-1];
	HyV[0] = HyV[1];

	for(int k=1; k<NzV; ++k)
		HzV[k] = ZV[k] - ZV[k-1];
	HzV[0] = HzV[1];

	// Узлы для компоненты W скорости
	for(int i=1; i<NxW-1; ++i)
		XW[i] = ( X[i-1] + X[i] ) / 2.0;
	XW[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XW[NxW-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyW-1; ++j)
		YW[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YW[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YW[NyW-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=0; k<NzW; ++k)
		ZW[k] = Z[k];

	for(int i=1; i<NxW; ++i)
		HxW[i] = XW[i] - XW[i-1];
	HxW[0] = HxW[1];

	for(int j=1; j<NyW; ++j)
		HyW[j] = YW[j] - YW[j-1];
	HyW[0] = HyW[1];

	for(int k=1; k<NzW; ++k)
		HzW[k] = ZW[k] - ZW[k-1];
	HzW[0] = HzW[1];

	// Узлы для давления
	for(int i=1; i<NxP-1; ++i)
		XP[i] = ( X[i-1] + X[i] ) / 2.0;
	XP[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XP[NxP-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyP-1; ++j)
		YP[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YP[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YP[NyP-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzP-1; ++k)
		ZP[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZP[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZP[NzP-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxP; ++i)
		HxP[i] = XP[i] - XP[i-1];
	HxP[0] = HxP[1];

	for(int j=1; j<NyP; ++j)
		HyP[j] = YP[j] - YP[j-1];
	HyP[0] = HyP[1];

	for(int k=1; k<NzP; ++k)
		HzP[k] = ZP[k] - ZP[k-1];
	HzP[0] = HzP[1];
}
//*************************************
void TCubeInsideNewDivWithPressureMaskGenerator::CreateMultyGrids()
{
	// Заполнение маски для U //

	for(int i=0; i<NxU; ++i)
		for(int j=0; j<NyU; ++j)
			for(int k=0; k<NzU; ++k)
			{
				MaskU[i][j][k].mask = TActualPoint;
				MaskU[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyU-1; ++j)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TNormalBorderPoint;
			//MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TEquationBorderPoint;
			MaskU[0][j][k].normal = new glTVector(-1,0,0);
			MaskU[NxU-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[i][0][k].mask = MaskU[i][NyU-1][k].mask = TPreDefinedBorderPoint;
			MaskU[i][0][k].normal =  new glTVector(0,-1,0);
			MaskU[i][NyU-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int j=1; j<NyU-1; ++j)
		{
			MaskU[i][j][0].mask = MaskU[i][j][NzU-1].mask = TPreDefinedBorderPoint;
			MaskU[i][j][0].normal = new glTVector(0,0,-1);
			MaskU[i][j][NzU-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxU; ++i)
	{
		MaskU[i][0][0].mask = TFictivePoint;
		MaskU[i][0][NzU-1].mask = TFictivePoint;
		MaskU[i][NyU-1][0].mask = TFictivePoint;
		MaskU[i][NyU-1][NzU-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyU; ++j)
	{
		MaskU[0][j][0].mask = TFictivePoint;
		MaskU[0][j][NzU-1].mask = TFictivePoint;
		MaskU[NxU-1][j][0].mask = TFictivePoint;
		MaskU[NxU-1][j][NzU-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzU; ++k)
	{
		MaskU[0][0][k].mask = TFictivePoint;
		MaskU[0][NyU-1][k].mask = TFictivePoint;
		MaskU[NxU-1][0][k].mask = TFictivePoint;
		MaskU[NxU-1][NyU-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1; i<Nx1+Nx2-2; ++i)
		for(int j=Ny1+1; j<Ny1+Ny2-2; ++j)
			for(int k=Nz1+1; k<Nz1+Nz2-2; ++k)
				MaskU[i][j][k].mask = TFictivePoint;

	// Грани препятствия
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskU[i][Ny1][k].mask = MaskU[i][Ny1+Ny2-2][k].mask = TPreDefinedBorderPoint;
			MaskU[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskU[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			MaskU[i][j][Nz1].mask = MaskU[i][j][Nz1+Nz2-2].mask = TPreDefinedBorderPoint;
			MaskU[i][j][Nz1].normal = new glTVector(0,0,1);
			MaskU[i][j][Nz1+Nz2-2].normal = new glTVector(0,0,-1);
		}

	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskU[Nx1-1][j][k].mask = MaskU[Nx1+Nx2-2][j][k].mask = TDefinedBorderPoint;
			MaskU[Nx1-1][j][k].normal = new glTVector(1,0,0);
			MaskU[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}

	// Заполнение маски для V //
	for(int i=0; i<NxV; ++i)
		for(int j=0; j<NyV; ++j)
			for(int k=0; k<NzV; ++k)
			{
				MaskV[i][j][k].mask = TActualPoint;
				MaskV[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyV-1; ++j)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[0][j][k].mask = MaskV[NxV-1][j][k].mask = TPreDefinedBorderPoint;
			MaskV[0][j][k].normal = new glTVector(-1,0,0);
			MaskV[NxV-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[i][0][k].mask = MaskV[i][NyV-1][k].mask = TDefinedBorderPoint;
			MaskV[i][0][k].normal =  new glTVector(0,-1,0);
			MaskV[i][NyV-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int j=1; j<NyV-1; ++j)
		{
			MaskV[i][j][0].mask = MaskV[i][j][NzV-1].mask = TPreDefinedBorderPoint;
			MaskV[i][j][0].normal = new glTVector(0,0,-1);
			MaskV[i][j][NzV-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxV; ++i)
	{
		MaskV[i][0][0].mask = TFictivePoint;
		MaskV[i][0][NzV-1].mask = TFictivePoint;
		MaskV[i][NyV-1][0].mask = TFictivePoint;
		MaskV[i][NyV-1][NzV-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyV; ++j)
	{
		MaskV[0][j][0].mask = TFictivePoint;
		MaskV[0][j][NzV-1].mask = TFictivePoint;
		MaskV[NxV-1][j][0].mask = TFictivePoint;
		MaskV[NxV-1][j][NzV-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzV; ++k)
	{
		MaskV[0][0][k].mask = TFictivePoint;
		MaskV[0][NyV-1][k].mask = TFictivePoint;
		MaskV[NxV-1][0][k].mask = TFictivePoint;
		MaskV[NxV-1][NyV-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1+1; i<Nx1+Nx2-2; ++i)
		for(int j=Ny1; j<Ny1+Ny2-2; ++j)
			for(int k=Nz1+1; k<Nz1+Nz2-2; ++k)
				MaskV[i][j][k].mask = TFictivePoint;
	
	// Грани препятствия
	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskV[Nx1][j][k].mask = MaskV[Nx1+Nx2-2][j][k].mask = TPreDefinedBorderPoint;
			MaskV[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskV[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		{
			MaskV[i][j][Nz1].mask = MaskV[i][j][Nz1+Nz2-2].mask = TPreDefinedBorderPoint;
			MaskV[i][j][Nz1].normal = new glTVector(0,0,1);
			MaskV[i][j][Nz1+Nz2-2].normal = new glTVector(0,0,-1);
		}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskV[i][Ny1-1][k].mask = MaskV[i][Ny1+Ny2-2][k].mask = TDefinedBorderPoint;
			MaskV[i][Ny1-1][k].normal = new glTVector(0,1,0);
			MaskV[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	// Заполнение маски для W //

	for(int i=0; i<NxW; ++i)
		for(int j=0; j<NyW; ++j)
			for(int k=0; k<NzW; ++k)
			{
				MaskW[i][j][k].mask = TActualPoint;
				MaskW[i][j][k].normal = 0;
			}
	// Границы области
	for(int j=1; j<NyW-1; ++j)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[0][j][k].mask = MaskW[NxW-1][j][k].mask = TPreDefinedBorderPoint; 
			MaskW[0][j][k].normal = new glTVector(-1,0,0);
			MaskW[NxW-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[i][0][k].mask = MaskW[i][NyW-1][k].mask = TPreDefinedBorderPoint;
			MaskW[i][0][k].normal =  new glTVector(0,-1,0);
			MaskW[i][NyW-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int j=1; j<NyW-1; ++j)
		{
			MaskW[i][j][0].mask = MaskW[i][j][NzW-1].mask = TDefinedBorderPoint;
			MaskW[i][j][0].normal = new glTVector(0,0,-1);
			MaskW[i][j][NzW-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxW; ++i)
	{
		MaskW[i][0][0].mask = TFictivePoint;
		MaskW[i][0][NzW-1].mask = TFictivePoint;
		MaskW[i][NyW-1][0].mask = TFictivePoint;
		MaskW[i][NyW-1][NzW-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyW; ++j)
	{
		MaskW[0][j][0].mask = TFictivePoint;
		MaskW[0][j][NzW-1].mask = TFictivePoint;
		MaskW[NxW-1][j][0].mask = TFictivePoint;
		MaskW[NxW-1][j][NzW-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzW; ++k)
	{
		MaskW[0][0][k].mask = TFictivePoint;
		MaskW[0][NyW-1][k].mask = TFictivePoint;
		MaskW[NxW-1][0][k].mask = TFictivePoint;
		MaskW[NxW-1][NyW-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1+1; i<Nx1+Nx2-2; ++i)
		for(int j=Ny1+1; j<Ny1+Ny2-2; ++j)
			for(int k=Nz1; k<Nz1+Nz2-2; ++k)
				MaskW[i][j][k].mask = TFictivePoint;

	// Грани препятствия
	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		for(int k=Nz1-1; k<Nz1+Nz2-1; ++k)
		{
			MaskW[Nx1][j][k].mask = MaskW[Nx1+Nx2-2][j][k].mask = TPreDefinedBorderPoint;
			MaskW[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskW[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int k=Nz1-1; k<Nz1+Nz2-1; ++k)
		{
			MaskW[i][Ny1][k].mask = MaskW[i][Ny1+Ny2-2][k].mask = TPreDefinedBorderPoint;
			MaskW[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskW[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			MaskW[i][j][Nz1-1].mask = MaskW[i][j][Nz1+Nz2-2].mask = TDefinedBorderPoint;
			MaskW[i][j][Nz1-1].normal = new glTVector(0,0,1);
			MaskW[i][j][Nz1+Nz2-2].normal = new glTVector(0,0,-1);
		}

	// Заполнение маски для P
	for(int i=0; i<NxP; ++i)
		for(int j=0; j<NyP; ++j)
			for(int k=0; k<NzP; ++k)
			{
				MaskP[i][j][k].mask = TActualPoint;
				MaskP[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyP-1; ++j)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreDefinedBorderPoint;
			MaskP[0][j][k].normal = new glTVector(-1,0,0);
			MaskP[NxP-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[i][0][k].mask = MaskP[i][NyP-1][k].mask = TFictivePoint;
			MaskP[i][0][k].normal =  new glTVector(0,-1,0);
			MaskP[i][NyP-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int j=1; j<NyP-1; ++j)
		{
			MaskP[i][j][0].mask = MaskP[i][j][NzP-1].mask = TFictivePoint;
			MaskP[i][j][0].normal = new glTVector(0,0,-1);
			MaskP[i][j][NzP-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxP; ++i)
	{
		MaskP[i][0][0].mask = TFictivePoint;
		MaskP[i][0][NzP-1].mask = TFictivePoint;
		MaskP[i][NyP-1][0].mask = TFictivePoint;
		MaskP[i][NyP-1][NzP-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyP; ++j)
	{
		MaskP[0][j][0].mask = TFictivePoint;
		MaskP[0][j][NzP-1].mask = TFictivePoint;
		MaskP[NxP-1][j][0].mask = TFictivePoint;
		MaskP[NxP-1][j][NzP-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzP; ++k)
	{
		MaskP[0][0][k].mask = TFictivePoint;
		MaskP[0][NyP-1][k].mask = TFictivePoint;
		MaskP[NxP-1][0][k].mask = TFictivePoint;
		MaskP[NxP-1][NyP-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия (куб) - нерасчетные точки
	for(int i=Nx1+1; i<Nx1+Nx2-2; ++i)
		for(int j=Ny1+1; j<Ny1+Ny2-2; ++j)
			for(int k=Nz1+1; k<Nz1+Nz2-2; ++k)
				MaskP[i][j][k].mask = TFictivePoint;
	
	// Грани препятствия
	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskP[Nx1][j][k].mask = MaskP[Nx1+Nx2-2][j][k].mask = TFictivePoint;
			MaskP[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskP[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
		
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int k=Nz1; k<Nz1+Nz2-1; ++k)
		{
			MaskP[i][Ny1][k].mask = MaskP[i][Ny1+Ny2-2][k].mask = TFictivePoint;
			MaskP[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskP[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			MaskP[i][j][Nz1].mask = MaskP[i][j][Nz1+Nz2-2].mask = TFictivePoint;
			MaskP[i][j][Nz1].normal = new glTVector(0,0,1);
			MaskP[i][j][Nz1+Nz2-2].normal = new glTVector(0,0,-1);
		}

	// Ребра препятствия
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	{
		MaskP[i][Ny1][Nz1].normal->set(1,1,0);
		MaskP[i][Ny1+Ny2-2][Nz1].normal->set(-1,1,0);
		MaskP[i][Ny1][Nz1+Nz2-2].normal->set(1,-1,0);
		MaskP[i][Ny1+Ny2-2][Nz1+Nz2-2].normal->set(-1,-1,0);
	}

	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	{
		MaskP[Nx1][j][Nz1].normal->set(1,0,1);
		MaskP[Nx1+Nx2-2][j][Nz1].normal->set(-1,0,1);
		MaskP[Nx1][j][Nz1+Nz2-2].normal->set(1,0,-1);
		MaskP[Nx1+Nx2-2][j][Nz1+Nz2-2].normal->set(-1,0,-1);
	}

	for(int k=Nz1; k<Nz1+Nz2-1; ++k)
	{
		MaskP[Nx1][Ny1][k].normal->set(0,1,1);
		MaskP[Nx1+Nx2-2][Ny1][k].normal->set(0,-1,1);
		MaskP[Nx1][Ny1+Ny2-2][k].normal->set(0,1,-1);
		MaskP[Nx1+Nx2-2][Ny1+Ny2-2][k].normal->set(0,-1,-1);
	}

	// Угловые точки препятствия
	MaskP[Nx1][Ny1][Nz1].normal->set(1,1,1);
	MaskP[Nx1+Nx2-2][Ny1][Nz1].normal->set(-1,1,1);
	MaskP[Nx1][Ny1][Nz1+Nz2-2].normal->set(1,1,-1);
	MaskP[Nx1+Nx2-2][Ny1][Nz1+Nz2-2].normal->set(-1,1,-1);
	MaskP[Nx1][Ny1+Ny2-2][Nz1].normal->set(1,-1,1);
	MaskP[Nx1+Nx2-2][Ny1+Ny2-2][Nz1].normal->set(-1,-1,1);
	MaskP[Nx1][Ny1+Ny2-2][Nz1+Nz2-2].normal->set(1,-1,-1);
	MaskP[Nx1+Nx2-2][Ny1+Ny2-2][Nz1+Nz2-2].normal->set(-1,-1,-1);


	//Mask modification (Begin)
	for (int funcNum=1;funcNum<5;funcNum++)
	{
		glTMaskWithNormals*** mask = NULL;
	    int Nx = -1;
	    int Ny = -1;
	    int Nz = -1;

		if (funcNum==1) 
		{
			mask = MaskU;
			Nx = NxU; Ny = NyU; Nz = NzU;
		} 	 
		else if (funcNum==2)
        {
			mask = MaskV;
			Nx = NxV; Ny = NyV; Nz = NzV;
		} 	
		else if (funcNum==3)
        {
			mask = MaskW;
			Nx = NxW; Ny = NyW; Nz = NzW;
		} 	
		else if (funcNum==4)
        {
			mask = MaskP;
			Nx = NxP; Ny = NyP; Nz = NzP;
		} 	
		else {throw "Error in Mask generator: Unknown function number.";}

	    for(int i=0; i<Nx; ++i)
		{
		  for(int j=0; j<Ny; ++j)
		  {
			for(int k=0; k<Nz; ++k)
			{
				if((mask[i][j][k]).mask == TPreDefinedBorderPoint)
				{
					int actualValuesCount = 0;
					for (int dir=1;dir<4;dir++)
					{
						for (int sign=0;sign<2;sign++)
						{
							int shift = -1 + sign*2;
							int i__ = i;
							int j__ = j;
							int k__ = k;
							if (dir==1) {i__=i+shift;}
							else if (dir==2) {j__=j+shift;}
							else if (dir==3) {k__=k+shift;}
							else {throw "Error in Mask generator: Unknown direction.";} 
							if ( (i__>=0) && (i__<=Nx-1) && (j__>=0) && (j__<=Ny-1) && (k__>=0) && (k__<=Nz-1) )
							{
						       if ( (mask[i__][j__][k__]).mask==TActualPoint ) 
								{ 
									actualValuesCount = actualValuesCount + 1;
								}
							}
						}//by -1 or 1
					}//by direction

					if (actualValuesCount<=0)
					{
						throw "Error in Mask generator: Actual points number around PreDefined is zero."; 
					}
					else if (actualValuesCount==1)
					{
						//It is OK!
					}
					else if (actualValuesCount==2)
					{
						(mask[i][j][k]).mask = TDefinedBorderPoint;
					}
					else
					{
                       throw "Error in Mask generator: Actual points number around PreDefined more then two.";
					}
				}// if PreDefined
			}// by i
		  }//by j
		}//by k
	}//by function number
   //Mask modification (End)
}
//*************************************
double TCubeInsideNewDivWithPressureMaskGenerator::getUCondition(int i, int j, int k, int n)
{	
	return 0;	
}
//*************************************
double TCubeInsideNewDivWithPressureMaskGenerator::getVCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TCubeInsideNewDivWithPressureMaskGenerator::getWCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TCubeInsideNewDivWithPressureMaskGenerator::getPCondition(int i, int j, int k, int n)
{
	if(i==0)
		return PInner;
	else if(i==NxP-1)
		return POuter;
	else
		return 0;
}
//*************************************


//#################################################################



//*************************************
void TDoubleTrapezeBottomWithPressureMaskGenerator::CreateGeometryNodes()
{
	// Узлы, шаги сетки
	for(int i=0; i<Nx; ++i)
		X[i] = LengthX*i/(Nx-1);
	for(int i=1; i<Nx; ++i)
		Hx[i] = X[i]-X[i-1];
	Hx[0] = Hx[1];

	for(int j=0; j<Ny; ++j)
		Y[j] = LengthY*j/(Ny-1);
	for(int j=1; j<Ny; ++j)
		Hy[j] = Y[j]-Y[j-1];
	Hy[0] = Hy[1];

	for(int k=0; k<Nz; ++k)
		Z[k] = LengthZ*k/(Nz-1);
	for(int k=1; k<Nz; ++k)
		Hz[k] = Z[k]-Z[k-1];
	Hy[0] = Hy[1];
}
//*************************************
void TDoubleTrapezeBottomWithPressureMaskGenerator::CreateGeometryGrid()
{
	// Сначала все точки - расчетные
	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz; ++k)
			{
				Mask[i][j][k].mask = TActualPoint;
				Mask[i][j][k].normal = 0;
			}

	// Левая и правая границы
	for(int j=1; j<Ny-1; ++j)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[Nx-1][j][k].mask = Mask[0][j][k].mask = TBorderPoint;
			Mask[0][j][k].normal = new glTVector(-1,0,0);
			Mask[Nx-1][j][k].normal = new glTVector(1,0,0);
		}

	// Ближняя и дальняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[i][Ny-1][k].mask = Mask[i][0][k].mask = TBorderPoint;
			Mask[i][0][k].normal = new glTVector(0,-1,0);
			Mask[i][Ny-1][k].normal = new glTVector(0,1,0);
		}

	// Нижняя и верхняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int j=1; j<Ny-1; ++j)
		{
			Mask[i][j][Nz-1].mask = Mask[i][j][0].mask = TBorderPoint;
			Mask[i][j][0].normal = new glTVector(0,0,-1);
			Mask[i][j][Nz-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<Nx; ++i)
	{
		Mask[i][0][0].mask = TFictivePoint;
		Mask[i][0][Nz-1].mask = TFictivePoint;
		Mask[i][Ny-1][0].mask = TFictivePoint;
		Mask[i][Ny-1][Nz-1].mask = TFictivePoint;
	}

	for(int j=0; j<Ny; ++j)
	{
		Mask[0][j][0].mask = TFictivePoint;
		Mask[0][j][Nz-1].mask = TFictivePoint;
		Mask[Nx-1][j][0].mask = TFictivePoint;
		Mask[Nx-1][j][Nz-1].mask = TFictivePoint;
	}

	for(int k=0; k<Nz; ++k)
	{
		Mask[0][0][k].mask = TFictivePoint;
		Mask[0][Ny-1][k].mask = TFictivePoint;
		Mask[Nx-1][0][k].mask = TFictivePoint;
		Mask[Nx-1][Ny-1][k].mask = TFictivePoint;
	}


	// Внутри препятствия нерасчетные точки. Поверхность.
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		{
			int
				ii = i-(Nx1-1),
				jj = j-(Ny1-1),
				q;

			if( ( ii <= jj ) && ( ii <= (Ny2-1)-jj ) && ( ii < Nz2 ) )
			{
				q = ii;
				Mask[i][j][Nz1-1+q].mask = TBorderPoint;
				Mask[i][j][Nz1-1+q].normal = new glTVector(1,0,-1);
			}
			else if( ( (Nx2-1)-ii <= jj ) && ( (Nx2-1)-ii <= (Ny2)-jj ) && ( ii>Nx2-Nz2 ) )
			{
				q = Nx2 - 1 - ii;
				Mask[i][j][Nz1-1+q].mask = TBorderPoint;
				Mask[i][j][Nz1-1+q].normal = new glTVector(-1,0,-1);
			}
			else if( ( ii > jj ) && ( (Nx2-1)-ii > jj ) && ( jj < Nz2 ) )
			{
				q = jj;
				Mask[i][j][Nz1-1+q].mask = TBorderPoint;
				Mask[i][j][Nz1-1+q].normal = new glTVector(0,1,-1);
			}
			else if( ( ii > (Ny2-1)-jj ) && ( (Nx2-1)-ii > (Ny2-1)-jj ) && ( jj>Ny2-Nz2 ) )
			{
				q = Ny2 - 1 - jj;
				Mask[i][j][Nz1-1+q].mask = TBorderPoint;
				Mask[i][j][Nz1-1+q].normal = new glTVector(0,-1,-1);
			}
			else
			{
				q = Nz2;
			}
			for(int k=0; k<Nz1-1+q; ++k)
				Mask[i][j][k].mask = TFictivePoint;
		}


	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
	{
		for(int k=1; k<Nz1; ++k)
		{
			Mask[i][Ny1-1][k].mask = Mask[i][Ny1+Ny2-2][k].mask = TPreNormalBorderPoint;
			Mask[i][Ny1-1][k].normal = new glTVector(0,1,0);
			Mask[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}
	}
		
	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
	{
		for(int k=1; k<Nz1; ++k)
		{
			Mask[Nx1-1][j][k].mask = Mask[Nx1+Nx2-2][j][k].mask = TPreNormalBorderPoint;
			Mask[Nx1-1][j][k].normal = new glTVector(1,0,0);
			Mask[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
	}

}
//*************************************
void TDoubleTrapezeBottomWithPressureMaskGenerator::CreateMultyNodes()
{
	// Узлы для компоненты U скорости
	for(int i=0; i<NxU; ++i)
		XU[i] = X[i];

	for(int j=1; j<NyU-1; ++j)
		YU[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YU[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YU[NyU-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzU-1; ++k)
		ZU[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZU[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZU[NzU-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxU; ++i)
		HxU[i] = XU[i] - XU[i-1];
	HxU[0] = HxU[1];

	for(int j=1; j<NyU; ++j)
		HyU[j] = YU[j] - YU[j-1];
	HyU[0] = HyU[1];

	for(int k=1; k<NzU; ++k)
		HzU[k] = ZU[k] - ZU[k-1];
	HzU[0] = HzU[1];

	// Узлы для компоненты V скорости
	for(int i=1; i<NxV-1; ++i)
		XV[i] = ( X[i-1] + X[i] ) / 2.0;
	XV[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XV[NxV-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=0; j<NyV; ++j)
		YV[j] = Y[j];

	for(int k=1; k<NzV-1; ++k)
		ZV[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZV[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZV[NzV-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxV; ++i)
		HxV[i] = XV[i] - XV[i-1];
	HxV[0] = HxV[1];

	for(int j=1; j<NyV; ++j)
		HyV[j] = YV[j] - YV[j-1];
	HyV[0] = HyV[1];

	for(int k=1; k<NzV; ++k)
		HzV[k] = ZV[k] - ZV[k-1];
	HzV[0] = HzV[1];

	// Узлы для компоненты W скорости
	for(int i=1; i<NxW-1; ++i)
		XW[i] = ( X[i-1] + X[i] ) / 2.0;
	XW[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XW[NxW-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyW-1; ++j)
		YW[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YW[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YW[NyW-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=0; k<NzW; ++k)
		ZW[k] = Z[k];

	for(int i=1; i<NxW; ++i)
		HxW[i] = XW[i] - XW[i-1];
	HxW[0] = HxW[1];

	for(int j=1; j<NyW; ++j)
		HyW[j] = YW[j] - YW[j-1];
	HyW[0] = HyW[1];

	for(int k=1; k<NzW; ++k)
		HzW[k] = ZW[k] - ZW[k-1];
	HzW[0] = HzW[1];

	// Узлы для давления
	for(int i=1; i<NxP-1; ++i)
		XP[i] = ( X[i-1] + X[i] ) / 2.0;
	XP[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XP[NxP-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyP-1; ++j)
		YP[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YP[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YP[NyP-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzP-1; ++k)
		ZP[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZP[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZP[NzP-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxP; ++i)
		HxP[i] = XP[i] - XP[i-1];
	HxP[0] = HxP[1];

	for(int j=1; j<NyP; ++j)
		HyP[j] = YP[j] - YP[j-1];
	HyP[0] = HyP[1];

	for(int k=1; k<NzP; ++k)
		HzP[k] = ZP[k] - ZP[k-1];
	HzP[0] = HzP[1];
}
//*************************************
void TDoubleTrapezeBottomWithPressureMaskGenerator::CreateMultyGrids()
{
	// Заполнение маски для U //

	for(int i=0; i<NxU; ++i)
		for(int j=0; j<NyU; ++j)
			for(int k=0; k<NzU; ++k)
			{
				MaskU[i][j][k].mask = TActualPoint;
				MaskU[i][j][k].normal = 0;
			}

	// Границы области
	//for(int j=1; j<NyU-1; ++j)
	//	for(int k=1; k<NzU-1; ++k)
	//	{
	//		MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TDefinedBorderPoint;
	//		MaskU[0][j][k].normal = new glTVector(-1,0,0);
	//		MaskU[NxU-1][j][k].normal = new glTVector(1,0,0);
	//	}

	for(int j=1; j<NyU-1; ++j)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TNormalBorderPoint;
			MaskU[0][j][k].normal = new glTVector(-1,0,0);
			MaskU[NxU-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[i][0][k].mask = MaskU[i][NyU-1][k].mask = TPreDefinedBorderPoint;
			MaskU[i][0][k].normal =  new glTVector(0,-1,0);
			MaskU[i][NyU-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int j=1; j<NyU-1; ++j)
		{
			MaskU[i][j][0].mask = MaskU[i][j][NzU-1].mask = TPreDefinedBorderPoint;
			MaskU[i][j][0].normal = new glTVector(0,0,-1);
			MaskU[i][j][NzU-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxU; ++i)
	{
		MaskU[i][0][0].mask = TFictivePoint;
		MaskU[i][0][NzU-1].mask = TFictivePoint;
		MaskU[i][NyU-1][0].mask = TFictivePoint;
		MaskU[i][NyU-1][NzU-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyU; ++j)
	{
		MaskU[0][j][0].mask = TFictivePoint;
		MaskU[0][j][NzU-1].mask = TFictivePoint;
		MaskU[NxU-1][j][0].mask = TFictivePoint;
		MaskU[NxU-1][j][NzU-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzU; ++k)
	{
		MaskU[0][0][k].mask = TFictivePoint;
		MaskU[0][NyU-1][k].mask = TFictivePoint;
		MaskU[NxU-1][0][k].mask = TFictivePoint;
		MaskU[NxU-1][NyU-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия нерасчетные точки. Поверхность.
	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			int
				ii = i-(Nx1-1),
				jj = j-(Ny1-1),
				q;

			if( ( ii < jj ) && ( ii <= (Ny2-1)-jj ) && ( ii <= Nz2-2 ) )
			{
				q = ii + 1;
				MaskU[i][j][Nz1-1+q].mask = TDefinedBorderPoint;
				MaskU[i][j][Nz1-1+q].normal = new glTVector(1,0,-1);
			}
			else if( ( (Nx2-1)-ii < jj ) && ( (Nx2-1)-ii <= (Ny2-1)-jj ) && ( ii>=Nx2-Nz2+1 ) )
			{
				q = Nx2 - ii;
				MaskU[i][j][Nz1-1+q].mask = TDefinedBorderPoint;
				MaskU[i][j][Nz1-1+q].normal = new glTVector(-1,0,-1);
			}
			else if( ( ii >= jj ) && ( (Nx2-1)-ii >= jj ) && ( jj <= Nz2-1 ) )
			{
				q = jj;
				MaskU[i][j][Nz1-1+q].mask = TPreDefinedBorderPoint;
				MaskU[i][j][Nz1-1+q].normal = new glTVector(0,1,-1);
			}
			else if( ( ii > (Ny2-1)-jj ) && ( (Nx2-1)-ii > (Ny2-1)-jj ) && ( jj>=Ny2-Nz2+1 ) )
			{
				q = Ny2 - jj;
				MaskU[i][j][Nz1-1+q].mask = TPreDefinedBorderPoint;
				MaskU[i][j][Nz1-1+q].normal = new glTVector(0,-1,-1);
			}
			else
			{
				q = Nz2+1;
			}
			for(int k=0; k<Nz1-1+q; ++k)
				MaskU[i][j][k].mask = TFictivePoint;
		}


	for(int i=Nx1-1; i<Nx1+Nx2-1; ++i)
	{
		for(int k=1; k<Nz1; ++k)
		{
			MaskU[i][Ny1][k].mask = MaskU[i][Ny1+Ny2-2][k].mask = TPreDefinedBorderPoint;
			MaskU[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskU[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}
	}
		
	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	{
		for(int k=1; k<Nz1; ++k)
		{
			MaskU[Nx1-1][j][k].mask = MaskU[Nx1+Nx2-2][j][k].mask = TDefinedBorderPoint;
			MaskU[Nx1-1][j][k].normal = new glTVector(1,0,0);
			MaskU[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
	}

	// Заполнение маски для V //

	for(int i=0; i<NxV; ++i)
		for(int j=0; j<NyV; ++j)
			for(int k=0; k<NzV; ++k)
			{
				MaskV[i][j][k].mask = TActualPoint;
				MaskV[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyV-1; ++j)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[0][j][k].mask = MaskV[NxV-1][j][k].mask = TPreDefinedBorderPoint;
			MaskV[0][j][k].normal = new glTVector(-1,0,0);
			MaskV[NxV-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[i][0][k].mask = MaskV[i][NyV-1][k].mask = TDefinedBorderPoint;
			MaskV[i][0][k].normal =  new glTVector(0,-1,0);
			MaskV[i][NyV-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int j=1; j<NyV-1; ++j)
		{
			MaskV[i][j][0].mask = MaskV[i][j][NzV-1].mask = TPreDefinedBorderPoint;
			MaskV[i][j][0].normal = new glTVector(0,0,-1);
			MaskV[i][j][NzV-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxV; ++i)
	{
		MaskV[i][0][0].mask = TFictivePoint;
		MaskV[i][0][NzV-1].mask = TFictivePoint;
		MaskV[i][NyV-1][0].mask = TFictivePoint;
		MaskV[i][NyV-1][NzV-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyV; ++j)
	{
		MaskV[0][j][0].mask = TFictivePoint;
		MaskV[0][j][NzV-1].mask = TFictivePoint;
		MaskV[NxV-1][j][0].mask = TFictivePoint;
		MaskV[NxV-1][j][NzV-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzV; ++k)
	{
		MaskV[0][0][k].mask = TFictivePoint;
		MaskV[0][NyV-1][k].mask = TFictivePoint;
		MaskV[NxV-1][0][k].mask = TFictivePoint;
		MaskV[NxV-1][NyV-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия нерасчетные точки. Поверхность.
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
		{
			int
				ii = i-(Nx1-1),
				jj = j-(Ny1-1),
				q;

			if( ( ii <= jj ) && ( ii <= (Ny2-1)-jj ) && ( ii <= Nz2-1 ) )
			{
				q = ii;
				MaskV[i][j][Nz1-1+q].mask = TPreDefinedBorderPoint;
				MaskV[i][j][Nz1-1+q].normal = new glTVector(1,0,-1);
			}
			else if( ( (Nx2-1)-ii < jj ) && ( (Nx2-1)-ii < (Ny2-1)-jj ) && ( ii>=Nx2-Nz2+1 ) )
			{
				q = Nx2 - ii;
				MaskV[i][j][Nz1-1+q].mask = TPreDefinedBorderPoint;
				MaskV[i][j][Nz1-1+q].normal = new glTVector(-1,0,-1);
			}
			else if( ( ii > jj ) && ( (Nx2-1)-ii >= jj ) && ( jj <= Nz2-2 ) )
			{
				q = jj+1;
				MaskV[i][j][Nz1-1+q].mask = TDefinedBorderPoint;
				MaskV[i][j][Nz1-1+q].normal = new glTVector(0,1,-1);
			}
			else if( ( ii > (Ny2-1)-jj ) && ( (Nx2-1)-ii >= (Ny2-1)-jj ) && ( jj>=Ny2-Nz2+1 ) )
			{
				q = Ny2 - jj;
				MaskV[i][j][Nz1-1+q].mask = TDefinedBorderPoint;
				MaskV[i][j][Nz1-1+q].normal = new glTVector(0,-1,-1);
			}
			else
			{
				q = Nz2+1;
			}
			for(int k=0; k<Nz1-1+q; ++k)
				MaskV[i][j][k].mask = TFictivePoint;
		}


	for(int j=Ny1-1; j<Ny1+Ny2-1; ++j)
	{
		for(int k=1; k<Nz1; ++k)
		{
			MaskV[Nx1][j][k].mask = MaskV[Nx1+Nx2-2][j][k].mask = TPreDefinedBorderPoint;
			MaskV[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskV[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
	}

	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	{
		for(int k=1; k<Nz1; ++k)
		{
			MaskV[i][Ny1-1][k].mask = MaskV[i][Ny1+Ny2-2][k].mask = TDefinedBorderPoint;
			MaskV[i][Ny1-1][k].normal = new glTVector(0,1,0);
			MaskV[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}
	}
		
	// Заполнение маски для W //

	for(int i=0; i<NxW; ++i)
		for(int j=0; j<NyW; ++j)
			for(int k=0; k<NzW; ++k)
			{
				MaskW[i][j][k].mask = TActualPoint;
				MaskW[i][j][k].normal = 0;
			}
	// Границы области
	for(int j=1; j<NyW-1; ++j)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[0][j][k].mask = MaskW[NxW-1][j][k].mask = TPreDefinedBorderPoint; 
			MaskW[0][j][k].normal = new glTVector(-1,0,0);
			MaskW[NxW-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[i][0][k].mask = MaskW[i][NyW-1][k].mask = TPreDefinedBorderPoint;
			MaskW[i][0][k].normal =  new glTVector(0,-1,0);
			MaskW[i][NyW-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int j=1; j<NyW-1; ++j)
		{
			MaskW[i][j][0].mask = MaskW[i][j][NzW-1].mask = TDefinedBorderPoint;
			MaskW[i][j][0].normal = new glTVector(0,0,-1);
			MaskW[i][j][NzW-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxW; ++i)
	{
		MaskW[i][0][0].mask = TFictivePoint;
		MaskW[i][0][NzW-1].mask = TFictivePoint;
		MaskW[i][NyW-1][0].mask = TFictivePoint;
		MaskW[i][NyW-1][NzW-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyW; ++j)
	{
		MaskW[0][j][0].mask = TFictivePoint;
		MaskW[0][j][NzW-1].mask = TFictivePoint;
		MaskW[NxW-1][j][0].mask = TFictivePoint;
		MaskW[NxW-1][j][NzW-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzW; ++k)
	{
		MaskW[0][0][k].mask = TFictivePoint;
		MaskW[0][NyW-1][k].mask = TFictivePoint;
		MaskW[NxW-1][0][k].mask = TFictivePoint;
		MaskW[NxW-1][NyW-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия нерасчетные точки. Поверхность.
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			int
				ii = i-(Nx1-1),
				jj = j-(Ny1-1),
				q;

			if( ( ii <= jj ) && ( ii <= (Ny2)-jj ) && ( ii < Nz2-1 ) )
			{
				q = ii;
				MaskW[i][j][Nz1-1+q].mask = TDefinedBorderPoint;
				MaskW[i][j][Nz1-1+q].normal = new glTVector(1,0,-1);
			}
			else if( ( (Nx2)-ii <= jj ) && ( (Nx2)-ii <= (Ny2)-jj ) && ( ii>Nx2-Nz2+1 ) )
			{
				q = Nx2 - ii;
				MaskW[i][j][Nz1-1+q].mask = TDefinedBorderPoint;
				MaskW[i][j][Nz1-1+q].normal = new glTVector(-1,0,-1);
			}
			else if( ( ii > jj ) && ( (Nx2)-ii > jj ) && ( jj < Nz2-1 ) )
			{
				q = jj;
				MaskW[i][j][Nz1-1+q].mask = TDefinedBorderPoint;
				MaskW[i][j][Nz1-1+q].normal = new glTVector(0,1,-1);
			}
			else if( ( ii > (Ny2)-jj ) && ( (Nx2)-ii > (Ny2)-jj ) && ( jj>Ny2-Nz2+1 ) )
			{
				q = Ny2 - jj;
				MaskW[i][j][Nz1-1+q].mask = TDefinedBorderPoint;
				MaskW[i][j][Nz1-1+q].normal = new glTVector(0,-1,-1);
			}
			else
			{
				q = Nz2;
			}
			for(int k=0; k<Nz1-1+q; ++k)
				MaskW[i][j][k].mask = TFictivePoint;
		}


	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	{
		for(int k=1; k<Nz1; ++k)
		{
			MaskW[i][Ny1][k].mask = MaskW[i][Ny1+Ny2-2][k].mask = TPreDefinedBorderPoint;
			MaskW[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskW[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}
	}
		
	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	{
		for(int k=1; k<Nz1; ++k)
		{
			MaskW[Nx1][j][k].mask = MaskW[Nx1+Nx2-2][j][k].mask = TPreDefinedBorderPoint;
			MaskW[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskW[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
	}


	// Заполнение маски для P
	for(int i=0; i<NxP; ++i)
		for(int j=0; j<NyP; ++j)
			for(int k=0; k<NzP; ++k)
			{
				MaskP[i][j][k].mask = TActualPoint;
				MaskP[i][j][k].normal = 0;
			}

	// Границы области
	//for(int j=1; j<NyP-1; ++j)
	//	for(int k=1; k<NzP-1; ++k)
	//	{
	//		MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreNormalBorderPoint;
	//		MaskP[0][j][k].normal = new glTVector(-1,0,0);
	//		MaskP[NxP-1][j][k].normal = new glTVector(1,0,0);
	//	}

	for(int j=1; j<NyP-1; ++j)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreDefinedBorderPoint; 
			MaskP[0][j][k].normal = new glTVector(-1,0,0);
			MaskP[NxP-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[i][0][k].mask = MaskP[i][NyP-1][k].mask = TFictivePoint;
			MaskP[i][0][k].normal =  new glTVector(0,-1,0);
			MaskP[i][NyP-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int j=1; j<NyP-1; ++j)
		{
			MaskP[i][j][0].mask = MaskP[i][j][NzP-1].mask = TFictivePoint;
			MaskP[i][j][0].normal = new glTVector(0,0,-1);
			MaskP[i][j][NzP-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxP; ++i)
	{
		MaskP[i][0][0].mask = TFictivePoint;
		MaskP[i][0][NzP-1].mask = TFictivePoint;
		MaskP[i][NyP-1][0].mask = TFictivePoint;
		MaskP[i][NyP-1][NzP-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyP; ++j)
	{
		MaskP[0][j][0].mask = TFictivePoint;
		MaskP[0][j][NzP-1].mask = TFictivePoint;
		MaskP[NxP-1][j][0].mask = TFictivePoint;
		MaskP[NxP-1][j][NzP-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzP; ++k)
	{
		MaskP[0][0][k].mask = TFictivePoint;
		MaskP[0][NyP-1][k].mask = TFictivePoint;
		MaskP[NxP-1][0][k].mask = TFictivePoint;
		MaskP[NxP-1][NyP-1][k].mask = TFictivePoint;
	}

	// Внутри препятствия нерасчетные точки. Поверхность.
	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
		for(int j=Ny1; j<Ny1+Ny2-1; ++j)
		{
			int
				ii = i-(Nx1-1),
				jj = j-(Ny1-1),
				q;

			if( ( ii <= jj ) && ( ii <= (Ny2)-jj ) && ( ii <= Nz2-1 ) )
			{
				q = ii;
				MaskP[i][j][Nz1-1+q].mask = TFictivePoint;
				MaskP[i][j][Nz1-1+q].normal = new glTVector(1,0,-1);
			}
			else if( ( (Nx2)-ii <= jj ) && ( (Nx2)-ii <= (Ny2)-jj ) && ( ii>=Nx2-Nz2+1 ) )
			{
				q = Nx2 - ii;
				MaskP[i][j][Nz1-1+q].mask = TFictivePoint;
				MaskP[i][j][Nz1-1+q].normal = new glTVector(-1,0,-1);
			}
			else if( ( ii > jj ) && ( (Nx2)-ii > jj ) && ( jj <= Nz2-1 ) )
			{
				q = jj;
				MaskP[i][j][Nz1-1+q].mask = TFictivePoint;
				MaskP[i][j][Nz1-1+q].normal = new glTVector(0,1,-1);
			}
			else if( ( ii > (Ny2)-jj ) && ( (Nx2)-ii > (Ny2)-jj ) && ( jj>=Ny2-Nz2+1 ) )
			{
				q = Ny2 - jj;
				MaskP[i][j][Nz1-1+q].mask = TFictivePoint;
				MaskP[i][j][Nz1-1+q].normal = new glTVector(0,-1,-1);
			}
			else
			{
				q = Nz2+1;
			}
			for(int k=0; k<Nz1-1+q; ++k)
				MaskP[i][j][k].mask = TFictivePoint;
		}


	for(int i=Nx1; i<Nx1+Nx2-1; ++i)
	{
		for(int k=1; k<Nz1; ++k)
		{
			MaskP[i][Ny1][k].mask = MaskP[i][Ny1+Ny2-2][k].mask = TFictivePoint;
			MaskP[i][Ny1][k].normal = new glTVector(0,1,0);
			MaskP[i][Ny1+Ny2-2][k].normal = new glTVector(0,-1,0);
		}
	}
		
	for(int j=Ny1; j<Ny1+Ny2-1; ++j)
	{
		for(int k=1; k<Nz1; ++k)
		{
			MaskP[Nx1][j][k].mask = MaskP[Nx1+Nx2-2][j][k].mask = TFictivePoint;
			MaskP[Nx1][j][k].normal = new glTVector(1,0,0);
			MaskP[Nx1+Nx2-2][j][k].normal = new glTVector(-1,0,0);
		}
	}


	{


        /*
		FILE *f = fopen("u.txt", "w");
		for(int k=0; k<NzU; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyU-1; j>=0; --j)
			{
				for(int i=0; i<NxU; ++i)
					fprintf(f, "%4d ", MaskU[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);


		f = fopen("v.txt", "w");
		for(int k=0; k<NzV; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyV-1; j>=0; --j)
			{
				for(int i=0; i<NxV; ++i)
					fprintf(f, "%4d ", MaskV[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);

		f = fopen("w.txt", "w");
		for(int k=0; k<NzW; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyW-1; j>=0; --j)
			{
				for(int i=0; i<NxW; ++i)
					fprintf(f, "%4d ", MaskW[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);

		f = fopen("p.txt", "w");
		for(int k=0; k<NzP; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyP-1; j>=0; --j)
			{
				for(int i=0; i<NxP; ++i)
					fprintf(f, "%4d ", MaskP[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);
        */

	}

}
//*************************************
double TDoubleTrapezeBottomWithPressureMaskGenerator::getUCondition(int i, int j, int k, int n)
{

	if(k==NzU-1)
	{
		return 1;
	}

	return 0;
}
//*************************************
double TDoubleTrapezeBottomWithPressureMaskGenerator::getVCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TDoubleTrapezeBottomWithPressureMaskGenerator::getWCondition(int i, int j, int k, int n)
{
	
	//if ( ((i==Nx1-1)||(i==Nx1))&&(j>Ny1-1)&&(j<Ny1+Ny2-1) )
	//{
	//	return
	//		-0.25;
	//}

	return 0;
}
//*************************************
double TDoubleTrapezeBottomWithPressureMaskGenerator::getPCondition(int i, int j, int k, int n)
{
	if(i==0)
		return PInner;
	else if(i==NxP-1)
		return POuter;
	else
		return 0;
}
//*************************************



//#################################################################



//*************************************
void TBackwardFacingStepMaskGenerator::CreateGeometryNodes()
{
	for(int i=0; i<Nx; ++i)
		X[i] = LengthX*i/(Nx-1);
	for(int i=1; i<Nx; ++i)
		Hx[i] = X[i]-X[i-1];
	Hx[0] = Hx[1];

	for(int j=0; j<Ny; ++j)
		Y[j] = LengthY*j/(Ny-1);
	for(int j=1; j<Ny; ++j)
		Hy[j] = Y[j]-Y[j-1];
	Hy[0] = Hy[1];

	for(int k=0; k<Nz; ++k)
		Z[k] = LengthZ*k/(Nz-1);
	for(int k=1; k<Nz; ++k)
		Hz[k] = Z[k]-Z[k-1];
	Hy[0] = Hy[1];
}
//*************************************
void TBackwardFacingStepMaskGenerator::CreateGeometryGrid()
{
	// Сначала все точки - расчетные
	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz; ++k)
			{
				Mask[i][j][k].mask = TActualPoint;
				Mask[i][j][k].normal = 0;
			}

	// Левая и правая границы
	for(int j=1; j<Ny-1; ++j)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[Nx-1][j][k].mask = Mask[0][j][k].mask = TBorderPoint;
			Mask[0][j][k].normal = new glTVector(-1,0,0);
			Mask[Nx-1][j][k].normal = new glTVector(1,0,0);
		}

	// Ближняя и дальняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int k=1; k<Nz-1; ++k)
		{
			Mask[i][Ny-1][k].mask = Mask[i][0][k].mask = TBorderPoint;
			Mask[i][0][k].normal = new glTVector(0,-1,0);
			Mask[i][Ny-1][k].normal = new glTVector(0,1,0);
		}

	// Нижняя и верхняя границы
	for(int i=1; i<Nx-1; ++i)
		for(int j=1; j<Ny-1; ++j)
		{
			Mask[i][j][Nz-1].mask = Mask[i][j][0].mask = TBorderPoint;
			Mask[i][j][0].normal = new glTVector(0,0,-1);
			Mask[i][j][Nz-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<Nx; ++i)
	{
		Mask[i][0][0].mask = TFictivePoint;
		Mask[i][0][Nz-1].mask = TFictivePoint;
		Mask[i][Ny-1][0].mask = TFictivePoint;
		Mask[i][Ny-1][Nz-1].mask = TFictivePoint;
	}

	for(int j=0; j<Ny; ++j)
	{
		Mask[0][j][0].mask = TFictivePoint;
		Mask[0][j][Nz-1].mask = TFictivePoint;
		Mask[Nx-1][j][0].mask = TFictivePoint;
		Mask[Nx-1][j][Nz-1].mask = TFictivePoint;
	}

	for(int k=0; k<Nz; ++k)
	{
		Mask[0][0][k].mask = TFictivePoint;
		Mask[0][Ny-1][k].mask = TFictivePoint;
		Mask[Nx-1][0][k].mask = TFictivePoint;
		Mask[Nx-1][Ny-1][k].mask = TFictivePoint;
	}

	// Уступ
	for(int i=0; i<Nx1; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz1; ++k)
				Mask[i][j][k].mask = TFictivePoint;

	for(int i=0; i<Nx1; ++i)
		for(int j=0; j<Ny; ++j)
			Mask[i][j][Nz1-1].mask = TBorderPoint;

	for(int j=0; j<Ny; ++j)
		for(int k=0; k<Nz1; ++k)
			Mask[Nx1-1][j][k].mask = TBorderPoint;

}
//*************************************
void TBackwardFacingStepMaskGenerator::CreateMultyNodes()
{
	// Узлы для компоненты U скорости
	for(int i=0; i<NxU; ++i)
		XU[i] = X[i];

	for(int j=1; j<NyU-1; ++j)
		YU[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YU[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YU[NyU-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzU-1; ++k)
		ZU[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZU[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZU[NzU-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxU; ++i)
		HxU[i] = XU[i] - XU[i-1];
	HxU[0] = HxU[1];

	for(int j=1; j<NyU; ++j)
		HyU[j] = YU[j] - YU[j-1];
	HyU[0] = HyU[1];

	for(int k=1; k<NzU; ++k)
		HzU[k] = ZU[k] - ZU[k-1];
	HzU[0] = HzU[1];

	// Узлы для компоненты V скорости
	for(int i=1; i<NxV-1; ++i)
		XV[i] = ( X[i-1] + X[i] ) / 2.0;
	XV[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XV[NxV-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=0; j<NyV; ++j)
		YV[j] = Y[j];

	for(int k=1; k<NzV-1; ++k)
		ZV[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZV[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZV[NzV-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxV; ++i)
		HxV[i] = XV[i] - XV[i-1];
	HxV[0] = HxV[1];

	for(int j=1; j<NyV; ++j)
		HyV[j] = YV[j] - YV[j-1];
	HyV[0] = HyV[1];

	for(int k=1; k<NzV; ++k)
		HzV[k] = ZV[k] - ZV[k-1];
	HzV[0] = HzV[1];

	// Узлы для компоненты W скорости
	for(int i=1; i<NxW-1; ++i)
		XW[i] = ( X[i-1] + X[i] ) / 2.0;
	XW[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XW[NxW-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyW-1; ++j)
		YW[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YW[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YW[NyW-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=0; k<NzW; ++k)
		ZW[k] = Z[k];

	for(int i=1; i<NxW; ++i)
		HxW[i] = XW[i] - XW[i-1];
	HxW[0] = HxW[1];

	for(int j=1; j<NyW; ++j)
		HyW[j] = YW[j] - YW[j-1];
	HyW[0] = HyW[1];

	for(int k=1; k<NzW; ++k)
		HzW[k] = ZW[k] - ZW[k-1];
	HzW[0] = HzW[1];

	// Узлы для давления
	for(int i=1; i<NxP-1; ++i)
		XP[i] = ( X[i-1] + X[i] ) / 2.0;
	XP[0] = X[0] - ( X[1] - X[0] ) / 2.0;
	XP[NxP-1] = X[Nx-1] + ( X[Nx-1] - X[Nx-2] ) / 2.0;

	for(int j=1; j<NyP-1; ++j)
		YP[j] = ( Y[j-1] + Y[j] ) / 2.0;
	YP[0] = Y[0] - ( Y[1] - Y[0] ) / 2.0;
	YP[NyP-1] = Y[Ny-1] + ( Y[Ny-1] - Y[Ny-2] ) / 2.0;

	for(int k=1; k<NzP-1; ++k)
		ZP[k] = ( Z[k-1] + Z[k] ) / 2.0;
	ZP[0] = Z[0] - ( Z[1] - Z[0] ) / 2.0;
	ZP[NzP-1] = Z[Nz-1] + ( Z[Nz-1] - Z[Nz-2] ) / 2.0;

	for(int i=1; i<NxP; ++i)
		HxP[i] = XP[i] - XP[i-1];
	HxP[0] = HxP[1];

	for(int j=1; j<NyP; ++j)
		HyP[j] = YP[j] - YP[j-1];
	HyP[0] = HyP[1];

	for(int k=1; k<NzP; ++k)
		HzP[k] = ZP[k] - ZP[k-1];
	HzP[0] = HzP[1];
}
//*************************************
void TBackwardFacingStepMaskGenerator::CreateMultyGrids()
{
	// Заполнение маски для U //

	for(int i=0; i<NxU; ++i)
		for(int j=0; j<NyU; ++j)
			for(int k=0; k<NzU; ++k)
			{
				MaskU[i][j][k].mask = TActualPoint;
				MaskU[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyU-1; ++j)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TNormalBorderPoint;
			//MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TEquationBorderPoint;
			MaskU[0][j][k].normal = new glTVector(-1,0,0);
			MaskU[NxU-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[i][0][k].mask = MaskU[i][NyU-1][k].mask = TPreDefinedBorderPoint;
			MaskU[i][0][k].normal =  new glTVector(0,-1,0);
			MaskU[i][NyU-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int j=1; j<NyU-1; ++j)
		{
			MaskU[i][j][0].mask = MaskU[i][j][NzU-1].mask = TPreDefinedBorderPoint;
			MaskU[i][j][0].normal = new glTVector(0,0,-1);
			MaskU[i][j][NzU-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxU; ++i)
	{
		MaskU[i][0][0].mask = TFictivePoint;
		MaskU[i][0][NzU-1].mask = TFictivePoint;
		MaskU[i][NyU-1][0].mask = TFictivePoint;
		MaskU[i][NyU-1][NzU-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyU; ++j)
	{
		MaskU[0][j][0].mask = TFictivePoint;
		MaskU[0][j][NzU-1].mask = TFictivePoint;
		MaskU[NxU-1][j][0].mask = TFictivePoint;
		MaskU[NxU-1][j][NzU-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzU; ++k)
	{
		MaskU[0][0][k].mask = TFictivePoint;
		MaskU[0][NyU-1][k].mask = TFictivePoint;
		MaskU[NxU-1][0][k].mask = TFictivePoint;
		MaskU[NxU-1][NyU-1][k].mask = TFictivePoint;
	}

	// Уступ
	for(int i=0; i<Nx1; ++i)
		for(int j=0; j<NyU; ++j)
			for(int k=0; k<Nz1; ++k)
				MaskU[i][j][k].mask = TFictivePoint;
	
	for(int i=1; i<Nx1; ++i)
		for(int j=1; j<NyU-1; ++j)
		{
			MaskU[i][j][Nz1-1].mask = TPreDefinedBorderPoint;
			MaskU[i][j][Nz1-1].normal = new glTVector(0,0,-1);
		}

	for(int j=0; j<NyU; ++j)
		for(int k=0; k<Nz1; ++k)
		{
			MaskU[Nx1-1][j][k].mask = TDefinedBorderPoint;
			MaskU[Nx1-1][j][k].normal = new glTVector(-1,0,0);
		}

	// Заполнение маски для V //

	for(int i=0; i<NxV; ++i)
		for(int j=0; j<NyV; ++j)
			for(int k=0; k<NzV; ++k)
			{
				MaskV[i][j][k].mask = TActualPoint;
				MaskV[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyV-1; ++j)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[0][j][k].mask = MaskV[NxV-1][j][k].mask = TPreDefinedBorderPoint;
			MaskV[0][j][k].normal = new glTVector(-1,0,0);
			MaskV[NxV-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[i][0][k].mask = MaskV[i][NyV-1][k].mask = TDefinedBorderPoint;
			MaskV[i][0][k].normal =  new glTVector(0,-1,0);
			MaskV[i][NyV-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int j=1; j<NyV-1; ++j)
		{
			MaskV[i][j][0].mask = MaskV[i][j][NzV-1].mask = TPreDefinedBorderPoint;
			MaskV[i][j][0].normal = new glTVector(0,0,-1);
			MaskV[i][j][NzV-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxV; ++i)
	{
		MaskV[i][0][0].mask = TFictivePoint;
		MaskV[i][0][NzV-1].mask = TFictivePoint;
		MaskV[i][NyV-1][0].mask = TFictivePoint;
		MaskV[i][NyV-1][NzV-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyV; ++j)
	{
		MaskV[0][j][0].mask = TFictivePoint;
		MaskV[0][j][NzV-1].mask = TFictivePoint;
		MaskV[NxV-1][j][0].mask = TFictivePoint;
		MaskV[NxV-1][j][NzV-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzV; ++k)
	{
		MaskV[0][0][k].mask = TFictivePoint;
		MaskV[0][NyV-1][k].mask = TFictivePoint;
		MaskV[NxV-1][0][k].mask = TFictivePoint;
		MaskV[NxV-1][NyV-1][k].mask = TFictivePoint;
	}

	// Уступ
	for(int i=0; i<Nx1; ++i)
		for(int j=0; j<NyV; ++j)
			for(int k=0; k<Nz1; ++k)
				MaskV[i][j][k].mask = TFictivePoint;
	
	for(int i=1; i<Nx1; ++i)
		for(int j=1; j<NyV-1; ++j)
		{
			MaskV[i][j][Nz1-1].mask = TPreDefinedBorderPoint;
			MaskV[i][j][Nz1-1].normal = new glTVector(0,0,-1);
		}

	for(int j=1; j<NyV-1; ++j)
		for(int k=1; k<Nz1; ++k)
		{
			MaskV[Nx1-1][j][k].mask = TPreDefinedBorderPoint;
			MaskV[Nx1-1][j][k].normal = new glTVector(-1,0,0);
		}

	for(int j=0; j<NyV; ++j)
		MaskV[Nx1-1][j][Nz1-1].normal = new glTVector(-1/sqrt(2.0),0,-1/sqrt(2.0));

	for(int i=0; i<Nx1; ++i)
	{
		MaskV[i][0][Nz1-1].mask = TDefinedBorderPoint;
		MaskV[i][0][Nz1-1].normal = new glTVector(0,-1,0);
		MaskV[i][NyV-1][Nz1-1].mask = TDefinedBorderPoint;
		MaskV[i][NyV-1][Nz1-1].normal = new glTVector(0,1,0);
	}

	// Заполнение маски для W //

	for(int i=0; i<NxW; ++i)
		for(int j=0; j<NyW; ++j)
			for(int k=0; k<NzW; ++k)
			{
				MaskW[i][j][k].mask = TActualPoint;
				MaskW[i][j][k].normal = 0;
			}
	// Границы области
	for(int j=1; j<NyW-1; ++j)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[0][j][k].mask = MaskW[NxW-1][j][k].mask = TPreDefinedBorderPoint; 
			MaskW[0][j][k].normal = new glTVector(-1,0,0);
			MaskW[NxW-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[i][0][k].mask = MaskW[i][NyW-1][k].mask = TPreDefinedBorderPoint;
			MaskW[i][0][k].normal =  new glTVector(0,-1,0);
			MaskW[i][NyW-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int j=1; j<NyW-1; ++j)
		{
			MaskW[i][j][0].mask = MaskW[i][j][NzW-1].mask = TDefinedBorderPoint;
			MaskW[i][j][0].normal = new glTVector(0,0,-1);
			MaskW[i][j][NzW-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxW; ++i)
	{
		MaskW[i][0][0].mask = TFictivePoint;
		MaskW[i][0][NzW-1].mask = TFictivePoint;
		MaskW[i][NyW-1][0].mask = TFictivePoint;
		MaskW[i][NyW-1][NzW-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyW; ++j)
	{
		MaskW[0][j][0].mask = TFictivePoint;
		MaskW[0][j][NzW-1].mask = TFictivePoint;
		MaskW[NxW-1][j][0].mask = TFictivePoint;
		MaskW[NxW-1][j][NzW-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzW; ++k)
	{
		MaskW[0][0][k].mask = TFictivePoint;
		MaskW[0][NyW-1][k].mask = TFictivePoint;
		MaskW[NxW-1][0][k].mask = TFictivePoint;
		MaskW[NxW-1][NyW-1][k].mask = TFictivePoint;
	}

	// Уступ
	for(int i=0; i<Nx1; ++i)
		for(int j=0; j<NyW; ++j)
			for(int k=0; k<Nz1; ++k)
				MaskW[i][j][k].mask = TFictivePoint;
	
	for(int j=1; j<NyW-1; ++j)
		for(int k=1; k<Nz1; ++k)
		{
			MaskW[Nx1-1][j][k].mask = TPreDefinedBorderPoint;
			MaskW[Nx1-1][j][k].normal = new glTVector(-1,0,0);
		}

	for(int i=1; i<Nx1; ++i)
		for(int j=1; j<NyW-1; ++j)
		{
			MaskW[i][j][Nz1-1].mask = TDefinedBorderPoint;
			MaskW[i][j][Nz1-1].normal = new glTVector(0,0,-1);
		}

	// Заполнение маски для P
	for(int i=0; i<NxP; ++i)
		for(int j=0; j<NyP; ++j)
			for(int k=0; k<NzP; ++k)
			{
				MaskP[i][j][k].mask = TActualPoint;
				MaskP[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyP-1; ++j)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreDefinedBorderPoint;
			MaskP[0][j][k].normal = new glTVector(-1,0,0);
			MaskP[NxP-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[i][0][k].mask = MaskP[i][NyP-1][k].mask = TFictivePoint;
			MaskP[i][0][k].normal =  new glTVector(0,-1,0);
			MaskP[i][NyP-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int j=1; j<NyP-1; ++j)
		{
			MaskP[i][j][0].mask = MaskP[i][j][NzP-1].mask = TFictivePoint;
			MaskP[i][j][0].normal = new glTVector(0,0,-1);
			MaskP[i][j][NzP-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxP; ++i)
	{
		MaskP[i][0][0].mask = TFictivePoint;
		MaskP[i][0][NzP-1].mask = TFictivePoint;
		MaskP[i][NyP-1][0].mask = TFictivePoint;
		MaskP[i][NyP-1][NzP-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyP; ++j)
	{
		MaskP[0][j][0].mask = TFictivePoint;
		MaskP[0][j][NzP-1].mask = TFictivePoint;
		MaskP[NxP-1][j][0].mask = TFictivePoint;
		MaskP[NxP-1][j][NzP-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzP; ++k)
	{
		MaskP[0][0][k].mask = TFictivePoint;
		MaskP[0][NyP-1][k].mask = TFictivePoint;
		MaskP[NxP-1][0][k].mask = TFictivePoint;
		MaskP[NxP-1][NyP-1][k].mask = TFictivePoint;
	}

	// Уступ
	for(int i=0; i<Nx1; ++i)
		for(int j=0; j<NyP; ++j)
			for(int k=0; k<Nz1; ++k)
				MaskP[i][j][k].mask = TFictivePoint;
	
	for(int i=0; i<Nx1; ++i)
		for(int j=0; j<NyP; ++j)
		{
			MaskP[i][j][Nz1-1].normal = new glTVector(0,0,-1);
		}

	for(int j=0; j<NyP; ++j)
		for(int k=0; k<Nz1; ++k)
		{
			MaskP[Nx1-1][j][k].normal = new glTVector(-1,0,0);
		}

	for(int j=0; j<NyP; ++j)
		MaskP[Nx1-1][j][Nz1-1].normal = new glTVector(-1/sqrt(2.0),0,-1/sqrt(2.0));

	for(int i=0; i<Nx1; ++i)
	{
		MaskP[i][0][Nz1-1].normal = new glTVector(0,-1/sqrt(2.0),-1/sqrt(2.0));
		MaskP[i][NyP-1][Nz1-1].normal = new glTVector(0,1/sqrt(2.0),-1/sqrt(2.0));
	}


	//Mask modification (Begin)
	for (int funcNum=1;funcNum<5;funcNum++)
	{
		glTMaskWithNormals*** mask = NULL;
	    int Nx = -1;
	    int Ny = -1;
	    int Nz = -1;

		if (funcNum==1) 
		{
			mask = MaskU;
			Nx = NxU; Ny = NyU; Nz = NzU;
		} 	 
		else if (funcNum==2)
        {
			mask = MaskV;
			Nx = NxV; Ny = NyV; Nz = NzV;
		} 	
		else if (funcNum==3)
        {
			mask = MaskW;
			Nx = NxW; Ny = NyW; Nz = NzW;
		} 	
		else if (funcNum==4)
        {
			mask = MaskP;
			Nx = NxP; Ny = NyP; Nz = NzP;
		} 	
		else {throw "Error in Mask generator: Unknown function number.";}

	    for(int i=0; i<Nx; ++i)
		{
		  for(int j=0; j<Ny; ++j)
		  {
			for(int k=0; k<Nz; ++k)
			{
				if((mask[i][j][k]).mask == TPreDefinedBorderPoint)
				{
					int actualValuesCount = 0;
					for (int dir=1;dir<4;dir++)
					{
						for (int sign=0;sign<2;sign++)
						{
							int shift = -1 + sign*2;
							int i__ = i;
							int j__ = j;
							int k__ = k;
							if (dir==1) {i__=i+shift;}
							else if (dir==2) {j__=j+shift;}
							else if (dir==3) {k__=k+shift;}
							else {throw "Error in Mask generator: Unknown direction.";} 
							if ( (i__>=0) && (i__<=Nx-1) && (j__>=0) && (j__<=Ny-1) && (k__>=0) && (k__<=Nz-1) )
							{
						       if ( (mask[i__][j__][k__]).mask==TActualPoint ) 
								{ 
									actualValuesCount = actualValuesCount + 1;
								}
							}
						}//by -1 or 1
					}//by direction

					if (actualValuesCount<=0)
					{
						throw "Error in Mask generator: Actual points number around PreDefined is zero."; 
					}
					else if (actualValuesCount==1)
					{
						//It is OK!
					}
					else if (actualValuesCount==2)
					{
						(mask[i][j][k]).mask = TDefinedBorderPoint;
					}
					else
					{
                       throw "Error in Mask generator: Actual points number around PreDefined more then two.";
					}
				}// if PreDefined
			}// by i
		  }//by j
		}//by k
	}//by function number
   //Mask modification (End)


/*
		FILE *f = fopen("u.txt", "w");
		for(int k=0; k<NzU; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyU-1; j>=0; --j)
			{
				for(int i=0; i<NxU; ++i)
					fprintf(f, "%4d ", MaskU[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);


		f = fopen("v.txt", "w");
		for(int k=0; k<NzV; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyV-1; j>=0; --j)
			{
				for(int i=0; i<NxV; ++i)
					fprintf(f, "%4d ", MaskV[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);

		f = fopen("w.txt", "w");
		for(int k=0; k<NzW; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyW-1; j>=0; --j)
			{
				for(int i=0; i<NxW; ++i)
					fprintf(f, "%4d ", MaskW[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);

		f = fopen("p.txt", "w");
		for(int k=0; k<NzP; ++k)
		{
			fprintf(f, "k = %d\n", k);
			for(int j=NyP-1; j>=0; --j)
			{
				for(int i=0; i<NxP; ++i)
					fprintf(f, "%4d ", MaskP[i][j][k].mask);
				fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		fclose(f);
*/
}
//*************************************
double TBackwardFacingStepMaskGenerator::getUCondition(int i, int j, int k, int n)
{
	//if(k==NzU-1)
		//return 10;
	return 0;	
}
//*************************************
double TBackwardFacingStepMaskGenerator::getVCondition(int i, int j, int k, int n)
{
	return 0;
}
//*************************************
double TBackwardFacingStepMaskGenerator::getWCondition(int i, int j, int k, int n)
{
	//int
	//	WaveCnt = 5,
	//	Nw1 = 7; // Nw1*WaveCnt - (WaveCnt-1) = Nx

	//if(k==NzW-1)
	//	return sin(2*M_PI*i/(Nw1-1));
	return 0;
}
//*************************************
double TBackwardFacingStepMaskGenerator::getPCondition(int i, int j, int k, int n)
{
	if(i==0)
		return PInner;
	else if(i==NxP-1)
		return POuter;
	else
		return 0;
}
//*************************************



double TSimpleWithCylinderPressureMaskGenerator::getPCondition(int i, int j, int k, int n)
{
	//printf("\nHEHEY!\n");
	//int time = (n / 100);
	//int time = (int)(n*fTimeStepValue);
    //if(i==0)
    if(i==0)
    {
		//printf(" >>Statement -> TSimpleWithCylinderPressureMaskGenerator::getPCondition i=%d, j=%d, k=%d, t=%d\n", i, j, k, n);
		//return PLeft*sin(M_PI*time) + 1;
		//return 2.0*sin(M_PI*time) + 1;
		//return sin(M_PI*time) + 2.;
		return PLeft;
    }
    //else if(i==NxP-1)
    else if(i==NxP-1)
    {
		//printf(" >>Statement -> TSimpleWithCylinderPressureMaskGenerator::getPCondition i=%d, j=%d, k=%d, t=%d\n", i, j, k, n);
        //return PRight*sin(M_PI*time - 0.25) + 1;
        //return 2.0*sin(M_PI*time - 0.25) + 1;d
		//return sin(M_PI*time - 0.25) + 1.;
		return PRight;
    }
    else
    {
        return 0;
    }
}

void TSimpleWithCylinderPressureMaskGenerator::CreateMultyGrids()
{
		// Заполнение маски для U //

	for(int i=0; i<NxU; ++i)
		for(int j=0; j<NyU; ++j)
			for(int k=0; k<NzU; ++k)
			{
				MaskU[i][j][k].mask = TActualPoint;
				MaskU[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyU-1; ++j)
		for(int k=1; k<NzU-1; ++k)
		{			
            /*if (pressureCondition(0, j, k, NxU, NyU, NzU))
            {
                MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TDefinedBorderPoint;
            }
            else
            {
                MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TNormalBorderPoint;
            }*/

			MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TDefinedBorderPoint;
			MaskU[0][j][k].normal = new glTVector(-1,0,0);
			MaskU[NxU-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int k=1; k<NzU-1; ++k)
		{
			MaskU[i][0][k].mask = MaskU[i][NyU-1][k].mask = TPreDefinedBorderPoint;
			MaskU[i][0][k].normal =  new glTVector(0,-1,0);
			MaskU[i][NyU-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxU-1; ++i)
		for(int j=1; j<NyU-1; ++j)
		{
			MaskU[i][j][0].mask = MaskU[i][j][NzU-1].mask = TPreDefinedBorderPoint;
			MaskU[i][j][0].normal = new glTVector(0,0,-1);
			MaskU[i][j][NzU-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxU; ++i)
	{
		MaskU[i][0][0].mask = TFictivePoint;
		MaskU[i][0][NzU-1].mask = TFictivePoint;
		MaskU[i][NyU-1][0].mask = TFictivePoint;
		MaskU[i][NyU-1][NzU-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyU; ++j)
	{
		MaskU[0][j][0].mask = TFictivePoint;
		MaskU[0][j][NzU-1].mask = TFictivePoint;
		MaskU[NxU-1][j][0].mask = TFictivePoint;
		MaskU[NxU-1][j][NzU-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzU; ++k)
	{
		MaskU[0][0][k].mask = TFictivePoint;
		MaskU[0][NyU-1][k].mask = TFictivePoint;
		MaskU[NxU-1][0][k].mask = TFictivePoint;
		MaskU[NxU-1][NyU-1][k].mask = TFictivePoint;
	}


	// Заполнение маски для V //

	for(int i=0; i<NxV; ++i)
		for(int j=0; j<NyV; ++j)
			for(int k=0; k<NzV; ++k)
			{
				MaskV[i][j][k].mask = TActualPoint;
				MaskV[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyV-1; ++j)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[0][j][k].mask = MaskV[NxV-1][j][k].mask = TPreDefinedBorderPoint;
			MaskV[0][j][k].normal = new glTVector(-1,0,0);
			MaskV[NxV-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int k=1; k<NzV-1; ++k)
		{
			MaskV[i][0][k].mask = MaskV[i][NyV-1][k].mask = TDefinedBorderPoint;
			MaskV[i][0][k].normal =  new glTVector(0,-1,0);
			MaskV[i][NyV-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxV-1; ++i)
		for(int j=1; j<NyV-1; ++j)
		{
			MaskV[i][j][0].mask = MaskV[i][j][NzV-1].mask = TPreDefinedBorderPoint;
			MaskV[i][j][0].normal = new glTVector(0,0,-1);
			MaskV[i][j][NzV-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxV; ++i)
	{
		MaskV[i][0][0].mask = TFictivePoint;
		MaskV[i][0][NzV-1].mask = TFictivePoint;
		MaskV[i][NyV-1][0].mask = TFictivePoint;
		MaskV[i][NyV-1][NzV-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyV; ++j)
	{
		MaskV[0][j][0].mask = TFictivePoint;
		MaskV[0][j][NzV-1].mask = TFictivePoint;
		MaskV[NxV-1][j][0].mask = TFictivePoint;
		MaskV[NxV-1][j][NzV-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzV; ++k)
	{
		MaskV[0][0][k].mask = TFictivePoint;
		MaskV[0][NyV-1][k].mask = TFictivePoint;
		MaskV[NxV-1][0][k].mask = TFictivePoint;
		MaskV[NxV-1][NyV-1][k].mask = TFictivePoint;
	}


	// Заполнение маски для W //

	for(int i=0; i<NxW; ++i)
		for(int j=0; j<NyW; ++j)
			for(int k=0; k<NzW; ++k)
			{
				MaskW[i][j][k].mask = TActualPoint;
				MaskW[i][j][k].normal = 0;
			}
	// Границы области
	for(int j=1; j<NyW-1; ++j)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[0][j][k].mask = MaskW[NxW-1][j][k].mask = TPreDefinedBorderPoint; 
			MaskW[0][j][k].normal = new glTVector(-1,0,0);
			MaskW[NxW-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int k=1; k<NzW-1; ++k)
		{
			MaskW[i][0][k].mask = MaskW[i][NyW-1][k].mask = TPreDefinedBorderPoint;
			MaskW[i][0][k].normal =  new glTVector(0,-1,0);
			MaskW[i][NyW-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxW-1; ++i)
		for(int j=1; j<NyW-1; ++j)
		{
			MaskW[i][j][0].mask = MaskW[i][j][NzW-1].mask = TDefinedBorderPoint;
			MaskW[i][j][0].normal = new glTVector(0,0,-1);
			MaskW[i][j][NzW-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxW; ++i)
	{
		MaskW[i][0][0].mask = TFictivePoint;
		MaskW[i][0][NzW-1].mask = TFictivePoint;
		MaskW[i][NyW-1][0].mask = TFictivePoint;
		MaskW[i][NyW-1][NzW-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyW; ++j)
	{
		MaskW[0][j][0].mask = TFictivePoint;
		MaskW[0][j][NzW-1].mask = TFictivePoint;
		MaskW[NxW-1][j][0].mask = TFictivePoint;
		MaskW[NxW-1][j][NzW-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzW; ++k)
	{
		MaskW[0][0][k].mask = TFictivePoint;
		MaskW[0][NyW-1][k].mask = TFictivePoint;
		MaskW[NxW-1][0][k].mask = TFictivePoint;
		MaskW[NxW-1][NyW-1][k].mask = TFictivePoint;
	}


	// Заполнение маски для P
	for(int i=0; i<NxP; ++i)
		for(int j=0; j<NyP; ++j)
			for(int k=0; k<NzP; ++k)
			{
				MaskP[i][j][k].mask = TActualPoint;
				MaskP[i][j][k].normal = 0;
			}

	// Границы области
	for(int j=1; j<NyP-1; ++j)
		for(int k=1; k<NzP-1; ++k)
		{
			 MaskP[0][j][k].mask = MaskP[NxP-1][j][k].mask = TPreDefinedBorderPoint;
            if (pressureCondition(0, j, k, NxP, NyP, NzP))
            {
                MaskP[0][j][k].mask =  TFictivePoint;
            }
            else
            {
                MaskP[0][j][k].mask = TPreDefinedBorderPoint;
            }

			MaskP[0][j][k].normal = new glTVector(-1,0,0);
			MaskP[NxP-1][j][k].normal = new glTVector(1,0,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int k=1; k<NzP-1; ++k)
		{
			MaskP[i][0][k].mask = MaskP[i][NyP-1][k].mask = TFictivePoint;
			MaskP[i][0][k].normal =  new glTVector(0,-1,0);
			MaskP[i][NyP-1][k].normal = new glTVector(0,1,0);
		}

	for(int i=1; i<NxP-1; ++i)
		for(int j=1; j<NyP-1; ++j)
		{
			MaskP[i][j][0].mask = MaskP[i][j][NzP-1].mask = TFictivePoint;
			MaskP[i][j][0].normal = new glTVector(0,0,-1);
			MaskP[i][j][NzP-1].normal = new glTVector(0,0,1);
		}

	// Фиктивные линии по ребрам
	for(int i=0; i<NxP; ++i)
	{
		MaskP[i][0][0].mask = TFictivePoint;
		MaskP[i][0][NzP-1].mask = TFictivePoint;
		MaskP[i][NyP-1][0].mask = TFictivePoint;
		MaskP[i][NyP-1][NzP-1].mask = TFictivePoint;
	}

	for(int j=0; j<NyP; ++j)
	{
		MaskP[0][j][0].mask = TFictivePoint;
		MaskP[0][j][NzP-1].mask = TFictivePoint;
		MaskP[NxP-1][j][0].mask = TFictivePoint;
		MaskP[NxP-1][j][NzP-1].mask = TFictivePoint;
	}

	for(int k=0; k<NzP; ++k)
	{
		MaskP[0][0][k].mask = TFictivePoint;
		MaskP[0][NyP-1][k].mask = TFictivePoint;
		MaskP[NxP-1][0][k].mask = TFictivePoint;
		MaskP[NxP-1][NyP-1][k].mask = TFictivePoint;
	}
}

//#################################################################