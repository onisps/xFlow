#ifndef StatementsH
#define StatementsH

#include "XService.h"
#include <functional>
#include <math.h>


//#################################################################


// Глобальные типы данных
struct glTMaskWithNormals
{
	glTVector *normal;
	int mask;
	glTMaskWithNormals(): normal(0), mask(0) {}
};


//#################################################################



// Абстрактный класс генерации масок
class TMaskGenerator
{
public:

	// Область решения
	double
		LengthX, LengthY, LengthZ, // размеры области
		*X, *Y, *Z, // массивы узлов сетки
		*Hx, *Hy, *Hz; // массивы шагов сетки glHx[i] = glX[i] - glX[i-1]
	
	// Координаты и шаги смещенных сеток
	double 
		*XU, *YU, *ZU,
		*HxU, *HyU, *HzU,

		*XV, *YV, *ZV,
		*HxV, *HyV, *HzV,

		*XW, *YW, *ZW,
		*HxW, *HyW, *HzW,

		*XP, *YP, *ZP,
		*HxP, *HyP, *HzP;

	// Размерности области решения
	int
		Nx, Ny, Nz;

	// Размерности смещенных сеток
	int
		NxU, NyU, NzU, 
		NxV, NyV, NzV,
		NxW, NyW, NzW,
		NxP, NyP, NzP;
	
	// Маски
	glTMaskWithNormals
		***Mask, // маска геометрии
		***MaskU, ***MaskV, ***MaskW, ***MaskP; // маски для смещенных сеток
	
	// Функции краевых условий
	virtual double getUCondition(int i, int j, int k, int n) = 0;
	virtual double getVCondition(int i, int j, int k, int n) = 0;
	virtual double getWCondition(int i, int j, int k, int n) = 0;
	virtual double getPCondition(int i, int j, int k, int n) = 0;

	// Деструктор
	~TMaskGenerator()
	{
		// Узлы сетки и шаги
		glFreeVector(X, Nx);
		glFreeVector(Y, Ny);
		glFreeVector(Z, Nz);
		glFreeVector(Hx, Nx);
		glFreeVector(Hy, Ny);
		glFreeVector(Hz, Nz);

		// Маска
		glFreeMatrix(Mask, Nx, Ny, Nz);

		// смещенные сетки
		glFreeVector(XU, Nx);
		glFreeVector(YU, Ny);
		glFreeVector(ZU, Nz);
		glFreeVector(HxU, Nx);
		glFreeVector(HyU, Ny);
		glFreeVector(HzU, Nz);

		glFreeVector(XV, Nx);
		glFreeVector(YV, Ny);
		glFreeVector(ZV, Nz);
		glFreeVector(HxV, Nx);
		glFreeVector(HyV, Ny);
		glFreeVector(HzV, Nz);

		glFreeVector(XW, Nx);
		glFreeVector(YW, Ny);
		glFreeVector(ZW, Nz);
		glFreeVector(HxW, Nx);
		glFreeVector(HyW, Ny);
		glFreeVector(HzW, Nz);

		glFreeVector(XP, Nx);
		glFreeVector(YP, Ny);
		glFreeVector(ZP, Nz);
		glFreeVector(HxP, Nx);
		glFreeVector(HyP, Ny);
		glFreeVector(HzP, Nz);

		glFreeMatrix(MaskU, Nx, Ny, Nz);
		glFreeMatrix(MaskV, Nx, Ny, Nz);
		glFreeMatrix(MaskW, Nx, Ny, Nz);
		glFreeMatrix(MaskP, Nx, Ny, Nz);
	}
protected:
	TMaskGenerator( 
		double lengthX, double lengthY, double lengthZ, 
		int nx, int ny, int nz
		):
			LengthX(lengthX), LengthY(lengthY), LengthZ(lengthZ), 
			Nx(nx), Ny(ny), Nz(nz) {};
	void initialize()
	{
		// Размеры разнесенных сеток
		NxU = NxV = NxW = NxP = Nx + 1;
		NyU = NyV = NyW = NyP = Ny + 1;
		NzU = NzV = NzW = NzP = Nz + 1;
		--NxU;
		--NyV;
		--NzW;

		InitMemory();
		CreateGeometryNodes();
		CreateGeometryGrid();
		CreateMultyNodes();
		CreateMultyGrids();
	}
	virtual void InitMemory();
	virtual void CreateGeometryNodes() = 0;
	virtual void CreateGeometryGrid() = 0;
	virtual void CreateMultyNodes() = 0;
	virtual void CreateMultyGrids() = 0;
};



//#################################################################



// Простой параллелепипед с давлением на левой и правой границах
class TSimpleWithPressureMaskGenerator: public TMaskGenerator
{
public:
	// Краевые условия
	double
		PLeft, PRight;

	// Функции краевых условий
	virtual double getUCondition(int i, int j, int k, int n);
	virtual double getVCondition(int i, int j, int k, int n);
	virtual double getWCondition(int i, int j, int k, int n);
	virtual double getPCondition(int i, int j, int k, int n);

	// Фабрика
	static TSimpleWithPressureMaskGenerator *CreateInstance(
		double lengthX, double lengthY, double lengthZ, 
		int nx, int ny, int nz,
		double pLeft, double pRight
		)
	{
		TSimpleWithPressureMaskGenerator *t = new TSimpleWithPressureMaskGenerator(lengthX, lengthY, lengthZ, nx, ny, nz, pLeft, pRight);
		t->initialize();
		return t;
	}
protected:
	TSimpleWithPressureMaskGenerator(
		double lengthX, double lengthY, double lengthZ, 
		int nx, int ny, int nz, 
		double pLeft, double pRight
		):
			TMaskGenerator(lengthX, lengthY, lengthZ, nx, ny, nz),
			PLeft(pLeft), PRight(pRight) {}
	
	virtual void CreateGeometryNodes();
	virtual void CreateGeometryGrid();
	virtual void CreateMultyNodes();
	virtual void CreateMultyGrids();
};


//#################################################################



// Параллелепипед с кубиком внутри
class TCubeInsideWithPressureMaskGenerator: public TMaskGenerator
{
public:
	// Краевые условия
	double
		PInner, POuter;

	// Вспомогательные переменные
	int // размеры области
		Nx1, Nx2, Nx3,
		Ny1, Ny2, Ny3,
		Nz1, Nz2, Nz3;

	// Функции краевых условий
	virtual double getUCondition(int i, int j, int k, int n);
	virtual double getVCondition(int i, int j, int k, int n);
	virtual double getWCondition(int i, int j, int k, int n);
	virtual double getPCondition(int i, int j, int k, int n);

	// Фабрика
	static TCubeInsideWithPressureMaskGenerator *CreateInstance(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2, int nx3,
		int ny1, int ny2, int ny3,
		int nz1, int nz2, int nz3,
		double pInner, double pOuter
		)
	{
		TCubeInsideWithPressureMaskGenerator *t = new TCubeInsideWithPressureMaskGenerator(lengthX, lengthY, lengthZ, nx1, nx2, nx3, ny1, ny2, ny3, nz1, nz2, nz3, pInner, pOuter);
		t->initialize();
		return t;
	}
protected:
	TCubeInsideWithPressureMaskGenerator(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2, int nx3,
		int ny1, int ny2, int ny3,
		int nz1, int nz2, int nz3,
		double pInner, double pOuter
		):
			TMaskGenerator(lengthX, lengthY, lengthZ, nx1+nx2+nx3-2, ny1+ny2+ny3-2, nz1+nz2+nz3-2),
			Nx1(nx1), Nx2(nx2), Nx3(nx3),
			Ny1(ny1), Ny2(ny2), Ny3(ny3),
			Nz1(nz1), Nz2(nz2), Nz3(nz3),
			PInner(pInner), POuter(pOuter) {}
	
	virtual void CreateGeometryNodes();
	virtual void CreateGeometryGrid();
	virtual void CreateMultyNodes();
	virtual void CreateMultyGrids();
};



//#################################################################



// Два концентрических параллелепипеда с давлением на нижней и верхней границах
class TTwoParallelepipedsWithPressureMaskGenerator: public TMaskGenerator
{
public:
	// Краевые условия
	double
		PInner, POuter;

	// Вспомогательные переменные
	int // размеры области
		Nx1, Nx2, Nx3, Nx4, Nx5, Nx6, Nx7,
		Ny1, Ny2, Ny3, Ny4, Ny5, Ny6, Ny7,
		Nz1, Nz2, Nz3, Nz4;

	// Функции краевых условий
	virtual double getUCondition(int i, int j, int k, int n);
	virtual double getVCondition(int i, int j, int k, int n);
	virtual double getWCondition(int i, int j, int k, int n);
	virtual double getPCondition(int i, int j, int k, int n);

	// Фабрика
	static TTwoParallelepipedsWithPressureMaskGenerator *CreateInstance(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2, int nx3, int nx4, int nx5, int nx6, int nx7,
		int ny1, int ny2, int ny3, int ny4, int ny5, int ny6, int ny7,
		int nz, int nz1, int nz2, int nz3, int nz4,
		double pInner, double pOuter
		)
	{
		TTwoParallelepipedsWithPressureMaskGenerator *t = new TTwoParallelepipedsWithPressureMaskGenerator(lengthX, lengthY, lengthZ, nx1, nx2, nx3, nx4, nx5, nx6, nx7, ny1, ny2, ny3, ny4, ny5, ny6, ny7, nz, nz1, nz2, nz3, nz4, pInner, pOuter);
		t->initialize();
		return t;
	}
protected:
	TTwoParallelepipedsWithPressureMaskGenerator(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2, int nx3, int nx4, int nx5, int nx6, int nx7,
		int ny1, int ny2, int ny3, int ny4, int ny5, int ny6, int ny7,
		int nz, int nz1, int nz2, int nz3, int nz4,
		double pInner, double pOuter
		):
			TMaskGenerator(lengthX, lengthY, lengthZ, nx1+nx2+nx3-2, ny1+ny2+ny3-2, nz),
			Nx1(nx1), Nx2(nx2), Nx3(nx3), Nx4(nx4), Nx5(nx5), Nx6(nx6), Nx7(nx7),
			Ny1(ny1), Ny2(ny2), Ny3(ny3), Ny4(ny4), Ny5(ny5), Ny6(ny6), Ny7(ny7),
			Nz1(nz1), Nz2(nz2), Nz3(nz3), Nz4(nz4), 
			PInner(pInner), POuter(pOuter) {}
	
	virtual void CreateGeometryNodes();
	virtual void CreateGeometryGrid();
	virtual void CreateMultyNodes();
	virtual void CreateMultyGrids();
};



//#################################################################



// Параллелепипед с кубиком на нижней грани
class TCubeBottomWithPressureMaskGenerator: public TMaskGenerator
{
public:
	// Краевые условия
	double
		PInner, POuter;

	// Вспомогательные переменные
	int // размеры области
		Nx1, Nx2, Nx3,
		Ny1, Ny2, Ny3,
		Nz1, Nz2;

	// Функции краевых условий
	virtual double getUCondition(int i, int j, int k, int n);
	virtual double getVCondition(int i, int j, int k, int n);
	virtual double getWCondition(int i, int j, int k, int n);
	virtual double getPCondition(int i, int j, int k, int n);

	// Фабрика
	static TCubeBottomWithPressureMaskGenerator *CreateInstance(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2, int nx3,
		int ny1, int ny2, int ny3,
		int nz1, int nz2,
		double pInner, double pOuter
		)
	{
		TCubeBottomWithPressureMaskGenerator *t = new TCubeBottomWithPressureMaskGenerator(lengthX, lengthY, lengthZ, nx1, nx2, nx3, ny1, ny2, ny3, nz1, nz2, pInner, pOuter);
		t->initialize();
		return t;
	}
protected:
	TCubeBottomWithPressureMaskGenerator(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2, int nx3,
		int ny1, int ny2, int ny3,
		int nz1, int nz2,
		double pInner, double pOuter
		):
			TMaskGenerator(lengthX, lengthY, lengthZ, nx1+nx2+nx3-2, ny1+ny2+ny3-2, nz1+nz2-1),
			Nx1(nx1), Nx2(nx2), Nx3(nx3),
			Ny1(ny1), Ny2(ny2), Ny3(ny3),
			Nz1(nz1), Nz2(nz2),
			PInner(pInner), POuter(pOuter) {}
	
	virtual void CreateGeometryNodes();
	virtual void CreateGeometryGrid();
	virtual void CreateMultyNodes();
	virtual void CreateMultyGrids();
};




//#################################################################



// Параллелепипед с трапецией на нижней грани
class TTrapezeBottomWithPressureMaskGenerator: public TMaskGenerator
{
public:
	// Краевые условия
	double
		PInner, POuter;

	// Вспомогательные переменные
	int // размеры области
		Nx1, Nx2, Nx3,
		Ny1, Ny2, Ny3,
		Nz1, Nz2;

	// Функции краевых условий
	virtual double getUCondition(int i, int j, int k, int n);
	virtual double getVCondition(int i, int j, int k, int n);
	virtual double getWCondition(int i, int j, int k, int n);
	virtual double getPCondition(int i, int j, int k, int n);

	// Фабрика
	static TTrapezeBottomWithPressureMaskGenerator *CreateInstance(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2, int nx3,
		int ny1, int ny2, int ny3,
		int nz1, int nz2,
		double pInner, double pOuter
		)
	{
		TTrapezeBottomWithPressureMaskGenerator *t = new TTrapezeBottomWithPressureMaskGenerator(lengthX, lengthY, lengthZ, nx1, nx2, nx3, ny1, ny2, ny3, nz1, nz2, pInner, pOuter);
		t->initialize();
		return t;
	}

	static TTrapezeBottomWithPressureMaskGenerator *CreateInstance(
		double lengthX, double lengthY, double lengthZ, 
		double bottomSide, double topSide, double straightPartHeight, double trapezoidPartHeight,
		double centerX, double centerY,
		double hX, double hY, double hZ,
		double pInner, double pOuter
		)
	{
		// Задаем угол 45 градусов и высоту области
		trapezoidPartHeight = (bottomSide-topSide)/2;
		lengthZ = straightPartHeight + trapezoidPartHeight;
		
		// Вычисляем расстояния
		double
			lx1 = centerX - bottomSide / 2,
			lx2 = bottomSide,
			lx3 = lengthX - centerX - bottomSide / 2,
			ly1 = centerY - bottomSide / 2,
			ly2 = bottomSide,
			ly3 = lengthY - centerY - bottomSide / 2,
			lz1 = straightPartHeight,
			lz2 = trapezoidPartHeight;
		
		// Вычисляем число узлов
		int
			nx1 = int( ceil( lx1 / hX ) + 0.1 )+1,
			nx2 = int( ceil( lx2 / hX ) + 0.1)+1,
			nx3 = int( ceil( lx3 / hX ) + 0.1 )+1,
			ny1 = int( ceil( ly1 / hY ) + 0.1 )+1,
			ny2 = int( ceil( ly2 / hY ) + 0.1 )+1,
			ny3 = int( ceil( ly3 / hY ) + 0.1 )+1,
			nz1 = int( ceil( lz1 / hZ ) + 0.1 )+1,
			nz2 = int( ceil( lz2 / hZ ) + 0.1 )+1;

		// Проверка
		if(ny2!=nx2)
			ny2 = nx2;

		// Конструктор
		TTrapezeBottomWithPressureMaskGenerator *t = new TTrapezeBottomWithPressureMaskGenerator(lengthX, lengthY, lengthZ, nx1, nx2, nx3, ny1, ny2, ny3, nz1, nz2, pInner, pOuter);
		t->initialize();
		return t;
	}

protected:
	TTrapezeBottomWithPressureMaskGenerator(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2, int nx3,
		int ny1, int ny2, int ny3,
		int nz1, int nz2,
		double pInner, double pOuter
		):
			TMaskGenerator(lengthX, lengthY, lengthZ, nx1+nx2+nx3-2, ny1+ny2+ny3-2, nz1+nz2-1),
			Nx1(nx1), Nx2(nx2), Nx3(nx3),
			Ny1(ny1), Ny2(ny2), Ny3(ny3),
			Nz1(nz1), Nz2(nz2),
			PInner(pInner), POuter(pOuter) {}
	
	virtual void CreateGeometryNodes();
	virtual void CreateGeometryGrid();
	virtual void CreateMultyNodes();
	virtual void CreateMultyGrids();
};




//#################################################################



// Параллелепипед с перевернутой трапецией на нижней грани
class TInvertedTrapezeBottomAlongXWithPressureMaskGenerator: public TMaskGenerator
{
public:
	// Краевые условия
	double
		PInner, POuter;

	// Вспомогательные переменные
	int // размеры области
		Nx1, Nx2, Nx3,
		Ny1, Ny2, Ny3;

	// Функции краевых условий
	virtual double getUCondition(int i, int j, int k, int n);
	virtual double getVCondition(int i, int j, int k, int n);
	virtual double getWCondition(int i, int j, int k, int n);
	virtual double getPCondition(int i, int j, int k, int n);

	// Фабрика
	static TInvertedTrapezeBottomAlongXWithPressureMaskGenerator *CreateInstance(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2, int nx3,
		int ny1, int ny2, int ny3,
		int nz,
		double pInner, double pOuter
		)
	{
		TInvertedTrapezeBottomAlongXWithPressureMaskGenerator *t = new TInvertedTrapezeBottomAlongXWithPressureMaskGenerator(lengthX, lengthY, lengthZ, nx1, nx2, nx3, ny1, ny2, ny3, nz, pInner, pOuter);
		t->initialize();
		return t;
	}

	static TInvertedTrapezeBottomAlongXWithPressureMaskGenerator *CreateInstance(
		double lengthX, double lengthY, double lengthZ, 
		double bargeWidth, double bottomLength, double topLength, double bargeHeight,
		double centerX, double centerY,
		double hX, double hY, double hZ,
		double pInner, double pOuter
		)
	{
		// Задаем угол 45 градусов и высоту области
		bargeHeight = (topLength-bottomLength)/2;
		lengthZ = bargeHeight;
		
		// Вычисляем расстояния
		double
			lx1 = centerX - topLength / 2,
			lx2 = topLength,
			lx3 = lengthX - centerX - topLength / 2,
			ly1 = centerY - bargeWidth / 2,
			ly2 = bargeWidth,
			ly3 = lengthY - centerY - bargeWidth / 2;
		
		// Вычисляем число узлов
		int
			nx1 = int( ceil( lx1 / hX ) + 0.1 )+1,
			nx2 = int( ceil( lx2 / hX ) + 0.1 )+1,
			nx3 = int( ceil( lx3 / hX ) + 0.1 )+1,
			ny1 = int( ceil( ly1 / hY ) + 0.1 )+1,
			ny2 = int( ceil( ly2 / hY ) + 0.1 )+1,
			ny3 = int( ceil( ly3 / hY ) + 0.1 )+1,
			nz = int( ceil( lengthZ / hZ ) + 0.1 )+1;

		// Конструктор
		TInvertedTrapezeBottomAlongXWithPressureMaskGenerator *t = new TInvertedTrapezeBottomAlongXWithPressureMaskGenerator(lengthX, lengthY, lengthZ, nx1, nx2, nx3, ny1, ny2, ny3, nz, pInner, pOuter);
		t->initialize();
		return t;
	}

protected:
	TInvertedTrapezeBottomAlongXWithPressureMaskGenerator(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2, int nx3,
		int ny1, int ny2, int ny3,
		int nz,
		double pInner, double pOuter
		):
			TMaskGenerator(lengthX, lengthY, lengthZ, nx1+nx2+nx3-2, ny1+ny2+ny3-2, nz),
			Nx1(nx1), Nx2(nx2), Nx3(nx3),
			Ny1(ny1), Ny2(ny2), Ny3(ny3),
			PInner(pInner), POuter(pOuter) {}
	
	virtual void CreateGeometryNodes();
	virtual void CreateGeometryGrid();
	virtual void CreateMultyNodes();
	virtual void CreateMultyGrids();
};




//#################################################################



// Параллелепипед с перевернутой трапецией на нижней грани
class TInvertedTrapezeBottomAlongYWithPressureMaskGenerator: public TMaskGenerator
{
public:
	// Краевые условия
	double
		PInner, POuter;

	// Вспомогательные переменные
	int // размеры области
		Nx1, Nx2, Nx3,
		Ny1, Ny2, Ny3;

	// Функции краевых условий
	virtual double getUCondition(int i, int j, int k, int n);
	virtual double getVCondition(int i, int j, int k, int n);
	virtual double getWCondition(int i, int j, int k, int n);
	virtual double getPCondition(int i, int j, int k, int n);

	// Фабрика
	static TInvertedTrapezeBottomAlongYWithPressureMaskGenerator *CreateInstance(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2, int nx3,
		int ny1, int ny2, int ny3,
		int nz,
		double pInner, double pOuter
		)
	{
		TInvertedTrapezeBottomAlongYWithPressureMaskGenerator *t = new TInvertedTrapezeBottomAlongYWithPressureMaskGenerator(lengthX, lengthY, lengthZ, nx1, nx2, nx3, ny1, ny2, ny3, nz, pInner, pOuter);
		t->initialize();
		return t;
	}

	static TInvertedTrapezeBottomAlongYWithPressureMaskGenerator *CreateInstance(
		double lengthX, double lengthY, double lengthZ, 
		double bargeWidth, double bottomLength, double topLength, double bargeHeight,
		double centerX, double centerY,
		double hX, double hY, double hZ,
		double pInner, double pOuter
		)
	{
		// Задаем угол 45 градусов и высоту области
		bargeHeight = (topLength-bottomLength)/2;
		lengthZ = bargeHeight;
		
		// Вычисляем расстояния
		double
			lx1 = centerX - bargeWidth / 2,
			lx2 = bargeWidth,
			lx3 = lengthX - centerX - bargeWidth / 2,
			ly1 = centerY - topLength / 2,
			ly2 = topLength,
			ly3 = lengthY - centerY - topLength / 2;
		
		// Вычисляем число узлов
		int
			nx1 = int( ceil( lx1 / hX ) + 0.1 )+1,
			nx2 = int( ceil( lx2 / hX ) + 0.1 )+1,
			nx3 = int( ceil( lx3 / hX ) + 0.1 )+1,
			ny1 = int( ceil( ly1 / hY ) + 0.1 )+1,
			ny2 = int( ceil( ly2 / hY ) + 0.1 )+1,
			ny3 = int( ceil( ly3 / hY ) + 0.1 )+1,
			nz = int( ceil( lengthZ / hZ ) + 0.1 )+1;

		// Конструктор
		TInvertedTrapezeBottomAlongYWithPressureMaskGenerator *t = new TInvertedTrapezeBottomAlongYWithPressureMaskGenerator(lengthX, lengthY, lengthZ, nx1, nx2, nx3, ny1, ny2, ny3, nz, pInner, pOuter);
		t->initialize();
		return t;
	}

protected:
	TInvertedTrapezeBottomAlongYWithPressureMaskGenerator(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2, int nx3,
		int ny1, int ny2, int ny3,
		int nz,
		double pInner, double pOuter
		):
			TMaskGenerator(lengthX, lengthY, lengthZ, nx1+nx2+nx3-2, ny1+ny2+ny3-2, nz),
			Nx1(nx1), Nx2(nx2), Nx3(nx3),
			Ny1(ny1), Ny2(ny2), Ny3(ny3),
			PInner(pInner), POuter(pOuter) {}
	
	virtual void CreateGeometryNodes();
	virtual void CreateGeometryGrid();
	virtual void CreateMultyNodes();
	virtual void CreateMultyGrids();
};


//#################################################################


// Параллелепипед с кубиком внутри (новая дивергенция)
class TCubeInsideNewDivWithPressureMaskGenerator: public TMaskGenerator
{
public:
	// Краевые условия
	double
		PInner, POuter;

	// Вспомогательные переменные
	int // размеры области
		Nx1, Nx2, Nx3,
		Ny1, Ny2, Ny3,
		Nz1, Nz2, Nz3;

	// Функции краевых условий
	virtual double getUCondition(int i, int j, int k, int n);
	virtual double getVCondition(int i, int j, int k, int n);
	virtual double getWCondition(int i, int j, int k, int n);
	virtual double getPCondition(int i, int j, int k, int n);

	// Фабрика
	static TCubeInsideNewDivWithPressureMaskGenerator *CreateInstance(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2, int nx3,
		int ny1, int ny2, int ny3,
		int nz1, int nz2, int nz3,
		double pInner, double pOuter
		)
	{
		TCubeInsideNewDivWithPressureMaskGenerator *t = new TCubeInsideNewDivWithPressureMaskGenerator(lengthX, lengthY, lengthZ, nx1, nx2, nx3, ny1, ny2, ny3, nz1, nz2, nz3, pInner, pOuter);
		t->initialize();
		return t;
	}
protected:
	TCubeInsideNewDivWithPressureMaskGenerator(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2, int nx3,
		int ny1, int ny2, int ny3,
		int nz1, int nz2, int nz3,
		double pInner, double pOuter
		):
			TMaskGenerator(lengthX, lengthY, lengthZ, nx1+nx2+nx3-2, ny1+ny2+ny3-2, nz1+nz2+nz3-2),
			Nx1(nx1), Nx2(nx2), Nx3(nx3),
			Ny1(ny1), Ny2(ny2), Ny3(ny3),
			Nz1(nz1), Nz2(nz2), Nz3(nz3),
			PInner(pInner), POuter(pOuter) {}
	
	virtual void CreateGeometryNodes();
	virtual void CreateGeometryGrid();
	virtual void CreateMultyNodes();
	virtual void CreateMultyGrids();
};


//#################################################################


// Параллелепипед со скошенным кубиком на нижней грани
class TDoubleTrapezeBottomWithPressureMaskGenerator: public TMaskGenerator
{
public:
	// Краевые условия
	double
		PInner, POuter;

	// Вспомогательные переменные
	int // размеры области
		Nx1, Nx2, Nx3,
		Ny1, Ny2, Ny3,
		Nz1, Nz2;

	// Функции краевых условий
	virtual double getUCondition(int i, int j, int k, int n);
	virtual double getVCondition(int i, int j, int k, int n);
	virtual double getWCondition(int i, int j, int k, int n);
	virtual double getPCondition(int i, int j, int k, int n);

	// Фабрика
	static TDoubleTrapezeBottomWithPressureMaskGenerator *CreateInstance(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2, int nx3,
		int ny1, int ny2, int ny3,
		int nz1, int nz2,
		double pInner, double pOuter
		)
	{
		TDoubleTrapezeBottomWithPressureMaskGenerator *t = new TDoubleTrapezeBottomWithPressureMaskGenerator(lengthX, lengthY, lengthZ, nx1, nx2, nx3, ny1, ny2, ny3, nz1, nz2, pInner, pOuter);
		t->initialize();
		return t;
	}

	static TDoubleTrapezeBottomWithPressureMaskGenerator *CreateInstance(
		double lengthX, double lengthY, double lengthZ, 
		double bottomSide, double topSide, double straightPartHeight, double trapezoidPartHeight,
		double centerX, double centerY,
		double hX, double hY, double hZ,
		double pInner, double pOuter
		)
	{
		// Задаем угол 45 градусов и высоту области
		trapezoidPartHeight = (bottomSide-topSide)/2;
		lengthZ = straightPartHeight + trapezoidPartHeight;
		
		// Вычисляем расстояния
		double
			lx1 = centerX - bottomSide / 2,
			lx2 = bottomSide,
			lx3 = lengthX - centerX - bottomSide / 2,
			ly1 = centerY - bottomSide / 2,
			ly2 = bottomSide,
			ly3 = lengthY - centerY - bottomSide / 2,
			lz1 = straightPartHeight,
			lz2 = trapezoidPartHeight;
		
		// Вычисляем число узлов
		int
			nx1 = int( ceil( lx1 / hX ) + 0.1 )+1,
			nx2 = int( ceil( lx2 / hX ) + 0.1)+1,
			nx3 = int( ceil( lx3 / hX ) + 0.1 )+1,
			ny1 = int( ceil( ly1 / hY ) + 0.1 )+1,
			ny2 = int( ceil( ly2 / hY ) + 0.1 )+1,
			ny3 = int( ceil( ly3 / hY ) + 0.1 )+1,
			nz1 = int( ceil( lz1 / hZ ) + 0.1 )+1,
			nz2 = int( ceil( lz2 / hZ ) + 0.1 )+1;

		// Проверка
		if(ny2!=nx2)
			ny2 = nx2;

		// Конструктор
		TDoubleTrapezeBottomWithPressureMaskGenerator *t = new TDoubleTrapezeBottomWithPressureMaskGenerator(lengthX, lengthY, lengthZ, nx1, nx2, nx3, ny1, ny2, ny3, nz1, nz2, pInner, pOuter);
		t->initialize();
		return t;
	}

protected:
	TDoubleTrapezeBottomWithPressureMaskGenerator(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2, int nx3,
		int ny1, int ny2, int ny3,
		int nz1, int nz2,
		double pInner, double pOuter
		):
			TMaskGenerator(lengthX, lengthY, lengthZ, nx1+nx2+nx3-2, ny1+ny2+ny3-2, nz1+nz2-1),
			Nx1(nx1), Nx2(nx2), Nx3(nx3),
			Ny1(ny1), Ny2(ny2), Ny3(ny3),
			Nz1(nz1), Nz2(nz2),
			PInner(pInner), POuter(pOuter) {}
	
	virtual void CreateGeometryNodes();
	virtual void CreateGeometryGrid();
	virtual void CreateMultyNodes();
	virtual void CreateMultyGrids();
};


//#################################################################



// Параллелепипед с уступом внутри (новая дивергенция)
class TBackwardFacingStepMaskGenerator: public TMaskGenerator
{
public:
	// Краевые условия
	double
		PInner, POuter;

	// Вспомогательные переменные
	int // размеры области
		Nx1, Nx2,
		Nz1, Nz2;

	// Функции краевых условий
	virtual double getUCondition(int i, int j, int k, int n);
	virtual double getVCondition(int i, int j, int k, int n);
	virtual double getWCondition(int i, int j, int k, int n);
	virtual double getPCondition(int i, int j, int k, int n);

	// Фабрика
	static TBackwardFacingStepMaskGenerator *CreateInstance(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2,
		int ny,
		int nz1, int nz2,
		double pInner, double pOuter
		)
	{
		TBackwardFacingStepMaskGenerator *t = new TBackwardFacingStepMaskGenerator(lengthX, lengthY, lengthZ, nx1, nx2, ny, nz1, nz2, pInner, pOuter);
		t->initialize();
		return t;
	}
protected:
	TBackwardFacingStepMaskGenerator(
		double lengthX, double lengthY, double lengthZ, 
		int nx1, int nx2,
		int ny,
		int nz1, int nz2,
		double pInner, double pOuter
		):
			TMaskGenerator(lengthX, lengthY, lengthZ, nx1+nx2-1, ny, nz1+nz2-1),
			Nx1(nx1), Nx2(nx2),
			Nz1(nz1), Nz2(nz2),
			PInner(pInner), POuter(pOuter) {}
	
	virtual void CreateGeometryNodes();
	virtual void CreateGeometryGrid();
	virtual void CreateMultyNodes();
	virtual void CreateMultyGrids();
};

//#################################################################


class TSimpleWithCylinderPressureMaskGenerator: public TSimpleWithPressureMaskGenerator
{
private:
    std::function<bool(int, int, int, int, int, int)> pressureCondition;

public:
    void CreateMultyGrids();
	double getPCondition(int i, int j, int k, int n);
	static TSimpleWithCylinderPressureMaskGenerator *CreateInstance(
		double lengthX, double lengthY, double lengthZ, 
		int nx, int ny, int nz,
		double pLeft, double pRight,
        std::function<bool(int, int, int, int, int, int)> pressureCondition
		)
	{
		TSimpleWithCylinderPressureMaskGenerator *t = new TSimpleWithCylinderPressureMaskGenerator(
                lengthX,
                lengthY,
                lengthZ,
                nx, ny, nz,
                pLeft, pRight,
                pressureCondition 
        );
		t->initialize();
		return t;
	}

	TSimpleWithCylinderPressureMaskGenerator(
		double lengthX, double lengthY, double lengthZ, 
		int nx, int ny, int nz, 
		double pLeft, double pRight,
        std::function<bool(int, int, int, int, int, int)> pressureCondition
		): TSimpleWithPressureMaskGenerator(lengthX, lengthY, lengthZ, nx, ny, nz, pLeft, pRight), pressureCondition(pressureCondition) {}
};

#endif
