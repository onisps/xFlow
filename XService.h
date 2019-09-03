#ifndef ServicesH
#define ServicesH


//#################################################################

// Работа с памятью
template <typename T> void glAllocVector(T *&arr, int n)
{
	arr = new T [n];
}
template <typename T> void glFreeVector(T *&arr, int n)
{
	delete[] arr;// arr = NULL;
}
template <typename T> void glAllocPlane(T **&arr, int n, int m)
{
	arr = new T* [n];
	for(int i=0; i<n; ++i)
	{
		arr[i] = new T [m];
	}
}
template <typename T> void glFreePlane(T **&arr, int n, int m)
{
	for(int i=0; i<n; ++i)
		delete [] arr[i];
	delete[] arr;// arr = NULL;
}
template <typename T> void glAllocMatrix(T ***&arr, int n, int m, int l)
 {
	arr = new T** [n];
	for(int i=0; i<n; ++i)
	{
		arr[i] = new T* [m];
		for(int j=0; j<m; ++j)
		{
			arr[i][j] = new T [l];
		}
	}
}
template <typename T> void glFreeMatrix(T ***&arr, int n, int m, int l)
{
	for(int i=0; i<n; ++i)
	{
		for(int j=0; j<m; ++j)
			delete [] arr[i][j];
		delete [] arr[i];
	}
	delete [] arr;
	//arr = NULL;
}


//#################################################################

// Глобальные типы данных
struct glTVector
{
	double x,y,z;
	glTVector() {}
	glTVector(double a, double b, double c): x(a), y(b), z(c) {}
	void set(double a, double b, double c)
	{
		x = a;
		y = b;
		z = c;
	}
};
struct glTMatrix
{
	double arr[3][3];
	glTMatrix() {}
	glTMatrix(double m11, double m12, double m13, double m21, double m22, double m23, double m31, double m32, double m33)
	{
		arr[0][0] = m11; arr[0][1] = m12; arr[0][2] = m13;
		arr[1][0] = m21; arr[1][1] = m22; arr[1][2] = m23;
		arr[2][0] = m31; arr[2][1] = m32; arr[2][2] = m33;
	}
};


//#################################################################


// Линейная алгебра в R3
void cross(glTVector a, glTVector b, glTVector &t);
double dot(glTVector a, glTVector b);
double len(glTVector a);
void normalize(glTVector &a);
void minus(glTVector a, glTVector b, glTVector &c);
double det(glTMatrix &a);
bool invers(glTMatrix a, glTMatrix &inv);
void mult(glTMatrix a, glTMatrix b, glTMatrix &s);
void mult(glTMatrix a, glTVector x, glTVector &s);
void mult(glTVector x, double t, glTVector &s);


//#################################################################

#endif