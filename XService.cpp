#pragma hdrstop

#include "XService.h"
#include <math.h>
#include <stdlib.h>

//#################################################################



//*************************************
void cross(glTVector a, glTVector b, glTVector &t)
{
	t.x = a.y*b.z - a.z*b.y;
	t.y = - (a.x*b.z - a.z*b.x);
	t.z = a.x*b.y - a.y*b.x;
}
//*************************************
double dot(glTVector a, glTVector b)
{
	return
		a.x*b.x+a.y*b.y+a.z*b.z;
}
//*************************************
double len(glTVector a)
{
	return
		sqrt(dot(a,a));
}
//*************************************
void normalize(glTVector &a)
{
	double l = len(a);
	a.x /= l;
	a.y /= l;
	a.z /= l;
}
//*************************************
void minus(glTVector a, glTVector b, glTVector &c)
{
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;
}
//*************************************
double det(glTMatrix &a)
{
	return
		a.arr[0][0]*( a.arr[1][1]*a.arr[2][2] - a.arr[1][2]*a.arr[2][1] ) -
		a.arr[0][1]*( a.arr[1][0]*a.arr[2][2] - a.arr[1][2]*a.arr[2][0] ) +
		a.arr[0][2]*( a.arr[1][0]*a.arr[2][1] - a.arr[1][1]*a.arr[2][0] );
}
//*************************************
bool invers(glTMatrix a, glTMatrix &inv)
{
	if( fabs(det(a))<1e-16 )
		return false;

	inv.arr[0][0] =   ( a.arr[1][1]*a.arr[2][2] - a.arr[1][2]*a.arr[2][1] ) / det(a);
	inv.arr[0][1] = - ( a.arr[0][1]*a.arr[2][2] - a.arr[0][2]*a.arr[2][1] ) / det(a);
	inv.arr[0][2] =   ( a.arr[0][1]*a.arr[1][2] - a.arr[0][2]*a.arr[1][1] ) / det(a);
	inv.arr[1][0] = - ( a.arr[1][0]*a.arr[2][2] - a.arr[1][2]*a.arr[2][0] ) / det(a);
	inv.arr[1][1] =   ( a.arr[0][0]*a.arr[2][2] - a.arr[0][2]*a.arr[2][0] ) / det(a);
	inv.arr[1][2] = - ( a.arr[0][0]*a.arr[1][2] - a.arr[0][2]*a.arr[1][0] ) / det(a);
	inv.arr[2][0] =   ( a.arr[1][0]*a.arr[2][1] - a.arr[1][1]*a.arr[2][0] ) / det(a);
	inv.arr[2][1] = - ( a.arr[0][0]*a.arr[2][1] - a.arr[0][1]*a.arr[2][0] ) / det(a);
	inv.arr[2][2] =   ( a.arr[0][0]*a.arr[1][1] - a.arr[0][1]*a.arr[1][0] ) / det(a);

	return true;
}
//*************************************
void mult(glTMatrix a, glTMatrix b, glTMatrix &s)
{
	for(int i=0; i<3; ++i)
		for(int j=0; j<3; ++j)
		{
			double c = 0;
			for(int k=0; k<3; ++k)
				c += a.arr[i][k] * b.arr[k][j];
			s.arr[i][j] = c;
		}
}
//*************************************
void mult(glTMatrix a, glTVector x, glTVector &s)
{
	s.x = a.arr[0][0]*x.x + a.arr[0][1]*x.y + a.arr[0][2]*x.z;
	s.y = a.arr[1][0]*x.x + a.arr[1][1]*x.y + a.arr[1][2]*x.z;
	s.z = a.arr[2][0]*x.x + a.arr[2][1]*x.y + a.arr[2][2]*x.z;
}
//*************************************
void mult(glTVector x, double t, glTVector &s)
{
	s.x = t*x.x;
	s.y = t*x.y;
	s.z = t*x.z;
}
//*************************************



//#################################################################