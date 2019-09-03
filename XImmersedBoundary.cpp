#include "XImmersedBoundary.h"
#include "XProblem.h"
#include "XUtils.h"
#include "stl_reader.h"
#include <cmath>
#include <functional>
#include <stdio.h>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <iostream>
using namespace std;

#define ISNAN(value) value != value
#define IN(begin, end, value) std::find(begin, end, value) != end
//#define M_PI acos(-1.)
struct xyz
{
    float x;
    float y;
    float z;

    xyz operator = (xyz b)
    {
        if (this == &b) {
            return *this;
        }
        x = b.x; y = b.y; z = b.z;
        return *this;
    }
    xyz& operator = (int b)
    {
        if (this->x == b && this->y == b && this->z == b) {
            return *this;
        }
        x = b; y = b; z = b;
        return *this;
    }
    //*
    bool operator==(   xyz& right)
    {
        if(x == right.x &&
            y == right.y &&
            z == right.z)
            return true;
            else return false;
    }
    //*/
};

xyz stepByCircle(xyz prevPoint, int delta)
{
    float betta  = delta * M_PI/180.0;
    xyz ll = prevPoint;
    ll.x = prevPoint.x*cos(betta) - prevPoint.y*sin(betta);
    ll.y = prevPoint.x*sin(betta) + prevPoint.y*cos(betta);
    return ll;
}

xyz stepByCircle(xyz prevPoint, xyz center, double delta)
{
    float betta  = delta * M_PI/180.0;
    xyz ll = prevPoint;
    ll.x = (prevPoint.x - center.x)*cos(betta) + (prevPoint.y-center.y)*sin(betta);
    ll.y = (prevPoint.x  -center.x)*sin(betta) - (prevPoint.y-center.y)*cos(betta);
    return ll;
}

bool isIntersect (xyz p1, xyz p2, xyz p3, xyz p4)
{
    //сначала расставим точки по порядку, т.е. чтобы было p1.x <= p2.x
    if (p2.x < p1.x) {

        xyz tmp = p1;
        p1 = p2;
        p2 = tmp;
    }
    //и p3.x <= p4.x
    if (p4.x < p3.x) {

        xyz tmp = p3;
        p3 = p4;
        p4 = tmp;
    }

    //проверим существование потенциального интервала для точки пересечения отрезков
    if (p2.x < p3.x) {
        return false; //ибо у отрезков нету взаимной абсциссы
    }

    //если оба отрезка вертикальные
    if((p1.x - p2.x == 0) && (p3.x - p4.x == 0)) {

        //если они лежат на одном X
        if(p1.x == p3.x) {

            //проверим пересекаются ли они, т.е. есть ли у них общий Y
            //для этого возьмём отрицание от случая, когда они НЕ пересекаются
            if (!((std::max(p1.y, p2.y) < std::min(p3.y, p4.y)) ||
                    (std::min(p1.y, p2.y) > std::max(p3.y, p4.y)))) {

                return true;
            }
        }

        return false;
    }

    //найдём коэффициенты уравнений, содержащих отрезки
    //f1(x) = A1*x + b1 = y
    //f2(x) = A2*x + b2 = y

    //если первый отрезок вертикальный
    if (p1.x - p2.x == 0) {

        //найдём Xa, Ya - точки пересечения двух прямых
        double Xa = p1.x;
        double A2 = (p3.y - p4.y) / (p3.x - p4.x);
        double b2 = p3.y - A2 * p3.x;
        double Ya = A2 * Xa + b2;

        if (p3.x <= Xa && p4.x >= Xa && std::min(p1.y, p2.y) <= Ya &&
                std::max(p1.y, p2.y) >= Ya) {

            return true;
        }

        return false;
    }

    //если второй отрезок вертикальный
    if (p3.x - p4.x == 0) {

        //найдём Xa, Ya - точки пересечения двух прямых
        double Xa = p3.x;
        double A1 = (p1.y - p2.y) / (p1.x - p2.x);
        double b1 = p1.y - A1 * p1.x;
        double Ya = A1 * Xa + b1;

        if (p1.x <= Xa && p2.x >= Xa && std::min(p3.y, p4.y) <= Ya &&
                std::max(p3.y, p4.y) >= Ya) {

            return true;
        }

        return false;
    }

    //оба отрезка невертикальные
    double A1 = (p1.y - p2.y) / (p1.x - p2.x);
    double A2 = (p3.y - p4.y) / (p3.x - p4.x);
    double b1 = p1.y - A1 * p1.x;
    double b2 = p3.y - A2 * p3.x;

    if (A1 == A2) {
        return false; //отрезки параллельны
    }

    //Xa - абсцисса точки пересечения двух прямых
    double Xa = (b2 - b1) / (A1 - A2);

    if ((Xa < std::max(p1.x, p3.x)) || (Xa > std::min( p2.x, p4.x))) {
        return false; //точка Xa находится вне пересечения проекций отрезков на ось X
    }
    else {
        return true;
    }
}

xyz findIntersection(xyz a, xyz b, xyz c, xyz d) //точки a и b концы первого отрезка  c и d второго
{
    xyz T;
    float x1 = a.x;
    float x2 = b.x;
    float x3 = c.x;
    float x4 = d.x;
    float y1 = a.y;
    float y2 = b.y;
    float y3 = c.y;
    float y4 = d.y;
    float avg_z = (a.z +b.z + c.z +d.z )/4;
    T.x = (((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)- (y1-y2)*(x3-x4)));
    T.y = (((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)- (y1-y2)*(x3-x4)));
    T.z=avg_z;
    return T;
}

void quickSort(xyz* arr, int left, int right, char flag_coord, char flag_dir)
{
    int i = left, j = right;
    xyz tmp;
    float pivot;
    if(flag_coord == 'z')
        pivot = arr[(left + right) / 2].z;
    else if(flag_coord == 'x')
        pivot = arr[(left + right) / 2].x;
    else if(flag_coord == 'y')
        pivot = arr[(left + right) / 2].y;
    else throw("quickSort:: error in direction!");
    /* partition */
    while (i <= j)
    {
        if(flag_coord == 'z')
        {
            if(flag_dir == 'u'){
                while (arr[i].z < pivot)
                    i++;
                while (arr[j].z > pivot)
                    j--;
                }else if(flag_dir == 'd'){
                    while (arr[i].z > pivot)
                        i++;
                    while (arr[j].z < pivot)
                        j--;
                }
            }
            else if(flag_coord == 'x')
            {
                if(flag_dir == 'u'){
                    while (arr[i].z < pivot)
                        i++;
                    while (arr[j].z > pivot)
                        j--;
                }else if(flag_dir == 'd'){
                    while (arr[i].z > pivot)
                        i++;
                    while (arr[j].z < pivot)
                        j--;
                }
            }
            else if(flag_coord == 'y')
            {
                if(flag_dir == 'u'){
                    while (arr[i].y < pivot)
                        i++;
                    while (arr[j].y > pivot)
                        j--;
                }else if(flag_dir == 'd'){
                    while (arr[i].y > pivot)
                        i++;
                    while (arr[j].y < pivot)
                        j--;
                }
            }else throw("quickSort:: error in direction!");

            if (i <= j) {
                tmp = arr[i];
                arr[i] = arr[j];
                arr[j] = tmp;
                i++;
                j--;
            }
        }
        if (left < j)   quickSort(arr, left, j,flag_coord,flag_dir);
        if (i < right)  quickSort(arr, i, right,flag_coord,flag_dir);


}

void quickSort(xyz*arr, float* angle_arr, int left, int right, char flag)
{
    int i = left, j = right;
    float angle_tmp;
    xyz tmp;
    float pivot = angle_arr[(left + right) / 2];
    /* partition */
    while (i <= j)
    {
        if(flag == 'u'){
            while (angle_arr[i] < pivot)
                i++;
            while (angle_arr[j] > pivot)
                j--;
        }else if(flag == 'd'){
            while (angle_arr[i] > pivot)
                i++;
            while (angle_arr[j] < pivot)
                j--;
        }else throw("quickSort(x,f,i,i,c):: error in direction!");
        if (i <= j) {
            angle_tmp = angle_arr[i];
            angle_arr[i] = angle_arr[j];
            angle_arr[j] = angle_tmp;
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
            i++;
            j--;
        }
    }
    if (left < j)   quickSort(arr, angle_arr, left, j,flag);
    if (i < right)  quickSort(arr, angle_arr, i, right,flag);
}

void findMassCenter(xyz* points, xyz& midpoint, int ibegin, int iend)
{
    float summ_x = 0;
    float summ_y = 0;
    float summ_z = 0;
    int count_point = iend - ibegin;
    //cout<<"fMCnum "<<count_point<<endl;
    for(auto i = ibegin;i<iend;i++)
    {
        summ_x+=points[i].x;
        summ_y+=points[i].y;
        summ_z+=points[i].z;
    }
    midpoint.x = summ_x/(1.0*count_point);
    midpoint.y = summ_y/(1.0*count_point);
    midpoint.z = summ_z/(1.0*count_point);
}

void findAngleInCircle(xyz* arr, xyz mp, float* angle_arr, int ibegin, int iend)
{
    for(auto i = ibegin; i<iend; i++)
        angle_arr[i]=atan((arr[i].y-mp.y)/(arr[i].x-mp.x));
}


void sortByRadius(xyz* arr, xyz midpoint, int ibegin, int iend)
{
    int num = iend-ibegin;
    //cout<<"num "<<num<<endl;
    xyz* x_plus  = new xyz [num];
    xyz* x_minus = new xyz [num];
    int count_x_plus = 0;
    int count_x_minus = 0;
    /*
    char fileName1[40];
    char fileName2[40];
    string_print(fileName1, 40, "%d-debug_x+.txt", glob_iter);
    string_print(fileName2, 40, "%d-debug_x-.txt", glob_iter);
    FILE *xp = fopen(fileName1,"w");
    FILE *xm = fopen(fileName2,"w");
    glob_iter++; //*/
    for(auto i = ibegin; i<iend; i++)
    {
        if(arr[i].x - midpoint.x>=0)
        {
            x_plus[count_x_plus] = arr[i];
            //fprintf(xp,"%f\n", x_plus[count_x_plus]);
            count_x_plus++;

        }else{
            x_minus[count_x_minus] = arr[i];
            //fprintf(xm,"%f\n", x_minus[count_x_minus]);
            count_x_minus++;
        }
    }
    float* angle_arr_x_plus = new float [count_x_plus];
    float* angle_arr_x_minus = new float [count_x_minus];

    findAngleInCircle(x_plus,midpoint,angle_arr_x_plus,0,count_x_plus);
    findAngleInCircle(x_minus,midpoint,angle_arr_x_minus,0,count_x_minus);
    quickSort(x_plus, angle_arr_x_plus, 0,count_x_plus-1,'d');
    quickSort(x_minus, angle_arr_x_minus, 0,count_x_minus-1,'d');

    for(auto i = ibegin; i<iend;i++)
    {
        //printf("ibeg > %4d | iend > %4d | cx+: %4d |  cx-: %4d | ib+cx+: %4d |   i: %4d | i-li-ib+cx+: %4d\n",ibegin,iend,count_x_plus,count_x_minus,ibegin+count_x_plus,i,i-count_x_plus-ibegin);
        if(i < ibegin+count_x_plus) arr[i] = x_plus[i-ibegin];
        else                        arr[i] = x_minus[i-count_x_plus-ibegin];
    }
    //fclose(xp);
    //fclose(xm);
    delete []x_plus;
    delete []x_minus;
    delete []angle_arr_x_plus;
    delete []angle_arr_x_minus;
}


TNode::TNode()
{
    x = y = z = 0;
    xPrev = yPrev = zPrev = -100000;
    xRef = yRef = zRef = 0;
    xForce = yForce = zForce = 0;
    xVel = yVel = zVel = 0;
    neighbors.next = NULL;
    neighbors.prev = NULL;
}

TNode::Getter TNode::Ox = {TNode::getX, TNode::getXRef, OX};
TNode::Getter TNode::Oy = {TNode::getY, TNode::getYRef, OY};
TNode::Getter TNode::Oz = {TNode::getZ, TNode::getZRef, OZ};

TImmersedBoundary::~TImmersedBoundary()
{
    //delete this->nodes;
}

TSphereBoundary* TSphereBoundary::Initialize()
{
    this->nodes = new TNode*[this->nodesCount];

    for (int theta = 0; theta < 10; theta++) {
        for (int phi = 0; phi < 10; phi++) {
            this->nodes[phi + theta*10] = new TNode();

            double x = this->xCenter + this->radius * cos(double(theta)/9.0 * 2 * M_PI) * sin(double(phi)/9.0 * M_PI);
            this->nodes[phi + theta*10]->x = x;
            this->nodes[phi + theta*10]->xRef = x;
            this->nodes[phi + theta*10]->xVel = 0;
            this->nodes[phi + theta*10]->xForce = 0;

            double y = this->yCenter + this->radius * sin(double(theta)/9.0 * 2 * M_PI) * sin(double(phi)/9.0 * M_PI);
            this->nodes[phi + theta*10]->y = y;
            this->nodes[phi + theta*10]->yRef = y;
            this->nodes[phi + theta*10]->yVel = 0;
            this->nodes[phi + theta*10]->yForce = 0;

            double z = this->zCenter + this->radius * cos(double(phi)/9.0 * M_PI);
            this->nodes[phi + theta*10]->z = z;
            this->nodes[phi + theta*10]->zRef = z;
            this->nodes[phi + theta*10]->zVel = 0;
            this->nodes[phi + theta*10]->zForce = 0;
        }
    }

    return this;
}

double TSphereBoundary::GetArea()
{
    return 4 * M_PI * this->radius * this->radius/ this->nodesCount;
}

TSimpleValve* TSimpleValve::Initialize()
{
    this->nodes = new TNode*[this->nodesCount];

    for(int n = 0; n < 10; ++n) {
        for (int i = 0; i < 10; i++) {
            this->nodes[i+n*10] = new TNode();

            double x = this->xCenter + this->radius * sin(-M_PI * 0.5 * (double)(9-i)/9);
            this->nodes[i+n*10]->x = x;
            this->nodes[i+n*10]->xRef = x;
            this->nodes[i+n*10]->xVel = 0;
            this->nodes[i+n*10]->xForce = 0;

            double y = 0.5 - this->radius * cos(-M_PI * 0.5 * (double)(9-i)/9);
            this->nodes[i+n*10]->y = y;
            this->nodes[i+n*10]->yRef = y;
            this->nodes[i+n*10]->yVel = 0;
            this->nodes[i+n*10]->yForce = 0;

            double z = 0.5 * n / 9;
            this->nodes[i+n*10]->z = z;
            this->nodes[i+n*10]->zRef = z;
            this->nodes[i+n*10]->zVel = 0;
            this->nodes[i+n*10]->zForce = 0;
        }
    }

    for(int n = 0; n < 10; ++n) {
        for (int i = 0; i < 10; i++) {
            double x = this->xCenter + this->radius * sin(-M_PI * 0.5 * (double)(9-i)/9);
            this->nodes[i+n*10+100]->x = x;
            this->nodes[i+n*10+100]->xRef = x;
            this->nodes[i+n*10+100]->xVel = 0;
            this->nodes[i+n*10+100]->xForce = 0;

            double y = 0.0 + this->radius * cos(-M_PI * 0.5 * (double)(9-i)/9);
            this->nodes[i+n*10+100]->y = y;
            this->nodes[i+n*10+100]->yRef = y;
            this->nodes[i+n*10+100]->yVel = 0;
            this->nodes[i+n*10+100]->yForce = 0;

            double z = 0.5 * n / 9 ;
            this->nodes[i+n*10+100]->z = z;
            this->nodes[i+n*10+100]->zRef = z;
            this->nodes[i+n*10+100]->zVel = 0;
            this->nodes[i+n*10+100]->zForce = 0;
        }
    }


    return this;
}

double TSimpleValve::GetArea()
{
    return 2 * M_PI * (this->radius / 4.0) * 0.5 / this->nodesCount;
}

double TCylinderNode::GetStiffnessX(int timestep)
{
    double coeff = 400;
    if (this->x < 0.2 || this->x > 0.8)
    {
        coeff = 400;
    }

    if(this->x > 0.4 && this->x < 0.6)
    {
        coeff = 20.0;
    }

    return (this->x - this->xRef) * coeff;
}

double TCylinderNode::GetStiffnessY(int timestep)
{
    double coeff = 400;
    if (this->x < 0.2 || this->x > 0.8)
    {
        coeff = 400;
    }

    TCylinderBoundary *boundary = (TCylinderBoundary *)this->boundary;

    if(
            this->x > 0.3 && this->x < 0.7 &&
            (this->x - 0.5) * (this->x - 0.5) + (this->z - 0.25) * (this->z - 0.25) < 0.4 &&
            this->y > boundary->yCenter && timestep > 1)
    {
        if (this->y > this->yRef)
        {
            coeff = 0.0000001;
        }
    }

    return (this->y - this->yRef) * coeff;
}

double TCylinderNode::GetStiffnessZ(int timestep)
{
    double coeff = 400;
    if (this->x < 0.2 || this->x > 0.8)
    {
        coeff = 400;
    }

    TCylinderBoundary *boundary = (TCylinderBoundary *)this->boundary;

    if(
            this->x > 0.3 && this->x < 0.7 &&
            (this->x - 0.5) * (this->x - 0.5) + (this->z - 0.25) * (this->z - 0.25) < 0.4 &&
            this->y > boundary->yCenter && timestep > 1)
    {
        if (this->y > this->yRef)
        {
            coeff = 0.000001;
        }
    }

    return (this->z - this->zRef) * coeff;
}

TCylinderBoundary* TCylinderBoundary::Initialize()
{
    this->nodes = new TNode*[this->nodesCount];

    for(int n = 0; n < 120; ++n) {
        for (int i = 0; i < 120; i++) {
            this->nodes[i+n*120] = new TNode();
            this->nodes[i+n*120]->boundary = this;

            if (n < 119) this->nodes[i+n*120]->neighbors.next = this->nodes[i+(n+1)*120];
            if (n > 0) this->nodes[i+n*120]->neighbors.prev = this->nodes[i+(n-1)*120];

            double x = float(n) / float(120);
            this->nodes[i+n*120]->x = x;
            this->nodes[i+n*120]->xRef = x;
            this->nodes[i+n*120]->xVel = 0;
            this->nodes[i+n*120]->xForce = 0;

            double y = this->yCenter - this->radius * cos(-M_PI * 2 * (double)(119-i)/119);
            this->nodes[i+n*120]->y = y;
            this->nodes[i+n*120]->yRef = y;
            this->nodes[i+n*120]->yVel = 0;
            this->nodes[i+n*120]->yForce = 0;

            double z = this->zCenter + this->radius * sin(-M_PI * 2 * (double)(119-i)/119);
            this->nodes[i+n*120]->z = z;
            this->nodes[i+n*120]->zRef = z;
            this->nodes[i+n*120]->zVel = 0;
            this->nodes[i+n*120]->zForce = 0;
        }
    }

    return this;
}

double TCylinderBoundary::GetArea()
{
    return 2 * M_PI * this->radius * 1 / this->nodesCount;
}

TDeformedCylinderBoundary* TDeformedCylinderBoundary::Initialize()
{
    this->nodes = new TNode*[this->nodesCount];
    double radius = this->radius;

    for(int n = 0; n < 120; ++n) {
        for (int i = 0; i < 120; i++) {
            this->nodes[i+n*120] = new TCylinderNode();
            this->nodes[i+n*120]->boundary = this;

            double x = float(n) / float(120);
            this->nodes[i+n*120]->x = x;
            this->nodes[i+n*120]->xRef = x;
            this->nodes[i+n*120]->xVel = 0;
            this->nodes[i+n*120]->xForce = 0;

            if (x > 0.4 && x < 0.6)
            {
                radius = (this->radius - deformation * (1 - fabs(0.5-x)/0.1) );
            }

            double y = this->yCenter - radius * cos(-M_PI * 2 * (double)(119-i)/119);
            this->nodes[i+n*120]->y = y;
            this->nodes[i+n*120]->yRef = y;
            this->nodes[i+n*120]->yVel = 0;
            this->nodes[i+n*120]->yForce = 0;

            double z = this->zCenter + radius * sin(-M_PI * 2 * (double)(119-i)/119);
            this->nodes[i+n*120]->z = z;
            this->nodes[i+n*120]->zRef = z;
            this->nodes[i+n*120]->zVel = 0;
            this->nodes[i+n*120]->zForce = 0;
        }
    }

    return this;
}

double TDeformedCylinderBoundary::GetArea()
{
    return 2 * M_PI * this->radius * 1 / this->nodesCount;
}

double TTargetNode::GetStretchingForce(TNode *neighbor, TNode::Getter getters)
{
    return 0;
}

double TTargetNode::GetBendingForce(TNode *next, TNode *prev, TNode::Getter getters)
{
    return 0;
}

double TTargetNode::GetTargetForce(TNode::Getter getters)
{
    return this->boundary->stiffness * (getters.coord(this) - getters.ref(this));
}

double TElasticNode::GetStretchingForce(TNode *neighbor, TNode::Getter getters)
{
    if (neighbor == NULL)
    {
        return 0;
    }

    double distance = TNode::distance_3d(this, neighbor);
    double distance_axis = getters.coord(this) - getters.coord(neighbor);
    double initial_distance = this->neighbors.initialDistanceNext;

    return this->boundary->stretchingStiffness * (distance - initial_distance) * distance_axis / distance;
}

double TElasticNode::GetBendingForce(TNode *next, TNode *prev, TNode::Getter getters)
{
    if (next == NULL || prev == NULL)
    {
        return 0;
    }

    double initial_curvature = this->neighbors.initialCurvatures[getters.type];

    return this->boundary->bendingStiffness * (initial_curvature - TNode::curvature(this, next, prev, getters));
}

double TElasticNode::GetTargetForce(TNode::Getter getters)
{
    return 0;
}

double TFixedElasticNode::GetStretchingForce(TNode *neighbor, TNode::Getter getters)
{
    if (neighbor == NULL)
    {
        return 0;
    }

    double distance = TNode::distance_3d(this, neighbor);
    double distance_axis = getters.coord(this) - getters.coord(neighbor);
    double initial_distance = this->neighbors.initialDistanceNext;

    return this->stretchingStiffness * (distance - initial_distance) * distance_axis / distance;
}

double TFixedElasticNode::GetBendingForce(TNode *next, TNode *prev, TNode::Getter getters)
{
    return 0;

    //if (next == NULL || prev == NULL)
    //{
        //return 0;
    //}

    //double initial_curvature = this->neighbors.initialCurvatures[getters.type];

    //return this->boundary->bendingStiffness * (initial_distance - TNode::curvature(this, next, prev, getters));
}

double TFixedElasticNode::GetTargetForce(TNode::Getter getters)
{
     return this->stretchingStiffness * (getters.coord(this) - getters.ref(this));
}



double TFixedStretchingNode::GetStretchingForce(TNode *neighbor, TNode::Getter getters)
{
    if (neighbor == NULL)
    {
        return 0;
    }

    double distance = TNode::distance_3d(this, neighbor);
    double distance_axis = getters.coord(this) - getters.coord(neighbor);
    double initial_distance = this->neighbors.initialDistanceNext;

    return this->boundary->stretchingStiffness * (distance - initial_distance) * distance_axis / distance;
}

double TFixedStretchingNode::GetBendingForce(TNode *next, TNode *prev, TNode::Getter getters)
{
    return 0;
}

double TFixedStretchingNode::GetTargetForce(TNode::Getter getters)
{
    return this->boundary->stiffness * (getters.coord(this) - getters.ref(this));
}



TCylinderElasticBoundary* TCylinderElasticBoundary::Initialize()
{
	//------------- NEW ---------------------------------------------------------------------
	int tube_type = 8;
	//bool make_shift = true;
	double shift_angle = 0.;
	double shift_step = shift_angle*M_PI / 180.;
	double shift = 0.0;

	this->yCenter = this->GetYCenter();// - 0.022;
	this->zCenter = this->GetZCenter();// - 0.022;
	
	if (tube_type == 1) //line chanal
	{
		this->step = 0.01; // 0.069;
		double iPart = 0.;
		//this->height = 1.0;

		modf(this->height / this->step, &iPart);
		int endIndex = (int)iPart;
		this->height_nodes = endIndex + 1;

		modf(this->radius*2.*M_PI / this->step, &iPart);
		endIndex = (int)iPart;
		this->radius_nodes = endIndex;
		this->nodesCount = this->height_nodes*this->radius_nodes;
		this->nodes = new TNode*[this->nodesCount];

	//*
		double begin_center = 0.5;

		//printf(" >> this->yCenter = %lf, this->zCenter = %lf\n",this->yCenter, this->zCenter);
	//*/
		/*double begin_center = 0.319;
		this->yCenter = 0.319;
		this->zCenter = 0.319;
		*/
		for (int n = 0; n < this->height_nodes; ++n) {
			for (int i = 0; i < this->radius_nodes; i++) {
				this->nodes[i + n*this->radius_nodes] = new TFixedElasticNode();
				//this->nodes[i + n * 120]->;
			}
		}


		double rc, stiffness_plus = this->stiffness / 30;
		double iRadius = (double) this->radius_nodes;
		double iHeight = (double) this->height_nodes - 1.;
		int cx = this->height_nodes / 2., cy = this->radius_nodes / 4., Rpip = 0.00; // this->radius_nodes / 10.,
		int ix, iy;

		for (int n = 0; n < this->height_nodes; ++n) {
			for (int i = 0; i < this->radius_nodes; i++) {
				double x = this->height*float(n) / iHeight;
				double y = this->yCenter + this->radius * cos(-M_PI * 2 * (iRadius - i) / iRadius); // +shift_step);
				double z = this->zCenter + this->radius * sin(-M_PI * 2 * (iRadius - i) / iRadius); // +shift_step);

				this->nodes[i + n*this->radius_nodes]->x = x;
				this->nodes[i + n*this->radius_nodes]->xRef = x;
				this->nodes[i + n*this->radius_nodes]->xVel = 0;
				this->nodes[i + n*this->radius_nodes]->xForce = 0;

				this->nodes[i + n*this->radius_nodes]->y = y;
				this->nodes[i + n*this->radius_nodes]->yRef = y;
				this->nodes[i + n*this->radius_nodes]->yVel = 0;
				this->nodes[i + n*this->radius_nodes]->yForce = 0;

				this->nodes[i + n*this->radius_nodes]->z = z;
				this->nodes[i + n*this->radius_nodes]->zRef = z;
				this->nodes[i + n*this->radius_nodes]->zVel = 0;
				this->nodes[i + n*this->radius_nodes]->zForce = 0;

				/*
				ix = abs(n - cx);
				iy = abs(i - cy);
				rc = sqrt((double)(ix*ix + iy*iy));
				if (rc <= Rpip)
					this->nodes[i + n * this->radius_nodes]->stretchingStiffness = stiffness_plus; //+(this->stiffness - stiffness_plus)*sqrt(rc);    //sin(rc*0.5*M_PI / Rpip);
				else
				*/
					this->nodes[i + n *this->radius_nodes]->stretchingStiffness = this->stiffness;


			}
		}
	}
	else if (tube_type == 2) // dyga
	{
		char filename[] = "line.msh";
		FILE* SplFile;

		SplFile = fopen(filename, "r"); // do exception !!!
		if (SplFile == NULL)
		{
			printf("No line.msh file.");
		}
		else
		{
			int nSpl;
			int tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;

			fscanf(SplFile, "%d %d %d %d %d", &tmp1, &tmp2, &tmp3, &nSpl, &tmp4);
			double *xSpl;
			xSpl = new double[nSpl];
			double *zSpl;
			zSpl = new double[nSpl];
			double dtmp1;

			for (int i = 0; i < nSpl; i++)
			{
				fscanf(SplFile, "%d %d %d %d %d %d %d %d", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6, &tmp7, &tmp8);
				fscanf(SplFile, "%lf %lf %lf", &xSpl[i], &zSpl[i], &dtmp1);
			}
			fclose(SplFile);

			//this->step = sqrtl(pow(xSpl[1] - xSpl[0], 2) + pow(zSpl[1] - zSpl[0], 2));
			//printf("\n==========================\n");
			//printf("this->step", this->step);
			double iPart = 0.;
			this->step = 0.01;
			this->height = 1.0;
			modf(this->height / this->step, &iPart);
			int endIndex = (int)iPart;
			this->height_nodes = endIndex + 1;


			//this->height_nodes = nSpl;
			modf(this->radius*2.*M_PI / this->step, &iPart);
			endIndex = (int)iPart;
			this->radius_nodes = endIndex;
			this->nodesCount = this->height_nodes*this->radius_nodes;
			this->nodes = new TNode*[this->nodesCount];
			this->yCenter = 0.4;
			//this->zCenter = 0.0;

			for (int n = 0; n < this->height_nodes; ++n) {
				for (int i = 0; i < this->radius_nodes; i++) {
					this->nodes[i + n*this->radius_nodes] = new TFixedElasticNode();
					//this->nodes[i + n * 120]->;
				}
			}

			double rc, stiffness_plus = this->stiffness / 100;
			double iRadius = (double) this->radius_nodes;
			double iHeight = (double) this->height_nodes - 1.;
			int cx = this->height_nodes / 2., cy = 3.*this->radius_nodes / 4., Rpip = this->radius_nodes / 15., ix, iy;
			int iSpl = 0;
			double Alpha, z_plus;

			for (int n = 0; n < this->height_nodes; ++n) {
				double x;
				if (n == 0)
				{
					x = xSpl[0];
					z_plus = zSpl[0];
					printf("-->z_b=%lf   ", z_plus);
				}
				else if (n == this->height_nodes - 1)
				{
					x = xSpl[nSpl - 1];
					z_plus = zSpl[nSpl - 1];
					printf("-->z_e=%lf   ", z_plus);
				}
				else // Finding intersection point
					//--------------------------------
				{
					x = this->height*float(n) / iHeight;
					while (!(xSpl[iSpl] <= x && x < xSpl[iSpl + 1]))
						iSpl++;
					Alpha = (x - xSpl[iSpl + 1]) / (xSpl[iSpl] - xSpl[iSpl + 1]);
					z_plus = zSpl[iSpl] * Alpha + (1. - Alpha)*zSpl[iSpl + 1];
					if (iSpl != 0) iSpl--;
				}
				for (int i = 0; i < this->radius_nodes; i++) {


					double y = this->yCenter + this->radius * cos(-M_PI * 2 * (iRadius - i) / iRadius + shift_step);
					double z = z_plus + this->radius * sin(-M_PI * 2 * (iRadius - i) / iRadius + shift_step);

					this->nodes[i + n*this->radius_nodes]->x = x;
					this->nodes[i + n*this->radius_nodes]->xRef = x;
					this->nodes[i + n*this->radius_nodes]->xVel = 0;
					this->nodes[i + n*this->radius_nodes]->xForce = 0;

					this->nodes[i + n*this->radius_nodes]->y = y;
					this->nodes[i + n*this->radius_nodes]->yRef = y;
					this->nodes[i + n*this->radius_nodes]->yVel = 0;
					this->nodes[i + n*this->radius_nodes]->yForce = 0;

					this->nodes[i + n*this->radius_nodes]->z = z;
					this->nodes[i + n*this->radius_nodes]->zRef = z;
					this->nodes[i + n*this->radius_nodes]->zVel = 0;
					this->nodes[i + n*this->radius_nodes]->zForce = 0;

					/*
					ix = abs(n - cx);
					iy = abs(i - cy);
					rc = sqrt((double)(ix*ix + iy*iy));
					if (rc <= Rpip)
					//this->nodes[i + n * 120]->stretchingStiffness = this->stiffness*rc*rc + sin(rc/Rpip)*stiffness_plus;
					this->nodes[i + n * this->radius_nodes]->stretchingStiffness = stiffness_plus; //+(this->stiffness - stiffness_plus)*sqrt(rc);    //sin(rc*0.5*M_PI / Rpip);
					else*/
					this->nodes[i + n *this->radius_nodes]->stretchingStiffness = this->stiffness;

					// Verison 1 (not good)
					/*
					if ((n>=49 && n<=69) && (i>=20 && i<=39))
					this->nodes[i + n * 120]->stretchingStiffness = this->stiffness/100;
					else
					this->nodes[i + n * 120]->stretchingStiffness=this->stiffness;
					*/
				}
				//this->zCenter += 0.005;
			}
			//this->zCenter -= 0.005;
			//printf("\nzCenter=%lf\n", this->zCenter);
		}
		delete[]xSpl;// xSpl = NULL;
		delete[]zSpl;// zSpl = NULL;
	}
	else if (tube_type == 3) //skrucheno
	{
		printf("\nBOUNDARY TYPE = %1d step=%lf\n", tube_type, shift_step);
		this->step = 0.01; // 0.069;
		double iPart = 0.;
		//this->height = 1.0;

		modf(this->height / this->step, &iPart);
		int endIndex = (int)iPart;
		this->height_nodes = endIndex + 1;

		modf(this->radius*2.*M_PI / this->step, &iPart);
		endIndex = (int)iPart;
		this->radius_nodes = endIndex;
		this->nodesCount = this->height_nodes*this->radius_nodes;
		this->nodes = new TNode*[this->nodesCount];

		//double begin_center = 0.4;
		this->yCenter = 0.522;
		this->zCenter = 0.524;
		double an_begin = 0.25;
		double an1 = 0.35, an2 = 0.65;
		double an_end = 0.75;
		double step_div = (an1 - an_begin) / this->step;
		double Alpha_step = 1. / step_div;
		double Alpha = 1.;
		bool Alpha_flag = true;

		for (int n = 0; n < this->height_nodes; ++n)
		{
			for (int i = 0; i < this->radius_nodes; i++)
			{
				this->nodes[i + n*this->radius_nodes] = new TFixedElasticNode();
				//this->nodes[i + n * 120]->;
			}
		}


		double rc, stiffness_plus = this->stiffness / 30;
		double iRadius = (double) this->radius_nodes;
		double iHeight = (double) this->height_nodes - 1.;
		//int cx = this->height_nodes / 2., cy = this->radius_nodes / 4., Rpip = this->radius_nodes / 10., ix, iy;

		for (int n = 0; n < this->height_nodes; ++n)
		{
			if (n != 0) shift += shift_step;
			printf("\nn=%d shift=%lf ", n, shift);

			for (int i = 0; i < this->radius_nodes; i++)
			{
				double x = this->height*float(n) / iHeight;
				double y = this->yCenter + this->radius * cos(-M_PI * 2 * (iRadius - i) / iRadius + shift);
				double z = this->zCenter + this->radius * sin(-M_PI * 2 * (iRadius - i) / iRadius + shift);
				if (i == 0)
				{
					printf("  --> %lf %lf", y, z);
					if ((n + 1) % 2 == 0) printf("\n");
				}

				this->nodes[i + n*this->radius_nodes]->x = x;
				this->nodes[i + n*this->radius_nodes]->xRef = x;
				this->nodes[i + n*this->radius_nodes]->xVel = 0;
				this->nodes[i + n*this->radius_nodes]->xForce = 0;

				this->nodes[i + n*this->radius_nodes]->y = y;
				this->nodes[i + n*this->radius_nodes]->yRef = y;
				this->nodes[i + n*this->radius_nodes]->yVel = 0;
				this->nodes[i + n*this->radius_nodes]->yForce = 0;

				this->nodes[i + n*this->radius_nodes]->z = z;
				this->nodes[i + n*this->radius_nodes]->zRef = z;
				this->nodes[i + n*this->radius_nodes]->zVel = 0;
				this->nodes[i + n*this->radius_nodes]->zForce = 0;


				//ix = abs(n - cx);
				//iy = abs(i - cy);
				//rc = sqrt((double)(ix*ix + iy*iy));
				if (x >= an_begin && x <= an1)
				{
					this->nodes[i + n * this->radius_nodes]->stretchingStiffness = Alpha*this->stiffness + (1. - Alpha)*stiffness_plus; //+(this->stiffness - stiffness_plus)*sqrt(rc);    //sin(rc*0.5*M_PI / Rpip);
					if (i == this->radius_nodes - 1)
					{
						Alpha -= Alpha_step;
						printf("-->%lf  ", this->nodes[i + n * this->radius_nodes]->stretchingStiffness);
					}
				}
				else if (x > an1 && x < an2)
				{
					this->nodes[i + n * this->radius_nodes]->stretchingStiffness = stiffness_plus;
					if (Alpha_flag && i == this->radius_nodes - 1)
					{
						Alpha_flag = false;
						Alpha = 0.;
						printf("-->%lf  ", this->nodes[i + n * this->radius_nodes]->stretchingStiffness);
					}
				}
				else if (x >= an2 && x <= an_end)
				{
					this->nodes[i + n * this->radius_nodes]->stretchingStiffness = Alpha*this->stiffness + (1. - Alpha)*stiffness_plus; //+(this->stiffness - stiffness_plus)*sqrt(rc);    //sin(rc*0.5*M_PI / Rpip);
					if (i == this->radius_nodes - 1)
					{
						Alpha += Alpha_step;
						printf("-->%lf  ", this->nodes[i + n * this->radius_nodes]->stretchingStiffness);
					}
				}
				else
					this->nodes[i + n *this->radius_nodes]->stretchingStiffness = this->stiffness;


			}
		}

		delete[]xSpl;// xSpl = NULL;
		delete[]ySpl;// ySpl = NULL;
		delete[]zSpl;// zSpl = NULL;
	}
	else if (tube_type == 4) // dyga
	{
		char filename[] = "line.msh";
		FILE* SplFile;

		SplFile=fopen(filename, "r"); // do exception !!!
		if (SplFile == NULL)
		{
			printf("No line.msh file.");
		}
		else
		{
			//int nSpl;
			int tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;

			fscanf(SplFile, "%d %d %d %d %d", &tmp1, &tmp2, &tmp3, &nSpl, &tmp4);
			//double *xSpl;
			xSpl = new double[nSpl];
			//double *ySpl;
			ySpl = new double[nSpl];
			//double *zSpl;
			zSpl = new double[nSpl];
			double dtmp1;

			for (int i =0 ; i <=nSpl-1 ; i++)
			{
				fscanf(SplFile, "%d %d %d %d %d %d %d %d", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6, &tmp7, &tmp8);
				fscanf(SplFile, "%lf %lf %lf", &xSpl[i], &zSpl[i], &dtmp1);
				//zSpl[i] += 1.0;
				ySpl[i] = 0.4;
				printf("\n %2d  %lf %lf %lf", i, xSpl[i], ySpl[i], zSpl[i]);
			}
			fclose(SplFile);

			//this->step = sqrtl(pow(xSpl[1] - xSpl[0], 2) + pow(zSpl[1] - zSpl[0], 2));
			//printf("\n==========================\n");
			//printf("this->step", this->step);
			double iPart = 0.;
			this->step = 0.01;
			this->height = 1.0;
			modf(this->height / this->step, &iPart);
			int endIndex = (int)iPart;
			this->height_nodes = endIndex + 1;


			//this->height_nodes = nSpl;
			modf(this->radius*2.*M_PI / this->step, &iPart);
			endIndex = (int)iPart;
			this->radius_nodes = endIndex;
			this->nodesCount = this->height_nodes*this->radius_nodes;
			this->nodes = new TNode*[this->nodesCount];
			//this->yCenter = 0.4;
			//this->zCenter = 0.0;

			for (int n = 0; n < this->height_nodes; ++n) {
				for (int i = 0; i < this->radius_nodes; i++) {
					this->nodes[i + n*this->radius_nodes] = new TFixedElasticNode();
					//this->nodes[i + n * 120]->;
				}
			}

			double rc, stiffness_plus = this->stiffness / 28;
			double iRadius = (double) this->radius_nodes;
			double iHeight = (double) this->height_nodes - 1.;
			double dcx =0.5;  //this->height_nodes / 2.;
			double dcy = 0.4;  //3.*this->radius_nodes / 4.;
			double dcz = 0.35;
			double dRmax = 0.18, dRmin = 0.9;
			double dRpip = dRmax;  //this->radius_nodes / 15.;

			int ix, iy;
			int iSpl = 0;
			double Alpha, y_plus, z_plus;

			for (int n = 0; n < this->height_nodes; ++n) {
				double x;
				if (n == 0)
				{
					x = xSpl[0];
					y_plus = ySpl[0];
					z_plus = zSpl[0];
					printf("-->z_b=%lf   ", z_plus);
				}
				else if (n == this->height_nodes - 1)
				{
					x = xSpl[nSpl - 1];
					y_plus = ySpl[nSpl - 1];
					z_plus = zSpl[nSpl - 1];
					printf("-->z_e=%lf   ", z_plus);
				}
				else // Finding intersection point
					//--------------------------------
				{
					x = this->height*float(n) / iHeight;
					while (!(xSpl[iSpl] <= x && x < xSpl[iSpl + 1]))
						iSpl++;
					Alpha = (x - xSpl[iSpl + 1]) / (xSpl[iSpl] - xSpl[iSpl + 1]);
 					y_plus = ySpl[iSpl] * Alpha + (1. - Alpha)*ySpl[iSpl + 1];
					z_plus = zSpl[iSpl] * Alpha + (1. - Alpha)*zSpl[iSpl + 1];
					if (iSpl != 0) iSpl--;
				}
				for (int i = 0; i < this->radius_nodes; i++) {


					double y = y_plus + this->radius * cos(-M_PI * 2 * (iRadius - i) / iRadius + shift_step);
					double z = z_plus + this->radius * sin(-M_PI * 2 * (iRadius - i) / iRadius + shift_step);

					this->nodes[i + n*this->radius_nodes]->x = x;
					this->nodes[i + n*this->radius_nodes]->xRef = x;
					this->nodes[i + n*this->radius_nodes]->xVel = 0;
					this->nodes[i + n*this->radius_nodes]->xForce = 0;

					this->nodes[i + n*this->radius_nodes]->y = y;
					this->nodes[i + n*this->radius_nodes]->yRef = y;
					this->nodes[i + n*this->radius_nodes]->yVel = 0;
					this->nodes[i + n*this->radius_nodes]->yForce = 0;

					this->nodes[i + n*this->radius_nodes]->z = z;
					this->nodes[i + n*this->radius_nodes]->zRef = z;
					this->nodes[i + n*this->radius_nodes]->zVel = 0;
					this->nodes[i + n*this->radius_nodes]->zForce = 0;


					//ix = x - dcx;
					//iy = abs(i - cy);
					rc = sqrt((x - dcx)*(x - dcx) + (y - dcy)*(y - dcy) + (z - dcz)*(z - dcz));

					if (rc <= dRmax)
					{
						//Count Alpha
						double Alpha_rad = (rc - dRmax) / (dRmin - dRmax);
						if (0.0 <= Alpha_rad && Alpha_rad <= 1.0)
						{
							this->nodes[i + n * this->radius_nodes]->stretchingStiffness = (1.0 - Alpha_rad)*this->stiffness + Alpha_rad*stiffness_plus;
						}
						else
							this->nodes[i + n * this->radius_nodes]->stretchingStiffness = stiffness_plus;

					}
					else
					this->nodes[i + n *this->radius_nodes]->stretchingStiffness = this->stiffness;
					if (i == 0) printf("s(%2d)=%lf  ", n, this->nodes[i + n *this->radius_nodes]->stretchingStiffness);

					// Verison 1 (not good)
					/*
					if ((n>=49 && n<=69) && (i>=20 && i<=39))
					this->nodes[i + n * 120]->stretchingStiffness = this->stiffness/100;
					else
					this->nodes[i + n * 120]->stretchingStiffness=this->stiffness;
					*/
				}
				//this->zCenter += 0.005;
			}
			//this->zCenter -= 0.005;
			//printf("\nzCenter=%lf\n", this->zCenter);
		}
		delete[]xSpl;// xSpl = NULL;
		delete[]ySpl;// ySpl = NULL;
		delete[]zSpl;// zSpl = NULL;
	}
	else if (tube_type == 5) //linear stiff change
	{
		this->step = 0.01; // 0.069;
		double iPart = 0.;
		//this->height = 1.0;

		modf(this->height / this->step, &iPart);
		int endIndex = (int)iPart;
		this->height_nodes = endIndex + 1;

		modf(this->radius*2.*M_PI / this->step, &iPart);
		endIndex = (int)iPart;
		this->radius_nodes = endIndex;
		this->nodesCount = this->height_nodes*this->radius_nodes;
		this->nodes = new TNode*[this->nodesCount];

		//*
		double begin_center = 0.5;
		//this->yCenter = 0.522;
		//this->zCenter = 0.524;
		//printf(" >> this->yCenter = %lf, this->zCenter = %lf\n", this->yCenter, this->zCenter);
		//*/
		/*double begin_center = 0.319;
		this->yCenter = 0.319;
		this->zCenter = 0.319;
		*/
		for (int n = 0; n < this->height_nodes; ++n) {
			for (int i = 0; i < this->radius_nodes; i++) {
				this->nodes[i + n*this->radius_nodes] = new TFixedElasticNode();
				//this->nodes[i + n * 120]->;
			}
		}


		double rc=0, stiffness_plus = this->stiffness / 30;
		double iRadius = (double) this->radius_nodes;  //InCircle = fBoundary->radius_nodes;
		double iHeight = (double) this->height_nodes - 1.; //Circles = fBoundary->height_nodes;
		int cx = this->height_nodes / 2., cy = this->radius_nodes / 4., Rpip = 0.00; // this->radius_nodes / 10.,
		int ix, iy;
		std::ofstream hout;
		hout.open("toexel.txt", std::ios::out);
		for (int n = 0; n < this->height_nodes; ++n) {
			for (int i = 0; i < this->radius_nodes; i++) {
				double x = this->height*float(n) / iHeight;
				double y = this->yCenter + this->radius * cos(-M_PI * 2 * (iRadius - i) / iRadius); // +shift_step);
				double z = this->zCenter + this->radius * sin(-M_PI * 2 * (iRadius - i) / iRadius); // +shift_step);

				this->nodes[i + n*this->radius_nodes]->x = x;
				this->nodes[i + n*this->radius_nodes]->xRef = x;
				this->nodes[i + n*this->radius_nodes]->xVel = 0;
				this->nodes[i + n*this->radius_nodes]->xForce = 0;

				this->nodes[i + n*this->radius_nodes]->y = y;
				this->nodes[i + n*this->radius_nodes]->yRef = y;
				this->nodes[i + n*this->radius_nodes]->yVel = 0;
				this->nodes[i + n*this->radius_nodes]->yForce = 0;

				this->nodes[i + n*this->radius_nodes]->z = z;
				this->nodes[i + n*this->radius_nodes]->zRef = z;
				this->nodes[i + n*this->radius_nodes]->zVel = 0;
				this->nodes[i + n*this->radius_nodes]->zForce = 0;

				/*
				ix = abs(n - cx);
				iy = abs(i - cy);
				rc = sqrt((double)(ix*ix + iy*iy));
				if (rc <= Rpip)
				this->nodes[i + n * this->radius_nodes]->stretchingStiffness = stiffness_plus; //+(this->stiffness - stiffness_plus)*sqrt(rc);    //sin(rc*0.5*M_PI / Rpip);
				else
				*/

				double vx1 = 1.*n / 3;
				double vx2 = 2.*n / 3;
				double x1 = -1;
				double x2 = 1;
				double left = 1.*this->height_nodes / 3.;
				double right = 2.*this->height_nodes / 3.;
				double c2 = this->stiffness;
				double c1 = 0.4*this->stiffness;
				double c12 = (c1 + c2) / 2;
				double k = (x1*vx1 - x2*vx1 - 2 * x1*vx2) / (vx1*(vx1 - vx2)); //kx+b
				double b = (x2*vx1 - x1*vx2) / (vx1 - vx2);
				double stepen = 1.0 / 3.0;
				//printf("left > %d right > %d height > %d\n",left, right, this->height);
				if (n < left)
					this->nodes[i + n *this->radius_nodes]->stretchingStiffness = c1;
				else if ((n >= left) && (n <= right))
				{
					//x[i] = pow((gr[i])- (left+right)/2,stepen)+c1;
					//cout<<" under > "<<(gr[i])- (left+right)/2<<" \tpow > "<<x[i]<<endl;
					this->nodes[i + n *this->radius_nodes]->stretchingStiffness = (c2 - c1)*((n - left) / (right - left)) + c1;
					//this->nodes[(int)left + n *this->radius_nodes]->stretchingStiffness = (this->nodes[(int)left - 1 + n *this->radius_nodes]->stretchingStiffness + this->nodes[(int)left + 1 + n *this->radius_nodes]->stretchingStiffness ) / 2;
					//this->nodes[(int)right + n *this->radius_nodes]->stretchingStiffness = (this->nodes[(int)right - 1 + n *this->radius_nodes]->stretchingStiffness + this->nodes[(int)right + 1 + n *this->radius_nodes]->stretchingStiffness) / 2;

					//cout<<x[i]<<endl;
				}
				else
					this->nodes[i + n *this->radius_nodes]->stretchingStiffness = c2;
				//left = 0;
				/*
				if((n >= left) && (n <= right))

					//printf("im in\n");
					this->nodes[i + n *this->radius_nodes]->stretchingStiffness = this->stiffness - (this->stiffness - 100)*sin(this->height*float(n) / iHeight * M_PI);
					//this->nodes[i + n *this->radius_nodes]->stretchingStiffness = this->stiffness + (n - left)*(n - right);
				else
					this->nodes[i + n *this->radius_nodes]->stretchingStiffness = this->stiffness;
				//*/

				hout << this->nodes[i + n *this->radius_nodes]->stretchingStiffness<<"\n";

			}
		}
		hout.close();
	}
	else if (tube_type == 6) //circled stiff on border
	{
		this->step = 0.01; // 0.069;
		double iPart = 0.;
		//this->height = 1.0;

		modf(this->height / this->step, &iPart);
		int endIndex = (int)iPart;
		this->height_nodes = endIndex + 1;

		modf(this->radius*2.*M_PI / this->step, &iPart);
		endIndex = (int)iPart;
		this->radius_nodes = endIndex;
		this->nodesCount = this->height_nodes*this->radius_nodes;
		this->nodes = new TNode*[this->nodesCount];

		//*
		double begin_center = 0.5;
		//this->yCenter = 0.522;
		//this->zCenter = 0.524;
		//printf(" >> this->yCenter = %lf, this->zCenter = %lf\n", this->yCenter, this->zCenter);
		//*/
		/*double begin_center = 0.319;
		this->yCenter = 0.319;
		this->zCenter = 0.319;
		*/
		for (int n = 0; n < this->height_nodes; ++n) {
			for (int i = 0; i < this->radius_nodes; i++) {
				this->nodes[i + n*this->radius_nodes] = new TFixedElasticNode();
				//this->nodes[i + n * 120]->;
			}
		}


		double rc = 0, stiffness_plus = this->stiffness / 30;
		double iRadius = (double) this->radius_nodes;  //InCircle = fBoundary->radius_nodes;
		double iHeight = (double) this->height_nodes - 1.; //Circles = fBoundary->height_nodes;
		int cx = this->height_nodes / 2., cy = this->radius_nodes / 4., Rpip = 0.00; // this->radius_nodes / 10.,
		int ix, iy;
		for (int n = 0; n < this->height_nodes; ++n) {
			for (int i = 0; i < this->radius_nodes; i++) {
				double x = this->height*float(n) / iHeight;
				double y = this->yCenter + this->radius * cos(-M_PI * 2 * (iRadius - i) / iRadius); // +shift_step);
				double z = this->zCenter + this->radius * sin(-M_PI * 2 * (iRadius - i) / iRadius); // +shift_step);

				this->nodes[i + n*this->radius_nodes]->x = x;
				this->nodes[i + n*this->radius_nodes]->xRef = x;
				this->nodes[i + n*this->radius_nodes]->xVel = 0;
				this->nodes[i + n*this->radius_nodes]->xForce = 0;

				this->nodes[i + n*this->radius_nodes]->y = y;
				this->nodes[i + n*this->radius_nodes]->yRef = y;
				this->nodes[i + n*this->radius_nodes]->yVel = 0;
				this->nodes[i + n*this->radius_nodes]->yForce = 0;

				this->nodes[i + n*this->radius_nodes]->z = z;
				this->nodes[i + n*this->radius_nodes]->zRef = z;
				this->nodes[i + n*this->radius_nodes]->zVel = 0;
				this->nodes[i + n*this->radius_nodes]->zForce = 0;

				/*
				ix = abs(n - cx);
				iy = abs(i - cy);
				rc = sqrt((double)(ix*ix + iy*iy));
				if (rc <= Rpip)
				this->nodes[i + n * this->radius_nodes]->stretchingStiffness = stiffness_plus; //+(this->stiffness - stiffness_plus)*sqrt(rc);    //sin(rc*0.5*M_PI / Rpip);
				else
				*/


				//printf("left > %d right > %d height > %d\n",left, right, this->height);

				int ind_mid_circles = this->height_nodes/2;
				int ind_mid_incircle = this->radius_nodes /2;

				double coord_mid_circles = ind_mid_circles*this->step;
				double coord_mid_incircle = ind_mid_incircle*this->step;

				double c1 = this->stiffness;
				double c2 = 100;//this->stiffness*0.0075;


				double rad = 25;
				//((j - M / 2.0)*(j - M / 2.0) + (k - K / 2.0)*(k - K / 2.0) < (double(M) / 10.1)*(double(M) / 10.1))
				if ((n - ind_mid_circles)*(n - ind_mid_circles) + (i - ind_mid_incircle)*(i - ind_mid_incircle) < rad*rad)
				{
					//linear pyramid
					//double r = sqrt((n - ind_mid_circles)*(n - ind_mid_circles) + (i - ind_mid_incircle)*(i - ind_mid_incircle));
					//this->nodes[i + n * this->radius_nodes]->stretchingStiffness = r/rad*c1+(rad-r)/rad*c2;

					//linear trapezoid
					double r = sqrt((n - ind_mid_circles)*(n - ind_mid_circles) + (i - ind_mid_incircle)*(i - ind_mid_incircle));
					if ((n - ind_mid_circles)*(n - ind_mid_circles) + (i - ind_mid_incircle)*(i - ind_mid_incircle) > (rad - 10)*(rad - 10))
					{
						this->nodes[i + n * this->radius_nodes]->stretchingStiffness = r / rad*c1 + (rad - r) / rad*c2;
					}
					else
						this->nodes[i + n * this->radius_nodes]->stretchingStiffness =  c2;

					//sin
					//double loc = sin(M_PI*(n-ind_mid_circles))*sin(M_PI*(i - ind_mid_incircle));
					//this->nodes[i + n * this->radius_nodes]->stretchingStiffness = loc*c1 + (1.0-loc)*c2;

				}
				else this->nodes[i + n * this->radius_nodes]->stretchingStiffness = c1;
			}
		}
		
	}
	else if (tube_type == 7)
	{
		this->step = 0.01; // 0.069;
		double iPart = 0.;
		//this->height = 1.0;

		modf(this->height / this->step, &iPart);
		int endIndex = (int)iPart;
		this->height_nodes = endIndex + 1;

		modf(this->radius*2.*M_PI / this->step, &iPart);
		endIndex = (int)iPart;
		this->radius_nodes = endIndex;
		this->nodesCount = this->height_nodes*this->radius_nodes;
		this->nodes = new TNode*[this->nodesCount];

		//*
		double begin_center = 0.5;
		//this->yCenter = 0.522;
		//this->zCenter = 0.524;
		//printf(" >> this->yCenter = %lf, this->zCenter = %lf\n", this->yCenter, this->zCenter);
		//*/
		/*double begin_center = 0.319;
		this->yCenter = 0.319;
		this->zCenter = 0.319;
		*/
		for (int n = 0; n < this->height_nodes; ++n) {
			for (int i = 0; i < this->radius_nodes; i++) {
				this->nodes[i + n*this->radius_nodes] = new TFixedElasticNode();
				//this->nodes[i + n * 120]->;
			}
		}


		double rc = 0, stiffness_plus = this->stiffness / 30;
		double iRadius = (double) this->radius_nodes;  //InCircle = fBoundary->radius_nodes;
		double iHeight = (double) this->height_nodes - 1.; //Circles = fBoundary->height_nodes;
		int cx = this->height_nodes / 2., cy = this->radius_nodes / 4., Rpip = 0.00; // this->radius_nodes / 10.,
		int ix, iy;

		for (int n = 0; n < this->height_nodes; ++n) {
			for (int i = 0; i < this->radius_nodes; i++) {
				double x = this->height*float(n) / iHeight;
				double y = this->yCenter + this->radius * cos(-M_PI * 2 * (iRadius - i) / iRadius); // +shift_step);
				double z = this->zCenter + this->radius * sin(-M_PI * 2 * (iRadius - i) / iRadius); // +shift_step);

				this->nodes[i + n*this->radius_nodes]->x = x;
				this->nodes[i + n*this->radius_nodes]->xRef = x;
				this->nodes[i + n*this->radius_nodes]->xVel = 0;
				this->nodes[i + n*this->radius_nodes]->xForce = 0;

				this->nodes[i + n*this->radius_nodes]->y = y;
				this->nodes[i + n*this->radius_nodes]->yRef = y;
				this->nodes[i + n*this->radius_nodes]->yVel = 0;
				this->nodes[i + n*this->radius_nodes]->yForce = 0;

				this->nodes[i + n*this->radius_nodes]->z = z;
				this->nodes[i + n*this->radius_nodes]->zRef = z;
				this->nodes[i + n*this->radius_nodes]->zVel = 0;
				this->nodes[i + n*this->radius_nodes]->zForce = 0;


				//double c1 = this->stiffness;
								
				this->nodes[i + n * this->radius_nodes]->stretchingStiffness = this->stiffness;

			}
		}
	}
	else if (tube_type == 8)
	{
		
		xyz* coord;
		xyz* coord_r;
		xyz* restructed_coord;
		xyz* output_coord;
		xyz* outline_cylinder;
		
		stl_reader::StlMesh <float, unsigned int> mesh ("geometry.stl");
		//stl_reader::StlMesh <float, unsigned int> mesh ("Shell_1.stl");
		int num_vert = mesh.num_vrts();
		int count_z = 0;
		double angle = 6;
		double len = 0, step_z = 0;
		int ieo;
		int inCircle = 360.0/angle;
		float lmax = -1, lmin = -1;
		char dir;
		xyz cp;
		cp = 0;
		{
			//read stl and write to struct
			cout<<"Read stl file"<<endl;
			coord = new xyz [num_vert];
			coord_r = new xyz [num_vert];

			const float* c;

			for(size_t ivert = 0; ivert < mesh.num_vrts(); ivert++)
			{
				c = mesh.vrt_coords(ivert);
				coord_r[ivert].x = c[0]; coord_r[ivert].y = c[1]; coord_r[ivert].z = c[2];
			}
		}
		//forming circles
		{
			//sorting by points circle
			double maxz=coord_r[0].z,minz=coord_r[0].z;
			double maxx=coord_r[0].x,minx=coord_r[0].x;
			double maxy=coord_r[0].y,miny=coord_r[0].y;
			float lmax;
			for(auto i = 0; i<num_vert;++i)
			{
				if(maxz < coord_r[i].z) maxz = coord_r[i].z;
				if(minz > coord_r[i].z) minz = coord_r[i].z;
				if(maxy < coord_r[i].y) maxy = coord_r[i].y;
				if(miny > coord_r[i].y) miny = coord_r[i].y;
				if(maxx < coord_r[i].x) maxx = coord_r[i].x;
				if(minx > coord_r[i].x) minx = coord_r[i].x;
			}
			//Define vessel orientation and orient them to Z-orientation\n
			if( (fabs(maxz - minz) > fabs(maxx - minx)) && (fabs(maxz - minz) > fabs(maxy - miny)))
			{
				printf("Shape is Z-oriended\n");
				dir = 'z';
				lmax = maxz; lmin = minz;
				for(int i=0;i<num_vert;i++)
				{
					coord[i].x=coord_r[i].x;
					coord[i].y=coord_r[i].y;
					coord[i].z=coord_r[i].z;
					//printf("c: (%6.3f, %6.3f, %6.3f) | c_r: (%6.3f, %6.3f, %6.3f)\n",coord[i].x,coord[i].y,coord[i].z,coord_r[i].x,coord_r[i].y,coord_r[i].z);
				}
				lmax = maxz;
			}else if ((fabs(maxy - miny) > fabs(maxx - minx)) && (fabs(maxy - miny) > fabs(maxz - minz)))
			{
				printf("Shape is Y-oriended\n");
				dir = 'y';
				lmax = maxy; lmin = miny;
				for(int i=0;i<num_vert;i++)
				{
					coord[i].x=coord_r[i].z;
					coord[i].y=coord_r[i].x;
					coord[i].z=coord_r[i].y;
					//printf("c: (%6.3f, %6.3f, %6.3f) | c_r: (%6.3f, %6.3f, %6.3f)\n",coord[i].x,coord[i].y,coord[i].z,coord_r[i].x,coord_r[i].y,coord_r[i].z);
				}
				lmax = maxy;
			}else if ((fabs(maxx - minx) > fabs(maxz - minz)) && (fabs(maxx - minx) > fabs(maxy - miny)))
			{
				printf("Shape is X-oriended\n");
				dir = 'x';
				lmax = maxx; lmin = minx;
				for(int i=0;i<num_vert;i++)
				{
					coord[i].x=coord_r[i].z;
					coord[i].y=coord_r[i].y;
					coord[i].z=coord_r[i].x;
					//printf("c: (%6.3f, %6.3f, %6.3f) | c_r: (%6.3f, %6.3f, %6.3f)\n",coord[i].x,coord[i].y,coord[i].z,coord_r[i].x,coord_r[i].y,coord_r[i].z);
				}
				lmax = maxx;
			}else throw("Unexpected error in shape");
			len = fabs(lmax - lmin);
			//*
			quickSort(coord, 0, num_vert - 1,'z','u');
			for(auto i = 1; i<num_vert;++i)
			{
				if(
					(step_z < fabs(coord[i-1].z - coord[i].z)) &&
					fabs(coord[i-1].z - coord[i].z) > 0
					) step_z = fabs(coord[i-1].z - coord[i].z);
			}
			count_z = std::ceil(len/step_z);//*/
			printf("max: %f | min: %f | len: %f\n",lmax, lmin, len);

			int iend = 0;
			int ibegin = 0;
			float curr_z = coord[0].z;

			for(ieo=0;curr_z < lmax;ieo++ ) //count of circles
			{
				curr_z = coord[ibegin].z;
				for(int i=0;i<num_vert;i++)
				{
					if( (coord[i].z < curr_z + (step_z/3 ))
					 && (coord[i].z > curr_z - (step_z/3 )))	iend=i;
				}
				ibegin=iend+1;
			}
			cout<<ieo<<" * "<< inCircle<<" = " <<ieo*inCircle<<endl;
			ibegin=iend=curr_z=0;
			restructed_coord = new xyz [ieo*inCircle];
			output_coord = new xyz [ieo*inCircle];
			outline_cylinder = new xyz [ieo*inCircle];

			//FILE* f1 = fopen("line.txt","w");
			for(int iieo = 0; iieo < ieo;iieo++ ) //restruction loop
			{
				curr_z = coord[ibegin].z;
				xyz midpoint;
				for(int i=0;i<num_vert;i++)
				{
					if( (coord[i].z <= curr_z + step_z/2 )
					 && (coord[i].z >= curr_z - step_z/2 ))	iend=i;
				}
				findMassCenter(coord, midpoint,ibegin,iend);
				sortByRadius(coord, midpoint,ibegin,iend);
				cp.x+=midpoint.x;
				cp.y+=midpoint.y;
				cp.z+=midpoint.z;

				//if(iieo == 0 || curr_z==lmax) printf("Mass Center %3d: x: %8.6f | y: %8.6f | z: %8.6f\n",iieo,midpoint.x,midpoint.y,midpoint.z);
				//fprintf(f1,"%f %f %f\n",midpoint.x,midpoint.y,midpoint.z);
				int loc_start_iter;
				int loc_end_iter;
				for(int ii1 = iieo*inCircle; ii1<(iieo+1)*inCircle;ii1++)   //loop by all points in restructed shape. Circles*inCircle
				{
					if(ii1 == iieo*inCircle )
					{
						outline_cylinder[ii1].x = 0.;
						outline_cylinder[ii1].y = 1.2*maxy;
						outline_cylinder[ii1].z = curr_z;
					}else
					{
						float betta = angle * M_PI/180.0;
						outline_cylinder[ii1].x = midpoint.x + (outline_cylinder[ii1-1].x - midpoint.x)*cos(betta) - (outline_cylinder[ii1-1].y - midpoint.y)*sin(betta);
						outline_cylinder[ii1].y = midpoint.y + (outline_cylinder[ii1-1].x - midpoint.x)*sin(betta) + (outline_cylinder[ii1-1].y - midpoint.y)*cos(betta);
						outline_cylinder[ii1].z = curr_z;
						}
				}
				for(int ii1 = iieo*inCircle; ii1<(iieo+1)*inCircle;ii1++)   //loop by all points in restructed shape. Circles*inCircle
				{
					for(int ci = ibegin; ci<iend; ci++) // loop by incircle of original shape
					{
						if(ci == iend-1)
						{
							loc_start_iter = ci;
							loc_end_iter = ibegin;
						}else{
							loc_start_iter=ci;
							loc_end_iter=ci+1;
						}
						if(isIntersect(midpoint,outline_cylinder[ii1],coord[loc_start_iter],coord[loc_end_iter]))
						{
							restructed_coord[ii1] = findIntersection(midpoint,outline_cylinder[ii1],coord[loc_start_iter],coord[loc_end_iter]);
							break;
						}
					}
				}
				ibegin = iend+1;
			}
			//fclose(f1);
		}
		//restore orientation
		{

			if(dir == 'x')
			{
				printf("Move back to X-orient\n");
				for(int i=0;i<ieo*inCircle;i++)
				{
					output_coord[i].x = restructed_coord[i].z;
					output_coord[i].y = restructed_coord[i].y;
					output_coord[i].z = restructed_coord[i].x;
				}
			}else if(dir == 'y')
			{
				printf("Move back to Y-orient\n");
				for(int i=0;i<ieo*inCircle;i++)
				{
					output_coord[i].x = restructed_coord[i].z;
					output_coord[i].y = restructed_coord[i].x;
					output_coord[i].z = restructed_coord[i].y;
				}
			}else{ printf("Already Z-orient\n"); }
						//*
			cp.x/=(ieo*inCircle);
			cp.y/=(ieo*inCircle);
			cp.z/=(ieo*inCircle);
			double tx, ty, tz;
			double tx1, ty1, tz1;				

			for(int i=0;i<ieo*inCircle;i++)
			{
				//float betta = 184 * M_PI/180.0;
				//float betta1 = -2 * M_PI/180.0;
				float betta = 180 * M_PI/180.0;
				float betta1 = 0 * M_PI/180.0;
				//tx = output_coord[i].x;
				//ty = cp.y + (output_coord[i].y - cp.y)*cos(betta) - (output_coord[i].z - cp.z)*sin(betta);
				//tz = cp.z + (output_coord[i].y - cp.y)*sin(betta) + (output_coord[i].z - cp.z)*cos(betta);
				tx1 = cp.x + (output_coord[i].x - cp.x)*cos(betta) + (output_coord[i].z - cp.z)*sin(betta);
				ty1 = output_coord[i].y;
				tz1 = cp.z - (output_coord[i].x - cp.x)*sin(betta) + (output_coord[i].z - cp.z)*cos(betta);
				
				tx = cp.x + (tx1 - cp.x)*cos(betta1) - (ty1 - cp.y)*sin(betta1);
				ty = cp.y + (tx1 - cp.x)*sin(betta1) + (ty1 - cp.y)*cos(betta1);
				tz = tz1;
				
				output_coord[i].x = tx;
				output_coord[i].y = ty;
				output_coord[i].z = tz;

			
				//output_coord[i].x += 4.13;
				//output_coord[i].y -= 0.05;
				//output_coord[i].z += 1;
				//output_coord[i].x += 0.009624+len;
				//output_coord[i].y += 0.63;
				//output_coord[i].z += 0.54;
			}//*/

			double min_x;
			int i_start, i_end;
			for(int i=0;i<ieo;i++)
			{
				min_x = output_coord[i*ieo].x;
				if(min_x > output_coord[i*inCircle].x)
				{
					min_x = output_coord[i*inCircle].x;
					i_end = ieo*inCircle;
					//i_end = ieo*inCircle + j;
				}
				i_start = i*inCircle - inCircle;
			}			
			printf("i_start and i_end is: (%6d %6d)\n",i_start, i_end);
			xyz midpoint;
			xyz midpoint2;
			findMassCenter(output_coord, midpoint,i_start,i_end);
			//findMassCenter(output_coord, midpoint,0,inCircle);
			printf("center of inflow to vessel before is: (%8.6f %8.6f %8.6f)\n",midpoint.x, midpoint.y, midpoint.z);
			//*
			double shift_by_x = len;
			double shift_by_y = this->GetYCenter() - midpoint.y;
			double shift_by_z = this->GetZCenter() - midpoint.z;
			printf("shifts is: (%8.6f %8.6f %8.6f)\n",shift_by_x, shift_by_y, shift_by_z);
			for(int i=0;i<ieo*inCircle;i++)
			{
				output_coord[i].x += shift_by_x;
				output_coord[i].y += shift_by_y;
				output_coord[i].z += shift_by_z;
			}//*/
			findMassCenter(output_coord, midpoint,i_start,i_end);
			//findMassCenter(output_coord, midpoint,ieo*inCircle,ieo*inCircle + inCircle);
			//findMassCenter(output_coord, midpoint,0,inCircle);
			printf("center of inflow to vessel after  is: (%8.6f %8.6f %8.6f)\n",midpoint.x, midpoint.y, midpoint.z);
		}

		//debug output to file
		if(0)
		{
			cout<<"Print to txt file"<<endl;
			FILE* h = fopen("restructed.txt","w");
			FILE* h2 = fopen("ouput.txt","w");
			FILE* h1 = fopen("outline_cylinder.txt","w");
			for(int circles = 0; circles<ieo;circles++)
			{
				for(int iC = 0; iC<inCircle;iC++)
				{
					fprintf(h,"%f %f %f\n",restructed_coord[circles*inCircle + iC].x,restructed_coord[circles*inCircle + iC].y,restructed_coord[circles*inCircle + iC].z);
					fprintf(h2,"%f %f %f\n",output_coord[circles*inCircle + iC].x,output_coord[circles*inCircle + iC].y,output_coord[circles*inCircle + iC].z);
					fprintf(h1,"%f %f %f\n",outline_cylinder[circles*inCircle + iC].x,outline_cylinder[circles*inCircle + iC].y,outline_cylinder[circles*inCircle + iC].z);
				}
			}
			fclose(h);
			fclose(h1);
			fclose(h2);
			FILE* f = fopen("test_coord.txt","w");
			fprintf(f,"%d\n",num_vert);
			for(int i = 0; i<num_vert;i++)
				fprintf(f,"%16.12f %16.12f %16.12f\n",coord[i].x,coord[i].y,coord[i].z);
			fclose(f);
		}
	
		delete []coord_r;
		delete []coord;
		delete []outline_cylinder;
		delete []restructed_coord;
		//*
		
		this->step = step_z; // 0.069;
		double iPart = 0.;
		//this->height = 1.0;

		modf(this->height / this->step, &iPart);
		int endIndex = (int)iPart;
		this->height_nodes = ieo; //Circles

		modf(this->radius*2.*M_PI / this->step, &iPart);
		endIndex = (int)iPart;
		this->radius_nodes = inCircle;
		this->nodesCount = this->height_nodes*this->radius_nodes;
		this->nodes = new TNode*[this->nodesCount];

		//*
		double begin_center = 0.5;

		for (int n = 0; n < this->height_nodes; ++n) {
			for (int i = 0; i < this->radius_nodes; i++) {
				this->nodes[i + n*this->radius_nodes] = new TFixedElasticNode();
				//this->nodes[i + n * 120]->;
			}
		}


		double rc = 0, stiffness_plus = this->stiffness / 30;
		double iRadius = (double) this->radius_nodes;  //InCircle = fBoundary->radius_nodes;
		double iHeight = (double) this->height_nodes - 1.; //Circles = fBoundary->height_nodes;
		//int cx = this->height_nodes / 2., cy = this->radius_nodes / 4., Rpip = 0.00; // this->radius_nodes / 10.,
		int ix, iy;

		for (int n = 0; n < this->height_nodes; ++n) {
			for (int i = 0; i < this->radius_nodes; i++) {
				double x = output_coord[i + n*this->radius_nodes].x;
				double y = output_coord[i + n*this->radius_nodes].y;
				double z = output_coord[i + n*this->radius_nodes].z;

				this->nodes[i + n*this->radius_nodes]->x = x;
				this->nodes[i + n*this->radius_nodes]->xRef = x;
				this->nodes[i + n*this->radius_nodes]->xVel = 0;
				this->nodes[i + n*this->radius_nodes]->xForce = 0;


				this->nodes[i + n*this->radius_nodes]->y = y;
				this->nodes[i + n*this->radius_nodes]->yRef = y;
				this->nodes[i + n*this->radius_nodes]->yVel = 0;
				this->nodes[i + n*this->radius_nodes]->yForce = 0;

				this->nodes[i + n*this->radius_nodes]->z = z;
				this->nodes[i + n*this->radius_nodes]->zRef = z;
				this->nodes[i + n*this->radius_nodes]->zVel = 0;
				this->nodes[i + n*this->radius_nodes]->zForce = 0;

				this->nodes[i + n*this->radius_nodes]->p = 0;
				//double c1 = this->stiffness;
								
				this->nodes[i + n * this->radius_nodes]->stretchingStiffness = this->stiffness;

			}
		}
	}
	/*
	for (int kk = 0; kk < 2; kk++)
	{
		for (int n = 0; n < this->height_nodes; ++n)
		{
			for (int i = 0; i < this->radius_nodes; i++)
			{
				if ((n - this->height_nodes / 2)*(n - this->height_nodes / 2) + (i - this->radius_nodes / 2)*(i - this->radius_nodes / 2) < 25*25)
					this->nodes[i + n * this->radius_nodes]->stretchingStiffness = (this->nodes[i + (n-	1) * this->radius_nodes]->stretchingStiffness + this->nodes[i + (n+1) * this->radius_nodes]->stretchingStiffness)/2;
			}
		}
	}
	//*/
	for (int n = 0; n < this->height_nodes; ++n) {
		for (int i = 0; i < this->radius_nodes; i++) {
			this->nodes[i + n*this->radius_nodes]->boundary = this;

			if (n <  this->height_nodes - 1)
			{
				this->nodes[i + n*this->radius_nodes]->neighbors.next = this->nodes[i + (n + 1) * this->radius_nodes];
				this->nodes[i + n*this->radius_nodes]->neighbors.initialDistanceNext =
					TNode::distance_3d(this->nodes[i + n * this->radius_nodes], this->nodes[i + n * this->radius_nodes]->neighbors.next);
			}
			if (n > 0)
			{
				this->nodes[i + n*this->radius_nodes]->neighbors.prev = this->nodes[i + (n - 1) * this->radius_nodes];
				this->nodes[i + n*this->radius_nodes]->neighbors.initialDistancePrev = TNode::distance_3d(this->nodes[i + n * this->radius_nodes], this->nodes[i + n * this->radius_nodes]->neighbors.prev);
			}

			if (this->nodes[i + n*this->radius_nodes]->neighbors.prev != NULL && this->nodes[i + n * this->radius_nodes]->neighbors.next != NULL)
			{
				this->nodes[i + n*this->radius_nodes]->neighbors.initialCurvatures[OX] = TNode::curvature(
					this->nodes[i + n*this->radius_nodes],
					this->nodes[i + n*this->radius_nodes]->neighbors.next,
					this->nodes[i + n*this->radius_nodes]->neighbors.prev,
					TNode::Ox
					);
				this->nodes[i + n*this->radius_nodes]->neighbors.initialCurvatures[OY] = TNode::curvature(
					this->nodes[i + n*this->radius_nodes],
					this->nodes[i + n*this->radius_nodes]->neighbors.next,
					this->nodes[i + n*this->radius_nodes]->neighbors.prev,
					TNode::Oy
					);
				this->nodes[i + n*this->radius_nodes]->neighbors.initialCurvatures[OZ] = TNode::curvature(
					this->nodes[i + n*this->radius_nodes],
					this->nodes[i + n*this->radius_nodes]->neighbors.next,
					this->nodes[i + n*this->radius_nodes]->neighbors.prev,
					TNode::Oz
					);
			}
		}
	}

	this->n_pathes = 6;
	this->max_path_step = 1000;
	this->out_path_step =30;

	this->X_pathes = new double* [this->n_pathes];
	this->Y_pathes = new double* [this->n_pathes];
	this->Z_pathes = new double* [this->n_pathes];

	for (int i = 0; i < this->n_pathes; i++)
	{
		this->X_pathes[i] = new double [this->max_path_step];
		this->Y_pathes[i] = new double [this->max_path_step];
		this->Z_pathes[i] = new double [this->max_path_step];
	}

	for (int i = 0; i < this->n_pathes; i++)
	for (int j = 0; j < this->max_path_step; j++)
	{
		this->X_pathes[i][j] = 0.0;
		this->Y_pathes[i][j] = 0.0;
		this->Z_pathes[i][j] = 0.0;
	}

	this->Y_pathes[0][0] = 0.4;
	this->Z_pathes[0][0] = 0.495;

	this->Y_pathes[1][0] = 0.448;
	this->Z_pathes[1][0] = 0.482;

	this->Y_pathes[2][0] = 0.482;
	this->Z_pathes[2][0] = 0.448;

	this->Y_pathes[3][0] = 0.495;
	this->Z_pathes[3][0] = 0.4;

	this->Y_pathes[4][0] = 0.467;
	this->Z_pathes[4][0] = 0.333;

	this->Y_pathes[5][0] = 0.4;
	this->Z_pathes[5][0] = 0.305;


    return this;
}

double TCylinderElasticBoundary::GetArea()
{
    return 2 * M_PI * this->radius * 1 / this->nodesCount;
}


