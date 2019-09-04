#ifndef IMMERSED_BOUNDARY_H
#define IMMERSED_BOUNDARY_H

#include <map>
#include <tuple>
#include <cmath>
#include <unordered_map>
#include <functional>

using std::map;
using std::tuple;
using std::unordered_map;

#define OX 1
#define OY 2
#define OZ 3

#define AXIS int

#define ELASTIC_NODE 0
#define FIXED_NODE 1

class TNode;
class TImmersedBoundary;

class TNode
{
    public:
        double x,y,z;
        double x_0,y_0,z_0;
        double xPrev,yPrev,zPrev;
        double xRef, yRef, zRef;
        double xForce, yForce, zForce;
        double xVel, yVel, zVel;
		double Nx, Ny, Nz;
		double p;
		double concentration;
		double vis;
		double ShearStressX,ShearStressY,ShearStressZ,ShearStress;
    	int type;
		double stretchingStiffness;

        TImmersedBoundary *boundary;
        struct Neighbor {
            TNode *next;
            double initialDistanceNext;
            TNode *prev;
            double initialDistancePrev;
            unordered_map<AXIS, double> initialCurvatures;
        } neighbors;

        // debug
        map<tuple<int, int, int>, TNode*> surround;

        TNode();
        virtual double GetStiffnessX(int timestep) {return 1;};
        virtual double GetStiffnessY(int timestep) {return 1;};
        virtual double GetStiffnessZ(int timestep) {return 1;};

        struct Getter {
            std::function<double(TNode*)> coord;
            std::function<double(TNode*)> ref;
            int type;
        };

        static double getX(TNode* node) {return node->x;};
        static double getY(TNode* node) {return node->y;};
        static double getZ(TNode* node) {return node->z;};        
		
		static double getX_0(TNode* node) {return node->x_0;};
        static double getY_0(TNode* node) {return node->y_0;};
        static double getZ_0(TNode* node) {return node->z_0;};

        static double getXRef(TNode* node) {return node->xRef;};
        static double getYRef(TNode* node) {return node->yRef;};
        static double getZRef(TNode* node) {return node->zRef;};

        static Getter Ox;
        static Getter Oy;
        static Getter Oz;

        static double distance_3d(TNode* P1, TNode* P2)
        {
            double DX = P1->x - P2->x;
            double DY = P1->y - P2->y;
            double DZ = P1->z - P2->z;

            return sqrt(DX * DX + DY * DY + DZ * DZ);
        }

        static double curvature(TNode* master, TNode* next, TNode* prev, TNode::Getter getters)
        {
            return getters.coord(next) - 2 * getters.coord(master) + getters.coord(prev);
        }

        virtual double GetStretchingForce(TNode *neighbor, TNode::Getter getter) {return 0;};
        virtual double GetBendingForce(TNode *next, TNode *prev, TNode::Getter getter) {return 0;};
        virtual double GetTargetForce(TNode::Getter getter) {return 0;};

        ~TNode() {};
        };

class TImmersedBoundary
{
    public:
        double stiffness;
        double stretchingStiffness = 10000;//10000;
        double bendingStiffness = 0;//10000;//20; //1800;//15000;
		double step=0.01;
        int nodesCount;
		int radius_nodes, height_nodes;
        TNode **nodes;
		int max_path_step, out_path_step;
		int n_pathes;
		double **X_pathes;
		double **Y_pathes;
		double **Z_pathes;

		double radius;
		int nSpl;
		double *xSpl, *ySpl, *zSpl;

        TImmersedBoundary(){};
        ~TImmersedBoundary();
        virtual double GetArea() = 0;
		virtual double GetYCenter() = 0;
		virtual double GetZCenter() = 0;
};

class TSphereBoundary: public TImmersedBoundary
{
    public:
        double radius;
        double xCenter;
        double yCenter;
        double zCenter;

        TSphereBoundary(): TImmersedBoundary(){};
        ~TSphereBoundary() {}; 
        TSphereBoundary* Initialize();
        TSphereBoundary* withStiffness(double stiffness) { this->stiffness = stiffness; return this; };
        TSphereBoundary* withNodesCount(int lNodesCount) { this->nodesCount = lNodesCount; return this; };
        TSphereBoundary* withRadius(double lRadius) { this->radius = lRadius; return this; };
        TSphereBoundary* withCenter(
            double lXCenter,
            double lYCenter,
            double lZCenter
        ) { this->xCenter = lXCenter; this->yCenter = lYCenter; this->zCenter = lZCenter; return this; };

        double GetArea();
		double GetYCenter() {return yCenter;};;
		double GetZCenter() {return yCenter;};;
};


class TSimpleValve: public TImmersedBoundary
{
    public:
        double radius;
        double xCenter;
        double height;

        TSimpleValve(): TImmersedBoundary(){};
        ~TSimpleValve() {}; 
        TSimpleValve* Initialize();
        TSimpleValve* withStiffness(double stiffness) { this->stiffness = stiffness; return this; };
        TSimpleValve* withNodesCount(int lNodesCount) { this->nodesCount = lNodesCount; return this; };
        TSimpleValve* withRadius(double lRadius) { this->radius = lRadius; return this; };
        TSimpleValve* withCenter(double lXCenter) { this->xCenter = lXCenter; return this;};
        TSimpleValve* withHeight(double lHeight) { this->height = lHeight; return this;};

        double GetArea();
		double GetYCenter() {};
		double GetZCenter() {};
};


class TCylinderBoundary: public TImmersedBoundary
{
    public:
        double radius;
        double yCenter;
        double zCenter;
        double height;

        TCylinderBoundary(): TImmersedBoundary(){};
        ~TCylinderBoundary() {};
        TCylinderBoundary* Initialize();
        TCylinderBoundary* withStiffness(double stiffness) { this->stiffness = stiffness; return this; };
        TCylinderBoundary* withNodesCount(int lNodesCount) { this->nodesCount = lNodesCount; return this; };
        TCylinderBoundary* withRadius(double lRadius) { this->radius = lRadius; return this; };
        TCylinderBoundary* withCenter(double lYCenter, double lZCenter) { this->yCenter = lYCenter; this->zCenter = lZCenter;return this;};

        double GetArea();
		double GetYCenter() {return yCenter;};;
		double GetZCenter() {return yCenter;};;
};

class TDeformedCylinderBoundary: public TImmersedBoundary
{
    public:
        double radius;
        double deformation;
        double yCenter;
        double zCenter;
        double height;

        TDeformedCylinderBoundary(): TImmersedBoundary(){};
        ~TDeformedCylinderBoundary() {};
        TDeformedCylinderBoundary* Initialize();
        TDeformedCylinderBoundary* withStiffness(double stiffness) { this->stiffness = stiffness; return this; };
        TDeformedCylinderBoundary* withNodesCount(int lNodesCount) { this->nodesCount = lNodesCount; return this; };
        TDeformedCylinderBoundary* withRadius(double lRadius) { this->radius = lRadius; return this; };
        TDeformedCylinderBoundary* withDeformation(double lDeformation) { this->deformation = lDeformation; return this; };
        TDeformedCylinderBoundary* withCenter(double lYCenter, double lZCenter) { this->yCenter = lYCenter; this->zCenter = lZCenter;return this;};

        double GetArea();
		double GetYCenter() {return yCenter;};
		double GetZCenter() {return yCenter;};
};

class TCylinderNode: public TNode
{
    public:
        TCylinderNode():TNode() {};
        ~TCylinderNode() {};
        double GetStiffnessX(int timestep);
        double GetStiffnessY(int timestep);
        double GetStiffnessZ(int timestep);
};


class TElasticNode: public TNode
{
    public:
        int type = ELASTIC_NODE;
        TElasticNode():TNode() {};
        ~TElasticNode() {};

        virtual double GetStretchingForce(TNode *neighbor, TNode::Getter getter);
        virtual double GetBendingForce(TNode *next, TNode *prev, TNode::Getter getter);
        virtual double GetTargetForce(TNode::Getter getter);
};

class TTargetNode: public TNode
{
    public:
        TTargetNode():TNode() {};
        ~TTargetNode() {};

        virtual double GetStretchingForce(TNode *neighbor, TNode::Getter getter);
        virtual double GetBendingForce(TNode *next, TNode *prev, TNode::Getter getter);
        virtual double GetTargetForce(TNode::Getter getter);
};

class TFixedElasticNode: public TNode
{
    public:
    	int type = FIXED_NODE;
		//double stretchingStiffness;

        TFixedElasticNode():TNode() {};
        ~TFixedElasticNode() {};

        virtual double GetStretchingForce(TNode *neighbor, TNode::Getter getter);
        virtual double GetBendingForce(TNode *next, TNode *prev, TNode::Getter getter);
        virtual double GetTargetForce(TNode::Getter getter);
};

class TFixedStretchingNode: public TNode
{
    public:
        TFixedStretchingNode():TNode() {};
        ~TFixedStretchingNode() {};

        virtual double GetStretchingForce(TNode *neighbor, TNode::Getter getter);
        virtual double GetBendingForce(TNode *next, TNode *prev, TNode::Getter getter);
        virtual double GetTargetForce(TNode::Getter getter);
};


class TCylinderElasticBoundary: public TImmersedBoundary
{
    public:
        //double radius;
        double yCenter;
        double zCenter;
        double height;
		double fRadius;
		//int nSpl;
		//double *xSpl, *ySpl, *zSpl;

        TCylinderElasticBoundary(): TImmersedBoundary(){};
        ~TCylinderElasticBoundary() {};
        TCylinderElasticBoundary* Initialize();
        TCylinderElasticBoundary* withStiffness(double stiffness) { this->stiffness = stiffness; return this; };
        TCylinderElasticBoundary* withNodesCount(int lNodesCount) { this->nodesCount = lNodesCount; return this; };
		TCylinderElasticBoundary* withRadius(double lRadius) { fRadius = lRadius; this->radius = lRadius; return this; };
		TCylinderElasticBoundary* withHeight(double lHeight) { this->height = lHeight; return this; };
		TCylinderElasticBoundary* withStep(double lStep)   { this->step = lStep; return this; };
        TCylinderElasticBoundary* withCenter(double lYCenter, double lZCenter) { this->yCenter = lYCenter; this->zCenter = lZCenter;return this;};
        double GetArea();
		double GetYCenter() {return yCenter;};
		double GetZCenter() {return zCenter;};
};

#endif
