
#include"tsunami.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define Error(a)   femError(a,__LINE__,__FILE__)
#define Warning(a) femWarning(a,  __LINE__, __FILE__)

typedef enum {FEM_TRIANGLE,FEM_EDGE} femElementType;

typedef struct {
    int *elem;
    double *X;
    double *Y;
    double *Z;                //       _____   _____   _   _   __   _       ___       ___  ___   _   
    int nElem;                //      |_   _| /  ___/ | | | | |  \ | |     /   |     /   |/   | | |
    int nNode;                //        | |   | |___  | | | | |   \| |    / /| |    / /|   /| | | |
    int nLocalNode;           //        | |   \___  \ | | | | | |\   |   / / | |   / / |__/ | | | | 
} femMesh;                    //        | |    ___| | | |_| | | | \  |  / /  | |  / /       | | | | 
                              //        |_|   /_____/ \_____/ |_|  \_| /_/   |_| /_/        |_| |_| 
typedef struct {
    int n;    
    const double *xsi;
    const double *eta;
    const double *weight;
} femIntegration;

typedef struct {
    int n;
    int order;
    void (*x2)(double *xsi, double *eta);
    void (*phi2)(double xsi, double eta, double *phi);
    void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
    void (*x1)(double *xsi);
    void (*phi1)(double xsi, double *phi);
} femDiscrete;

                                //..........
typedef struct {                //......................................
   int elem[2];                 //........................................................
   int node[2];                 //...............................................................
}femEdge;                       //.................:///////+oosyhhysyo++/:..............................
                                //.............-/+/-`       `-/+ooooooshhhhso/-..................................
typedef struct {                //...........:++.`      `-/ossyyhho--/++++oooos+-....................................
   femMesh *mesh;               //.........:o/`       .:+ooosyhhy:.....-+o/////+s:....................................
   femEdge *edges;              //.......-+/`       -/+++osyhhhy-........-y//////s:.......................................
   int nEdge;                   //....../+.       ./++++oyhhhhh-..........-ossssss+..........................................
   int nBoundary;               //....-o:       ./++++oyhhhhhh+........-:::///++oo++//:--......................................
} femEdges;                     //...:o.      `:+++++shhhhhhhd.....-///:..``....-/+syhhddhys+/-.................................
                                //../+`     `:+++++syhhhhhhhhh..-/+:.`      `-/osyhhho++oosyysss/-...............................
typedef struct {                //-o/     `:+oooosyhhhhhhhhhhh-/+.`      `-/osssyhhs:..-:++++++++o+.................................
   femMesh *mesh;               //o.  `.-----/+shdddhhhhhhhhhd/.       .:+ooosyhhho.......-/o/////+s..................................
   femEdges *edges;             //  `-:.  `-/+oshhhsooshhhhhs.       .:++++oyhhhho.......`../++ooooss.................................
   femIntegration *rule1d;      // -:`  .:/+osyyyhhhyoosdhh:       `:/+++oyhhhhhy........``..-://///:.................................
   femIntegration *rule2d;      ///.  `:+++/:/ohdhhhhhhhds`      `:/++++shhhhhhh/..........`..........................................
   femDiscrete *space;          //` `-++:.`-+shdhhsohhhd+      `-/++++syhhhhhhhd......................................................
   int size;                    //`-/o+. -/syhdhhhhhhhh:     `-/++++syhhhhhhhhhd.`````````````````````````````````````````.:.````-.```
   double *E;                   //+oo+`./syhhhmhhhhhds.    `-/++++syhhhhhhhhhhhd-````````````````````````````````````````.//-:----+```
   double *U;                   //so--//:/shhhhdhhhh/    `-/++++syhhhhhhhhhhhhhho```````````````````````````````````````-oyyyyhhyhy-``
   double *V;                   //+os+../shhhhydddo`   `-/++++syhhhhhhhhhhhhhhhhh-```````````````````````````````````.:oyyyyyyyyyyyy-`
   double *FE;                  //hy/./syhdhhhhh+.   `-/++++oyhhhhhhhhhhhhhhhhhhhy.```````````````````````````````.-+syyyyyyyyyyyyyyy:
   double *FU;                  //o:+syhhhhdhs/`   `:+oooosyhhhhhhhhhhhhhhhhhhhhhhy:``````````````````````````.-/+syyyyyyyyyyyyyyyyyyy
   double *FV;                  //hhhhhhhyo:.   `-+osssyyhhhhhhhhhhhhhhhhhhhhhhhhhhh+.```````````````````..:/osyyyyyyyyyyyyyyyyyyyyyyy
   double omega;                //oo+/:-`  `.:/osyyhhhhhhhhhhhhhhhyysyyhhhhhhhhhhhhhhy+.``````````...:/+ossyssysssoossssyyssssssssssss
   double gamma;                //----:/oshdmdhddhhhhhhhhhhhhhhhhyssyssyhhhhhhhhhhhhhhhhs:..-::/+oosysssssssssss+ooooossssssssssssssss
   double g;                    //hhhhhhhhhyyhhhhhhhhhhhhhhhhhhhhhsssssyhhhhhhhhhhhhhhhhhyyhhyssssssssssssssssss+ooooossssssssssssssss
   double R;                    //hhhhhhhyyhhhhhhhhhhhhhhhhhhhhhhhhyyhhhhhhhhhhhhhhhhhhyyyhhhhhhhyysssssssssssssssssssssssssssssssssso
   double dt;                   //hhhhhhhyhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhyhhhhhhhhhhhhhhhhhyyyysssssssssssssssssssssyyyy
} femTsunami;                   //hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
                                //hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
                                //+osooso+++++++++++++++++++++++++++++++++++++++++++++++++++++++++oo++++++++++++++++++++++++++++++++++


femMesh             *femMeshRead(const char *filename);
void                 femMeshWrite(const femMesh* myMesh, const char *filename);
void                 femMeshFree(femMesh *theMesh);

femIntegration      *femIntegrationCreate(int n, femElementType type);
void                 femIntegrationFree(femIntegration *theRule);

femEdges*            femEdgesCreate(femMesh *theMesh);
void                 femEdgesFree(femEdges *theEdges);
int                  femEdgesCompare(const void *edgeOne, const void *edgeTwo);
void                 femEdgesMap(femEdges *theEdges, int index, int map[2][2]);

femDiscrete*         femDiscreteCreate(int n, femElementType type);
void                 femDiscreteFree(femDiscrete* mySpace);
void                 femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta);
void                 femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi);
void                 femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);
void                 femDiscreteXsi1(femDiscrete* mySpace, double *xsi);
void                 femDiscretePhi1(femDiscrete* mySpace, double xsi, double *phi);
double               femDiscreteInterpolate(double *phi, double *U, int *map, int n);

void                 femError(char *text, int line, char *file);
void                 femWarning(char *text, int line, char *file);

femTsunami          *femTsunamiCreate(const char *meshFileName, double dt);
void						femTsunamiInit(femTsunami *myTsunami);
void                 femTsunamiCompute(femTsunami* myTsunami);
void                 femTsunamiAddIntegralsElement(femTsunami *myTsunami);
void                 femTsunamiAddIntegralsEdges(femTsunami *myTsunami);
void                 femTsunamiMultiplyInverseMatrix(femTsunami *myTsunami);
void                 femTsunamiFree(femTsunami *myTsunami);

void                 convertTo2D(femTsunami *myTsunami);

static const double _gaussEdge2Xsi[2]     = { 0.577350269189626,-0.577350269189626 };
static const double _gaussEdge2Weight[2]  = { 1.000000000000000, 1.000000000000000 };

static const double _gaussTri3Xsi[3]      = { 0.166666666666667, 0.666666666666667, 0.166666666666667 };
static const double _gaussTri3Eta[3]      = { 0.166666666666667, 0.166666666666667, 0.666666666666667 };
static const double _gaussTri3Weight[3]   = { 0.166666666666667, 0.166666666666667, 0.166666666666667 };

void tsunamiCompute(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName)
{
	femTsunami *myTsunami = femTsunamiCreate(meshFileName, dt);
	int n;
	for(n=0;n<=nmax;++n){
		femTsunamiCompute(myTsunami);
		if(n%sub == 0){
			tsunamiWriteFile(baseResultName, n, myTsunami->U, myTsunami->V, myTsunami->E, myTsunami->mesh->nElem);
		}
	}
	femTsunamiFree(myTsunami);
}

////// FEMTSUNAMI //////

femTsunami *femTsunamiCreate(const char *meshFileName, double dt){
	femTsunami *myTsunami = malloc(sizeof(femTsunami));

	myTsunami->omega = 2*M_PI / 86400.0;
	myTsunami->gamma = pow(10,-7);
	myTsunami->g = 9.81;
	myTsunami->R = 6371220.0;
	myTsunami->dt = dt;

	myTsunami->mesh = femMeshRead(meshFileName);
	myTsunami->edges = femEdgesCreate(myTsunami->mesh);
	myTsunami->rule1d = femIntegrationCreate(2, FEM_EDGE);
	myTsunami->rule2d = femIntegrationCreate(3, FEM_TRIANGLE);
	myTsunami->space = femDiscreteCreate(3, FEM_TRIANGLE);

	int size = myTsunami->mesh->nElem * 3 + 1;
	myTsunami->size = size;
	myTsunami->E = malloc(sizeof(double)*size);
	myTsunami->U = malloc(sizeof(double)*size);
	myTsunami->V = malloc(sizeof(double)*size);
	myTsunami->FE = malloc(sizeof(double)*size);
	myTsunami->FU = malloc(sizeof(double)*size);
	myTsunami->FV = malloc(sizeof(double)*size);


	int i,j;
	for(i=0	;i<myTsunami->mesh->nElem;++i){
		int *mapCoord = &(myTsunami->mesh->elem[3*i]);
		for(j=0;j<3;++j){
			myTsunami->E[3*i+j] = tsunamiInitialConditionOkada(myTsunami->mesh->X[mapCoord[j]],myTsunami->mesh->Y[mapCoord[j]]);
			myTsunami->U[3*i+j] = 0.0;
			myTsunami->V[3*i+j] = 0.0;
		}
	}
	convertTo2D(myTsunami);
	
	return myTsunami;
}

void femTsunamiFree(femTsunami *myTsunami){
	free(myTsunami->E);
	free(myTsunami->U);
	free(myTsunami->V);
	free(myTsunami->FE);
	free(myTsunami->FU);
	free(myTsunami->FV);
	femMeshFree(myTsunami->mesh);
	femEdgesFree(myTsunami->edges);
	femIntegrationFree(myTsunami->rule1d);
	femIntegrationFree(myTsunami->rule2d);
	femDiscreteFree(myTsunami->space);
	free(myTsunami);
}

void femTsunamiCompute(femTsunami* myTsunami){
	int size=myTsunami->size, i;
	double *E = myTsunami->E;
	double *U = myTsunami->U;
	double *V = myTsunami->V;
	double *FE = myTsunami->FE;
	double *FU = myTsunami->FU;
	double *FV = myTsunami->FV;
	double dt = myTsunami->dt;
	for(i=0;i<size;++i){
		FE[i] = 0.0;
		FU[i] = 0.0;
		FV[i] = 0.0;
	}
	femTsunamiAddIntegralsElement(myTsunami);
	femTsunamiAddIntegralsEdges(myTsunami);
	femTsunamiMultiplyInverseMatrix(myTsunami);
	for(i=0;i<size;++i){
		E[i] += dt * FE[i];
		U[i] += dt * FU[i];
		V[i] += dt * FV[i];
	}
}

void femTsunamiAddIntegralsElement(femTsunami *myTsunami){
	double *BE = myTsunami->FE;
	double *BU = myTsunami->FU;
	double *BV = myTsunami->FV;
	double *E = myTsunami->E;
	double *U = myTsunami->U;
	double *V = myTsunami->V;
	double *Y = myTsunami->mesh->Y;
	double *X = myTsunami->mesh->X;
	double *Z = myTsunami->mesh->Z;

	femIntegration *theRule = myTsunami->rule2d;
	femDiscrete *theSpace = myTsunami->space;
	double Xloc[3], Yloc[3], phi[3], dphidxsi[3], dphideta[3], dphidx[3], dphidy[3];
	double xsi, eta, weight, jac;
	double x,y,e,u,v,h;
	int i, j, k, iElem, mapElem[3];

	double  g = myTsunami->g;
	double  gamma = myTsunami->gamma;
	double  R = myTsunami->R;
	double  omega = myTsunami->omega;

	for(iElem=0;iElem<myTsunami->mesh->nElem;++iElem){
		int *mapCoord = &(myTsunami->mesh->elem[3*iElem]);
		for(j=0;j<3;++j){
			mapElem[j] = 3*iElem+j;
			Xloc[j] = myTsunami->mesh->X[mapCoord[j]];
			Yloc[j] = myTsunami->mesh->Y[mapCoord[j]];
		}
		for(k=0;k<theRule->n;++k){
			xsi = theRule->xsi[k];
			eta = theRule->eta[k];
			weight = theRule->weight[k];
			femDiscretePhi2(theSpace,xsi,eta, phi);
			femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);
			double dxdxsi = 0.0;
			double dxdeta = 0.0;
			double dydxsi = 0.0;
			double dydeta = 0.0;
			for(i=0;i<3;++i){
				dxdxsi += Xloc[i]*dphidxsi[i];
				dxdeta += Xloc[i]*dphideta[i];
				dydxsi += Yloc[i]*dphidxsi[i];
				dydeta += Yloc[i]*dphideta[i];
			}
			jac = fabs(dxdxsi*dydeta - dydxsi*dxdeta);
			for(i=0;i<3;++i){
				dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
				dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
			}

			y = femDiscreteInterpolate(phi,Y,mapCoord,3);
			x = femDiscreteInterpolate(phi,X,mapCoord,3);
			h = femDiscreteInterpolate(phi,Z,mapCoord,3);
			e = femDiscreteInterpolate(phi,E,mapElem,3);
			u = femDiscreteInterpolate(phi,U,mapElem,3);
			v = femDiscreteInterpolate(phi,V,mapElem,3);
			
			double lat = h/R;//((4*R*R-x*x-y*y))/(4*R*R+x*x+y*y); == sin(theta)
			const double coriolis = 2*omega*lat;
			const double sphere = ((4*R*R+x*x+y*y)/(4*R*R));

			for(i=0;i<3;++i){
				BE[mapElem[i]] += ((dphidx[i]*h*u+dphidy[i]*h*v)*sphere + (phi[i]*(h*(x*u+y*v)/(R*R))))            * jac*weight;
				BU[mapElem[i]] += ((phi[i]*(coriolis*v-gamma*u)+dphidx[i]*g*e*sphere) + (phi[i]*g*x*e/(2*R*R)))    * jac*weight;
				BV[mapElem[i]] += ((phi[i]*(-coriolis*u-gamma*v) + dphidy[i]*g*e*sphere) + (phi[i]*g*y*e/(2*R*R))) * jac*weight;
			}
		}
	}
}

void femTsunamiAddIntegralsEdges(femTsunami *myTsunami){
	double *BE = myTsunami->FE;
	double *BU = myTsunami->FU;
	double *BV = myTsunami->FV;
	double *E = myTsunami->E;
	double *U = myTsunami->U;
	double *V = myTsunami->V;
	femIntegration *theRule = myTsunami->rule1d;
	femDiscrete *theSpace = myTsunami->space;
	int size = myTsunami->size;
	double R = myTsunami->R;

	double xEdge[2], yEdge[2], hEdge[2], phiEdge[2];
	double xsi, weight, jac;
	double eL,eR,uL,uR,vL,vR,unL,unR,qe,qu,qv;
	double x,y,h,sphere;
	int i,j,k,iEdge,mapEdge[2][2];

	double  g = myTsunami->g;

	for(iEdge=0;iEdge<myTsunami->edges->nEdge;++iEdge){
		femEdgesMap(myTsunami->edges,iEdge,mapEdge);
		for(j=0;j<2;++j){
			int node = myTsunami->edges->edges[iEdge].node[j];
			xEdge[j] = myTsunami->mesh->X[node];
			yEdge[j] = myTsunami->mesh->Y[node];
			hEdge[j] = myTsunami->mesh->Z[node];
		}
		int boundary = (mapEdge[1][0] == size-1);

		double dx = xEdge[1] - xEdge[0];
		double dy = yEdge[1] - yEdge[0];
		double norm = sqrt(dx*dx+dy*dy);
		double nx = dy/norm;
		double ny = -dx/norm;
		jac = norm / 2.0;
		for(k=0;k<2;++k){
			xsi = theRule->xsi[k];
			weight = theRule->weight[k];
			femDiscretePhi1(theSpace, xsi, phiEdge);
			
			x=0.0;
			y=0.0;
			h=0.0;
			for(i=0;i<2;++i){
				x += xEdge[i]*phiEdge[i];
				y += yEdge[i]*phiEdge[i];
				h += hEdge[i]*phiEdge[i];
			}
			eL = femDiscreteInterpolate(phiEdge,E,mapEdge[0],2);
			eR = boundary ? eL : femDiscreteInterpolate(phiEdge,E,mapEdge[1],2);
			uL = femDiscreteInterpolate(phiEdge, U, mapEdge[0], 2);
			uR = femDiscreteInterpolate(phiEdge, U, mapEdge[1], 2);
			vL = femDiscreteInterpolate(phiEdge, V, mapEdge[0], 2);
			vR = femDiscreteInterpolate(phiEdge, V, mapEdge[1], 2);
			unL = uL*nx + vL*ny;
			unR = boundary ? -unL : uR*nx + vR*ny;

			sphere = ((4*R*R+x*x+y*y)/(4*R*R));

			qe = 0.5*h*   ((unL+unR) + sqrt(g/h)*(eL - eR))*sphere;
			qu = 0.5*g*nx*((eL + eR) + sqrt(h/g)*(unL-unR))*sphere;
			qv = 0.5*g*ny*((eL + eR) + sqrt(h/g)*(unL-unR))*sphere;
			
			for(i=0;i<2;++i){
				BE[mapEdge[0][i]] -= qe*phiEdge[i]*jac*weight;
				BU[mapEdge[0][i]] -= qu*phiEdge[i]*jac*weight;
				BV[mapEdge[0][i]] -= qv*phiEdge[i]*jac*weight;
				BE[mapEdge[1][i]] += qe*phiEdge[i]*jac*weight;
				BU[mapEdge[1][i]] += qu*phiEdge[i]*jac*weight;
				BV[mapEdge[1][i]] += qv*phiEdge[i]*jac*weight;
			}
		}
	}	
}

void femTsunamiMultiplyInverseMatrix(femTsunami *myTsunami){
	double *BE = myTsunami->FE;
	double *BU = myTsunami->FU;
	double *BV = myTsunami->FV;
	double invA[3][3] = {{18.0,-6.0,-6.0},{-6.0,18.0,-6.0},{-6.0,-6.0,18.0}};
	double BEloc[3],BVloc[3],BUloc[3];

	double Xloc[3],Yloc[3],jac;
	int i,j,iElem,mapElem[3];

	for(iElem=0;iElem<myTsunami->mesh->nElem;++iElem){
		int *mapCoord = &(myTsunami->mesh->elem[3*iElem]);
		for(j=0;j<3;++j){
			mapElem[j] = 3*iElem + j;
			Xloc[j] = myTsunami->mesh->X[mapCoord[j]];
			Yloc[j] = myTsunami->mesh->Y[mapCoord[j]];
		}
		jac = (Xloc[1]-Xloc[0])*(Yloc[2]-Yloc[0]) - (Yloc[1]-Yloc[0])*(Xloc[2]-Xloc[0]);
		for(i=0;i<3;++i){
			BEloc[i] = BE[mapElem[i]]; BE[mapElem[i]] = 0.0;
			BUloc[i] = BU[mapElem[i]]; BU[mapElem[i]] = 0.0;
			BVloc[i] = BV[mapElem[i]]; BV[mapElem[i]] = 0.0;
		}
		for(i=0;i<3;++i){
			for(j=0;j<3;++j){
				BE[mapElem[i]] += invA[i][j] * BEloc[j] / jac;
				BU[mapElem[i]] += invA[i][j] * BUloc[j] / jac;
				BV[mapElem[i]] += invA[i][j] * BVloc[j] / jac;
			}
		}
	}
}

void convertTo2D(femTsunami *myTsunami){
	femMesh *theMesh = myTsunami->mesh;
	double R = myTsunami->R;
	double *X = theMesh->X;
	double *Y = theMesh->Y;
	double *Z = theMesh->Z;
	int i;
	for(i=0;i<theMesh->nNode;++i){
		X[i] = 2*R*X[i] / (R+Z[i]);
		Y[i] = 2*R*Y[i] / (R+Z[i]);
	}
	for(i=0;i<theMesh->nNode;++i){
		if(Z[i] < 100){
			Z[i] = 100;
		}
	}
}

////// FEMINTEGRATION //////
femIntegration *femIntegrationCreate(int n, femElementType type)
{
    femIntegration *theRule = malloc(sizeof(femIntegration));
    if (type == FEM_TRIANGLE && n == 3) {
        theRule->n      = 3;
        theRule->xsi    = _gaussTri3Xsi;
        theRule->eta    = _gaussTri3Eta;
        theRule->weight = _gaussTri3Weight; }
    else if (type == FEM_EDGE && n == 2) {
        theRule->n      = 2;
        theRule->xsi    = _gaussEdge2Xsi;
        theRule->eta    = NULL;
        theRule->weight = _gaussEdge2Weight; }

    else Error("Cannot create such an integration rule !");
    return theRule; 
}

void femIntegrationFree(femIntegration *theRule)
{
    free(theRule);
}



////// FEMMESH //////

femMesh *femMeshRead(const char *filename)
{
    femMesh *theMesh = malloc(sizeof(femMesh));

    int i,trash,*elem;
    
    FILE* file = fopen(filename,"r");
    if (file == NULL) Error("No mesh file !");

    if(!fscanf(file, "Number of nodes %d \n", &theMesh->nNode))
		Error("Problem reading mesh ... no node");
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Z = malloc(sizeof(double)*theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        if(!fscanf(file,"%d : %le %le %le \n",&trash,&theMesh->X[i],&theMesh->Y[i],&theMesh->Z[i]))
				Error("Cannot read mesh .... node i");
	 }
    
    char str[256]; 
	 if(!fgets(str, sizeof(str), file))
		  Error("hahahahahahah");;
    if (!strncmp(str,"Number of triangles",19))  { 
        sscanf(str,"Number of triangles %d \n", &theMesh->nElem);
        theMesh->elem = malloc(sizeof(int)*3*theMesh->nElem);
        theMesh->nLocalNode = 3;
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*3]);
            if(!fscanf(file,"%d : %d %d %d\n", &trash,&elem[0],&elem[1],&elem[2]))
					Error("Cannot read mesh .... node i"); }}
    fclose(file);
    return theMesh;
}

void femMeshFree(femMesh *theMesh)
{
    free(theMesh->X);
    free(theMesh->Y);
    free(theMesh->Z);
    free(theMesh->elem);
    free(theMesh);
}

////// FEMEDGES ///////

femEdges *femEdgesCreate(femMesh *theMesh)
{
    femEdges *theEdges = malloc(sizeof(femEdges));
    int nLoc = theMesh->nLocalNode;
    int i,j,n = theMesh->nElem * nLoc;
    femEdge* edges = malloc(n * sizeof(femEdge));
    theEdges->mesh  = theMesh;
    theEdges->edges = edges;
    theEdges->nEdge = n;
    theEdges->nBoundary = n;
    
    for (i = 0; i < theMesh->nElem; i++) {
        int *elem = &(theMesh->elem[i*nLoc]);
        for (j = 0; j < nLoc; j++) {
            int id = i * nLoc + j;
            edges[id].elem[0] = i;
            edges[id].elem[1] = -1;
            edges[id].node[0] = elem[j];
            edges[id].node[1] = elem[(j + 1) % nLoc]; }}

    qsort(theEdges->edges, theEdges->nEdge, sizeof(femEdge), femEdgesCompare);

    int index = 0;          
    int nBoundary = 0;
    
    for (i=0; i < theEdges->nEdge; i++) {
      if (i == theEdges->nEdge - 1 || femEdgesCompare(&edges[i],&edges[i+1]) != 0) {
              edges[index] = edges[i];
              nBoundary++; }
      else {  edges[index] = edges[i];
              edges[index].elem[1] = edges[i+1].elem[0];
              i = i+1;}
      index++; }
      
    theEdges->edges = realloc(edges, index * sizeof(femEdge));
    theEdges->nEdge = index;
    theEdges->nBoundary = nBoundary;
    return theEdges;
}

void femEdgesFree(femEdges *theEdges)
{
    free(theEdges->edges);
    free(theEdges);
}

int femEdgesCompare(const void *edgeOne, const void *edgeTwo)
{
    int *nodeOne = ((femEdge*) edgeOne)->node;
    int *nodeTwo = ((femEdge*) edgeTwo)->node;  
    int  diffMin = fmin(nodeOne[0],nodeOne[1]) - fmin(nodeTwo[0],nodeTwo[1]);
    int  diffMax = fmax(nodeOne[0],nodeOne[1]) - fmax(nodeTwo[0],nodeTwo[1]);
    
    if (diffMin < 0)    return  1;
    if (diffMin > 0)    return -1;
    if (diffMax < 0)    return  1;
    if (diffMax > 0)    return -1; 
                        return  0;
}

void femEdgesMap(femEdges *theEdges, int index, int map[2][2])
{
    int i,j,k;
    int n = theEdges->mesh->nLocalNode;
    
    for (j=0; j < 2; ++j) {
        int node = theEdges->edges[index].node[j];
        for (k=0; k < 2; k++) {
            int elem = theEdges->edges[index].elem[k];
            map[k][j] = (theEdges->mesh->nElem)*n;
            if (elem >= 0) {
                for (i=0; i < n; i++) {
                    if (theEdges->mesh->elem[elem*n + i] == node) {
                        map[k][j] = elem*n + i;  }}}}}
}

////// FEMDISCRETE //////

void _p1c0_x(double *xsi, double *eta) {
    xsi[0] =  0.0;  eta[0] =  0.0;
    xsi[1] =  1.0;  eta[1] =  0.0;
    xsi[2] =  0.0;  eta[2] =  1.0;}

void _p1c0_phi(double xsi, double eta, double *phi){
    phi[0] = 1 - xsi - eta;  
    phi[1] = xsi;
    phi[2] = eta;}

void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta){
    dphidxsi[0] = -1.0;  
    dphidxsi[1] =  1.0;
    dphidxsi[2] =  0.0;
    dphideta[0] = -1.0;  
    dphideta[1] =  0.0;
    dphideta[2] =  1.0;}


void _1c0_x(double *xsi) {
    xsi[0] = -1.0; 
    xsi[1] =  1.0; }


void _1c0_phi(double xsi, double *phi) {   
    phi[0] = (1.0 - xsi)/2.0;
    phi[1] = (1.0 + xsi)/2.0;}

femDiscrete *femDiscreteCreate(int n, femElementType type)
{
    femDiscrete *theSpace = malloc(sizeof(femDiscrete)); 
    if (type == FEM_TRIANGLE && n == 3) {
        theSpace->n       = 3;
        theSpace->order   = 1;
        theSpace->x2      = _p1c0_x;
        theSpace->phi2    = _p1c0_phi;
        theSpace->dphi2dx = _p1c0_dphidx;
        theSpace->x1      = _1c0_x;
        theSpace->phi1    = _1c0_phi;}
    else Error("Cannot create such a discrete space !");
    return theSpace; 
}

void femDiscreteFree(femDiscrete *theSpace)
{
    free(theSpace);
}

void femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta)
{
    mySpace->x2(xsi,eta);
}

void femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi)
{
    mySpace->phi2(xsi,eta,phi);
}

void femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta)
{
    mySpace->dphi2dx(xsi,eta,dphidxsi,dphideta);
}

void femDiscreteXsi1(femDiscrete* mySpace, double *xsi)
{
    mySpace->x1(xsi);
}

void femDiscretePhi1(femDiscrete* mySpace, double xsi, double *phi)
{
    mySpace->phi1(xsi,phi);
}

double femDiscreteInterpolate(double *phi, double *U, int *map, int n)
{
    double u = 0.0; int i;
    for (i=0; i <n; i++)
        u += phi[i]*U[map[i]];
    return u;
}

////// ERROR //////

void femError(char *text, int line, char *file){
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);                                                 
}

void femWarning(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek\n");                                              
}




