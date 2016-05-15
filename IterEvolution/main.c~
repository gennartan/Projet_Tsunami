/**
Author : Vincent Legat
Description : Interface graphique pour le projet Tsunami du cours d'élément fini de l'EPL
Modifié par : Antoine Gennart
DATE : 14 mai 2016
**/ 

# include "tsunami.h"

# define GLFW_INCLUDE_GLU
# include <GLFW/glfw3.h>
# define NOPLOT_MODE 1
# define PLOT_MODE   2


int main(int argc, char *argv[])
{   
// Les seuls et uniques paramètres à modifier si on veut changer ce que le programme fait !!!!!
   int theProgramMode = NOPLOT_MODE; 		// PLOT_MODE || NOPLOT_MODE
	int nmax = 25;                 		// nombre d'itération maximale > nsub
	int nsub = 25;                    		// Nombre d'itération entre chaque enregistrement de parametres
	double dt = 5;						  			// timeStep
	const char *name = "tsunamiMedium";		// Nom du fichier pour la lecture





	char meshfile[256];
	char data[256];
	const char *mesh = "%s.txt";
	const char *database = "data/%s-%04.0f";
	sprintf(meshfile,mesh,name);
	sprintf(data,database,name, dt);


	//const char *meshfile = "tsunamiSmall.txt"; // Pour calculer un nouveau tsunami (PLOT_MODE)
	//const char *data = "data/tsunamiSmall"; // Tusnami à lire (NOPLOT_MODE)

//
//  Options possibles
//    -noplot : version non graphique avec creation de fichiers de resultats
//    -plot : lecture des fichiers et animation
//    (par défaut) : version graphique interactive
//

//
//  Options graphiques
//    B = montre la bathymétrie
//    E = montre l'elevation
//    V = zoom sur le japon (3 vues possibles)
//
//

    while(argc--) {
        if (strcmp( *argv, "-noplot" ) == 0)    theProgramMode = NOPLOT_MODE;
        if (strcmp( *argv++, "-plot" ) == 0)    theProgramMode = PLOT_MODE; }

    if (theProgramMode == NOPLOT_MODE) {
        printf("Computing the tsunami : be patient  :-) \n");
        tsunamiCompute(dt,nmax,nsub,meshfile,data);
			theProgramMode = PLOT_MODE;
        }
        
    if (theProgramMode == PLOT_MODE) {
        int nElem,nNode,i,j,index,trash,*elem;
        double *X,*Y,*H,*E,*U,*V;
        double value = 0;
        int width,height;
            
        double t;
        double t0 = 0;
        double R = 6371220;
        double BathMax = 9368;
        GLfloat colors[9], coord[9];
        int theZoom = 0;
        int thePlot = -1;
       
        FILE* file = fopen(meshfile,"r");
        if (file == NULL) {
            printf("Error : cannot open mesh file :-) \n");
            exit(0); }
        if (!fscanf(file, "Number of nodes %d \n", &nNode)) exit(-1);
        X = malloc(sizeof(double)*nNode);
        Y = malloc(sizeof(double)*nNode);
        H = malloc(sizeof(double)*nNode);
        for (i = 0; i < nNode; i++) 
            if (!fscanf(file,"%d : %le %le %le  \n",&trash,&X[i],&Y[i],&H[i])) exit(-1); 
        if (!fscanf(file, "Number of triangles %d \n", &nElem)) exit(-1);  
        elem = malloc(sizeof(int)*3*nElem);
        E = malloc(sizeof(double)*3*nElem);
        U = malloc(sizeof(double)*3*nElem);
        V = malloc(sizeof(double)*3*nElem);
        for (i = 0; i < nElem; i++)     
            if (!fscanf(file,"%d : %d %d %d \n", &trash,&elem[i*3],&elem[i*3+1],&elem[i*3+2])) exit(-1); 
        fclose(file);
        
    
        glfwInit();   	
        GLFWwindow* window = glfwCreateWindow(480,480,"MECA1120 : Tsunami",NULL,NULL);    
        glfwMakeContextCurrent(window);
        glShadeModel(GL_SMOOTH);      
        glfwSwapInterval(1);
    
        GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 0.0 };
        GLfloat mat_shininess[] = { 50.0 };
        GLfloat light_position[] = { 8.0, 8.0, 8.0, 0.0 };
        
        glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
        glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
        glLightfv(GL_LIGHT0, GL_POSITION, light_position);
        GLfloat light_radiance[] = {1.0, 1.0, 1.0, 1.0};
    
        glLightfv(GL_LIGHT0, GL_DIFFUSE, light_radiance);
        glLightfv(GL_LIGHT0, GL_SPECULAR, light_radiance);
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        glDepthFunc(GL_LEQUAL);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_COLOR_MATERIAL);
        glEnable(GL_NORMALIZE);	
    
        int frame=0;
        do {
            t = glfwGetTime();                  
            frame++;
				frame = frame%(nmax/nsub+1);
				
    
            if (thePlot == 1) {            
                char filename[256];
                const char *basename = "%s-%08d.txt";
                sprintf(filename, basename, data, frame * nsub);    
          	    if (access(filename, F_OK)) {
                    t0 = t; 
							}   
                else {
                    printf("===  Reading local file %s %d %f \n",filename,frame,t);
                    tsunamiReadFile(data,frame*nsub,U,V,E,nElem); }}
       
            glfwGetFramebufferSize(window,&width,&height);
            height = height > 0 ? height : 1;
            glViewport( 0, 0, width, height );

            glClearColor( 0.9f, 0.9f, 0.8f, 0.0f );
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            gluPerspective(65.0f,(GLfloat)width/(GLfloat)height,1.0f,100.0f);

            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            gluLookAt(0.0f,1.0f,0.0f,0.0f, 20.0f, 0.0f,0.0f,0.0f,1.0f);  
            if (theZoom == 0)  glTranslatef(0.0f,14.0f,0.0f);     // Vue terre complete
            if (theZoom == 1)  glTranslatef(0.0f,9.0f,-3.0f);     // Vue zoom moyen
            if (theZoom == 2)  glTranslatef(1.0f,7.0f,-4.0f);     // Vue zoom rapproche
            
            double alpha = 389;
            double beta = 0;
            glRotatef(0.3f*(GLfloat)alpha + (GLfloat)beta*10.0f,0.0f,0.0f,1.0f);
            
            // A commenter pour supprimer la sphere interieure
            // Conseille pour les maillages grossiers :-)
            
            GLUquadricObj *quadratic = gluNewQuadric();         
            gluQuadricNormals(quadratic, GLU_SMOOTH); 
            glColor3f(1.0,1.0,1.0);
            gluSphere(quadratic,5.95,400,200);
            
     
            for (i=0; i < nElem; ++i) {
                for (j=0; j < 3; ++j) {
                    index = elem[3*i+j];
                    if (thePlot == -1)  value = H[index]/BathMax;
                    if (thePlot ==  1)  value = E[3*i+j]*10;
                    if (value < 0) value = 0;
                    if (value > 1) value = 1; 
                    colors[j*3+0] = 3.5*(value)*(value);
                    colors[j*3+1] = (1-value)*(value)*3.5;
                    colors[j*3+2] = (1-value)*(1-value);
                    double x = X[index]; 
                    double y = Y[index]; 
                    double Factor = (4*R*R + x*x + y*y)*(R/6);
                    coord[j*3+0] = 4*R*R * x / Factor;
                    coord[j*3+1] = 4*R*R * y / Factor;
                    coord[j*3+2] = (4*R*R - x*x - y*y)*R / Factor;  } 
      
                glEnableClientState(GL_VERTEX_ARRAY);
                glEnableClientState(GL_COLOR_ARRAY);
                glEnableClientState(GL_NORMAL_ARRAY);
                glVertexPointer(3, GL_FLOAT, 0, coord);
                glNormalPointer(GL_FLOAT, 0, coord);
                glColorPointer(3, GL_FLOAT, 0, colors);
                glDrawArrays(GL_TRIANGLES, 0, 3);
                glDisableClientState(GL_NORMAL_ARRAY);    
                glDisableClientState(GL_COLOR_ARRAY);
                glDisableClientState(GL_VERTEX_ARRAY);      
                                 
                glColor3f(0.0, 0.0, 0.0);
                glEnableClientState(GL_VERTEX_ARRAY);
                for (j=0; j < 9; ++j)
                     coord[j] = coord[j] * 1.001;
                glVertexPointer(3, GL_FLOAT, 0, coord);
                glDrawArrays(GL_LINE_LOOP, 0, 3);
                glDisableClientState(GL_VERTEX_ARRAY); }
                
            if (glfwGetKey(window,'V') == GLFW_PRESS)   theZoom = (theZoom + 1) % 3; 
            if (glfwGetKey(window,'E') == GLFW_PRESS)   {thePlot = 1; t0 = t; 
				}
            if (glfwGetKey(window,'B') == GLFW_PRESS)   thePlot = -1;
 
            glfwSwapBuffers(window);
            glfwPollEvents();

        
        }  while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
    
        glfwTerminate();
			free(E); free(U); free(V); free(H); free(X); free(Y); free(elem);
        exit( EXIT_SUCCESS ); }
}

