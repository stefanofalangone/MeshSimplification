#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#include <GL/glu.h>
#endif
#include<math.h>
#include <stdbool.h>
#include<stdio.h>
#include <stdlib.h>
#include <time.h>
#include<string.h>
#include"lib/meshSimplification.h"


double r=0.5;
double angleInRadiants=0.0;
double j=1.0, i=1.0, k1=1.0, k2=5.0;
int counter=0;
double colorStep=0.0;
double ratio=1.0;
double angle=0.0;
double angleHours=0.0;
int timeMinutes=0;
double xLookingPoint=0.0;
double zLookingPoint= 10.0;
double yLookingPoint=0.0;

GLvoid drawScene( GLvoid ){
  glMatrixMode(GL_PROJECTION);
  //glClearColor( 1.0, 1.0, 1.0, 1.0 );
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glEnable( GL_DEPTH_TEST );
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glShadeModel( GL_FLAT );
  if(counter==0){
    gluPerspective(45, 1, 0.1, 15.0);
    counter++;
  }

  glMatrixMode(GL_MODELVIEW);
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glPushMatrix();
    gluLookAt( zLookingPoint*sin(angle), yLookingPoint, zLookingPoint*cos(angle), 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 );
    glPushMatrix();
      drawMesh(mesh);
    glPopMatrix();
  glPopMatrix();
  glFlush();
  glLoadIdentity();
//glDisable(GL_DEPTH_TEST);
}


double toRadians (double angle) {
  return angle * (M_PI / 180);
}

void keyboard(unsigned char key, int x, int y){
switch (key) {
  case 'a':
    angle=angle-0.1;
    glutPostRedisplay();
    break;
  case 'd':
    angle=angle+0.1;
    glutPostRedisplay();
    break;
  case 'w':
    yLookingPoint+=0.1;
    glutPostRedisplay();
    //}
    break;
  case 's':
      yLookingPoint-=0.1;
      glutPostRedisplay();
    break;
  case 'p':
      zLookingPoint-=0.5;
      glutPostRedisplay();
    break;
  case 'l':
      zLookingPoint+=0.4;
      glutPostRedisplay();
    break;
  case 'z':
    zLookingPoint=zLookingPoint*(-1);
    glutPostRedisplay();
    break;
  case '-':
    zLookingPoint-=-1;
    glutPostRedisplay();
    break;
  case '+':
    zLookingPoint-=1;
    glutPostRedisplay();
    break;
  case 'k':
    simplification(mesh, STEP_EDGES_SIMPLIFIED);
    glutPostRedisplay();
    break;
  case 27: /* ESC */
  writeOutput(mesh);
  exit(0);
  }
}


void testMath(){
    double a[3]={1,2,3};
    double b[3]={1,2,3};
    printf("%f\n", dotProduct(a,b,3));
    double *vet1=sumVectors(a,b,3);
    printf("vectors summed: %f %f %f\n", vet1[0], vet1[1], vet1[2]);
    double **A=(double**)malloc(sizeof(double*)*3);
    for(int i=0;i<3;i++) A[i]=(double*)malloc(sizeof(double)*3);
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) A[i][j]=1.0;
    }
    A[2][0]=5.0; A[0][2]=2.0;
    printf("det should be -4 and is %f\n",det3x3(A));
    double v1[]={1, 2, 1};
    double v2[]={1, 0, 1};
    double *crossP=crossProduct3dim(v1,v2);
    printf("cross product should be 2 0 -2, is: %f %f %f\n", crossP[0], crossP[1], crossP[2]);
    b[2]=2.5;
    vet1=solveLinearSystem(A, b, 3);
    printf("solutions to system are:\n %f %f %f (should be 0.125, 2.8, -1)\n",vet1[0], vet1[1], vet1[2]);
    A[0][0]=1; A[0][1]=1;
    A[1][0]=0; A[1][1]=-1;
    b[0]=2; b[1]=1;
    vet1=solveLinearSystem(A, b, 2);
    printf("solutions to system are:\n %f %f (should be 0.125, 2.8, -1)\n",vet1[0], vet1[1]);
    A=allocMatrix(3, 3);
    A[0][0]=1; A[0][1]=1;
    A[1][0]=0; A[1][1]=-1;
    A[2][2]=-2.2;
    vet1=matvet(A, 3, 3, b, 3);
    printf("solutions to matvet are:\n %f %f %f (should be 3.000000 -1.000000 -5.500000)\n",vet1[0], vet1[1], vet1[2]);
    A=inverseMatrix3x3(A);
    printf("inverse matrix\n");
    printMatrix(A, 3, 3);

    A=allocMatrix(2, 3);
    A[0][0]=1; A[0][1]=1;
    A[1][0]=0; A[1][1]=-1; A[1][2]=-0.15;
    double **Z=allocMatrix(3, 1);
    Z[0][0]=Z[1][0]=-3;
    A=matmat(A, 2, 3, Z, 3, 1);
    printf("matrix product is: (should be -6 3)\n ");
    printMatrix(A, 2, 1);
    A=allocMatrix(4,4);
    for(int i=0; i<4; i++){
      for(int j=0; j<4; j++) A[i][j]=double_rand(-10.0, 10.0);
    }
    printf("starting 4x4 inverse testing\n");
    printMatrix(A, 4, 4);
    A=inverseMatrix4x4(A);
    printMatrix(A, 4,4);
}

void testMeshBuilding(Mesh *m){
  printf("number of vertices=%d, triangles=%d\n", m->numV, m->numT);
  printf("vertices are:\n");
  for(int i=0;i<m->numV;i++){
    printf("n: %d, (%f, %f, %f)\n", m->v[i]->n, m->v[i]->x, m->v[i]->y, m->v[i]->z);
    List *curr=m->v[i]->triangles;
    printf("Triangles\n");
    while(curr!=NULL){
      printf("-> %d", ((Triangle*)curr->value)->n);
      curr=curr->next;
    }
    printf("\n");
  }
}

void testMeshSimplification(Mesh *m){
  for(int i=0; i<m->numT-1; i++){
    Triangle *t=m->t[i];
    if(t->v1==t->v2 || t->v2==t->v3 || t->v1==t->v3) printf("FATAL ERROR IN BUILDING MESH\n");
  }
  for(int i=0; i<m->numE; i++) m->calculateSolutions(m, i);
  printf("considered %d arcs, found null solution for %d\n", m->numE, numberOfNullSolutions);
  printf("cases before third constraint with: 0 previous constraints %d, 1 previous constr %d, 2 previous constr %d\n", cases[0], cases[1], cases[2]);
  quicksort(m, 0, m->numE-1);
  //simplification(m, 6*STEP_EDGES_SIMPLIFIED);
  printf("deleted, number of triangles was %d\n", m->numT);
}

//first argument path of file, second method chosen, third error distance, fourth number STEP_EDGES_SIMPLIFIED (0 means set it to default),
int main( int argc, char** argv ){
  srand((unsigned int)time(NULL));
  //testMath();
  strcpy(PATH_FILE, "input/bunny_simple2.off");
  //strcpy(PATH_FILE, "input/bunny_simple2.off.out");
  strcpy(PATH_FILE, "input/cow.off");
  //strcpy(PATH_FILE, "input/cow.off.out");
  //strcpy(PATH_FILE, "input/homer.off");
  //strcpy(PATH_FILE, "input/torus_simple.off");
  //strcpy(PATH_FILE, "input/bunny.off");
  //strcpy(PATH_FILE, "input/hand.ply");
  strcpy(PATH_FILE, "input/dragon.ply");
  int methodChosen=atoi(argv[2]); //0 is calculateSolutionsLength, 1 is calculateSolutionsLindstrom //2 is Garland
  mesh=initialize(methodChosen);
  //testMeshBuilding(mesh);
  MAX_DISTANCE_ALLOWED=atof(argv[3]);
  if(atof(argv[4])>1) STEP_EDGES_SIMPLIFIED=atof(argv[4]);
  STOPPING_PERCENTAGE=5.0;
  printf("before simplification\n");
  testMeshSimplification(mesh);
  glutInit(&argc, argv);
  glutInitDisplayMode( GLUT_RGBA | GLUT_DEPTH );
  glutInitWindowPosition( 300, 100 );
  glutInitWindowSize ( 600, 600 );
  glutCreateWindow ( argv[0] );
  glutDisplayFunc( drawScene );
  glutKeyboardFunc( keyboard );
  glutMainLoop ();
  return 0;
}
