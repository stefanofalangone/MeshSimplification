#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include"linearAlgebra_lib.h"
#include"list.h"
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#include <GL/glu.h>
#endif


char PATH_FILE[300];

struct triangle{
  int n;
  int v1, v2, v3;
  int isDeleted;
};
typedef struct triangle Triangle;



struct edge{
  int n;
  int v1, v2;
  double cost;
  int isDeleted;
};
typedef struct edge Edge;

struct vertex{
  int n;
  double x,y,z; //coordinates
  List *triangles; //list of triangles in which this vertex is present
  List *e;
};
typedef struct vertex Vertex;

struct vector{
  double x,y,z; //coordinates
};
typedef struct vector Vector;

struct constraint{
  double A[3][3];
  double b[3];
};
typedef struct constraint Constraint;


/*fv Data contains matrix H, vector c and scalar k used in quadratic function of volume to calculate edge cost */
struct fv_data{
  double **H;
  double c[3];
  double k;
};
typedef struct fv_data Fv_data;



struct mesh{
  int numV, numE, numT;
  int dimE;
  Vertex **v;
  Edge **e;
  Triangle **t;
  Vector **normals; //normals to the faces
  Vector **solutions; //points calculated for each edge that could be substituted
  void (*calculateSolutions)(struct mesh *, int);
};
typedef struct mesh Mesh;


Mesh *mesh;
int numberOfNullSolutions=0;
int cases[4], trianglesDeleted=0, trianglesDeletedWithLength=0;
int STEP_EDGES_SIMPLIFIED=0;
double MAX_DISTANCE_ALLOWED=0.1;
double STOPPING_PERCENTAGE;

double double_rand( double min, double max )
{
    double scale = rand() / (double) RAND_MAX; /* [0, 1.0] */
    return min + scale * ( max - min );      /* [min, max] */
}

Triangle *newTriangle(int n, int a, int b, int c){
  Triangle *t=(Triangle*)malloc(sizeof(Triangle));
  t->n=n; t->v1=a; t->v2=b; t->v3=c; t->isDeleted=0;
  return t;
}


Constraint *newConstraint(){
  Constraint *ret=(Constraint*)malloc(sizeof(Constraint));
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      ret->A[i][j]=0.0;
    }
    ret->b[i]=0.0;
  }
  return ret;
}

Fv_data *newFv(){
  Fv_data *ret=(Fv_data*)malloc(sizeof(Fv_data));
  return ret;
}

Edge *newEdge(int n, int a, int b){
  Edge *e=(Edge*)malloc(sizeof(Edge));
  e->n=n; e->v1=a; e->v2=b; e->cost=99999.0; e->isDeleted=0;
  return e;
}

int insertEdge(Mesh *m, int counter, int v1, int v2){
  List *edges=m->v[v1]->e;
  Edge *currE=NULL;
  int edgeIsPresent=0;
  Edge *e;
  while(edges!=NULL && edgeIsPresent==0){
    currE=(Edge*)edges->value;
    if((currE->v1==v1 && currE->v2==v2) || (currE->v1==v2 && currE->v2==v1)){
      edgeIsPresent=1;
    }
    edges=edges->next;
  }

  if(edgeIsPresent==0){
    if(counter+1==m->dimE){
      //enlarge array for E and for solutions
      int newDim=m->dimE*2;
      Edge **newE=(Edge**)malloc(sizeof(Edge*)*newDim);
      Vector **newSolutions=(Vector**)malloc(sizeof(Vector*)*newDim);
      for(int i=0; i<m->dimE; i++){
        newE[i]=m->e[i];
      }
      free(m->e);
      free(m->solutions);
      m->e=newE;
      m->solutions=newSolutions;
      m->dimE=newDim;
    }
    e=newEdge(counter, v1, v2);
    m->e[counter]=e;
    m->v[ v1 ]->e=addList( m->v[ v1 ]->e , e );
    m->v[ v2 ]->e=addList( m->v[ v2 ]->e , e );
    m->numE++;
    return 1;
  }else{
    e=currE;
    return 0;
  }


}

Vertex *newVertex(int n, double a, double b, double c){
  Vertex *v=(Vertex*)malloc(sizeof(Vertex));
  v->n=n; v->x=a; v->y=b; v->z=c; v->triangles=NULL; v->e=NULL;
  return v;
}

Vector *newVector(double a, double b, double c){
  Vector *v=(Vector*)malloc(sizeof(Vector));
  v->x=a; v->y=b; v->z=c;
  return v;
}







void readVertices(Mesh *m, FILE *fp, int n){
  double x,y,z;
  for(int i=0; i<n; i++){
    fscanf(fp, "%lf %lf %lf", &x, &y, &z);
    m->v[i]=newVertex(i, x, y, z);

    if(strcmp(PATH_FILE, "input/hand.ply")==0) fscanf(fp, "%lf %lf %lf", &x, &y, &z); //this is needed only for v1 v2 v3 n1 n2 n3 file to waste the normals
  }
  printf("read vertices\n");
}

void insertNormal(Mesh *m, int i, int v1, int v2, int v3){
  double x=0, y=0, z=0;
  int v[3]; v[0]=v1; v[1]=v2; v[2]=v3;
  for(int j=0; j<2; j++){ //Newell method to calculate normal for a triangle
    x=x + (m->v[v[j]]->y - m->v[v[j+1]]->y)*(m->v[v[j]]->z + m->v[v[j+1]]->z);
    y=y + (m->v[v[j]]->z - m->v[v[j+1]]->z)*(m->v[v[j]]->x + m->v[v[j+1]]->x);
    z=z + (m->v[v[j]]->x - m->v[v[j+1]]->x)*(m->v[v[j]]->y + m->v[v[j+1]]->y);
  }
  x=x + (m->v[v[2]]->y - m->v[v[0]]->y)*(m->v[v[2]]->z + m->v[v[0]]->z);
  y=y + (m->v[v[2]]->z - m->v[v[0]]->z)*(m->v[v[2]]->x + m->v[v[0]]->x);
  z=z + (m->v[v[2]]->x - m->v[v[0]]->x)*(m->v[v[2]]->y + m->v[v[0]]->y);
  m->normals[i]=newVector(x, y, z);
}

void readTriangles(Mesh *m, FILE *fp, int n){
  int dim, v1, v2, v3, counter=0;;
  for(int i=0; i<n; i++){
    fscanf(fp, "%d %d %d %d", &dim, &v1, &v2, &v3);
    if(dim!=3){
      //ERROR: mesh simplification works with triangular polygons!\n@@@@@@@@@@@@\n");
      return;
    }
    m->t[i]=newTriangle(i, v1, v2, v3);
    m->v[v1]->triangles=addList(m->v[v1]->triangles, m->t[i]);
    m->v[v2]->triangles=addList(m->v[v2]->triangles, m->t[i]);
    m->v[v3]->triangles=addList(m->v[v3]->triangles, m->t[i]);
    insertNormal(m, i, v1, v2, v3);
    counter=counter+insertEdge(m, counter, v1, v2);
    counter=counter+insertEdge(m, counter, v2, v3);
    counter=counter+insertEdge(m, counter, v3, v1);
  }
}


void drawMesh(Mesh *m){
  double col=0.3;
  int trianglesIgnored=0;
  glEnable(GL_NORMALIZE);
  for(int i=0;i<m->numT;i++){
    if( m->t[i]->isDeleted==0 ){
        glBegin(GL_POLYGON);
          glColor3f( 0.1, 0.0+col, 0.0+col );
          col=col+0.1;
          glNormal3f(m->normals[i]->x,m->normals[i]->y, m->normals[i]->z);
          Triangle *curr=m->t[i];
          if(m->v[ curr->v1 ]->z<=-9.9 && m->v[ curr->v1 ]->z>=-10.1) printf("ERROR IN DRAWING at triangle %d %d %d, value of z %lf\n", curr->v1, curr->v2, curr->v3, m->v[ curr->v1 ]->z);
          glVertex3f(m->v[ curr->v1 ]->x, m->v[ curr->v1 ]->y, m->v[ curr->v1 ]->z );
          glVertex3f(m->v[ curr->v2 ]->x, m->v[ curr->v2 ]->y, m->v[ curr->v2 ]->z );
          glVertex3f(m->v[ curr->v3 ]->x, m->v[ curr->v3 ]->y, m->v[ curr->v3 ]->z );
        glEnd();
      }else{
        trianglesIgnored++;
      }
    }
    printf("deleted %d triangles, remaining %d, percentage deleted %lf\n", trianglesIgnored, m->numT - trianglesIgnored, ((double)trianglesIgnored/m->numT)*100 );
}

void swapEdges(Mesh *m, int i, int j){
  Edge *tmpE;
  Vector *tmpV;
  int tmpPos;
  tmpE = m->e[i]; m->e[i] = m->e[j]; m->e[j] = tmpE;
  tmpPos = m->e[i]->n; m->e[i]->n = m->e[j]->n; m->e[j]->n = tmpPos;
  tmpV = m->solutions[i]; m->solutions[i] = m->solutions[j]; m->solutions[j] = tmpV;
}

int leftArray (Mesh *m, int pos){
    if((2*pos)>m->numE) return -1;
    else return (2*pos-1);
}

int rightArray (Mesh *m, int pos){
    if((2*pos+1)>m->numE) return -1;
    else return (2*pos);
}

int fatherArray (Mesh *m, int pos){
    int heapsize=m->numE;
    if((pos)>heapsize || pos==1) return -1;
    else{
        pos=pos-1;
        return (pos/2);
    }
}
void heapify(Mesh *m, int i){ //1 is the head, first element of the heap
    if(i>=m->numE) return;
    int min;
    int newindex;
    int local=i-1;
    int l=leftArray(m, i);
    int r=rightArray(m, i);
    if(l!=-1 && m->e[ local ]->cost > m->e[ l ]->cost){
        min=l;
        newindex=2*i;
    }else min=local;
    if(r!=-1 && m->e[ r ]->cost < m->e[ min ]->cost  ){
        min=r;
        newindex=2*i + 1;
    }
    if(min!=local){
        swapEdges(m, local, min);
        heapify(m, newindex);
    }
}



void quicksort(Mesh *m, int left, int right) {
  if(left >= right) return;
  int i = left, j = right;
  double  pivot = m->e[i]->cost;
  for(;;) {
    while(m->e[i]->cost < pivot) i++;
    while(pivot < m->e[j]->cost) j--;
    if(i >= j) break;
    swapEdges(m, i, j);
    i++; j--;
  }

  quicksort(m, left, i-1);
  quicksort(m, j+1, right);
}
