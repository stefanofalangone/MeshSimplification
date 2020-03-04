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




int checkConstrains(Mesh *m, Constraint *c, int n_constraints){ //1 means accepted constraint, 0 refuted
  double a1[3], a2[3], a3[3];
  double left, right, **matrix;
  switch(n_constraints){
    case 1:
      if(c->A[0][0]<=0.000001 && c->A[0][1]<=0.000001 && c->A[0][2]<=0.000001 ) return 0;
      else return 1;
    case 2:
      a1[0]=c->A[0][0]; a1[1]=c->A[0][1]; a1[2]=c->A[0][2];
      a2[0]=c->A[1][0]; a2[1]=c->A[1][1]; a2[2]=c->A[1][2];
      /*left=dotProduct(a1, a2, 3);
      left=left*left;
      right=norm(a1, 3)*norm(a2, 3)*cosTwoVectors(a1,a2,3);
      right=right*right;
      return (left<right);*/
      //now my version of this constraint: (it seems to work better)
      double cosine=cosTwoVectors(a1,a2,3);
      return (cosine*cosine < 0.97);
    case 3:
    a1[0]=c->A[0][0]; a1[1]=c->A[0][1]; a1[2]=c->A[0][2];
    a2[0]=c->A[1][0]; a2[1]=c->A[1][1]; a2[2]=c->A[1][2];
    a3[0]=c->A[2][0]; a3[1]=c->A[2][1]; a3[2]=c->A[2][2];
    double *crossProduct=crossProduct3dim(a1, a2);
    left=dotProduct( crossProduct , a3, 3 );
    left=left*left;
    right=norm( crossProduct , 3 ) * norm(a3, 3) * sinTwoVectors3dim(crossProduct, a3);
    right=right*right;
    free(crossProduct);
    return (left > right);
    //now my version of constraint
    /*matrix=allocMatrix(3,3);
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++) matrix[i][j]=c->A[i][j];
    }
    double det=det3x3(matrix);
    deallocMatrix(matrix, 3, 3);
    if(fabs(det)>=0.00001) return 1;
    else return 0; */
  }
  return 0;
}

int isBoundaryEdge(Mesh *m, Edge *e){
  int v1=e->v1;
  int v2=e->v2;
  List *intersection=intersectList(m->v[ v1 ]->triangles, m->v[ v2 ]->triangles );
  int sizeOfList=sizeList(intersection);
  deallocList(intersection);
  return (sizeOfList==1);
  return 0;
}

int edgesAreEqual(Edge *e1, Edge *e2){
  int n=0;
  if(e1->v1==e2->v1 || e1->v1==e2->v2) n++;
  if(e1->v2==e2->v1 || e1->v2==e2->v2) n++;
  if(n==2) return 1;
  else return 0;
}

//takes a  triangle, does the external product of v0 x v1, v1 x v2, v2 x v0 and sums to newC
//furthermore t=t + discriminant of matrix [v0, v1, v2]
void sumTriangleComponents(Mesh *m, double *newC, double *t,  Triangle *t_i){
  double a[3], b[3], c[3];
  double **M=allocMatrix(3,3); //auxiliary matrix to call determinant
  a[0]=m->v[ t_i->v1 ]->x; a[1]=m->v[ t_i->v1 ]->y; a[2]=m->v[ t_i->v1 ]->z;
  b[0]=m->v[ t_i->v2 ]->x; b[1]=m->v[ t_i->v2 ]->y; b[2]=m->v[ t_i->v2 ]->z;
  c[0]=m->v[ t_i->v3 ]->x; c[1]=m->v[ t_i->v3 ]->y; c[2]=m->v[ t_i->v3 ]->z;
  double **vertices=allocMatrix(3,3);
  vertices[0][0]=m->v[ t_i->v1 ]->x; vertices[0][1]=m->v[ t_i->v1 ]->y; vertices[0][2]=m->v[ t_i->v1 ]->z;
  vertices[1][0]=m->v[ t_i->v2 ]->x; vertices[1][1]=m->v[ t_i->v2 ]->y; vertices[1][2]=m->v[ t_i->v2 ]->z;
  vertices[2][0]=m->v[ t_i->v3 ]->x; vertices[2][1]=m->v[ t_i->v3 ]->y; vertices[2][2]=m->v[ t_i->v3 ]->z;
  double *normal=normalToTriangle(vertices);

  double a1=normal[0], b1=normal[1], c1=normal[2];
  for(int i=0;i < 3; i++){
    M[0][i]=a[i]; M[1][i]=b[i]; M[2][i]=c[i];
  }
  /*Transposing doesn't change the determinant, but to be faithful towards alghoritm described in the paper this characteristic is preserved */
  double *firstPoint=allocVector(3);
  for(int i=0; i<3; i++) firstPoint[i]=vertices[0][i];
  double d=-1*dotProduct(normal, firstPoint, 3);
  newC[0]+=a1; newC[1]+=b1; newC[2]+=c1;
  double **transposed=transpose(M, 3, 3);
  *t=*t-det3x3(transposed);
  //*t=*t+d;
  deallocMatrix(M, 3, 3); deallocMatrix(transposed, 3, 3); deallocMatrix(vertices, 3,3); free(firstPoint); free(normal);
}

double *normalToT(Mesh *m, Triangle *t_i){
  double a[3], b[3], c[3];
  double *newC=allocVector(3);
  a[0]=m->v[ t_i->v1 ]->x; a[1]=m->v[ t_i->v1 ]->y; a[2]=m->v[ t_i->v1 ]->z;
  b[0]=m->v[ t_i->v2 ]->x; b[1]=m->v[ t_i->v2 ]->y; b[2]=m->v[ t_i->v2 ]->z;
  c[0]=m->v[ t_i->v3 ]->x; c[1]=m->v[ t_i->v3 ]->y; c[2]=m->v[ t_i->v3 ]->z;
  double *v0x1=crossProduct3dim( a, b );
  double *v1x2=crossProduct3dim( b, c );
  double *v2x0=crossProduct3dim( c, a );
  newC[0]=newC[0]+v0x1[0]+v1x2[0]+v2x0[0];
  newC[1]=newC[1]+v0x1[1]+v1x2[1]+v2x0[1];
  newC[2]=newC[2]+v0x1[2]+v1x2[2]+v2x0[2];
  free(v0x1); free(v1x2); free(v2x0);
  return newC;
}

int addConstraintIfIndependent(Mesh *m, Constraint *constr, int n_constraints_accepted, double *newC, double b){
  for(int i=0; i<3 ; i++) constr->A[n_constraints_accepted][i]=newC[i];
  constr->b[n_constraints_accepted]=b;
  if(checkConstrains(m, constr, n_constraints_accepted+1)==1){
    n_constraints_accepted++; //new constraint is setted
  } else {
    for(int i=0; i<3 ; i++) constr->A[n_constraints_accepted][i]=0.0;
    constr->b[n_constraints_accepted]=0.0;
  }
  return n_constraints_accepted;
}

void printConstr(Constraint *constr){
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++) printf("%lf ", constr->A[i][j]);
    printf("\n");
  }
  printf("b: ");
  for(int i=0; i<3; i++) printf("%lf ", constr->b[i]);
  printf("\n");
}

int quadraticOpt(Mesh *m, Constraint *constr, int n_constraints_accepted, double **A){ //expects A to be 4 x 4
  double **Q=allocMatrix(3 - n_constraints_accepted, 3);
  /*Build Q, three rows orthogonal to each other and to vectors already present. If n_constraints_accepted=0, then Q is already the null matrix and shall have 1 on the diagonal */
    if(n_constraints_accepted==0){
      for(int i=0; i<3; i++) Q[i][i]=100000.0;
    }
    double *curr=allocVector(3);
    for(int i=0; i<3; i++) curr[i]=constr->A[0][i];
    double multiply=10000000.0;
    if(n_constraints_accepted==1){
      free(Q[0]); free(Q[1]);
      Q[0]=orthogonalVectorToVector(curr);
      Q[1]=orthogonalVectorTo2Vectors(curr, Q[0]);
      for(int j=0; j<3; j++) Q[0][j]*=multiply;
      for(int j=0; j<3; j++) Q[1][j]*=multiply;
    }
    if(n_constraints_accepted==2){
       double *curr2=allocVector(3);
       for(int i=0; i<3; i++){
         curr[i]=constr->A[0][i]; curr2[i]=constr->A[1][i];
       }
       free(Q[0]);
       Q[0]=orthogonalVectorTo2Vectors(curr, curr2);
       for(int j=0; j<3; j++) Q[0][j]*=multiply;
       free(curr2);
    }

    double **reducedA=allocMatrix(3,3);
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++) reducedA[i][j]=A[i][j];
    }
    double *b=allocVector(3);
    b[0]=A[0][3]; b[1]=A[1][3]; b[2]=A[2][3];
    double *Qb=matvet(Q, 3 - n_constraints_accepted, 3, b, 3); //shall change sign when considered as term noted
    double **QA=matmat(Q, 3 - n_constraints_accepted, 3, reducedA, 3, 3);
    int n_iterations= 3 - n_constraints_accepted;
    for(int i=0; i<n_iterations; i++){
      n_constraints_accepted=addConstraintIfIndependent(m, constr, n_constraints_accepted, QA[i], Qb[i]);
    }

  deallocMatrix(Q, 3-n_constraints_accepted, 3); deallocMatrix(reducedA, 3, 3); free(curr); free(Qb); deallocMatrix(QA, 3 - n_constraints_accepted, 3); free(b);
  return n_constraints_accepted;
}

void calculateBoundaryEdge(double **H, double *e1, double *e2){ //H is 4x4, e1 e2 have dim=3
  double **x_e1=allocMatrix(4,4), **e1_x;
  x_e1[0][1]=-e1[2]; x_e1[1][0]=e1[2];
  x_e1[0][2]=e1[1]; x_e1[2][0]=-e1[1];
  x_e1[1][2]=-e1[0]; x_e1[2][1]=e1[0];
  e1_x=transpose(x_e1, 3, 3);
  double **p=matmat(x_e1, 3, 3, e1_x, 3, 3);
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++) H[i][j]+=p[i][j];
  }
  double *crossP=crossProduct3dim(e1, e2);
  for(int j=0; j<3; j++) H[j][3]+=crossP[j];
  free(crossP);
  H[3][3]=0.5*dotProduct(e2, e2, 3);
  deallocMatrix(x_e1, 3, 3); deallocMatrix(e1_x, 3, 3); deallocMatrix(p, 3, 3);
}
// n_constraints_accepted is the number of constraints checked and setted
int getAllConstraints(Mesh *m, Constraint *constr, Fv_data *fv, Fv_data *fv2, int n_constraints_accepted, int edge){
  double *newC=malloc(sizeof(double)*3), *t=(double*)malloc(sizeof(double));
  Vertex *a1, *b1;
  newC[0]=0.0; newC[1]=0.0; newC[2]=0.0;
  *t=0.0;
  int v1=m->e[edge]->v1;
  int v2=m->e[edge]->v2;
  int i, j, counter=0;
  List *l1, *l2, *curr;
  double a[3], b[3], c[3];
  double *e1=malloc(sizeof(double)*3), *e2=malloc(sizeof(double)*3), *e3=malloc(sizeof(double)*3);
  double *v1_e1=malloc(sizeof(double)*3);
  double *v2_e1=malloc(sizeof(double)*3);
  double *crossP;
  Edge *currE;
  Triangle *t_i;

    double **vertices, **transposedVertices, **H, det, k;
  switch(n_constraints_accepted){
    case 0: //Constraint VOLUME PRESERVATION
      l1=m->v[v1]->triangles; l2=m->v[v2]->triangles;
      curr=l1;
      while(curr!=NULL){
        t_i=(Triangle*)curr->value;
        curr=curr->next;
        sumTriangleComponents(m, newC, t,  t_i);
      }
      curr=l2;
      while(curr!=NULL){
        if(contains(l1, curr->value)==0){
          t_i=curr->value;
          sumTriangleComponents(m, newC, t,  t_i);
        }
        curr=curr->next;
      }
      for(int i=0; i<3; i++){
        newC[i]=newC[i]/6;
      }
      *t=1*(*t)/6;
      n_constraints_accepted=addConstraintIfIndependent(m, constr, n_constraints_accepted, newC, *t);
      *t=0.0;

    case 1: //BOUNDARY PRESERVATION
      for(int i=0; i<3; i++){
        e1[i]=0.0; e2[i]=0.0; e3[i]=0.0;
      }
      l1=m->v[ v1 ]->e;
      l2=m->v[ v2 ]->e;
      while(l2!=NULL){
        if(isBoundaryEdge(m, (Edge*)l2->value) && l2->value!=m->e[ edge ]){
          counter++;
          currE=(Edge*)l2->value;
            v1_e1[0]=m->v[ currE->v1 ]->x; v1_e1[1]=m->v[ currE->v1 ]->y; v1_e1[2]=m->v[ currE->v1 ]->z;
            v2_e1[0]=m->v[ currE->v2 ]->x; v2_e1[1]=m->v[ currE->v2 ]->y; v2_e1[2]=m->v[ currE->v2 ]->z;
          for(int i=0; i<3; i++){
            e1[i]=e1[i] + v2_e1[i] - v1_e1[i];
          }
          crossP=crossProduct3dim(v2_e1, v1_e1);
          for(int i=0; i<3; i++) e2[i]= e2[i] + crossP[i];
          free(crossP);
        }
        l2=l2->next;
      }

      while(l1!=NULL){
        if(isBoundaryEdge(m, (Edge*)l1->value) && l1->value!=m->e[ edge ]){
          counter++;
          currE=(Edge*)l1->value;
            v1_e1[0]=m->v[ currE->v1 ]->x; v1_e1[1]=m->v[ currE->v1 ]->y; v1_e1[2]=m->v[ currE->v1 ]->z;
            v2_e1[0]=m->v[ currE->v2 ]->x; v2_e1[1]=m->v[ currE->v2 ]->y; v2_e1[2]=m->v[ currE->v2 ]->z;
          //  ////printf("here\n");
          for(int i=0; i<3; i++){
            e1[i]=e1[i] + v2_e1[i] - v1_e1[i];
          }
          crossP=crossProduct3dim(v2_e1, v1_e1);
          for(int i=0; i<3; i++) e2[i]= e2[i] + crossP[i];
          free(crossP);
        }
        l1=l1->next;
      }
      if(counter>0){
        v1_e1[0]=m->v[ m->e[edge]->v1 ]->x; v1_e1[1]=m->v[ m->e[edge]->v1 ]->y; v1_e1[2]=m->v[ m->e[edge]->v1 ]->z;
        v2_e1[0]=m->v[ m->e[edge]->v2 ]->x; v2_e1[1]=m->v[ m->e[edge]->v2 ]->y; v2_e1[2]=m->v[ m->e[edge]->v2 ]->z;
        for(int i=0; i<3; i++) e1[i]=e1[i] + v2_e1[i] - v1_e1[i];
        crossP=crossProduct3dim(v2_e1, v1_e1);
        for(int i=0; i<3; i++) e2[i]= e2[i] + crossP[i];
        free(crossP);
        e3=crossProduct3dim(e1, e2);
        newC=scalarVector( dotProduct(e1, e1, 3), e3, 3 );
        n_constraints_accepted=addConstraintIfIndependent(m, constr, n_constraints_accepted, newC, -1*dotProduct(e3, e3, 3));
        newC=crossProduct3dim(e1, e3);
        n_constraints_accepted=addConstraintIfIndependent(m, constr, n_constraints_accepted, newC, 0.0);
      }
    case 2: //volume optimization
      l1=m->v[v1]->triangles; l2=m->v[v2]->triangles;
      curr=l1;
      vertices=allocMatrix(3,3);
      H=allocMatrix(4,4);
      k=0;
      int counter=0;
      double a, b, c, d, **transposed, *firstPoint=allocVector(3);
      Triangle *t_i;
      while(curr!=NULL){
          t_i=(Triangle*)curr->value;
          curr=curr->next;
          vertices[0][0]=m->v[ t_i->v1 ]->x; vertices[0][1]=m->v[ t_i->v1 ]->y; vertices[0][2]=m->v[ t_i->v1 ]->z;
          vertices[1][0]=m->v[ t_i->v2 ]->x; vertices[1][1]=m->v[ t_i->v2 ]->y; vertices[1][2]=m->v[ t_i->v2 ]->z;
          vertices[2][0]=m->v[ t_i->v3 ]->x; vertices[2][1]=m->v[ t_i->v3 ]->y; vertices[2][2]=m->v[ t_i->v3 ]->z;
          double *normal=normalToTriangle(vertices);
          a=normal[0]; b=normal[1], c=normal[2];
          transposed=transpose(vertices, 3, 3);
          d=-1*det3x3(transposed);

          for(int i=0; i<3; i++) firstPoint[i]=vertices[0][i];
          a=normal[0]; b=normal[1], c=normal[2];
          H[0][0]+=a*a; H[0][1]+=a*b; H[0][2]+=a*c; H[0][3]+=a*d;
          H[1][0]+=b*a; H[1][1]+=b*b; H[1][2]+=b*c; H[1][3]+=b*d;
          H[2][0]+=c*a; H[2][1]+=c*b; H[2][2]+=c*c; H[2][3]+=c*d;
          H[3][0]+=d*a; H[3][1]+=d*b; H[3][2]+=d*c; H[3][3]+=d*d;
          free(normal);
          deallocMatrix(transposed, 3, 3);

      }
      curr=l2;
      while(curr!=NULL){
        t_i=(Triangle*)curr->value;
        if(contains(l1, curr->value)==0){
          vertices[0][0]=m->v[ t_i->v1 ]->x; vertices[0][1]=m->v[ t_i->v1 ]->y; vertices[0][2]=m->v[ t_i->v1 ]->z;
          vertices[1][0]=m->v[ t_i->v2 ]->x; vertices[1][1]=m->v[ t_i->v2 ]->y; vertices[1][2]=m->v[ t_i->v2 ]->z;
          vertices[2][0]=m->v[ t_i->v3 ]->x; vertices[2][1]=m->v[ t_i->v3 ]->y; vertices[2][2]=m->v[ t_i->v3 ]->z;
          double *normal=normalToTriangle(vertices);
          a=normal[0]; b=normal[1], c=normal[2];
          transposed=transpose(vertices, 3, 3);
          d=-1*det3x3(transposed);


          for(int i=0; i<3; i++) firstPoint[i]=vertices[0][i];
          a=normal[0]; b=normal[1], c=normal[2];
          H[0][0]+=a*a; H[0][1]+=a*b; H[0][2]+=a*c; H[0][3]+=a*d;
          H[1][0]+=b*a; H[1][1]+=b*b; H[1][2]+=b*c; H[1][3]+=b*d;
          H[2][0]+=c*a; H[2][1]+=c*b; H[2][2]+=c*c; H[2][3]+=c*d;
          H[3][0]+=d*a; H[3][1]+=d*b; H[3][2]+=d*c; H[3][3]+=d*d;
          free(normal);
          deallocMatrix(transposed, 3, 3);
        }
        curr=curr->next;
      }

      double multiply=1;
      for(int i=0; i<4; i++){
        for(int j=0; j<4; j++) H[i][j]=fv->H[i][j]=H[i][j]/18;
      }
      n_constraints_accepted=quadraticOpt(m, constr, n_constraints_accepted, H);
      cases[n_constraints_accepted]++;
      deallocMatrix(H, 4, 4); deallocMatrix(vertices, 3, 3); free(firstPoint);
      if(n_constraints_accepted==3) break;

    case 3: //BOUNDARY OPTIMIZATION
      for(int i=0; i<3; i++){
        e1[i]=0.0; e2[i]=0.0; e3[i]=0.0;
      }
      l1=m->v[ v1 ]->e;
      l2=m->v[ v2 ]->e;
      H=allocMatrix(4,4);
      int isFirst=1;
      while(l2!=NULL){
        if(isBoundaryEdge(m, (Edge*)l2->value) && l2->value!=m->e[ edge ]){
          for(int i=0; i<3; i++) e1[i]=0.0;
          for(int i=0; i<3; i++) e2[i]=0.0;
          counter++;
          currE=(Edge*)l2->value;
            v1_e1[0]=m->v[ currE->v1 ]->x; v1_e1[1]=m->v[ currE->v1 ]->y; v1_e1[2]=m->v[ currE->v1 ]->z;
            v2_e1[0]=m->v[ currE->v2 ]->x; v2_e1[1]=m->v[ currE->v2 ]->y; v2_e1[2]=m->v[ currE->v2 ]->z;

          for(int i=0; i<3; i++) e1[i] = v2_e1[i] - v1_e1[i];
          crossP=crossProduct3dim(v2_e1, v1_e1);
          for(int i=0; i<3; i++) e2[i] =  crossP[i];
          free(crossP);
          calculateBoundaryEdge(H, e1, e2);
        }
        l2=l2->next;
        if(l2==NULL && isFirst==1){
          isFirst=0;
          l2=l1;
        }
      }


      if(counter>0){
          v1_e1[0]=m->v[ m->e[edge]->v1 ]->x; v1_e1[1]=m->v[ m->e[edge]->v1 ]->y; v1_e1[2]=m->v[ m->e[edge]->v1 ]->z;
          v2_e1[0]=m->v[ m->e[edge]->v2 ]->x; v2_e1[1]=m->v[ m->e[edge]->v2 ]->y; v2_e1[2]=m->v[ m->e[edge]->v2 ]->z;
          for(int i=0; i<3; i++) e1[i] = v2_e1[i] - v1_e1[i];
          crossP=crossProduct3dim(v2_e1, v1_e1);
          for(int i=0; i<3; i++) e2[i] = crossP[i];
          free(crossP);
          calculateBoundaryEdge(H, e1, e2);
          for(int i=0; i<4; i++){
            for(int j=0; j<4; j++) H[i][j]=fv2->H[i][j]=H[i][j]/2;
          }
          n_constraints_accepted=quadraticOpt(m, constr, n_constraints_accepted, H);
          cases[n_constraints_accepted]++;
      }
      deallocMatrix(H, 4, 4);
      if(n_constraints_accepted==3) break;

    case 4:
    //TRIANGLE SHAPE
      a1=m->v[ v1 ];
      b1=m->v[ v2 ];
      Edge *currE;
      Vertex *currV;
      int firstIteration=1;
      double newc[3];
      for(int i=0; i<3; i++) newc[i]=fv->c[i]=0;
      H=allocMatrix(4,4); k=0;
      curr=a1->e;
      counter=0;
      while(curr!=NULL){
        currE=(Edge*)curr->value;
        curr=curr->next;
        if(firstIteration==1 && curr==NULL){
          firstIteration=0;
          curr=b1->e;
        }
        if(edgesAreEqual(currE, m->e[ edge ])==0){
          if(currE->v1==v1 || currE->v1==v2) currV=m->v[ currE->v2 ];
          else currV=m->v[ currE->v1 ];
          counter++;
          newc[0]=newc[0]-currV->x; newc[1]=newc[1]-currV->y; newc[2]=newc[2]-currV->z;
          k=k+dotProduct(newc, newc, 3);
        }

      }
      for(int i=0; i<3; i++) H[i][i]=2*counter;
      for(int i=0; i<3; i++) H[i][3]=2*newc[i];
      int prev=n_constraints_accepted;
      n_constraints_accepted=quadraticOpt(m, constr, n_constraints_accepted, H);
      cases[n_constraints_accepted]++;
      deallocMatrix(H, 4, 4);
      break;
  }
  free(e1); free(e2); free(e3);
  free(newC); free(t); free(v1_e1); free(v2_e1);
    return n_constraints_accepted;
}

void calculateSolutionsLength(Mesh *m, int n){
  double *v1=(double*)malloc(sizeof(double)*3);
  double *v2=(double*)malloc(sizeof(double)*3);
  Vertex *a=m->v[ m->e[n]->v1 ], *b=m->v[ m->e[n]->v2 ];
  v1[0]=a->x; v1[1]=a->y; v1[2]=a->z;
  v2[0]=b->x; v2[1]=b->y; v2[2]=b->z;
  double distance=fabs(distance3d(v1, v2));
  m->e[n]->cost=distance;
  m->solutions[n]=newVector(  (v1[0]+v2[0])/2 , (v1[1]+v2[1])/2 , (v1[2]+v2[2])/2 );
  free(v1); free(v2);
}


void calculateSolutionsGarland(Mesh *m, int n){
  double **quadratics=allocMatrix(4, 4);
  double *terms=allocVector(4), *firstPoint=allocVector(3);
  terms[0]=0.0; terms[1]=0.0; terms[2]=0.0; terms[3]=1.0;
  List *l1=m->v[ m->e[n]->v1  ]->triangles, *l2=m->v[ m->e[n]->v2  ]->triangles;
  List *curr=l1;
  double **vertices=allocMatrix(3,3);
  double a, b, c, d;
  Triangle *t_i;
  while(curr!=NULL){
      t_i=(Triangle*)curr->value;
      curr=curr->next;
      vertices[0][0]=m->v[ t_i->v1 ]->x; vertices[0][1]=m->v[ t_i->v1 ]->y; vertices[0][2]=m->v[ t_i->v1 ]->z;
      vertices[1][0]=m->v[ t_i->v2 ]->x; vertices[1][1]=m->v[ t_i->v2 ]->y; vertices[1][2]=m->v[ t_i->v2 ]->z;
      vertices[2][0]=m->v[ t_i->v3 ]->x; vertices[2][1]=m->v[ t_i->v3 ]->y; vertices[2][2]=m->v[ t_i->v3 ]->z;
      double *normal=normalUnit(vertices);
      for(int i=0; i<3; i++) firstPoint[i]=vertices[0][i];
      a=normal[0]; b=normal[1], c=normal[2];
      d=-1*dotProduct(normal, firstPoint, 3);
      quadratics[0][0]+=a*a; quadratics[0][1]+=a*b; quadratics[0][2]+=a*c; quadratics[0][3]+=a*d;
      quadratics[1][0]+=b*a; quadratics[1][1]+=b*b; quadratics[1][2]+=b*c; quadratics[1][3]+=b*d;
      quadratics[2][0]+=c*a; quadratics[2][1]+=c*b; quadratics[2][2]+=c*c; quadratics[2][3]+=c*d;
      quadratics[3][0]+=d*a; quadratics[3][1]+=d*b; quadratics[3][2]+=d*c; quadratics[3][3]+=d*d;
      free(normal);

  }
  curr=l2;
  while(curr!=NULL){
    t_i=(Triangle*)curr->value;
    if(contains(l1, curr->value)==0){
      vertices[0][0]=m->v[ t_i->v1 ]->x; vertices[0][1]=m->v[ t_i->v1 ]->y; vertices[0][2]=m->v[ t_i->v1 ]->z;
      vertices[1][0]=m->v[ t_i->v2 ]->x; vertices[1][1]=m->v[ t_i->v2 ]->y; vertices[1][2]=m->v[ t_i->v2 ]->z;
      vertices[2][0]=m->v[ t_i->v3 ]->x; vertices[2][1]=m->v[ t_i->v3 ]->y; vertices[2][2]=m->v[ t_i->v3 ]->z;
      double *normal=normalUnit(vertices);
      for(int i=0; i<3; i++) firstPoint[i]=vertices[0][i];
      a=normal[0]; b=normal[1], c=normal[2];
      d=-1*dotProduct(normal, firstPoint, 3);
      quadratics[0][0]+=a*a; quadratics[0][1]+=a*b; quadratics[0][2]+=a*c; quadratics[0][3]+=a*d;
      quadratics[1][0]+=b*a; quadratics[1][1]+=b*b; quadratics[1][2]+=b*c; quadratics[1][3]+=b*d;
      quadratics[2][0]+=c*a; quadratics[2][1]+=c*b; quadratics[2][2]+=c*c; quadratics[2][3]+=c*d;
      quadratics[3][0]+=d*a; quadratics[3][1]+=d*b; quadratics[3][2]+=d*c; quadratics[3][3]+=d*d;
      free(normal);
    }
    curr=curr->next;
  }

  //calculate quadratics matrix summing for Garland, it shall sum for every triangle

  double **garland=allocMatrix(4,4);
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++) garland[i][j]=quadratics[i][j];
  }
  garland[3][0]=0.0; garland[3][1]=0.0; garland[3][2]=0.0; garland[3][3]=1.0;
  double *solution=solveLinearSystem(garland, terms, 4);
  if(solution!=NULL){
    m->solutions[n]=newVector( solution[0] , solution[1] , solution[2] );
    double *tmpVector=matvet(quadratics, 4, 4, solution, 4);
    m->e[n]->cost=dotProduct( solution , tmpVector , 4);
    free(solution); free(tmpVector);
  }else{
    m->solutions[n]=newVector(9999, 9999, 9999); numberOfNullSolutions++;
    m->e[n]->cost=9999999.99  ;
  }
  deallocMatrix(garland, 4, 4); deallocMatrix(quadratics, 4, 4); free(firstPoint); free(terms);
  deallocMatrix(vertices, 3, 3);
}

void calculateSolutionsLindstrom(Mesh *m, int n){
  Constraint *c=(Constraint*)malloc(sizeof(Constraint));
  Fv_data *fv=newFv();
  Fv_data *fv2=newFv();
  fv->H=allocMatrix(4,4);
  fv2->H=allocMatrix(4,4);
  double *solution, *vv1, vv2, cv;
  double **currA=allocMatrix(3, 3), *currB=allocVector(3);
  int counter=0;
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++) c->A[i][j]=0.0;
  }
    int n_constr=getAllConstraints(m, c, fv, fv2, 0, n);
    if(n_constr!=3){
      m->solutions[n]=newVector(9999, 9999, 9999); numberOfNullSolutions++;
      m->e[n]->cost=9999999.99  ;
      return;
    }
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++) currA[i][j]=c->A[i][j];

    }
    for(int i=0; i<3; i++)  currB[i]=-1*c->b[i];
    solution=solveLinearSystem(currA , currB , 3);
    if(solution==NULL){
      m->solutions[n]=newVector(9999, 9999, 9999);  numberOfNullSolutions++;
      m->e[n]->cost=9999999.99  ;
    }else{
      m->solutions[n]=newVector(solution[0], solution[1], solution[2]);
      double *currSol=allocVector(4);
      for(int i=0; i<3; i++) currSol[i]=solution[i];
      currSol[3]=1.0;
      vv1=matvet(fv->H, 4, 4, currSol, 4);
      vv2=dotProduct(vv1, currSol, 4);
      m->e[n]->cost=0.5*(vv2);
      free(vv1);
      vv1=matvet(fv2->H, 4, 4, currSol, 4);
      vv2=dotProduct(vv1, currSol, 4);
      double *pointV1=allocVector(3), *pointV2=allocVector(3);
      pointV1[0]=m->v[ m->e[n]->v1  ]->x; pointV1[1]=m->v[ m->e[n]->v1  ]->y; pointV1[2]=m->v[ m->e[n]->v1  ]->z;
      pointV2[0]=m->v[ m->e[n]->v2  ]->x; pointV2[1]=m->v[ m->e[n]->v2  ]->y; pointV2[2]=m->v[ m->e[n]->v2  ]->z;
      double len=distance3d(pointV1, pointV2);
      m->e[n]->cost+=0.5*(vv2)*len*len;
      free(vv1); free(currSol); free(pointV1); free(pointV2);
      free(solution);

    }
    counter++;
    if(counter%20000==0) printf("counter %d \n", counter);
  deallocMatrix(currA, 3, 3); free(currB);
  deallocMatrix(fv->H, 4, 4);
  deallocMatrix(fv2->H, 4, 4);
  free(fv); free(c); free(fv2);
}

Mesh *newMesh(int nV, int nT, int choice){
  Mesh *m=(Mesh*)malloc(sizeof(Mesh));
  m->numV=nV; m->numT=nT;
  m->v=(Vertex**)malloc(sizeof(Vertex*)*nV);
  m->numE=0; m->dimE=2048;
  m->e=(Edge**)malloc(sizeof(Edge*)*2048);
  m->solutions=(Vector**)malloc(sizeof(Vector*)*2048);
  m->t=(Triangle**)malloc(sizeof(Triangle*)*nT);
  m->normals=(Vector**)malloc(sizeof(Vector*)*nT);
  switch(choice){
    case 0:
      m->calculateSolutions=&calculateSolutionsLength;
      break;
    case 1:
      m->calculateSolutions=&calculateSolutionsLindstrom;
      break;
    case 2:
      m->calculateSolutions=&calculateSolutionsGarland;
      break;

  }
  return m;
}

Mesh *initialize(int choice){
  for(int i=0; i<4; i++) cases[i]=0;
  if(PATH_FILE==NULL){
    printf("error: path file is not defined!\n");
    return NULL;
  }
  FILE *fp=fopen(PATH_FILE, "r");
  if(fp==NULL) printf("ERROR IN INPUT FILE, MISSING\n");
  int nV, nT;
  fscanf(fp, "%d %d", &nV, &nT); //expected first row of file as number of vertex and triangles
  Mesh *ret=newMesh(nV, nT, choice);
  readVertices(ret, fp, nV);
  readTriangles(ret, fp, nT);
  STEP_EDGES_SIMPLIFIED=ret->numT/12 ;
  return ret;
}


int containsVertex(Triangle *t, int v){
  if(t->v1==v || t->v2==v || t->v3==v) return 1;
  return 0;
}

Edge *getLinkingEdge(Vertex *v1, Vertex *v2){
  List *curr=v1->e;
  Edge *edge;
  while(curr!=NULL){
    edge=(Edge*)curr->value;
    curr=curr->next;
    if(contains(v2->e, edge)>0){
      if(edge->isDeleted==1) printf("FATAL ERROR, linking edge is deleted\n");
      return edge;
    }
  }
  printf("error, no linking edge from v1 and v2\n");
  return NULL;
}

Edge *getLinkingEdgeByInt(Mesh *m, int v1, int v2){
  List *curr=m->v[v1]->e;
  Edge *edge;
  while(curr!=NULL){
    edge=(Edge*)curr->value;
    curr=curr->next;
    if(edge->v1==v2 || edge->v2==v2){
      if(edge->isDeleted==1) printf("FATAL ERROR, linking edge is deleted\n");
      return edge;
    }
  }
  printf("error, no linking edge from v1 and v2\n");
  return NULL;
}

int getThirdVertex( Triangle *t , int v1 , int v2){
  if((t->v1==v1) || (t->v1==v2)){
    if((t->v2==v1) || (t->v2==v2)){
      return t->v3;
    }else{
      return t->v2;
    }
  }else{
    return t->v1;
  }
}

int triangleContainsTwoVertices(Triangle *t, int v1, int v2){
  int n=0;
  if((t->v1==v1) || (t->v1==v2)) n++;
  if((t->v2==v1) || (t->v2==v2)) n++;
  if((t->v3==v1) || (t->v3==v2)) n++;
  if(n==3){
    printf("FATAL ERROR, triangle with 3 vertices has two equals %d %d %d, called upon %d %d\n", t->v1, t->v2, t->v3, v1, v2);
  }
  return (n==2);
}

void printTriangleList(List *l){
  if(l!=NULL){
    Triangle *currT;
    List *curr=l;
    printf("triangle: ");
    while(curr!=NULL){
      currT=(Triangle*)curr->value;
      curr=curr->next;
      printf(" {%d} %d %d %d -> ", currT->isDeleted, currT->v1, currT->v2, currT->v3);
    }
  }

  printf(" end\n");
}

void printEdgeList(List *l){
  Edge *currE;
  List *curr=l;
  printf("edge: ");
  while(curr!=NULL){
    currE=(Edge*)curr->value;
    curr=curr->next;
    printf("edge {%d} %d %d -> ", currE->isDeleted, currE->v1, currE->v2);
  }
  printf(" end\n");
}

int trianglesAreEqual(Triangle *t1, Triangle *t2){
  int n=0;
  if(t1->v1==t2->v1 || t1->v1==t2->v2 || t1->v1==t2->v3) n++;
  if(t1->v2==t2->v1 || t1->v2==t2->v2 || t1->v2==t2->v3) n++;
  if(t1->v3==t2->v1 || t1->v3==t2->v2 || t1->v3==t2->v3) n++;
  if(n==3) return 1;
  else return 0;
}

int checkTriangleIsPresent(List *l1, int v1, int v2, int v3){
  Triangle *currT;
  Triangle *t=newTriangle(10, v1, v2, v3);
  while(l1!=NULL){
    currT=(Triangle*)l1->value;
    l1=l1->next;
    if(trianglesAreEqual(currT, t)==1){
      //found triangles are equal
      free(t);
      return 1;
    }
  }
  free(t);
  return 0;
}

int checkTriangleListIntegrity(List *l1){ //0 if integrity is violated, 1 otherwhise
  Triangle *currT, *currT2;
  List *curr1=l1, *curr2=l1;
  int violation=0;
  while(curr1!=NULL){
    currT=(Triangle*)curr1->value;
    while(curr2!=NULL){
      currT2=(Triangle*)curr2->value;
      if(trianglesAreEqual(currT, currT2)==1 && currT!=currT2){//means a triangle is duplicated
        violation=1;
      }
      curr2=curr2->next;
    }
    curr2=l1;
    curr1=curr1->next;
  }
  if (violation==0) return 1;
  else return 0;
}

int checkEdgeListIntegrity(List *l1){ //0 if integrity is violated, 1 otherwhise
  Edge *currE, *currE2;
  List *curr1=l1, *curr2=l1;
  int violation=0;
  while(curr1!=NULL){
    currE=(Edge*)curr1->value;
    while(curr2!=NULL){
      currE2=(Edge*)curr2->value;
      if(edgesAreEqual(currE, currE2)==1 && currE!=currE2){//means a triangle is duplicated
        printf("violation on %d %d and %d %d\n", currE->v1, currE->v2, currE2->v1, currE2->v2);
        violation=1;
      }
      curr2=curr2->next;
    }
    curr2=l1;
    curr1=curr1->next;
  }
  if (violation==0) return 1;
  else return 0;
}

void deleteTriangle(Mesh *m, Triangle *t){
  t->isDeleted=1;
  m->v[ t->v1 ]->triangles=removeList(m->v[ t->v1 ]->triangles, t);
  m->v[ t->v2 ]->triangles=removeList(m->v[ t->v2 ]->triangles, t);
  m->v[ t->v3 ]->triangles=removeList(m->v[ t->v3 ]->triangles, t);
}

int edgeFormedByTwoIntegers(Edge *e, int v1, int v2){
  int c=0;
  if(e->v1==v1 || e->v2==v1) c++;
  if(e->v1==v2 || e->v2==v2) c++;
  if(c==2) return 1;
  else return 0;
}


int numberOccurrencesEdge(List *l, int v1, int v2){
  int ret=0;
  Edge *currE;
  while(l!=NULL){
    currE=(Edge*)l->value;
    if(edgeFormedByTwoIntegers(currE, v1, v2)==1) ret++;
    l=l->next;
  }
  return ret;
}

void deleteEdge(Mesh *m, Edge *e){
  Vertex *v1=m->v[ e->v1 ];
  Vertex *v2=m->v[ e->v2 ];
  v1->e=removeList(v1->e, e);
  v2->e=removeList(v2->e, e);
  e->isDeleted=1;
  e->cost=9999999.99;
  heapify(m, e->n+1);
}

void deleteEdgeDuplicates(Mesh *m,Triangle *t, int a, int b, int c){
  Vertex *v1=m->v[ a ], *v2=m->v[ b ], *v3=m->v[ c ];
  List *curr;
  Edge *currE;
  int found=0;
  curr=v2->e;
  List *curr2;
  Edge *newCurr;
  int check=0;
  while(curr!=NULL && found==0){
    currE=(Edge*)curr->value;
    curr=curr->next;
    if(edgeFormedByTwoIntegers(currE, a, b)==1 && numberOccurrencesEdge(v1->e, a, b)>1){
      v1->e=removeList(v1->e, currE);
      v2->e=removeList(v2->e, currE);
      currE->isDeleted=1;
      currE->cost=9999999.99  ;
      heapify(m, currE->n+1);
      found=1;
      curr2=v2->e;
      while(curr2!=NULL){
        newCurr=(Edge*)curr2->value;
        curr2=curr2->next;
        if((edgeFormedByTwoIntegers(newCurr, a, b)==1) && newCurr!=currE) check++;
      }
      if(check==0) printf("ERROR, edged deleted has no duplicate!\n");
    }
  }



  found=0;
  curr=v3->e;
  check=0;
  while(curr!=NULL && found==0){
    currE=(Edge*)curr->value;
    curr=curr->next;
    if((edgeFormedByTwoIntegers(currE, a, c)==1)&& numberOccurrencesEdge(v1->e, a, c)>1){
      v3->e=removeList(v3->e, currE);
      v1->e=removeList(v1->e, currE);
      currE->isDeleted=1;
      currE->cost=9999999.99;
      curr2=v3->e;
      while(curr2!=NULL){
        newCurr=(Edge*)curr2->value;
        curr2=curr2->next;
        if((edgeFormedByTwoIntegers(newCurr, a, c)==1)&& newCurr!=currE) check++;
      }
      if(check==0) printf("ERROR, edged deleted has no duplicate!\n");
      found=1;
    }
  }


}


/*Takes t1 and t2 which have v3 in common, e2 and e3 should collapse on one edge. v1 belongs to t1, v2 to t2. e2=(v1, v3) e3=(v2, v3) */
void attachTriangles(Mesh *m, int v1, int v2, int v3, Triangle *t1, Triangle *t2, Edge *e2, Edge *e3){
  /*Change coordinates of t2 so that it references e2 instead of e3 */
  printf("t2 is %d %d %d,v1 v2 v3 %d %d %d \n", t2->v1, t2->v2, t2->v3, v1, v2, v3);
  if(t2->v1 == v2) t2->v1=v1;
  else if(t2->v2 == v2) t2->v2=v1;
  else if(t2->v3 == v2) t2->v3=v1;
  else printf("FATAL ERROR in attaching\n");
  /*v3 indexes edges, among them there is e3, it shall be changed to e2 (e2 is already inserted, so it has only to remove e3) */
  printf("ret attach triangles;\n");
}

int checkMeshIntegrity(Mesh *m, int n){
  int violation=0;
  printf("BEGIN CHECK MESH INTEGRITY on vertex %d\n", n);
  for(int i=0; i<m->numV-1; i++){
    Vertex *v=m->v[i];
    Triangle *currT;
    Edge *currE;
    List *curr;
    curr=v->triangles;
    while(curr!=NULL){
      currT=(Triangle*)curr->value;
      if(currT->v1==n || currT->v2==n || currT->v3==n || currT->isDeleted==1){
        violation++;
        printf("\n\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ERROR detected in vertex %d, with triangle %d %d %d\n", i, currT->v1, currT->v2, currT->v3);
      }
      curr=curr->next;
    }

    curr=v->e;
    while(curr!=NULL){
      currE=(Edge*)curr->value;
      if(currE->v1==n || currE->v2==n || currE->isDeleted==1){
          violation++;
         printf("\n\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ERROR detected in vertex %d, with edge %d %d \n", i, currE->v1, currE->v2);
      }
      curr=curr->next;
    }
  }
  return violation;
}

void checkEdgeErrors(Mesh *m){
  List *curr;
  Edge *currE;
  for(int i=0; i<m->numV-1; i++){
    curr=m->v[i]->e;
    while(curr!=NULL){
      currE=(Edge*)curr->value;
      if(numberOccurrencesEdge(m->v[i]->e, currE->v1, currE->v2)>1){
        printf("vertex %d; ERROR EDGE DUPLICATED %d %d, number of occurrences %d \n", i, currE->v1, currE->v2, numberOccurrencesEdge(m->v[i]->e, currE->v1, currE->v2));
        printEdgeList(m->v[i]->e);
      }
      curr=curr->next;
    }
  }

}

void updateCost(Mesh *m, int i){
  Vertex *v=m->v[i];
  List *curr=v->e;
  Edge *currE;
  while(curr!=NULL){
    currE=(Edge*)curr->value;
    curr=curr->next;
    m->calculateSolutions(m, currE->n);
    heapify(m, currE->n+1);
  }
}

void printTriangle(Triangle *t){
  printf("%d %d %d \n", t->v1, t->v2, t->v3);
}

void writeOutput(Mesh *m){
  int nT=m->numT - trianglesDeleted;
  int counter=0;
  char outFile[304];
  for(int i=0; i<300; i++) outFile[i]=PATH_FILE[i];
  strcat(outFile, ".out");
  FILE *fp=fopen(outFile, "w");
  printf("tried to write on %s\n", outFile);
  if(fp==NULL) printf("ERROR IN INPUT FILE, MISSING\n");
  /*Every vertex shall have a mapping from old int pointer to new logic pointer represented by counter. If before was 10.000 and new is 200, then new vertex is simply index[10000] */
  int *index=malloc(sizeof(int)*m->numV);
  for( int i=0; i<m->numV ; i++ ){
    if(m->v[i]->n!=-1){
      index[i]=counter;
      counter++;
    }else index[i]=-1;
  }
  fprintf(fp, "%d %d\n", counter, nT);
  for(int i=0; i<m->numV; i++){
    if(m->v[i]->n!=-1){
      fprintf(fp, "%lf %lf %lf\n", m->v[i]->x, m->v[i]->y, m->v[i]->z);
    }
  }
  int polygon_n_size=3;
  for(int i=0; i<m->numT; i++){
    if(m->t[i]->isDeleted==0){
      if(index[ m->t[i]->v1 ]==-1 || index[ m->t[i]->v2 ]==-1 || index[ m->t[i]->v3 ]==-1) printf("FATAL ERROR ON INDEX %d\n", i);
      fprintf(fp, "%d %d %d %d\n", polygon_n_size, index[ m->t[i]->v1 ], index[ m->t[i]->v2 ], index[ m->t[i]->v3 ]);
    }
  }


//  printf("was numT %d, deleted %d, vertices %d \n", m->numT, trianglesDeleted, counter);

  free(index);
}


/*Deletes first n edges from the mesh according to the cost */
int simplification(Mesh *m, int n){
  printf("begin simplification\n");
  int start=0, deleted, i; //later this can be saved maybe to decide which is the index of first non deleted edge
  List *curr, *curr1, *last;
  printf("BEFORE: triangles deleted %d\n", trianglesDeleted);
  for(i=0; i<m->numV -1; i++){
    if(checkTriangleListIntegrity(m->v[i]->triangles)==0) printf("ERROR: integrity list violated\n");
  }
  int numOfComputations=0, casesIsZero=0;
  while(n>0){

    if(trianglesDeleted>=(m->numT/100*(100-STOPPING_PERCENTAGE))){
      printf("TRUE\n\n");
      break;
    }

    for(i=0; i<m->numE ; i++){
      if((m->e[i]->cost)<=9999990.99 && m->e[i]->isDeleted==0){

        numOfComputations++;
        if(i==0) casesIsZero++;
        //first find t1 and t2 upper and lower triangles of edge e to act reduction (it might be even one if edge is boundary)
        Vertex *v1=m->v[ m->e[i]->v1 ], *v2=m->v[ m->e[i]->v2 ];
        if(m->e[i]->v1==m->e[i]->v2){
          printf("equal vertices in arc, returning from iteration %d\n", i);
          return 0;
        }

        List *commonTriangles=NULL;
        commonTriangles=intersectList(v1->triangles, v2->triangles);
        Triangle *currentT;
        //if t1 is a common triangle of e1, t2 and t3 are the external triangles that shall now get attached from collapse of t1
        curr=commonTriangles;
        Vertex *v3;
        Edge *e2, *e3;
        int numV3;
        while(curr!=NULL){
          Triangle *thisT=(Triangle*)curr->value;
          numV3=getThirdVertex((Triangle*)curr->value , m->e[i]->v1 , m->e[i]->v2);
          v3=m->v[ numV3 ] ;
          if(thisT->isDeleted==1) printf("ERROR, trying to delete a triangle already deleted\n");
          else thisT->isDeleted=1;
          e2=getLinkingEdge(v1, v3);
          e3=getLinkingEdge(v2, v3);
          v1->triangles=removeList(v1->triangles, thisT);
          v2->triangles=removeList(v2->triangles, thisT);
          v3->triangles=removeList(v3->triangles, thisT);
          trianglesDeleted=trianglesDeleted+1;
            e3->isDeleted=1;
            e3->v1=-1; e3->v2=-1; e3->cost=9999999.99  ;
            heapify(m, e3->n+1);
            v2->e=removeList(v2->e, e3);
            v3->e=removeList(v3->e, e3);
          curr=curr->next;
        }

        //Delete e1 from v1 and v2
        v1->e=removeList(v1->e, m->e[i]);
        v2->e=removeList(v2->e, m->e[i]);

        //Edges from v2 shall be passed to v1. v2 has no more e2 and e3, so it should be a simple merge moving the edges
        curr=v2->e;
        Edge* currE;
        while(curr!=NULL){
          currE=(Edge*)curr->value;
          curr=curr->next;
          if(currE->v1==m->e[i]->v2) currE->v1=m->e[i]->v1;
          else if(currE->v2==m->e[i]->v2) currE->v2=m->e[i]->v1;
          else printf("FATAL ERROR: moving edge from v2 that does not contain v2. Edge is %d %d\n", currE->v1, currE->v2);
          v1->e=addList(v1->e, currE);
          v2->e=removeList(v2->e, currE);
        }

        //every triangle that was in v2 shall have changed coordinates to v1
        Triangle *currT;
        curr=v2->triangles;
        while(curr!=NULL){
          currT=(Triangle*)curr->value;
          curr=curr->next;
          if( currT->v1==m->e[i]->v2 ) currT->v1=m->e[i]->v1;
          else if( currT->v2==m->e[i]->v2 ) currT->v2=m->e[i]->v1;
          else if( currT->v3==m->e[i]->v2 ) currT->v3=m->e[i]->v1;
          else printf("FATAL ERROR\n");
          v1->triangles=addList(v1->triangles, currT);

          if (checkTriangleListIntegrity(v1->triangles)==0){
            trianglesDeleted=trianglesDeleted+1;
            if(currT->v1==m->e[i]->v1) deleteEdgeDuplicates(m, currT, m->e[i]->v1, currT->v2, currT->v3);
            else if(currT->v2==m->e[i]->v1) deleteEdgeDuplicates(m, currT, m->e[i]->v1, currT->v1, currT->v3);
            else if(currT->v3==m->e[i]->v1) deleteEdgeDuplicates(m, currT, m->e[i]->v1, currT->v1, currT->v2);
            else printf("FATAL ERROR in finding common vertex\n");
            deleteTriangle(m, currT);
          }
          v2->triangles=removeList(v2->triangles, currT);
        }

        curr=v1->e;
        while(curr!=NULL){
          currE=(Edge*)curr->value;
          curr=curr->next;
          if(numberOccurrencesEdge(v1->e, currE->v1, currE->v2)>1){
            deleteEdge(m, currE);
          }
        }


        double *oldSol=allocVector(3), *newSol=allocVector(3);
        oldSol[0]=v1->x; oldSol[1]=v1->y; oldSol[2]=v1->z;
        newSol[0]=m->solutions[i]->x; newSol[1]=m->solutions[i]->y; newSol[2]=m->solutions[i]->z;
        v2->n=-1;
        if(distance3d(oldSol, newSol) < MAX_DISTANCE_ALLOWED){
          v1->x=m->solutions[i]->x;
          v1->y=m->solutions[i]->y;
          v1->z=m->solutions[i]->z;
        }else{
          trianglesDeletedWithLength+=2; /*It's an extimation, might be more. It is only used for analytics */
          calculateSolutionsLength(m, i);
        }
        free(oldSol); free(newSol);
        m->e[i]->isDeleted=1;
        updateCost(m, m->e[i]->v1);
         m->e[i]->v1=-1; m->e[i]->v2=-1; m->e[i]->cost=9999999.99  ;
        heapify(m, i+1);
        i=m->numE;
      } else break;
    }
    n--;
  }
  printf("AFTER: triangles deleted %d\n", trianglesDeleted);
  printf("initial numT is %d, trianglesDeleted %d trianglesDeletedWithLength %d \n numE is %d, num of computations is %d, cases is zero %d, STEP_EDGES_SIMPLIFIED = %d \n", m->numT, trianglesDeleted, trianglesDeletedWithLength, m->numE, numOfComputations, casesIsZero, STEP_EDGES_SIMPLIFIED);
  return 0;
}
