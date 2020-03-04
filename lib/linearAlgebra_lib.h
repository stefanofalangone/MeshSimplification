#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double dotProduct(double *a, double *b, int dim){
  double ret=0.0;
  for(int i=0; i<dim; i++){
    ret=ret+a[i]*b[i];
  }
  return ret;
}

double *orthogonalVectorToVector(double *a){ //a shall be 3dim. Vector calculated with scalar product imposed zero
  double *ret=(double*)malloc(sizeof(double)*3);
  if(fabs(a[2])>=0.00001){
    //printf("case 1\n");
    ret[0]=1; ret[1]=1; ret[2]=(-a[0]-a[1])/a[2];
  }else if(fabs(a[1])>=0.00001){
    //printf("case 2\n");
    ret[0]=1; ret[2]=1; ret[1]=(-a[0]-a[2])/a[1];
  }else if(fabs(a[0])>=0.00001){
    //printf("case 3\n");
    ret[1]=1; ret[2]=1; ret[0]=(-a[1]-a[2])/a[0];
  }else{
    //printf("case 4\n");
    ret[0]=1; ret[1]=1; ret[2]=1;
  }
  for(int i=0; i<3; i++) ret[i]=ret[i]*10000;
  return ret;
}

/*double *orthogonalVectorToVector(double *a){ //a shall be 3dim. Vector calculated with scalar product imposed zero
  double *ret=(double*)malloc(sizeof(double)*3);
  if(fabs(a[2])>=0.00001){
    //printf("case 1\n");
    ret[0]=1; ret[1]=1; ret[2]=0;
  }else if(fabs(a[1])>=0.00001){
    //printf("case 2\n");
    ret[0]=1; ret[2]=1; ret[1]=0;
  }else if(fabs(a[0])>=0.00001){
    //printf("case 3\n");
    ret[1]=1; ret[2]=1; ret[0]=0;
  }else{
    //printf("case 4\n");
    ret[0]=0; ret[1]=1; ret[2]=1;
  }
  return ret;
}*/

double *allocVector(int dim){
  double *result = (double*)malloc( sizeof(double)*dim );
  for(int j=0; j<dim; j++) result[j]=0.0;
  return result;
}



void deallocMatrix(double **M, int rows, int cols){
  if(M==NULL) return;
  for(int i=0; i<rows; i++){
    free(M[i]);
  }
  free(M);
}


int vectorIsNull(double *a, int dim){
  for(int i=0; i<dim; i++){
    if(a[i]>=0.0001) return 0;
  }
  return 1;
}





double **allocMatrix(int row, int col){
  double **result = (double**)malloc( sizeof(double*)*row );
  for(int i=0; i<row; i++){
    result[i]=(double*)malloc(sizeof(double)*col);
    for(int j=0; j<col; j++) result[i][j]=0.0;
  }
  return result;
}

double *scalarVector(double s, double *v, int dim){
  double *ret=(double*)malloc(sizeof(double)*dim);
  for(int i=0; i<dim; i++){
    ret[i]=s*v[i];
  }
  return ret;
}

double *sumVectors(double *a, double *b, int dim){
  double *ret=(double*)malloc(sizeof(double)*dim);
  for(int i=0; i<dim; i++){
    ret[i]=a[i]+b[i];
  }
  return ret;
}

double *subtractVectors(double *a, double *b, int dim){
  double *ret=(double*)malloc(sizeof(double)*dim);
  for(int i=0; i<dim; i++){
    ret[i]=a[i]-b[i];
  }
  return ret;
}

double distance3d(double *old, double *new){
  double *difference=subtractVectors(old, new, 3 );
  double ret=sqrt( dotProduct( difference, difference, 3  )   );
  free(difference);
  return ret;
}

double det2x2(double **A){
  return (A[0][0]*A[1][1] - A[0][1]*A[1][0]);
}

/*double det3x3(double **A){ //optimizable with sarrus's rule
  double **temp=allocMatrix(2,2);
  double ret=0.0;
  temp[0][0]=A[1][1]; temp[0][1]=A[1][2]; temp[1][0]=A[2][1]; temp[1][1]=A[2][2];
  ret=A[0][0]*det2x2(temp);
  temp[0][0]=A[1][0]; temp[0][1]=A[1][2]; temp[1][0]=A[2][0]; temp[1][1]=A[2][2];
  ret=ret - A[0][1]*det2x2(temp);
  temp[0][0]=A[1][0]; temp[0][1]=A[1][1]; temp[1][0]=A[2][0]; temp[1][1]=A[2][1];
  ret=ret + A[0][2]*det2x2(temp);
  free(temp[0]); free(temp[1]); free(temp);
  return ret;
}*/

double det3x3(double **a){ //optimizable with sarrus's rule
  double x=0.0, y=0.0;
  for(int i=0;i<3;i++)
  {
    x = x + a[0][i] * a[1][(i+1)%3] * a[2][(i+2)%3] ;
    y = y + a[2][i] * a[1][(i+1)%3] * a[0][(i+2)%3] ;
  }
  return (x-y);
}




double **matrixComplement(double **A, int n, int m, int a, int b){
 /*A the matrix is (n, m). We want the algebric complement to element (a,b) */
 double **ret=allocMatrix(n-1, m-1);
 int retI=0;
 for(int i=0; i<n; i++){
  int retJ=0;
  if (i==a)
   continue;
  /*Implicit else*/
  for(int j=0; j<m; j++){
   if (j==b)
    continue;
   /*Implicit else*/
   ret[retI][retJ]=A[i][j];
   retJ++;
  }
  retI++;
 }
 return ret;
}

double det4x4(double **A){
  double **tmp;
  double ret=0.0;
  for(int i=0; i<4; i++){
    tmp=matrixComplement(A, 4, 4, 0, i);
    if(i%2==0) ret=ret + A[0][i]*det3x3(tmp);
    else ret=ret - A[0][i]*det3x3(tmp);
    deallocMatrix(tmp, 3, 3);
  }
  return ret;
}



double **transpose(double **A, int m, int n){
    //takes A m x n
    double **B=(double **)malloc(sizeof(double*)*n);
    for(int i=0; i<n; i++) B[i]=(double *)malloc(sizeof(double)*m);
    for(int i=0; i<m; i++){
      for(int j=0; j<n; j++) B[j][i]=A[i][j];
    }
    return B;
}

double *crossProduct3dim(double *a, double *b){
  double *ret=(double*)malloc(sizeof(double)*3);
  ret[0]=a[1]*b[2]-a[2]*b[1];
  ret[1]=a[2]*b[0]-a[0]*b[2];
  ret[2]=a[0]*b[1]-a[1]*b[0];
  return ret;
}

double norm(double *v, int dim){
  double ret=0.0;
  for(int i=0; i<dim; i++){
    ret=ret + v[i]*v[i];
  }
  ret=sqrt(ret);
  return ret;
}

double cosTwoVectors(double *v1, double *v2, int dim){
  return (dotProduct(v1, v2, dim))/(norm(v1, dim)*norm(v2, dim));
}

double sinTwoVectors3dim(double *v1, double *v2){
  double *crossP=crossProduct3dim(v1, v2);
  double res=norm(crossP, 3)/(norm(v1, 3)*norm(v2, 3));
  free(crossP);
  return res;
}


double *solveLinearSystem(double **A, double *b, int dim){ //solves up to 4x4 linear systems
  double *ret=(double*)malloc(sizeof(double)*dim);
  for(int i=0; i<dim ; i++) ret[i]=0.0;
  double detA;
  if(dim==2) detA=det2x2(A);
  else if(dim==3) detA=det3x3(A);
  else if(dim==4) detA=det4x4(A);
  if(fabs(detA)<=0.0001){
    return NULL;
  }
  double currDet;
  double col[4];
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      col[j]=A[j][i];
      A[j][i]=b[j];
    }

    if(dim==3) currDet=det3x3(A);
    else if(dim==2) currDet=det2x2(A);
    else if(dim==4) currDet=det4x4(A);
    ret[i]=currDet/detA;
      for(int j=0; j<dim; j++) A[j][i]=col[j];
  }
  return ret;
}



double **matmat( double **m1, int row1, int col1, double **m2, int row2, int col2 ){
  int i, j, m;
  if(col1!=row2){
    printf("mat product not doable, error\n");
    return NULL;
  }
  double **result = allocMatrix(row1, col2);
  /* Perform multiplication algorithm */
  for( m = 0; m < row1 ; m++ ) {
    for( i = 0; i < col2 ; i++ ) {
      for( j = 0; j < row2 ; j++ ) {

        result[ m ][ i ] += m1[ m ][ j ] *
                                      m2[ j ][ i ];
      }
    }
  }
  return result;
}

double *matvet(double **A, int rows, int cols, double *c, int dim){
  if(cols!=dim){
    printf("mat product not doable, error\n");
    return NULL;
  }
  double *ret=(double*)malloc(sizeof(double)*rows);
  for(int i=0; i<rows; i++){
    ret[i]=0.0;
    for(int j=0; j<cols; j++) ret[i]=ret[i] + A[i][j]*c[j];
  }
  return ret;
}

void printMatrix(double **A, int rows, int cols){
  for(int i=0; i<rows; i++){
      for(int j=0; j<cols; j++) printf("%lf ", A[i][j]);
      printf("\n");
    }
}

void printVector(double *v, int dim){
  for(int i=0; i<dim; i++){
      printf("%lf ", v[i]);
    }
    printf("\n");
}

double **inverseMatrix3x3(double **A){ //assumes matrix is squared and 3x3
  double det=det3x3(A);
  if(fabs(det)<0.00001){
    //matrix is not invertible, det=0
    return NULL;
  }
  double **res=allocMatrix(3, 3);
  double **tmp=allocMatrix(2,2);
  int nextx[]={0, 0, 1, 1}, nexty[]={0, 1, 0, 1}, next;
  double val;
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      //calculate algebraic complement/det and put it in trasposed position
      next=0;
      for(int k=0; k<3; k++){
        for(int m=0; m<3; m++){
          if(k!=i && m!=j){
            tmp[ nextx[next] ][ nexty[next] ]=A[k][m];
            next++;
          }
        }
      }
      if((i+j)%2==0) val=(det2x2(tmp)/det);
      else val=(double)(-1.0*det2x2(tmp)/det);
      res[j][i]=val;
    }
  }
  deallocMatrix(tmp, 2, 2);
  return res;
}

double **inverseMatrix4x4(double **A){ //assumes matrix is squared and 3x3
  double det=det4x4(A);
  //printf("matrix given\n");
  //printMatrix(A, 3, 3);
  if(fabs(det)<0.0001){
    //printf("matrix is not invertible, det=0\n");
    return NULL;
  }
  double **res=allocMatrix(4, 4);
  double **tmp;
  double val;
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      //calculate algebraic complement/det and put it in trasposed position
      tmp=matrixComplement(A, 4, 4, i, j);
      //printf("algebric complement of %d %d\n", i, j);
      //printMatrix(tmp, 2, 2);
      if((i+j)%2==0) val=(det3x3(tmp)/det);
      else val=(double)(-1.0*det3x3(tmp)/det);
      res[j][i]=val;
      free(tmp);
    }
  }
  //deallocMatrix(tmp, 2, 2);
  return res;
}

double *normalToTriangle(double **v){
  double x=0, y=0, z=0, *res=malloc(sizeof(double)*3);
  for(int j=0; j<2; j++){ //Newell method to calculate normal for a triangle
    x=x + (v[j][1] - v[j+1][1])*(v[j][2] + v[j+1][2]);
    y=y + (v[j][2] - v[j+1][2])*(v[j][0] + v[j+1][0]);
    z=z + (v[j][0] - v[j+1][0])*(v[j][1] + v[j+1][1]);
  }
  x=x + (v[2][1] - v[0][1])*(v[2][2] + v[0][2]);
  y=y + (v[2][2] - v[0][2])*(v[2][0] + v[0][0]);
  z=z + (v[2][0] - v[0][0])*(v[2][1] + v[0][1]);
  //printf("normals calcolated %lf %lf %lf\n", x, y, z);
  res[0]=x; res[1]=y; res[2]=z;
  return res;
}

double *normalUnit(double **v){
  double *ret=normalToTriangle(v);
  double length=norm(ret, 3);
  for(int i=0; i<3; i++) ret[i]=ret[i]/length;
  return ret;
}


/*double *orthogonalVectorTo2Vectors(double *a, double *b){
  double *ret=(double*)malloc(sizeof(double)*3);
  if(fabs(a[2])>=0.00001 && fabs(b[2])>=0.00001){
    //printf("case 1\n");
    ret[0]=1; ret[1]=1; ret[2]=0;
  }else if(fabs(a[1])>=0.00001 && fabs(b[1])>=0.00001){
    //printf("case 2\n");
    ret[0]=1; ret[2]=1; ret[1]=0;
  }else if(fabs(a[0])>=0.00001 && fabs(b[0])>=0.00001){
    //printf("case 3\n");
    ret[1]=1; ret[2]=1; ret[0]=0;
  }else{
    //printf("case 4\n");
    ret=crossProduct3dim(a, b);
  }
  return ret;
}*/

double *orthogonalVectorTo2Vectors(double *a, double *b){
  double *ret=crossProduct3dim(a, b);
  for(int i=0; i<3; i++) ret[i]*=10000;
  return ret;
}
