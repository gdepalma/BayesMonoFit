#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>

static double mspline(double x, int k, int i, double *T){

  double d1, d2, v;
  

  if(x<T[i]) return 0;
  if(x>=T[i+k]) return 0;
  
  if(k==1) return 1/(T[i+1]-T[i]);
    
  d1 = x-T[i];
  d2 = T[i+k]-x;
 
  v = d1*mspline(x, k-1, i, T) + d2*mspline(x, k-1, i+1, T);
  v = v*k;
  v = v/((k-1)*(T[i+k]-T[i]));
 
  return v;
}

static double ispline(double x, int k, int i,double *T,int lknots){

    double v;
    int m,j=0;

    do{
      j=j+1;
    }while(T[j]<=x && j<lknots);
    j=j-1;
    
    if (j<i) return 0;
    if (j-k+1>i) return 1;
 
    v = 0;
    for (m=i; m<=j; m++)
        v =v+ mspline(x, k+1, m, T)*(T[m+k+1]-T[m])/(k+1);
    

   return v;
}

void getIspline(double *grid, int *lgrid, double *knotseq, double *mat,int *numBases) {
  int k=3;
  int i,j,lknots,idx=0;      
  
  lknots=*numBases+5;
    
  for(i=1; i<=*numBases; i++)     
    for(j=0; j<*lgrid; j++){
        mat[idx]=1-ispline(grid[j],k,i,knotseq,lknots);
        idx=idx+1;
     }
  
}

void getytrue(double *icoefs, double *knotseq, double *xtrue, double *ytrue, int *numBases, int *lxtrue) {
  int k=3;
  int i,j,lknots;  
  
  lknots=*numBases+5;
    
  for(i=0; i<*lxtrue; i++){    
    for(j=1; j<=*numBases; j++){
        ytrue[i]=ytrue[i]+(1-ispline(xtrue[i],k,j,knotseq,lknots)*icoefs[j]);
     }
     ytrue[i]=ytrue[i]+icoefs[0];
   }  
}

void isMonotone(double *y, int *nobs, int *mono){

  int i=0;

  for(i=0; i<*nobs-1; i++)
    if (y[i]<y[i+1]){
      *mono=0;
      break;
    }
  

}




