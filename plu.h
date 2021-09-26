#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
using namespace std;
#define B_BLOCK 32


#pragma once

class pLU {

public:
  int num_procs = 1;
  void init(){
    this->num_procs = omp_get_num_procs();
  }
  void close();
  
  void trsm(double **A, double **L, double **U, int n, int block_start_i, int u_block_start_j){
    // we have L00 and want to compute U01 using A=L00*U01
    // cout<<"starting trsm"<<block_start_i<<","<< u_block_start_j<<"\n";
    int matrix_size = min(n-block_start_i,  B_BLOCK);
    int num_cols_u = min(n-u_block_start_j,  B_BLOCK); 

    for( int  u_idx=0;u_idx<num_cols_u;u_idx++){
      for( int i=0; i<matrix_size;i++){
        double this_item = A[block_start_i+i][u_block_start_j+u_idx] ;
        for( int j=0;j<i;j++){
          this_item -= L[block_start_i+ i][block_start_i+j] * U[block_start_i+j][u_block_start_j+u_idx];
        } 
        // this_item = this_item/L[block_start_i+ i][block_start_i+i]; 
        U[block_start_i+i][u_block_start_j+u_idx]=this_item;
      }
    }
    // cout<<"finished trsm"<<block_start_i<<","<< u_block_start_j<<"\n";
  }

  void trsm_t(double **A, double **L, double **U, int n, int l_block_start_i, int block_start_j){
    // we have U00 and want to compute L10 using A=L01*U00
    // cout<<"starting trsm_t"<<l_block_start_i<<","<< block_start_j<<"\n";
    int matrix_size = min(n-block_start_j,  B_BLOCK);
    int num_rows_l = min(n-l_block_start_i,  B_BLOCK); 

    for( int  l_idx=0;l_idx<num_rows_l;l_idx++){
      for( int j=0; j<matrix_size;j++){
        double this_item = A[l_block_start_i+l_idx][block_start_j+j] ;
        for( int i=0;i<j;i++){
          this_item -= U[block_start_j+i][block_start_j+j] * L[l_block_start_i+ l_idx][block_start_j+i]  ;
        } 
        this_item = this_item/U[block_start_j+ j][block_start_j+j]; 
        // cout<<"div"<<U[block_start_j+ j][block_start_j+j]<<"\n"; 
        L[l_block_start_i+l_idx][block_start_j+j] = this_item;
      }
    }
    // cout<<"finished trsm_t"<<l_block_start_i<<","<< block_start_j<<"\n";
  }


void seq_sub_block_multiply(double **A, double **L, double **U, int n,int block_start_orig=0,int block_start_L=0, int block_start_U=0){
  int p = min(n-block_start_L,  B_BLOCK);
  int q = min(n-block_start_U,  B_BLOCK);
  
  for (int i = 0; i < p; i++) {
        for (int j = 0; j < q; j++) {
            double this_ele = 0;
            for (int k = 0; k < B_BLOCK; k++) {
                this_ele += L[i+block_start_L][k+block_start_orig] * U[k+block_start_orig][j+block_start_U];
            }
            A[i+block_start_L][j+block_start_U]-=this_ele;
         }
    }
}

  
  void lu_p_crout( double **A, int *P,double **L, double **U, int n_orig,int block_start=0){  // parallelized crout Algo 2
    int n=min(n_orig-block_start,B_BLOCK);
    // Overwrite A on U
    for( int i=0;i<n;i++){
      for( int j=0;j<n;j++){
        U[i+block_start][j+block_start]=A[i+block_start][j+block_start];
      }
    }
    
    // Apply crout on U that contains elements of A
    for(int k=0; k<n;k++){
      
      // Pivoting
      int max_pivot_idx=k;
      for(int i=k+1;i<n;i++){
        if(U[block_start+i][block_start+k]>U[block_start+max_pivot_idx][block_start+k]){
          max_pivot_idx=i;
        }
      }
      
      // Now we want to replace row(max_pivot_idx) and row(block_start+k)
      int temp_idx = P[block_start+max_pivot_idx];
      P[block_start+max_pivot_idx] = P[block_start+k];
      P[block_start+k] = temp_idx;
      for( int j =block_start;j<block_start+n;j++){
        double temp_u = U[block_start+k][j];
        U[block_start+k][j] = U[block_start+max_pivot_idx][j];
        U[block_start+max_pivot_idx][j]=temp_u;
      }
      // now we apply same operations on U01
      for( int j =block_start+n;j<n_orig;j++){
        double temp_u = A[block_start+k][j];
        A[block_start+k][j] = A[block_start+max_pivot_idx][j];
        A[block_start+max_pivot_idx][j]=temp_u;
      } 

      // now we apply same operations on L
      for( int j =0;j<block_start;j++){
        double temp_u = L[block_start+k][j];
        L[block_start+k][j] = L[block_start+max_pivot_idx][j];
        L[block_start+max_pivot_idx][j]=temp_u;
      } 

      // divide by Ukk
      for( int i = k+1; i<n;i++){
        U[i+block_start][k+block_start]=U[i+block_start][k+block_start]/U[k+block_start][k+block_start];
      }
      // #pragma omp parallel for 
        for( int j=k+1;j<n;j++){
          for( int i =k+1;i<n;i++){
            U[i+block_start][j+block_start]=U[i+block_start][j+block_start]-U[i+block_start][k+block_start]*U[k+block_start][j+block_start];
          }
        }
    }

    // Move lower triangular matrix to L and set the Lower triangle of U 0 
    for( int i=0;i<n;i++){
      for( int j=0;j<n;j++){
        if (i<j){
          L[i+block_start][j+block_start]=0.0;
        }
        else if (i==j)
        {
          L[i+block_start][j+block_start]=1.0;
        }
        else{
          L[i+block_start][j+block_start]=U[i+block_start][j+block_start];
          U[i+block_start][j+block_start]=0.0;
        }
      }
    }
  }

void lu(  double **A_orig, double **L, double **U, int n){
    
    double **A=new double*[n];
    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
        for(int j=0;j<n;j++){
            A[i][j] = A_orig[i][j];
        }
    }

    int P[n] ;
    for( int i=0;i<n;i++){
      P[i]=i;
    }
    // if (n<=B_BLOCK){
    //   lu_p_crout(A, P, L, U, n,0);
    //   return;
    // }
    
    int block_start = 0;
    
    
  // #pragma omp parallel num_threads(8)
  // #pragma omp master
  {
      while(block_start<n){
        lu_p_crout( A,P, L, U,n,block_start);
        
        #pragma omp parallel num_threads(8)
        #pragma omp master
        {
          int block_start_j = block_start+B_BLOCK;
          while (block_start_j<n){
              //Task creation
              
              #pragma omp task depend(out:A[block_start][block_start_j])
                trsm(A, L, U, n, block_start, block_start_j);// computes U01

              #pragma omp task depend(out:A[block_start_j][block_start])
                trsm_t(A, L,U,n, block_start_j,block_start); // computes L10
              block_start_j+=B_BLOCK;
          }
          // #pragma omp taskwait 

          // parallel block multiplication A' =  A - L10 U01
          for(int sub_block_start_i = block_start+B_BLOCK;sub_block_start_i<n;sub_block_start_i+=B_BLOCK){
            for(int sub_block_start_j = block_start+B_BLOCK;sub_block_start_j<n;sub_block_start_j+=B_BLOCK){
              #pragma omp task depend(in:A[sub_block_start_i][block_start],A[block_start][sub_block_start_j])
                seq_sub_block_multiply(A,L,U,n,block_start,sub_block_start_i,sub_block_start_j);
            }
          }
          #pragma omp taskwait  
        }// end of master
        
        block_start+=B_BLOCK;
      }
    }
  // Delete The extra A 

  // reorder L
  for( int i=0;i<n;i++){
      // cout<<"P["<<i<<"]="<<P[i]<<'\n';
      while(P[i]!=i){
        //swap row(P[i]) row(i)
        for(int j=0;j<n;j++){
          double temp_u = L[i][j];
          L[i][j] = L[P[i]][j];
          L[P[i]][j]=temp_u;
        }
        int h = P[i];
        P[i]=P[P[i]];
        P[h]=h;
      }
    }
  // for( int i=0;i<n;i++){
  //     cout<<"P["<<i<<"]="<<P[i]<<'\n';
  // }
  }
};
