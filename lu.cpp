#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <omp.h>
#include <fstream>
#include <iomanip>
#include "plu.h"
using namespace std;
// # define DEBUG 
// # define CHECK_CORRECTNESS
int main(){
    
    pLU p= pLU();
    double **A;
    double **L, **U;
    int n;
    ifstream infile;
    infile.open("input_dir/2900.txt");
    infile>>n;
    cout<<n<<"\n";
    A=new double*[n];
    L=new double*[n];
    U=new double*[n];
    // cout<<"getting input\n";
    
    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
        L[i] = new double[n];
        U[i] = new double[n];
        for(int j=0;j<n;j++){
            infile>>A[i][j];
            U[i][j]=0;
            L[i][j]=0;
        }
    }
    
    // const double **A2 = *(&A);
    p.lu(A,L,U,n);
    // cout<<"printing output\n";
    
    // cout<<"L\n";
    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < n; j++){
    //         cout << setw(3) << L[i][j] << " ";}
    //      cout << endl;
    // }
    
    // cout<<"U\n";
    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < n; j++){
    //         cout << setw(3) << U[i][j] << " ";
    //         }
    //      cout << endl;
    // }
#ifdef DEBUG    
    cout<<"A orig\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cout << setw(3) << A[i][j] << " ";
         cout << endl;
    }
    cout<<"A new\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = 0;
            for (int k = 0; k < n; k++) {
                A[i][j] += L[i][k] * U[k][j];
            }
            if(A[i][j]<0.001)
                A[i][j]=0.0;
            cout <<  setw(3) <<A[i][j] << " ";
        }
        cout << endl;
    }
#endif

#ifdef CHECK_CORRECTNESS
    cout<<"Checking correctness.\n";
    int errors_found =0;
    double **A_new =new double*[n];
    for (int i = 0; i < n; i++) {
        A_new[i] = new double[n];
        for (int j = 0; j < n; j++) {
            A_new[i][j] = 0;
            for (int k = 0; k < n; k++) {
                A_new[i][j] += L[i][k] * U[k][j];
            }
            if(fabs(A_new[i][j]-A[i][j])>0.01){
                cout<<"Found error"<<A_new[i][j]<<"is actually"<<A[i][j]<<endl;
                errors_found+=1;
                }
        }
    }
    cout<<"errors_found="<<errors_found<<endl;
#endif
}

