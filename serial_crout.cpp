void lu_serial( double **A, double **L, double **U, int n){
    
      for (int i = 0; i < n; i++) {
          for (int j = 0; j < n; j++) {
            if (j < i)
              L[j][i] = 0;
            else {
                L[j][i] = A[j][i];
                for (int k = 0; k < i; k++) {
                  L[j][i] = L[j][i] - L[j][k] * U[k][i];
                }
            }
          }
          for (int j = 0; j < n; j++) {
            if (j < i)
            U[i][j] = 0;
            else if (j == i)
            U[i][j] = 1;
            else {
                U[i][j] = A[i][j] / L[i][i];
                for (int k = 0; k < i; k++) {
                  U[i][j] = U[i][j] - ((L[i][k] * U[k][j]) / L[i][i]);
                }
            }
          }
      }
    // for(int i =0;i<n;i++){
      
      
      
      // // calculate u
      
      // for( int j=i;j<n;j++){
      //   int new_item=A[i][j];
      //   for( int k=0;k<i;k++){
      //       new_item -= ( L[i][k] * U[k][j] );
      //   }
      //   U[i][j]= new_item ;
      // }
      
      // // calculate l
      // // here i behaves as second index
      // for( int i_n=i;i_n<n;i_n++){
      //   if (i_n==i)
      //     L[i][i]=1;
      //   else{
      //     int new_item=A[i_n][i];
      //     for( int j=0;j<i;j++){
      //          new_item -= (L[i_n][j]*U[j][i]);
      //     }
      //     L[i_n][i]=new_item/U[i][i];
      //   }
      // }
    // }
  } 
