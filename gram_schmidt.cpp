#include <sstream>
#include <fstream>
#include <iostream>
#include <NTL/ZZ_p.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/mat_ZZ_p.h>
using namespace std;
using namespace NTL;

// 次元
#define dimention 3

/// @brief Gram-Schmidt直行化を行う関数(正規化はしてない)
/// @param base 基底行列
/// @param GSO_matrix グラムシュミット行列
/// @param GSO_mu_matrix グラムシュミット係数行列
void GramSchmidt(mat_RR base, mat_RR GSO_matrix, mat_RR GSO_mu_matrix){
  
   transpose(base,base);
  
   for(int i = 0;i < dimention; i++){
      GSO_matrix[i] = base[i];
      
      for(int j = 0;j < i ;j++){
         RR product, square_bj, mu;
         product = base[i] * GSO_matrix[j];
         square_bj = GSO_matrix[j] * GSO_matrix[j];
         // cout << product << endl;
         div(mu,product,square_bj);
         GSO_mu_matrix[i][j] = mu;
         GSO_matrix[i] = GSO_matrix[i] - mu * GSO_matrix[j];
      }
   }
   transpose(GSO_matrix,GSO_matrix);
   cout << "グラムシュミット行列" << endl;
   for(int i=0;i<dimention;i++){
      for(int j=0;j<dimention;j++){
         cout << GSO_matrix[i][j]<< "\t";
      } 
      cout << "" << endl;  
   }
}



int main()
{
   // 基底行列
   mat_RR base;
   // グラムシュミット直行化行列
   mat_RR GSO_matrix;
   // グラムシュミット直行化係数行列
   mat_RR GSO_mu_matrix;

   // 次元の設定
   base.SetDims(dimention,dimention);
   GSO_matrix.SetDims(dimention,dimention);
   GSO_mu_matrix.SetDims(dimention,dimention);

   // 入力ファイル
   // ifstream base_file;
   // string filename = "base.txt";
   // base_file.open(filename, ios::in);
   // string reading_line_buffer;
   // while(getline(base_file, reading_line_buffer)){
      
   // }
   /*

   この例だと

     1      1/√6    -1/√3
     1     -1/√6     1/√3
     0       2/√6     1/√3

   になるはず。
   */ 
   base[0][0] = 1;
   base[0][1] = 0;
   base[0][2] = 1;

   base[1][0] = 1;
   base[1][1] = -1;
   base[1][2] = 1;

   base[2][0] = 0;
   base[2][1] = 1;
   base[2][2] = 1;
   // cout << base[0] * base[0] << endl;

   // 基底行列の確認
   // for(int i=0;i<dimention;i++){
   //    for(int j=0;j<dimention;j++){
   //       cout << base[i][j] << "\t";
   //    } 
   //    cout << "" << endl;
   // }

   GramSchmidt(base,GSO_matrix,GSO_mu_matrix);

   // cout << "グラムシュミット行列" << endl;
   // for(int i=0;i<dimention;i++){
   //    for(int j=0;j<dimention;j++){
   //       cout << GSO_matrix[i][j]<< "\t";
   //    } 
   //    cout << "" << endl;  
   // }
   // cout << "グラムシュミット係数行列" << endl;
   // for(int i=0;i<dimention;i++){
   //    for(int j=0;j<dimention;j++){
   //       cout << GSO_mu_matrix[i][j]<< "\t";
   //    } 
   //    cout << "" << endl;
   // }

}
NTL_CLIENT