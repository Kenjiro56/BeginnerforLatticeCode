#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

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
void GramSchmidt(mat_RR& base, mat_RR& GSO_matrix, mat_RR& GSO_mu_matrix){
  
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

}

vector<string> split(string& input, char delimiter){
   istringstream stream(input);
   string field;
   vector<string> result;
   while (getline(stream, field, delimiter)) {
      result.push_back(field);
   }
   return result;
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
   ifstream ifs("d=" + to_string(dimention) +"baseMatrix.csv");
   string line;
   while (getline(ifs, line)) {
   vector<string> strvec = split(line, ',');
      for (int i=0; i<strvec.size();i++){
         printf("%5d\n", stoi(strvec.at(i)));
      }    
   }
  
  /*
  確かめるための小さな次元の既定行列*/
   base[0][0] = 1;
   base[0][1] = 0;
   base[0][2] = 1;
   
   base[1][0] = 1;
   base[1][1] = -1;
   base[1][2] = 1;

   base[2][0] = 0;
   base[2][1] = 1;
   base[2][2] = 1;

   transpose(base,base);
   
   // cout << base[0] * base[0] << endl;

   // 基底行列の確認
   for(int i=0;i<dimention;i++){
      for(int j=0;j<dimention;j++){
         cout << base[i][j] << "\t";
      } 
      cout << "" << endl;
   }

   GramSchmidt(base,GSO_matrix,GSO_mu_matrix);

   cout << "グラムシュミット行列" << endl;
   for(int i=0;i<dimention;i++){
      for(int j=0;j<dimention;j++){
         cout << GSO_matrix[i][j]<< "\t";
      } 
      cout << "" << endl;  
   }
   cout << "グラムシュミット係数行列" << endl;
   for(int i=0;i<dimention;i++){
      for(int j=0;j<dimention;j++){
         cout << GSO_mu_matrix[i][j]<< "\t";
      } 
      cout << "" << endl;
   }

}
NTL_CLIENT