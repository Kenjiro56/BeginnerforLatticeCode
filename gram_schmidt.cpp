// 出入力ライブラリ
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// ディレクトリ作成ライブラリ
#include <sys/stat.h>

// NTLライブラリ
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
  transpose(base,base);
  // これが正解？
  ident(GSO_mu_matrix,dimention);
   for(int i = 0;i < dimention; i++){
      GSO_matrix[i] = base[i];
      
      for(int j = 0;j < i ;j++){
         RR product, square_bj, mu;
         product = base[i] * GSO_matrix[j];
         square_bj = GSO_matrix[j] * GSO_matrix[j];
         div(mu,product,square_bj);
         GSO_mu_matrix[i][j] = mu;
         GSO_matrix[i] = GSO_matrix[i] - mu * GSO_matrix[j];
      }
   }
   transpose(GSO_matrix,GSO_matrix);

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
   string base_file_path = "./BaseMatrix/d=" + to_string(dimention) + "baseMatrix.csv";
   ifstream file(base_file_path);
   string line;

   if(file.fail()){
      cout << "Failed to open file." << endl;
      return -1;
   }else{
      cout << "open successfully" << endl;
   }

   int col = 0;
   while(getline(file,line)){
      string value;
      istringstream stream(line);
      int row = 0;
      while(getline(stream,value,',')){
         base[row][col] = stoi(value);
         row++;
      }
      col++;
   }
  
   transpose(base,base);

   // // 基底行列の確認
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



   // ファイルの出力
   string folder_name = ("./GSOMatrix/d=" + to_string(dimention));
   if(mkdir(folder_name.c_str(),0777) == 0){
      cout << "ディレクトリ"<< folder_name << "を作成しました" << endl;
   }

   // 出力ファイル
   string GSOMatrix_csv_path = folder_name + "/GSOMatrix.csv";
   string GSOCoeff_csv_path = folder_name + "/GSOCoeff.csv";
   ofstream Matrix_csv(GSOMatrix_csv_path);
   ofstream Coeff_csv(GSOCoeff_csv_path);

   // csv書き込み
   for(int i=0;i<dimention;i++){
      for(int j=0;j<dimention;j++){
         if(j != 0){
            Matrix_csv << ",";
            Coeff_csv << ",";
         }
         Matrix_csv << GSO_matrix[i][j];
         Coeff_csv << GSO_mu_matrix[i][j];
      } 
      Matrix_csv << endl;
      Coeff_csv << endl;
   }
}
NTL_CLIENT