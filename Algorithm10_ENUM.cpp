/*
Algorithm10 ENUM : 格子上の最短ベクトルの数え上げを行うアルゴリズム
Input :
    - n次元格子Lの基底のGSO係数
    - GSOベクトルの2乗ノルム
    - 数え上げ上界列
Output :
    - (3.5)式を満たす格子ベクトルの係数ベクトル

方針 :
    ステップ1. k = nにおいて(3.5)式を満たすよう整数を一つ決める.
    ステップ2. k = n-1において(3.5)式を満たすよう整数を一つ決める.
    　　　　　 もし存在しない場合ステップ1に戻る.
    ステップ3. ステップ1,2を深さ優先探索で回し、n個の整数の組（最短ベクトルを成す格子ベクトルの係数ベクトル）を見つける

 　　※上界の決め方：
        Minkowskiの第一定理を用いて設定。短い格子ベクトルを見つけるたびに入力する数え上げ上界列を小さく設定すれば良い。
*/

// 出入力ライブラリ
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// ディレクトリ作成ライブラリ
#include <sys/stat.h>

// NTLライブラリ
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>
using namespace std;
using namespace NTL;

#define dimension 3 // 次元数

/// @brief Gram-Schmidt直行化を行う関数(正規化はしてない)
/// @param base 基底行列
/// @param GSO_matrix グラムシュミット行列
/// @param GSO_mu_matrix グラムシュミット係数行列
void GramSchmidt(mat_RR &base, mat_RR &GSO_matrix, mat_RR &GSO_mu_matrix)
{
    transpose(base, base);
    ident(GSO_mu_matrix, dimension);
    for (int i = 0; i < dimension; i++)
    {
        GSO_matrix[i] = base[i];

        for (int j = 0; j < i; j++)
        {
            RR product, square_bj, mu;
            product = base[i] * GSO_matrix[j];
            square_bj = GSO_matrix[j] * GSO_matrix[j];
            div(mu, product, square_bj);
            GSO_mu_matrix[i][j] = mu;
            GSO_matrix[i] = GSO_matrix[i] - mu * GSO_matrix[j];
        }
    }
    transpose(GSO_matrix, GSO_matrix);
}

/// @brief 数え上げ上界列を求める関数
/// @param upper_bound 数え上げ上界列の二乗を格納する配列
/// @param epsilon パラメータ（基本0.99）
/// @param B （GSOベクトルの二乗ノルム）
void cal_upper_bound(vec_RR &upper_bound, RR epsilon, vec_RR B)
{
    upper_bound[dimension - 1] = epsilon * B[0];
    for (int i = dimension - 2; i >= 0; i--)
    {
        upper_bound[i] = i * upper_bound[dimension - 1] / dimension;
    }
}

/// @brief 格子上の最短ベクトルの数え上げを行う関数
/// @param base N次元基底ベクトル
/// @param GSO_matrix 基底ベクトルのGSO係数
/// @param GSO_norm GSOベクトルの二乗ノルム
/// @param upper_bound 数え上げ上界列
/// @return 格子ベクトルの係数ベクトル
vec_RR ENUM(mat_RR &base, mat_RR &GSO_matrix, vec_RR &GSO_norm, vec_RR &upper_bound)
{
    ///////// Initialize ////////////
    // mu,vのシグマを格納する用の多次元配列(初期値は全て0)
    mat_RR sigma;
    sigma.SetDims(dimension + 1, dimension);
    clear(sigma);

    // こいつ何者だ
    vec_RR r;
    r.SetLength(dimension + 1);
    for (int i = 0; i < dimension + 1; i++)
        r.at(i) = i;

    // 直交射影を代入する用の配列（初期値は全て0）
    vec_RR rho;
    rho.SetLength(dimension + 1);
    clear(rho);

    // 格子ベクトルの係数ベクトル
    vec_RR v;
    v.SetLength(dimension + 1);
    clear(v);
    v[0] = 1;

    // 負のmu,vのシグマを格納する用の多次元配列(初期値は全て0)
    vec_RR c;
    c.SetLength(dimension);
    clear(c);

    // 重み!バッファが大きい順に深さ優先探索
    vec_RR w;
    w.SetLength(dimension);
    clear(w);

    // v(i) != 0 となる最大の1 <= i <= n

    int last_nonzero = 1;

    int k = 1;
    while (true)
    {
        RR vk_ck;
        mul(vk_ck, sqr(v[k] - c[k]), GSO_norm[k]);
        rho[k] = rho[k + 1] + vk_ck;

        if (rho[k] <= upper_bound[dimension + 1 - k])
        {
            if (k == 1)
            {
                return v;
            }
            else
            {
                k = k - 1;
                r[k - 1] = max(r[k - 1], r[k]);
                for (int i = conv<int>(r[k]) - 1; i >= k + 1; i--)
                {
                    sigma[i][k] = sigma[i + 1][k] + GSO_matrix[i][k] * v[i];
                } // end for
                c[k] = -sigma[k + 1][k];
                v[k] = round(c[k]);
                w[k] = 1;
            } //  end if
        }
        else
        {
            k = k + 1;
            if (k == dimension + 1)
            {
                cout << "Not Exists" << endl;
                return v;
            } // end if
            r[k - 1] = k;
            if (k >= last_nonzero)
            {
                last_nonzero = k;
                v[k] = v[k] - 1;
            }
            else
            {
                if (v[k] > c[k])
                {
                    v[k] = v[k] - w[k];
                }
                else
                {
                    v[k] = v[k] + w[k];
                } // end if
            }     // end if
            w[k] = w[k] + 1;
        } // end if
    }     // end if
} // end while

int main()
{
    // 基底行列
    mat_RR base;
    // グラムシュミット直行化行列
    mat_RR GSO_matrix;
    // グラムシュミット直行化係数行列
    mat_RR GSO_mu_matrix;

    vec_RR B;
    RR epsilon = conv<RR>(0.99);
    vec_RR upper_bound;
    vec_RR v;

    // 次元の設定
    base.SetDims(dimension, dimension);
    GSO_matrix.SetDims(dimension, dimension);
    GSO_mu_matrix.SetDims(dimension, dimension);
    B.SetLength(dimension);
    upper_bound.SetLength(dimension);
    v.SetLength(dimension);
    // 入力ファイル
    string base_file_path = "./BaseMatrix/d=" + to_string(dimension) + "baseMatrix.csv";
    ifstream file(base_file_path);
    string line;

    if (file.fail())
    {
        cout << "Failed to open file." << endl;
        return -1;
    }
    else
    {
        cout << "open successfully" << endl;
    }

    int col = 0;
    while (getline(file, line))
    {
        string value;
        istringstream stream(line);
        int row = 0;
        while (getline(stream, value, ','))
        {
            base[row][col] = stoi(value);
            row++;
        }
        col++;
    }
    transpose(base, base);
    GramSchmidt(base, GSO_matrix, GSO_mu_matrix);

    for (int i = 0; i < dimension; i++)
    {
        B[i] = GSO_matrix[i] * GSO_matrix[i];
    }
    cal_upper_bound(upper_bound, epsilon, B);
    v = ENUM(base, GSO_matrix, B, upper_bound);

    mat_ZZ BKZ_base = conv<mat_ZZ>(base);

    BKZ_FP(BKZ_base, 0.99, dimension);
    for (int i = 0; i < v.length(); i++)
    {
        cout << v[i] << " ";
    }
    return 0;
}
NTL_CLIENT