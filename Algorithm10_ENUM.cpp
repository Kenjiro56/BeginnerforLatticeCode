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
        upper_bound[i] = (i + 1) * upper_bound[dimension - 1] / dimension;
        // cout << i << "に" << upper_bound[i] << "を入れるよ" << endl;
    }
}

/// @brief 格子上の最短ベクトルの数え上げを行う関数
/// @param base N次元基底ベクトル
/// @param GSO_matrix 基底ベクトルのGSO係数
/// @param GSO_norm GSOベクトルの二乗ノルム
/// @param upper_bound 数え上げ上界列
/// @return 格子ベクトルの係数ベクトル
vec_RR ENUM(mat_RR &base, mat_RR &GSO_mu_matrix, vec_RR &GSO_norm, vec_RR &upper_bound)
{
    // transpose(GSO_matrix, GSO_matrix);
    // cout << "グラムシュミット行列" << endl;
    // cout << GSO_matrix << "\t";
    // cout << "" << endl;

    // cout << "B" << endl;
    // cout << GSO_norm << "\t";
    // cout << "" << endl;

    // cout << "R^2" << endl;
    // cout << upper_bound << "\t";
    // cout << "" << endl;

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
    v.SetLength(dimension);
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
        mul(vk_ck, sqr(v[k - 1] - c[k - 1]), GSO_norm[k - 1]);
        rho[k - 1] = rho[k] + vk_ck;
        // cout << rho[k - 1] << endl;
        // cout << upper_bound[dimension - k] << endl;
        // cout << "k:" << k << endl;
        if (rho[k - 1] <= upper_bound[dimension - k])
        {
            if (k == 1)
            {
                cout << "1:お返しします" << endl;
                return v;
            }
            else
            {
                cout << "2:" << endl;
                k = k - 1;
                r[k - 1] = max(r[k - 1], r[k]);
                // cout << conv<int>(r[k - 1]) << "から " << k + 1 << "まで" << endl;

                for (int i = conv<int>(r[k - 1]); i >= k + 1; i--)
                {
                    // cout << "格納してるにょん" << endl;
                    sigma[i - 1][k - 1] = sigma[i][k - 1] + GSO_mu_matrix[i - 1][k - 1] * v[k - 1];
                } // end for
                c[k - 1] = -sigma[k][k - 1];
                v[k - 1] = round(c[k - 1]);
                w[k - 1] = 1;
            } //  end if
        }
        else
        {
            k = k + 1;
            if (k == dimension + 1)
            {
                cout << "3:Not Exists" << endl;
                return v;
            } // end if
            r[k - 1] = k;
            cout << "4:" << endl;
            if (k >= last_nonzero)
            {
                last_nonzero = k;
                v[k - 1] = v[k - 1] + 1;
            }
            else
            {
                if (v[k - 1] > c[k - 1])
                {
                    v[k - 1] = v[k - 1] - w[k - 1];
                }
                else
                {
                    v[k - 1] = v[k - 1] + w[k - 1];
                } // end if
                w[k - 1] = w[k - 1] + 1;
            } // end if
        }     // end if
        cout << "k = " << k << endl;
        cout << "sigma = " << endl;
        cout << sigma << endl;
        cout << "r = " << r << endl;
        cout << "rho = " << rho << endl;
        cout << "v = " << v << endl;
        cout << "c = " << c << endl;
        cout << "w = " << w << endl;
        cout << "last_nonzero = " << last_nonzero << endl;
        cout << "----------------------------------------" << endl;
    } // end while
} // end function
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
    // 配列出力
    // for (int i = 0; i < dimension; i++)
    // {
    //     for (int j = 0; j < dimension; j++)
    //     {
    //         cout << base[i][j] << " ";
    //     }
    //     cout << "" << endl;
    // }
    // こいつ使ってない
    // mat_RR mu;
    // mu.SetDims(dimension, dimension);

    // 下の関数で代替（ntlのライブラリ）
    // GramSchmidt(base, GSO_matrix, GSO_mu_matrix);
    // 行列が転置するから
    transpose(base, base);
    ComputeGS(conv<mat_ZZ>(base), GSO_mu_matrix, B);
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            if (i == j)
            {
                GSO_mu_matrix[i][j] = 1;
            }
        }
    }
    /////// 入力パラメータの確認 //////////////////
    cout << "------------------入力パラメータ----------------" << endl;
    cout << "base" << endl;
    cout << base << endl;
    cout << "mu" << endl;
    cout << GSO_mu_matrix << endl;
    cout << "B" << endl;
    cout << B << endl;

    // 上の関数で代替できる
    // for (int i = 0; i < dimension; i++)
    // {
    //     transpose(GSO_matrix, GSO_matrix);
    //     //  これどっちも同じ
    //     InnerProduct(B[i], GSO_matrix[i], GSO_matrix[i]);
    //     // B[i] = GSO_matrix[i] * GSO_matrix[i];
    // }
    cal_upper_bound(upper_bound, epsilon, B);
    cout << "R^2" << endl;
    cout << upper_bound << endl;
    cout << "----------------------------------------------" << endl;
    // 数え上げ上界がうまく計算できているか確認
    // for (int i = 0; i < dimension; i++)
    // {
    //     cout << upper_bound.at(i) << endl;
    // }
    // cout << GSO_matrix << endl;
    // cout << B << endl;
    // cout << upper_bound << endl;
    v = ENUM(base, GSO_mu_matrix, B, upper_bound);

    // mat_ZZ BKZ_base = conv<mat_ZZ>(base);

    // BKZ_FP(BKZ_base, 0.99, dimension);
    cout << "最終出力結果" << endl;
    cout << v << endl;
    cout << "----------------------------------" << endl;
    cout << "最終出力結果によるSVP" << endl;
    vec_RR myself_SVP;
    myself_SVP.SetLength(dimension);
    clear(myself_SVP);
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            // cout << myself_SVP[i] << " + " << v[j] << " * " << base[j][i] << endl;
            myself_SVP[i] += v[j] * base[j][i];
        }
        // cout << "---" << endl;
    }
    cout << myself_SVP << endl;

    cout << "ライブラリによる出力" << endl;
    mat_ZZ lib_SVP = conv<mat_ZZ>(base);
    LLL_FP(lib_SVP);
    cout << lib_SVP[0] << endl;
    return 0;
}
NTL_CLIENT