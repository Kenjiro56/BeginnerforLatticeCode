// 出入力ライブラリ
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// NTLライブラリ
#include <NTL/ZZ_p.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/LLL.h> // LLLを使うためのライブラリ
using namespace std;
using namespace NTL;
#define dimension 4

int main()
{

    // 基底行列
    mat_ZZ base;

    // 次元の設定
    base.SetDims(dimension, dimension);

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

    // 浮動小数点LLL
    LLL_FP(base);

    // 出力
    cout << base << endl;
}
NTL_CLIENT