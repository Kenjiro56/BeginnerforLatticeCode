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
#define mindimenstion 40
#define maxdimension 50
#define trytime 100

int main()
{
    for (int i = mindimenstion; i <= maxdimension; i++)
    {
        // 基底行列
        mat_ZZ base;

        // 次元の設定
        base.SetDims(i, i);

        // 入力ファイル
        // string base_file_path = "./BaseMatrix/d=" + to_string(i) + "baseMatrix.csv";
        string base_file_path = "./BaseMatrix/d=" + to_string(i) + "baseMatrix.txt";
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
                cout << value << endl;
                cout << "uiui" << endl;
                base[row][col] = stoi(value);
                row++;
            }
            col++;
        }
        transpose(base, base);

        int successtime = 0;
        // 浮動小数点LLL
        mat_ZZ LLLbase = base, BKZbase = base;
        for (int j = 0; j < trytime; j++)
        {
            LLL_FP(LLLbase);
            BKZ_FP(BKZbase, 0.99, i);
            if (LLLbase == BKZbase)
            {
                successtime++;
            }
        }
        cout << "dimsntion:" << i << "," << successtime / trytime << endl;
    }
}
NTL_CLIENT