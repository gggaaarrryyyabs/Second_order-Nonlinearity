#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <iostream>
#include <functional>
#include <cstdint>

using namespace std;

// g:=x3*x2*x1 + x4 + x4*x1 + x4*x2 + x4*x2*x1 + x4*x3*x1 + x4*x3*x2 +x5 + x5*x1 + x5*x2*x1 + x5*x3 + x5*x3*x1 + x5*x3*x2 + x5*x4 +x5*x4*x1 + x5*x4*x2 + x5*x4*x3;

void WalshTransform(int *tt, const int N)
{
    for (int i = 0; i < N; i++)
    {
        int gap = 1 << (N - i - 1);
        for (int j = 0; j < (2 << i); j = j + 2)
        {
            for (int l = 0; l < gap; l++)
            {
                tt[gap * j + l] = tt[gap * j + l] + tt[gap * (j + 1) + l];
                tt[gap * (j + 1) + l] = tt[gap * j + l] - 2 * tt[gap * (j + 1) + l];
            }
        }
    }
}

int hw(int n)
{
    int sum = 0;
    while (n != 0)
    {
        if (n & 1)
            sum++;
        n = n >> 1;
    }
    return sum;
}

void FMEI(int *tt, int n_var)
{
    // float will cause error since the value of influence is small
    double influence = 0;
    float min_entropy = 0;

    int max = 0;
    for (auto i = 0; i < 1 << n_var; i++)
        tt[i] = 1 - 2 * tt[i];

    WalshTransform(tt, n_var);

    for (auto i = 0; i < 1 << n_var; i++)
    {
        influence += hw(i) * pow((double)tt[i] / (1 << n_var), 2);
        if (abs(tt[i]) > max)
            max = abs(tt[i]);
    }
    cout << "max walsh spectrum:" << max << endl;

    min_entropy = 2 * n_var - 2 * log2((float)max);

    // min_entropy is 12
    cout << "Min Entropy: " << min_entropy << endl;

    // Influence should be 9/4 * 15/8 = 2.8125
    cout << "Influence: " << influence << endl;

    cout << "FMEI: " << min_entropy / influence << endl;
}

int main()
{
    int n_var = 30;

    // g:=x3*x2*x1 + x4 + x4*x1 + x4*x2 + x4*x2*x1 + x4*x3*x1 + x4*x3*x2 +x5 + x5*x1 + x5*x2*x1 + x5*x3 + x5*x3*x1 + x5*x3*x2 + x5*x4 +x5*x4*x1 + x5*x4*x2 + x5*x4*x3;
    int tt_g[32];
    for (auto i = 0; i < 32; i++)
    {
        int x[5];
        for (auto j = 0; j < 5; j++)
            x[j] = (i >> j) & 1;

        tt_g[16 * x[0] + 8 * x[1] + 4 * x[2] + 2 * x[3] + x[4]] = (x[4] * x[3] * x[2]) ^ x[1] ^ (x[1] * x[4]) ^ (x[1] * x[3]) ^ (x[1] * x[3] * x[4]) ^ (x[1] * x[2] * x[4]) ^ (x[1] * x[2] * x[3]) ^ x[0] ^ (x[0] * x[4]) ^ (x[0] * x[3] * x[4]) ^ (x[0] * x[2]) ^ (x[0] * x[2] * x[4]) ^ (x[0] * x[2] * x[3]) ^ (x[0] * x[1]) ^ (x[0] * x[1] * x[4]) ^ (x[0] * x[1] * x[3]) ^ (x[0] * x[1] * x[2]);
    }

    // g||g_inverse = g(x1,x2,x3,x4,x5)||g(x1+1,x2+1,x3+1,x4+1,x5+1)
    int tt_g_g_inverse[64];
    for (auto i = 0; i < 32; i++)
    {
        tt_g_g_inverse[i] = tt_g[i];
        tt_g_g_inverse[i + 32] = tt_g[31 - i];
    }

    // g||g_inverse(g(x1),g(x2),g(x3),g(x4),g(x5),g(x6))
    int *tt = new int[1 << 30];
    for (auto i = 0; i < 1 << 30; i++)
    {
        int i_binary[30];
        for (auto j = 0; j < 30; j++)
            i_binary[j] = (i >> j) & 1;

        int index_final[6];
        for (auto j = 0; j < 6; j++)
            index_final[j] = 16 * i_binary[j * 5] + 8 * i_binary[j * 5 + 1] + 4 * i_binary[j * 5 + 2] + 2 * i_binary[j * 5 + 3] + 1 * i_binary[j * 5 + 4];

        tt[i] = tt_g_g_inverse[32 * tt_g[index_final[0]] + 16 * tt_g[index_final[1]] + 8 * tt_g[index_final[2]] + 4 * tt_g[index_final[3]] + 2 * tt_g[index_final[4]] + tt_g[index_final[5]]];
    }

    FMEI(tt, 30);

    vector<int> walsh_spectrum = {0};
    for (auto i = 0; i < 1 << 30; i++)
    {
        if (find(walsh_spectrum.begin(), walsh_spectrum.end(), abs(tt[i])) == walsh_spectrum.end())
            walsh_spectrum.push_back(abs(tt[i]));
    }

    cout << "walsh spectrum size: " << walsh_spectrum.size() << endl;

    cout << "walsh spectrum: ";
    for (auto &i : walsh_spectrum)
        cout << i << " ";
    cout << endl;

    vector<set<int>> walsh_spectrum_hw;

    for (auto i = 0; i < walsh_spectrum.size(); i++)
        walsh_spectrum_hw.push_back(set<int>());

    for (auto j = 0; j < 1 << 30; j++)
    {
        for (auto i = 0; i < walsh_spectrum.size(); i++)
        {
            if (abs(tt[j]) == walsh_spectrum[i])
            {
                walsh_spectrum_hw[i].insert(hw(j));
            }
        }
    }

    cout << "walsh spectrum hw size: " << walsh_spectrum_hw.size() << endl;

    for (auto i = 0; i < walsh_spectrum_hw.size(); i++)
    {
        for (auto &j : walsh_spectrum_hw[i])
        {
            cout << j << " ";
        }
        cout << endl;
    }

    delete[] tt;
}

// 0: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
// 16777216: 2 3 4 6 7 10
// 65536: 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 26 27 30
