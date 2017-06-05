#include <bits/stdc++.h>
using namespace std;

const int N = 100;

int d[N][N];

int main() {
    int n = 5;
    d[0][N + 0] = 1;
    for (int i = 0; i < n; i++) {
        for (int j = -i; j <= i; j++) {
            cerr << d[i][N + j] << " ";
            d[i + 1][N + j - 1] += d[i][N + j];
            d[i + 1][N + j] += 2 * d[i][N + j];
            d[i + 1][N + j + 1] += d[i][N + j];
        }
        cerr << endl;
    }
}
