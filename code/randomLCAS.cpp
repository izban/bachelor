#include <bits/stdc++.h>
using namespace std;
const int INF = (int)1e9;

string randomString(int n) {
    string res = "";
    for (int i = 0; i < n; i++) res += (char)('a' + rand() % 2);
    return res;
}

vector<pair<int, int> > abelianSubstrings(string s) {
    int n = s.size();
    vector<pair<int, int> > res(n + 1, make_pair(INF, -INF));
    res[0] = make_pair(0, 0);
    for (int i = 0; i < (int)s.length(); i++) {
        int cur = 0;
        for (int j = i; j < (int)s.length(); j++) {
            cur += s[j] == 'b';
            res[j - i + 1].first = min(res[j - i + 1].first, cur);
            res[j - i + 1].second = max(res[j - i + 1].second, cur);
        }
    }
    return res;
}

int lcas(string a, string b) {
    auto ca = abelianSubstrings(a);
    auto cb = abelianSubstrings(b);
    //int n = a.length();
    //cerr << ca[n].first << " " << ca[n].second << " " << cb[n].first << " " << cb[n].second << endl;
    for (int i = min(a.length(), b.length()); i >= 0; i--) {
        if (max(ca[i].first, cb[i].first) <= min(ca[i].second, cb[i].second)) {
            return i;
        }
    }
    assert(0);
}

pair<int, int> getMinMaxSubstr(string &a, int len) {
    pair<int, int> res = make_pair(INF, -INF);
    if (len > (int)a.length()) return res;
    int cur = 0;
    for (int i = 0; i < len - 1; i++) cur += a[i] == 'b';
    for (int i = 0; i + len <= (int)a.length(); i++) {
        cur += a[i + len - 1] == 'b';
        res.first = min(res.first, cur);
        res.second = max(res.second, cur);
        cur -= a[i] == 'b';
    }
    return res;
}

int jumps;
int lcasFast(string a, string b) {
    int len = a.length();
    jumps = 0;
    while (1) {
        pair<int, int> ca = getMinMaxSubstr(a, len);
        pair<int, int> cb = getMinMaxSubstr(b, len);
        if (ca > cb) swap(ca, cb);
        if (max(ca.first, cb.first) <= min(ca.second, cb.second)) {
            return len;
        }
        len -= (cb.first - ca.second + 1) / 2;
        jumps++;
    }
    assert(0);
}

int main() {
    //const int N = 10000;
    for (int N = 500; N <= 10000; N += 500) {
    long long sum = 0;
    long long sumJumps = 0;
    const int ITERATIONS = 5e3;
    for (int it = 0; it < ITERATIONS; it++) {
        srand(it);
        string a = randomString(N);
        string b = randomString(N);
        int len = lcasFast(a, b);
        sum += len;
        sumJumps += jumps;
        //cout << it << " " << len << "/" << N << ", jumps: " << jumps << ", avLCAS: " << (double)sum / ((it + 1)) << ", av jumps: " << (double)sumJumps / (it + 1) << endl;
    } 
    cout << N << " " << sum / ITERATIONS << endl;
    }
}
