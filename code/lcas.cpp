#include <bits/stdc++.h>
using namespace std;
const int INF = (int)1e9;

namespace slow {
    int solve(vector<int> a, vector<int> b) {
        assert(a.size() == b.size());
        int n = (int)a.size();
        int k = max(*max_element(a.begin(), a.end()), *max_element(b.begin(), b.end())) + 1;
        int ans = 0;
        for (int len = 1; len <= n; len++) {
            vector<vector<int> > aa;
            for (int i = 0; i + len <= n; i++) {
                vector<int> cur(k);
                for (int j = 0; j < len; j++) {
                    cur[a[i + j]]++;
                }
                aa.push_back(cur);
            }
            sort(aa.begin(), aa.end());
            for (int i = 0; i + len <= n; i++) {
                vector<int> cur(k);
                for (int j = 0; j < len; j++) {
                    cur[b[i + j]]++;
                }
                if (binary_search(aa.begin(), aa.end(), cur)) {
                    ans = max(ans, len);
                }
            }
        }
        return ans;
    }
}

namespace attabi {
    int solve(vector<int> a, vector<int> b) {
        assert(a.size() == b.size());
        int n = (int)a.size();
        int k = max(*max_element(a.begin(), a.end()), *max_element(b.begin(), b.end())) + 1;
        int ans = 0;
        for (int len = 1; len <= n; len++) {
            vector<vector<int> > aa, bb;
            for (int i = 0; i + len <= n; i++) {
                vector<int> cur(k);
                for (int j = 0; j < len; j++) {
                    cur[a[i + j]]++;
                }
                aa.push_back(cur);
            }
            int brd = aa.size();
            for (int i = 0; i + len <= n; i++) {
                vector<int> cur(k);
                for (int j = 0; j < len; j++) {
                    cur[b[i + j]]++;
                }
                aa.push_back(cur);
            }
            int m = aa.size();
            vector<int> c(m);
            vector<int> p(m);
            for (int i = 0; i < m; i++) p[i] = i;
            int cl = 1;
            for (int i = k - 1; i >= 0; i--) {
                vector<int> cnt(n + 1);
                for (int j = 0; j < m; j++) cnt[aa[j][i]]++;
                for (int j = 1; j <= n; j++) cnt[j] += cnt[j - 1];
                vector<int> pn = p;
                for (int j = m - 1; j >= 0; j--) p[--cnt[aa[pn[j]][i]]] = pn[j];

                cl = 1;
                vector<int> cn(m);
                cn[p[0]] = cl - 1;
                for (int j = 1; j < m; j++) {
                    cl += c[p[j]] != c[p[j - 1]] || aa[p[j]][i] != aa[p[j - 1]][i];
                    cn[p[j]] = cl - 1;
                }
                c = cn;
            }
            vector<int> mask(cl);
            for (int i = 0; i < m; i++) {
                if (i < brd) mask[c[i]] |= 1;
                else mask[c[i]] |= 2;
            }
            for (int i = 0; i < cl; i++) if (mask[i] == 3) ans = max(ans, len);
        }
        return ans;
    }
}

namespace fast {
    int curNode;

    struct node {
        node *ch[2];
        node *fat;
        node *pr;
        vector<pair<pair<node*, int>, pair<int, int> > > ppr;
        int which, time;
        int h;
        int pos, val; // for leafs
        int level;
        vector<pair<int, int> > vct;
        int created;

        vector<vector<pair<int, int> > > hson;

        node() {}

        void updPar(node *par, int son) {
            if (!ppr.empty()) {
                ppr.back().second.second = curNode;
            }
            ppr.push_back(make_pair(make_pair(par, son), make_pair(curNode, INF)));
            pr = par;
        }
    };

    template<typename T>
    struct clearingArray {
        int n;
        int time;
        vector<int> lastCleared;
        vector<T> data;

        clearingArray(int nn) {
            n = nn;
            time = 0;
            lastCleared.assign(n, 0);
            data.assign(n, T());
        }

        T& operator[](int i) {
            assert(0 <= i && i < n);
            if (lastCleared[i] < time) {
                lastCleared[i] = time;
                data[i] = T();
            }
            return data[i];
        }

        void clear() {
            time++;
        }
    };

    const int MAXN = 1 << 17;
    node nodes[MAXN];

    node* newNode(node *l, node *r) {
        assert(curNode < MAXN);
        node *t = &nodes[curNode++];
        t->ch[0] = l;
        t->ch[1] = r;
        l->updPar(t, 0);
        r->updPar(t, 1);
        t->fat = 0;
        t->which = -1;
        t->h = -1;
        t->pos = t->val = 0;
        //assert(l->level == r->level);
        //level = l->level + 1;
        t->level = max(l->level, r->level) + 1;
        t->pr = 0;
        t->time = 0;
        t->created = curNode;
        t->hson.resize(2);
        return t;
    }
    node* newNode(int npos, int nval) {
        assert(curNode < MAXN);
        node *t = &nodes[curNode++];
        t->ch[0] = 0;
        t->ch[1] = 0;
        t->fat = 0;
        t->which = -1;
        t->h = -1;
        t->pos = npos;
        t->val = nval;
        t->level = 0;
        t->pr = 0;
        t->time = 0;
        t->created = curNode;
        t->hson.resize(2);
        return t;
    }

    int solve(vector<int> a, vector<int> b) {
        assert(a.size() == b.size());
        int n = a.size();
        int k = max(*max_element(a.begin(), a.end()), *max_element(b.begin(), b.end())) + 1;
        if (k == 1) return n;
        //int kk = 1;
        //while (kk < k) kk <<= 1;
        //k = kk;

        int ans = 0;
        double sum = 0;
        auto gettime = [&]() {
            return (double)clock() / (double)CLOCKS_PER_SEC;
        };
        for (int len = 1; len <= n; len++) {
            double start = gettime();
            curNode = 0;

            vector<node*> last(k);

            function<node*(int, int)> build = [&](int tl, int tr) {
                if (tl == tr) {
                    auto t = newNode(tl, 0);
                    last[tl] = t;
                    return t;
                }
                int tm = (tl + tr) >> 1;
                auto bl = build(tl, tm);
                auto br = build(tm + 1, tr);
                auto t = newNode(bl, br);
                bl->ppr.back().second.first = br->ppr.back().second.first = 0;
                return t;
            };
            auto croot = build(0, k - 1);
            for (int i = 0; i < curNode; i++) nodes[i].created = 0;
            int stime = 0;
            auto add = [&](int x, int dx) {
                auto prev = last[x];
                auto cur = newNode(x, last[x]->val + dx);
                stime = cur->created;
                last[x] = cur;
                while (1) {
                    if (prev->pr == 0) return cur;
                    auto par = prev->pr;
                    node *nch[2] = {par->ch[0], par->ch[1]};
                    int nwhich = -1;
                    if (prev == par->ch[0]) nwhich = 0;
                    else if (prev == par->ch[1]) nwhich = 1;
                    else if (prev == par->fat) nwhich = par->which;
                    else assert(0);
                    if (par->which == -1) {
                        par->which = nwhich;
                        par->fat = cur;
                        par->time = curNode;
                        cur->updPar(par, nwhich);
                        cur->ppr.back().second.first = stime;
                        return croot;
                    } else {
                        nch[par->which] = par->fat;
                        nch[nwhich] = cur;
                        auto npar = newNode(nch[0], nch[1]);
                        npar->created = stime;
                        nch[0]->ppr.back().second.first = nch[1]->ppr.back().second.first = stime;
                        cur->updPar(npar, nwhich);
                        cur->ppr.back().second.first = stime;
                        prev = par;
                        cur = npar;
                    }
                }
            };
            int frst = -1;
            vector<int> from, to;
            vector<int> badTimes;
            for (vector<int> s : {a, b}) {
                for (int i = 0; i < len; i++) {
                    badTimes.push_back(stime);
                    croot = add(s[i], +1);
                }
                from.push_back(stime);
                for (int i = 1; i + len <= n; i++) {
                    croot = add(s[i - 1], -1);
                    badTimes.push_back(stime);
                    croot = add(s[i + len - 1], +1);
                }
                to.push_back(stime);
                for (int i = n - len; i < n; i++) {
                    croot = add(s[i], -1);
                    badTimes.push_back(stime);
                }
                if (frst == -1) frst = curNode;
            }
            vector<char> badTime(curNode + 1);
            for (int x : badTimes) badTime[x] = 1;

            int height = croot->level + 1;
            vector<vector<node*> > vert(height);
            for (int i = 0; i < curNode; i++) {
                vert[nodes[i].level].push_back(&nodes[i]);
            }

            auto passToParent = [&](node *v) {
                for (int i = 0; i + 1 < (int)v->ppr.size(); i++) v->ppr[i].second.second = min(v->ppr[i].second.second, v->ppr[i + 1].second.first);
                int i = 0;
                for (int curP = 0; curP < (int)v->ppr.size(); curP++) {
                    while (i < (int)v->vct.size() && v->vct[i].first <= v->ppr[curP].second.second) {
                        auto o = v->vct[i];
                        o.first = max(o.first, v->ppr[curP].second.first);
                        v->ppr[curP].first.first->hson[v->ppr[curP].first.second].push_back(o);
                        i++;
                    }
                    i--;
                }
                v->hson[0].clear();
                v->hson[1].clear();
                v->vct.clear();
                v->ppr.clear();
            };
            {
                vector<int> mx(k);
                for (auto v : vert[0]) {
                    mx[v->pos] = max(mx[v->pos], v->val + 1);
                }
                vector<int> st(k);
                for (int i = 1; i < k; i++) {
                    st[i] = st[i - 1] + mx[i - 1];
                }

                for (auto v : vert[0]) {
                    v->h = st[v->pos] + v->val;
                    v->vct = vector<pair<int, int> >(1, make_pair(v->created, v->h));
                    passToParent(v);
                }
            }
            sum += gettime() - start;
            clearingArray<int> arr(curNode);
            for (int lvl = 1; lvl < height; lvl++) {
                vector<vector<pair<int, int*> > > calc(curNode);
                for (auto v : vert[lvl]) {
                    for (int i = 0; i < 2; i++) {
                        vector<vector<pair<int, int> > > vv(1);
                        for (int j = 0; j < (int)v->hson[i].size(); j++) {
                            if (j > 0 && v->hson[i][j].first < v->hson[i][j - 1].first) {
                                vv.push_back(vector<pair<int, int> >());
                            }
                            vv.back().push_back(v->hson[i][j]);
                        }
                        assert(1 <= vv.size() && vv.size() <= 2);
                        if (vv.size() == 2) {
                            vector<pair<int, int> > nv(vv[0].size() + vv[1].size());
                            merge(vv[0].begin(), vv[0].end(), vv[1].begin(), vv[1].end(), nv.begin());
                            vv.resize(1);
                            vv[0] = nv;
                        }
                        v->hson[i] = vv[0];
                        if (v->hson[i].empty()) {
                            assert(!v->hson[i].empty());
                        }
                    }
                    int ci = 0, cj = 0;
                    vector<pair<int, int> > oo;
                    while (1) {
                        int ctime = max(v->hson[0][ci].first, v->hson[1][cj].first);
                        if (ctime >= v->created) {
                            v->vct.push_back(make_pair(ctime, -1));
                            oo.push_back({ci, cj});
                        } else {
                            ctime = ctime;
                        }
                        if (ci + 1 < (int)v->hson[0].size() && (cj + 1 == (int)v->hson[1].size() || v->hson[0][ci + 1].first <= v->hson[1][cj + 1].first)) {
                            ci++;
                        } else if (cj + 1 < (int)v->hson[1].size() && (ci + 1 == (int)v->hson[0].size() || v->hson[0][ci + 1].first > v->hson[1][cj + 1].first)) {
                            cj++;
                        } else break;
                    }
                    for (int i = 0; i < (int)v->vct.size(); i++) {
                        calc[v->hson[0][oo[i].first].second].push_back(make_pair(v->hson[1][oo[i].second].second, &v->vct[i].second));
                    }
                }

                int curState = 0;
                for (int i = 0; i < curNode; i++) {
                    arr.clear();
                    for (auto o : calc[i]) {
                        if (arr[o.first] == 0) {
                            arr[o.first] = ++curState;
                        }
                        *o.second = arr[o.first];
                    }
                }
                if (lvl + 1 < height) {
                    for (auto v : vert[lvl]) {
                        passToParent(v);
                    }
                } else {
                    arr.clear();
                    for (auto v : vert[lvl]) {
                        for (auto o : v->vct) if (!badTime[o.first]) {
                            if (o.first >= from[0] && o.first <= to[0]) {
                                arr[o.second] |= 1;
                            }
                            if (o.first >= from[1] && o.first <= to[1]) {
                                arr[o.second] |= 2;
                            }
                        }
                        passToParent(v);
                    }
                    for (int i = 0; i < curNode; i++) {
                        if (arr[i] == 3) {
                            ans = max(ans, len);
                        }
                    }
                }
            }
        }
        return ans;
    }
}

auto randString = [&](int n, int k) {
    vector<int> res;
    for (int i = 0; i < n; i++) {
        res.push_back(rand() % k);
        //res.push_back(k - 1 - rand() % 3);
    }
    return res;
};

void stress() {
    for (int it = 100;; it++) {
        cerr << it << endl;
        srand(it);

        //int n = rand() % 1000 + 1;
        //int k = rand() % n + 1;
        int n = 1600;
        int k = n;
        vector<int> a = randString(n, k);
        vector<int> b = randString(n, k);

        auto print = [&]() {
            cerr << n << endl;
            for (int i = 0; i < n; i++) cerr << a[i] << " \n"[i + 1 == n];
            for (int i = 0; i < n; i++) cerr << b[i] << " \n"[i + 1 == n];
        };


        double st1 = clock() / (double)CLOCKS_PER_SEC;
        int ans1 = fast::solve(a, b);
        double fn1 = clock() / (double)CLOCKS_PER_SEC;

        double st2 = clock() / (double)CLOCKS_PER_SEC;
        int ans2 = ans1;//attabi::solve(a, b);
        double fn2 = clock() / (double)CLOCKS_PER_SEC;

        double st3 = clock() / (double)CLOCKS_PER_SEC;
        int ans3 = slow::solve(a, b);
        double fn3 = clock() / (double)CLOCKS_PER_SEC;

        if (1 && (ans1 != ans2 || ans1 != ans3)) {
            cerr << ans1 << " instead of " << ans2 << endl;
            print();
            exit(0);
        }
        cerr << "fast: " << fn1 - st1 << endl;
        cerr << "attabi: " << fn2 - st2 << endl;
        cerr << "slow: " << fn3 - st3 << endl;
    }
}

template <typename T>
void testTime(T f, int N) {
    for (int n = N; n <= N; n += 100) {
        int k = min(n, (int)(log(n) * 50));
        cerr << n << " " << k << endl;
        auto a = randString(n, k);
        auto b = randString(n, k);
        
        double st = clock() / (double)CLOCKS_PER_SEC;
        f(a, b);
        double fn = clock() / (double)CLOCKS_PER_SEC;
        cout << n << " " << fn - st << endl;
    } 
}

int main() {
    //stress();
    for (int N = 600; N <= 1500; N += 100) { 
        testTime([&](vector<int> a, vector<int> b) {
            fast::solve(a, b);
        }, N);
        testTime([&](vector<int> a, vector<int> b) {
            attabi::solve(a, b);
        }, N);
    }
    if (0) for (int N = 600; N <= 1500; N += 100) { 
    }
    return 0;
}
