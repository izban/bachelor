#include <sstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cassert>
#include <ctime>
#include <map>
#include <math.h>
#include <cstdio>
#include <set>
#include <deque>
#include <memory.h>
#include <queue>

//#pragma comment(linker, "/STACK:64000000")
typedef long long ll;

using namespace std;

namespace FFT
{
	const int maxBase = 21;
	const int maxN = 1 << maxBase;
	typedef double dbl;
#define forn(i, n) for (int i = 0; i < (n); i++)
#define pw(n) (1LL << (n))

	struct num
	{
		dbl x, y;
		num() {}
		num(dbl xx, dbl yy) : x(xx), y(yy) {}
		num(dbl alp) : x(cos(alp)), y(sin(alp)) {}
	};

	inline num operator + (num a, num b) { return num(a.x + b.x, a.y + b.y); }
	inline num operator - (num a, num b) { return num(a.x - b.x, a.y - b.y); }
	inline num operator * (num a, num b) { return num(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }
	inline num conj(num a) { return num(a.x, -a.y); }

	const dbl PI = acos(-1);

	num root[maxN];
	int rev[maxN];
	bool rootsPrepared = false;

	void prepRoots()
	{
		if (rootsPrepared) return;
		rootsPrepared = true;
		root[1] = num(1, 0);
		for (int k = 1; k < maxBase; ++k)
		{
			num x(2 * PI / pw(k + 1));
			for (int i = pw(k - 1); i < pw(k); ++i)
			{
				root[2 * i] = root[i];
				root[2 * i + 1] = root[i] * x;
			}
		}
	}

	int base, N;

	int lastRevN = -1;
	void prepRev()
	{
		if (lastRevN == N) return;
		lastRevN = N;
		forn(i, N) rev[i] = (rev[i >> 1] >> 1) + ((i & 1) << (base - 1));
	}

	void fft(num *a, num *f)
	{
		forn(i, N) f[i] = a[rev[i]];
		for (int k = 1; k < N; k <<= 1) for (int i = 0; i < N; i += 2 * k) forn(j, k)
		{
			num z = f[i + j + k] * root[j + k];
			f[i + j + k] = f[i + j] - z;
			f[i + j] = f[i + j] + z;
		}
	}

	num a[maxN], b[maxN], f[maxN], g[maxN];
	ll A[maxN], B[maxN], C[maxN];

	void _multMod(int mod)
	{
		forn(i, N)
		{
			int x = A[i] % mod;
			a[i] = num(x & (pw(15) - 1), x >> 15);
		}
		forn(i, N)
		{
			int x = B[i] % mod;
			b[i] = num(x & (pw(15) - 1), x >> 15);
		}
		fft(a, f);
		fft(b, g);

		forn(i, N)
		{
			int j = (N - i) & (N - 1);
			num a1 = (f[i] + conj(f[j])) * num(0.5, 0);
			num a2 = (f[i] - conj(f[j])) * num(0, -0.5);
			num b1 = (g[i] + conj(g[j])) * num(0.5 / N, 0);
			num b2 = (g[i] - conj(g[j])) * num(0, -0.5 / N);
			a[j] = a1 * b1 + a2 * b2 * num(0, 1);
			b[j] = a1 * b2 + a2 * b1;
		}

		fft(a, f);
		fft(b, g);

		forn(i, N)
		{
			ll aa = f[i].x + 0.5;
			ll bb = g[i].x + 0.5;
			ll cc = f[i].y + 0.5;
			C[i] = (aa + bb % mod * pw(15) + cc % mod * pw(30)) % mod;
		}
	}

	void prepAB(int n1, int n2)
	{
		base = 1;
		N = 2;
		while (N < n1 + n2) base++, N <<= 1;

		for (int i = n1; i < N; ++i) A[i] = 0;
		for (int i = n2; i < N; ++i) B[i] = 0;

		prepRoots();
		prepRev();
	}

	void mult(int n1, int n2)
	{
		prepAB(n1, n2);
		forn(i, N) a[i] = num(A[i], B[i]);
		fft(a, f);
		forn(i, N)
		{
			int j = (N - i)  & (N - 1);
			a[i] = (f[j] * f[j] - conj(f[i] * f[i])) * num(0, -0.25 / N);
		}
		fft(a, f);
		forn(i, N) C[i] = (ll)round(f[i].x);
	}


	void multMod(int n1, int n2, int mod)
	{
		prepAB(n1, n2);
		_multMod(mod);
	}

	int D[maxN];

	void multLL(int n1, int n2)
	{
		prepAB(n1, n2);

		int mod1 = 1.5e9;
		int mod2 = mod1 + 1;

		_multMod(mod1);

		forn(i, N) D[i] = C[i];

		_multMod(mod2);

		forn(i, N)
		{
			C[i] = D[i] + (C[i] - D[i] + (ll)mod2) * (ll)mod1 % mod2 * mod1;
		}
	}
	// HOW TO USE ::
	// -- set correct maxBase
	// -- use mult(n1, n2), multMod(n1, n2, mod) and multLL(n1, n2)
	// -- input  : A[], B[]
	// -- output : C[]
}


struct pt {
	int x, y;

	pt() {}
	pt(int _x, int _y) : x(_x), y(_y) {}

	pt operator+ (const pt &p) const {
		return pt(x + p.x, y + p.y);
	}

	bool operator< (const pt &p) const {
		if (x != p.x) return x < p.x;
		return y < p.y;
	}

	bool operator> (const pt &p) const {
		if (x != p.x) return x > p.x;
		return y > p.y;
	}

	bool operator== (const pt &p) const {
		return x == p.x && y == p.y;
	}
};

const int MX = 1e7;
const int PW = 1 << 22;
char isPr[MX];
vector<int> pr;
bool done = 0;
int bits[PW];

void precalc() {
	if (done) return;
	done = 1;
	for (int i = 0; i < MX; i++) isPr[i] = 1;
	for (int i = 2; i < MX; i++) {
		if (!isPr[i]) continue;
		pr.push_back(i);
		for (int j = i + i; j < MX; j += i) isPr[j] = 0;
	}
	bits[0] = 0;
	for (int i = 1; i < PW; i++) bits[i] = bits[i / 2] + i % 2;
}

int popcount(unsigned long long x) {
	return bits[x % PW] + bits[x / PW % PW] + bits[x / PW / PW];
}

int Rand(int x) {
	return (rand() + (rand() << 15)) % x;
}

ll slowsolve(string a) {
	int n = a.length();
	vector<int> s(n + 1);
	for (int i = 0; i < n; i++) s[i + 1] = s[i] + (a[i] == 'b');

	ll ans = 0;
	for (int i = 0; i < n; i++) {
		for (int d = 1; i + 2 * d <= n; d++) {
			ans += s[i + d] * 2 == s[i] + s[i + 2 * d];
		}
	}
	return ans;
}

const int SQ = 10;
struct my_bitset {
	int n;
	vector<ll> a;

	my_bitset(int _n) {
		this->n = _n;
		a.resize((n + 63) / 64);
	}

	void change(int id) {
		a[id / 64] ^= 1LL << (id % 64);
	}

	my_bitset operator& (const my_bitset &rhs) const {
		my_bitset res(n);
		for (int i = 0; i < (int)a.size(); i++) res.a[i] = a[i] & rhs.a[i];
		return move(res);
	}

	int count() const {
		int ans = 0;
		for (ll x : a) ans += popcount(x);
		return ans;
	}

	int get(int id) const {
		return (a[id / 64] >> (id % 64)) & 1;
	}
};

void BSGLemma(int n, int m, vector<pt> G, double alpha, vector<int> &ra, vector<int> &rb) {
	//cerr << "BSG" << endl;
	vector<my_bitset> e(n, my_bitset(m));
	vector<vector<int> > g(m);
	vector<int> deg(n);
	for (auto p : G) deg[p.x]++;
	for (auto p : G) g[p.y].push_back(p.x);
	for (auto p : G) e[p.x].change(p.y);
	vector<int> p(m);
	for (int i = 0; i < m; i++) p[i] = i;
	random_shuffle(p.begin(), p.end());
	vector<int> A0;

	for (int i = 0; i < n; i++) if (deg[i] >= alpha * m / 2) A0.push_back(i);
	for (int jj : p) {
		vector<int> A;
		for (int x : A0) if (e[x].get(jj)) A.push_back(x);
		vector<pt> bad;
		for (int i = 0; i < (int)A.size(); i++) {
			for (int j = i + 1; j < (int)A.size(); j++) {
				if ((e[A[i]] & e[A[j]]).count() <= alpha * alpha * alpha * m / 2048) {
					bad.push_back(pt(A[i], A[j]));
				}
			}
		}
		if (A.size() < alpha * n / 4) continue;
		if (bad.size() <= alpha * alpha * n * A.size() / 256) {
			vector<int> deg2(n);
			for (pt pp : bad) deg2[pp.x]++, deg2[pp.y]++;
			vector<int> aP, bP;
			for (int x : A) if (deg2[x] <= alpha * alpha * A.size() / 64) aP.push_back(x);
			vector<int> inA(n);
			for (int x : A) inA[x] = 1;
			vector<int> deg3(m);
			for (auto pp : G) if (inA[pp.x]) deg3[pp.y]++;
			for (int i = 0; i < m; i++) if (deg3[i] >= alpha * aP.size() / 4) bP.push_back(i);
			ra = aP;
			rb = bP;
			return;
		}
	}
	assert(0);
}

void BSGCorollary(vector<pt> A, vector<pt> B, vector<pt> C, double alpha, vector<vector<int> > &ra, vector<vector<int> > &rb, vector<pt> &r) {
	vector<pt > G;
	int n = A.size();
	int m = B.size();
	for (int i = 0; i < n; i++) {
		int k = 0;
		for (int j = 0; j < m; j++) {
			while (k < (int)C.size() && A[i] + B[j] > C[k]) k++;
			if (k < (int)C.size() && A[i] + B[j] == C[k]) G.push_back(pt(i, j));
		}
	}
	double need = A.size() * B.size() * alpha;
	while (G.size() > need) {
		ra.push_back(vector<int>());
		rb.push_back(vector<int>());
		BSGLemma(n, m, G, G.size() * 1.0 / (A.size() * B.size()), ra.back(), rb.back()); // TODO
		vector<char> ca(n), cb(m);
		for (int x : ra.back()) ca[x] = 1;
		for (int x : rb.back()) cb[x] = 1;
		vector<pt> nG;
		for (auto p : G) if (!ca[p.x] || !cb[p.y]) nG.push_back(p);
		assert(G.size() > nG.size());
		G = nG;
	}
	r = G;
	return;
}

ll FFTLemma(vector<pt> a, vector<pt> b, vector<pt> t, vector<pt> c) {
	ll ans = 0;
	//cerr << "FFT" << endl;
	//return -1;
	/*for (int i = 0; i < (int)a.size(); i++) {
		int k = 0;
		for (int j = 0; j < (int)b.size(); j++) {
			while (k < (int)c.size() && a[i] + b[j] > c[k]) k++;
			ans += k < (int)c.size() && a[i] + b[j] == c[k];
		}
	}
	return ans;*/
	int U = 0;
	for (auto p : a) U = max(U, max(p.x, p.y));
	for (auto p : b) U = max(U, max(p.x, p.y));
	for	(auto p : t) U = max(U, max(p.x, p.y));
	U++;
	vector<ll> aa, bb, cc;
	for (auto p : a) aa.push_back(1LL * U * p.x + p.y);
	for (auto p : b) bb.push_back(1LL * U * p.x + p.y);
	for (auto p : t) cc.push_back(1LL * U * p.x + p.y);
	//sort(cc.begin(), cc.end());

	int N = max((int)(t.size() * 1), 8); // TODO LOL
	int cnt = 0;
	vector<int> distinguished(t.size(), -1);
	vector<int> primes;
	vector<vector<ll> > product;

	while (cnt < (int)t.size()) {
		N *= 1.5;
		//N = min(N, (int)t.size() * 10);
		//N = max(N + 10, (int)(N * 1.01));
		vector<int> p;
		for (int i = 0; i < 5; i++) p.push_back(*lower_bound(pr.begin(), pr.end(), Rand(N - N / 2) + N / 2));
		int mx = -1, pp = -1;
		for (int i = 0; i < (int)p.size(); i++) {
			vector<int> curcnt(p[i]);
			for (ll x : cc) curcnt[x % p[i]]++;
			int cur = 0;
			for (int j = 0; j < (int)t.size(); j++) cur += distinguished[j] == -1 && curcnt[cc[j] % p[i]] == 1;
			if (mx < cur) {
				mx = cur;
				pp = p[i];
			}
		}
		if (mx < ((int)t.size() - cnt) * 0.01) {
			//N *= 2;
			//continue;
		}
		if (mx > 0) {
			cnt += mx;
			vector<int> curcnt(pp);
			for (ll x : cc) curcnt[x % pp]++;
			for (int j = 0; j < (int)t.size(); j++) if (distinguished[j] == -1 && curcnt[cc[j] % pp] == 1) distinguished[j] = primes.size();
			for (int i = 0; i < pp; i++) FFT::A[i] = FFT::B[i] = 0;
			for (int i = 0; i < (int)aa.size(); i++) FFT::A[aa[i] % pp]++;
			for (int i = 0; i < (int)bb.size(); i++) FFT::B[bb[i] % pp]++;
			FFT::mult(pp, pp);
			product.push_back(vector<ll>(pp));
			for (int i = 0; i < pp + pp - 1; i++) product.back()[i % pp] += FFT::C[i];
			primes.push_back(pp);
		}
	}
	for (int ii = 0; ii < (int)c.size(); ii++) {
		ll hsh = 1LL * c[ii].x * U + c[ii].y;
		int i = lower_bound(cc.begin(), cc.end(), hsh) - cc.begin();
		ans += product[distinguished[i]][cc[i] % primes[distinguished[i]]];
	}
	return ans;
}

ll solve3Sum(vector<pt > A, vector<pt > B, vector<pt > C, int h) {
	int n = !A.empty() && !B.empty() && !C.empty() ? max({ A.back().x, A.back().y, B.back().x, B.back().y, C.back().y }) : 0;
	if (n <= SQ * 1 || h > 0) {
		int ans = 0;
		for (int i = 0; i < (int)A.size(); i++) {
			int k = 0;
			int k2 = 0;
			for (int j = 0; j < (int)B.size(); j++) {
				while (k < (int)C.size() && A[i] + B[j] > C[k]) {
					k++;
				}
				while (k2 < (int)C.size() && !(A[i] + B[j] < C[k2])) k2++;
				ans += k2 - k;
			}
		}
		return ans;
	}
	int l = pow(n, 0.0707);
	l = max(l, 2);
	if (l % 2) l++;
	double alpha = pow(n, -0.1313) / 3;

	auto decompose = [&](vector<pt > a) {
		vector<vector<pt > > res;
		for (int mask = 0; mask < 4; mask++) {
			vector<pt > cur;
			for (int i = 0; i < (int)a.size(); i++) {
				int cmask = (int)(a[i].x % l * 2 >= l) + ((int)(a[i].y % l * 2 >= l) * 2);
				if (cmask == mask) {
					cur.push_back(a[i]);
				}
			}
			res.push_back(cur);
		}
		return res;
	};
	vector<vector<pt > > ca = decompose(A);
	vector<vector<pt > > cb = decompose(B);
	int ans = 0;
	for (int ma = 0; ma < 4; ma++) {
		for (int mb = 0; mb < 4; mb++) {
			vector<pt > a = ca[ma];
			vector<pt > b = cb[mb];
			vector<pt > c = C;
			if (ma & 1) {
				for (auto &x : a) x.x += l / 2;
				for (auto &x : c) x.x += l / 2;
			}
			if (ma & 2) {
				for (auto &x : a) x.y += l / 2;
				for (auto &x : c) x.y += l / 2;
			}
			if (mb & 1) {
				for (auto &x : b) x.x += l / 2;
				for (auto &x : c) x.x += l / 2;
			}
			if (mb & 2) {
				for (auto &x : b) x.y += l / 2;
				for (auto &x : c) x.y += l / 2;
			}

			if (a.size() * b.size() <= 1 * SQ * SQ) {
				int res = 0;
				for (auto aa : a) {
					for (auto bb : b) {
						res += upper_bound(c.begin(), c.end(), aa + bb) -
							   lower_bound(c.begin(), c.end(), aa + bb);
					}
				}
				ans += res;
				continue;
			}
			auto getCells = [&](vector<pt > aa) {
				vector<pt > res;
				for (int i = 0; i < (int)aa.size(); i++) {
					res.push_back(pt(aa[i].x / l, aa[i].y / l));
				}
				sort(res.begin(), res.end());
				res.resize(unique(res.begin(), res.end()) - res.begin());
				return res;
			};
			vector<pt > cellA, cellB, cellC;
			cellA = getCells(a);
			cellB = getCells(b);
			cellC = getCells(c);
			vector<vector<int> > ra, rb;
			vector<pt > r;
			BSGCorollary(cellA, cellB, cellC, alpha, ra, rb, r);
			vector<vector<pt> > bA(cellA.size()), bB(cellB.size()), bC(cellC.size());
			auto fill = [&](const vector<pt> &aa, const vector<pt> &bb, vector<vector<pt> > &cc) {
				int j = 0;
				for (auto x : aa) {
					if (pt(x.x / l, x.y / l) > bb[j]) j++;
					cc[j].push_back(x);
				}
			};
			fill(a, cellA, bA);
			fill(b, cellB, bB);
			fill(c, cellC, bC);
			int k = 0;
			for (int i = 0; i < (int)r.size(); i++) {
				if (i == 0 || r[i - 1].x != r[i].x) k = 0;
				while (cellA[r[i].x] + cellB[r[i].y] > cellC[k]) k++;
				assert(k < (int)cellC.size() && cellA[r[i].x] + cellB[r[i].y] == cellC[k]);
				vector<pt> na = bA[r[i].x], nb = bB[r[i].y], nc = bC[k];
				for (pt &p : na) p.x %= l, p.y %= l;
				for (pt &p : nb) p.x %= l, p.y %= l;
				for (pt &p : nc) p.x %= l, p.y %= l;
				ans += solve3Sum(na, nb, nc, h + 1);
			}
			for (int i = 0; i < (int)ra.size(); i++) {
				vector<pt> t;
				for (int iA : ra[i]) {
					for (int iB : rb[i]) {
						t.push_back(cellA[iA] + cellB[iB]);
					}
				}
				vector<pt> pA, pB;
				for (int iA : ra[i]) for (pt p : bA[iA]) pA.push_back(p);
				for (int iB : rb[i]) for (pt p : bB[iB]) pB.push_back(p);
				sort(t.begin(), t.end());
				t.resize(unique(t.begin(), t.end()) - t.begin());
				vector<pt> tt;
				for (pt p : t) {
					for (int ii = 0; ii < l; ii++) {
						for (int j = 0; j < l; j++) {
							tt.push_back(pt(p.x * l + ii, p.y * l + j));
						}
					}
				}
				sort(tt.begin(), tt.end());
				vector<pt> pC;
				for (pt p : t) {
					int id = lower_bound(cellC.begin(), cellC.end(), p) - cellC.begin();
					if (id != (int)cellC.size() && cellC[id] == p) {
						for (pt pp : bC[id]) {
							pC.push_back(pp);
						}
					}
				}
				ans += FFTLemma(pA, pB, tt, pC);
			}
		}
	}
	return ans;
}

ll fastsolve(string a) {
	precalc();
	int n = a.length();
	vector<int> s(n + 1);
	for (int i = 0; i < n; i++) s[i + 1] = s[i] + (a[i] == 'b');

	vector<pt > A;
	vector<pt > C;
	for (int i = 0; i <= n; i++) {
		A.push_back(pt(i - s[i], s[i]));
		C.push_back(pt(2 * (i - s[i]), 2 * s[i]));
	}
	
	ll ans = solve3Sum(A, A, C, 0);
	return (ans - (n + 1)) / 2;
	/*ll ans = 0;
	for (int i = 0; i < n; i++) {
		int k = 0;
		for (int j = i + 1; j <= n; j++) {
			pair<int, int> o = make_pair(A[i].first + A[j].first, A[i].second + A[j].second);
			while (k <= n && C[k] < o) k++;
			ans += k <= n && C[k] == o;
		}
	}
	return ans;*/
}

void stress() {
	cerr.precision(3);
	cerr << fixed;
	for (int it = 19;; it++) {
		srand(it);
		
		//int n = (int)(4e4);
		int n = 1000;
		string s = "";
		for (int i = 0; i < n; i++) s += (char)('a' + rand() % 2);
		//for (int i = 0; i < n; i++) s += (char)('a');

		double st, fn;

		st = clock() / (double)CLOCKS_PER_SEC;
		ll ans2 = slowsolve(s);
		fn = clock() / (double)CLOCKS_PER_SEC;
		double t1 = fn - st;

		//cerr << it << " " << ans2 << " " << fn - st << endl;

		st = clock() / (double)CLOCKS_PER_SEC;
		ll ans1 = fastsolve(s);
		fn = clock() / (double)CLOCKS_PER_SEC;
		double t2 = fn - st;

		//cerr << ans1 << " " << fn - st << endl;
		if (ans1 != ans2) {
			cerr << ans1 << " instead of " << ans2 << endl;
			cout << s << endl;
			slowsolve(s);
			fastsolve(s);
			exit(0);
		}
		
		cerr << "it: " << it << ", ans: " << ans1 << ", time_slow: " << t1 << ", time_fast: " << t2 << endl;
	}
}

void testSpeedEasy() {
    cout.precision(3);
    cout << fixed;
    const int ITERS = 30;
    for (int n = 300; n <= 5000; n += 300) {
        double sum = 0;
        for (int it = 0; it < ITERS; it++) {
            cerr << "n: " << n << ", iteration " << it << " of " << ITERS << endl;
            string s(n, 'a');
            double st = clock() / (double)CLOCKS_PER_SEC;
            fastsolve(s);
            double fn = clock() / (double)CLOCKS_PER_SEC;
            sum += fn - st;
        } 
        cout << n << " " << sum / ITERS << endl;
    }
}

void testSpeedRandom() {
    cout.precision(3);
    cout << fixed;
    const int ITERS = 20;
    for (int n = 200; n <= 3000; n += 200) {
        double sum = 0;
        for (int it = 0; it < ITERS; it++) {
            cerr << "n: " << n << ", iteration " << it << " of " << ITERS << endl;
            string s;
            for (int i = 0; i < n; i++) s += (char)('a' + rand() % 2);
            double st = clock() / (double)CLOCKS_PER_SEC;
            fastsolve(s);
            double fn = clock() / (double)CLOCKS_PER_SEC;
            sum += fn - st;
        } 
        cout << n << " " << sum / ITERS << endl;
    }
}

int main() {
	freopen("input.txt", "r", stdin);
    testSpeedRandom();
    //testSpeedEasy();
	//stress();

	return 0;
}
