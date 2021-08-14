#include <bits/stdc++.h>
using namespace std;

template <typename T>
T modinv(T a, T MOD){
    T res = 1, p = a;
    if (p >= MOD){
        p %= MOD;
    }
    while (p > 1){
        res = (long long)res * (MOD - MOD / p) % MOD;
        p = MOD % p;
    }
    return res;
}

template <typename T>
T garner(vector<T> a, vector<T> m, T MOD){
    int n = a.size();
    vector<T> c(n, 0);
    T b = 1, res = 0, t, tmp;
    for (int i=0;i<n;i++){
        t = a[i];
        while (t >= m[i]){
            t -= m[i];
        }
        for (int j=0;j<i;j++){
            t -= c[j];
            while (t >= m[i]){
                t -= m[i];
            }
            while (t < 0){
                t += m[i];
            }
            t = (long long)t * modinv(m[j], m[i]) % m[i];
        }
        if (MOD){
            tmp = (long long)b * t % MOD;
            res = (res + tmp >= MOD ? res + tmp - MOD : res + tmp);
            b = (long long)b * m[i] % MOD;
        } else{
            res = res + b * t;
            b = b * m[i];
        }
        c[i] = t;
    }
    return res;
}

template <int m>
class modint{
    using mint = modint<m>;
    private:
    static constexpr int mod() {return m;}
    static constexpr unsigned int umod() {return m;}
    unsigned int x = 0;
    public:
    modint() {}
    template <class T>
    modint(T a){
        long long z = (long long)(a % mod());
        z = (z >= 0 ? z : z + umod());
        x = (unsigned int)z;
    }
    modint(bool a){
        x = (uint32_t)a;
    }
    unsigned int val(){
        return x;
    }
    mint inv() const{
        mint res = 1;
        int n = x;
        int mm = mod();
        while (n > 1){
            res *= (mm - mm / n);
            n = mm % n;
        }
        return res;
    }
    friend ostream& operator<<(ostream& os, const mint& rhs){
        os << rhs.x;
        return os;
    }
    friend istream& operator>>(istream& ist, mint& rhs){
        unsigned int s;
        ist >> s;
        rhs = s;
        return (ist);
    }
    bool operator==(const mint& rhs){
        return x == rhs.x;
    }
    template <typename T>
    bool operator==(const T rhs){
        mint res = mint(rhs);
        return x == res.x;
    }
    bool operator!=(const mint& rhs){
        return x != rhs.x;
    }
    template <typename T>
    bool operator!=(const T rhs){
        mint res = mint(rhs);
        return x != res.x;
    }
    mint operator+(){
        mint res = mint(*this);
        return res;
    }
    mint operator-(){
        mint res = mint(*this);
        if (res.x){
            res.x = umod() - res.x;
        }
        return res;
    }
    mint& operator++(){
        x++;
        if (x == umod()) x = 0;
        return *this;
    }
    mint& operator--(){
        if (x == 0) x = umod();
        x--;
        return *this;
    }
    mint operator++(int){
        mint res = *this;
        ++*this;
        return res;
    }
    mint operator--(int){
        mint res = *this;
        --*this;
        return res;
    }
    mint& operator+=(const mint& rhs){
        x = (x + rhs.x < umod() ? x + rhs.x : x + rhs.x - umod());
        return *this;
    }
    mint& operator-=(const mint& rhs){
        x = (x < rhs.x ? x + umod() - rhs.x : x - rhs.x);
        return *this;
    }
    mint& operator*=(const mint& rhs){
        x = (unsigned int)((long long)x * rhs.x % umod());
        return *this;
    }
    mint& operator/=(const mint& rhs){
        *this *= rhs.inv();
        return *this;
    }
    mint operator+(const mint& rhs){
        mint res = mint(*this);
        res += rhs;
        return res;
    }
    mint operator-(const mint& rhs){
        mint res = mint(*this);
        res -= rhs;
        return res;
    }
    mint operator*(const mint& rhs){
        mint res = mint(*this);
        res *= rhs;
        return res;
    }
    mint operator/(const mint& rhs){
        mint res = mint(*this);
        res /= rhs;
        return res;
    }
};

template <int m, typename T>
modint<m> modpow(modint<m> x, T n){
    modint<m> y = 1;
    modint<m> b = x;
    while (n){
        if (n & 1){
            y *= b;
        }
        b *= b;
        n >>= 1;
    }
    return y;
}

int getPrimitiveRoot(int MOD, int depth){
    int i = 0, t;
    bool flg;
    while (i < MOD){
        t = i;
        flg = true;
        for (int i=0;i<depth;i++){
            if (t == 1){
                flg = false;
                break;
            }
            t = (long long)t * t % MOD;
        }
        if (!flg || t != 1){
            i++;
            continue;
        }
        return i;
    }
    return -1;
}

template <int m>
vector<modint<m>> ntt(vector<modint<m>> f, bool isInverse=false){
    int n = f.size(), b = 1, d = 0;
    while (b < n){
        b <<= 1;
        d++;
    }
    vector<modint<m>> res(b);
    int pos = 0, dd = 0, tmp;
    for (int i=0;i<b;i++){
        pos ^= (b-(1<<(d-dd)));
        if (pos < n){
            res[i] = f[pos];
        }
        tmp = (i & (i + 1)) ^ (i + 1);
        dd = 0;
        while (tmp){
            dd++;
            tmp >>= 1;
        }
    }
    int bm = 0;
    tmp = m - 1;
    while ((tmp + 1) & 1){
        bm++;
        tmp >>= 1;
    }
    modint<m> root = getPrimitiveRoot(m, bm);
    if (!isInverse){
        root = root.inv();
    }
    vector<modint<m>> rootList(bm+1, 0);
    rootList[bm] = root;
    for (int i=bm;i>0;i--){
        rootList[i-1] = rootList[i];
        rootList[i-1] *= rootList[i];
    }
    modint<m> grow = 1, seed = 0, u, v;
    int offset = 1;
    for (int i=0;i<d;i++){
        grow = 1;
        seed = rootList[i+1];
        offset = 1 << i;
        for (int k=0;k<offset;k++){
            for (int j=k;j<b;j+=1<<(i+1)){
                u = res[j];
                v = res[j+offset] * grow;
                res[j] = u + v;
                res[j+offset] = u - v;
            }
            grow *= seed;
        }
    }
    if (isInverse){
        u = b;
        u = u.inv();
        for (int i=0;i<b;i++){
            res[i] *= u;
        }
    }
    return res;
}

template <int m>
vector<modint<m>> convolveOne(vector<modint<m>> f, vector<modint<m>> g){
    int nf = (int)f.size(), ng = (int)g.size();
    int b = 1;
    while (b < nf + ng - 1){
        b <<= 1;
    }
    vector<modint<m>> x(b);
    vector<modint<m>> y(b);
    for (int i=0;i<nf;i++){
        x[i] = f[i];
    }
    for (int i=0;i<ng;i++){
        y[i] = g[i];
    }
    x = ntt<m>(x);
    y = ntt<m>(y);
    for (int i=0;i<b;i++){
        x[i] *= y[i];
    }
    x = ntt<m>(x, true);
    x.resize(nf+ng-1);
    return x;
}

template <int MOD, int m0 = 645922817, int m1 = 924844033, int m2 = 998244353>
vector<int> convolve(vector<int> f, vector<int> g){
    int n = (int)f.size(), m = (int)g.size();
    if (MOD == m0 || MOD == m1 || MOD == m2){
        vector<modint<MOD>> x(n);
        vector<modint<MOD>> y(m);
        for (int i=0;i<n;i++){
            x[i] = f[i];
        }
        for (int i=0;i<m;i++){
            y[i] = g[i];
        }
        x = convolveOne<MOD>(x, y);
        vector<int> res(n+m-1);
        for (int i=0;i<n+m-1;i++){
            res[i] = x[i].val();
        }
        return res;
    }
    vector<modint<m0>> x0(n), y0(m);
    vector<modint<m1>> x1(n), y1(m);
    vector<modint<m2>> x2(n), y2(m);
    for (int i=0;i<n;i++){
        x0[i] = f[i];
        x1[i] = f[i];
        x2[i] = f[i];
    }
    for (int i=0;i<m;i++){
        y0[i] = g[i];
        y1[i] = g[i];
        y2[i] = g[i];
    }
    x0 = convolveOne<m0>(x0, y0);
    x1 = convolveOne<m1>(x1, y1);
    x2 = convolveOne<m2>(x2, y2);
    vector<int> al(3);
    vector<int> ml = {m0, m1, m2};
    vector<int> res(n+m-1);
    for (int i=0;i<n+m-1;i++){
        al[0] = x0[i].val();
        al[1] = x1[i].val();
        al[2] = x2[i].val();
        res[i] = garner(al, ml, MOD);
    }
    return res;
}

template <typename T>
void printVector(vector<T> a){
    int n = (int)a.size();
    if (n <= 0){
        return;
    }
    cout << a[0];
    for (int i=1;i<n;i++){
        cout << " " << a[i];
    }
    cout << endl;
}

constexpr const int mod = 998244353;
using mint = modint<mod>;
