#include <bits/stdc++.h>
using namespace std;

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

template <typename T, typename S, typename U>
T modpow(T n, S m, U p){
    T y = 1;
    T b = n;
    while (m){
        if (m & 1){
            y = (long long)y * b % p;
        }
        b = (long long)b * b % p;
        m >>= 1;
    }
    return y;
}

int modsqrt(int a, int m){
    if (m == 2){
        return a;
    }
    if (a == 0){
        return 0;
    }
    if (modpow(a, m>>1, m) == m - 1){
        return -1;
    }
    int z = 1;
    while (modpow(z, m>>1, m) == 1){
        z++;
    }
    int q = m - 1, s = 0;
    while (~q & 1){
        q >>= 1;
        s++;
    }
    z = modpow(z, q, m);
    int r = modpow(a, (q+1)>>1, m);
    int t = modpow(a, q, m);
    int b, i, c;
    while (t - 1){
        c = t;
        i = 0;
        while (c - 1){
            i++;
            c = (long long)c * c % m;
        }
        b = modpow(z, 1<<(s-i-1), m);
        z = (long long)b * b % m;
        s = i;
        t = (long long)t * z % m;
        r = (long long)r * b % m;
    }
    if (r > m - r){
        r = m - r;
    }
    return r;
}

template <const int MOD>
class FormalPowerSeries{
    using T = modint<MOD>;
    using FPS = FormalPowerSeries<MOD>;
    protected:
    T primitiveRoot = 0;
    T inverseRoot = 0;
    vector<T> primitiveBaseList;
    vector<T> inverseBaseList;
    int baseSize = 0;
    template <typename Tp>
    static T modpow(T a, Tp b){
        T y = 1;
        T bas = a;
        while (b){
            if (b & 1){
                y *= bas;
            }
            bas *= bas;
            b >>= 1;
        }
        return y;
    }
    void initialize(){
        if (MOD == 998244353){
            inverseRoot = (T)31;
            primitiveRoot = inverseRoot.inv();
            baseSize = 23;
            primitiveBaseList.resize(baseSize+1);
            inverseBaseList.resize(baseSize+1);
        }
        if (primitiveRoot == 0){
            int q = MOD-1;
            while (!(q & 1)){
                q >>= 1;
                baseSize++;
            }
            primitiveBaseList.resize(baseSize+1);
            inverseBaseList.resize(baseSize+1);
            while (true){
                bool flg = true;
                T tmp = inverseRoot;
                for (int i=0;i<baseSize;i++){
                    if (tmp == 1){
                        flg = false;
                        break;
                    }
                    tmp *= tmp;
                }
                if (tmp != 1){
                    flg = false;
                }
                if (flg){
                    break;
                }
                inverseRoot++;
            }
            primitiveRoot = modpow(inverseRoot, MOD-2);
        }
        primitiveBaseList[baseSize] = primitiveRoot;
        inverseBaseList[baseSize] = inverseRoot;
        for (int i=baseSize;i>0;i--){
            primitiveBaseList[i-1] = primitiveBaseList[i] * primitiveBaseList[i];
            inverseBaseList[i-1] = inverseBaseList[i] * inverseBaseList[i];
        }
    }
    vector<T> getNumberTheoremTransform(vector<T> v){
        int n = v.size();
        int m = 1;
        int depth = 0;
        while (m < n){
            depth++;
            m <<= 1;
        }
        int pos = 0, tmp = 0, d;
        vector<T> res(m, 0);
        for (int i=0;i<m-1;i++){
            if (pos < n){
                res[i] = v[pos];
            }
            tmp = (i+1)&-(i+1);
            d = 0;
            while (tmp){
                d++;
                tmp >>= 1;
            }
            pos ^= n - (1<<(depth - d));
        }
        if (pos < n){
            res[m-1] = v[pos];
        }
        T grow = 1, seed = 0;
        int offset = 0;
        T x, y;
        for (int i=0;i<depth;i++){
            grow = 1;
            seed = primitiveBaseList[i+1];
            offset = 1 << i;
            for (int k=0;k<offset;k++){
                for (int j=k;j<n;j+=1<<(i+1)){
                    x = res[j];
                    y = res[j+offset] * grow;
                    res[j] = x + y;
                    res[j+offset] = x - y;
                }
                grow *= seed;
            }
        }
        return res;
    }
    vector<T> getInverseNumberTheoremTransform(vector<T> v){
        int n = v.size();
        int m = 1;
        int depth = 0;
        while (m < n){
            depth++;
            m <<= 1;
        }
        int pos = 0, tmp = 0, d;
        vector<T> res(m, 0);
        for (int i=0;i<m-1;i++){
            if (pos < n){
                res[i] = v[pos];
            }
            tmp = (i+1)&-(i+1);
            d = 0;
            while (tmp){
                d++;
                tmp >>= 1;
            }
            pos ^= n - (1<<(depth - d));
        }
        if (pos < n){
            res[m-1] = v[pos];
        }
        T grow = 1, seed = 0;
        int offset = 0;
        T x, y;
        for (int i=0;i<depth;i++){
            grow = 1;
            seed = inverseBaseList[i+1];
            offset = 1 << i;
            for (int k=0;k<offset;k++){
                for (int j=k;j<n;j+=1<<(i+1)){
                    x = res[j];
                    y = res[j+offset] * grow;
                    res[j] = x + y;
                    res[j+offset] = x - y;
                }
                grow *= seed;
            }
        }
        T mInverse = modpow(m, MOD-2);
        for (int i=0;i<m;i++){
            res[i] *= mInverse;
        }
        return res;
    }
    vector<T> convolve(vector<T> f, vector<T> g){
        int n = f.size(), m = g.size();
        int b = 1;
        while (b < n + m){
            b <<= 1;
        }
        vector<T> x(b, 0);
        vector<T> y(b, 0);
        for (int i=0;i<n;i++){
            x[i] = f[i];
        }
        for (int i=0;i<m;i++){
            y[i] = g[i];
        }
        x = getNumberTheoremTransform(x);
        y = getNumberTheoremTransform(y);
        for (int i=0;i<b;i++){
            x[i] *= y[i];
        }
        x = getInverseNumberTheoremTransform(x);
        return x;
    }
    void regularize(){
        int t = degree;
        for (int i=degree;i>=0;i--){
            if (f[i] != 0){
                break;
            }
            t--;
        }
        if (t == degree){
            return;
        }
        f.resize(t+1);
        degree = t;
    }
    public:
    int degree = -1;
    vector<T> f;
    FormalPowerSeries(){
        initialize();
    }
    FormalPowerSeries(const int n, bool reg=true){
        degree = 0;
        f.resize(1);
        f[0] = (T)n;
        initialize();
        if (reg){
            regularize();
        }
    }
    FormalPowerSeries(const T n, bool reg=true){
        degree = 0;
        f.resize(1);
        f[0] = n;
        initialize();
        if (reg){
            regularize();
        }
    }
    FormalPowerSeries(const vector<T> g, bool reg=true){
        degree = g.size() - 1;
        f.resize(degree+1);
        for (int i=0;i<=degree;i++){
            f[i] = g[i];
        }
        initialize();
        if (reg){
            regularize();
        }
    }
    FormalPowerSeries(const FPS& fps, bool reg=true){
        degree = fps.degree;
        f.resize(degree+1);
        for (int i=0;i<=degree;i++){
            f[i] = fps.f[i];
        }
        initialize();
        if (reg){
            regularize();
        }
    }
    void resize(const int size){
        degree = size - 1;
        f.resize(size);
    }
    void addManually(const int i, const T a){
        if (degree < i){
            degree = i;
            f.resize(degree+1);
        }
        f[i] = a;
    }
    T operator[](const int i) const{
        return f[i];
    }
    T substitute(const T x) const{
        T res = 0;
        T t = 1;
        for (int i=0;i<=degree;i++){
            res += f[i] * t;
            t *= x;
        }
        return t;
    }
    T substitute(const int x) const{
        return substitute(T(x));
    }
    friend ostream& operator<<(ostream& os, FPS& fps){
        if (fps.degree == -1){
            os << 0;
            return os;
        }
        for (int i=0;i<=fps.degree;i++){
            if (i){
                os << " ";
            }
            os << fps.f[i];
        }
        return os;
    }
    friend istream& operator>>(istream& ist, FPS& fps){
        for (int i=0;i<=fps.degree;i++){
            T s;
            ist >> s;
            fps.f[i] = s;
        }
        return ist;
    }
    bool operator==(const FPS& fps) const{
        int lt = degree, rt = fps.degree;
        for (int i=degree;i>=0;i--){
            if (f[i] != (T)0){
                break;
            }
            lt--;
        }
        for (int i=fps.degree;i>=0;i--){
            if (fps.f[i] != (T)0){
                break;
            }
            rt--;
        }
        if (lt != rt){
            return false;
        }
        for (int i=0;i<=lt;i++){
            if (f[i] != (T)fps.f[i]){
                return false;
            }
        }
        return true;
    }
    bool operator!=(const FPS& fps) const{
        return !(*this == fps);
    }
    FPS rev(){
        FPS res = FPS(*this);
        for (int i=0;2*i<=degree;i++){
            T u = res.f[i];
            res.f[i] = res.f[degree-i];
            res.f[degree-i] = u;
        }
        return res;
    }
    FPS operator+() const{
        FPS res = FPS(*this);
        return res;
    }
    FPS operator-() const{
        FPS res = FPS(*this);
        for (int i=0;i<=degree;i++){
            res.f[i] = -res.f[i];
        }
        return res;
    }
    FPS& operator++(){
        if (degree == -1){
            degree = 0;
            f.resize(1);
            f[0] = 0;
        }
        f[0]++;
        regularize();
        return *this;
    }
    FPS& operator--(){
        if (degree == -1){
            degree = 0;
            f.resize(1);
            f[0] = 0;
        }
        f[0]--;
        regularize();
        return *this;
    }
    FPS operator++(int){
        FPS res = FPS(*this);
        ++*this;
        return res;
    }
    FPS operator--(int){
        FPS res = FPS(*this);
        --*this;
        return res;
    }
    FPS& operator+=(const FPS& rhs){
        if (degree < rhs.degree){
            f.resize(rhs.degree+1);
            degree = rhs.degree;
        }
        for (int i=0;i<=rhs.degree;i++){
            f[i] += rhs.f[i];
        }
        regularize();
        return *this;
    }
    FPS& operator+=(const int rhs){
        if (degree == -1){
            f.resize(1);
            degree = 0;
            f[0] = 0;
        }
        f[0] += (T)rhs;
        return *this;
    }
    FPS& operator+=(const T rhs){
        if (degree == -1){
            f.resize(1);
            degree = 0;
            f[0] = 0;
        }
        f[0] += rhs;
        return *this;
    }
    FPS& operator-=(const FPS& rhs){
        if (degree < rhs.degree){
            f.resize(rhs.degree+1);
            degree = rhs.degree;
        }
        for (int i=0;i<=rhs.degree;i++){
            f[i] -= rhs.f[i];
        }
        regularize();
        return *this;
    }
    FPS& operator-=(const int rhs){
        if (degree == -1){
            f.resize(1);
            degree = 0;
            f[0] = 0;
        }
        f[0] -= (T)rhs;
        return *this;
    }
    FPS& operator-=(const T rhs){
        if (degree == -1){
            f.resize(1);
            degree = 0;
            f[0] = 0;
        }
        f[0] -= rhs;
        return *this;
    }
    FPS& operator*=(const FPS& rhs){
        degree += rhs.degree;
        if (degree == -2){
            degree = -1;
        }
        f = convolve(f, rhs.f);
        regularize();
        return *this;
    }
    FPS& operator*=(const int rhs){
        for (int i=0;i<=degree;i++){
            f[i] *= (T)rhs;
        }
        return *this;
    }
    FPS& operator*=(const T rhs){
        for (int i=0;i<=degree;i++){
            f[i] *= rhs;
        }
        return *this;
    }
    template <typename Tp>
    FPS operator+(const Tp rhs){
        FPS res = FPS(*this);
        res += rhs;
        return res;
    }
    template <typename Tp>
    FPS operator-(const Tp rhs){
        FPS res = FPS(*this);
        res -= rhs;
        return res;
    }
    template <typename Tp>
    FPS operator*(const Tp rhs){
        FPS res = FPS(*this);
        res *= rhs;
        return res;
    }
    FPS& shiftMultiply(int n){
        degree += n;
        f.resize(degree+1);
        for (int i=degree;i>=n;i--){
            f[i] = f[i-n];
        }
        for (int i=0;i<n;i++){
            f[i] = 0;
        }
        regularize();
        return *this;
    }
    FPS shiftMultiplied(int n){
        FPS res = FPS(*this);
        res.shiftMultiply(n);
        return res;
    }
    FPS& shiftDivide(int n){
        if (degree < n){
            degree = -1;
            f.resize(0);
            return *this;
        }
        for (int i=n;i<=degree;i++){
            f[i-n] = f[i];
        }
        degree -= n;
        f.resize(degree+1);
        regularize();
        return *this;
    }
    FPS shiftDivided(int n){
        FPS res = FPS(*this);
        res.shiftDivide(n);
        return res;
    }
    FPS& differentiate(){
        if (degree == -1){
            return *this;
        }
        for (int i=0;i<degree;i++){
            f[i] = (T)(i+1)*f[i+1];
        }
        f[degree] = 0;
        regularize();
        return *this;
    }
    FPS differentiated(){
        FPS res = FPS(*this);
        res.differentiate();
        return res;
    }
    FPS& integrate(){
        degree++;
        f.push_back(0);
        vector<T> invl(degree+1, 1);
        for (int i=2;i<=degree;i++){
            invl[i] = (T)(MOD / i) * (-invl[MOD%i]);
        }
        for (int i=degree;i>0;i--){
            f[i] = invl[i] * f[i-1];
        }
        f[0] = 0;
        regularize();
        return *this;
    }
    FPS integrated(){
        FPS res = FPS(*this);
        res.integrate();
        return res;
    }
    FPS& inverse(){
        T a = 0;
        if (degree == -1){
            degree = -2;
            return *this;
        }
        a = modpow(f[0], MOD-2);
        FPS res(a), twoFactor(2), tmp(f[0]), h;
        int b = 1;
        while (b < degree+1){
            for (int i=0;i<b;i++){
                if (i+b >= degree+1){
                    break;
                }
                tmp.f.push_back(f[i+b]);
                tmp.degree++;
            }
            b <<= 1;
            h = twoFactor - tmp * res;
            h.resize(b);
            res *= h;
        }
        res.resize(degree+1);
        *this = res;
        return *this;
    }
    FPS inversed(){
        FPS res = FPS(*this);
        res.inverse();
        return res;
    }
    FPS& operator/=(const FPS& rhs){
        if (degree < rhs.degree){
            *this = FPS(0);
            return *this;
        }
        FPS tmp = rhs;
        FPS fr = rev();
        FPS gr = tmp.rev();
        fr.resize(degree - rhs.degree + 1);
        gr.resize(degree - rhs.degree + 1);
        fr *= gr.inverse();
        fr.resize(degree - rhs.degree + 1);
        *this = fr.rev();
        regularize();
        return *this;
    }
    FPS& operator/=(const T rhs){
        for (int i=0;i<=degree;i++){
            f[i] /= rhs;
        }
        return *this;
    }
    template<typename Tp>
    FPS operator/(const Tp rhs){
        FPS res = FPS(*this);
        res /= rhs;
        return res;
    }
    FPS& operator%=(const FPS& rhs){
        *this -= (*this / rhs) * rhs;
        return *this;
    }
    FPS operator%(const FPS rhs){
        FPS res = FPS(*this);
        res %= rhs;
        return res;
    }
    FPS& convertTaylorShift(T c){
        if (degree <= 0){
            return *this;
        }
        vector<T> fact(degree+1, 1);
        vector<T> invf(degree+1, 1);
        for (int i=0;i<degree;i++){
            fact[i+1] = (fact[i] * (i+1));
        }
        invf[degree] = fact[degree].inv();
        for (int i=degree;i>0;i--){
            invf[i-1] = (invf[i] * i);
        }
        vector<T> x(degree+1, 0);
        vector<T> y(degree+1, 0);
        T t = 1;
        for (int i=0;i<=degree;i++){
            x[degree-i] = f[i];
            x[degree-i] *= fact[i];
            y[i] = t;
            y[i] *= invf[i];
            t *= c;
        }
        x = convolve(x, y);
        for (int i=0;i<=degree;i++){
            f[i] = invf[i];
            f[i] *= x[degree-i];
        }
        regularize();
        return *this;
    }
    FPS getTaylorShift(T c){
        FPS res = FPS(*this);
        res.convertTaylorShift(c);
        return res;
    }
    vector<T> getMultipointEvaluation(const vector<T>& p){
        int n = p.size();
        int b = n;
        if (b & (b - 1)){
            while (b & (b - 1)){
                b &= b - 1;
            }
            b <<= 1;
        }
        vector<FPS> fs(b<<1, 1);
        T t;
        for (int i=0;i<n;i++){
            t = p[i];
            fs[i+b] = FPS({-t, 1});
        }
        for (int i=b-1;i>=1;i--){
            fs[i] = fs[1|(i<<1)] * fs[i<<1];
        }
        fs[1] = *this % fs[1];
        for (int i=2;i<2*b;i++){
            fs[i] = fs[i>>1] % fs[i];
        }
        vector<T> res(n);
        for (int i=0;i<n;i++){
            res[i] = fs[i+b].f[0];
        }
        return res;
    }
    FPS& convertSqrt(){
        if (degree == -1){
            return *this;
        }
        int btm = 0;
        for (int i=0;i<=degree;i++){
            if (f[i] == 0){
                btm++;
            } else{
                break;
            }
        }
        if (btm & 1){
            degree = -2;
            return *this;
        }
        int z = degree - (btm >> 1);
        FPS res = FPS(*this, false).shiftDivided(btm);
        int b = 1;
        int t = modsqrt(res.f[0].val(), MOD);
        if (t == -1){
            degree = -2;
            return *this;
        }
        int bt = 1;
        while (bt <= z){
            bt <<= 1;
        }
        T x = t;
        FPS g(x), h(res.f[0]);
        T i2 = T(2).inv();
        while (b <= bt){
            for (int i=0;i<b;i++){
                if (i+b>res.degree){
                    break;
                }
                h.addManually(i+b, res.f[i+b]);
            }
            g += g.inversed() * h;
            b <<= 1;
            g *= i2;
            g.resize(b);
        }
        *this = g.shiftMultiply(btm>>1);
        return *this;
    }
    FPS getSqrt(){
        FPS res = FPS(*this);
        res.convertSqrt();
        return res;
    }
};

template <int m>
FormalPowerSeries<m> logarithm(FormalPowerSeries<m> fps){
    int n = fps.degree+1;
    FormalPowerSeries<m> res;
    res = fps.inversed() * fps.differentiated();
    res.resize(n);
    res.integrate();
    return res;
}

template <int m>
FormalPowerSeries<m> exponential(FormalPowerSeries<m> fps, int degreeRequire=-1){
    int resSize = 1;
    if (fps.degree < degreeRequire){
        fps.resize(degreeRequire+1);
    }
    int x = fps.degree + 1;
    if (x & (x - 1)){
        while (x & (x - 1)){
            x &= x - 1;
        }
        x <<= 1;
    }
    FormalPowerSeries<m> res(1), g(1), h(2), q, r, tmp;
    if (fps.degree == -1){
        tmp = 0;
    } else{
        tmp = fps.f[0];
    }
    while (2*resSize <= x){
        q = res * g;
        q.resize(resSize);
        g *= (h - q);
        g.resize(resSize);
        q = tmp.differentiated();
        r = q + g * (res.differentiated() - res * q);
        r.resize(2*resSize);
        r.integrate();
        for (int i=0;i<resSize;i++){
            if (i+resSize > fps.degree){
                break;
            }
            tmp.addManually(i+resSize, fps.f[i+resSize]);
        }
        q = tmp - r;
        q *= res;
        res += q;
        resSize <<= 1;
    }
    return res;
}

template <int m, typename T>
FormalPowerSeries<m> power(FormalPowerSeries<m> fps, T x){
    int t = -1;
    int d = fps.degree;
    modint<m> a = 0;
    for (int i=0;i<=d;i++){
        if (fps.f[i] != 0){
            t = i;
            a = fps.f[i];
            break;
        }
    }
    if (t == -1){
        if (x == 0){
            return FormalPowerSeries<m>(1);
        }
        return FormalPowerSeries<m>(0);
    }
    FormalPowerSeries<m> gps = fps;
    if ((long long) x * t > d){
        return FormalPowerSeries<m>(0);
    }
    gps.shiftDivide(t);
    for (int i=0;i<=gps.degree;i++){
        gps.f[i] /= a;
    }
    gps.resize(d+1);
    gps = logarithm(gps);
    for (int i=0;i<=d;i++){
        gps.f[i] *= x;
    }
    gps = exponential(gps, d);
    gps.shiftMultiply(x * t);
    gps.resize(d+1);
    a = modpow<m>(a, (long long)x);
    for (int i=0;i<=d;i++){
        gps.f[i] *= a;
    }
    return gps;
}

template <int m>
FormalPowerSeries<m> interpolate(vector<modint<m>> x, vector<modint<m>> y){
    int n = x.size();
    int b = 1;
    while (b < n){
        b <<= 1;
    }
    vector<FormalPowerSeries<m>> v(b << 1, 1);
    vector<FormalPowerSeries<m>> v2(b<<1, 0);
    for (int i=0;i<n;i++){
        v[i+b].resize(2);
        v[i+b].f[0] = -x[i];
        v[i+b].f[1] = 1;
    }
    for (int i=b-1;i>0;i--){
        v[i] = v[i<<1] * v[(i<<1)|1];
    }
    FormalPowerSeries<m> F = v[1];
    FormalPowerSeries<m> dF = F.differentiated();
    v2[1] = dF % v[1];
    for (int i=2;i<2*b;i++){
        v2[i] = v2[i>>1] % v[i];
    }
    for (int i=0;i<n;i++){
        v2[i+b].f[0] = y[i] / v2[i+b].f[0];
    }
    for (int i=b-1;i>0;i--){
        v2[i] = v2[i<<1]*v[(i<<1)|1] + v[i<<1]*v2[(i<<1)|1];
    }
    return v2[1];
}

template <int m>
vector<modint<m>> getBernoulliNumberTable(int n){
    vector<modint<m>> fact(n+2, 1);
    for (int i=1;i<=n+1;i++){
        fact[i] = (fact[i-1] * i);
    }
    vector<modint<m>> invf(n+1, 1);
    invf[n] = fact[n+1].inv();
    for (int i=n;i>0;i--){
        invf[i-1] = invf[i] * (i+1);
    }
    FormalPowerSeries<m> r(invf);
    r.inverse();
    vector<modint<m>> res(n+1);
    for (int i=0;i<=n;i++){
        res[i] = r.f[i];
        res[i] *= fact[i];
    }
    return res;
}

template <int m>
vector<modint<m>> getPartitionNumberTable(int n){
    vector<modint<m>> invn(n+1, 1);
    for (int i=2;i<=n;i++){
        invn[i] = -(invn[m%i]);
        invn[i] *= (m / i);
    }
    FormalPowerSeries<m> f;
    f.resize(n+1);
    for (int i=1;i<=n;i++){
        for (int j=1;i*j<=n;j++){
            f.f[i*j] += invn[j];
        }
    }
    f = exponential<m>(f);
    for (int i=0;i<=n;i++){
        invn[i] = f.f[i];
    }
    return invn;
}

template <int m>
vector<modint<m>> getStirlingNumberOfFirstTable(int n){
    FormalPowerSeries<m> g = 1;
    FormalPowerSeries<m> h({0, 1});
    int b = 1;
    int t = 0;
    while (t < n){
        if (b & n){
            g *= h.getTaylorShift(-modint<m>(t));
            t += b;
        }
        h *= h.getTaylorShift(-modint<m>(b));
        b <<= 1;
    }
    vector<modint<m>> res(n+1);
    for (int i=0;i<=n;i++){
        res[i] = g.f[i];
    }
    return res;
}

template <int m>
vector<modint<m>> getStirlingNumberOfSecondTable(int n){
    vector<modint<m>> a(n+1);
    vector<modint<m>> b(n+1);
    vector<modint<m>> invn(n+2, 1);
    for (int i=2;i<=n;i++){
        invn[i] = (-invn[m%i])*(m/i);
    }
    modint<m> t = 1, sign = 1;
    for (int i=0;i<=n;i++){
        a[i] = modpow<m>(i, n) * t;
        b[i] = sign * t;
        if (sign == 1){
            sign = m - 1;
        } else{
            sign = 1;
        }
        t *= invn[i+1];
    }
    FormalPowerSeries<m> f(a), g(b);
    f *= g;
    for (int i=0;i<=n;i++){
        a[i] = f.f[i];
    }
    return a;
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
