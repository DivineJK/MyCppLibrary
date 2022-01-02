#include <bits/stdc++.h>
using namespace std;

template <typename T>
void printVector(vector<T> a){
    int n = (int)a.size();
    if (n <= 0){
        cout << endl;
        return;
    }
    cout << a[0];
    for (int i=1;i<n;i++){
        cout << " " << a[i];
    }
    cout << endl;
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
    inline static bool isInitialized = false;
    inline static T primitiveRoot = 0;
    inline static T inverseRoot = 0;
    inline static vector<T> primitiveBaseList = {};
    inline static vector<T> inverseBaseList = {};
    inline static vector<T> cumulativeBase = {};
    inline static int baseSize = 0;
    bool isValid = true;
    int degree = -1;
    vector<T> polynomial;
    protected:
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
    static void initialize(){
        if (isInitialized){
            return;
        }
        isInitialized = true;
        if (MOD == 998244353){
            inverseRoot = (T)31;
            primitiveRoot = inverseRoot.inv();
            baseSize = 23;
            primitiveBaseList.resize(baseSize+1);
            inverseBaseList.resize(baseSize+1);
            cumulativeBase.resize(baseSize-1);
        }
        if (primitiveRoot == 0){
            int q = MOD-1;
            while (!(q & 1)){
                q >>= 1;
                baseSize++;
            }
            primitiveBaseList.resize(baseSize+1);
            inverseBaseList.resize(baseSize+1);
            cumulativeBase.resize(baseSize-1);
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
        T ie = 1;
        for (int i=0;i<baseSize-1;i++){
            cumulativeBase[i] = ie * primitiveBaseList[i+2];
            ie *= inverseBaseList[i+2];
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
        T grow = 1;
        int p = 1, t = 1 << (depth - 1), offset = 0, w, z;
        T x, y;
        for (int i=0;i<depth;i++){
            grow = 1;
            for (int j=0;j<p;j++){
                offset = j << (depth - i);
                for (int k=0;k<t;k++){
                    x = v[k+offset];
                    y = v[k+offset+t] * grow;
                    v[k+offset] = x + y;
                    v[k+offset+t] = x - y;
                }
                w = (j+1)^(j&(j+1)), z = -1;
                while (w){
                    z++;
                    w >>= 1;
                }
                grow *= cumulativeBase[z];
            }
            p <<= 1;
            t >>= 1;
        }
        return v;
    }
    vector<T> getInverseNumberTheoremTransform(vector<T> v){
        int n = v.size();
        int m = 1;
        int depth = 0;
        while (m < n){
            depth++;
            m <<= 1;
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
                    x = v[j];
                    y = v[j+offset] * grow;
                    v[j] = x + y;
                    v[j+offset] = x - y;
                }
                grow *= seed;
            }
        }
        T mInverse = modpow(m, MOD-2);
        for (int i=0;i<m;i++){
            v[i] *= mInverse;
        }
        return v;
    }
    vector<T> convolve(vector<T> f, vector<T> g){
        int n = f.size(), m = g.size();
        if (n == 0 || m == 0){
            return {};
        }
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
        x.resize(n+m-1);
        return x;
    }
    void regularize(){
        int t = degree;
        for (int i=degree;i>=0;i--){
            if (polynomial[i] != 0){
                break;
            }
            t--;
        }
        if (t == degree){
            return;
        }
        polynomial.resize(t+1);
        degree = t;
    }
    public:
    FormalPowerSeries(){
        isValid = true;
        initialize();
    }
    FormalPowerSeries(const int n){
        initialize();
        degree = 0;
        polynomial.resize(1);
        polynomial[0] = (T)n;
        regularize();
    }
    FormalPowerSeries(const T n){
        initialize();
        degree = 0;
        polynomial.resize(1);
        polynomial[0] = n;
        regularize();
    }
    FormalPowerSeries(const vector<int> g){
        initialize();
        degree = g.size() - 1;
        polynomial.resize(degree+1);
        for (int i=0;i<=degree;i++){
            polynomial[i] = (T)g[i];
        }
        regularize();
    }
    FormalPowerSeries(const vector<T> g){
        initialize();
        degree = g.size() - 1;
        polynomial.resize(degree+1);
        for (int i=0;i<=degree;i++){
            polynomial[i] = g[i];
        }
        regularize();
    }
    FormalPowerSeries(const FPS& fps){
        initialize();
        degree = fps.degree;
        polynomial.resize(degree+1);
        for (int i=0;i<=degree;i++){
            polynomial[i] = fps.polynomial[i];
        }
        regularize();
    }
    void setIsValid(const bool b){
        if (!b){
            setDegree(-1);
        }
        isValid = b;
    }
    bool getIsValid() const{
        return isValid;
    }
    void setDegree(const int aDegree){
        degree = aDegree;
        polynomial.resize(degree+1);
    }
    int getDegree() const{
        return degree;
    }
    void setCoefficient(const int i, const T a){
        if (degree < i){
            degree = i;
            polynomial.resize(degree+1);
        }
        polynomial[i] = a;
    }
    T getCoefficient(const int i) const{
        if (degree < i){
            return 0;
        }
        return polynomial[i];
    }
    T substitute(const T x) const{
        T res = 0;
        T t = 1;
        for (int i=0;i<=degree;i++){
            res += polynomial[i] * t;
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
            os << fps.polynomial[i];
        }
        return os;
    }
    friend istream& operator>>(istream& ist, FPS& fps){
        for (int i=0;i<=fps.degree;i++){
            T s;
            ist >> s;
            fps.polynomial[i] = s;
        }
        fps.regularize();
        return ist;
    }
    bool operator==(const FPS& fps) const{
        int lt = degree, rt = fps.degree;
        for (int i=degree;i>=0;i--){
            if (polynomial[i] != (T)0){
                break;
            }
            lt--;
        }
        for (int i=fps.degree;i>=0;i--){
            if (fps.polynomial[i] != (T)0){
                break;
            }
            rt--;
        }
        if (lt != rt){
            return false;
        }
        for (int i=0;i<=lt;i++){
            if (polynomial[i] != (T)fps.polynomial[i]){
                return false;
            }
        }
        return true;
    }
    bool operator!=(const FPS& fps) const{
        return !(*this == fps);
    }
    FPS& convertToReversedPolynomial(){
        for (int i=0;2*i<=degree;i++){
            T u = polynomial[i];
            polynomial[i] = polynomial[degree-i];
            polynomial[degree-i] = u;
        }
        return *this;
    }
    FPS getReversedPolynomial(){
        FPS res = FPS(*this);
        res.convertToReversedPolynomial();
        return res;
    }
    FPS operator+() const{
        FPS res = FPS(*this);
        return res;
    }
    FPS operator-() const{
        FPS res = FPS(*this);
        for (int i=0;i<=degree;i++){
            res.polynomial[i] = -res.polynomial[i];
        }
        return res;
    }
    FPS& operator++(){
        if (degree == -1){
            setDegree(0);
        }
        polynomial[0]++;
        regularize();
        return *this;
    }
    FPS& operator--(){
        if (degree == -1){
            setDegree(0);
        }
        polynomial[0]--;
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
            polynomial.resize(rhs.degree+1);
            degree = rhs.degree;
        }
        for (int i=0;i<=rhs.degree;i++){
            polynomial[i] += rhs.polynomial[i];
        }
        regularize();
        return *this;
    }
    FPS& operator+=(const int rhs){
        if (degree == -1){
            setDegree(0);
        }
        polynomial[0] += (T)rhs;
        regularize();
        return *this;
    }
    FPS& operator+=(const T rhs){
        if (degree == -1){
            setDegree(0);
        }
        polynomial[0] += rhs;
        regularize();
        return *this;
    }
    FPS& operator-=(const FPS& rhs){
        if (degree < rhs.degree){
            polynomial.resize(rhs.degree+1);
            degree = rhs.degree;
        }
        for (int i=0;i<=rhs.degree;i++){
            polynomial[i] -= rhs.polynomial[i];
        }
        regularize();
        return *this;
    }
    FPS& operator-=(const int rhs){
        if (degree == -1){
            setDegree(0);
        }
        polynomial[0] -= (T)rhs;
        regularize();
        return *this;
    }
    FPS& operator-=(const T rhs){
        if (degree == -1){
            setDegree(0);
        }
        polynomial[0] -= rhs;
        return *this;
    }
    FPS& operator*=(const FPS& rhs){
        if (degree == -1 || rhs.degree == -1){
            setDegree(-1);
            return *this;
        }
        degree += rhs.degree;
        polynomial = convolve(polynomial, rhs.polynomial);
        regularize();
        return *this;
    }
    FPS& operator*=(const int rhs){
        for (int i=0;i<=degree;i++){
            polynomial[i] *= (T)rhs;
        }
        return *this;
    }
    FPS& operator*=(const T rhs){
        for (int i=0;i<=degree;i++){
            polynomial[i] *= rhs;
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
    FPS& operator<<=(int n){
        degree += n;
        polynomial.resize(degree+1);
        for (int i=degree;i>=n;i--){
            polynomial[i] = polynomial[i-n];
        }
        for (int i=0;i<n;i++){
            polynomial[i] = 0;
        }
        regularize();
        return *this;
    }
    FPS operator<<(int n){
        FPS res = FPS(*this);
        res <<= n;
        return res;
    }
    FPS& operator>>=(int n){
        if (degree < n){
            setDegree(-1);
            return *this;
        }
        for (int i=n;i<=degree;i++){
            polynomial[i-n] = polynomial[i];
        }
        degree -= n;
        polynomial.resize(degree+1);
        regularize();
        return *this;
    }
    FPS operator>>(int n){
        FPS res = FPS(*this);
        res >>= n;
        return res;
    }
    FPS& differentiate(){
        if (degree == -1){
            return *this;
        }
        for (int i=0;i<degree;i++){
            polynomial[i] = (T)(i+1)*polynomial[i+1];
        }
        polynomial[degree] = 0;
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
        polynomial.push_back(0);
        vector<T> invl(degree+1, 1);
        for (int i=2;i<=degree;i++){
            invl[i] = (T)(MOD / i) * (-invl[MOD%i]);
        }
        for (int i=degree;i>0;i--){
            polynomial[i] = invl[i] * polynomial[i-1];
        }
        polynomial[0] = 0;
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
        a = modpow(polynomial[0], MOD-2);
        FPS res(a), twoFactor(2), tmp(polynomial[0]), h;
        int b = 1;
        while (b < degree+1){
            for (int i=0;i<b;i++){
                if (i+b >= degree+1){
                    break;
                }
                tmp.polynomial.push_back(polynomial[i+b]);
                tmp.degree++;
            }
            b <<= 1;
            h = twoFactor - tmp * res;
            h.setDegree(b-1);
            res *= h;
        }
        *this = res;
        regularize();
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
        FPS fr = getReversedPolynomial();
        FPS gr = tmp.getReversedPolynomial();
        fr.setDegree(degree - rhs.degree);
        gr.setDegree(degree - rhs.degree);
        fr *= gr.inverse();
        fr.setDegree(degree - rhs.degree);
        *this = fr.getReversedPolynomial();
        regularize();
        return *this;
    }
    FPS& operator/=(const T rhs){
        for (int i=0;i<=degree;i++){
            polynomial[i] /= rhs;
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
            x[degree-i] = polynomial[i];
            x[degree-i] *= fact[i];
            y[i] = t;
            y[i] *= invf[i];
            t *= c;
        }
        x = convolve(x, y);
        for (int i=0;i<=degree;i++){
            polynomial[i] = invf[i];
            polynomial[i] *= x[degree-i];
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
        vector<T> tmpt = {0, 1};
        vector<FPS> fs(b<<1, 1);
        T t;
        for (int i=0;i<n;i++){
            tmpt[0] = p[i];
            tmpt[0] *= -1;
            fs[i+b] = FPS(tmpt);
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
            if (fs[i+b].degree < 0){
                res[i] = 0;
                continue;
            }
            res[i] = fs[i+b].polynomial[0];
        }
        return res;
    }
    FPS& convertSqrt(){
        if (degree == -1){
            return *this;
        }
        int btm = 0;
        for (int i=0;i<=degree;i++){
            if (polynomial[i] == 0){
                btm++;
            } else{
                break;
            }
        }
        if (btm & 1){
            setIsValid(false);
            return *this;
        }
        int z = degree - (btm >> 1);
        FPS res = FPS(*this) >> btm;
        res.setDegree(degree - btm);
        int b = 1;
        int t = modsqrt(res.polynomial[0].val(), MOD);
        if (t == -1){
            setIsValid(false);
            return *this;
        }
        int bt = 1;
        while (bt <= z){
            bt <<= 1;
        }
        T x = t;
        FPS g(x), h(res.polynomial[0]);
        T i2 = T(2).inv();
        while (b <= bt){
            for (int i=0;i<b;i++){
                if (i+b>res.degree){
                    break;
                }
                h.setCoefficient(i+b, res.polynomial[i+b]);
            }
            g += g.inversed() * h;
            b <<= 1;
            g *= i2;
            g.setDegree(b-1);
        }
        *this = g << (btm>>1);
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
    int n = fps.getDegree();
    FormalPowerSeries<m> res;
    res = fps.inversed() * fps.differentiated();
    res.setDegree(n);
    res.integrate();
    return res;
}

template <int m>
FormalPowerSeries<m> exponential(FormalPowerSeries<m> fps, int degreeRequire=-1){
    int resSize = 1;
    if (fps.getDegree() < degreeRequire){
        fps.setDegree(degreeRequire);
    }
    int x = fps.getDegree() + 1;
    if (x & (x - 1)){
        while (x & (x - 1)){
            x &= x - 1;
        }
        x <<= 1;
    }
    FormalPowerSeries<m> res(1), g(1), h(2), q, r, tmp;
    if (fps.getDegree() == -1){
        tmp = 0;
    } else{
        tmp = fps.getCoefficient(0);
    }
    while (2*resSize <= x){
        q = res * g;
        q.setDegree(resSize-1);
        g *= (h - q);
        g.setDegree(resSize-1);
        q = tmp.differentiated();
        r = q + g * (res.differentiated() - res * q);
        r.setDegree(2*resSize-1);
        r.integrate();
        for (int i=0;i<resSize;i++){
            if (i+resSize > fps.getDegree()){
                break;
            }
            tmp.setCoefficient(i+resSize, fps.getCoefficient(i+resSize));
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
    int d = fps.getDegree();
    modint<m> a = 0;
    for (int i=0;i<=d;i++){
        if (fps.getCoefficient(i) != 0){
            t = i;
            a = fps.getCoefficient(i);
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
    gps >>= t;
    for (int i=0;i<=gps.getDegree();i++){
        gps.setCoefficient(i, gps.getCoefficient(i) / a);
    }
    gps.setDegree(d);
    gps = logarithm(gps);
    for (int i=0;i<=d;i++){
        gps.setCoefficient(i, gps.getCoefficient(i) * x);
    }
    gps = exponential(gps, d);
    gps <<= (x * t);
    gps.setDegree(d);
    a = modpow<m>(a, (long long)x);
    for (int i=0;i<=d;i++){
        gps.setCoefficient(i, gps.getCoefficient(i) * a);
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
        v[i+b] = FormalPowerSeries<m>(vector<modint<m>>({-x[i], 1}));
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
        v2[i+b].setCoefficient(0, y[i] / v2[i+b].getCoefficient(0));
    }
    for (int i=b-1;i>0;i--){
        v2[i] = v2[i<<1]*v[(i<<1)|1] + v[i<<1]*v2[(i<<1)|1];
    }
    return v2[1];
}

template <int m>
FormalPowerSeries<m> compose(FormalPowerSeries<m> lhs, FormalPowerSeries<m> rhs, int degreeRequire=-1){
    int n = (degreeRequire == -1 ? lhs.getDegree() : degreeRequire);
    if (n <= 0){
        return lhs;
    }
    int b = 1, c = 0;
    while (b <= n){
        b <<= 1;
        c++;
    }
    int l = 0, r = n + 2;
    int k = n / 2 + 1;
    while (r - l > 1){
        if ((long long)c * k * k <= n + 1){
            l = k;
        } else{
            r = k;
        }
        k = (l + r) / 2;
    }
    if (c * k * k <= n){
        k++;
    }
    FormalPowerSeries<m> p = rhs, q = rhs >> k;
    p.setDegree(k-1);
    FormalPowerSeries<m> tmpsq = p;
    FormalPowerSeries<m> bas = 1;
    if (tmpsq.getDegree() == -1){
        FormalPowerSeries<m> res = 0;
        for (int i=0;i<=(n+k-1)/k;i++){
            if (i <= lhs.getDegree()){
                tmpsq = bas * FormalPowerSeries<m>(lhs.getCoefficient(i));
                tmpsq <<= (i*k);
                res += tmpsq;
                if (res.getDegree() > n){
                    res.setDegree(n);
                }
                bas *= q;
                if (bas.getDegree() > n){
                    bas.setDegree(n);
                }
            }
        }
        return res;
    }
    vector<FormalPowerSeries<m>> v(b<<1, 0);
    for (int i=0;i<=n;i++){
        v[i+b] = FormalPowerSeries<m>(lhs.getCoefficient(i));
    }
    for (int i=b-1;i>0;i--){
        v[i] = v[i<<1] + v[(i<<1)|1] * tmpsq;
        if (v[i].getDegree() > n){
            v[i].setDegree(n);
        }
        if ((i & (i - 1)) == 0){
            tmpsq *= tmpsq;
            if (tmpsq.getDegree() > n){
                tmpsq.setDegree(n);
            }
        }
    }
    vector<modint<m>> invf((n+k-1)/k + 1, 1);
    for (int i=2;i<=(n+k-1)/k;i++){
        invf[i] = -(invf[m%i] * (m/i));
    }
    for (int i=0;i<(n+k-1)/k;i++){
        invf[i+1] *= invf[i];
    }
    FormalPowerSeries<m> t = v[1];
    FormalPowerSeries<m> ip = p.differentiated();
    int lz = 0;
    for (int i=0;i<=ip.getDegree();i++){
        if (ip.getCoefficient(i) != 0){
            break;
        }
        lz++;
    }
    ip >>= lz;
    ip.setDegree(n);
    ip.inverse();
    bas = 1;
    FormalPowerSeries<m> res = v[1];
    t.differentiate();
    t *= ip;
    t >>= lz;
    bas *= q;
    FormalPowerSeries<m> X = ip * p.differentiated();
    for (int i=1;i<=(n+lz+k-1)/k;i++){
        tmpsq = t * bas;
        tmpsq <<= (i*k);
        if (tmpsq.getDegree() > n){
            tmpsq.setDegree(n);
        }
        res += tmpsq * invf[i];
        bas *= q;
        if (bas.getDegree() > n){
            bas.setDegree(n);
        }
        if (res.getDegree() > n){
            res.setDegree(n);
        }
        t.differentiate();
        t *= ip;
        t >>= lz;
        if (t.getDegree() > n+1){
            t.setDegree(n+1);
        }
    }
    return res;
}

template <int m>
FormalPowerSeries<m> naiveCompose(FormalPowerSeries<m> lhs, FormalPowerSeries<m> rhs){
    int n = lhs.getDegree();
    FormalPowerSeries<m> g(1);
    FormalPowerSeries<m> res(0);
    for (int i=0;i<=n;i++){
        res += g * lhs.getCoefficient(i);
        g *= rhs;
        g.setDegree(n);
    }
    return res;
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
        res[i] = r.getCoefficient(i);
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
    f.setDegree(n);
    for (int i=1;i<=n;i++){
        for (int j=1;i*j<=n;j++){
            f.setCoefficient(i*j, f.getCoefficient(i*j) + invn[j]);
        }
    }
    f = exponential<m>(f);
    for (int i=0;i<=n;i++){
        invn[i] = f.getCoefficient(i);
    }
    return invn;
}

template <int m>
vector<modint<m>> getStirlingNumberOfFirstTable(int n){
    FormalPowerSeries<m> g = 1;
    FormalPowerSeries<m> h(vector<modint<m>>({0, 1}));
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
        res[i] = g.getCoefficient(i);
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
        a[i] = f.getCoefficient(i);
    }
    return a;
}

template <int m>
vector<modint<m>> getSubsetSum(vector<modint<m>> v){
    int t = v.size() - 1;
    vector<modint<m>> invn(t+1, 1);
    for (int i=2;i<=t;i++){
        invn[i] = -invn[m%i]*(m/i);
    }
    modint<m> offset = modpow(2, v[0].val(), m), sign = 1;
    vector<modint<m>> res(t+1);
    for (int i=1;i<=t;i++){
        sign = 1;
        for (int j=1;i*j<=t;j++){
            res[i*j] += sign * invn[j] * v[i];
            sign *= m - 1;
        }
    }
    FormalPowerSeries<m> fps(res);
    fps = exponential(fps);
    fps *= offset;
    for (int i=0;i<=t;i++){
        res[i] = fps.getCoefficient(i);
    }
    return res;
}

constexpr const int mod = 924844033;
using mint = modint<mod>;
using FPS = FormalPowerSeries<mod>;
