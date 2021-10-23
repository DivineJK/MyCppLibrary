#include <bits/stdc++.h>
using namespace std;

template <typename T>
void writeVector(vector<T> v){
    if (v.size() <= 0){
        return;
    }
    cout << v[0];
    for (int i=1;i<v.size();i++){
        cout << " " << v[i];
    }
    cout << "\n";
}

#pragma mark - ModInt

// ModInt
class RuntimeModInt{
    protected:
    static int& mod() {
        static int _m = 0;
        return _m;
    }
    public:
    unsigned int x = 0;
    static void setModulo(const int m) { mod() = m; }
    static int getModulo() {return mod();}
    const unsigned int& getValue() const{
        return x;
    }
    RuntimeModInt(){}
    RuntimeModInt(const RuntimeModInt& a){
        x = a.x;
    }
    template <typename T>
    RuntimeModInt(T a){
        long long z = (long long)(a % mod());
        z = (z >= 0 ? z : z + mod());
        x = (unsigned int)z;
    }
    RuntimeModInt(bool a){
        x = (uint32_t)a;
    }
    RuntimeModInt inv() const{
        RuntimeModInt res = 1;
        int p = x;
        int mm = mod();
        while (p > 1){
            res *= (mm - mm / p);
            p = mm % p;
        }
        return res;
    }
    bool operator==(const RuntimeModInt& rhs) const{
        return x == rhs.x;
    }
    template <typename T>
    bool operator==(const T rhs) const{
        RuntimeModInt res = RuntimeModInt(rhs);
        return x == res.x;
    }
    bool operator!=(const RuntimeModInt& rhs) const{
        return x != rhs.x;
    }
    template <typename T>
    bool operator!=(const T rhs) const{
        RuntimeModInt res = RuntimeModInt(rhs);
        return x != res.x;
    }
    friend ostream& operator<<(ostream& os, const RuntimeModInt& rhs){
        os << rhs.x;
        return os;
    }
    friend istream& operator>>(istream& ist, RuntimeModInt& rhs){
        unsigned int s;
        ist >> s;
        rhs = s;
        return (ist);
    }
    RuntimeModInt operator+() const{
        RuntimeModInt res = RuntimeModInt(*this);
        return res;
    }
    RuntimeModInt operator-() const{
        RuntimeModInt res = RuntimeModInt(*this);
        if (res.x){
            res.x = mod() - res.x;
        }
        return res;
    }
    RuntimeModInt& operator++(){
        x++;
        if (x == mod()){
            x = 0;
        }
        return *this;
    }
    RuntimeModInt& operator--(){
        if (x == 0){
            x = mod();
        }
        x--;
        return *this;
    }
    RuntimeModInt operator++(int){
        RuntimeModInt res = RuntimeModInt(*this);
        ++*this;
        return res;
    }
    RuntimeModInt operator--(int){
        RuntimeModInt res = RuntimeModInt(*this);
        --*this;
        return res;
    }
    RuntimeModInt& operator+=(const RuntimeModInt& rhs){
        x = (x + rhs.x >= mod() ? x + rhs.x - mod() : x + rhs.x);
        return *this;
    }
    RuntimeModInt& operator-=(const RuntimeModInt& rhs){
        x = (x < rhs.x ? x - rhs.x + mod() : x - rhs.x);
        return *this;
    }
    RuntimeModInt& operator*=(const RuntimeModInt& rhs){
        x = (unsigned int)((long long)x * rhs.x % mod());
        return *this;
    }
    RuntimeModInt& operator/=(const RuntimeModInt& rhs){
        x *= rhs.inv().getValue();
        return *this;
    }
    RuntimeModInt operator+(const RuntimeModInt& rhs) const{
        RuntimeModInt res = RuntimeModInt(*this);
        res += rhs;
        return res;
    }
    RuntimeModInt operator-(const RuntimeModInt& rhs) const{
        RuntimeModInt res = RuntimeModInt(*this);
        res -= rhs;
        return res;
    }
    RuntimeModInt operator*(const RuntimeModInt& rhs) const{
        RuntimeModInt res = RuntimeModInt(*this);
        res *= rhs;
        return res;
    }
    RuntimeModInt operator/(const RuntimeModInt& rhs) const{
        RuntimeModInt res = RuntimeModInt(*this);
        res /= rhs;
        return res;
    }
};

template <typename T>
RuntimeModInt modpow(RuntimeModInt n, T p){
    RuntimeModInt y = 1, bas = n;
    while (p){
        if (p & 1){
            y *= bas;
        }
        bas *= bas;
        p >>= 1;
    }
    return y;
}

template <typename T, typename U>
RuntimeModInt modpow(U n, T p){
    return modpow(RuntimeModInt(n), p);
}

int modinv(int a, int m){
    int res = 1;
    int p = a % m;
    while (p > 1){
        res = (long long)res * ((long long)m - m / p) % m;
        p = m % p;
    }
    return res;
}

int modsqrt(RuntimeModInt a){
    int m = RuntimeModInt::getModulo();
    if (m == 2){
        return a.getValue();
    }
    if (a == 0){
        return 0;
    }
    if (modpow(a, m>>1) == -1){
        return -1;
    }
    RuntimeModInt z = 1;
    while (modpow(z, m>>1) == 1){
        z++;
    }
    int q = m - 1;
    RuntimeModInt s = 0;
    while (~q & 1){
        q >>= 1;
        s++;
    }
    z = modpow(z, q);
    RuntimeModInt r = modpow(a, (q+1)>>1);
    RuntimeModInt t = modpow(a, q);
    RuntimeModInt b, i, c;
    while (t != 1){
        c = t;
        i = 0;
        while (c != 1){
            i++;
            c *= c;
        }
        b = modpow(z, 1<<(s.getValue()-i.getValue()-1));
        z = b * b;
        s = i;
        t *= z;
        r *= b;
    }
    if (r.getValue() > m - r.getValue()){
        r = -r;
    }
    return r.getValue();
}

#pragma mark - FormalPowerSeries

static const int DEFAULT_MODULO_COUNT = 3;
static const vector<int> DEFAULT_BASE_THRESHOLD = {21, 22, 23};
static const vector<int> DEFAULT_MODULO_LIST = {924844033, 943718401, 998244353};
static const vector<int> DEFAULT_PRIMITIVE_ROOT = {508059997, 803508810, 128805723};
static const vector<int> DEFAULT_INVERSE_ROOT = {3597, 175, 31};

class FormalPowerSeriesArbitraryPrimeModulo{
    using FPS = FormalPowerSeriesArbitraryPrimeModulo;
    protected:
    static int moduloCount;
    static vector<int> baseThreshold;
    static vector<int> moduloList;
    static vector<int> primitiveRoot;
    static vector<int> inverseRoot;
    static vector<vector<int>> primitiveBaseList;
    static vector<vector<int>> inverseBaseList;
    static vector<vector<int>> cumulativeBaseList;
    int degree = -1;
    vector<RuntimeModInt> polynomial = {};
    bool isValid = true;
    protected:
    RuntimeModInt solveModularEquation(vector<int> values){
        assert(values.size() == moduloCount);
        int t[moduloCount];
        RuntimeModInt b = 1;
        int nv;
        RuntimeModInt res = 0;
        for (int i=0;i<moduloCount;i++){
            nv = values[i] % moduloList[i];
            for (int j=0;j<i;j++){
                nv += moduloList[i] - t[j];
                nv = (long long)nv * modinv(moduloList[j], moduloList[i]) % moduloList[i];
            }
            res += b * nv;
            b *= moduloList[i];
            t[i] = nv;
        }
        return res;
    }
    vector<RuntimeModInt> getNumberTheoremTransform(vector<RuntimeModInt> f){
        int currentMod = RuntimeModInt::getModulo();
        int d = -1;
        for (int i=0;i<moduloCount;i++){
            if (currentMod == moduloList[i]){
                d = i;
                break;
            }
        }
        if (d == -1){
            assert(false);
        }
        int n = f.size();
        int m = 1;
        int depth = 0;
        while (m < n){
            m <<= 1;
            depth++;
        }
        RuntimeModInt grow = 1;
        int p = 1, t = 1 << (depth - 1), offset = 0, w, z;
        RuntimeModInt x, y;
        for (int i=0;i<depth;i++){
            grow = 1;
            for (int j=0;j<p;j++){
                offset = j << (depth - i);
                for (int k=0;k<t;k++){
                    x = f[k+offset];
                    y = f[k+offset+t] * grow;
                    f[k+offset] = x + y;
                    f[k+offset+t] = x - y;
                }
                w = (j+1)^(j&(j+1)), z = -1;
                while (w){
                    w >>= 1;
                    z++;
                }
                grow *= cumulativeBaseList[d][z];
            }
            p <<= 1;
            t >>= 1;
        }
        return f;
    }
    vector<RuntimeModInt> getInverseNumberTheoremTransform(vector<RuntimeModInt> f){
        int currentMod = RuntimeModInt::getModulo();
        int d = -1;
        for (int i=0;i<moduloCount;i++){
            if (currentMod == moduloList[i]){
                d = i;
                break;
            }
        }
        if (d == -1){
            assert(false);
        }
        int n = f.size();
        int m = 1;
        int depth = 0;
        while (m < n){
            m <<= 1;
            depth++;
        }
        RuntimeModInt grow = 1, seed = 0, x, y;
        int offset = 0;
        for (int i=0;i<depth;i++){
            seed = inverseBaseList[d][i+1];
            grow = 1;
            offset = 1 << i;
            for (int k=0;k<offset;k++){
                for (int j=k;j<n;j+=1<<(i+1)){
                    x = f[j];
                    y = f[j+offset] * grow;
                    f[j] = x + y;
                    f[j+offset] = x - y;
                }
                grow *= seed;
            }
        }
        RuntimeModInt mInv = m;
        mInv = mInv.inv();
        for (int i=0;i<m;i++){
            f[i] *= mInv;
        }
        return f;
    }
    vector<RuntimeModInt> convolveOne(vector<RuntimeModInt> a, vector<RuntimeModInt> b){
        int n = a.size(), m = b.size();
        int bt = 1;
        while (bt < n + m){
            bt <<= 1;
        }
        vector<RuntimeModInt> x(bt), y(bt);
        if (n == 0 && m == 0){
            return {};
        }
        if ((long long)m * n < 100000){
            for (int i=0;i<n;i++){
                for (int j=0;j<m;j++){
                    x[i+j] += a[i] * b[j];
                }
            }
            return x;
        }
        for (int i=0;i<n;i++){
            x[i] = a[i];
        }
        for (int i=0;i<m;i++){
            y[i] = b[i];
        }
        x = getNumberTheoremTransform(x);
        y = getNumberTheoremTransform(y);
        for (int i=0;i<bt;i++){
            x[i] *= y[i];
        }
        x = getInverseNumberTheoremTransform(x);
        x.resize(n+m-1);
        return x;
    }
    vector<RuntimeModInt> convolve(vector<RuntimeModInt> a, vector<RuntimeModInt> b){
        int tmod = RuntimeModInt::getModulo();
        for (int i=0;i<moduloCount;i++){
            if (tmod == moduloList[i]){
                return convolveOne(a, b);
            }
        }
        int n = a.size(), m = b.size();
        if (n == 0 && m == 0){
            return {};
        }
        if ((long long)m * n < 100000){
            vector<RuntimeModInt> x(n+m-1);
            for (int i=0;i<n;i++){
                for (int j=0;j<m;j++){
                    x[i+j] += a[i] * b[j];
                }
            }
            return x;
        }
        vector<vector<int>> conv(moduloCount, vector<int>(n+m-1, 0));
        for (int i=0;i<moduloCount;i++){
            RuntimeModInt::setModulo(moduloList[i]);
            vector<RuntimeModInt> v = convolveOne(a, b);
            for (int j=0;j<n+m-1;j++){
                conv[i][j] = v[j].getValue();
            }
        }
        vector<int> temp(moduloCount);
        RuntimeModInt::setModulo(tmod);
        vector<RuntimeModInt> res(n+m-1, 0);
        for (int i=0;i<n+m-1;i++){
            for (int j=0;j<moduloCount;j++){
                temp[j] = conv[j][i];
            }
            res[i] = solveModularEquation(temp);
        }
        return res;
    }
    FPS& convertToMonomialMultiplication(int deg){
        regularize();
        if (deg == 0){
            return *this;
        }
        if (getDegree() == -1){
            return *this;
        }
        int d = getDegree();
        setDegree(d+deg);
        for (int i=d;i>=0;i--){
            polynomial[i+deg] = polynomial[i];
        }
        for (int i=0;i<deg;i++){
            polynomial[i] = 0;
        }
        return *this;
    }
    FPS getMonomialMultiplication(int deg){
        FPS res = FPS(*this);
        res.convertToMonomialMultiplication(deg);
        return res;
    }
    FPS& convertToMonomialDivision(int deg){
        regularize();
        if (deg > getDegree()){
            setDegree(-1);
            return *this;
        }
        int d = getDegree();
        for (int i=deg;i<=d;i--){
            polynomial[i-deg] = polynomial[i];
        }
        setDegree(d-deg);
        regularize();
        return *this;
    }
    FPS getMonomialDivision(int deg){
        FPS res = FPS(*this);
        res.convertToMonomialDivision(deg);
        return res;
    }
    public:
    static void initialize(){
        moduloCount = DEFAULT_MODULO_COUNT;
        moduloList = DEFAULT_MODULO_LIST;
        baseThreshold = DEFAULT_BASE_THRESHOLD;
        primitiveRoot = DEFAULT_PRIMITIVE_ROOT;
        inverseRoot = DEFAULT_INVERSE_ROOT;
        primitiveBaseList.resize(moduloCount);
        inverseBaseList.resize(moduloCount);
        primitiveBaseList.resize(moduloCount);
        cumulativeBaseList.resize(moduloCount);
        for (int i=0;i<moduloCount;i++){
            primitiveBaseList[i].resize(baseThreshold[i]+1);
            inverseBaseList[i].resize(baseThreshold[i]+1);
            primitiveBaseList[i][baseThreshold[i]] = primitiveRoot[i];
            inverseBaseList[i][baseThreshold[i]] = inverseRoot[i];
            for (int j=baseThreshold[i];j>0;j--){
                primitiveBaseList[i][j-1] = (long long)primitiveBaseList[i][j] * primitiveBaseList[i][j] % moduloList[i];
                inverseBaseList[i][j-1] = (long long)inverseBaseList[i][j] * inverseBaseList[i][j] % moduloList[i];
            }
        }
        for (int i=0;i<moduloCount;i++){
            cumulativeBaseList[i].resize(baseThreshold[i]-1);
            int ie = 1;
            for (int j=0;j<baseThreshold[i]-1;j++){
                cumulativeBaseList[i][j] = (long long)ie * primitiveBaseList[i][j+2] % moduloList[i];
                ie = (long long)ie * inverseBaseList[i][j+2] % moduloList[i];
            }
        }
    }
    void regularize(){
        int t = -1, n = polynomial.size() - 1;
        for (int i=n;i>=0;i--){
            if (polynomial[i] != 0){
                t = i;
                break;
            }
        }
        if (t != n){
            polynomial.resize(t+1);
        }
    }
    FormalPowerSeriesArbitraryPrimeModulo(){
    }
    FormalPowerSeriesArbitraryPrimeModulo(const vector<RuntimeModInt> a){
        setDegree(a.size()-1);
        for (int i=0;i<=degree;i++){
            polynomial[i] = a[i];
        }
        regularize();
    }
    FormalPowerSeriesArbitraryPrimeModulo(const RuntimeModInt a){
        setCoefficient(0, a);
        regularize();
    }
    FormalPowerSeriesArbitraryPrimeModulo(const FPS& f){
        setIsValid(f.getIsValid());
        setDegree(f.getDegree());
        for (int i=0;i<=degree;i++){
            polynomial[i] = f.polynomial[i];
        }
        regularize();
    }
    FPS& operator=(const FPS& rhs){
        isValid = rhs.getIsValid();
        degree = rhs.getDegree();
        polynomial = rhs.polynomial;
        return *this;
    }
    void setIsValid(bool b){
        isValid = b;
    }
    bool getIsValid() const{
        return isValid;
    }
    int getDegree() const{
        return degree;
    }
    void setDegree(int d){
        degree = d;
        polynomial.resize(d+1);
    }
    RuntimeModInt getCoefficient(int pos) const{
        if (pos > degree || pos < 0){
            return 0;
        }
        return polynomial[pos];
    }
    void setCoefficient(int pos, RuntimeModInt newCoeff){
        if (pos < 0){
            return;
        }
        if (pos > degree){
            if (newCoeff == 0){
                return;
            }
            setDegree(pos);
        }
        polynomial[pos] = newCoeff;
    }
    template <typename TYPE>
    void setCoefficient(int pos, TYPE newCoeff){
        RuntimeModInt v = RuntimeModInt(newCoeff);
        setCoefficient(pos, v);
    }
    FPS& convertToMonomialShift(int deg){
        if (deg >= 0){
            convertToMonomialMultiplication(deg);
        } else{
            convertToMonomialDivision(-deg);
        }
        return *this;
    }
    FPS getMonomialShift(int deg){
        FPS res = FPS(*this);
        res.convertToMonomialShift(deg);
        return res;
    }
    FPS& convertToMonomialModulo(int deg){
        setDegree(deg-1);
        regularize();
        return *this;
    }
    FPS getMonomialModulo(int deg){
        FPS res = FPS(*this);
        res.convertToMonomialModulo(deg);
        return res;
    }
    friend ostream& operator<<(ostream& os, FPS& fps){
        int fd = fps.getDegree();
        if (!fps.getIsValid()){
            os << -1;
            return os;
        }
        if (fd == -1){
            os << 0;
            return os;
        }
        for (int i=0;i<=fd;i++){
            if (i){
                os << " ";
            }
            os << fps.getCoefficient(i);
        }
        return os;
    }
    friend istream& operator>>(istream& ist, FPS& fps){
        int fd = fps.getDegree();
        for (int i=0;i<=fd;i++){
            RuntimeModInt s;
            ist >> s;
            fps.setCoefficient(i, s);
        }
        return ist;
    }
    bool operator==(const FPS& rhs) const{
        int n = getDegree(), m = rhs.getDegree();
        int s, t;
        if (n > m){
            s = n;
            t = m;
        } else{
            s = m;
            t = n;
        }
        for (int i=0;i<=t;i++){
            if (getCoefficient(i) != rhs.getCoefficient(i)){
                return false;
            }
        }
        if (s == n){
            for (int i=t+1;i<=s;i++){
                if (getCoefficient(i) != 0){
                    return false;
                }
            }
        } else{
            for (int i=t+1;i<=s;i++){
                if (rhs.getCoefficient(i) != 0){
                    return false;
                }
            }
        }
        return true;
    }
    bool operator!=(const FPS& rhs) const{
        return !(*this == rhs);
    }
    FPS operator+(){
        FPS res = FPS(*this);
        return res;
    }
    FPS operator-(){
        FPS res = FPS(*this);
        for (int i=0;i<=res.degree;i++){
            res.setCoefficient(i, -res.getCoefficient(i));
        }
        return res;
    }
    FPS& operator++(){
        if (getDegree() == -1){
            setDegree(0);
            setCoefficient(0, 1);
            return *this;
        }
        polynomial[0]++;
        return *this;
    }
    FPS operator++(int){
        FPS res = FPS(*this);
        ++*this;
        return res;
    }
    FPS& operator--(){
        if (polynomial.size() == 0){
            setDegree(0);
            setCoefficient(0, -1);
            return *this;
        }
        polynomial[0]--;
        return *this;
    }
    FPS operator--(int){
        FPS res = FPS(*this);
        --*this;
        return res;
    }
    FPS& operator+=(const FPS& rhs){
        int d = rhs.getDegree();
        if (d > getDegree()){
            setDegree(d);
        }
        for (int i=0;i<=d;i++){
            polynomial[i] += rhs.getCoefficient(i);
        }
        regularize();
        return *this;
    }
    FPS& operator+=(const int rhs){
        if (degree == -1){
            setDegree(0);
            polynomial[0] += RuntimeModInt(rhs);
            regularize();
            return *this;
        }
        polynomial[0] += RuntimeModInt(rhs);
        regularize();
        return *this;
    }
    FPS& operator+=(const RuntimeModInt rhs){
        if (degree == -1){
            setDegree(0);
            polynomial[0] += rhs;
            regularize();
            return *this;
        }
        polynomial[0] += rhs;
        regularize();
        return *this;
    }
    FPS& operator-=(const FPS& rhs){
        int d = rhs.getDegree();
        if (d > getDegree()){
            setDegree(d);
        }
        for (int i=0;i<=d;i++){
            polynomial[i] -= rhs.getCoefficient(i);
        }
        regularize();
        return *this;
    }
    FPS& operator-=(const int rhs){
        if (degree == -1){
            setDegree(0);
            polynomial[0] -= RuntimeModInt(rhs);
            regularize();
            return *this;
        }
        polynomial[0] += RuntimeModInt(rhs);
        regularize();
        return *this;
    }
    FPS& operator-=(const RuntimeModInt rhs){
        if (degree == -1){
            setDegree(0);
            polynomial[0] += rhs;
            regularize();
            return *this;
        }
        polynomial[0] += rhs;
        regularize();
        return *this;
    }
    FPS& operator*=(const FPS& rhs){
        degree += rhs.getDegree();
        if (degree == -2){
            degree = -1;
        }
        polynomial = convolve(polynomial, rhs.polynomial);
        regularize();
        return *this;
    }
    FPS& operator*=(const int rhs){
        for (int i=0;i<=degree;i++){
            polynomial[i] *= (RuntimeModInt)rhs;
        }
        return *this;
    }
    FPS& operator*=(const RuntimeModInt rhs){
        for (int i=0;i<=degree;i++){
            polynomial[i] *= rhs;
        }
        return *this;
    }
    template <typename TYPE>
    FPS operator+(const TYPE rhs){
        FPS res = FPS(*this);
        res += rhs;
        return res;
    }
    template <typename TYPE>
    FPS operator-(const TYPE rhs){
        FPS res = FPS(*this);
        res -= rhs;
        return res;
    }
    template <typename TYPE>
    FPS operator*(const TYPE rhs){
        FPS res = FPS(*this);
        res *= rhs;
        return res;
    }
    FPS& differentiate(){
        if (degree == -1){
            return *this;
        }
        for (int i=0;i<degree;i++){
            polynomial[i] = polynomial[i+1] * (i + 1);
        }
        setDegree(degree-1);
        regularize();
        return *this;
    }
    FPS differentiated(){
        FPS res = FPS(*this);
        res.differentiate();
        return res;
    }
    FPS& integrate(){
        if (degree == -1){
            return *this;
        }
        vector<RuntimeModInt> invn(degree+2, 1);
        int p = RuntimeModInt::getModulo();
        for (int i=2;i<=degree+1;i++){
            invn[i] = invn[p%i] * (p - p / i);
        }
        setDegree(degree+1);
        for (int i=degree;i>=1;i--){
            polynomial[i] = polynomial[i-1] * invn[i];
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
    FPS& convertToReversedPolynomial(){
        regularize();
        int d = getDegree();
        if (d == -1){
            return *this;
        }
        RuntimeModInt t;
        for (int i=0;2*i<=d;i++){
            t = polynomial[i];
            polynomial[i] = polynomial[d-i];
            polynomial[d-i] = t;
        }
        regularize();
        return *this;
    }
    FPS getReversedPolynomial(){
        FPS res = FPS(*this);
        res.convertToReversedPolynomial();
        return res;
    }
    FPS& convertToInversedPolynomial(){
        if (getDegree() == -1){
            setIsValid(false);
            return *this;
        }
        if (polynomial[0] == 0){
            setIsValid(false);
            return *this;
        }
        FPS res(polynomial[0].inv());
        FPS two(2), tmp(polynomial[0]), h;
        int b = 1;
        while (b <= getDegree()){
            for (int i=0;i<b;i++){
                if (i+b > getDegree()){
                    break;
                }
                tmp.setCoefficient(i+b, polynomial[i+b]);
            }
            b <<= 1;
            h = two - tmp * res;
            h.setDegree(b-1);
            res *= h;
        }
        *this = res;
        regularize();
        return *this;
    }
    FPS getInversedPolynomial(){
        FPS res = FPS(*this);
        res.convertToInversedPolynomial();
        return res;
    }
    FPS& operator/=(const FPS& rhs){
        int rd = rhs.getDegree();
        regularize();
        if (rd == -1){
            setIsValid(false);
            return *this;
        }
        if (rhs.getCoefficient(rd) == 0){
            setIsValid(false);
            return *this;
        }
        if (rd > getDegree()){
            setDegree(-1);
            return *this;
        }
        FPS tmp = rhs;
        FPS fr = getReversedPolynomial();
        FPS rr = tmp.getReversedPolynomial();
        rr.setDegree(getDegree() - rd);
        fr.setDegree(getDegree() - rd);
        fr *= rr.convertToInversedPolynomial();
        fr.convertToMonomialModulo(getDegree() - rd + 1);
        *this = fr.convertToReversedPolynomial();
        return *this;
    }
    FPS& operator/=(const int rhs){
        if (rhs == 0){
            setIsValid(false);
            return *this;
        }
        for (int i=0;i<=getDegree();i++){
            polynomial[i] /= (RuntimeModInt)rhs;
        }
        return *this;
    }
    FPS& operator/=(const RuntimeModInt rhs){
        if (rhs == 0){
            setIsValid(false);
            return *this;
        }
        for (int i=0;i<=getDegree();i++){
            polynomial[i] /= rhs;
        }
        return *this;
    }
    FPS& operator%=(const FPS& rhs){
        int rd = rhs.getDegree();
        FPS res = FPS(*this);
        res /= rhs;
        *this -= res * rhs;
        convertToMonomialModulo(rd);
        regularize();
        return *this;
    }
    FPS& operator%=(const int rhs){
        int rd = 0;
        FPS res = FPS(*this);
        FPS r = FPS(rhs);
        res /= r;
        *this -= r * res;
        convertToMonomialModulo(rd);
        regularize();
        return *this;
    }
    FPS& operator%=(const RuntimeModInt rhs){
        FPS res = FPS(*this);
        int rd = 0;
        FPS r = FPS(rhs);
        res /= r;
        *this -= r * res;
        convertToMonomialModulo(rd);
        regularize();
        return *this;
    }
    template <typename TYPE>
    FPS operator/(const TYPE rhs){
        FPS res = FPS(*this);
        res /= rhs;
        return res;
    }
    template <typename TYPE>
    FPS operator%(const TYPE rhs){
        FPS res = FPS(*this);
        res %= rhs;
        return res;
    }
};

FormalPowerSeriesArbitraryPrimeModulo logarithm(FormalPowerSeriesArbitraryPrimeModulo fps){
    int deg = fps.getDegree();
    FormalPowerSeriesArbitraryPrimeModulo res = fps.differentiated();
    res *= fps.getInversedPolynomial();
    res.setDegree(deg-1);
    res.integrate();
    return res;
}

constexpr const int GLOBAL_MODULO = (int)1e9 + 7;
using FPS = FormalPowerSeriesArbitraryPrimeModulo;
using mint = RuntimeModInt;
int FPS::moduloCount;
vector<int> FPS::moduloList;
vector<int> FPS::baseThreshold;
vector<int> FPS::primitiveRoot;
vector<int> FPS::inverseRoot;
vector<vector<int>> FPS::primitiveBaseList;
vector<vector<int>> FPS::inverseBaseList;
vector<vector<int>> FPS::cumulativeBaseList;

int main(void){
    // !!SET UP FOR FORMAL POWER SERIES BEGIN!! (DO NOT EDIT THIS ZONE!)
    RuntimeModInt::setModulo(GLOBAL_MODULO);
    FPS::initialize();
    // SET UP FOR FORMAL POWER SERIES END
    int N;
    cin >> N;
    vector<RuntimeModInt> A(N);
    for (int i=0;i<N;i++){
        cin >> A[i];
    }
    FPS f(A);
    f = logarithm(f);
    f.setDegree(N-1);
    cout << f << endl;
    return 0;
}
