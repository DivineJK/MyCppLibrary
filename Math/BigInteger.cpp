#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <deque>
#include <bitset>
#include <string>
#include <tuple>
#include <complex>
#include <stdexcept>

using namespace std;

class BigInteger{
    using bint = BigInteger;
    private:
    static const int NTT_MOD = 998244353; // 998244353 = 119 * 2^23 + 1
    static const int PRIMITIVE_ROOT = 31;
    static const int INVERSE_PRIMITIVE_ROOT = 128805723;
    static const int MAX_NTT_DEPTH = 23;
    static const int MAX_INTEGER_DIGITS = 1 << (MAX_NTT_DEPTH - 1);
    public:
    int8_t sign;
    vector<int8_t> num;
    template <typename T>
    static void printVector(vector<T> v){
        if (v.size() <= 0){
            cout << "[]" << endl;
            return;
        }
        cout << "[" << v[0]-48+48;
        for (int i=1;i<v.size();i++){
            cout << ", " << v[i]-48+48;
        }
        cout << "]" <<endl;
    }
    static int modpow(int n, int m, int modulo){
        int y = 1, bas = n;
        while (m){
            if (m & 1){
                y = (long long)y * bas % modulo;
            }
            bas = (long long)bas * bas % modulo;
            m >>= 1;
        }
        return y;
    }
    static vector<int> ntt(vector<int> f, bool inverse=false){
        int n = f.size();
        int bin_top = 1, depth = 0;
        while (bin_top < n){
            bin_top <<= 1;
            depth++;
        }
        vector<int> res(bin_top, 0);
        int pnt = 0, tmp, d;
        for (int i=0;i<bin_top;i++){
            if (pnt < n){
                res[i] = f[pnt];
            }
            tmp = (-(i+1))&(i+1);
            d = 0;
            while (tmp){
                d++;
                tmp >>= 1;
            }
            pnt ^= (bin_top - (bin_top >> d));
        }
        vector<int> bases(MAX_NTT_DEPTH+1, 0);
        bases[MAX_NTT_DEPTH] = (inverse ? PRIMITIVE_ROOT : INVERSE_PRIMITIVE_ROOT);
        for (int i=MAX_NTT_DEPTH;i>0;i--){
            bases[i-1] = (long long)bases[i] * bases[i] % NTT_MOD;
        }
        int u, v, grow, seed, offset;
        for (int i=0;i<depth;i++){
            grow = 1;
            seed = bases[i+1];
            offset = 1 << i;
            for (int k=0;k<offset;k++){
                for (int j=k;j<n;j+=1<<(i+1)){
                    u = res[j];
                    v = (long long)res[j+offset] * grow % NTT_MOD;
                    res[j] = (u+v>=NTT_MOD ? u+v-NTT_MOD : u+v);
                    res[j+offset] = (u-v<0 ? u-v+NTT_MOD : u-v);
                }
                grow = (long long)grow * seed % NTT_MOD;
            }
        }
        if (inverse){
            int ir = modpow(bin_top, NTT_MOD-2, NTT_MOD);
            for (int i=0;i<bin_top;i++){
                res[i] = (long long)res[i] * ir % NTT_MOD;
            }
        }
        return res;
    }
    static vector<int> convolute(vector<int> f, vector<int> g){
        int n = f.size(), m = g.size();
        int bin_top = 1;
        while (bin_top < n + m){
            bin_top <<= 1;
        }
        vector<int> x(bin_top, 0);
        vector<int> y(bin_top, 0);
        for (int i=0;i<n;i++){
            x[i] = f[i];
        }
        for (int i=0;i<m;i++){
            y[i] = g[i];
        }
        x = ntt(x);
        y = ntt(y);
        for (int i=0;i<bin_top;i++){
            x[i] = (long long)x[i] * y[i] % NTT_MOD;
        }
        x = ntt(x, true);
        return x;
    }
    BigInteger(){
        sign = 1;
        num = {};
    }
    BigInteger(string a){
        int n = a.size();
        sign = 1;
        num = vector<int8_t>(0, 0);
        for (int i=n;i>0;i--){
            if ('0' <= a[i-1] && a[i-1] <= '9'){
                num.push_back(a[i-1] - '0');
            } else if (i == 1 && a[i-1] == '-'){
                sign = -1;
            }
        }
        int r = n;
        while (r > 0 && num[r-1] == 0){
            r--;
        }
        if (r == 0){
            sign = 1;
        }
        num.resize(r);
    }
    template <typename T>
    BigInteger(T a){
        num = vector<int8_t>(0, 0);
        sign = 1;
        if (a < 0){
            sign = -1;
            a = -a;
        }
        while (a){
            num.push_back(a%10);
            a /= 10;
        }
    }
    void regularize(){
        int n = num.size();
        while (n > 0 && num[n-1] == 0){
            n--;
        }
        if (n == 0){
            sign = 1;
        }
        num.resize(n);
    }
    string bint_to_string() const{
        int n = num.size();
        if (sign == 1 && n == 0){
            return "0";
        }
        string res = "";
        if (sign == -1){
            res = "-";
        }
        for (int i=n;i>0;i--){
            res = res + to_string(num[i-1]);
        }
        return res;
    }
    friend ostream& operator<<(ostream& os, const bint& rhs){
        os << rhs.bint_to_string();
        return os;
    }
    friend istream& operator>>(istream& ist, bint& rhs){
        string s;
        ist >> s;
        rhs = bint(s);
        return (ist);
    }
    bool operator==(const bint& rhs) const{
        int n = num.size(), m = rhs.num.size();
        if (n != m || sign != rhs.sign){
            return false;
        }
        int l = n;
        while (l > 0 && num[l-1] == rhs.num[l-1]){
            l--;
        }
        if (l == 0){
            return true;
        }
        return false;
    }
    template <typename T>
    bool operator==(const T rhs) const{
        bint x = rhs;
        return *this == x;
    }
    bool operator<(const bint& rhs) const{
        int n = num.size(), m = rhs.num.size();
        if (sign != rhs.sign){
            if (sign == 1){
                return false;
            }
            return true;
        }
        if (n < m){
            if (sign == -1){
                return false;
            }
            return true;
        } else if (n > m){
            if (sign == -1){
                return true;
            }
            return false;
        } else{
            int l = n;
            while (l > 0 && num[l-1] == rhs.num[l-1]){
                l--;
            }
            if (l == 0){
                return false;
            }
            if (num[l-1] < rhs.num[l-1]){
                if (sign == 1){
                    return true;
                }
                return false;
            } else{
                if (sign == 1){
                    return false;
                }
                return true;
            }
        }
    }
    template <typename T>
    bool operator<(const T rhs) const{
        bint x = rhs;
        return *this < x;
    }
    bool operator!=(const bint& rhs) const{
        return !(bint(*this) == rhs);
    }
    bool operator<=(const bint& rhs) const{
        return (bint(*this) < rhs) || (bint(*this) == rhs);
    }
    bool operator>(const bint& rhs) const{
        return !(bint(*this) <= rhs);
    }
    bool operator>=(const bint& rhs) const{
        return !(bint(*this) < rhs);
    }
    bint& operator++(){
        if (sign == 1){
            num.push_back(0);
            num[0]++;
            int pnt = 0;
            while (num[pnt] >= 10){
                num[pnt] -= 10;
                pnt++;
                num[pnt]++;
            }
        } else{
            num.push_back(0);
            num[0]--;
            int pnt = 0;
            while (num[pnt] < 0){
                num[pnt] += 10;
                pnt++;
                num[pnt]--;
            }
        }
        regularize();
        return *this;
    }
    bint& operator--(){
        if (num.size()){
            sign *= -1;
        }
        ++*this;
        sign *= -1;
        regularize();
        return *this;
    }
    bint operator++(int){
        bint result = *this;
        ++*this;
        return result;
    }
    bint operator--(int){
        bint result = *this;
        --*this;
        return result;
    }
    bint operator-(){
        bint result = *this;
        result.sign *= -1;
        result.regularize();
        return result;
    }
    bint& operator+=(const bint& rhs){
        int n = num.size(), m = rhs.num.size();
        if (sign == rhs.sign){
            vector<int> res(max(n, m)+1, 0);
            for (int i=0;i<max(n, m);i++){
                if (i < n){
                    res[i] += num[i];
                }
                if (i < m){
                    res[i] += rhs.num[i];
                }
                if (res[i] >= 10){
                    res[i+1]++;
                    res[i] -= 10;
                }
            }
            num = vector<int8_t>(res.size(), 0);
            for (int i=0;i<(int)res.size();i++){
                num[i] = (int8_t)res[i];
            }
        } else{
            vector<int> larger(max(n, m)+1, 0), smaller(max(n, m)+1, 0);
            bool isLeftLarger;
            if (n > m){
                isLeftLarger = true;
            } else if (n < m){
                isLeftLarger = false;
            } else{
                int pnt = n;
                while (pnt > 0 && num[pnt-1] == rhs.num[pnt-1]){
                    pnt--;
                }
                if (pnt > 0 && num[pnt-1] < rhs.num[pnt-1]){
                    isLeftLarger = false;
                } else{
                    isLeftLarger = true;
                }
            }
            if (isLeftLarger){
                for (int i=0;i<n;i++){
                    larger[i] = num[i];
                }
                for (int i=0;i<m;i++){
                    smaller[i] = rhs.num[i];
                }
            } else{
                sign *= -1;
                for (int i=0;i<m;i++){
                    larger[i] = rhs.num[i];
                }
                for (int i=0;i<n;i++){
                    smaller[i] = num[i];
                }
            }
            for (int i=0;i<max(n, m);i++){
                larger[i] -= smaller[i];
                if (larger[i] < 0){
                    larger[i] += 10;
                    larger[i+1] -= 1;
                }
            }
            num = vector<int8_t>(max(n, m)+1);
            for (int i=0;i<max(n, m)+1;i++){
                num[i] = (int8_t)larger[i];
            }
        }
        regularize();
        return *this;
    }
    template <typename T>
    bint& operator+=(const T rhs){
        bint x = rhs;
        *this += x;
        return *this;
    }
    bint& operator-=(const bint& rhs){
        sign *= -1;
        *this += rhs;
        *this = -*this;
        regularize();
        return *this;
    }
    template <typename T>
    bint& operator-=(const T rhs){
        bint x = rhs;
        *this -= x;
        return *this;
    }
    bint& operator*=(const bint& rhs){
        int n = num.size(), m = rhs.num.size();
        vector<int> left(n, 0);
        vector<int> right(m, 0);
        for (int i=0;i<n;i++){
            left[i] = (int)num[i];
        }
        for (int i=0;i<m;i++){
            right[i] = (int)rhs.num[i];
        }
        vector<int> res = convolute(left, right);
        for (int i=0;i<(int)res.size()-1;i++){
            res[i+1] += res[i] / 10;
            res[i] %= 10;
        }
        num = vector<int8_t>(res.size(), 0);
        for (int i=0;i<(int)res.size();i++){
            num[i] = (int8_t)res[i];
        }
        sign = sign*rhs.sign;
        regularize();
        return *this;
    }
    template <typename T>
    bint& operator*=(const T rhs){
        bint x = rhs;
        *this *= x;
        return *this;
    }
    bint& operator/=(const bint& rhs){
        try{
            if (rhs == 0){
                throw 0;
            }
        } catch(int zde){
            throw runtime_error("Zero Division Exception");
        }
        if (*this == 0){
            return *this;
        }
        int e = 16, loopTime = 0;
        int f = num.size(), d = rhs.num.size();
        while ((e - 1) * (1 << loopTime) + d + 1 <= f){
            loopTime++;
        }
        vector<int> a(num.size(), 0), b(rhs.num.size(), 0);
        for (int i=0;i<(int)num.size();i++){
            a[i] = (int)num[i];
        }
        for (int i=0;i<(int)rhs.num.size();i++){
            b[i] = (int)rhs.num[i];
        }
        vector<int> rem(d+2, 0), tmp(d+2, 0);
        d--;
        rem[d] = 1;
        vector<int> res(e+1, 0);
        bool isDivisible = false, isRemLarge, isZero;
        for (int i=0;i<=e;i++){
            for (int j=9;j>=0;j--){
                for (int k=0;k<=d;k++){
                    tmp[k] = j * b[k];
                }
                for (int k=d+1;k<=d+2;k++){
                    tmp[k] = 0;
                }
                for (int k=0;k<=d+2;k++){
                    tmp[k+1] += tmp[k] / 10;
                    tmp[k] %= 10;
                }
                isRemLarge = true;
                for (int k=d+2;k>=0;k--){
                    if (rem[k] > tmp[k]){
                        isRemLarge = true;
                        break;
                    } else if (rem[k] < tmp[k]){
                        isRemLarge = false;
                        break;
                    }
                }
                if (isRemLarge){
                    res[e-i] = j;
                    break;
                }
            }
            for (int j=0;j<d+2;j++){
                rem[j] -= tmp[j];
                if (rem[j] < 0){
                    rem[j] += 10;
                    rem[j+1]--;
                }
            }
            rem[d+2] -= tmp[d+2];
            isZero = true;
            for (int j=0;j<=d+2;j++){
                if (rem[j]){
                    isZero = false;
                    break;
                }
            }
            if (isZero){
                isDivisible = true;
            }
            for (int j=d+2;j>0;j--){
                rem[j] = rem[j-1];
            }
            rem[0] = 0;
        }
        int t = d + e;
        if (!isDivisible){
            vector<int> two(t+1, 0);
            for (int i=0;i<loopTime;i++){
                two[t] = 2;
                tmp = convolute(res, b);
                for (int j=0;j<t;j++){
                    tmp[j+1] += tmp[j] / 10;
                    tmp[j] %= 10;
                }
                tmp.resize(t+1);
                for (int j=0;j<=t;j++){
                    tmp[j] = two[j] - tmp[j];
                    if (tmp[j] < 0){
                        tmp[j] += 10;
                        tmp[j+1]++;
                    }
                }
                res = convolute(res, tmp);
                for (int j=1;j<(int)res.size();j++){
                    res[j] += res[j-1] / 10;
                    res[j-1] %= 10;
                }
                two[t] = 0;
                t <<= 1;
                two.resize(t+1);
            }
        }
        res = convolute(a, res);
        for (int i=1;i<(int)res.size();i++){
            res[i] += res[i-1] / 10;
            res[i-1] %= 10;
        }
        vector<int> p(res.size()-t, 0);
        for (int i=t;i<(int)res.size();i++){
            p[i-t] = res[i];
        }
        p[0]++;
        int posCarry = 0;
        while (p[posCarry] >= 10){
            p[posCarry] -= 10;
            posCarry++;
            p[posCarry]++;
        }
        vector<int> q = convolute(p, b);
        for (int i=1;i<(int)q.size();i++){
            q[i] += q[i-1] / 10;
            q[i-1] %= 10;
        }
        bool isALarge = true;
        int qSize = q.size(), aSize = a.size();
        int ln, rn;
        for (int i=max(qSize, aSize);i>0;i--){
            ln = 0;
            rn = 0;
            if (i-1 < qSize){
                rn = q[i-1];
            }
            if (i-1 < aSize){
                ln = a[i-1];
            }
            if (ln < rn){
                isALarge = false;
                break;
            } else if (ln > rn){
                break;
            }
        }
        if (!(isALarge)){
            p[0]--;
            posCarry = 0;
            while (p[posCarry] < 0){
                p[posCarry] += 10;
                posCarry++;
                p[posCarry]--;
            }
        }
        q = convolute(p, b);
        for (int i=1;i<(int)q.size();i++){
            q[i] += q[i-1] / 10;
            q[i] %= 10;
        }
        q.resize(a.size());
        bool isRemZero = true;
        for (int i=0;i<(int)a.size();i++){
            if (q[i] != a[i]){
                isRemZero = false;
                break;
            }
        }
        num.clear();
        num.resize(p.size());
        for (int i=0;i<(int)p.size();i++){
            num[i] = (int8_t)p[i];
        }
        if (!(isRemZero) && sign < 0){
            num.push_back(0);
            num[0]++;
        }
        if (rhs.sign * sign < 0){
            sign = -1;
        } else{
            sign = 1;
        }
        regularize();
        return *this;
    }
    template <typename T>
    bint& operator/=(const T rhs){
        bint x = rhs;
        *this /= x;
        return *this;
    }
    bint& operator%=(const bint& rhs){
        try{
            if (rhs == 0){
                throw 0;
            }
        } catch(int zde){
            throw runtime_error("Zero Division Exception");
        }
        bint result = *this;
        result /= rhs;
        result *= rhs;
        *this -= result;
        return *this;
    }
    template <typename T>
    bint& operator%=(const T rhs){
        bint x = rhs;
        *this %= x;
        return *this;
    }
    template <typename T>
    bint operator+(const T rhs){
        bint result = bint(*this);
        result += rhs;
        return result;
    }
    template <typename T>
    bint operator-(const T rhs){
        bint result = bint(*this);
        result -= rhs;
        return result;
    }
    template <typename T>
    bint operator*(const T rhs){
        bint result = bint(*this);
        result *= rhs;
        return result;
    }
    template <typename T>
    bint operator/(const T rhs){
        bint result = bint(*this);
        result /= rhs;
        return result;
    }
    template <typename T>
    bint operator%(const T rhs){
        bint result = bint(*this);
        result %= rhs;
        return result;
    }
};
BigInteger gcd(BigInteger a, BigInteger b){
    BigInteger x = a, y = b, tmp;
    while (y != 0){
        tmp = x;
        x = y;
        y = tmp % y;
    }
    return x;
}
BigInteger lcm(BigInteger a, BigInteger b){
    return a * b / gcd(a, b);
}
BigInteger factorial(int n){
    BigInteger f = 1;
    for (int i=1;i<=n;i++){
        f *= i;
    }
    return f;
}
BigInteger abs(BigInteger a){
    if (a >= 0){
        return a;
    } else{
        return -a;
    }
}
BigInteger pow(BigInteger a, BigInteger b){
    BigInteger y = 1, bas = a;
    while (b != 0){
        if (b % 2 == 1){
            y *= bas;
        }
        bas *= bas;
        b /= 2;
    }
    return y;
}
int main(void){
    // Your code here!
    BigInteger A, B;
    cin >> A >> B;
    BigInteger R = A % B;
    cout << R << endl;
    return 0;
}
