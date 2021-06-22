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
#include <cassert>
using namespace std;

class Quaternion{
    public:
    double r;
    double i;
    double j;
    double k;
    Quaternion(){
        r = 0.0f;
        i = 0.0f;
        j = 0.0f;
        k = 0.0f;
    }
    template <typename T>
    Quaternion(complex<T> x){
        r = (double)x.real;
        i = (double)x.imag;
        j = 0.0f;
        k = 0.0f;
    }
    template <typename T>
    Quaternion(T x){
        r = (double)x;
        i = 0.0f;
        j = 0.0f;
        k = 0.0f;
    }
    bool operator==(const Quaternion rhs) const{
        return (r == rhs.r && i == rhs.i && j == rhs.j && k == rhs.k);
    }
    template <typename T>
    bool operator==(const T rhs) const{
        return (r == rhs && i == 0.0f && j == 0.0f && k == 0.0f);
    }
    template <typename T>
    bool operator!=(const T rhs) const{
        return !(*this == rhs);
    }
    friend ostream& operator<<(ostream& os, const Quaternion q){
        os << "quaternion(" << q.r << ", " << q.i << ", " << q.j << ", " << q.k << ")";
        return os;
    }
    Quaternion conjugated() const{
        Quaternion res = *this;
        res.i = -i;
        res.j = -j;
        res.k = -k;
        return res;
    }
    Quaternion& conjugate(){
        i = -i;
        j = -j;
        k = -k;
        return *this;
    }
    double norm2() const{
        return r*r + i*i + j*j + k*k;
    }
    double norm() const{
        return sqrt((*this).norm2());
    }
    Quaternion operator+() const{
        Quaternion res = *this;
        return res;
    }
    Quaternion operator-() const{
        Quaternion res = *this;
        res.r = -res.r;
        res.i = -res.i;
        res.j = -res.j;
        res.k = -res.k;
        return res;
    }
    Quaternion& operator+=(const Quaternion rhs){
        r += rhs.r;
        i += rhs.i;
        j += rhs.j;
        k += rhs.k;
        return *this;
    }
    template <typename T>
    Quaternion& operator+=(T x){
        r += x;
        return *this;
    }
    Quaternion& operator-=(const Quaternion rhs){
        r -= rhs.r;
        i -= rhs.i;
        j -= rhs.j;
        k -= rhs.k;
        return *this;
    }
    template <typename T>
    Quaternion& operator-=(T x){
        r -= x;
        return *this;
    }
    Quaternion& operator*=(const Quaternion rhs){
        double nr = r, ni = i, nj = j, nk = k;
        r = nr*rhs.r - ni*rhs.i - nj*rhs.j - nk*rhs.k;
        i = nr*rhs.i + ni*rhs.r + nj*rhs.k - nk*rhs.j;
        j = nr*rhs.j - ni*rhs.k + nj*rhs.r + nk*rhs.i;
        k = nr*rhs.k + ni*rhs.j - nj*rhs.i + nk*rhs.r;
        return *this;
    }
    template <typename T>
    Quaternion& operator*=(T x){
        r *= x;
        i *= x;
        j *= x;
        k *= x;
        return *this;
    }
    Quaternion& operator/=(const Quaternion rhs){
        assert(rhs.norm2() > 0.0f);
        *this *= rhs.conjugated();
        r /= rhs.norm2();
        i /= rhs.norm2();
        j /= rhs.norm2();
        k /= rhs.norm2();
        return *this;
    }
    template <typename T>
    Quaternion& operator/=(T x){
        assert(x != 0.0f);
        r /= x;
        i /= x;
        j /= x;
        k /= x;
        return *this;
    }
    template <typename T>
    Quaternion operator+(T x) const{
        Quaternion res = *this;
        res += x;
        return res;
    }
    template <typename T>
    Quaternion operator-(T x) const{
        Quaternion res = *this;
        res -= x;
        return res;
    }
    template <typename T>
    Quaternion operator*(T x) const{
        Quaternion res = *this;
        res *= x;
        return res;
    }
    template <typename T>
    Quaternion operator/(T x) const{
        Quaternion res = *this;
        res /= x;
        return res;
    }
};

double abs(Quaternion q){
    return q.norm();
}
