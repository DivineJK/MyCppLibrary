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

template <typename T>
class Rational{
    public:
    T numer;
    T denom;
    protected:
    static T gcd(T a, T b){
        T x = a, y = b;
        T tmp = x;
        while (y){
            tmp = x;
            x = y;
            y = tmp % y;
        }
        return x;
    }
    void regularize(){
        T g = gcd(abs(numer), abs(denom));
        assert(g != 0);
        numer /= g;
        denom /= g;
    }
    public:
    Rational(){
        numer = 0;
        denom = 1;
    }
    Rational(T a){
        numer = a;
        denom = 1;
    }
    Rational(pair<T, T> a){
        numer = a.first;
        denom = a.second;
        if (denom < 0){
            numer = -numer;
            denom = -denom;
        }
        regularize();
    }
    Rational(T a, T b){
        numer = a;
        denom = b;
        if (denom < 0){
            numer = -numer;
            denom = -denom;
        }
        regularize();
    }
    Rational(Rational<T> a, Rational<T> b){
        numer = a.numer * b.denom;
        denom = a.denom * b.numer;
        if (denom < 0){
            numer = -numer;
            denom = -denom;
        }
        regularize();
    }
    Rational(Rational<T> a, T b){
        numer = a.numer;
        denom = a.denom * b;
        if (denom < 0){
            numer = -numer;
            denom = -denom;
        }
        regularize();
    }
    Rational(T a, Rational<T> b){
        numer = a * b.denom;
        denom = b.numer;
        if (denom < 0){
            numer = -numer;
            denom = -denom;
        }
        regularize();
    }
    friend ostream& operator<<(ostream& os, const Rational<T> rhs){
        if (rhs.denom == 0){
            if (rhs.numer > 0){
                os << "inf";
            } else{
                os << "inf";
            }
            return os;
        }
        if (rhs.denom == 1){
            os << rhs.numer;
        } else{
            os << rhs.numer << "/" << rhs.denom;
        }
        return os;
    }
    bool operator==(const Rational<T> rhs) const{
        return (numer == rhs.numer) && (denom == rhs.denom);
    }
    bool operator!=(const Rational<T> rhs) const{
        return !(*this == rhs);
    }
    bool operator>(const Rational<T> rhs) const{
        return (numer * rhs.denom - denom * rhs.numer) > 0;
    }
    bool operator>=(const Rational<T> rhs) const{
        return (*this > rhs) || (*this == rhs);
    }
    bool operator<(const Rational<T> rhs) const{
        return !(*this >= rhs);
    }
    bool operator<=(const Rational<T> rhs) const{
        return !(*this > rhs);
    }
    Rational<T>& operator++(){
        numer += denom;
        return *this;
    }
    Rational<T>& operator--(){
        numer -= denom;
        return *this;
    }
    Rational<T> operator++(int){
        Rational<T> res = *this;
        ++*this;
        return res;
    }
    Rational<T> operator--(int){
        Rational<T> res = *this;
        --*this;
        return res;
    }
    Rational<T> operator+(){
        Rational<T> res = *this;
        return res;
    }
    Rational<T> operator-(){
        Rational<T> res = *this;
        res.numer = -res.numer;
        return res;
    }
    Rational<T>& operator+=(const Rational<T> rhs){
        T n = numer * rhs.denom + denom * rhs.numer;
        T d = denom * rhs.denom;
        assert(n != 0 || d != 0);
        numer = n;
        denom = d;
        regularize();
        return *this;
    }
    Rational<T>& operator-=(const Rational<T> rhs){
        T n = numer * rhs.denom - denom * rhs.numer;
        T d = denom * rhs.denom;
        assert(n != 0 || d != 0);
        numer = n;
        denom = d;
        regularize();
        return *this;
    }
    Rational<T>& operator*=(const Rational<T> rhs){
        T n = numer * rhs.numer;
        T d = denom * rhs.denom;
        assert(n != 0 || d != 0);
        numer = n;
        denom = d;
        regularize();
        return *this;
    }
    Rational<T>& operator/=(const Rational<T> rhs){
        T n = numer * rhs.denom;
        T d = denom * rhs.numer;
        assert(n != 0 || d != 0);
        numer = n;
        denom = d;
        regularize();
        return *this;
    }
    template <typename TYPE>
    Rational<T>& operator+=(const TYPE rhs){
        Rational<T> r = rhs;
        T n = numer * r.denom + denom * r.numer;
        T d = denom * r.denom;
        assert(n != 0 || d != 0);
        numer = n;
        denom = d;
        regularize();
        return *this;
    }
    template <typename TYPE>
    Rational<T>& operator-=(const TYPE rhs){
        Rational<T> r = rhs;
        T n = numer * r.denom - denom * r.numer;
        T d = denom * r.denom;
        assert(n != 0 || d != 0);
        numer = n;
        denom = d;
        regularize();
        return *this;
    }
    template <typename TYPE>
    Rational<T>& operator*=(const TYPE rhs){
        Rational<T> r = rhs;
        T n = numer * r.numer;
        T d = denom * r.denom;
        assert(n != 0 || d != 0);
        numer = n;
        denom = d;
        regularize();
        return *this;
    }
    template <typename TYPE>
    Rational<T>& operator/=(const TYPE rhs){
        Rational<T> r = rhs;
        T n = numer * r.denom;
        T d = denom * r.numer;
        assert(n != 0 || d != 0);
        numer = n;
        denom = d;
        regularize();
        return *this;
    }
    template <typename TYPE>
    Rational<T> operator+(const TYPE rhs) const{
        Rational<T> res = *this;
        res += rhs;
        return res;
    }
    template <typename TYPE>
    Rational<T> operator-(const TYPE rhs) const{
        Rational<T> res = *this;
        res -= rhs;
        return res;
    }
    template <typename TYPE>
    Rational<T> operator*(const TYPE rhs) const{
        Rational<T> res = *this;
        res *= rhs;
        return res;
    }
    template <typename TYPE>
    Rational<T> operator/(const TYPE rhs) const{
        Rational<T> res = *this;
        res /= rhs;
        return res;
    }
};
