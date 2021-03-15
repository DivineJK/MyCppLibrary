#include <bits/stdc++.h>
#include <complex>
using namespace std;

const double pi = acos(-1);

template<typename T>
vector<complex<long double>> simply_fft(vector<T> f, bool inverse=false){
    unsigned int n = f.size();
    unsigned int bin_top = 1, depth = 0;
    while (bin_top < n){
        bin_top <<= 1;
        depth++;
    }
    unsigned int pos = 0, tmp = 0, d = 0;
    vector<complex<long double>> res(bin_top, 0);
    for (int i=0;i<n-1;i++){
        if (pos < n) res[i] = f[pos];
        tmp = (i+1)&-(i+1);
        d = 0;
        while (tmp){
            d++;
            tmp >>= 1;
        }
        pos ^= (n-(1<<(depth-d)));
    }
    if (pos < n) res[n-1] = f[pos];
    complex<long double> grow, seed, u, v;
    unsigned int offset = 1, left = 2;
    for (int i=0;i<depth;i++){
        grow = 1.0+0.0i;
        seed = (inverse ? cos(2.0*pi/left) + sin(2.0*pi/left)*(1.0i) : cos(2.0*pi/left) - sin(2.0*pi/left)*(1.0i));
        offset = 1 << i;
        for (int k=0;k<offset;k++){
            for (int j=k;j<n;j+=1<<(i+1)){
                u = res[j];
                v = res[j+offset] * grow;
                res[j] = u + v;
                res[j+offset] = u - v;
            }
            grow *= seed;
        }
        left <<= 1;
    }
    if (inverse){
        for (int i=0;i<bin_top;i++){
            res[i] /= bin_top;
        }
    }
    return res;
}
template <typename T>
vector<complex<long double>> bluestein_fft(vector<T> f, bool inverse=false){
    unsigned int n = f.size(), bin_top = 1;
    while (bin_top < 2*n+1) bin_top <<= 1;
    vector<complex<long double>> a(bin_top, 0), b(bin_top, 0), c(bin_top, 0);
    int tmp;
    for (int i=0;i<2*n-1;i++){
        tmp = (i-n+1)*(i-n+1);
        b[i] = (inverse ? cos(pi*tmp/n) + sin(pi*tmp/n)*(1.0i) : cos(pi*tmp/n) - sin(pi*tmp/n)*(1.0i));
        c[i] = (inverse ? cos(pi*tmp/n) - sin(pi*tmp/n)*(1.0i) : cos(pi*tmp/n) + sin(pi*tmp/n)*(1.0i));
    }
    c = simply_fft(c);
    for (int i=0;i<n;i++){
        a[i] = f[i];
        a[i] *= b[i+n-1];
    }
    a = simply_fft(a);
    for (int i=0;i<bin_top;i++){
        a[i] *= c[i];
    }
    a = simply_fft(a, true);
    vector<complex<long double>> res(n);
    for (int i=0;i<n;i++){
        res[i] = a[i+n-1] * b[i+n-1];
    }
    if (inverse){
        for (int i=0;i<n;i++){
            res[i] /= n;
        }
    }
    return res;
}
