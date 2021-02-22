#include <bits/stdc++.h>
using namespace std;

const int mod_free = 998244353;
const int root_free = 128805723;
const int invr_free = 31;
int safety_divide(int a, int b, bool mod=false){
    int r, m;
    r = a / b;
    m = a % b;
    if (m < 0){
        if (b > 0){
            m += b;
            r -= 1;
        }
        else{
            m -= b;
            r += 1;
        }
    }
    if (mod){
        return m;
    }
    return r;
}
template<typename T>
void vec_out(vector<T> v, int m=1<<30){
    int n = min((int)v.size(), m);
    printf("%d", v[0]);
    for (int i=1;i<n;i++){
        printf(" %d", v[i]);
    }
    printf("\n");
}
int inved(int a, int MOD=mod_free){
    int x = 1, y = 0, u = 0, v = 1, k = a, l = MOD;
    int tmp1, tmp2, tmp3;
    while (l){
        tmp1 = x;
        tmp2 = y;
        x = u;
        y = v;
        u = (long long)tmp1 - u * safety_divide(k, l);
        v = (long long)tmp2 - v * safety_divide(k, l);
        tmp3 = k;
        k = l;
        l = safety_divide(tmp3, l, true);
    }
    if (x < 0){x+=MOD;}
    return x;
}
vector<int> ntt_free(vector<int> f, int n, int root=root_free){
    int depth = 0;
    int bin_top = 1;
    while (bin_top < n){
        bin_top <<= 1;
        depth++;
    }
    int pos = 0, tmp = 0, d = 0;
    vector<int> res(bin_top, 0);
    for (int i=0;i<n-1;i++){
        res[i] = f[pos];
        tmp = (i+1)&-(i+1);
        d = 0;
        while (tmp){
            d++;
            tmp >>= 1;
        }
        pos ^= (n-(1<<(depth-d)));
    }
    int MOD = mod_free;
    res[n-1] = f[pos];
    vector<int> base_list(24);
    base_list[23] = root;
    for (int i=23;i>0;i--){
        base_list[i-1] = (long long)base_list[i]*base_list[i]%MOD;
    }
    int grow = 1, seed = 1, offset = 1;
    int u, v;
    for (int i=0;i<depth;i++){
        grow = 1;
        seed = base_list[i+1];
        offset = 1 << i;
        for (int k=0;k<offset;k++){
            for (int j=k;j<n;j+=1<<(i+1)){
                u = res[j];
                v = (int)((long long)res[j+offset] * grow % MOD);
                res[j] = (u+v>=MOD ? u+v-MOD : u+v);
                res[j+offset] = (u-v<0 ? u-v+MOD : u-v);
            }
            grow = (int)((long long)grow * seed % MOD);
        }
    }
    return res;
}
vector<int> inv_ntt_free(vector<int> f, int n){
    vector<int> res = ntt_free(f, n, invr_free);
    int MOD = mod_free;
    int x = inved(n);
    for (int i=0;i<n;i++){
        res[i] = (int)((long long)res[i] * x % MOD);
    }
    return res;
}
vector<int> convolute_one(vector<int> f, vector<int> g){
    int n = f.size();
    int m = g.size();
    int MOD = mod_free;
    int bin_top = 1;
    while (bin_top < n + m){
        bin_top <<= 1;
    }
    vector<int> x(bin_top);
    vector<int> y(bin_top);
    for (int i=0;i<n;i++){
        x[i] = f[i];
    }
    for (int i=0;i<m;i++){
        y[i] = g[i];
    }
    x = ntt_free(x, bin_top);
    y = ntt_free(y, bin_top);
    for (int i=0;i<bin_top;i++){
        x[i] = (int)((long long)x[i] * y[i] % MOD);
    }
    return inv_ntt_free(x, bin_top);
}
vector<int> inverse(vector<int> f, int MOD=mod_free){
    int n = f.size();
    vector<int> g(1, inved(f[0]));
    int bin_top = 1;
    while (bin_top < n){
        bin_top <<= 1;
    }
    vector<int> t(bin_top);
    for (int i=0;i<n;i++){
        t[i] = f[i];
    }
    vector<int> x(1, f[0]);
    vector<int> tmp;
    vector<int> h;
    int m = 1;
    while (m < n){
        for (int i=0;i<m;i++){x.push_back(t[i+m]);}
        m <<= 1;
        h = convolute_one(x, g);
        h.resize(m);
        h[0] = (2 - h[0]);
        if (h[0] < 0){h[0] += MOD;}
        for (int i=1;i<m;i++){
            h[i] = (h[i] ? MOD - h[i] : 0);
        }
        g = convolute_one(g, h);
        g.resize(m);
    }
    return g;
}
vector<int> differentiate(vector<int> f, int MOD=mod_free){
    int n = f.size();
    vector<int> res(n);
    for (int i=1;i<n;i++){
        res[i-1] = (long long)f[i] * i % MOD;
    }
    return res;
}
vector<int> integrate(vector<int> f, int MOD=mod_free){
    int n = f.size();
    vector<int> res(n);
    vector<int> invn(n+1, 1);
    for (int i=0;i<n-1;i++){
        res[i+1] = (long long)f[i] * invn[i+1] % MOD;
        invn[i+2] = (long long)(MOD - invn[MOD%(i+2)]) * (MOD/(i+2)) % MOD;
    }
    return res;
}
vector<int> logarithm(vector<int> f, int MOD=mod_free){
    int n = f.size();
    vector<int> res = differentiate(f);
    res = convolute_one(res, inverse(f));
    res.resize(n);
    res = integrate(res);
    return res;
}
vector<int> exponential(vector<int> f, int MOD=mod_free){
    int n = f.size();
    int bin_top = 1;
    while (bin_top < n){
        bin_top <<= 1;
    }
    vector<int> t(bin_top);
    for (int i=0;i<n;i++){
        t[i] = f[i];
    }
    vector<int> g(1, 1);
    vector<int> p(1, 1);
    int m = 1;
    vector<int> h, w, q, r, x;
    while (2*m <= bin_top){
        h = convolute_one(g, p);
        h.resize(m);
        h[0] = (2 - h[0] < 0 ? 2 - h[0] + MOD : 2 - h[0]);
        for (int i=1;i<m;i++){
            h[i] = (h[i] ? MOD - h[i] : 0);
        }
        g = convolute_one(g, h);
        g.resize(m);
        q.resize(m);
        for (int i=0;i<m;i++){
            q[i] = t[i];
        }
        q = differentiate(q);
        r = convolute_one(p, q);
        r.resize(2*m);
        for (int i=1;i<m;i++){
            r[i-1] = (long long)i * p[i] % MOD - r[i-1];
            r[i-1] = (r[i-1] < 0 ? r[i-1] + MOD : r[i-1]);
        }
        for (int i=m;i<2*m;i++){
            r[i-1] = (r[i-1] ? MOD - r[i-1] : 0);
        }
        w = convolute_one(g, r);
        w.resize(2*m);
        for (int i=0;i<m;i++){
            w[i] = (w[i] + q[i] >= MOD ? w[i] + q[i] - MOD : w[i] + q[i]);
        }
        w = integrate(w);
        for (int i=0;i<2*m;i++){
            w[i] = (t[i] > w[i] ? t[i] - w[i] : t[i] - w[i] + MOD);
        }
        x = convolute_one(p, w);
        for (int i=0;i<m;i++){
            p[i] = (p[i] + x[i] >= MOD ? p[i] + x[i] - MOD : p[i] + x[i]);
            p.push_back(x[i+m]);
        }
        m <<= 1;
    }
    return p;
}
vector<int> pow_of_fps(vector<int> f, int m, int MOD=mod_free){
    int n = f.size();
    int nonzero_pos = 0;
    while (f[nonzero_pos] == 0){
        nonzero_pos++;
        if (nonzero_pos == n){
            break;
        }
    }
    vector<int> res(n, 0);
    if ((long long)nonzero_pos*m >= n){
        return res;
    }
    int v = inved(f[nonzero_pos]);
    int b = f[nonzero_pos];
    for (int i=nonzero_pos;i<n;i++){
        res[i-nonzero_pos] = (long long)v * f[i] % MOD;
    }
    res = logarithm(res);
    for (int i=0;i<n;i++){
        res[i] = (long long)m * res[i] % MOD;
    }
    int c = 1, tmp = nonzero_pos * m;
    while (m){
        if (m & 1){
            c = (long long)c * b % MOD;
        }
        b = (long long) b * b % MOD;
        m >>= 1;
    }
    res = exponential(res);
    for (int i=0;i<n;i++){
        res[i] = (long long)res[i] * c % MOD;
    }
    for (int i=n;i>tmp;i--){
        res[i-1] = res[i-tmp-1];
    }
    for (int i=0;i<tmp;i++){
        res[i] = 0;
    }
    return res;
}
vector<int> polynomial_taylor_shift(vector<int> f, int c, int MOD=mod_free){
    int n = f.size();
    vector<int> invn(n+1, 1);
    int fact = 1;
    int invf = 1;
    vector<int> invf_list(n+1, 1);
    int tmp = 1;
    vector<int> x(n);
    vector<int> y(n);
    for (int i=0;i<n;i++){
        x[n-i-1] = (long long)f[i] * fact % MOD;
        y[i] = (long long)tmp * invf % MOD;
        if (i){
            invn[i+1] = (long long)(MOD - invn[MOD%(i+1)]) * (MOD/(i+1)) % MOD;
        }
        fact = (long long)fact * (i+1) % MOD;
        invf = (long long)invf * invn[i+1] % MOD;
        invf_list[i+1] = invf;
        if (c >= 0){
            tmp = (long long)tmp * c % MOD;
        }
        else{
            tmp = (long long)tmp * (-c) % MOD;
            tmp = (tmp ? MOD - tmp : 0);
        }
    }
    x = convolute_one(x, y);
    x.resize(n);
    for (int i=0;2*i<n;i++){
        tmp = x[i];
        x[i] = x[n-i-1];
        x[n-i-1] = tmp;
    }
    for (int i=0;i<n;i++){
        x[i] = (long long)x[i] * invf_list[i] % MOD;
    }
    return x;
}
int main(void){
    // Your code here!
    int N, M;
    cin >> N >> M;
    vector<int> A(N);
    for (int i=0;i<N;i++){
        cin >> A[i];
    }
    vector<int> X = polynomial_taylor_shift(A, M);
    vec_out(X, N);
    return 0;
}
