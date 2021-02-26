#include <bits/stdc++.h>
using namespace std;

int mod_hand = 998244353;
int root_hand = 128805723;
int invr_hand = 31;
long long safety_divide(long long a, long long b, bool mod=false){
    long long r, m;
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
int modpow(int n, int m, int mod){
    int y = 1, bas = n;
    while (m){
        if (m & 1){
            y = (long long)y * bas % mod;
        }
        bas = (long long)bas * bas % mod;
        m >>= 1;
    }
    return y;
}
bool IsPrime(int n){
    if (n == 1) return false;
    int p = 2;
    while (p * p <= n){
        if (n % p == 0) return false;
        p++;
    }
    return true;
}
template<typename T>
void vec_out(vector<T> v, int l=0, int r=1<<30){
    int n = min((int)v.size(), r);
    if (l >= n){
        printf("\n");
        return;
    }
    printf("%d", v[l]);
    for (int i=l+1;i<n;i++){
        printf(" %d", v[i]);
    }
    printf("\n");
}
int sqrt_mod(int p, int y){
    if (y == 0) return 0;
    if (p == 2) return y & 1;
    int c = modpow(y, (p - 1) / 2, p);
    if (c == p - 1) return -1;
    int q = p - 1;
    int m = 0, b;
    while (q % 2 == 0){
        q >>= 1;
        m++;
    }
    int z = 1;
    while (modpow(z, (p-1)/2, p) == 1){
        z++;
    }
    int t = modpow(y, q, p), r = modpow(y, (q+1)/2, p);
    c = modpow(z, q, p);
    int cnt;
    while (t - 1){
        cnt = 0, b = t;
        while (b - 1){
            b = (long long)b * b % p;
            cnt++;
        }
        b = modpow(c, 1<<(m-cnt-1), p);
        m = cnt, c = (long long)b * b % p, t = (long long) t * c % p, r = (long long)r * b % p;
    }
    return r;
}
int inved(long long a, long long MOD=mod_hand){
    long long x = 1, y = 0, u = 0, v = 1, k = a, l = MOD;
    long long tmp1, tmp2, tmp3;
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
vector<int> ntt_hand(vector<int> f, int n, int root=root_hand){
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
    int MOD = mod_hand;
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
                v = (long long)res[j+offset] * grow % MOD;
                res[j] = (u+v>=MOD ? u+v-MOD : u+v);
                res[j+offset] = (u-v<0 ? u-v+MOD : u-v);
            }
            grow = (long long)grow * seed % MOD;
        }
    }
    return res;
}
vector<int> inv_ntt_hand(vector<int> f, int n){
    vector<int> res = ntt_hand(f, n, invr_hand);
    int MOD = mod_hand;
    int x = inved(n);
    for (int i=0;i<n;i++){
        res[i] = (long long)res[i] * x % MOD;
    }
    return res;
}
vector<int> convolute_one(vector<int> f, vector<int> g){
    int n = f.size();
    int m = g.size();
    int MOD = mod_hand;
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
    x = ntt_hand(x, bin_top);
    y = ntt_hand(y, bin_top);
    for (int i=0;i<bin_top;i++){
        x[i] = (long long)x[i] * y[i] % MOD;
    }
    return inv_ntt_hand(x, bin_top);
}
vector<int> ntt_free(vector<int> f, int n, int MOD, int root){
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
                v = (long long)res[j+offset] * grow % MOD;
                res[j] = (u+v>=MOD ? u+v-MOD : u+v);
                res[j+offset] = (u-v<0 ? u-v+MOD : u-v);
            }
            grow = (long long)grow * seed % MOD;
        }
    }
    return res;
}
vector<int> inv_ntt_free(vector<int> f, int n, int MOD, int invr_free){
    vector<int> res = ntt_free(f, n, MOD, invr_free);
    int x = inved(n, MOD);
    for (int i=0;i<n;i++){
        res[i] = (int)((long long)res[i] * x % MOD);
    }
    return res;
}
vector<int> convolute_general(vector<int> f, vector<int> g, int MOD, int convolute_rank=3, int binary_rank=23, int modulo_limit=(int)1e8){
    int n = f.size();
    int m = g.size();
    int bin_top = 1;
    while (bin_top < n + m){
        bin_top <<= 1;
    }
    vector<int> res(bin_top, 0);
    int s = 0;
    if ((long long)n * m <= 524288){
        for (int i=0;i<n;i++){
            for (int j=0;j<m;j++){
                s = (long long)f[i] * g[j] % MOD;
                res[i+j] = (res[i+j] + s >= MOD ? res[i+j] + s - MOD : res[i+j] + s);
            }
        }
        return res;
    }
    vector<int> x(bin_top);
    vector<int> y(bin_top);
    for (int i=0;i<n;i++){
        x[i] = f[i];
    }
    for (int i=0;i<m;i++){
        y[i] = g[i];
    }
    vector<int> mod_list(convolute_rank, 0);
    vector<int> primitive_root_list(convolute_rank, 0);
    int mod_cnt = 0, r = 0;
    bool flg, fflg;
    int b = (1 << binary_rank);
    int odd_part = (modulo_limit + b - 1) / b;
    if (odd_part & ~1) odd_part++;
    while (mod_cnt < convolute_rank){
        if (IsPrime(odd_part*b+1)){
            r = 1;
            flg = true;
            while (flg){
                fflg = true;
                for (int i=0;i<binary_rank;i++){
                    if (modpow(r, 1<<i, odd_part*b+1) == 1) fflg = false;
                }
                if (modpow(r, 1<<binary_rank, odd_part*b+1) != 1) fflg = false;
                if (fflg){
                    mod_list[mod_cnt] = odd_part*b+1;
                    primitive_root_list[mod_cnt] = r;
                    flg = false;
                    mod_cnt++;
                }
                r++;
            }
        }
        odd_part++;
    }
    vector<vector<int>> res_x(convolute_rank, vector<int>(bin_top, 0));
    vector<vector<int>> res_y(convolute_rank, vector<int>(bin_top, 0));
    for (int i=0;i<convolute_rank;i++){
        res_x[i] = ntt_free(x, bin_top, mod_list[i], inved(primitive_root_list[i], mod_list[i]));
        res_y[i] = ntt_free(y, bin_top, mod_list[i], inved(primitive_root_list[i], mod_list[i]));
        for (int j=0;j<bin_top;j++){
            res_x[i][j] = (long long)res_x[i][j] * res_y[i][j] % mod_list[i];
        }
        res_x[i] = inv_ntt_free(res_x[i], bin_top, mod_list[i], primitive_root_list[i]);
    }
    for (int i=1;i<convolute_rank;i++){
        for (int j=0;j<bin_top;j++){
            for (int k=0;k<i;k++){
                res_x[i][j] = (res_x[i][j] - res_x[k][j] >= 0 ? res_x[i][j] - res_x[k][j] : res_x[i][j] - res_x[k][j] + mod_list[i]);
                res_x[i][j] = (long long)res_x[i][j] * inved(mod_list[k], mod_list[i]) % mod_list[i];
            }
        }
    }
    int tmp = 1;
    for (int i=0;i<convolute_rank;i++){
        for (int j=0;j<bin_top;j++){
            s = (long long)tmp * res_x[i][j] % MOD;
            res[j] = (res[j] + s >= MOD ? res[j] + s - MOD : res[j] + s);
            
        }
        tmp = (long long) tmp * mod_list[i] % MOD;
    }
    return res;
}
vector<int> inverse(vector<int> f){
    int n = f.size();
    int MOD = mod_hand;
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
vector<int> differentiate(vector<int> f){
    int n = f.size();
    vector<int> res(n);
    int MOD = mod_hand;
    for (int i=1;i<n;i++){
        res[i-1] = (long long)f[i] * i % MOD;
    }
    return res;
}
vector<int> integrate(vector<int> f){
    int n = f.size();
    vector<int> res(n);
    int MOD = mod_hand;
    vector<int> invn(n+1, 1);
    for (int i=0;i<n-1;i++){
        res[i+1] = (long long)f[i] * invn[i+1] % MOD;
        invn[i+2] = (long long)(MOD - invn[MOD%(i+2)]) * (MOD/(i+2)) % MOD;
    }
    return res;
}
vector<int> logarithm(vector<int> f){
    int n = f.size();
    int MOD = mod_hand;
    vector<int> res = differentiate(f);
    res = convolute_one(res, inverse(f));
    res.resize(n);
    res = integrate(res);
    return res;
}
vector<int> exponential(vector<int> f){
    int n = f.size();
    int MOD = mod_hand;
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
vector<int> pow_of_fps(vector<int> f, int m){
    int n = f.size();
    int MOD = mod_hand;
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
vector<int> sqrt_of_fps(vector<int> f){
    int n = f.size();
    int MOD = mod_hand;
    vector<int> g = {1};
    vector<int> res(n, 0);
    int nonzero_pos = 0;
    while (f[nonzero_pos] == 0){
        nonzero_pos++;
        if (nonzero_pos == n) break;
    }
    if (nonzero_pos == n) return res;
    if (nonzero_pos & 1) return {-1};
    int s = sqrt_mod(MOD, f[nonzero_pos]);
    if (s == -1) return {-1};
    g[0] = s;
    int bin_top = 1;
    while (bin_top < n){
        bin_top <<= 1;
    }
    vector<int> t(bin_top, 0);
    for (int i=nonzero_pos;i<n;i++){
        t[i-nonzero_pos] = f[i];
    }
    int m = 1, i2 = inved(2);
    vector<int> u(1, t[0]), v;
    while (m <= bin_top){
        v = convolute_one(inverse(g), u);
        for (int i=0;i<m;i++){
            u.push_back(t[i+m]);
        }
        m <<= 1;
        g.resize(m);
        for (int i=0;i<m;i++){
            g[i] = (v[i] + g[i] >= MOD ? v[i] + g[i] - MOD : v[i] + g[i]);
            g[i] = (long long)g[i] * i2 % MOD;
        }
    }
    for (int i=n-nonzero_pos/2;i>0;i--){
        g[i+nonzero_pos/2-1] = g[i-1];
    }
    for (int i=0;i<nonzero_pos/2;i++){
        g[i] = 0;
    }
    return g;
}
vector<int> polynomial_taylor_shift(vector<int> f, int c){
    int n = f.size();
    int MOD = mod_hand;
    vector<int> invn(n+1, 1);
    int fact = 1;
    int invf = 1;
    vector<int> invf_list(n+1, 1);
    int tmp = 1;
    vector<int> x(n);
    vector<int> y(n);
    c = safety_divide(c, MOD, true);
    for (int i=0;i<n;i++){
        x[n-i-1] = (long long)f[i] * fact % MOD;
        y[i] = (long long)tmp * invf % MOD;
        if (i){
            invn[i+1] = (long long)(MOD - invn[MOD%(i+1)]) * (MOD/(i+1)) % MOD;
        }
        fact = (long long)fact * (i+1) % MOD;
        invf = (long long)invf * invn[i+1] % MOD;
        invf_list[i+1] = invf;
        tmp = (long long)tmp * c % MOD;
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
vector<int> subset_sum(vector<int> f, int t){
    vector<int> num_cnt(t+1, 0);
    int MOD = mod_hand;
    for (int i: f) num_cnt[i]++;
    vector<int> invn(t+1, 1);
    for (int i=2;i<=t;i++){
        invn[i] = (long long)(MOD - invn[MOD%i]) * (MOD/i) % MOD;
    }
    vector<int> g(t+1, 0);
    int tmp, sign;
    for (int i=1;i<=t;i++){
        if (num_cnt[i] == 0) continue;
        sign = 1;
        for (int j=1;i*j<=t;j++){
            tmp = (long long)invn[j] * sign % MOD;
            tmp = (long long)tmp * num_cnt[i] % MOD;
            g[i*j] = (g[i*j] + tmp >= MOD ? g[i*j] + tmp - MOD : g[i*j] + tmp);
            sign = MOD - sign;
        }
    }
    return exponential(g);
}
vector<int> first_stirling(int n){
    int MOD = mod_hand;
    vector<int> g = {1};
    vector<int> h = {0, 1};
    int m = 0;
    int b = 1;
    while (m < n){
        if (b & n){
            g = convolute_one(g, polynomial_taylor_shift(h, MOD - m));
            m += b;
            g.resize(m+1);
        }
        h = convolute_one(h, polynomial_taylor_shift(h, MOD - b));
        b <<= 1;
        h.resize(b+1);
    }
    return g;
}
vector<int> second_stirling(int n){
    int MOD = mod_hand;
    vector<int> g(n+1, 0);
    vector<int> h(n+1, 0);
    vector<int> invn(n+3, 1);
    int tmp = 1, sign = 1;
    for (int i=0;i<=n;i++){
        invn[i+2] = (long long)(MOD - invn[MOD%(i+2)]) * (MOD/(i+2)) % MOD;
        g[i] = (long long)sign * tmp % MOD;
        h[i] = (long long)modpow(i, n, MOD) * tmp % MOD;
        tmp = (long long)tmp * invn[i+1] % MOD;
        sign = MOD - sign;
    }
    g = convolute_one(g, h);
    g.resize(n+1);
    return g;
}
vector<int> partition_number(int n){
    int MOD = mod_hand;
    vector<int> invn(n+2, 1);
    for (int i=2;i<n;i++){
        invn[i] = (long long)(MOD - invn[MOD%i]) * (MOD/i) % MOD;
    }
    vector<int> g(n+1, 0);
    for (int i=1;i<=n;i++){
        for (int j=1;i*j<=n;j++){
            g[i*j] = (g[i*j] + invn[j] >= MOD ? g[i*j] + invn[j] - MOD : g[i*j] + invn[j]);
        }
    }
    g = exponential(g);
    return g;
}
int main(void){
    // Your code here!
    int N;
    scanf("%d", &N);
    vector<int> A(N);
    for (int i=0;i<N;i++){
        scanf("%d", &A[i]);
    }
    vector<int> X = exponential(A);
    vec_out(X, 0, N);
    return 0;
}
