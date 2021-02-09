#include <bits/stdc++.h>
using namespace std;
struct discrete_log{
    int prev_num;
    int prev_mod;
    int prev_inv;
    int prev_sqrt;
    unordered_map<int, int> baby_dict;
    discrete_log(){
        prev_num = -1;
        prev_mod = -1;
        prev_inv = -1;
        prev_sqrt = -1;
        baby_dict = {};
    }
    template<typename T_res, typename T1, typename T2>
    T_res safety_divide(T1 a, T2 b, bool mod=false){
        T_res r, m;
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
    int modpow(int n, int m, int modulo=0){
        long long y = (long long)1;
        long long bas = (long long)n;
        while (m){
            if (m & 1){
                y *= bas;
                if (modulo){
                    y %= modulo;
                }
            }
            bas *= bas;
            if (modulo){
                bas %= modulo;
            }
            m >>= 1;
        }
        return (int)y;
    }
    int gcd(int a, int b){
        int tmp;
        while (b){
            tmp = a;
            a = b;
            b = safety_divide<int>(tmp, b, true);
        }
        return a;
    }
    pair<int, int> extgcd(int a, int b, int c){
        if (b == 0){
            if (a == 0){
                if (c == 0){
                    return make_pair(0, 0);
                }
                return make_pair(-1, -1);
            }
            if (safety_divide<int>(c, a, true)){
                return make_pair(-1, -1);
            }
            return make_pair(safety_divide<int>(c, a), 0);
        }
        if (a < 0){
            a = -a;
            b = -b;
            c = -c;
        }
        int g = gcd(a, b);
        if (safety_divide<int>(c, g, true)){
            return make_pair(-1, -1);
        }
        a = safety_divide<int>(a, g);
        b = safety_divide<int>(b, g);
        c = safety_divide<int>(c, g);
        int tmp1, tmp2, tmp3;
        int x = 1, y = 0, u = 0, v = 1, k = a, l = b;
        while (l){
            tmp1 = x;
            tmp2 = y;
            tmp3 = k;
            x = u;
            y = v;
            u = tmp1 - u * safety_divide<int>(k, l);
            v = tmp2 - v * safety_divide<int>(k, l);
            k = l;
            l = safety_divide<int>(tmp3, l, true);
        }
        k = safety_divide<int>((long long)x*c, b);
        x = safety_divide<int>((long long)x*c, b, true);
        y *= c;
        y += k * a;
        return make_pair(x, y);
    }
    int ext_inved(int a, int c, int modulo){
        return safety_divide<int>(extgcd(a, modulo, c).first, modulo, true);
    }
    vector<pair<int, int>> factorization(int n){
        vector<pair<int, int>> res;
        int a = n;
        int p = 3;
        int cnt = 0;
        while (a % 2 == 0){
            a >>= 1;
            cnt += 1;
        }
        if (cnt){
            res.push_back(make_pair(2, cnt));
        }
        while (a != 1){
            cnt = 0;
            while (a % p == 0){
                cnt += 1;
                a /= p;
            }
            if (cnt){
                res.push_back(make_pair(p, cnt));
            }
            p += 2;
            if (p * p > n && a != 1){
                res.push_back(make_pair(a, 1));
                break;
            }
        }
        return res;
    }
    vector<int> divisors(int n){
        vector<int> res;
        int p = 1;
        int cnt = 0;
        while (p * p <= n){
            if (n % p == 0){
                res.push_back(p);
                cnt += 1;
            }
            p += 1;
        }
        for (int i=cnt;i>0;i--){
            if (res[i-1]*res[i-1] < n){
                res.push_back(n/res[i-1]);
            }
        }
        return res;
    }
    int totient(int n){
        int a = n;
        int res = n;
        int p = 2;
        while (a != 1){
            if (a % p == 0){
                res /= p;
                res *= p-1;
                while (a % p == 0){
                    a /= p;
                }
            }
            p += 1;
            if (p * p > n && a != 1){
                res /= a;
                res *= a-1;
                break;
            }
        }
        return res;
    }
    int generalized_bsgs(int n, int x, int y){
        if (x == 0){
            if (y == 0){
                if (n == 1){
                    return 0;
                }
                return 1;
            }
            if (y == 1){
                if (n == 1){
                    return -1;
                }
                return 0;
            }
            return -1;
        }
        vector<pair<int, int>> fc = factorization(x);
        int Mp = n;
        int tail = 0;
        int cnt = 0;
        for (pair<int, int> p: fc){
            cnt = 0;
            while (Mp % p.first == 0){
                cnt += 1;
                Mp /= p.first;
            }
            if (tail < (p.second+cnt-1)/p.second){
                tail = (p.second+cnt-1)/p.second;
            }
        }
        int fMp = totient(Mp);
        vector<int> div_fMp = divisors(fMp);
        int bas = 1;
        for (int i=0;i<tail;i++){
            if (y == bas){
                return i;
            }
            bas = (int)((long long)bas*x%n);
        }
        if (y % gcd(bas, n)){
            return -1;
        }
        int loop = fMp;
        for (int k: div_fMp){
            if ((long long)bas*modpow(x, k, n)%n == (long long)bas){
                loop = k;
                break;
            }
        }
        int b = prev_inv;
        int m = prev_sqrt;
        int e = modpow(x, loop, n);
        if (n != prev_mod || x != prev_num){
            baby_dict = {};
            prev_num = x;
            prev_mod = n;
            int l = 0, r = loop;
            m = (l + r) / 2;
            while (r - l > 1){
                if (((long long)m * m) <= loop){
                    l = m;
                }
                else{
                    r = m;
                }
                m = (l + r) / 2;
            }
            if (m * m < loop){
                m += 1;
            }
            prev_sqrt = m;
            b = modpow(ext_inved(x, e, n), m, n);
            prev_inv = b;
            int f = bas;
            for (int i=0;i<m;i++){
                baby_dict[f] = i;
                f = (int)((long long)f*x%n);
            }
        }
        int g = y;
        for (int i=0;i<m;i++){
            if (baby_dict.count(g)){
                return i*m+baby_dict[g]+tail;
            }
            g = (int)((long long)g*b%n);
        }
        return -1;
    }
};
