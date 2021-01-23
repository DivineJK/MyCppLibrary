#include <bits/stdc++.h>
using namespace std;
template <typename T1, typename T2, typename T3>
T1 safety_divide(T2 a, T3 b, bool mod=false){
    T1 r, m;
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
template <typename T_out, typename T_bas, typename T_pow, typename T_mod>
T_out modpow(T_bas n, T_pow m, T_mod modulo=0){
    T_out y = 1;
    T_pow tmp = m;
    T_out bas = n;
    while (tmp){
        if (tmp & 1){
            y *= bas;
            if (modulo){
                y = safety_divide<T_out>(y, modulo, true);
            }
        }
        bas *= bas;
        if (modulo){
            bas = safety_divide<T_out>(bas, modulo, true);
        }
        tmp >>= 1;
    }
    return y;
}
