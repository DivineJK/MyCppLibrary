template <int m>
class modint{
    using mint = modint;
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
    mint inv() const{
        mint res;
        long long y = 1, z = 0, u = 0, v = 1, k = x, l = mod();
        long long tmp1, tmp2, tmp3;
        while (l){
            tmp1 = y;
            tmp2 = z;
            tmp3 = k;
            y = u;
            z = v;
            u = tmp1 - u * (k / l);
            v = tmp2 - v * (k / l);
            k = l;
            l = tmp3 % l;
            if (tmp3 < 0 && l){
                if (k < 0){
                    l += k;
                    u -= y;
                    v -= z;
                } else{
                    l -= k;
                    u += y;
                    v += z;
                }
            }
        }
        res = y;
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
    bool operator!=(const mint& rhs){
        return x != rhs.x;
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
