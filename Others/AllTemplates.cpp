// modint
// DijkstraWithPathRestoring
// Dijkstra
// UnionFind
// LazySegmentTree
// FormalPowerSeries

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

// modint
template <int m>
class modint{
    using mint = modint<m>;
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
    unsigned int val(){
        return x;
    }
    mint inv() const{
        mint res = 1;
        int n = x;
        int mm = mod();
        while (n > 1){
            res *= (mm - mm / n);
            n = mm % n;
        }
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
    template <typename T>
    bool operator==(const T rhs){
        mint res = mint(rhs);
        return x == res.x;
    }
    bool operator!=(const mint& rhs){
        return x != rhs.x;
    }
    template <typename T>
    bool operator!=(const T rhs){
        mint res = mint(rhs);
        return x != res.x;
    }
    mint operator+(){
        mint res = mint(*this);
        return res;
    }
    mint operator-(){
        mint res = mint(*this);
        if (res.x){
            res.x = umod() - res.x;
        }
        return res;
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

// Path restorable Dijkstra's algorithm
class DijkstraWithPathRestoring{
    private:
    bool isDirty = false;
    public:
    static constexpr long long inf = (long long)1e18 + 18;
    static constexpr long long threshold = (long long)1e18;
    int vertexCount = 0;
    vector<vector<pair<int, long long>>> graph;
    vector<long long> cost;
    vector<int> prev;
    DijkstraWithPathRestoring(int n){
        vertexCount = n;
        graph.resize(vertexCount);
        cost.resize(vertexCount);
        prev.resize(vertexCount);
        for (int i=0;i<vertexCount;i++){
            graph[i] = vector<pair<int, long long>>(0);
            cost[i] = inf;
            prev[i] = -1;
        }
    }
    void resetAll(){
        for (int i=0;i<vertexCount;i++){
            graph[i] = vector<pair<int, long long>>(0);
            cost[i] = inf;
        }
        isDirty = false;
    }
    void resetOnlyCost(){
        for (int i=0;i<vertexCount;i++){
            cost[i] = inf;
        }
        isDirty = false;
    }
    void connectEdge(int a, int b, long long c, bool undirected=true){
        graph[a].push_back(make_pair(b, c));
        if (undirected){
            graph[b].push_back(make_pair(a, c));
        }
        isDirty = false;
    }
    void getMinimumCost(int start){
        if (isDirty){
            for (int i=0;i<vertexCount;i++){
                cost[i] = inf;
            }
        }
        cost[start] = 0;
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<pair<long long, int>>> pq;
        pq.push(make_pair(0, start));
        int pq_cnt = 1;
        pair<long long, int> b;
        int r, i;
        long long v;
        while (pq_cnt){
            b = pq.top();
            pq.pop();
            pq_cnt--;
            r = b.second;
            if (cost[r] < b.first){
                continue;
            }
            for (pair<int, long long> p: graph[r]){
                i = p.first;
                v = p.second;
                if (cost[i] > cost[r] + v){
                    cost[i] = cost[r] + v;
                    pq.push(make_pair(cost[i], i));
                    pq_cnt++;
                    prev[i] = r;
                }
            }
        }
        isDirty = true;
    }
    vector<int> restorePath(int terminal){
        if (cost[terminal] > threshold){
            return {-1};
        }
        vector<int> res = {terminal};
        terminal = prev[terminal];
        int c = 1, tmp;
        while (terminal != -1){
            res.push_back(terminal);
            terminal = prev[terminal];
            c++;
        }
        for (int i=0;i<c/2;i++){
            tmp = res[i];
            res[i] = res[c-i-1];
            res[c-i-1] = tmp;
        }
        return res;
    }
};

// Dijkstra's algorithm without path restoring
class Dijkstra{
    private:
    bool isDirty = false;
    public:
    static constexpr long long inf = (long long)1e18 + 18;
    static constexpr long long threshold = (long long)1e18;
    int vertexCount = 0;
    vector<vector<pair<int, long long>>> graph;
    vector<long long> cost;
    Dijkstra(int n){
        vertexCount = n;
        graph.resize(vertexCount);
        cost.resize(vertexCount);
        for (int i=0;i<vertexCount;i++){
            graph[i] = vector<pair<int, long long>>(0);
            cost[i] = inf;
        }
    }
    void resetAll(){
        for (int i=0;i<vertexCount;i++){
            graph[i] = vector<pair<int, long long>>(0);
            cost[i] = inf;
        }
        isDirty = false;
    }
    void resetOnlyCost(){
        for (int i=0;i<vertexCount;i++){
            cost[i] = inf;
        }
        isDirty = false;
    }
    void connectEdge(int a, int b, long long c, bool undirected=true){
        graph[a].push_back(make_pair(b, c));
        if (undirected){
            graph[b].push_back(make_pair(a, c));
        }
        isDirty = false;
    }
    void getMinimumCost(int start){
        if (isDirty){
            for (int i=0;i<vertexCount;i++){
                cost[i] = inf;
            }
        }
        cost[start] = 0;
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<pair<long long, int>>> pq;
        pq.push(make_pair(0, start));
        int pq_cnt = 1;
        pair<long long, int> b;
        int r, i;
        long long v;
        while (pq_cnt){
            b = pq.top();
            pq.pop();
            pq_cnt--;
            r = b.second;
            if (cost[r] < b.first){
                continue;
            }
            for (pair<int, long long> p: graph[r]){
                i = p.first;
                v = p.second;
                if (cost[i] > cost[r] + v){
                    cost[i] = cost[r] + v;
                    pq.push(make_pair(cost[i], i));
                    pq_cnt++;
                }
            }
        }
        isDirty = true;
    }
};

// Union-Find (Disjoint Set Union)
class UnionFind{
    public:
    int num;
    vector<int> parent;
    vector<int> size;
    UnionFind(int n){
        num = n;
        parent = vector<int>(num, 0);
        for (int i=0;i<num;i++){
            parent[i] = i;
        }
        size = vector<int>(num, 1);
    }
    int getRoot(int x){
        int t = x;
        while (t != parent[t]){
            t = parent[t];
        }
        int s;
        while (x != parent[x]){
            s = x;
            x = parent[x];
            parent[s] = t;
        }
        return t;
    }
    bool isSame(int x, int y){
        return getRoot(x) == getRoot(y);
    }
    void unite(int x, int y){
        x = getRoot(x);
        y = getRoot(y);
        if (x == y){
            return;
        }
        if (size[x] < size[y]){
            int t = x;
            x = y;
            y = t;
        }
        size[x] += size[y];
        parent[y] = x;
    }
};

// Lazy Segment Tree
template <typename MONOID_TYPE, typename UPDATE_TYPE>
class LazySegmentTree{
    protected:
    MONOID_TYPE binomialOperator(MONOID_TYPE lhs, MONOID_TYPE rhs){}
    MONOID_TYPE updateOperator(int width, UPDATE_TYPE lhs, MONOID_TYPE rhs){}
    UPDATE_TYPE lazyOperator(UPDATE_TYPE lhs, UPDATE_TYPE rhs){}
    public:
    int treeSize = 1;
    int depth = 0;
    MONOID_TYPE monoId;
    UPDATE_TYPE updateId;
    vector<MONOID_TYPE> segmentTree;
    vector<UPDATE_TYPE> updateTree;
    vector<int> zoneTree;
    LazySegmentTree(int n, MONOID_TYPE identity, UPDATE_TYPE upd_id, vector<MONOID_TYPE> initial = {}){
        monoId = identity;
        updateId = upd_id;
        treeSize = 1;
        depth = 0;
        while (treeSize < n){
            treeSize <<= 1;
            depth++;
        }
        segmentTree = vector<MONOID_TYPE>(2*treeSize, identity);
        int m = (treeSize < initial.size() ? treeSize : initial.size());
        for (int i=0;i<m;i++){
            segmentTree[i+treeSize] = initial[i];
        }
        zoneTree = vector<int>(2*treeSize, 1);
        for (int i=treeSize;i>0;i--){
            segmentTree[i-1] = binomialOperator(segmentTree[2*i-2], segmentTree[2*i-1]);
            zoneTree[i-1] = zoneTree[2*i-2] + zoneTree[2*i-1];
        }
        updateTree = vector<UPDATE_TYPE>(2*treeSize, upd_id);
    }
    ~LazySegmentTree(){}
    // get lower zone tied with [l, r)
    vector<int> getLower(int l, int r){
        vector<int> F1(0, 0);
        vector<int> F2(0, 0);
        int L = l + treeSize, R = r + treeSize;
        int lsize = 0, rsize = 0;
        while (L < R){
            if (L & 1){
                F1.push_back(L);
                L++;
                lsize++;
            }
            if (R & 1){
                F2.push_back(R-1);
                R--;
                rsize++;
            }
            L >>= 1;
            R >>= 1;
        }
        vector<int> F(0, 0);
        for (int i=0;i<lsize;i++){
            F.push_back(F1[i]);
        }
        for (int i=rsize;i>0;i--){
            F.push_back(F2[i-1]);
        }
        return F;
    }
    // get path tied with [l, r)
    vector<int> getPath(int l, int r){
        vector<int> F(0, 0);
        set<int> Q;
        int L = l + treeSize, R = r + treeSize;
        while (!(L&1)){
            L >>= 1;
        }
        while (!(R&1)){
            R >>= 1;
        }
        int tmp = L >> 1;
        while (tmp){
            F.push_back(tmp);
            Q.insert(tmp);
            tmp >>= 1;
        }
        tmp = R >> 1;
        while (tmp){
            if (Q.count(tmp)){
                break;
            }
            F.push_back(tmp);
            Q.insert(tmp);
            tmp >>= 1;
        }
        sort(F.begin(), F.end());
        return F;
    }
    // update interval [l, r) by x
    void update(int l, int r, UPDATE_TYPE x){
        vector<int> lazyF = getLower(l, r);
        vector<int> pathF = getPath(l, r);
        int ps = pathF.size(), ls = lazyF.size();
        int k, tmp;
        for (int i=0;i<ps;i++){
            k = pathF[i];
            tmp = k << 1;
            updateTree[tmp] = lazyOperator(updateTree[k], updateTree[tmp]);
            updateTree[tmp^1] = lazyOperator(updateTree[k], updateTree[tmp^1]);
            segmentTree[k] = updateOperator(zoneTree[k], updateTree[k], segmentTree[k]);
            updateTree[k] = updateId;
        }
        for (int i=0;i<ls;i++){
            k = lazyF[i];
            updateTree[k] = lazyOperator(x, updateTree[k]);
        }
        MONOID_TYPE left, right;
        for (int i=ps;i>0;i--){
            k = pathF[i-1];
            tmp = k << 1;
            left = updateOperator(zoneTree[tmp], updateTree[tmp], segmentTree[tmp]);
            right = updateOperator(zoneTree[tmp], updateTree[tmp^1], segmentTree[tmp^1]);
            segmentTree[k] = binomialOperator(left, right);
        }
    }
    MONOID_TYPE getSegmentValue(int l, int r){
        vector<int> lazyF = getLower(l, r);
        vector<int> pathF = getPath(l, r);
        MONOID_TYPE S = monoId;
        int ps = pathF.size(), ls = lazyF.size();
        int k, tmp;
        for (int i=0;i<ps;i++){
            k = pathF[i];
            tmp = k << 1;
            updateTree[tmp] = lazyOperator(updateTree[k], updateTree[tmp]);
            updateTree[tmp^1] = lazyOperator(updateTree[k], updateTree[tmp^1]);
            segmentTree[k] = updateOperator(zoneTree[k], updateTree[k], segmentTree[k]);
            updateTree[k] = updateId;
        }
        for (int i=0;i<ls;i++){
            k = lazyF[i];
            S = binomialOperator(S, updateOperator(zoneTree[k], updateTree[k], segmentTree[k]));
        }
        return S;
    }
};

// x^n mod m
template <int m, typename T>
modint<m> modpow(modint<m> x, T n){
    modint<m> y = 1;
    modint<m> b = x;
    while (n){
        if (n & 1){
            y *= b;
        }
        b *= b;
        n >>= 1;
    }
    return y;
}

// Formal Power Series
template <const int MOD>
class FormalPowerSeries{
    using T = modint<MOD>;
    using FPS = FormalPowerSeries<MOD>;
    protected:
    T primitiveRoot = 0;
    T inverseRoot = 0;
    vector<T> primitiveBaseList;
    vector<T> inverseBaseList;
    int baseSize = 0;
    template <typename Tp>
    static T modpow(T a, Tp b){
        T y = 1;
        T bas = a;
        while (b){
            if (b & 1){
                y *= bas;
            }
            bas *= bas;
            b >>= 1;
        }
        return y;
    }
    void initialize(){
        int q = MOD-1;
        while (!(q & 1)){
            q >>= 1;
            baseSize++;
        }
        primitiveBaseList.resize(baseSize+1);
        inverseBaseList.resize(baseSize+1);
        while (true){
            bool flg = true;
            T tmp = inverseRoot;
            for (int i=0;i<baseSize;i++){
                if (tmp == 1){
                    flg = false;
                    break;
                }
                tmp *= tmp;
            }
            if (tmp != 1){
                flg = false;
            }
            if (flg){
                break;
            }
            inverseRoot++;
        }
        primitiveRoot = modpow(inverseRoot, MOD-2);
        primitiveBaseList[baseSize] = primitiveRoot;
        inverseBaseList[baseSize] = inverseRoot;
        for (int i=baseSize;i>0;i--){
            primitiveBaseList[i-1] = primitiveBaseList[i] * primitiveBaseList[i];
            inverseBaseList[i-1] = inverseBaseList[i] * inverseBaseList[i];
        }
    }
    vector<T> getNumberTheoremTransform(vector<T> v){
        int n = v.size();
        int m = 1;
        int depth = 0;
        while (m < n){
            depth++;
            m <<= 1;
        }
        int pos = 0, tmp = 0, d;
        vector<T> res(m, 0);
        for (int i=0;i<m-1;i++){
            if (pos < n){
                res[i] = v[pos];
            }
            tmp = (i+1)&-(i+1);
            d = 0;
            while (tmp){
                d++;
                tmp >>= 1;
            }
            pos ^= n - (1<<(depth - d));
        }
        if (pos < n){
            res[m-1] = v[pos];
        }
        T grow = 1, seed = 0;
        int offset = 0;
        T x, y;
        for (int i=0;i<depth;i++){
            grow = 1;
            seed = primitiveBaseList[i+1];
            offset = 1 << i;
            for (int k=0;k<offset;k++){
                for (int j=k;j<n;j+=1<<(i+1)){
                    x = res[j];
                    y = res[j+offset] * grow;
                    res[j] = x + y;
                    res[j+offset] = x - y;
                }
                grow *= seed;
            }
        }
        return res;
    }
    vector<T> getInverseNumberTheoremTransform(vector<T> v){
        int n = v.size();
        int m = 1;
        int depth = 0;
        while (m < n){
            depth++;
            m <<= 1;
        }
        int pos = 0, tmp = 0, d;
        vector<T> res(m, 0);
        for (int i=0;i<m-1;i++){
            if (pos < n){
                res[i] = v[pos];
            }
            tmp = (i+1)&-(i+1);
            d = 0;
            while (tmp){
                d++;
                tmp >>= 1;
            }
            pos ^= n - (1<<(depth - d));
        }
        if (pos < n){
            res[m-1] = v[pos];
        }
        T grow = 1, seed = 0;
        int offset = 0;
        T x, y;
        for (int i=0;i<depth;i++){
            grow = 1;
            seed = inverseBaseList[i+1];
            offset = 1 << i;
            for (int k=0;k<offset;k++){
                for (int j=k;j<n;j+=1<<(i+1)){
                    x = res[j];
                    y = res[j+offset] * grow;
                    res[j] = x + y;
                    res[j+offset] = x - y;
                }
                grow *= seed;
            }
        }
        T mInverse = modpow(m, MOD-2);
        for (int i=0;i<m;i++){
            res[i] *= mInverse;
        }
        return res;
    }
    vector<T> convolve(vector<T> f, vector<T> g){
        int n = f.size(), m = g.size();
        int b = 1;
        while (b < n + m){
            b <<= 1;
        }
        vector<T> x(b, 0);
        vector<T> y(b, 0);
        for (int i=0;i<n;i++){
            x[i] = f[i];
        }
        for (int i=0;i<m;i++){
            y[i] = g[i];
        }
        x = getNumberTheoremTransform(x);
        y = getNumberTheoremTransform(y);
        for (int i=0;i<b;i++){
            x[i] *= y[i];
        }
        x = getInverseNumberTheoremTransform(x);
        return x;
    }
    void regularize(){
        int t = degree;
        for (int i=degree;i>=0;i--){
            if (f[i] != 0){
                break;
            }
            t--;
        }
        if (t == degree){
            return;
        }
        f.resize(t+1);
        degree = t;
    }
    public:
    int degree = -1;
    vector<T> f;
    FormalPowerSeries(){
        initialize();
    }
    FormalPowerSeries(const int n){
        degree = 0;
        f.resize(1);
        f[0] = (T)n;
        initialize();
        regularize();
    }
    FormalPowerSeries(const T n){
        degree = 0;
        f.resize(1);
        f[0] = n;
        initialize();
        regularize();
    }
    FormalPowerSeries(const vector<T> g){
        degree = g.size() - 1;
        f.resize(degree+1);
        for (int i=0;i<=degree;i++){
            f[i] = g[i];
        }
        initialize();
        regularize();
    }
    FormalPowerSeries(const FPS& fps){
        degree = fps.degree;
        f.resize(degree+1);
        for (int i=0;i<=degree;i++){
            f[i] = fps.f[i];
        }
        initialize();
        regularize();
    }
    void resize(const int size){
        degree = size - 1;
        f.resize(size);
    }
    void addManually(const int i, const T a){
        if (degree < i){
            degree = i;
            f.resize(degree+1);
        }
        f[i] = a;
    }
    T operator[](const int i) const{
        return f[i];
    }
    friend ostream& operator<<(ostream& os, FPS& fps){
        if (fps.degree == -1){
            os << 0;
            return os;
        }
        for (int i=0;i<=fps.degree;i++){
            if (i){
                os << " ";
            }
            os << fps.f[i];
        }
        return os;
    }
    friend istream& operator>>(istream& ist, FPS& fps){
        for (int i=0;i<=fps.degree;i++){
            T s;
            ist >> s;
            fps.f[i] = s;
        }
        return ist;
    }
    bool operator==(const FPS& fps) const{
        int lt = degree, rt = fps.degree;
        for (int i=degree;i>=0;i--){
            if (f[i] != (T)0){
                break;
            }
            lt--;
        }
        for (int i=fps.degree;i>=0;i--){
            if (fps.f[i] != (T)0){
                break;
            }
            rt--;
        }
        if (lt != rt){
            return false;
        }
        for (int i=0;i<=lt;i++){
            if (f[i] != (T)fps.f[i]){
                return false;
            }
        }
        return true;
    }
    bool operator!=(const FPS& fps) const{
        return !(*this == fps);
    }
    FPS operator+() const{
        FPS res = FPS(*this);
        return res;
    }
    FPS operator-() const{
        FPS res = FPS(*this);
        for (int i=0;i<=degree;i++){
            res.f[i] = -res.f[i];
        }
        return res;
    }
    FPS& operator++(){
        if (degree == -1){
            degree = 0;
            f.resize(1);
            f[0] = 0;
        }
        f[0]++;
        regularize();
        return *this;
    }
    FPS& operator--(){
        if (degree == -1){
            degree = 0;
            f.resize(1);
            f[0] = 0;
        }
        f[0]--;
        regularize();
        return *this;
    }
    FPS operator++(int){
        FPS res = FPS(*this);
        ++*this;
        return res;
    }
    FPS operator--(int){
        FPS res = FPS(*this);
        --*this;
        return res;
    }
    FPS& operator+=(const FPS& rhs){
        if (degree < rhs.degree){
            f.resize(rhs.degree+1);
            degree = rhs.degree;
        }
        for (int i=0;i<=rhs.degree;i++){
            f[i] += rhs.f[i];
        }
        regularize();
        return *this;
    }
    FPS& operator+=(const int rhs){
        if (degree == -1){
            f.resize(1);
            degree = 0;
            f[0] = 0;
        }
        f[0] += (T)rhs;
        regularize();
        return *this;
    }
    FPS& operator+=(const T rhs){
        if (degree == -1){
            f.resize(1);
            degree = 0;
            f[0] = 0;
        }
        f[0] += rhs;
        regularize();
        return *this;
    }
    FPS& operator-=(const FPS& rhs){
        if (degree < rhs.degree){
            f.resize(rhs.degree+1);
            degree = rhs.degree;
        }
        for (int i=0;i<=rhs.degree;i++){
            f[i] -= rhs.f[i];
        }
        regularize();
        return *this;
    }
    FPS& operator-=(const int rhs){
        if (degree == -1){
            f.resize(1);
            degree = 0;
            f[0] = 0;
        }
        f[0] -= (T)rhs;
        regularize();
        return *this;
    }
    FPS& operator-=(const T rhs){
        if (degree == -1){
            f.resize(1);
            degree = 0;
            f[0] = 0;
        }
        f[0] -= rhs;
        regularize();
        return *this;
    }
    FPS& operator*=(const FPS& rhs){
        degree += rhs.degree;
        if (degree == -2){
            degree = -1;
        }
        f = convolve(f, rhs.f);
        regularize();
        return *this;
    }
    FPS& operator*=(const int rhs){
        for (int i=0;i<=degree;i++){
            f[i] *= (T)rhs;
        }
        regularize();
        return *this;
    }
    FPS& operator*=(const T rhs){
        for (int i=0;i<=degree;i++){
            f[i] *= rhs;
        }
        regularize();
        return *this;
    }
    template <typename Tp>
    FPS operator+(const Tp rhs){
        FPS res = FPS(*this);
        res += rhs;
        return res;
    }
    template <typename Tp>
    FPS operator-(const Tp rhs){
        FPS res = FPS(*this);
        res -= rhs;
        return res;
    }
    template <typename Tp>
    FPS operator*(const Tp rhs){
        FPS res = FPS(*this);
        res *= rhs;
        return res;
    }
    FPS& shiftMultiply(int n){
        degree += n;
        f.resize(degree+1);
        for (int i=degree;i>=n;i--){
            f[i] = f[i-n];
        }
        for (int i=0;i<n;i++){
            f[i] = 0;
        }
        return *this;
    }
    FPS shiftMultiplied(int n){
        FPS res = FPS(*this);
        res.shiftMultiply(n);
        return res;
    }
    FPS& shiftDivide(int n){
        if (degree < n){
            degree = -1;
            f.resize(0);
            return *this;
        }
        for (int i=n;i<=degree;i++){
            f[i-n] = f[i];
        }
        degree -= n;
        f.resize(degree+1);
        return *this;
    }
    FPS shiftDivided(int n){
        FPS res = FPS(*this);
        res.shiftDivide(n);
        return res;
    }
    FPS& differentiate(){
        if (degree == -1){
            return *this;
        }
        for (int i=0;i<degree;i++){
            f[i] = (T)(i+1)*f[i+1];
        }
        f[degree] = 0;
        regularize();
        return *this;
    }
    FPS differentiated(){
        FPS res = FPS(*this);
        res.differentiate();
        return res;
    }
    FPS& integrate(){
        degree++;
        f.push_back(0);
        vector<T> invl(degree+1, 1);
        for (int i=2;i<=degree;i++){
            invl[i] = (T)(MOD / i) * (-invl[MOD%i]);
        }
        for (int i=degree;i>0;i--){
            f[i] = invl[i] * f[i-1];
        }
        f[0] = 0;
        regularize();
        return *this;
    }
    FPS integrated(){
        FPS res = FPS(*this);
        res.integrate();
        return res;
    }
    FPS& inverse(){
        T a = 0;
        if (degree == -1){
            return *this;
        }
        a = modpow(f[0], MOD-2);
        FPS res(a), twoFactor(2), tmp(f[0]), h;
        int b = 1;
        while (b < degree+1){
            for (int i=0;i<b;i++){
                if (i+b >= degree+1){
                    break;
                }
                tmp.f.push_back(f[i+b]);
                tmp.degree++;
            }
            b <<= 1;
            h = twoFactor - tmp * res;
            h.resize(b);
            res *= h;
        }
        res.resize(degree+1);
        *this = res;
        return *this;
    }
    FPS inversed(){
        FPS res = FPS(*this);
        res.inverse();
        return res;
    }
    FPS& convertTaylorShift(T c){
        if (degree <= 0){
            return *this;
        }
        vector<T> fact(degree+1, 1);
        vector<T> invf(degree+1, 1);
        for (int i=0;i<degree;i++){
            fact[i+1] = (fact[i] * (i+1));
        }
        invf[degree] = fact[degree].inv();
        for (int i=degree;i>0;i--){
            invf[i-1] = (invf[i] * i);
        }
        vector<T> x(degree+1, 0);
        vector<T> y(degree+1, 0);
        T t = 1;
        for (int i=0;i<=degree;i++){
            x[degree-i] = f[i];
            x[degree-i] *= fact[i];
            y[i] = t;
            y[i] *= invf[i];
            t *= c;
        }
        x = convolve(x, y);
        for (int i=0;i<=degree;i++){
            f[i] = invf[i];
            f[i] *= x[degree-i];
        }
        return *this;
    }
    FPS getTaylorShift(T c){
        FPS res = FPS(*this);
        res.convertTaylorShift(c);
        return res;
    }
};

template <int m>
FormalPowerSeries<m> logarithm(FormalPowerSeries<m> fps){
    int n = fps.degree+1;
    FormalPowerSeries<m> res;
    res = fps.inversed() * fps.differentiated();
    res.resize(n);
    res.integrate();
    return res;
}

template <int m>
FormalPowerSeries<m> exponential(FormalPowerSeries<m> fps, int degreeRequire=-1){
    int resSize = 1;
    if (fps.degree < degreeRequire){
        fps.resize(degreeRequire+1);
    }
    int x = fps.degree + 1;
    if (x & (x - 1)){
        while (x & (x - 1)){
            x &= x - 1;
        }
        x <<= 1;
    }
    FormalPowerSeries<m> res(1), g(1), h(2), q, r, tmp;
    if (fps.degree == -1){
        tmp = 0;
    } else{
        tmp = fps.f[0];
    }
    while (2*resSize <= x){
        q = res * g;
        q.resize(resSize);
        g *= (h - q);
        g.resize(resSize);
        q = tmp.differentiated();
        r = q + g * (res.differentiated() - res * q);
        r.resize(2*resSize);
        r.integrate();
        for (int i=0;i<resSize;i++){
            if (i+resSize > fps.degree){
                break;
            }
            tmp.addManually(i+resSize, fps.f[i+resSize]);
        }
        q = tmp - r;
        q *= res;
        res += q;
        resSize <<= 1;
    }
    return res;
}

template <int m, typename T>
FormalPowerSeries<m> power(FormalPowerSeries<m> fps, T x){
    int t = -1;
    int d = fps.degree;
    modint<m> a = 0;
    for (int i=0;i<=d;i++){
        if (fps.f[i] != 0){
            t = i;
            a = fps.f[i];
            break;
        }
    }
    if (t == -1){
        if (x == 0){
            return FormalPowerSeries<m>(1);
        }
        return FormalPowerSeries<m>(0);
    }
    FormalPowerSeries<m> gps = fps;
    if ((long long) x * t > d){
        return FormalPowerSeries<m>(0);
    }
    gps.shiftDivide(t);
    for (int i=0;i<=gps.degree;i++){
        gps.f[i] /= a;
    }
    gps.resize(d+1);
    gps = logarithm(gps);
    for (int i=0;i<=d;i++){
        gps.f[i] *= x;
    }
    gps = exponential(gps, d);
    gps.shiftMultiply(x * t);
    gps.resize(d+1);
    a = modpow<m>(a, (long long)x);
    for (int i=0;i<=d;i++){
        gps.f[i] *= a;
    }
    return gps;
}

template <int m>
vector<modint<m>> getBernoulliTable(int n){
    vector<modint<m>> fact(n+2, 1);
    for (int i=1;i<=n+1;i++){
        fact[i] = (fact[i-1] * i);
    }
    vector<modint<m>> invf(n+1, 1);
    invf[n] = fact[n+1].inv();
    for (int i=n;i>0;i--){
        invf[i-1] = invf[i] * (i+1);
    }
    FormalPowerSeries<m> r(invf);
    r.inverse();
    vector<modint<m>> res(n+1);
    for (int i=0;i<=n;i++){
        res[i] = r.f[i];
        res[i] *= fact[i];
    }
    return res;
}

template <int m>
vector<modint<m>> getPartitionNumberTable(int n){
    vector<modint<m>> invn(n+1, 1);
    for (int i=2;i<=n;i++){
        invn[i] = -(invn[m%i]);
        invn[i] *= (m / i);
    }
    FormalPowerSeries<m> f;
    f.resize(n+1);
    for (int i=1;i<=n;i++){
        for (int j=1;i*j<=n;j++){
            f.f[i*j] += invn[j];
        }
    }
    f = exponential<m>(f);
    for (int i=0;i<=n;i++){
        invn[i] = f.f[i];
    }
    return invn;
}

template <typename T>
void printVector(vector<T> a){
    int n = (int)a.size();
    if (n <= 0){
        return;
    }
    cout << a[0];
    for (int i=1;i<n;i++){
        cout << " " << a[i];
    }
    cout << endl;
}
