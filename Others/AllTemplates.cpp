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
