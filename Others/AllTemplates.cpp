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
