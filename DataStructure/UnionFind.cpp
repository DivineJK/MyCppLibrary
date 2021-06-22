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

class UnionFind{
    public:
    int n = 0;
    vector<int> par;
    vector<int> size;
    UnionFind(int m){
        n = m;
        par = vector<int>(n, 0);
        size = vector<int>(n, 1);
        for (int i=0;i<n;i++){
            par[i] = i;
        }
    }
    int root(int x){
        if (x == par[x]){
            return x;
        }
        par[x] = root(par[x]);
        return par[x];
    }
    bool same(int x, int y){
        return (root(x) == root(y));
    }
    void unite(int x, int y){
        x = root(x);
        y = root(y);
        if (x == y){
            return;
        }
        if (size[x] < size[y]){
            int tmp = x;
            x = y;
            y = tmp;
        }
        par[y] = x;
        size[x] += size[y];
    }
};
