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
