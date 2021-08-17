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

template <typename TYPE>
class BinaryIndexedTree{
    private:
    TYPE binaryOperator(TYPE lhs, TYPE rhs){
        return lhs + rhs;
    }
    TYPE invert(TYPE lhs){
        return -lhs;
    }
    public:
    int size = 0;
    TYPE identity = 0;
    vector<TYPE> tree;
    BinaryIndexedTree(int n){
        size = n;
        tree.resize(size+1);
    }
    BinaryIndexedTree(vector<TYPE> v){
        size = (int)v.size();
        tree.resize(size+1);
        for (int i=0;i<size;i++){
            int p = i + 1;
            while (p <= size){
                tree[p] = binaryOperator(tree[p], v[i]);
                p += p - (p & (p - 1));
            }
        }
    }
    void updateOneElement(int pos, TYPE val){
        int x = pos + 1;
        while (x <= size){
            tree[x] = binaryOperator(tree[x], val);
            x += x - (x & (x - 1));
        }
    }
    // get [l, r) sum
    TYPE getSegment(int l, int r){
        TYPE resl = identity, resr = identity;
        int ll = l;
        int rr = r;
        while (rr){
            resr = binaryOperator(resr, tree[rr]);
            rr &= (rr - 1);
        }
        while (ll){
            resl = binaryOperator(resl, tree[ll]);
            ll &= (ll - 1);
        }
        return binaryOperator(invert(resl), resr);
    }
    TYPE getOneElement(int p){
        return getSegment(p, p+1);
    }
};
