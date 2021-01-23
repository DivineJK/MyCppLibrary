#include <bits/stdc++.h>
using namespace std;
template <typename T>
void test_out_vec(vector<T> a){
    unsigned n = a.size();
    cout << "[";
    for (int i=0;i<n;i++){
        if (i){
            cout << ", ";
        }
        cout << a[i];
    }
    cout << "]" << endl;
}
