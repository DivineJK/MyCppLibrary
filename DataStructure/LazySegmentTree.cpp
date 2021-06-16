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
