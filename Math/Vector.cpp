template <typename T>
class Vector{
    public:
    int vecSize = 0;
    vector<T> a;
    Vector(){}
    Vector(int n, T initVal = 0){
        vecSize = n;
        a = vector<T>(n, initVal);
    }
    Vector(vector<T> b){
        vecSize = (int)b.size();
        a = vector<T>(vecSize);
        for (int i=0;i<vecSize;i++){
            a[i] = b[i];
        }
    }
    Vector(const Vector<T>& b){
        vecSize = b.size();
        a = vector<T>(vecSize);
        for (int i=0;i<vecSize;i++){
            a[i] = b[i];
        }
    }
    int size() const{
        return vecSize;
    }
    friend ostream& operator<<(ostream& os, Vector& rhs){
        for (int i=0;i<rhs.size();i++){
            if (i){
                os << " ";
            }
            os << rhs[i];
        }
        return os;
    }
    friend istream& operator>>(istream& ist, Vector& rhs){
        for (int i=0;i<rhs.size();i++){
            T s;
            ist >> s;
            rhs[i] = s;
        }
        return (ist);
    }
    T operator[](int i){
        return a[i];
    }
    Vector<T> operator+(){
        Vector<T> res(vecSize);
        for (int i=0;i<vecSize;i++){
            res[i] = a[i];
        }
        return res;
    }
    Vector<T> operator-(){
        Vector<T> res(vecSize);
        for (int i=0;i<vecSize;i++){
            res[i] = -a[i];
        }
        return res;
    }
    Vector<T>& operator+=(Vector<T>& rhs){
        assert(rhs.size() == vecSize);
        for (int i=0;i<vecSize;i++){
            a[i] += rhs[i];
        }
        return *this;
    }
    Vector<T>& operator-=(Vector<T>& rhs){
        assert(rhs.size() == vecSize);
        for (int i=0;i<vecSize;i++){
            a[i] -= rhs[i];
        }
        return *this;
    }
    Vector<T>& operator*=(T x){
        for (int i=0;i<vecSize;i++){
            a[i] *= x;
        }
        return *this;
    }
    Vector<T>& operator/=(T x){
        for (int i=0;i<vecSize;i++){
            a[i] /= x;
        }
        return *this;
    }
};

template <typename T>
int len(Vector<T>& vec){
    return vec.size();
}

template <typename T>
T dot(Vector<T>& vec1, Vector<T>& vec2){
    assert(vec1.size() == vec2.size());
    T S = 0;
    int n = vec1.size();
    for (int i=0;i<n;i++){
        S += vec1[i] * vec2[i];
    }
    return S;
}
