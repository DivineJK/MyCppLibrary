template <class T>
class Matrix{
    public:
    int h, w;
    vector<vector<T>> a;
    Matrix(int n = 0, int m = 0, T x = 0){
        h = n;
        w = m;
        a = vector<vector<T>>(h, vector<T>(w, x));
    }
    Matrix(vector<vector<T>> x){
        int h = x.size(), w = x[0].size();
        a = vector<vector<T>>(h, vector<T>(w));
        for (int i=0;i<h;i++){
            assert(w == x[i].size());
            for (int j=0;j<w;j++){
                a[i][j] = x[i][j];
            }
        }
    }
    ~Matrix(){}
    vector<T>& operator[](const int i){
        return a[i];
    }
    int size() const{
        return h;
    }
    friend ostream& operator<<(ostream& os, Matrix& rhs){
        for (int i=0;i<rhs.size();i++){
            if (i){
                os << "\n";
            }
            for (int j=0;j<rhs[0].size();j++){
                if (j){
                    os << " ";
                }
                os << rhs[i][j];
            }
        }
        return os;
    }
    friend istream& operator>>(istream& ist, Matrix& rhs){
        for (int i=0;i<rhs.size();i++){
            for (int j=0;j<rhs[0].size();j++){
                T s;
                ist >> s;
                rhs[i][j] = s;
            }
        }
        return (ist);
    }
    bool operator==(Matrix& rhs){
        if (h != rhs.size() || w != rhs[0].size()){
            return false;
        }
        for (int i=0;i<h;i++){
            for (int j=0;j<w;j++){
                if (a[i][j] != rhs[i][j]){
                    return false;
                }
            }
        }
        return true;
    }
    bool operator!=(Matrix& rhs){
        return !(this == rhs);
    }
    Matrix<T>& operator+=(Matrix& rhs){
        assert(h == rhs.size() && w == rhs[0].size());
        for (int i=0;i<h;i++){
            for (int j=0;j<w;j++){
                a[i][j] += rhs[i][j];
            }
        }
        return *this;
    }
    Matrix<T> operator+(Matrix& rhs){
        Matrix<T> res = Matrix(*this);
        res += rhs;
        return res;
    }
    Matrix<T>& operator-=(Matrix& rhs){
        assert(h == rhs.size() && w == rhs[0].size());
        for (int i=0;i<h;i++){
            for (int j=0;j<w;j++){
                a[i][j] -= rhs[i][j];
            }
        }
        return *this;
    }
    Matrix<T> operator-(Matrix& rhs){
        Matrix<T> res = Matrix(*this);
        res -= rhs;
        return res;
    }
    Matrix<T>& operator*=(Matrix& rhs){
        assert(w == rhs.size());
        vector<vector<T>> res(h, vector<T>(rhs[0].size(), 0));
        for (int i=0;i<h;i++){
            for (int j=0;j<rhs[0].size();j++){
                for (int k=0;k<w;k++){
                    res[i][j] += a[i][k] * rhs[k][j];
                }
            }
        }
        a = res;
        return *this;
    }
    Matrix<T>& operator*=(T rhs){
        for (int i=0;i<h;i++){
            for (int j=0;j<w;j++){
                a[i][j] *= rhs;
            }
        }
    }
    template <class TYPE>
    Matrix<T> operator*(TYPE rhs){
        Matrix<T> res = Matrix(*this);
        res *= rhs;
        return res;
    }
};
