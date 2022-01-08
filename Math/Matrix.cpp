template <class T>
class Matrix{
    protected:
    int row = 0;
    int column = 0;
    vector<vector<T>> a;
    public:
    static Matrix<T> getIdentityMatrix(int size){
        Matrix<T> ret(size, size);
        for (int i=0;i<size;i++){
            ret[i][i] = (T)1;
        }
        return ret;
    }
    Matrix(int n = 0, int m = 0, T x = 0){
        row = n;
        column = m;
        a = vector<vector<T>>(row, vector<T>(column));
        for (int i=0;i<row;i++){
            for (int j=0;j<column;j++){
                a[i][j] = x;
            }
        }
    }
    Matrix(vector<vector<T>> x){
        int row = x.size(), column = x[0].size();
        a = vector<vector<T>>(row, vector<T>(column));
        for (int i=0;i<row;i++){
            assert(column == x[i].size());
            for (int j=0;j<column;j++){
                a[i][j] = x[i][j];
            }
        }
    }
    Matrix(T x);
    ~Matrix(){}
    vector<T>& operator[](const int i){
        return a[i];
    }
    int getRow() const{
        return row;
    }
    void setRow(int aRow){
        if (row == aRow){
            return;
        }
        a.resize(aRow);
        row = aRow;
    }
    int getColumn() const{
        return column;
    }
    void setColumn(int aColumn){
        if (column == aColumn){
            return;
        }
        for (int i=0;i<row;i++){
            a[i].resize(aColumn);
        }
        column = aColumn;
    }
    friend ostream& operator<<(ostream& os, Matrix& rhs){
        for (int i=0;i<rhs.getRow();i++){
            if (i){
                os << "\n";
            }
            for (int j=0;j<rhs.getColumn();j++){
                if (j){
                    os << " ";
                }
                os << rhs[i][j];
            }
        }
        return os;
    }
    friend istream& operator>>(istream& ist, Matrix& rhs){
        for (int i=0;i<rhs.getRow();i++){
            for (int j=0;j<rhs.getColumn();j++){
                T s;
                ist >> s;
                rhs[i][j] = s;
            }
        }
        return (ist);
    }
    bool operator==(Matrix& rhs){
        if (row != rhs.getRow() || column != rhs.getColumn()){
            return false;
        }
        for (int i=0;i<row;i++){
            for (int j=0;j<column;j++){
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
        assert(row == rhs.size() && column == rhs[0].size());
        for (int i=0;i<row;i++){
            for (int j=0;j<column;j++){
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
        assert(row == rhs.size() && column == rhs[0].size());
        for (int i=0;i<row;i++){
            for (int j=0;j<column;j++){
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
        assert(column == rhs.getRow());
        vector<vector<T>> res(row, vector<T>(rhs.getColumn(), 0));
        for (int i=0;i<row;i++){
            for (int j=0;j<rhs.getColumn();j++){
                for (int k=0;k<column;k++){
                    res[i][j] += a[i][k] * rhs[k][j];
                }
            }
        }
        a = res;
        return *this;
    }
    Matrix<T>& operator*=(T rhs){
        for (int i=0;i<row;i++){
            for (int j=0;j<column;j++){
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

template <class T, typename Tp>
Matrix<T> pow(Matrix<T> a, Tp n){
    assert(a.getRow() == a.getColumn());
    int m = a.getRow();
    Matrix<T> y = Matrix<T>::getIdentityMatrix(m);
    Matrix<T> b = a;
    while (n){
        if (n & 1){
            y *= b;
        }
        b *= b;
        n >>= 1;
    }
    return y;
}
