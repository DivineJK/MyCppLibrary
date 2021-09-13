template <typename T>
class LongStack{
    protected:
    bool isError = false;
    public:
    int maximumSize = 1 << 30;
    int back = -1;
    int valueStackSize = 0;
    vector<T> valueStack;
    vector<int> cumulativeStack = {0};
    LongStack(int aMaximumSize){
        maximumSize = aMaximumSize;
    }
    bool setIsErrorStack(bool v){
        isError = v;
        return v;
    }
    bool getIsErrorStack(){
        return isError;
    }
    int size(){
        return cumulativeStack[back+1];
    }
    bool topSub(){
        setIsErrorStack(cumulativeStack[back+1] == 0);
        return getIsErrorStack();
    }
    T top(){
        return valueStack[back];
    }
    void push(T value, int num){
        if (maximumSize < (long long)num + cumulativeStack[back+1]){
            setIsErrorStack(true);
            return;
        }
        if (back == -1){
            if (valueStackSize){
                valueStack[0] = value;
                cumulativeStack[1] = num;
            } else{
                valueStack.push_back(value);
                cumulativeStack.push_back(num);
                valueStackSize++;
            }
            back = 0;
            return;
        }
        if (value == valueStack[back]){
            cumulativeStack[back+1] += num;
        } else{
            if (back + 1 == valueStackSize){
                back++;
                valueStack.push_back(value);
                cumulativeStack.push_back(cumulativeStack[back]+num);
                valueStackSize++;
            } else{
                back++;
                valueStack[back] = value;
                cumulativeStack[back+1] = cumulativeStack[back]+num;
            }
        }
    }
    void pop(int num){
        if (cumulativeStack[back+1] < num){
            setIsErrorStack(true);
            return;
        }
        int s = cumulativeStack[back+1];
        if (s == num){
            back = -1;
            return;
        }
        int l = 0, r = back+2;
        int d = (back+2)>>1;
        while (r - l > 1){
            if (cumulativeStack[d] <= s-num){
                l = d;
            } else{
                r = d;
            }
            d = (l + r) >> 1;
        }
        if (cumulativeStack[d] == s-num){
            back = d - 1;
        } else{
            back = d;
            cumulativeStack[back+1] = s - num;
        }
        return;
    }
};
