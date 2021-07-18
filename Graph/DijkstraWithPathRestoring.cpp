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
