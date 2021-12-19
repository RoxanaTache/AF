#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <stack>
#include <tuple>
#include <queue>
#include <list>
 
using namespace std;
 
class Graph {
    struct nodeStruct {
        int node1, node2, cost;
        bool operator()(nodeStruct const& n1, nodeStruct const& n2) { return n1.cost > n2.cost; }
    };
    vector<list<nodeStruct>> adjacent;
    vector<list<nodeStruct>> transposed();
    void dfs(int &current, vector<bool>&visited, stack<int> &order);
    void bfs(int &start, vector<int> &costs, vector<bool> &visited, deque<int> &queue, int &diameter);
    void _biconnected(int &node, int parent, vector<bool> &visited, vector<int> &level, vector<int> &minLevel, stack<int> &s, vector<vector<int>> &components, vector<pair<int,int>> &criticalEdges);
    void _hardConnected(int &node, vector<bool> &visited, vector<vector<int>> &components, vector<list<nodeStruct>> &t);
    int findParent(int node, vector<int> &parent);
    void _union(vector<int> &parents, vector<int> &height, const int &parent1, const int &parent2);
    bool _hasPath(int &current, int &target, vector<bool> &visited);
public:
    Graph(vector<tuple<int, int, int>> &data, int &nrNodes, bool oriented);
    Graph(int &nrNodes) { adjacent.resize(nrNodes); }
    friend ostream& operator<< (ostream& os, Graph graph) {
        os << graph.adjacent.size() << " nodes\n";
        for(int i = 0; i < graph.adjacent.size(); i++) {
            os << "node " << i + 1 << ": ";
            for(nodeStruct j: graph.adjacent[i])
                os << "(" << j.node2 + 1 << ", " << j.cost << ") ";
            os << "\n";
        }
        return os;
    }
    int connected();
    pair<vector<int>, vector<bool>> costs(int start);
    stack<int> topologicalSort();
    pair<vector<vector<int>>, vector<pair<int,int>>> biconnected();
    vector<vector<int>> hardConnected();
    vector<int> dijkstra(int start);
    pair<vector<int>, bool> bellmanFord(int start);
    pair<vector<nodeStruct>, int> minimumTreeKruskall();
    bool havelHakimi(deque<int> degrees);
    int diameter();
    void insertEdge(int node1, int node2, int cost, bool oriented = true);
    bool hasPath(int current, int target);
    vector<vector<int>> costMatrix();
    vector<vector<int>> royfloyd();
    int hamilton();
    int getCost(int src, int dest) {
        for(auto i: adjacent[src])
            if(i.node2 == dest)
                return i.cost;
    }
};
 
Graph :: Graph(vector<tuple<int, int, int>> &data, int &nrNodes, bool oriented) {
    adjacent.resize(nrNodes);
    for(auto[node1, node2, cost]: data) {
        adjacent[node1].push_back(nodeStruct({node1, node2, cost}));
        if(!oriented)
            adjacent[node2].push_back(nodeStruct({node2, node1, cost}));
    }
}
vector<list<Graph :: nodeStruct>> Graph :: transposed() {
    vector<list<nodeStruct>> ret(adjacent.size());
    for(int i = 0; i < adjacent.size(); i++)
        for(nodeStruct node: adjacent[i])
            ret[node.node2].push_back(nodeStruct({i, node.cost}));
    return ret;
}
void Graph :: dfs(int &current, vector<bool>&visited, stack<int> &order) {
    visited[current] = true;
    for(auto i: adjacent[current])
        if(!visited[i.node2])
            dfs(i.node2, visited, order);
    order.push(current);
}
void Graph :: bfs(int &start, vector<int> &costs, vector<bool> &visited, deque<int> &queue, int &diameter) {
    queue.push_back(start);
    visited[start] = true;
    costs[start] = 0;
    while(queue.empty() != 1) {
        const int current = queue.front();
        for(auto i: adjacent[current])
            if(!visited[i.node2]) {
                costs[i.node2] += (costs[current] + 1);
                diameter = costs[i.node2];
                visited[i.node2] = true;
                start = i.node2;
                queue.push_back(i.node2);
            }
        queue.pop_front();
    }
}
void Graph :: _biconnected(int &node, int parent, vector<bool> &visited, vector<int> &level, vector<int> &minLevel, stack<int> &s, vector<vector<int>> &components, vector<pair<int,int>> &criticalEdges) {
    visited[node] = true;
    minLevel[node] = level[node] = level[parent] + 1;
    s.push(node);
    for (auto x: adjacent[node])
        if (x.node2 != parent) {
            if (visited[x.node2])
                minLevel[node] = min(minLevel[node], level[x.node2]);
            else {
                _biconnected(x.node2, node, visited, level, minLevel, s, components, criticalEdges);
                minLevel[node] = min(minLevel[node], minLevel[x.node2]);
                if(level[node] < minLevel[x.node2])
                    criticalEdges.emplace_back(node, x.node2);
                if (minLevel[x.node2] >= level[node]) {
                    components.resize(components.size()+1);
                    while (s.top() != x.node2) {
                        components[components.size()-1].push_back(s.top());
                        s.pop();
                    }
                    components[components.size()-1].push_back(x.node2);
                    s.pop();
                    components[components.size()-1].push_back(node);
                }
            }
        }
}
void Graph :: _hardConnected(int &node, vector<bool> &visited, vector<vector<int>> &components, vector<list<nodeStruct>> &t) {
    visited[node] = false;
    components[components.size() - 1].push_back(node);
    for(auto j: t[node])
        if(visited[j.node2])
            _hardConnected(j.node2, visited, components, t);
}
int Graph :: findParent(int node, vector<int> &parent) {
    while(parent[node] != 0)
        node = parent[node];
    return node;
}
void Graph :: _union(vector<int> &parents, vector<int> &height, const int &parent1, const int &parent2){
    if(height[parent1] > height[parent2])
        parents[parent2] = parent1;
    else {
        parents[parent1] = parent2;
        if(height[parent1] == height[parent2])
            height[parent2]++;
    }
}
bool Graph :: _hasPath(int &current, int &target, vector<bool> &visited) {
    if(visited[current])
        return false;
    visited[current] = true;
    if(current == target)
        return true;
    for(auto i: adjacent[current])
        if(_hasPath(i.node2, target, visited))
            return true;
    return false;
}
int Graph :: connected() {
    vector<bool> visited(adjacent.size());
    int nr = 0;
    for(int i = 0; i < adjacent.size(); i++)
        if(!visited[i]) {
            stack<int> _;
            nr++;
            dfs(i, visited, _);
        }
    return nr;
}
pair<vector<int>, vector<bool>> Graph :: costs(int start) {
    vector<int> costs(adjacent.size());
    vector<bool> visited(adjacent.size());
    deque<int> queue;
    int _;
    bfs(start, costs, visited, queue, _);
    return make_pair(costs, visited);
}
stack<int> Graph :: topologicalSort() {
    vector<bool> visited(adjacent.size());
    stack<int> order;
    for(int i = 0; i < adjacent.size(); i++)
        if(!visited[i])
            dfs(i, visited, order);
    return order;
}
pair<vector<vector<int>>, vector<pair<int,int>>> Graph :: biconnected(){
    stack<int> s;
    vector<int> level(adjacent.size()), minLevel(adjacent.size());
    vector<bool> visited(adjacent.size());
    vector<vector<int>> components;
    vector<pair<int,int>> criticalEdges;
    for (int i = 0; i < adjacent.size(); i++)
        if (visited[i] == 0)
            _biconnected(i, 0, visited, level, minLevel, s, components, criticalEdges);
    return make_pair(components, criticalEdges);
}
vector<vector<int>> Graph :: hardConnected() {
    stack<int> s;
    vector<bool> visited(adjacent.size());
    vector<vector<int>> components;
    vector<list<nodeStruct>> t = transposed();
    for(int i = 0; i < adjacent.size(); i++)
        if(!visited[i])
            dfs(i, visited, s);
    while(!s.empty()){
        if(visited[s.top()]){
            components.resize(components.size() + 1);
            _hardConnected(s.top(), visited, components, t);
        }
        s.pop();
    }
    return components;
}
vector<int> Graph :: dijkstra(int start) {
    vector<int> visited(adjacent.size()), distance(adjacent.size(), -1);
    priority_queue<nodeStruct, vector<nodeStruct>, nodeStruct> costs;
    costs.push({0, start, 0});
    distance[start] = 0;
    while(costs.empty() != 1) {
        int node = costs.top().node2;
        costs.pop();
        if(!visited[node])
            for(auto i: adjacent[node]){
                if(!visited[i.node2])
                    if(distance[i.node2] == -1 || distance[i.node2] > i.cost + distance[node]){
                        distance[i.node2] = i.cost + distance[node];
                        costs.push({0, i.node2, distance[i.node2]});
                    }
            }
        visited[node] = 1;
    }
    return distance;
}
pair<vector<int>, bool> Graph :: bellmanFord(int start) {
    const int inf = 250001;
    vector<int> visited(adjacent.size()), distance(adjacent.size(), inf);
    priority_queue<nodeStruct, vector<nodeStruct>, nodeStruct> costs;
    costs.push({start, 0});
    distance[start] = 0;
    while(costs.empty() != 1) {
        int node = costs.top().node2;
        costs.pop();
        for(auto i: adjacent[node]){
            if(distance[i.node2] == inf || distance[i.node2] > i.cost + distance[node]){
                distance[i.node2] = i.cost + distance[node];
                costs.push({i.node2, distance[i.node2]});
                visited[node]++;
                if(visited[i.node2] >= adjacent.size())
                    return make_pair(distance, 1);
            }
        }
        visited[node]++;
    }
    return make_pair(distance, 0);
}
pair<vector<Graph :: nodeStruct>, int> Graph :: minimumTreeKruskall() {
    vector<int> parents(adjacent.size()), height(adjacent.size());
    vector<nodeStruct> minimumSearchTree, sortedEdges;
    int cost = 0, nrEdges = 0;
    for(int i = 0; i < adjacent.size(); i++)
        for(auto j: adjacent[i])
            sortedEdges.push_back(j);
    sort(sortedEdges.begin(), sortedEdges.end(), nodeStruct());
    for(int i = sortedEdges.size() - 1; 0 <= i; i--){
        const nodeStruct node = sortedEdges[i];
        const int parent1 = findParent(node.node1, parents);
        const int parent2 = findParent(node.node2, parents);
        if(parent1 != parent2) {
            cost += node.cost;
            nrEdges++;
            minimumSearchTree.push_back(node);
            _union(parents, height, parent1, parent2);
            if(nrEdges == adjacent.size() - 1)
                return make_pair(minimumSearchTree, cost);
        }
    }
    return make_pair(minimumSearchTree, cost);
}
int Graph :: diameter() {
    vector<int> costs(adjacent.size());
    vector<bool> visited(adjacent.size());
    deque<int> queue;
    int current = 0, diameter;
    bfs(current, costs, visited, queue, diameter);
    fill(costs.begin(), costs.end(), 0);
    fill(visited.begin(), visited.end(), 0);
    bfs(current, costs, visited, queue, diameter);
    return diameter;
}
void Graph :: insertEdge(int node1, int node2, int cost, bool oriented) {
    adjacent[node1].emplace_back(nodeStruct{node1, node2, cost});
    if(!oriented)
        adjacent[node2].emplace_back(nodeStruct{node2, node1, cost});
}
bool Graph :: hasPath(int current, int target) {
    vector<bool> visited(adjacent.size());
    return _hasPath(current, target, visited);
}
bool Graph :: havelHakimi(deque<int> degrees) {
    int sum = 0;
    for(auto i: degrees) {
        sum += i;
        if(i >= degrees.size())
            return 0;
    }
    while(1) {
        sort(degrees.begin(), degrees.end(),greater<>());
        if(degrees[0] == 0)
            return 1;
        degrees.pop_front();
        for(int i = 0; i < degrees.size(); i++) {
            degrees[i]--;
            if(degrees[i] < 0)
                return 0;
        }
    }
}
vector<vector<int>> Graph::costMatrix() {
    vector<vector<int>> matrix;
    for(int i = 0; i < adjacent.size(); i++) {
        vector<int> line(adjacent.size());
        for(auto edge: adjacent[i])
            line[edge.node2] = edge.cost;
        matrix.push_back(line);
    }
    return matrix;
}
vector<vector<int>> Graph::royfloyd() {
    vector<vector<int>> matrix = costMatrix();
    for(int i = 0; i < adjacent.size(); i++)
        for(int j = 0; j < adjacent.size(); j++)
            for(int k = 0; k < adjacent.size(); k++)
                if (matrix[k][i]                                    // exista drum intre i si k
                    && matrix[i][j]                                 // exista drum intre i si j
                    && (matrix[k][j] > matrix[k][i] + matrix[i][j]  // daca drumul direct intre k si j este mai lung decat drumul ocolitor (prin i)
                        || !matrix[k][j]) && k != j)                //      sau nu exista deja drum
                    matrix[k][j] = matrix[k][i] + matrix[i][j];     // cel mai mic drum intre k si j devine k - i - j
    return matrix;
}
int Graph :: hamilton() {
    vector<int> perm;
    vector<long long> mins;
    for(int i = 0; i < adjacent.size(); i++)
        perm.push_back(i);
    while(next_permutation(perm.begin(), perm.end())) {
        long long s = 0;
        int ok = 1;
        for(int i = 0; i < perm.size() - 1; i++)
            if(!hasPath(perm[i], perm[i + 1])) {
                ok = 0;
                break;
            }
        if(!hasPath(perm[perm.size() - 1], perm[0]))
            ok = 0;
        if(ok) {
            for(int i = 0; i < perm.size() - 1; i++)
                s += getCost(perm[i], perm[i + 1]);
            s += getCost(perm[perm.size() - 1], perm[0]);
            mins.push_back(s);
        }
    }
    long long min;
    if(!mins.empty()) {
        min = mins[0];
        for(auto i: mins)
            if(min > i)
                min = i;
    }
    else
        min = -1;
    return min;
}
 
int main() {

    return 0;
}