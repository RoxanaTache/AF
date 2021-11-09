#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <deque>
#include <stack>
using namespace std;

class Graf {
    vector<list<int>> noduri, transpus;
    void dfs(int nod, vector<int> &vizitate) {
        vizitate[nod] = 1;
        for(auto i: noduri[nod])
            if(!vizitate[i])
                dfs(i, vizitate);
    }
    void bfs(int start, vector<int> &costuri, vector<bool> &vizitate,deque<int> &coada) {
        coada.push_back(start);
        vizitate[start] = true;
        costuri[start] = 0;
        while(coada.empty() != 1) {
            int k = coada.front();
            for(auto i: noduri[k])
                if(!vizitate[i]) {
                    costuri[i] += (costuri[k] + 1);
                    vizitate[i] = true;
                    coada.push_back(i);
                }
            coada.pop_front();
        }
    }
    void dfsTopo(int nod, vector<bool> &vizitate, stack<int> &rezultat) {
        vizitate[nod] = true;
        for(auto i: noduri[nod])
            if(!vizitate[i])
                dfsTopo(i, vizitate, rezultat);
        rezultat.push(nod);
    }
    void dfsBiconex(int nod, int tata, int &nrSolutii, vector<bool> &vizitate, vector<int> &nivel, vector<int> &nma, stack<int> &stiva, vector<int> solutie[], vector<pair<int,int>> &muchiiCritice) {
        vizitate[nod] = 1;
        nma[nod] = nivel[nod] = nivel[tata] + 1;
        stiva.push(nod);
        for (auto x: noduri[nod]) {
            if (x != tata){
                if (vizitate[x])
                    nma[nod] = min(nma[nod], nivel[x]);
                else{
                    dfsBiconex(x, nod, nrSolutii, vizitate, nivel, nma, stiva, solutie, muchiiCritice);
                    nma[nod] = min(nma[nod], nma[x]);

                    //muchii critice
                    if(nivel[nod] < nma[x])
                        muchiiCritice.push_back(make_pair(nod, x));

                    //componente biconexe
                    if (nma[x] >= nivel[nod]) {
                        nrSolutii++;
                        while (stiva.top() != x) {
                            solutie[nrSolutii].push_back(stiva.top());
                            stiva.pop();
                        }
                        solutie[nrSolutii].push_back(x);
                        stiva.pop();
                        solutie[nrSolutii].push_back(nod);
                    }
                }
            }
        }
    }
    void dfsConex(int nod, vector<int> &vizitate, stack<int> &s) {
        vizitate[nod] = 1;
        for(auto i: noduri[nod])
            if(!vizitate[i])
                dfsConex(i, vizitate, s);
        s.push(nod);
    }
    void dfsTranspus(int nod, int contor, vector<int> &vizitate, vector<int> *solutie) {
        vizitate[nod] = 2;
        solutie[contor].push_back(nod);
        for(auto j: transpus[nod])
            if(vizitate[j] == 1)
                dfsTranspus(j, contor, vizitate, solutie);
    }

public:
    Graf(vector<list<int>> _noduri, vector<list<int>> _transpuse) : noduri(_noduri), transpus(_transpuse) {}
    friend ostream& operator<< (ostream& os, Graf graf) {
        os << graf.noduri.size() << " nodes\n";
        for(int i = 0; i < graf.noduri.size(); i++) {
            os << "node " << i << ": ";
            for(int j: graf.noduri[i])
                os << j << " ";
            os << "\n";
        }
        return os;
    }
    int nrConexe() {
        vector<int> vizitate(noduri.size());
        int nr = 0;
        for(int i = 0; i < noduri.size(); i++)
            if(!vizitate[i]) {
                nr++;
                dfs(i, vizitate);
            }
        return nr;
    }
    void costuriBFS(int nod, ostream& os = cout) {
        vector<int> costuri(noduri.size());
        deque<int> coada;
        vector<bool> vizitate(noduri.size());
        bfs(nod, costuri, vizitate, coada);
        for (int i = 0; i < noduri.size(); i++) {
            if(vizitate[i])
                os << costuri[i] << " ";
            else os << -1 << " ";
        }
    }
    void sortareTopologica(ostream& os = cout){
        vector<bool> vizitate(noduri.size());
        stack<int> rezultat;
        int nr = 0;
        for(int i = 0; i < noduri.size(); i++)
            if(!vizitate[i])
                dfsTopo(i, vizitate, rezultat);
        while(rezultat.empty() != 1){
            os << rezultat.top() + 1 << " ";
            rezultat.pop();
        }
    }
    void componenteBiconexe(ostream& os = cout){
        vector<bool> vizitate(noduri.size());
        vector<int> nivel(noduri.size()), nma(noduri.size()),  solutie[noduri.size()];
        stack<int> stiva;
        vector<pair<int,int>> muchiiCritice;
        int nrSolutii = 0;
        for (int i = 0; i < noduri.size(); i++)
            if (vizitate[i] == 0)
                dfsBiconex(i, 0, nrSolutii, vizitate, nivel, nma, stiva, solutie, muchiiCritice);

        //Componentele biconexe
//        os << nrSolutii << "\n";
//        for (int i = 1; i <= nrSolutii; i++) {
//            for (auto j: solutie[i])
//                os << j + 1 << " ";
//            os << "\n";
//        }

//        //Muchiile critice
//        for(int i = 0; i < muchiiCritice.size(); i++)
//            cout << "(" << muchiiCritice[i].first + 1 << ","<< muchiiCritice[i].second + 1 << ")" << " ";
    }
    void componenteTareConexe(ostream& os = cout) {
        vector<int> vizitate(noduri.size()), solutie[noduri.size()];
        stack<int> s;
        int contor = 0;
        for(int i = 0; i < noduri.size(); i++)
            if(!vizitate[i])
                dfsConex(i, vizitate, s);
        while(s.empty() != 1){
            int nod = s.top();
            if(vizitate[nod] == 1){
                contor++;
                dfsTranspus(nod, contor, vizitate, solutie);
            }
            s.pop();
        }
        os << contor << "\n";
        for(int i = 1; i <= contor; i++){
            for(auto j: solutie[i])
                os << j + 1 << " ";
            os << "\n";
        }
    }
};

int main() {
    ifstream in("file");
    ofstream out("ctc.out");
    int noduri, muchii, start;
    in >> noduri >> muchii;
//    in >> start; //doar ptr bfs
    vector<list<int>> aux(noduri), aux1(noduri);
    for(int i = 0; i < muchii; i++) {
        int x1, x2;
        in >> x1 >> x2;
        aux[x1 - 1].push_back(x2 - 1);
//        aux[x2 - 1].push_back(x1 - 1);  // daca e neorientat
//        aux1[x2 - 1].push_back(x1 - 1); //ptr graful transpus
    }
    Graf graf(aux, aux1);
//    out << graf.nrConexe(); //ptr DFS
//    graf.costuriBFS(start - 1, out); //ptr BFS
//    graf.sortareTopologica(out);
//    graf.componenteBiconexe(out); // Componente biconexe si muchii critice
//    graf.componenteTareConexe(out);
    return 0;
}