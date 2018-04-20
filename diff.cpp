#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <bitset>
#include <cstdlib>
#include <ctime>
#include <string>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <set>

using namespace std;
double start_time;

double double_rand(){
    return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
}


class Chromosome {
    public:
        vector<int> genes;
        int size;
        Chromosome(int size){
            this->size = size;
            for(int i=0; i<size; i++){
                this->genes.push_back(rand()%2);
            } 
        }
        Chromosome(vector<int> genes){
            this->size = genes.size();
            this->genes = genes;
        }
        int fitness (const vector<pair<pair<int, int>, int> >& edge) {
            int e_size = edge.size();
            int ret = 0;
            for(int i=0; i<e_size; i++){
                int v1 = edge[i].first.first;
                int v2 = edge[i].first.second;
                if(this->genes[v1] != this->genes[v2]){
                    ret += edge[i].second;
                }
            }
            return ret;
        }
        void print (){
            for(int i=0; i<this->size; i++){
                cout << this->genes[i];
            }
            cout << endl; 
        }

        friend ostream& operator<<(ostream& os, const Chromosome& c);
        friend int operator -(Chromosome& a, const Chromosome& b);
        bool operator <( const Chromosome &b ) const{
            for(int i=0; i<this->size; i++){
                if((this->genes)[i]<b[i]){
                    return true;
                }
            }
            return false;
        }
        bool operator ==( const Chromosome &b ) const{
            bool equals = true;
            for(int i=0; i<this->size; i++){
                if((this->genes)[i]!=b[i]){
                    equals =false;
                }
            }
            if(equals == true) return true;
            
            equals = true;
            for(int i=0; i<this->size; i++){
                if((this->genes)[i]==b[i]){
                    equals = false;
                }
            }
            return equals;
        }

        const int& operator[] (size_t index) const{
            return this->genes[index];
        }
        int& operator[] (size_t index){
            return this->genes[index];
        }
}; 
ostream& operator<<(ostream& os, const Chromosome& c)  
{  
    for(int i=0; i<c.size; i++){
        os << c.genes[i];
    }
    return os;  
}  
int operator -(Chromosome& a, const Chromosome& b){
    int diff = 0;
    for(int i=0; i<a.size; i++){
        diff += a[i]!=b[i];
    }
    return diff;
}
bool compareChromosome(Chromosome a, Chromosome b, vector<pair<pair<int, int>, int> > edge){
    return a.fitness(edge) > b.fitness(edge);
}
class ChromosomeSorter{
    vector<pair<pair<int, int>, int> > edge;
    public:
        ChromosomeSorter(vector<pair<pair<int, int>, int> > edge){this->edge=edge;}
        bool operator()(Chromosome a, Chromosome b) const{
            return compareChromosome(a, b, this->edge);
        }
};


int main(){
    string s;
    vector<Chromosome> cs;
    for(int i=0; i<2; i++){
        getline(cin, s);
        cout << s << endl;
        vector<int> gene;
        for(int i=0; i<s.size(); i++){
            gene.push_back(s[i]-'0');
        }
        cs.push_back(Chromosome(gene));
    }

    cout << cs[0]-cs[1] << endl;
}