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

class Crossover {
    public:
        static Chromosome one_point(const Chromosome& c1, const Chromosome& c2){
            int cut_point = rand()%(c1.size-1);
            vector<int> genes;
            for(int i=0; i<c1.size; i++){
                if(i<=cut_point){
                    genes.push_back(c1[i]);
                }else{
                    genes.push_back(c2[i]);
                }
            }
            return Chromosome(genes);
        }

};

class Mutation {
    public:
        static void typical(Chromosome & c, int max_gen, double start_time){
            double mutation_param = 0.05;
            double p = mutation_param*pow(M_E, -((clock()-start_time)/CLOCKS_PER_SEC)/(max_gen/5));
            // cout << p <<endl;
            for(int i=0; i<c.size; i++){
                double _s = double_rand();
                if(p>_s){
                    c[i]+=1;
                    c[i]%=2;
                }
            }
        }
};

class Selection {
    public:
        static double shareF_i(vector<int>& diffs, double f_i){
            sort(diffs.begin(), diffs.end());
            int max_diff = diffs[diffs.size()-1];
            
            //linear
            double diff_share_sum = 0;
            int diffs_size = diffs.size();
            for(int i=0; i<diffs_size; i++){
                diff_share_sum += 1.0-pow((diffs[i]/max_diff),2);
            }
            return f_i/diff_share_sum;
        }
        static pair<int, int> Roulette(const int k, vector<Chromosome>& c_set, const vector<pair<pair<int, int>, int> >& edge, bool sharing){
            vector<int> fit_origin;
            vector<double> fitness;
            double fit_sum = 0;

            int c_set_size = c_set.size();
            int best = -1, worst = -1;

            for(int i=0; i<c_set_size; i++){
                int curr_fit = c_set[i].fitness(edge);
                fit_origin.push_back(curr_fit);
                if(best==-1 || curr_fit>best){
                    best = curr_fit;
                }
                if(worst==-1 || curr_fit<worst){
                    worst = curr_fit;
                }
            }
            for(int i=0; i<c_set_size; i++){
                double curr_fit = (fit_origin[i]-worst)+(double)(best-worst)/(k-1);
                fitness.push_back(curr_fit);
            }

            if(sharing){
                vector<int> diffs;
                for(int i=0; i<c_set_size; i++){
                    diffs.clear();
                    for(int j=0; j<c_set_size; j++){
                        if(i!=j){
                            diffs.push_back(c_set[i]-c_set[j]);
                        }
                    }
                    fitness[i] = Selection::shareF_i(diffs, fitness[i]); //F_is
                }
            }

            for(int i=0; i<c_set_size; i++){
                fit_sum += fitness[i];
            }

            int selected1 = -1, selected2 = -1;
            if(fit_sum!=0){
                while(selected1==-1 || selected2==-1){
                    double point = double_rand()*fit_sum;
                    double point_sum = 0;
                    for(int i=0; i<c_set_size; i++){
                        point_sum += fitness[i];
                        if(point < point_sum){
                            if(selected1==-1){
                                selected1 = i;
                                break;
                            }else if(selected2==-1 && i!=selected1){
                                selected2 = i;
                                break;
                            }else{
                                break;
                            }
                        }
                    }
                }    
            }
            

            return make_pair(selected1, selected2);
        }
        static pair<int, int> Tournament(const int k, vector<Chromosome>& c_set, const vector<pair<pair<int, int>, int> >& edge){
            assert(log((int)c_set.size())/log(2)>=k);
            //select_param t
            double t = 0.6;
            //select 2^k chromosome
            set<int> selected_idx;
            set<int> res;
            int selected1=-1, selected2=-1;
            while(selected1==-1 || selected2==-1){
                while(selected_idx.size()<(int)pow(2,k)){
                    int choice = rand()%(int)c_set.size();
                    selected_idx.insert(choice);
                }
                while(res.size()!=1){
                    while(!selected_idx.empty()){
                        res.clear();
                        int first = *selected_idx.begin();
                        selected_idx.erase(selected_idx.begin());
                        int second = *selected_idx.begin();
                        selected_idx.erase(selected_idx.begin());
                        if(c_set[first].fitness(edge) < c_set[second].fitness(edge)){
                            swap(first, second);
                        }
                        double p = double_rand();
                        if(t>p){
                            res.insert(first);
                        }else{
                            res.insert(second);
                        }
                    }
                    selected_idx = res;
                }
                assert(res.size()==1);
                if(selected1==-1)
                    selected1 = *res.begin();
                else{
                    if(*res.begin()!=selected1){
                        selected2 = *res.begin();
                    }
                }
                res.clear();
            }
            
            return make_pair(selected1, selected2);
        }
};

class Replace {
    public:
        static void Genitor(vector<Chromosome> & c_set, vector<Chromosome> solutions, const vector<pair<pair<int, int>, int> >& edge){
            int c_set_size = c_set.size();
            int sol_size = solutions.size();
            sort(c_set.begin(), c_set.end(), ChromosomeSorter(edge));
            for(int i=0, j=0; i<sol_size; i++){
                if(find(c_set.begin(), c_set.end(), solutions[i])!=c_set.end()) continue;
                c_set[j+c_set_size-sol_size] = solutions[i];
                j++;
            }
        }
        static void Preselection(vector<Chromosome> & c_set, pair<int, int> parents, Chromosome solution, const vector<pair<pair<int, int>, int> >& edge){
            double r = 0.6;
            double p = double_rand();
            int selection;
            int better = parents.first, worse = parents.second;
            if(c_set[better].fitness(edge) < c_set[worse].fitness(edge)){
                swap(better, worse);
            }
            if(p<r){
                selection = better;
            }else{
                selection = worse;
            }
            c_set[selection] = solution;
        }
        static void Crowd(vector<Chromosome> & c_set, vector<Chromosome>& solutions, const vector<pair<pair<int, int>, int> >& edge){
            set<int>choices;
            int sol_size = solutions.size();
            int c_set_size = c_set.size();
            for(int j=0; j<sol_size; j++){
                int min_diff;
                int min_diff_idx = -1;
                for(int i=0; i<c_set_size; i++){
                    if(find(choices.begin(), choices.end(), i)!=choices.end()) continue;
                    if(min_diff_idx == -1 || c_set[i]-solutions[j] < min_diff){
                        min_diff = c_set[i]-solutions[j];
                        min_diff_idx = i;
                    }
                }
                if(min_diff_idx!=-1 && min_diff!=0 && solutions[j].fitness(edge) > c_set[min_diff_idx].fitness(edge)){
                    c_set[min_diff_idx] = solutions[j];
                    choices.insert(min_diff_idx);   
                }
            }
        }
};



int main(){
    srand((unsigned int)time(NULL));
    start_time = clock();
    int v, e;
    cin >> v>> e;
    vector<pair<pair<int, int>, int> > edge;

    for(int i=0; i<e; i++){
        int v1, v2, w;
        cin >> v1 >> v2 >> w;
        edge.push_back(make_pair(make_pair(v1-1, v2-1), w)); // translate to 0 based index
    }

    //Generate initial solutions
    int population_size = 16;
    vector<Chromosome> c_set; //current population
    for(int i=0; i<population_size; i++){
        Chromosome new_chro = Chromosome(v);
        c_set.push_back(new_chro);
    }

    do{
        double tick = (clock()-start_time)/CLOCKS_PER_SEC;
        
        vector<Chromosome> next_gen;
        int num_next_gen = (population_size/2)*pow(M_E, -tick/40)+(population_size/2);

        
        int init_tour_val = 1;
        int init_tour_k;
        vector<pair<int, int> > parents;
        for(int i=0; i<num_next_gen; i++){
            /* roulette selection */
            int select_pressure = 3;
            pair<int, int> parent = Selection::Roulette(select_pressure, c_set, edge, false);
            
            if(parent.first==-1 || parent.second==-1) break;
            while(c_set[parent.first]-c_set[parent.second]<=2){ // performance key point #1
                parent = Selection::Roulette(select_pressure, c_set, edge, false);
            }
            parents.push_back(parent);
            
            Chromosome crossovered =  Crossover::one_point(c_set[parent.first], c_set[parent.second]);
            Mutation::typical(crossovered, 180, start_time);
            next_gen.push_back(crossovered);
            
        }
        if((int)next_gen.size()==0){
            break;
        }
        
            
        // Remove Duplicated
        set<Chromosome> non_duplicated_sols;
        for(int i=0; i<next_gen.size(); i++){
            non_duplicated_sols.insert(next_gen[i]);
        }
        next_gen.clear();
        for(set<Chromosome>::iterator ci=non_duplicated_sols.begin(); ci!=non_duplicated_sols.end(); ci++){    
            next_gen.push_back(*ci);
        }

        if(tick<45){
            Replace::Crowd(c_set, next_gen, edge);
        }else {
            Replace::Genitor(c_set, next_gen, edge);
        }

        sort(c_set.begin(), c_set.end(), ChromosomeSorter(edge));
        
        int max_fit = c_set[0].fitness(edge);
        
        max_fit = c_set[0].fitness(edge);
    }while((clock()-start_time)/CLOCKS_PER_SEC < 180);

    int max_fit = -1;
    int max_fit_idx = -1;
    for(int i=0; i<c_set.size(); i++){
        if(max_fit==-1 || c_set[i].fitness(edge) > max_fit){
            max_fit = c_set[i].fitness(edge);
            max_fit_idx = i;
        }
    }
    
    for(int i=0; i<c_set[max_fit_idx].size; i++){
        if(c_set[max_fit_idx][i]==0){
            cout << i+1 << " ";
        }
    }
    cout << endl;
}