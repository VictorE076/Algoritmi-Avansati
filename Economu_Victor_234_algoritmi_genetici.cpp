//
/// Economu Victor, grupa 234
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <tuple>
#include <cstring>
#include <string>
#include <cmath>
#include <stack>
#include <unordered_set>

using namespace std;

#define MAX_R 32767

ifstream in("intrare.txt");
ofstream out("Evolutie.txt");

////////// !!! Global Input Variables
int population_size;
double Start, End; // capetele intervalului [Start, End]
double a, b, c;
int precision;
double crossover_rate;
double mutation_rate;
int num_generations;
//////////

void read_global_input()
{
    in >> population_size;
    in >> Start >> End;
    in >> a >> b >> c;
    in >> precision;
    in >> crossover_rate;
    in >> mutation_rate;
    in >> num_generations;

    in.close();
}

// Transformam numarul zecimal "index" in binar, avand "num_bits" numar total de biti 
string toBinary(int index, int num_bits)
{
    stack<bool> biti;
    while(num_bits > 0 && index != 0)
    {
        biti.push(index % 2);
        index = index >> 1;
        num_bits--;
    }

    // Restul bitilor de la inceput vor fi "0" daca "num_bits" este mare
    while((num_bits--) > 0) 
    {
        biti.push(0);
    }
    
    // Construim sirul de biti ce reprezinta numarul "index"
    string binary = "";
    while(!biti.empty())
    {
        const int bit = biti.top();
        binary += ('0' + bit);
        biti.pop();
    }

    return binary; 
}

// Functia de fitness f(x) = a*(x^2) + b*x + c
double fitness(const double x)
{
    return (double)(a * x * x + b * x + c);
}

// Generarea Populatiei initiale 
vector<pair<double, string>> generatePopulation(const int num_bits, const double discret)
{
    vector<pair<double, string>> population;
    for(int i = 0; i < population_size; i++)
    {
        double num_unitInterval = (double)rand() / MAX_R; // random double between [0, 1]
        double cromozom = (double)(Start + ((double)num_unitInterval * (End - Start))); // cromozom = random value between [Start, End]
        int index = (int)((cromozom - Start) / discret);
        population.push_back(make_pair(cromozom, toBinary(index, num_bits)));
    }

    return population;
}

/// 
double SumFitness(const vector<pair<double, string>>& population)
{
    double sum = 0;
    for (const pair<double, string>& cromozom : population)
    {
        sum += (double)(fitness(cromozom.first));
    }

    return sum;
}
///

// Probabilitatile de selectie pt. fiecare cromozom
vector<double> SelectionProbabilities(const vector<pair<double, string>>& population)
{
    double sumFitness = SumFitness(population);
    vector<double> selectionProbabilities;
    for (const pair<double, string>& cromozom : population)
    {
        selectionProbabilities.push_back((double)(fitness(cromozom.first) / sumFitness));
    }

    return selectionProbabilities;
}

// Probabilitatile cumulate care dau intervalele pt. selectie
vector<double> intervals(const vector<double>& selectionProb)
{
    vector<double> intervalsProb;
    double interval = 0;
    intervalsProb.push_back(0);
    for(const double prob : selectionProb)
    {
        interval += (double)(prob);
        intervalsProb.push_back(interval);
    }

    return intervalsProb;
}

int BinarySearch_Intervals(const double u, const vector<double>& iProb, const int sizeI, const int left, const int right)
{
    if(left >= sizeI - 1 || left > right)
        return -1;

    if(left == right)
    {
        if(iProb[left] <= u && u < iProb[left + 1])
            return (left + 1);
        else
            return -1;
    }
        
    int m = left + ((right - left) >> 1);

    if(iProb[m] <= u && u < iProb[m + 1])
        return (m + 1);
    else if(iProb[m + 1] <= u)
        return BinarySearch_Intervals(u, iProb, sizeI, m + 1, right);
    else
        return BinarySearch_Intervals(u, iProb, sizeI, left, m);

}

// Evidentierea procesului de selectie
vector<pair<double, int>> SelectionProjection(const vector<double>& intervaleProb)
{
    const int sizeI = intervaleProb.size();
    vector<pair<double, int>> selectionProj;
    for(int i = 0; i < population_size; i++)
    {
        double u = (double)rand() / (MAX_R + 1); // random double between [0, 1)
        int cromozom = BinarySearch_Intervals(u, intervaleProb, sizeI, 0, sizeI - 1);
        selectionProj.push_back(make_pair(u, cromozom));
    }
    
    return selectionProj;
}

// Evidentierea cromozomilor care participa la recombinare (in ordinea din "SelectionProjection(...)")
vector<pair<double, string>> SelectionProjectionPopulation(const vector<pair<double, string>>& population, const vector<pair<double, int>>& selectionProj)
{
    vector<pair<double, string>> selectionProjPopulation;
    for(const pair<double, int>& selectie : selectionProj)
    {
        const int cromozom = selectie.second;
        selectionProjPopulation.push_back(population[cromozom - 1]);
    }

    return selectionProjPopulation;
}

// Selectie de cromozomi
vector<pair<double, string>> ProbCrossover(const vector<pair<double, string>>& selectionProjPopulation)
{
    vector<pair<double, string>> recombinari;
    for(const pair<double, string>& selectie : selectionProjPopulation)
    {
        double u = (double)rand() / MAX_R; // random double between [0, 1]
        recombinari.push_back(make_pair(u, selectie.second));
    }

    return recombinari;
} 

int selectRandomCromozom(vector<int>& cromozomiRecombinare)
{
    const int sz_v = cromozomiRecombinare.size();
    int index = abs(rand()) % sz_v;
    int crom = cromozomiRecombinare[index];
    swap(cromozomiRecombinare[index], cromozomiRecombinare[sz_v - 1]);
    cromozomiRecombinare.pop_back();
    return crom;
}

// Generare perechi + Punctul de rupere generat (in urma combinarii)
void Perechi(const int k, const int num_bits, vector<int>& CromozomiRecombinare, vector<pair<double, string>>& selectionProjPopulation)
{
    // !!! Unul dintre cromozomi nu va fi modificat (daca avem un numar impar de cromozomi care trebuie pusi in perechi de cate 2)
    while(CromozomiRecombinare.size() >= 2) 
    {
        const int crom1 = selectRandomCromozom(CromozomiRecombinare);
        const int crom2 = selectRandomCromozom(CromozomiRecombinare);

        int punctRupere = (int)rand() % num_bits;

        
        if(CromozomiRecombinare.size() != 1)
        {
            if(!k)
            {
                out << "Recombinare dintre cromozomul " << crom1 + 1 << " cu cromozomul " << crom2 + 1 << ":\n";
                out << selectionProjPopulation[crom1].second << " " << selectionProjPopulation[crom2].second << " punct  " << punctRupere << endl;
            }
        
            // Incrucisare de biti (2 cromozomi)
            for(int i = 0; i < punctRupere && i < num_bits; i++)
            {
                char c1 = selectionProjPopulation[crom1].second[i];
                char c2 = selectionProjPopulation[crom2].second[i];
                selectionProjPopulation[crom1].second[i] = c2;
                selectionProjPopulation[crom2].second[i] = c1;
            }

            if(!k)
                out << "Rezultat   " << selectionProjPopulation[crom1].second << " " << selectionProjPopulation[crom2].second << endl;
        }
        else
        {
            const int crom3 = selectRandomCromozom(CromozomiRecombinare);
            
            if(!k)
            {
                out << "Recombinare dintre cromozomul " << crom1 + 1 << " cu cromozomul " << crom2 + 1 << " si cu cromozomul " << crom3 + 1 << ":\n";
                out << selectionProjPopulation[crom1].second << " " << selectionProjPopulation[crom2].second << " " << selectionProjPopulation[crom3].second << " punct  " << punctRupere << endl;
            }
        
            // Incrucisare de biti (3 cromozomi)
            for(int i = 0; i < punctRupere && i < num_bits; i++)
            {
                char c1 = selectionProjPopulation[crom1].second[i];
                char c2 = selectionProjPopulation[crom2].second[i];
                char c3 = selectionProjPopulation[crom3].second[i];
                selectionProjPopulation[crom1].second[i] = c3;
                selectionProjPopulation[crom2].second[i] = c1;
                selectionProjPopulation[crom3].second[i] = c2;
            }

            if(!k)
                out << "Rezultat   " << selectionProjPopulation[crom1].second << " " << selectionProjPopulation[crom2].second << " " << selectionProjPopulation[crom3].second << endl;
        }

    }
}


int toInt(const string& bi)
{
    const int sz = bi.size();
    int pw = 1;
    int index = 0;
    for(int i = sz - 1; i >= 0; i--)
    {
        index += (bi[i] - '0') * pw;
        pw = pw << 1;
    }

    return index;
}

// Cromozomii rezultati dupa recombinare
void CromozomiRez(const double discret, vector<pair<double, string>>& selectionProjPopulation)
{
    const int sz_s = selectionProjPopulation.size();
    for(int i = 0; i < sz_s; i++)
    {
        string binar = selectionProjPopulation[i].second;
        selectionProjPopulation[i].first = (double)(Start + toInt(binar) * discret);
    }
}

// Populatie dupa mutatia aleatoare
vector<int> MutatiiRand_Cromozomi(const double discret, vector<pair<double, string>>& selectionProjPopulation)
{
    vector<int> cromozomi_cu_mutatii;
    const int sz_s = selectionProjPopulation.size();
    for(int i = 0; i < sz_s; i++)
    {
        const int lenBinar = selectionProjPopulation[i].second.size();
        bool flag = true;
        for(int j = 0; j < lenBinar; j++)
        {
            double u = (double)rand() / MAX_R; // random double between [0, 1]
            if(u < mutation_rate)
            {
                if(flag)
                {
                    cromozomi_cu_mutatii.push_back(i + 1);
                    flag = false;
                }
                
                if(selectionProjPopulation[i].second[j] == '0')
                    selectionProjPopulation[i].second[j] = '1';
                else
                    selectionProjPopulation[i].second[j] = '0';
            }
        }

        string binar = selectionProjPopulation[i].second;
        selectionProjPopulation[i].first = (double)(Start + toInt(binar) * discret);
    }

    return cromozomi_cu_mutatii;
}

// Valoarea maxima si cea medie a fitness-ului populatiei
pair<double, double> Max_Medie_Fitness(const vector<pair<double, string>>& selectionProjPopulation)
{
    double max_Fitness = 0;
    double mean_Fitness = 0;

    for(const pair<double, string>& cromozom : selectionProjPopulation)    
    {
        if(fitness(cromozom.first) > max_Fitness)
            max_Fitness = fitness(cromozom.first);

        mean_Fitness += (double)(fitness(cromozom.first));
    }
    mean_Fitness = (double)(mean_Fitness / selectionProjPopulation.size());

    return make_pair(max_Fitness, mean_Fitness);
}

int main() 
{
    srand(time(0)); // initializare seed pt. generarea de numere aleatoare

    read_global_input();
    
    int num_bits = ceil(log2((End - Start) * pow(10, precision)));
    double discret = (double)((End - Start) / pow(2, num_bits));


    ///// Populatia initiala:
    out << "Populatia initiala" << endl;
    vector<pair<double, string>> population = generatePopulation(num_bits, discret);
    for(int i = 0; i < population_size; i++)
    {
        out << (i + 1) << ": " << population[i].second << " x= " << population[i].first << " f= " << fitness(population[i].first) << endl;
    }
    out << endl;
    /////

    for(int k = 0; k < num_generations; k++)
    {
        ///// Probabilitati selectie:
        vector<double> selectionProb = SelectionProbabilities(population);
        
        if(!k)
        {
            out << "Probabilitati selectie" << endl;
            for(int i = 0; i < population_size; i++)
            {
                out << "cromozom  " << (i + 1) << " probabilitate " << selectionProb[i] << endl;
            }
            out << endl;
        }
        /////

        ///// Probabilitatile cumulate care dau intervalele pt. selectie:
        vector<double> intervaleProb = intervals(selectionProb);
        
        if(!k)
        {
            out << "Intervale probabilitati selectie" << endl;
            for(const double interval : intervaleProb)
            {
                out << interval << " ";
            }
            out << endl << endl;
        }
        /////
        
        ///// Procesul de selectie:
        vector<pair<double, int>> selectionProj = SelectionProjection(intervaleProb);
        
        if(!k)
        {
            for(const pair<double, int>& selectie : selectionProj)
            {
                out << "u= " << (double)(selectie.first) << " selectam cromozomul " << selectie.second << endl; 
            }
            out << endl;
        }
        /////

        ///// Evidentierea cromozomilor care participa la recombinare:
        vector<pair<double, string>> selectionProjPopulation = SelectionProjectionPopulation(population, selectionProj);
        
        if(!k)
        {
            out << "Dupa selectie:" << endl;
            for(int i = 0; i < population_size; i++)
            {
                out << (i + 1) << ": " << selectionProjPopulation[i].second << " x= " << selectionProjPopulation[i].first << " f= " << fitness(selectionProjPopulation[i].first) << endl;
            }
            out << endl;
        }
        /////

        ///// Participa la Recombinare - Incrucisare (Probabilitate):
        
        vector<int> cromozomiRecombinare;
        vector<pair<double, string>> recombinari = ProbCrossover(selectionProjPopulation);
        
        if(!k)
            out << "Probabilitatea de incrucisare " << crossover_rate << endl;
            
        for(int i = 0; i < population_size; i++)
        {
            if(!k)
                out << (i + 1) << ": " << recombinari[i].second << " u=" << recombinari[i].first;
            
            if(recombinari[i].first < crossover_rate)
            {
                if(!k)
                    out << "<" << crossover_rate << " participa";
                    
                cromozomiRecombinare.push_back(i);
            }

            if(!k)        
                out << endl;
        }

        if(!k)
            out << endl;
        /////

        ///// Generare perechi + Punctul de rupere generat (in urma combinarii): 
        Perechi(k, num_bits, cromozomiRecombinare, selectionProjPopulation);
        /////

        ///// Cromozomii rezultati dupa recombinare:
        CromozomiRez(discret, selectionProjPopulation);
        
        if(!k)
        {
            out << "\nDupa recombinare:" << endl;
            for(int i = 0; i < population_size; i++)
            {
                out << (i + 1) << ": " << selectionProjPopulation[i].second << " x= " << selectionProjPopulation[i].first << " f= " << fitness(selectionProjPopulation[i].first) << endl;
            }
            out << endl;
        }
        /////

        ///// Populatie dupa mutatia aleatoare:
        vector<int> cromozomi_cu_mutatii = MutatiiRand_Cromozomi(discret, selectionProjPopulation);
        
        if(!k)
        {
            out << "Probabilitate de mutatie pentru fiecare gena " << mutation_rate << endl;
            out << "Au fost modificati cromozomii:\n";
            for(const int cromozomi : cromozomi_cu_mutatii)
            {
                out << cromozomi << endl;
            }
            //
            out << "\nDupa mutatie:" << endl;
            for(int i = 0; i < population_size; i++)
            {
                out << (i + 1) << ": " << selectionProjPopulation[i].second << " x= " << selectionProjPopulation[i].first << " f= " << fitness(selectionProjPopulation[i].first) << endl;
            }
            out << endl;
        }
        ///// 

        ///// Valoarea maxima si cea medie a fitness-ului populatiei:
        if(!k)
            out << "Evolutia maximului\t\tEvolutia mediei\n";
        
        pair<double, double> max_medie_fitness = Max_Medie_Fitness(selectionProjPopulation);
        out << max_medie_fitness.first << "\t\t\t\t\t" << max_medie_fitness.second << endl; 
        /////
        
        population = selectionProjPopulation;
    }
    
    out.close();
    return 0;
}

