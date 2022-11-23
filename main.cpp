#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include "MathFunctions.h"
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include <numeric>

std::random_device rd;
std::mt19937 e2(rd());
std::uniform_real_distribution<float> dist(0, 1);


#define ELITISM 25

using std::vector;

int population_size = 0;
int dimensions = 0;
int dimension_length = 0;
int N = 0;
float translate_space = 0;
float a = 0;
float b = 0;

int precision = 0;
int math_func = 0;
float minimum_evaluation;
int T;
float mutation_probability, crossover_probability;
float selection_pressure = 1;


struct Cromozom
{
    vector<float> real_values;
    float fitness, evaluation, accumulated;
    vector<bool> bits;

    Cromozom()
    {
        real_values = vector<float>(dimensions);
        bits = vector<bool>(N);
    }
    
    Cromozom(const vector<bool>& _bits)
    {
        real_values = vector<float>(dimensions);
        bits = vector<bool>(_bits);

        for (int j = 0; j < dimensions; j++)
		{
		    int x_real = 0;
            for (int k = dimension_length * j; k < dimension_length * (j+1); k++)
            {
                x_real *= 2;
                x_real += bits[k];
            }
            real_values[j] = a + x_real * translate_space;
		}

    }
};

vector<Cromozom> population;


void init_interval() {
    if (math_func == DE_JONG || math_func == RASTRIGIN)
        a = -5.12, b = 5.12;
    else if (math_func == MICHALEWICS)
        a = 0, b = M_PI;
    else
        a = -500, b = 500;

}

void generate_population_bits()
{
    auto random_bits = vector<bool>(N);
    for (int i = 0; i < population_size; i++)
    {
       for(int j=0;j< N;j++)
    		random_bits[j] = dist(e2) > 0.5 ? 1 : 0;
       population[i] = Cromozom(random_bits);
    }
}


void evaluate_fitness()
{
    minimum_evaluation = FLT_MAX;
    float maximum_evaluation = FLT_MIN;
    for (int i = 0; i < population.size(); i++)
    {
        if (math_func == DE_JONG)
            population[i].evaluation = MathFunctions::de_jong(population[i].real_values);
        else if(math_func == RASTRIGIN)
            population[i].evaluation = MathFunctions::Rastrigin(population[i].real_values);
        else if(math_func == MICHALEWICS)
            population[i].evaluation = MathFunctions::Michalewics(population[i].real_values);
        else 
            population[i].evaluation = MathFunctions::Schwefel(population[i].real_values);

        minimum_evaluation = std::min(population[i].evaluation, minimum_evaluation);
        maximum_evaluation = std::max(population[i].evaluation, maximum_evaluation);
    }
    float total_fitness = 0;

    for (int i = 0; i < population.size(); i++)
    {
        population[i].fitness = 1. / pow(population[i].evaluation - minimum_evaluation + 0.000001, 10);
        total_fitness += population[i].fitness;
    }

    population[0].accumulated = population[0].fitness;
    for (int i = 1; i < population.size(); i++)
    {
        population[i].accumulated = population[i - 1].accumulated + population[i].fitness/total_fitness;
    }
}

void mutation()
{
    for (int i = 0; i < population_size; i++)
    {
        for (int j = 0; j < dimension_length; j++)
            if (dist(e2) <= mutation_probability)
            {
                population[i].bits[j] = !population[i].bits[j];
                population.push_back(Cromozom(population[i].bits));
                population[i].bits[j] = !population[i].bits[j];
            }
    }
}

void crossover()
{
    int initial_size = population.size();
    bool found = false;
    int prevIndex = 0;
    for (int i = 0; i < initial_size; i++)
    {
        if (dist(e2) <= crossover_probability)
        {
            if (found == false)
            {
                found = !found;
                prevIndex = i;
            }
            else {
                found = !found;
                int splitOne = rand() % N, splitTwo = rand() % N;
                auto candidateOne = vector<bool>(population[prevIndex].bits);
                auto candidateTwo = vector<bool>(population[i].bits);
               
                if (splitOne > splitTwo)
                    std::swap(splitOne, splitTwo);
                for (int i = splitOne; i <= splitTwo; i++)
                {
                    bool x = candidateOne[i]; 
                    candidateOne[i] = candidateTwo[i];
                    candidateTwo[i] = x;
                }

                population.push_back(Cromozom(candidateOne));
                population.push_back(Cromozom(candidateTwo));
 
            }
        }
    }
}


void selection()
{
    auto new_population = vector<Cromozom>(population_size);
    std::sort(population.begin(), population.end(), [](Cromozom a, Cromozom b) {
        return a.fitness > b.fitness;
        });
    int index_new = ELITISM;
    for (int i = ELITISM; i < population_size; i++)
    {
        float rnd = dist(e2);
        for(int j = 0; j < population.size();j++)
            if (population[j].accumulated <= rnd < population[j+1].accumulated)
            {
                new_population[index_new++] = population[i];
                break;
            }
    }
    for (int i = 0; i < ELITISM; i++)
        new_population[i] = population[i];

    population = new_population;
    if (new_population.size() != population_size)
        throw std::invalid_argument("Noua populatie trebuie sa aiba un nr identic de indivizi.");
}



float genetic_algoritm()
{
    generate_population_bits();
    
    for (int i = 0; i < T; i++)
    {
        minimum_evaluation = FLT_MAX;
        evaluate_fitness();
        selection();
        mutation();
        crossover();

    }
        std::cout << "Minimul este " << minimum_evaluation <<"\n";
   return minimum_evaluation;
}


float init_genetic(int _dimensions, int _math_func, int _precision, int _population_size, int _T, 
    float _crossover_probability, float _mutation_probability = -1)
{
    dimensions = _dimensions;
    math_func = _math_func;
    precision = _precision;
    population_size = _population_size;
    T = _T;
    crossover_probability = _crossover_probability;

    init_interval();
 
    dimension_length = (int)ceil(log2((b - a) * pow(10, precision)));
    N = dimension_length * dimensions;
    if(_mutation_probability == -1)
        mutation_probability = 1./N;

    translate_space = (b - a) / (pow(2, dimension_length) - 1);
    population = vector<Cromozom>(population_size, Cromozom());
    minimum_evaluation = FLT_MAX;

    return genetic_algoritm();
}

#define RUNS 3 
auto dims = vector<int>{ 5, 10, 30 };

void GA()
{
    for (int i = 0; i < dims.size(); i++)
    {
        char file_rastrigin[100];
        char file_michalewics[100];
        char file_dejong[100];
        char file_schewefel[100];

        sprintf(file_rastrigin, "ga_results/rastrigin_dim_%d.txt", dims[i]);
        sprintf(file_michalewics, "ga_results/michalewics_dim_%d.txt", dims[i]);
        sprintf(file_dejong, "ga_results/dejong_dim_%d.txt", dims[i]);
        sprintf(file_schewefel, "ga_results/schwefel_dim_%d.txt", dims[i]);

        std::ofstream rastrigin(file_rastrigin);
        std::ofstream michalewics(file_michalewics);
        std::ofstream dejong(file_dejong);
        std::ofstream schwefel(file_schewefel);


        for (int j = 0; j < RUNS; j++)
        {
            std::cout << "\n------------------------------------------- Dimensiuni "<<dims[i] << "\n";

			std::cout << "\Rastrigin\n";
			rastrigin << init_genetic(dims[i], RASTRIGIN, 5, 200, 2000, 0.4) << "\n";

			std::cout << "\nMichalewics\n";
			michalewics << init_genetic(dims[i], MICHALEWICS, 5, 200, 2000, 0.4) << "\n ";

			std::cout << "\nDeJong\n";
            dejong << init_genetic(dims[i], DE_JONG, 5, 200, 2000, 0.4) << "\n";

			std::cout << "\nSchwefel\n";
			schwefel << init_genetic(dims[i], SCHWEFEL, 5, 200, 2000, 0.8) << "\n";
        }
    }
}


vector<bool> random_cadidate_SA(int _size)
{
    auto candidate = vector<bool>(_size);
    for (int i = 0; i < _size; i++)
        candidate[i] = dist(e2) > 0.5 ? 1 : 0;
    return candidate;
}

vector<float> decode_to_real_SA(vector<bool>& candidate, int LDs, int Ds, double _translate_space, double As)
{
    auto x_real = vector<float>(Ds);
    for (auto jj = 0; jj < Ds; jj++)
    {
        auto x_bits = vector<bool>(LDs);
        copy(candidate.begin() + LDs * jj, candidate.begin() + LDs * (jj + 1), x_bits.begin());
        auto x_int = 0;
        for (auto ii = 0; ii < LDs; ii++)
        {
            x_int *= 2;
            x_int += x_bits[ii];
        }
        x_real[jj] = a + x_int * _translate_space;
    }
    return x_real;
}

vector<float> simulated_annealing(int static_dimenstion_GA, int _math_func_GA)
{
    int As = 0;
    int Bs = 1;
 
    int Ds = 4;
    int PRECISIONs = 2;
    int LDs = (int)ceil(log2((Bs - As) * pow(10, PRECISIONs)));
    int Ns = LDs * Ds;
    double translate_space_SA = (Bs - As) / (pow(2, LDs) - 1);
    float temperature = 100;
    int max_termination = 50;

    auto genetic_func = [LDs, Ds, translate_space_SA, As, static_dimenstion_GA, _math_func_GA](vector<bool>& SA_bits) -> float {
        auto paramsSA = decode_to_real_SA(SA_bits, LDs, Ds, translate_space_SA, As);
/*      params[0] = population_size
        params[1] = T
        params[2] = CP
        params[3] = MP
 */     
        int nr_runs = 1;
        auto results = vector<float>();
        for (int i = 0; i < nr_runs; i++)
        {
            //auto x = init_genetic(static_dimenstion_GA, _math_func_GA, 2,
            //    int(params[0] * 200), int(params[1] * 200), 0.1);
            std::cout << paramsSA[0] << " "<< paramsSA[1] << " " << paramsSA[2]<<" "<<paramsSA[3]<< " ";
            auto x = init_genetic(5, DE_JONG, 5, 200, 2000, 0.4);
           results.push_back(x);
        }
        return std::reduce(results.begin(), results.end()) / nr_runs;
    };
    auto best_bits = random_cadidate_SA(Ns);
    float best_fx = genetic_func(best_bits);

    auto candidate_now_bits = random_cadidate_SA(Ns);
    float candidate_now_fx = genetic_func(candidate_now_bits);

    while (temperature > 0.001)
    {
        for (int ii = 0; ii < max_termination; ii++)
        {
            auto candidate_next_bits = vector<bool>(candidate_now_bits);
            int random_index = rand() % candidate_next_bits.size();
            candidate_next_bits[random_index] = !candidate_next_bits[random_index];

            float candidate_next_fx = genetic_func(candidate_next_bits);

            if (candidate_next_fx < candidate_now_fx)
            {
                candidate_now_bits = vector<bool>(candidate_next_bits);
                candidate_now_fx = candidate_next_fx;
                if (candidate_now_fx < best_fx)
                {
                    best_bits = vector<bool>(candidate_now_bits);
                    best_fx = candidate_now_fx;
                }
            }
            else if (dist(e2) < exp(-abs(candidate_next_fx - candidate_now_fx) / temperature))
            {
                candidate_now_bits = vector<bool>(candidate_next_bits);
                candidate_now_fx = candidate_next_fx;
            }
        }
        std::cout << "Temperature is " << temperature << " and actual best solution is " << best_fx << '\n';
        temperature *= 0.99;
    }
    auto best_optimization = decode_to_real_SA(best_bits, LDs, Ds, translate_space_SA, As);

    std::cout << "Best Fx: " << best_fx << "\n";

    std::cout << "Best pop size: " << best_optimization[0] << "\n";
    std::cout << "Best nr of generatios (T): " << best_optimization[1] << "\n";
    std::cout << "Best crossover prop (CP): " << best_optimization[2] << "\n";
    std::cout << "Best mutation prop (MP): " << best_optimization[3] << "\n";

    return best_optimization;
}


void META_GA()
{
    int dimensions_ga = 5;
    auto meta_ga_result = simulated_annealing(5, MICHALEWICS);
}

int main()
{
    srand((unsigned)time(NULL));

    int which;
    std::cout << "Alege care varianta de algoritm folosesti (0 sau 1)\n";
    std::cin >> which;
    if (which == 0) GA();
    else META_GA();
}