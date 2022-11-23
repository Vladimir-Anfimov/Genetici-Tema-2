#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include "MathFunctions.h"
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include <numeric>
#include <iomanip>

std::random_device rd;
std::mt19937 e2(rd());
std::uniform_real_distribution<double> dist(0, 1);

double unif_random() {
    return (double)rand() / ((double)RAND_MAX + 1);
}

using std::vector;

int population_size = 0;
int dimensions = 0;
int dimension_length = 0;
int N = 0;
double translate_space = 0;
double a = 0;
double b = 0;

int precision = 0;
int math_func = 0;
double minimum_evaluation;
int T;
double mutation_probability, crossover_probability;
double selection_pressure = 1;
const int ELITISM = 25;
double best_solution = FLT_MAX;


struct Cromozom
{
    vector<double> real_values;
    double fitness, evaluation, accumulated;
    vector<bool> bits;

    Cromozom()
    {
        real_values = vector<double>(dimensions);
        bits = vector<bool>(N);
    }
    
    Cromozom(const vector<bool>& _bits)
    {
        real_values = vector<double>(dimensions);
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
        for (int j = 0; j < N; j++)
            random_bits[j] = rand() % 2;
       population[i] = Cromozom(random_bits);
    }
}


void evaluate_fitness()
{
    minimum_evaluation = FLT_MAX;
    double maximum_evaluation = FLT_MIN;
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
    best_solution = minimum_evaluation;

    double total_fitness = 0;
    
    for (int i = 0; i < population.size(); i++)
    {
        population[i].fitness = (maximum_evaluation - population[i].evaluation) / (maximum_evaluation - minimum_evaluation + 0.00001);
        total_fitness += population[i].fitness;
    }
    population[0].accumulated = population[0].fitness/total_fitness;
    for (int i = 1; i < population.size(); i++)
    {
        population[i].accumulated = population[i - 1].accumulated + population[i].fitness/total_fitness;
    }
    population[population.size() - 1].accumulated = 1;

     std::sort(population.begin(), population.end(), [](Cromozom a, Cromozom b) {
        return a.evaluation < b.evaluation;
     });

}

void mutation()
{
    for (int i = 0; i < population_size; i++)
    {
        for (int j = ELITISM; j < dimension_length; j++)
            if (unif_random() <= mutation_probability)
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
    for (int i = ELITISM; i < initial_size; i++)
    {
        if (unif_random() <= crossover_probability)
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
   
    for (int i = 0; i < ELITISM; i++)
        new_population[i] = population[i];
    
    int index_new = ELITISM;
    for (int i = ELITISM; i < population_size; i++)
    {
        double rnd = unif_random();
        for(int j = 0; j < population.size();j++)
            if (rnd < population[j].accumulated)
            {
                new_population[index_new++] = population[i];
                break;
            }
    }
   
    population = new_population;
    if (new_population.size() != population_size)
        throw std::invalid_argument("Noua populatie trebuie sa aiba un nr identic de indivizi.");
}


double genetic_algoritm()
{
    generate_population_bits();
    
    for (int i = 0; i < T; i++)
    {
        minimum_evaluation = FLT_MAX;
        evaluate_fitness();
        selection();
        evaluate_fitness();
        mutation();
        crossover();
    }
        std::cout << "Minimul este " << minimum_evaluation <<"\n";
        std::cout << "Evaluarea minima gasita in toate generatiile: " << best_solution << "\n\n";
   return minimum_evaluation;
}


double init_genetic(int _dimensions, int _math_func, int _precision, int _population_size, int _T, 
    double _crossover_probability, double _mutation_probability = -1)
{
    dimensions = _dimensions;
    math_func = _math_func;
    precision = _precision;
    population_size = _population_size;
    T = _T;
    crossover_probability = _crossover_probability;

    init_interval();
 
    best_solution = FLT_MAX;
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
auto dims = vector<int>{ 30 };

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


void META_GA()
{
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