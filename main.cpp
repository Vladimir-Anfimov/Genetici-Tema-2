#include <iostream>
#include "MathFunctions.h"
#include <vector>
#include <random>

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
        population[i].fitness = pow((maximum_evaluation - population[i].evaluation) / (maximum_evaluation - minimum_evaluation + 0.0000001) + 1, selection_pressure);
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



void genetic_algoritm()
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

}


void init_genetic(int _dimensions, int _math_func, int _precision, int _population_size, int _T, 
    float _crossover_probability)
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
    mutation_probability = 1./N;
    translate_space = (b - a) / (pow(2, dimension_length) - 1);
    population = vector<Cromozom>(population_size, Cromozom());
    minimum_evaluation = FLT_MAX;

    genetic_algoritm();
}

int main()
{
    srand((unsigned)time(NULL));

    std::cout << "Michalewics\n";
    init_genetic(30, MICHALEWICS, 5, 200, 2000, 0.4);

    std::cout << "\nRastrigin\n";
    init_genetic(30, RASTRIGIN, 5, 200, 2000, 0.4);

    std::cout << "\nDeJong\n";
    init_genetic(30, DE_JONG, 5, 200, 2000, 0.4);

    std::cout << "\nSCHWEFEL\n";
    init_genetic(30, SCHWEFEL, 5, 200, 2000, 0.4);






}