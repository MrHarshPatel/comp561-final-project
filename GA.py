'''
Sebastian Ballesteros
This is a genetic algorithm applied to the consensus sequence problem.
'''

'''
Algorithm Definitions:
population - array of candidate consensus sequences
fitness - array of fitness values for each member of the population
best_cost - best numerical cost found so far
optimal_consensus - consensus sequence with the best cost found so far
'''

############################## IMPORTS  ########################################

import pandas
import random
import math
import sys

########################## GLOBAL CONSTANTS ###################################

POPULATION_SIZE = 4
NUM_SEQUENCES = 4
GENERATIONS = 0

########################### GLOBAL VARIABLES ##################################

input_sequences = ['ACCTG','ACGAG','CCGCG','CCTGT']
population = []
fitness = []
best_cost = 100000000000 #very large number
optimal_consensus = ''

###########################  HELPER FUNCTIONS ##################################

def calculate_edit_cost(candidate_consensus, sequence):
    cost = 0
    for i in range(0, len(sequence)):
        # if nucleotides are the same in both sequences
        if candidate_consensus[i] == sequence[i]:
            cost += 0
        # for now let's just assume 1 if they are not the same
        elif candidate_consensus[i] != sequence[i]:
            cost += 1
    return cost

def calculate_consensus_score(candidate_consensus, sequences):
    consensus_score = 0
    for i in range(POPULATION_SIZE):
        # Calculate the edit cost for each sequence in the the population
        consensus_score += calculate_edit_cost(candidate_consensus, sequences[i])
    return consensus_score

def swap(sequence,i,j):
    list_sequence = list(sequence)
    list_sequence[i], list_sequence[j] = list_sequence[j], list_sequence[i]
    return ''.join(list_sequence)

def normalize_fitness(fitness):
    total_sum = 0
    for i in range(POPULATION_SIZE):
        total_sum += fitness[i]
    for i in range(POPULATION_SIZE):
        fitness[i] = (fitness[i] / total_sum)

def assign_fitness():
    global best_cost
    global optimal_consensus
    global fitness
    fitness = []
    for i in range(POPULATION_SIZE):
        #calcualte the edit cost of a particular chromosome
        element_cost = calculate_consensus_score(population[i], input_sequences)
        #keep record of element cost and the optimal consensus while creating population
        if element_cost < best_cost:
            best_cost = element_cost
            optimal_consensus = population[i]
        #since we want to have a greater fitness value for a smaller score
        #we have to modify it, i.e the greater the score the smaller the fitness value
        fitness.append( 1 / float(element_cost))
    normalize_fitness(fitness)

def choose_parent(population, probabilities):
    i = 0
    random_number = random.uniform(0,1)
    #generate an index based on its probabilty(fitness)
    while(random_number > 0):
        random_number = random_number - probabilities[i]
        i += 1
    i -= 1
    return population[i]

##just for progress bar
def startProgress(title):
    global progress_x
    sys.stdout.write(title + ": [" + "-"*40 + "]" + chr(8)*41)
    sys.stdout.flush()
    progress_x = 0

def progress(x):
    global progress_x
    x = int(x * 40 // 100)
    sys.stdout.write("#" * (x - progress_x))
    sys.stdout.flush()
    progress_x = x

def endProgress():
    sys.stdout.write("#" * (40 - progress_x) + "]\n")
    sys.stdout.flush()

###############################################################################

'''
1. INITIAL POPULATION.
    The initial population is the input strings
    For every element in the population, its fitness score its stored in the
    fitness array in the same index
'''

def populate():
    global population
    #initialize the population
    for i in range(NUM_SEQUENCES):
        population.append(input_sequences[i])

################################################################################

'''
2. REPRODUCTION
    We must fit the chromosomes to reproduce offspring according to
    their fitness values. So we are going to have a new population.
    In addition, we will calculate the parents through a function, which
    is going to pick a chromosome based on its fitness value. After producing
    offspring we are going to apply crossover and mutation for each child.
'''

def reproduce():
    global population
    global fitness
    new_population = []
    for i in range(POPULATION_SIZE):
        ##choose the parents according to their fitness value
        father = choose_parent(population,fitness)
        mother = choose_parent(population,fitness)
        child = crossover(father, mother)
        child = mutate(child)
        new_population.append(child)
    population = new_population


###############################################################################

'''
3. CROSSOVER
    Crossover is implemented the following way: father will share he's DNA from
    0 to a random number between 0 and the length of the sequence, mother then
    complement with her remaining part
'''

def crossover(father, mother):
    father_end = math.floor(random.uniform(0, len(input_sequences[0])-1))
    #create a new chromosome with the father's random DNA part
    child = father[0:father_end+1]
    #complement the new chromosomes' DNA with the mother's random DNA part
    for i in range(father_end+1, len(input_sequences[0])):
        child += mother[i]
    return child

###############################################################################

'''
4.MUTATION
    The implementation of mutation is really easy, just swap two random
    generated nucleotides in the sequence to produce the new chromosome.
'''

def mutate(sequence):
    #randomly pick two items to swap
    index_1 = random.randint(0, len(sequence)-1)
    index_2 = random.randint(0, len(sequence)-1)
    return swap(sequence, index_1, index_2)


###############################################################################

"""
* MAIN PROGRAM
    To start the Genetic Algorithm we need to start with an initial POPULATION
    then we repeat the reproduction (along with crossover and mutation) of the
    population according to the number of cycles (which is given)
"""

def main():
    global GENERATIONS
    GENERATIONS = int(input("Number of generations desired (1000 recommended) : "))
    startProgress('Generations')
    #set up the population
    populate()
    #repeat the algorithm until the desired generations have been reached
    for i in range(GENERATIONS):
        assign_fitness()
        reproduce();
        progress( i/GENERATIONS*100 )
    #end algorithm and print the result
    endProgress()
    print('Smaller cost found: {}'.format(best_cost))
    print('Consensus sequence: {}'.format(optimal_consensus))


if __name__ == "__main__":
    main()
