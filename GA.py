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
from weighted_levenshtein import lev, osa, dam_lev
from data.data import Data
import pandas
import random
import math
import sys
import numpy as np
########################### GLOBAL VARIABLES ##################################

base_list = ['A', 'C', 'T', 'G']
gap_char = '-'

# input_sequences = ['ACCTG','ACGAG','CCGCG','CCTGT']
# input_sequences = ['CGTAA','TACA'] # Want consensus to be AGC
orchid_data = Data('./data/ls_orchid.fasta')
input_sequences = list(orchid_data.unaligned_sequences.values())
delete_costs = np.ones(128, dtype=np.float64)

for base in base_list:
    # Make deletion costs 0.1 (causes consensus sequences to prioritize being longer).
    delete_costs[ord(base)] = 0.1
longest_sequence_length = len(max(input_sequences))

population = []
fitness = []
best_cost = 100000000000 #very large number
optimal_consensus = ''

########################## GLOBAL CONSTANTS ###################################

POPULATION_SIZE = len(input_sequences)
NUM_SEQUENCES =  len(input_sequences)
GENERATIONS = 0

###########################  HELPER FUNCTIONS ##################################

def calculate_edit_cost(candidate_consensus, sequence):
    # Return the levenshtein_distance between two strings.
    return lev(candidate_consensus, sequence, delete_costs=delete_costs)

def calculate_consensus_score(candidate_consensus, sequences):
    consensus_score = 0
    for i in range(len(sequences)):
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
    for i in range(len(input_sequences)):
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
    for i in range(len(input_sequences)):
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
        father, mother = random.choices(population, fitness, k=2)
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
    if(len(father) > len(mother)):
        # Add stars to beginning of mother.
        mother = '*' * (len(father) - len(mother)) + mother
    elif(len(mother) > len(father)):
        # Add stars to end of father.
        father = father + '*' * (len(mother) - len(father))

    father_end = math.floor(random.uniform(0, len(father)-1))
    #create a new chromosome with the father's random DNA part
    child = father[0:father_end+1].strip('*')
    #complement the new chromosomes' DNA with the mother's random DNA part
    child += mother[(father_end+1):].strip('*')

    return child

###############################################################################

'''
4.MUTATION
    The implementation of mutation is really easy, just swap two random
    generated nucleotides in the sequence to produce the new chromosome.
'''

def mutate(sequence):
    # IF within 1/population_size chance then we mutate.
    if(math.floor(random.uniform(0, 1/(1-POPULATION_SIZE)))):
        chance = 1/(len(sequence))
        new_sequence = ''
        # For every nucleotide there is a 1/sequence_length chance of mutating.
        for c in sequence:
            if(chance != 1 and math.floor(random.uniform(0, 1/(1-chance))) == 1):
                # If we have a mutation we can have substitution (most likely)
                # 80%, insertion at 5% and deletion at 5%
                if(math.floor(random.uniform(0, 1/(1-0.9))) == 1):
                    # Substitution
                    new_sequence += random.choice(base_list)
                elif(math.floor(random.uniform(0, 1/(1-0.5))) == 1):
                    # Insertion
                    # If size is greater than 1.5x the largest, then don't add!
                    if(len(sequence) > 1.5 * longest_sequence_length):
                        new_sequence += c
                        continue
                    new_sequence += (c + random.choice(base_list))
                else:
                    # Deletion
                    continue
            else:
                new_sequence += c

    return new_sequence


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
        reproduce()
        progress( i/GENERATIONS*100 )
    #end algorithm and print the result
    endProgress()
    print('Smaller cost found: {}'.format(best_cost))
    print('Consensus sequence: {}'.format(optimal_consensus))
    if(optimal_consensus in input_sequences):
        print("Optimal was one of the input_sequences...")


if __name__ == "__main__":
    main()
