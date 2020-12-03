'''
Sebastian Ballesteros 
This is a genetic algorithm used applied on a population of 10 citiies, each with an x and y
coordinate. The goal of the algorithm is to solve the famous travelling salesman problem.
'''

'''
Algorithm Definitions:
population - array of individuals tuples, each with (x,y) coordinate
order - array of randomly generate numbers which indicates  order to visit the cities
fitness - array of fitness values for each member of the population
best_distance - best numerical distance found so far
optimal_order - order for visiting cities which yields the best distance
'''

############################## IMPORTS  ########################################
import pandas
import random
import math
import sys
########################## GLOBAL CONSTANTS ###################################
POPULATION_SIZE = 0
NUM_CITIES = 10
GENERATIONS = 0
########################### GLOBAL VARIABLES ##################################
data_set = pandas.read_csv('cities.csv') #dataset is a 2D array with [coordinate][city]
cities = []
population = []
order = [1,2,3,4,5,6,7,8,9,0]
fitness = []
best_distance = 100000000000 #very large number
optimal_order = []
###########################  HELPER FUNCTIONS ##################################

def calculate_distance(order, cities):
    total_distance = 0
    for x in range(NUM_CITIES-1):
        first_city = cities[order[x]]
        second_city = cities [order[x+1]]
        distance = math.hypot(second_city[0] - first_city[0],  second_city[1] - first_city[1])
        total_distance += distance
    return total_distance

def randomize(order, num):
    for x in range (num):
        #get 2 random numbers (indeces) for swapping
        index_1 = random.randint(0,NUM_CITIES-1)
        index_2 = random.randint(0,NUM_CITIES-1)
        swap(order, index_1, index_2)
    return order

def swap(order,i,j):
    temp = order[i]
    order[i] = order[j]
    order[j] = temp

def normalize_fitness(fitness):
    total_sum = 0
    for i in range(POPULATION_SIZE):
        total_sum += fitness[i]
    for i in range(POPULATION_SIZE):
        fitness[i] = (fitness[i] / total_sum)

def assign_fitness():
    global best_distance
    global optimal_order
    global fitness
    fitness = []
    for i in range(POPULATION_SIZE):
        #calcualte the distance between of a particular chromosome
        element_distance = calculate_distance(population[i], cities)
        #keep record of best distance and the optimal order while creating population
        if element_distance < best_distance:
            best_distance = element_distance
            optimal_order = population[i]
        #since we want to have a greater fitness value for a shorter distance
        #we have to modify it, i.e the longer the distance the smaller the fitness value
        fitness.append( 1 / float(calculate_distance(population[i], cities)))
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
    Read CSV file with the dataset information needed to create the cities
    The population will consist of multiple arrays filled with random sequence
    of numbers i.e [1,2,4,10, etc..] which is the order on which the cities are
    going to be visited in
    For every element in the population, its fitness score its stored in the
    fitness array in the same index
'''

def populate():
    global population
    #store the cities from the dataset
    for i in range(NUM_CITIES):
        cities.append( (data_set['x'][i], data_set['y'][i]) ) #append the tupple x,y
        pass
    #initialize the population
    for i in range(POPULATION_SIZE):
        new_order = randomize(order, 5)
        population.append( new_order.copy() )

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
        mutate(child)
        new_population.append(child)
    population = new_population


###############################################################################

'''
3. CROSSOVER
    Crossover is implemented the following way: first by generating two random numbers
    between 0 and 9 which will be the indeces taken from the father as [start, end],
    later we will perfom crossover by taking the cities remained from the mother as long
    as they don't already exist in the child's chromosome.
'''

def crossover(father, mother):
    father_start = math.floor(random.uniform(0, NUM_CITIES-1))
    father_end = math.floor(random.uniform(father_start+1, NUM_CITIES-1))
    #create a new chromosome with the father's random DNA part
    child = father[father_start:father_end+1]
    #complement the new chromosomes' DNA with the mother's random DNA part
    for i in range(NUM_CITIES):
        #we need to make sure the DNA information is not already in the chromosome
        if (mother[i] not in child):
            child.append(mother[i])
    return child

###############################################################################

'''
4.MUTATION
    The implementation of mutation is really easy, just swap two random
    generated cities to produce the new chromosome.
'''

def mutate(order):
    #randomly pick two items to swap
    index_1 = random.randint(0,NUM_CITIES-1)
    index_2 = random.randint(0,NUM_CITIES-1)
    swap(order, index_1, index_2)

###############################################################################

"""
* MAIN PROGRAM
    To start the Genetic Algorithm we need to start with an initial POPULATION
    then we repeat the reproduction (along with crossover and mutation) of the
    population according to the number of cycles (which is given)
"""

def main():
    global GENERATIONS
    global POPULATION_SIZE
    GENERATIONS = int(input("Number of generations desired (1000 recommended) : "))
    POPULATION_SIZE = int(input("Population size desired (200 recommended) : "))
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
    print('Best distance found: {}'.format(best_distance))
    print('Order in which cities are visited: {}'.format(optimal_order))


if __name__ == "__main__":
    main()