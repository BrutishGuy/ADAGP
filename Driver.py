from GA import *
import time

__author__ = 'Victor Gueorguiev'


def main():
    print("HP model training using Genetic Algorithms\nRun commencing...\n...")

    # parameters, remember to set these!!
    population_size = 500
    monomer_length = 150
    number_h = 75
    generations = 3000
    iterations = 1
    
    start = time.time()
    algo = GA(population_size, monomer_length, number_h)
    best_conformation = algo.run(generations)   
    print(best_conformation.get_representation())
    print(best_conformation.get_fitness())   
    for iteration in range(1,iterations):     
        algo = GA(population_size, monomer_length, number_h)
    
        # run the algorithm and return the best monomer which is then written to file
        current_conformation = algo.run(generations)
        print(current_conformation.get_representation())
        print(current_conformation.get_fitness())
        if best_conformation.get_fitness() > current_conformation.get_fitness():
            best_conformation = current_conformation
    
    end = time.time()
    time_elapsed = end-start
    print("Time elapsed for " + str(iterations) + " runs/trials with " + str(generations) + ": " + str(time_elapsed))
    
    f = open('HPMoleculeResult.txt', 'w')
    f.write('Molecule   Direction\n')
    molecule = best_conformation.get_molecule()
    best_conformation = best_conformation.get_representation()
    for i in range(monomer_length):
        if i == 0:
            f.write(molecule[0] + '\n')
        elif i < monomer_length-1:
            f.write(molecule[i] + '\t' + best_conformation[i] + '\n')
        else:
            f.write(molecule[i] + '\t' + best_conformation[i])
    f.close()

if __name__ == "__main__":
    main()
