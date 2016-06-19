from HP import *
import random
import math

__author__ = 'Victor Gueorguiev - GRGVIC001'


class GA:
    # size of tournament for parental selection: F. Custodio et. al. has this set to 4
    tournament_size = 4
    # initial probabilities for the operators uniformly distributed
    prob = [1 / 6.0, 1 / 6.0, 1 / 6.0, 1 / 6.0, 1 / 6.0, 1 / 6.0]
    # progress maintained at least at 1.0 to ensure that if no progress exists, operators return to uniformly distributed
    progress = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    # number of calls to an operator's method - many calls wth no progress are penalised after some threshold defined below
    calls = [0, 0, 0, 0, 0, 0]

    # these two lists keep track of which method created the child and parental structures. Values range from 0 to 5
    # for the six operators; a -1 means no operator created it (initial population member)
    parental = [-1, -1]
    child = [-1]

    # function_calls = 0
    # energy_value = -34

    def __init__(self, pop_size, monomer_len, number_h):

        # initial empty population with pop. size and monomer length and number of H monomers initialized
        self.population = []
        self.population_size = pop_size
        self.monomer_length = monomer_len
        self.number_of_h = number_h # this is used in computing the path later on

        # initialize the HP members conformations randomly
        for i in range(self.population_size):
            member = HP(self.monomer_length)
            member.generate_hp_direction(self.monomer_length)
            member.my_operators = list(self.child)
            self.population.append(member)

            ''' uncomment this to check the representations of the initial population '''
            '''print("member "+str(i+1)+" born")
            print(member.get_representation())
            print(member.get_fitness())'''

    def run(self, generations):
        random.seed()
        for i in range(generations):
            print(i) # print current generation

            cumm = [] # cumulative probability distribution for the operators to be used in conditional statements
            for k in range(1, 7):
                cumm.append(sum(self.prob[0:k]))

            for j in range(10):
                # breed new children
                rn = random.random() # get some random number
                # crossover operators below
                if rn < cumm[0]:
                    children = self.__multi_px() # MPX crossover operator
                    self.child[0] = 0
                    self.calls[0] += 1
                elif rn < cumm[1]:
                    children = self.__two_px() # 2PX crossover operator
                    self.child[0] = 1
                    self.calls[1] += 1
                # mutation operators below
                elif rn < cumm[2]:
                    children = self.__mutation_selection()
                    self.__lp(children[0]) # LP mutation operator
                    self.child[0] = 2
                    self.calls[2] += 1
                elif rn < cumm[3]:
                    children = self.__mutation_selection()
                    self.__lm(children[0]) # LM mutation operator
                    self.child[0] = 3
                    self.calls[3] += 1
                elif rn < cumm[4]:
                    children = self.__mutation_selection()
                    self.__smut(children[0]) # SMUT mutation operator
                    self.child[0] = 4
                    self.calls[4] += 1
                else:
                    children = self.__mutation_selection()
                    self.__emut(children[0]) # EMUT mutation operator
                    self.child[0] = 5
                    self.calls[5] += 4

                # comment below is to measure function calls
                ''' if self.child[0] == 5:
                        self.function_calls += 4
                    else:
                        self.function_calls += 1'''
                # add new children to the set of all molecules
                ''' uncomment below to use either the DME insertion or SGA 'normal' insertion methods '''
                # self.__insert(children)
                self.__DME_insert(children)
                self.__check_calls() # check against the calls list to see if any have exceeded threshold and penalize operator if necessary
            self.adjust_prob() # every 10 members created adjust the probabilities as F. Custodio does in his paper

            # comment below is for testing number of function alls required to reach a specified energy value
            ''' best = self.__find_max()
                if best.get_fitness() >= self.energy_value:
                    print(self.function_calls)
                    break'''
        # comment code below is to measure number of monomers of a certain energy in final population
        ''' result = []
            counter = 0
            for i in range(self.energy_value + 10): # surplus energy added for graph completeness
                for member in self.population:
                    if member.get_fitness == i:
                        counter += 1
                result.append(counter)
                counter = 0
            print(result)'''
        return self.__find_max() # finally return the maximum member in the population

    # use the tournament to retrieve two parents for crossover
    def __crossover_selection(self):
        return self.__tournament(2)

    # use the tournament method to retrieve one parent for mutation
    def __mutation_selection(self):
        return self.__tournament(1)

    # select parents of certain specified number to generate new generation
    def __tournament(self, no_candidates):
        result = [] # end result
        candidates = [] # candidates for tournament from our population
        i = 0
        while i != self.tournament_size:
            candidate = self.population[random.randint(0, len(self.population) - 1)] # pick a random member of the pop.
            if not (candidate in candidates):
                candidates.append(candidate) # make sure new candidate has not yet been selected
                i += 1

        candidates = sorted(candidates, key=lambda hp: hp.get_fitness()) # sort the candidates by order of decreasing fitness

        for j in range(no_candidates):
            # choose the best candidate very often, else choose a random one from the remainder
            if random.random() < 0.9:
                parent = candidates[0]
                candidates.remove(candidates[0])
            else:
                candidate = random.randint(1, self.tournament_size - 2)
                parent = candidates[candidate]
                candidates.remove(candidates[candidate])
            self.parental[j] = parent.my_operators[0] # set the parental operators for the current child using the parent's specified operators

            parent = parent.get_representation() # we only use the parents sequence of directions
            result.append(list(parent))

        return result

    def __two_px(self):
        parents = self.__crossover_selection() # get a pair of parents
        cross_point = random.randint(1, self.monomer_length - 1) # choose a random crossover point
        # generate two children as per the usual way of crossing over two genomes
        child1 = list(parents[0][: cross_point] + parents[1][cross_point:])
        child2 = list(parents[1][: cross_point] + parents[0][cross_point:])
        return [child1, child2]

    def __multi_px(self):
        cut_points = [] # the set of crossover points
        c = int(0.1 * self.monomer_length) # the number of crossover points is a function of the monomer length
        for i in range(c):
            # ensure that the crossover points come after one another for list indexing purposes
            cut_points.append(random.randint(1 + i, self.monomer_length))
        cut_points.sort() # ensure that the crossover points come after one another for list indexing purposes
        child1 = []
        child2 = []
        parents = self.__crossover_selection() # get a set of parents
        prior_point = 0

        for cut_point in cut_points:
            chromosome = random.randint(0, 1) # randomly assign which parent donates which segment to a given child
            # give the segment to one child and the other to the second child
            child1 += list(parents[chromosome][prior_point: cut_point])
            child2 += list(parents[1 - chromosome][prior_point: cut_point])
            prior_point = cut_point
        # add final tail bit of the chromosome
        chromosome = random.randint(0, 1)
        child1 += list(parents[chromosome][prior_point:])
        child2 += list(parents[1 - chromosome][prior_point:])
        return [child1, child2]

    # swaps two random bases anywhere on the chain
    def __lp(self, molecule):
        base1 = random.randint(0, self.monomer_length - 1) # choose a random base to swap
        base2 = random.randint(0, self.monomer_length - 1) # choose the base to swap with

        # perform the two way swap with a temporary variable
        temp_direction = molecule[base1]
        molecule[base1] = molecule[base2]
        molecule[base2] = temp_direction

        return molecule

    # swaps two consecutive bases
    def __lm(self, molecule):
        base1 = random.randint(0, self.monomer_length - 1) # choose a random base to perform a swap with
        base2 = (base1 + 1) % self.monomer_length # get the next base in the chain to swap, and if it is the end piece swap with the head of the list

        # perform two way swap
        temp_direction = molecule[base2]
        molecule[base2] = molecule[base1]
        molecule[base1] = temp_direction
        return molecule

    # randomly permutes a segment starting at some random base. segment is size 2 to 7 permutations
    def __smut(self, molecule):
        base1 = random.randint(1, self.monomer_length - 1 - 7) # choose a random base to start at
        segment_len = random.randint(2, 7) # size of segment to randomly permute
        for i in range(segment_len):
            new_dir = random.randint(0, 5)
            molecule[base1 + i] = new_dir # select a new direction and assign it
        return molecule

    # randomly change the direction of a random base until a better energy is found, otherwise leave as is
    def __emut(self, molecule):
        # generate a new HP molecule with just that one random direction changed and compare it with the current molecule
        base1 = random.randint(1, self.monomer_length - 1)
        maximum = HP(self.monomer_length)
        maximum.assign_hp(list(molecule))
        maximum.calc_fitness()
        for i in range(5):
            temp = HP(self.monomer_length)
            molecule[base1] = HP.directions.get(i)
            temp.assign_hp(list(molecule))
            if HP.check_collisions(molecule[1:]):
                temp.calc_fitness()
                if temp.get_fitness() < maximum.get_fitness():
                    maximum.assign_hp(list(molecule))
                    maximum.calc_fitness()
        molecule = maximum.get_representation()
        return molecule

    # search the population for the least fit individual
    def __find_min(self):
        fitness = -10000
        worst = HP(self.monomer_length)
        for member in self.population:
            if member.get_fitness() > fitness:
                fitness = member.get_fitness()
                worst = member
        return worst

    # as above but finds the best
    def __find_max(self):
        fitness = 0
        best = HP(self.monomer_length)
        for member in self.population:
            if member.get_fitness() < fitness:
                fitness = member.get_fitness()
                best = member
        return best

    def __insert(self, children):
        for child in children:
            child_hp = HP(self.monomer_length)
            child_hp.assign_hp(list(child))
            child_hp.my_operators = list(self.child)
            child_hp.parent_operators = list(self.parental)
            min_member = self.__find_min()
            if child_hp.get_fitness() < min_member.get_fitness():
                self.adjust_progress(child_hp)
                self.population.remove(min_member)
                self.population.append(child_hp)

        # reset the operators that created the child to default
        self.child = [-1]
        self.parental = [-1, -1]

    def __DME_insert(self, children):
        for child in children:
            child_hp = HP(self.monomer_length)
            child_hp.assign_hp(list(child))
            child_hp.my_operators = list(self.child)
            child_hp.parent_operators = list(self.parental)

            min_member = self.__find_max()

            if child_hp.get_fitness() < min_member.get_fitness():
                self.adjust_progress(child_hp)

            old_member = self.__find_min_DME(child_hp)
            if child_hp.get_fitness() < old_member.get_fitness():
                self.population.remove(old_member)
                self.population.append(child_hp)
            elif child_hp.get_fitness() == old_member.get_fitness():
                if rn.random() < 0.5:
                    self.population.remove(old_member)
                    self.population.append(child_hp)

        # reset the operators that created the child to default
        self.child = [-1]
        self.parental = [-1, -1]

    def __find_min_DME(self, child):
        # get a list of all candidates sorted by their distance, using DME, from the current child
        candidates = sorted(self.population, key=lambda hp: hp.DME(child))
        return candidates[0] # return the closest individual

    def adjust_progress(self, child):
        max_member = self.__find_max()

        if child.get_fitness() < max_member.get_fitness():
            hh_difference = math.fabs(child.get_fitness() - max_member.get_fitness())
            child_operator = child.my_operators
            parent_operators = child.parent_operators

            if child_operator[0] == 5:
                self.progress[child_operator[0]] += 0.25 * hh_difference
            else:
                self.progress[child_operator[0]] += hh_difference

            self.calls[child_operator[0]] = 0

            if parent_operators != [-1, -1]:
                if parent_operators[0] == 5 and parent_operators[0] != -1:
                    self.progress[parent_operators[0]] += 0.25 * 0.25 * hh_difference # give parents only  a quarter
                elif parent_operators[0] != -1:
                    self.progress[parent_operators[0]] += 0.25 * hh_difference
                self.calls[parent_operators[0]] = 0 # set the calls to the operator back to zero
                if parent_operators[1] == 5 and parent_operators[1] != -1:
                    self.progress[parent_operators[1]] += 0.25 * 0.25 * hh_difference
                elif parent_operators[1] != -1:
                    self.progress[parent_operators[1]] += 0.25 * hh_difference
                self.calls[parent_operators[1]] = 0 # set back to zero

    # penalises a operators progress if there are too many calls and no progress is made
    def __check_calls(self):
        for i in range(len(self.calls)):
            if self.calls[i] > 1000:
                if self.progress[i] >= 2.0:
                    self.progress[i] -= 1.0
                else:
                    self.progress[i] = 1.0
                self.calls[i] = 0

    def adjust_prob(self):
        # probabilities for operator application cannot fall below this threshold
        min_threshold = 0.05

        # calculate the total progress for all operators for scaling purposes
        total = self.progress[0] + self.progress[1] + self.progress[2] + self.progress[3] + self.progress[4] + self.progress[5]

        # compute the new probability of an operator as the fraction of total progress
        if total > 0:
            for i in range(len(self.progress)):
                self.prob[i] = self.progress[i] / total

        # check that no operator prob. is below threshold; if one exists, subtract the difference from the leading
        # probability to ensure that they still sum to one.
        for i in range(len(self.prob)):
            if self.prob[i] < min_threshold:
                difference = min_threshold - self.prob[i]
                self.prob[i] += difference
                self.prob[self.prob.index(max(self.prob))] -= difference

        ''' uncomment this to track how the probabilities evolve after each update '''
        # print(self.prob, self.progress)
