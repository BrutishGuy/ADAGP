import random as rn
import math
import numpy as np
import Data as dat

__author__ = 'Victor Gueorguiev - GRGVIC001'


class HP(object):
    # default values for molecules and directions listed in dictionary data structure for ease of access
    molecules = {0: 'H', 1: 'P'}
    directions = {0: 'U', 1: 'D', 2: 'L', 3: 'R', 4: 'B', 5: 'F'}

    molecule = dat.molecules150[0] # molecule taken from Data.py, it is the actual HP sequence of monomers

    # the operators that created this HP molecule and those that created it's parents
    my_operators = []
    parent_operators = []

    ''' initialize basic variables'''
    def __init__(self, n):
        # initialize basic variables
        self.size = n
        self.representation = [''] * self.size  # represents the shape of the molecule as directions followed by monomers
        self.fitness = 0
        self.path = []  # path stores the Cartesian coordinate path followed by the molecule for later calculations
        self.mags = []  # stores the magnitudes of the coordinates on the path above for ease of calculations

    ''' assigns a shape to a molecule and performs some calculations for later usage'''
    def assign_hp(self, shape):
        self.representation = shape # assign the shape of the molecule in terms of directions to follow
        self.repair()  # repair the shape if necessary
        # calculate paths, magnitudes and fitness for later use
        self.calc_path()
        self.calc_mags()
        self.calc_fitness()

    ''' as with assign_hp() above but generates a sequence of directions to use instead of using assigned/ calculated ones '''
    def generate_hp_direction(self, n):
        shape = [''] * self.size
        rn.seed()

        for i in range(n):
            if self.molecule[i] != self.molecules.get(0):
                self.molecule[i] = self.molecules.get(1)

            direction = rn.randint(0, 5)
            shape[i] = self.directions.get(direction)
        # make first direction empty, this is so that it fits with representation while keeping same length list
        shape[0] = ''

        self.representation = shape
        self.repair()
        self.calc_path()
        self.calc_mags()
        self.calc_fitness()

    ''' generates a random molecule as well as the sequence of directions to follow. This one is not used in the actual algorithm but is included for interests sake '''
    def generate_hp_full(self, h):
        hp = [''] * self.size
        shape = [''] * self.size
        result = []
        rn.seed()
        n = self.size

        for i in range(h):
            location = rn.randint(0, n - 1)
            while hp[location] == self.molecules.get(0):
                location = rn.randint(0, n - 1)
            hp[location] = self.molecules.get(0)

        for i in range(n):
            if hp[i] != self.molecules.get(0):
                hp[i] = self.molecules.get(1)

            direction = rn.randint(0, 5)
            shape[i] = self.directions.get(direction)
        # make first direction empty, this is so that it fits with representation while keeping same length list
        shape[0] = ''
        result.append(hp)
        result.append(shape)
        self.representation = shape
        self.molecule = hp
        print(result)
        # repair mechanism
        self.repair()
        self.calc_path()
        self.calc_mags()
        self.calc_fitness()

    ''' generates the molecules path and stores the cartesian coordinate at each point, checking to see if the molecule
     has traversed that point already'''
    @staticmethod
    def check_collisions(sequence):
        coordinate = [0, 0, 0]
        coordinates = [list(coordinate)]

        # updates the coordinate, adds it to the sequence if not already there, else a collision exists and returns false
        for monomer in sequence:
            if monomer == 'U':
                coordinate[2] += 1
            elif monomer == 'D':
                coordinate[2] += -1
            elif monomer == 'F':
                coordinate[0] += 1
            elif monomer == 'B':
                coordinate[0] += -1
            elif monomer == 'L':
                coordinate[1] += 1
            elif monomer == 'R':
                coordinate[1] += -1
            if coordinate in coordinates:
                return False
            else:
                coordinates.append(list(coordinate))
        return True   # no collisions exist

    ''' return the fitness of the molecule once it is calculated '''
    def get_fitness(self):
        return self.fitness

    ''' return the shape of the molecule '''
    def get_representation(self):
        return self.representation

    ''' return the HP sequence of monomers taken from the data set '''
    def get_molecule(self):
        return self.molecule

    ''' calculate the magnitude of a 3D Cartesian point '''
    @staticmethod
    def __mag(coordinate):
        magnitude = 0
        for coord in coordinate:
            magnitude += coord ** 2.0
        return math.sqrt(magnitude)

    ''' calculate the distance between any two 3D Cartesian points as vectors '''
    @staticmethod
    def __dist(vec1, vec2):
        result = [0, 0, 0]
        for i in range(3):
            result[i] = vec1[i] - vec2[i]
        return result

    def DME(self, other):
        other_path = np.array(other.get_mags())

        this_path = np.array(self.mags)

        distance = (np.subtract(this_path, other_path))
        distance = distance * distance
        distance = np.sum(distance)

        distance = np.sqrt(distance / (self.size * (self.size - 1) / 2.0)) # computes the distance as per formula in the paper using numpy
        return distance

    def repair(self):
        rn.seed()
        done = False
        for i in range(1, self.size):
            tries = 0
            while not self.check_collisions(self.representation[1:i + 1]): # whenever a collision is detected the direction is randomly permuted until there is no collision
                tries += 1
                direction = rn.randint(0, 5)
                self.representation[i] = self.directions.get(direction)
                if tries == 50:
                    self.fitness = 1 # if too many tries are unsuccessful set the fitness to be invalid ( a value of 1 is later taken to mean bad)
                    done = True
                    break
            if done:
                done = False
                break
        if self.representation[0] != '':
            self.representation[0] = ''

    @staticmethod
    def __update_coordinate(monomer, coordinate):
        if monomer == 'U':
            coordinate[2] += 1
        elif monomer == 'D':
            coordinate[2] += -1
        elif monomer == 'F':
            coordinate[0] += 1
        elif monomer == 'B':
            coordinate[0] += -1
        elif monomer == 'L':
            coordinate[1] += 1
        elif monomer == 'R':
            coordinate[1] += -1

    def calc_fitness(self):
        # check if this is a valid conformation; if not, the fitness is zero, the value of one is set by
        # repair mechanism to differentiate between an invalid and default (new) molecule
        if self.fitness == 1:
            self.fitness = 0
            return

        fitness = 0

        coordinate = [0, 0, 0]
        coordinates = list(self.path)
        monomer = self.representation[1]
        next_coordinate = [0, 0, 0]
        self.__update_coordinate(monomer, next_coordinate)
        next_monomer = self.representation[2]
        # check adjacent monomers  if they are H monomers, and ensure that they are not on the path  as those don't count
        for i in range(1, len(self.representation)):
            prev_coordinate = list(coordinate)
            if i > 1:
                self.__update_coordinate(monomer, coordinate)

            monomer = self.representation[i]

            if 0 < i < len(self.representation) - 1:
                next_monomer = self.representation[i + 1]

                if self.molecule[i - 1] == 'H':
                    for j in range(3):
                        test = list(coordinate)
                        test[j] -= 1
                        if test in coordinates and test != next_coordinate and test != prev_coordinate:
                            fitness -= 1
                        test = list(coordinate)
                        test[j] += 1
                        if test in coordinates and test != next_coordinate and test != prev_coordinate:
                            fitness -= 1
                    coordinates.remove(coordinate)
                self.__update_coordinate(next_monomer, next_coordinate)
        self.fitness = fitness

    ''' calculate the path of the H molecules to calculate the fitness'''
    def calc_path(self):
        self.path = []
        coordinate = [0, 0, 0]
        self.path.append(list(coordinate))

        for i in range(1, len(self.representation)):
            self.__update_coordinate(self.representation[i], coordinate)
            if self.molecule[i] == 'H':
                self.path.append(list(coordinate))

    def calc_mags(self):
        for coord in self.path:
            self.mags.append(self.__mag(coord))

    def get_mags(self):
        return self.mags

    def get_path(self):
        return self.path

    def add_parent_operator(self, operator):
        self.parent_operators.append(operator)
