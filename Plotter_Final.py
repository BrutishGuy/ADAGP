#program to visualize an HP molecule stored in a specified file
#victor gueorguiev
#01-10-2015

#relevant imports
import numpy as np
from visual import *

#NB   #NB
filename = '200_30HH.txt'

#function to translate file entries int ovector coordinates
def get_translation(direction):
    d = direction[0]
    if d == 'U':
        return vector(0, step_size, 0)
    elif d == 'D':
        return vector(0, -step_size, 0)
    elif d == 'R':
        return vector(step_size, 0, 0)
    elif d == 'L':
        return vector(-step_size, 0, 0)
    elif d == 'F':
        return vector(0, 0, step_size)
    else:
        return vector(0, 0, -step_size)

#function to return the color of the H or P monomer
def get_color(molecule):
    m = molecule[0]
    if m == 'H':
        return color.blue
    else:
        return color.white

#set some defaults for the display window
scene = display(title='HP Model Plot', x = 0, y = 0, width = 800, height = 800, center = (0, 0, 0), background = color.white)

f = open(filename, 'r')
header = f.readline() #remove header of file
N=50
step_size = 3
#data from file stored in array for later use
molecule_data = []
direction_data = []
#3D objects that will appear on screen stored in an array
molecules = []
links = []
#variables to help calculate positions of 3D objects
current_pos = vector(0, 0, -step_size)
next_pos = vector(0, 0, 0)
prev_pos = vector(0, 0, -step_size)
#list of position coordinates for later use
coordinates = []
#counter
i = 0
prev_color = color.blue

for line in f:
    columns = line.split('\t')
    print(current_pos,columns[0])
    if i == 0:
        molecule_data.append(columns[0])
        molecules.append(sphere(pos = current_pos, radius = 0.6, color = get_color(columns[0])))
        if columns[0][0] == 'H':
            coordinates.append(vector(current_pos))
    else:
        prev_pos = vector(current_pos)
        current_pos += get_translation(columns[1])
        molecule_data.append(columns[0])
        molecules.append(sphere(pos = current_pos, radius = 0.6, color = get_color(columns[0])))
        direction_data.append(columns[1])
        if columns[0][0] == 'H':
            coordinates.append(vector(current_pos))
    links.append(cylinder(pos = prev_pos, axis = current_pos - prev_pos, radius = 0.15, color = color.white))
    i += 1

print(direction_data)
#reset certain variables to determine links for hydrophobic core/monomers
coordinate = vector(0, 0, -step_size)
prev_coordinate = vector(0, 0, -step_size)
next_coordinate = vector(0, 0, -step_size) + get_translation(direction_data[0])
print(next_coordinate)
next_monomer = molecule_data[1]
counter = 0
for i in range(len(molecule_data)):
    monomer = molecule_data[i]
    prev_coordinate = vector(coordinate)

    if molecule_data[i][0] == 'H':
        for j in range(3):
            test = vector(coordinate)
            test[j] -= step_size
            if test in coordinates and test != next_coordinate and test != prev_coordinate:
                links.append(helix(pos = coordinate, axis = test - coordinate, radius = 0.15, color = color.blue, coils = 20, thickness = 0.15/5))
                counter += 1
            test = vector(coordinate)
            test[j] += step_size
            if test in coordinates and test != next_coordinate and test != prev_coordinate:
                links.append(helix(pos = coordinate, axis = test - coordinate, radius = 0.15, color = color.blue, coils = 20, thickness = 0.15/5))
                counter += 1
        coordinates.remove(coordinate)
    
    if 0 <= i < len(molecule_data) - 2:
        coordinate += get_translation(direction_data[i])
        next_monomer = molecule_data[i+1]
        next_coordinate += get_translation(direction_data[i+1])
    elif i < len(molecule_data) - 1:
        coordinate += get_translation(direction_data[i])
        next_monomer = molecule_data[i+1]
        next_coordinate = vector(coordinate)
print(counter)
