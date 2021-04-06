"""
Specification:
N = 108
Lx = Ly = Lz = 18.0 A
epsilon = 0.238 kcal/mol
sigma = 3.4 A
initial constraint rij >= 3.4 A

Problem Statement:
(1) generation of initial configuration
(2a) calculating the negative gradient of the total potential energy function
(2b) total potential energy minimization
(3) sample initial velocities from the Maxwell-Boltzmann velocity distribution
(4) Solve the Hamilton's equations of motion and generate trajectories of atomic 
coordinates and velocities
(5) calculate mean square displacement, velocity correlation function, van Hove
correlation function, and dynamic structure factor
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import random

# Defining some constants
Lx = 18
Ly = 18
Lz = 18
N = 108
epsilon = 0.238
sigma = 3.4 # Diameter of the particle

# Coordinates of all the particles
coords = []

def compute_radius(x1, y1, z1, x2, y2, z2):
    return math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2))

def generate_random_config():
    print("Generating initial configuration...")
    for i in range(N):
        got = 1
        while got:
            x = random.random() * Lx
            y = random.random() * Ly
            z = random.random() * Lz
            got = 0
            for i in range(len(coords)):
                # Check if overlaps with any particles or not.
                if compute_radius(coords[i][0], coords[i][1], coords[i][2], x, y, z) >= 3.4 and\
                compute_radius(coords[i][0] - Lx, coords[i][1] , coords[i][2], x, y, z) >= 3.4 and\
                compute_radius(coords[i][0] + Lx, coords[i][1] , coords[i][2], x, y, z) >= 3.4 and\
                compute_radius(coords[i][0], coords[i][1] - Ly, coords[i][2], x, y, z) >= 3.4 and\
                compute_radius(coords[i][0], coords[i][1] + Ly, coords[i][2], x, y, z) >= 3.4 and\
                compute_radius(coords[i][0], coords[i][1], coords[i][2] - Lz, x, y, z) >= 3.4 and\
                compute_radius(coords[i][0], coords[i][1], coords[i][2] + Lz, x, y, z) >= 3.4:
                    continue
                else:
                    got = 1
                    break
            if got == 0:
                coords.append([x, y, z])

def create_xyz_file(filename):
    print("Creating xyz file of coordinates...")
    f = open(filename, "w")
    f.write(str(N) + "\n")
    f.write("System of Argon Gas\n")
    for i in range(N):
        f.write("C " + str(coords[i][0]) + " " + str(coords[i][1]) + " " + str(coords[i][2]) + "\n")
    f.close()

if __name__ == "__main__":
    print("Lennard Jones System")
    generate_random_config()
    create_xyz_file("abc.xyz")

    
    exit()