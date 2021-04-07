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

def compute_pairwise_potential(r):
    return 4 * epsilon * (math.pow((sigma * 1.0 / r), 12) - math.pow((sigma * 1.0 / r), 6))

def minimum(a, b, c, d, e, f, g):
    return min(a, min(b, min(c, min(d, min(e, min(f, g))))))

def min_image_distance(x1, y1, z1, x2, y2, z2):
    return minimum(compute_radius(x1, y1, z1, x2, y2, z2),\
    compute_radius(x1 - Lx, y1 , z1, x2, y2, z2),\
    compute_radius(x1 + Lx, y1 , z1, x2, y2, z2),\
    compute_radius(x1, y1 - Ly, z1, x2, y2, z2),\
    compute_radius(x1, y1 + Ly, z1, x2, y2, z2),\
    compute_radius(x1, y1, z1 - Lz, x2, y2, z2),\
    compute_radius(x1, y1, z1 + Lz, x2, y2, z2))

def compute_total_potential_energy():
    res = 0.0
    for i in range(N):
        for j in range(i + 1, N):
            res = res + compute_pairwise_potential(min_image_distance(coords[j][0], coords[j][1], coords[j][2], coords[i][0], coords[i][1], coords[i][2]))
    return res

def get_minimum_image_coords(x1, y1, z1, x2, y2, z2):
    min_d = 1000000000000000
    coords = []
    r = compute_radius(x1, y1, z1, x2, y2, z2)
    if r < min_d:
        min_d = r
        coords = [x1, y1, z1]
    r = compute_radius(x1 - Lx, y1 , z1, x2, y2, z2)
    if r < min_d:
        min_d = r
        coords = [x1 - Lx, y1, z1]
    r = compute_radius(x1 + Lx, y1 , z1, x2, y2, z2)
    if r < min_d:
        min_d = r
        coords = [x1 + Lx, y1, z1]
    r = compute_radius(x1, y1 - Ly, z1, x2, y2, z2)
    if r < min_d:
        min_d = r
        coords = [x1, y1 - Ly, z1]
    r = compute_radius(x1, y1 + Ly, z1, x2, y2, z2)
    if r < min_d:
        min_d = r
        coords = [x1, y1 + Ly, z1]
    r = compute_radius(x1, y1, z1 - Lz, x2, y2, z2)
    if r < min_d:
        min_d = r
        coords = [x1, y1, z1 - Lz]
    r = compute_radius(x1, y1, z1 + Lz, x2, y2, z2)
    if r < min_d:
        min_d = r
        coords = [x1, y1, z1 + Lz]
    return coords

def compute_gradient():
    grad_xs = []
    grad_ys = []
    grad_zs = []
    for i in range(N):
        gradx = 0.0
        grady = 0.0
        gradz = 0.0
        for j in range(i + 1, N):
            new_coords = get_minimum_image_coords(coords[j][0], coords[j][1], coords[j][2], coords[i][0], coords[i][1], coords[i][2])    
            r = compute_radius(coords[i][0], coords[i][1], coords[i][2], new_coords[0], new_coords[1], new_coords[2])
            gradx = gradx + 24 * epsilon * math.pow((sigma * 1.0 / r), 6) * (1.0 / (r * r)) * \
            (2 * math.pow((sigma * 1.0 / r), 6) - 1) * (coords[i][0] - new_coords[0])

            grady = grady + 24 * epsilon * math.pow((sigma * 1.0 / r), 6) * (1.0 / (r * r)) * \
            (2 * math.pow((sigma * 1.0 / r), 6) - 1) * (coords[i][1] - new_coords[1])

            gradz = gradz + 24 * epsilon * math.pow((sigma * 1.0 / r), 6) * (1.0 / (r * r)) * \
            (2 * math.pow((sigma * 1.0 / r), 6) - 1) * (coords[i][2] - new_coords[2])
        grad_xs.append(gradx)
        grad_ys.append(grady)
        grad_zs.append(gradz)   
    return [grad_xs, grad_ys, grad_zs]

def minimize_potential():
    print("Minimizing the potential...")
    # Current config which needs to be optimised...
    learning_rate = 0.01
    prev_potential = compute_total_potential_energy()
    itr = 1
    while True:
        grads = compute_gradient()
        for i in range(N):
            coords[i][0] = coords[i][0] + learning_rate * grads[0][i]
            coords[i][1] = coords[i][1] + learning_rate * grads[1][i]
            coords[i][2] = coords[i][2] + learning_rate * grads[2][i]
        curr_potential = compute_total_potential_energy()
        print("Total LJ Potential Energy after iteration", itr, "=", curr_potential)
        itr = itr + 1
        if abs(curr_potential - prev_potential) <= 1e-3:
            break
        prev_potential = curr_potential
    
if __name__ == "__main__":
    print("Lennard Jones System")
    generate_random_config()
    create_xyz_file("abc.xyz")
    print("Current total potential energy =", compute_total_potential_energy())
    minimize_potential()
    print("Total Potential Energy after Energy Minimization step =", compute_total_potential_energy())
    exit()