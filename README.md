# Lennard Jones System

## Link to repo: https://github.com/swetanjal/Lennard-Jones-System

## Specification of the system:

- Number of partciles (N) = 108
- Dimensions of box along X axis(Lx) = 18.0 A
- Dimensions of box along Y axis(Ly) = 18.0 A
- Dimensions of box along Z axis(Lz) = 18.0 A
- epsilon = 0.238 kcal/mol
- sigma = 3.4 A
- Initial constraint rij >= 3.4 A

## Problem Statement:

1. Generation of initial configuration
2. Calculating the negative gradient of the total potential energy function
3. Total potential energy minimization
4. Sample initial velocities from the Maxwell-Boltzmann velocity distribution
5. Solve the Hamilton's equations of motion and generate trajectories of atomic 
coordinates and velocities
6. Calculate mean square displacement, velocity correlation function, van Hove
correlation function, and dynamic structure factor

## Instructions to run code:

    # Changing into the src directory
    cd src
    
    python3 main.py [--minimize-potential] [--initial-config filename] [--initial-velocity filename] [--frames n_frames] [--temperature value_in_K]

    # Example Usage
    python3 main.py --minimize-potential --frames 250
    
When option minimize-potential is set to True, step 3 i.e total potential energy minimization step takes place and the configuration is saved after that in corresponding pdb files.

Use initial-config flag to load coordinates of particles from pdb file.

Use initial-velocity flag to load velocity of particles from pdb file.

### Instructions to generate mean square displacement plot:

    python3 msd.py ../outputs/argon_coordinates.pdb

### Instructions to generate velocity autocorrelation plot:

    python3 velocity_corr.py ../outputs/argon_velocities.pdb

### Instructions to generate Vann Hove correlation plot:

    python3 vann_hove.py ../outputs/argon_coordinates.pdb

### Instructions to generate Dynamic Structure Factor

    ./computeStructureFactor -i sample_s_k_full.in

    python3 plot_structure_factor.py ../outputs/traj_out.txt

### Instructions to generate Dynamic Structure Factor(Double Fourier Transform of Van Hove Correlation Fn.)

    python3 structure_factor.py ../outputs/argon_coordinates.pdb

## Description of files generated in outputs:

1. argon_coordinates.pdb : A PDB file containing the coordinates of trajectory of all N atoms.
2. argon_velocities.pdb : A PDB file containing the velocity trajectories of all N atoms.
3. initial_config_coords.xyz : A file in xyz format which contains coordinates of N atoms, generated randomly adhering to the constraints given in the problem statement. This generated file can be used to load initial coordinate configuration into the system.
4. initial_config_velocities.xyz: A file in xyz format which contains initial velocities of each of the N atoms sampled from Maxwell Boltzmann Velocity Distribution. This file can be used to load the inititial velocity configuration into the system.
5. minimized_config_coords.xyz : A file in xyz format which contains the coordinates of N atoms after the Potential Energy Minimization step.
6. traj.xyz : A file in xyz format which contains the coordinate trajectory of each of the N atoms.

## README: https://docs.google.com/document/d/1Ktl6CnwiAqsZnHV3XtD3a-qb3JYfnMJgpoVlE54ZoxU/edit?usp=sharing