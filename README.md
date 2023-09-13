# Masterarbeit "Computing lines on tropical cubic surface"

This repository contains the code belonging to the master thesis "Computing lines on tropical cubic surfaces" written by Lena Weis 
and supervised by Prof. Dr. Michael Joswig (TU Berlin) and Dr. Marta Panizzut (MPI Mathematics in Sciences). It was published in order to 
make the results of the thesis FAIR, i.e. findable, accessible, interoperable and reusable. 

## Problem Definition
The aim of this thesis was to compute the number of lines on a surface of one of five different combintorial types of smooth tropical cubic surface. 
They are
- #37
- #842554
- #2091566
- #5054117
- #12369387

## Reusing the results
The Schläfli fans can be found [here](Computations/Schlaefli_fans).

## Reproducing the results
### 1. Reproduce the Schläfli walls of a triangulation
   First include the secondary cone and Motifs of the triangulation of interest. They can be found [here](Computations/Code/info_triangulation).
   Include the algorithm to compute the [Schläfli walls](Computations/Code/schlaefliwalls.jl) are [print](Computations/Code/print_schlaefli_walls.jl) them with this. 

### 2. Verify the number of lines of examplary surfaces
   Include this [file](Computations/Code/compute_number_of_lines.jl) and procede as in the following example
   ```
   Which triangulation?
   2091566
   What are the coefficients of f defining the surface? (Write as x0 x1 x2...)
   19 -4 -13 -12 2 -7 -11 3 -10 12 -5 -28 -27 -6 -10 8 8 -10 7 21
   The corresponding surface has this many visible motifs of motifs 3A, 3B, 3C, 3D, 3E, 3H, i.e. lines:
   21
   ```
   Note that you can only check this for one of the five triangulations above. For triangulation #842554 we have included one [cone of the Schläfli fan](Computations/Code/cones_842554)
   such that the number of lines on the surfaces dual to the points in each cone is exactly 14, 15, 17, 18, 19, 20 or 21.

### 3. Check the remaining statement as the number of lines and the Schläfli pillars
   !MISSING
   
