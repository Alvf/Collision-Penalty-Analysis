# A Unified Analysis of Collision Penalty Energies - Codebase

Hello! This repository demonstrates the analysis of various collision penalties as described in the SCA 2023 paper [A Unified Analysis of Collision-Penalty Energies](https://www.alvin.pizza/unified-analysis-penalty-energies), by Alvin Shi and Theodore Kim.

The simulator's basic structure is based off of HOBAK, from the SIGGRAPH 2022 Course *[Dynamic Deformables: Implementation and Production Practicalities (Now With Code!)](http://www.tkim.graphics/DYNAMIC_DEFORMABLES/)*, by Theodore Kim and David Eberle.

## Running The Simulation
Run one of `make linux` or `make mac` (depending on your system) to compile dependencies and run the hammock simulation in `projects/simulateHammock`. 
Controls (much like the HOBAK interface) are as follows: 

 - `a` will start and stop the animation.
 - `q` will quit and write out a QuickTime movie of the animation. You need to have [FFMPEG](https://ffmpeg.org/)
installed for the movie write to succeed. Note that JSON scene descriptions are not available.
 - `Q` or `Esc` will quit without saving anything.

## Matlab Verification
Matlab scripts for testing the energies described in the paper are located in `Matlab_Verification.` Running any of the `Verify_[energy/length]` notebooks will validate analytic expressions against numerical differences.

## C++ Unit Tests
You can run unit tests for the C++ energy implementations like so:

    cd projects/unitTests
    make depend
    make
    cd ../../bin
    ./unitTests

We put various vertex-face/edge-edge collision frame/length/energy arrangements through a finite-differencing comparison for gradients and Hessians. 

## Collision Eigenanalysis
The files located in `src/Collision` contain C++ implementations of general energy eigenpairs and
secondary term handling. Feel free to swap out different frames, 1D energies, and lengths to get different qualitative results by modifying/making your own `setCollision()`-like functions in `src/Scenes/SIMULATION_SCENE.h`.

## FAQ
### Where is include_top.mk?
It gets made when `make linux` or `make mac` are run!
