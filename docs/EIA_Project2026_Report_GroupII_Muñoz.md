## Initialization, I/O & Analysis Workflow
**Contributor: Itxaso Muñoz-Aldalur (itxasoma)**

My contribution started from the initialization side of the project. I developed the `initial_conf` module to generate the starting conformation of the polyethylene chain, first at the only-carbon model level and then, when requested, with explicit hydrogens added for visualization and all-atom runs.

In parallel, I implemented the `io` module to keep the workflow simple and reproducible. This module reads the simulation parameters from `input.dat`, applies default values, strips comments, analizes user options such as `n_carbons`, `n_steps`, `explicit_h`, `conf_type`, and `rng_seed`, and writes the generated structures and trajectories in XYZ format.

## Initial Configuration Design
A central part of my work was the construction of different initial geometries for the polymer chain. The builder is based on internal coordinates, with fixed bond length, fixed bond angle, a selectable dihedral angle, and an overlap check to avoid unphysical self-intersections during chain creation.

I implemented several configuration modes with different physical purposes: fully planar all-trans, fully random dihedrals, a simple RIS-like discrete model, a spring/helix-like case (proposed by ai-murphy), and a perturbed trans case intended as a slight out-of-plane variation of the planar geometry. I also revised the hydrogen placement so that internal CH$_2$ hydrogens are no longer coplanar with the C-C-C backbone, but instead are placed in a more chemist-like arrangement, out-of-plane, consistent with a tetrahedral local geometry.

## Driver Improvements & Scientific Output
I also edited the serial driver program to make the production runs easier to analyze. In particular, I implemented an annealing schedule in the main Monte Carlo loop to study the approach to low-temperature equilibrium, and I added systematic string trimming in the output naming so that result files, trajectories, and plots carry a clean identifier of the simulation conditions.

In the same driver, I introduced `cpu_time` calls to monitor the execution time of the serial code. This gives a basic performance reference for later comparison with the MPI-parallel version, which is an important step to evaluate the actual benefit of parallelization.

On the analysis side, I improved the Python plotting script with a more scientific visual style and extended it to process all runs automatically. Besides energy, structural observables, and torsion distributions, I added an integrated autocorrelation-time estimate and an equilibration detector so that the torsion histograms are built only from the equilibrated production region. This follows standard practice based on autocorrelation analysis and automated equilibration detection aimed at maximizing the effectively uncorrelated sample count, based on [choderalab](https://www.choderalab.org/publications/2015/6/30/a-simple-method-for-automated-equilibration-detection-in-molecular-simulations).

## Parallelization
The current MPI contribution is intentionally minimal and focused on independent replicas rather than on parallelizing the energy evaluation itself. Starting from `main_serial.f90`, I produced a `main_parallel` version in which three MPI ranks are launched simultaneously, each one running the same Monte Carlo protocol with the same thermodynamic parameters but with a different initial configuration: `conf_type = 1`, `4`, and `5`.

This implementation preserves the structure of the serial code as much as possible. The only essential modifications are the MPI initialization/finalization, the rank-dependent assignment of the initial configuration and random seed, and the rank-labelled output files so that the three replicas can be compared directly without overwriting each other.