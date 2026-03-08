============================================================
MC SIMULATION TESTING & VALIDATION
============================================================

PROJECT INTENT:
This simulation serves as a comprehensive integration test 
to validate the interplay between the Monte Carlo engine, 
the simulated annealing protocol, and the custom energy 
forcefield (Lennard-Jones + Torsion). The goal was to 
verify the implementation of all mechanisms.

SYSTEM PARAMETERS:
- Atoms: 250 Carbon chain 
- Total Steps: 100.000.000 (10^8)
- Max Rotation (delta): 0.20 radians

ANNEALING PROTOCOL:
- Initial Temperature: 1000.0 K
- Final Temperature: 0.1 K
- Cooling Schedule: Linear decrement per step

KEY RESULTS:
- Final State: Compact globule (Energy < 0)
- Initial Rg: ~27.3 Å  -->  Final Rg: ~9.5 Å
- Integration Status: SUCCESSful. 
- Competition between torsional stiffness and 
  Van der Waals attraction
============================================================
