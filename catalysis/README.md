## Scripts for Catalysis Projects

----------
### About

`combine_slab_gas.py`: 

- Places gas from `gas.vasp` into slab and its box from `slab.vasp`, and intializes velocities for MD.
- Velocity initialization:
    - Slab: from Maxwell-Boltzmann at `slab_temp` (Kelvin)
    - Gas:  can separate vibrational, rotational, and translational velocities using `gas_vibration_temp`, `gas_rotation_temp`, and `gas_translation_temp`, all drawn from Maxwell-Boltzmann
- Gas placed `separation` (Angstroms) above slab at a random position after a random rotation to change orientation
- Finally, the gas's z-velocity is drawn from the Maxwell-Flux distribution to simulate gas hitting the slab surface
    - Details from [this old paper](https://theory.cm.utexas.edu/henkelman/pubs/sharia14_074706.pdf)
- Outputs `slab_gas.vasp`
