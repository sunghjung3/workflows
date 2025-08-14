import numpy as np
from ase.io import read, write
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
from ase.units import kB, Ang, fs
from ase import Atoms
import warnings
from collections import Counter

#================================================================
### Gas velocity initialization helper functions

# Sampling from Maxwell Flux Distribution
def mf_inverse_cdf(u, c):
    """
    Calculates the random velocity from a uniform random number u
    using the inverse CDF of Maxwell Flux.

    MF PDF: (mv)/(kT) * exp(-(mv**2)/(2kT))
    MF CDF: 1 - exp(-(mv**2)/(2kT))
    MF CDF^-1 for velocity: sqrt(-(2kT)/m ln(1-u))
    
    Args:
        u (float or array): Uniform random number(s) in (0, 1).
        c (float or array): kB*T/m (in eV/amu)
        
    Returns:
        float or array: The random sample(s) from your distribution.
    """
    return np.sqrt(-2*c * np.log(1 - u))


def sample_rotational_velocity(atoms, temperature_K):
    """
    Sample rotational velocities for a molecule from the Maxwell-Boltzmann
    distribution at a given temperature.

    Args:
        atoms (ase.Atoms): Atoms object representing the molecule.
        temperature_K (float): Temperature in Kelvin for the Maxwell-Boltzmann distribution.
    Returns:
        numpy.ndarray: The sampled rotational velocities for the molecule.
    """
    # Find the principal moments of inertia and principal axes basis vectors
    Ip, basis = atoms.get_moments_of_inertia(vectors=True)
    omega = np.zeros(3)

    # Sample rotational velocities for each principal axis
    for i in range(3):
        if Ip[i] < 1e-3:  # Avoid division by zero
            warnings.warn(f"Principal moment of inertia {Ip[i]} is too small, setting omega[{i}] to 0.")
            omega[i] = 0.0
            continue
        mb_std = np.sqrt(kB * temperature_K / Ip[i])  # sqrt(eV/(amu*Ang^2))
        omega[i] = np.random.normal(0, mb_std)
    
    # Convert back to Cartesian coordinates
    return omega @ basis

#================================================================
### Other helper functions

def check_poscar_element_groups(filename):
    vasp_keyword = ['vasp', 'poscar', 'contcar']
    if any(keyword in filename.lower() for keyword in vasp_keyword):
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Get the 6th line (index 5 in a 0-indexed list)
        sixth_line = lines[5]

        # Strip leading/trailing whitespace and split the line by spaces
        words = sixth_line.strip().split()

        # Use collections.Counter to count occurrences of each word
        word_counts = Counter(words)
        
        # Find words that appear more than once
        duplicates = {word: count for word, count in word_counts.items() if count > 1}

        if duplicates:
            print()
            print(f"***SAME ELEMENT FOUND IN MULTIPLE GROUPS IN {filename}***\n"
                   "Manual postprocessing is required to separate them again.")
            print()
#================================================================


def combine_slab_and_gas_with_velocities(
    slab_file="slab.vasp", 
    gas_file="gas.vasp", 
    output_file="slab_gas.vasp", 
    separation=5.0,
    slab_temp=300.0,
    gas_vibration_temp=3000.0,
    gas_rotation_temp=300.0,
    gas_translational_temp=300.0
):
    """
    Reads slab and gas files, initializes their velocities for an NVE
    bombardment simulation, combines them, and writes a new VASP file.

    - Slab atoms get Maxwell-Boltzmann velocities at slab_temp.
    - Gas molecule is randomly rotated.
    - Gas atoms get vibrational, rotational, and translational velocities drawn
        from the Maxwell-Boltzmann distribution at gas_..._temp.
    - Gas center-of-mass velocity in z is drawn from a Maxwell-Flux distribution
      at gas_translational_temp and added to the internal velocities.
    """
    try:
        slab = read(slab_file)
        gas = read(gas_file)
        print(f"Successfully read '{slab_file}' and '{gas_file}'.")
    except FileNotFoundError as e:
        print(f"Error: Could not find a required file. {e}")
        return

    # --- 1. Initialize Slab Velocities ---
    MaxwellBoltzmannDistribution(slab, temperature_K=slab_temp)
    print(f"Initialized slab velocities to {slab_temp} K.")

    # --- 2. Initialize Gas Velocities and Orientation ---
    total_gas_mass = np.sum(gas.get_masses())  # amu
    if len(gas) == 1:
        warnings.warn("Gas molecule has only one atom.")
        MaxwellBoltzmannDistribution(gas, temperature_K=gas_translational_temp)
        print(f"Initialized gas translational velocities to {gas_translational_temp} K.")
    else:
        # Apply a random rotation to the gas molecule
        gas.rotate(np.random.rand() * 360, 'x', rotate_cell=False)
        gas.rotate(np.random.rand() * 360, 'y', rotate_cell=False)
        gas.rotate(np.random.rand() * 360, 'z', rotate_cell=False)
        print("Applied random rotation to the gas molecule.")

        # Set gas vibrational velocities at vibrational temperature
        MaxwellBoltzmannDistribution(gas, temperature_K=gas_vibration_temp)
        print(f"Initialized gas vibrational velocities to {gas_vibration_temp} K.")

        # Zero center-of-mass (COM) linear and angular momenta for gas
        Stationary(gas, preserve_temperature=False)
        ZeroRotation(gas, preserve_temperature=False)
        print("Zeroed gas center-of-mass momenta.")

        # Set gas rotational velocities
        omega = sample_rotational_velocity(gas, temperature_K=gas_rotation_temp)
        pos_offset = gas.get_positions() - gas.get_center_of_mass()
        gas.set_velocities(gas.get_velocities() + np.cross(omega, pos_offset))
        print(f"Set gas rotational velocities at {gas_rotation_temp} K.")

        # Set gas translational velocities
        mb_std = np.sqrt(kB * gas_translational_temp / total_gas_mass)  # sqrt(eV/amu)
        gas.set_velocities(gas.get_velocities() + np.random.normal(0, mb_std, 3))
        print(f"Set gas translational velocities at {gas_translational_temp} K.")


    # --- 3. Set Gas Center-of-Mass (COM) Velocity from Maxwell-Flux ---
    # Draw from Maxwell-Flux distribution
    kBT_over_m = kB * gas_translational_temp / total_gas_mass  # eV/amu (unit of v**2)
    com_vz_target = mf_inverse_cdf( np.random.rand(), kBT_over_m )  # unit of (eV/amu)**0.5
    com_vz_target *= -1  # negative z direction to collide with slab
    print(f"Target COM velocity in z (Maxwell-Flux at {gas_translational_temp} K): {com_vz_target.round(3)} (eV/amu)**0.5")

    # --- 4. Adjust Atomic Velocities in z to Match Target COM Velocity ---
    # Get the current COM velocity from the internal thermal motion
    gas_momenta = gas.get_momenta()
    com_vz_current = (np.sum(gas_momenta, axis=0) / total_gas_mass)[2]
    
    # Calculate the z velocity correction needed per atom
    velocity_correction = com_vz_target - com_vz_current
    
    # Add this correction to each atom in the gas molecule
    current_velocities = gas.get_velocities()
    current_velocities[:, 2] += velocity_correction
    gas.set_velocities(current_velocities)
    print("Adjusted atomic velocities to match target COM velocity.")

    # --- 5. Position the Gas Molecule Above the Slab ---
    max_z_slab = slab.get_positions()[:, 2].max()
    min_z_gas = gas.get_positions()[:, 2].min()
    z_shift = (max_z_slab + separation) - min_z_gas
    
    gas_com = gas.get_center_of_mass()
    xy_shift_direct_coord = np.random.rand(3)
    xy_shift_direct_coord[2] = 0
    xy_shift = xy_shift_direct_coord @ slab.cell - gas_com

    gas.translate([xy_shift[0], xy_shift[1], z_shift])
    print(f"Placed gas molecule {separation:.1f} Å above the slab.")

    # --- 6. Combine and Write ---
    combined_system = slab + gas
    combined_system.set_cell(slab.get_cell())
    combined_system.set_pbc(slab.get_pbc())

    write(output_file, combined_system, format='vasp', vasp5=True)
    print(f"Successfully wrote combined system with velocities to '{output_file}'.")

    # Give warning if the orignal slab and/or gas file had an element in multiple groups.
    check_poscar_element_groups(gas_file)
    check_poscar_element_groups(slab_file)


if __name__ == "__main__":
    combine_slab_and_gas_with_velocities()
