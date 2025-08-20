import os
import argparse
from ase.io import read
from ase.ga.data import DataConnection
import matplotlib.pyplot as plt
import numpy as np

def monitor_ga_runs(ion_type, lower_bound, upper_bound, save_path=None):
    """
    Monitors the convergence of genetic algorithm runs by plotting the
    lowest relative energy found versus the number of candidates evaluated.

    Args:
        ion_type (str): The chemical symbol of the intercalated ion (e.g., 'Li').
        lower_bound (int): The starting number of ions to analyze.
        upper_bound (int): The ending number of ions to analyze.
        save_path (str, optional): Path to save the plot image. If None,
                                   the plot is displayed interactively.
    """
    print("--- Starting GA Convergence Analysis ---")
    
    # --- Set up the plot ---
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(12, 8))

    # --- Loop through each ion number ---
    for i in range(lower_bound, upper_bound + 1):
        # Construct the path to the database file
        run_dir = f"{ion_type}-{i}"
        db_path = os.path.join(run_dir, 'gadb.db')
        
        print(f"\nAnalyzing run for {ion_type}{i} using database: {db_path}")

        if not os.path.isfile(db_path):
            print(f"  [Warning] Database not found: {db_path}. Skipping.")
            continue

        try:
            # Connect to the GA database
            da = DataConnection(db_path)
            all_candidates = da.get_all_relaxed_candidates()
        except Exception as e:
            print(f"  [Error] Could not read database {db_path}: {e}")
            continue

        if not all_candidates:
            print("  No relaxed candidates found in the database. Skipping.")
            continue

        # --- Extract confid and energy, then sort by confid ---
        candidate_data = []
        for cand in all_candidates:
            try:
                confid = cand.info['confid']
                energy = cand.get_potential_energy()
                candidate_data.append((confid, energy))
            except (KeyError, AttributeError):
                print(f"  [Warning] Skipping a candidate with missing 'confid' or energy.")
        
        if not candidate_data:
            print("  No valid candidates with confid and energy found. Skipping.")
            continue
            
        # Sort candidates based on their ID number
        candidate_data.sort(key=lambda x: x[0])

        # --- Find overall minimum energy for this run to set the 0 eV reference ---
        all_energies = [data[1] for data in candidate_data]
        min_energy_for_run = min(all_energies)

        # --- Process energies and find the running minimum relative energy ---
        plot_indices = []
        plot_relative_energies = []
        lowest_energy_so_far = float('inf')

        for index, (confid, energy) in enumerate(candidate_data):
            if energy < lowest_energy_so_far:
                lowest_energy_so_far = energy
            
            # Calculate the relative energy with respect to the run's best energy
            relative_energy = lowest_energy_so_far - min_energy_for_run
            
            # Use the sorted index (plus 1) for the x-axis
            plot_indices.append(index + 1)
            plot_relative_energies.append(relative_energy)

        # --- Add this run's data to the plot ---
        if plot_indices:
            label = f'{ion_type}{i} ({len(plot_indices)} candidates)'
            ax.plot(plot_indices, plot_relative_energies, marker='o', 
                    linestyle='-', markersize=4, label=label)
            print(f"  Successfully processed and plotted {len(plot_indices)} candidates.")
            print(f"  Lowest absolute energy found: {min_energy_for_run:.4f} eV")

    # --- Finalize and show/save the plot ---
    ax.set_title(f'GA Convergence for {ion_type} Intercalation', fontsize=16)
    ax.set_xlabel('Number of Candidates Evaluated', fontsize=12)
    ax.set_ylabel('Lowest Relative Energy Found (eV)', fontsize=12)
    
    # Only add legend if there are plotted lines
    if ax.get_legend_handles_labels()[0]:
        ax.legend(title="Ion Intercalation Count")
    else:
        print("\nNo data was successfully plotted.")

    fig.tight_layout()

    if save_path:
        try:
            plt.savefig(save_path, dpi=300)
            print(f"\nPlot successfully saved to {save_path}")
        except Exception as e:
            print(f"\n[Error] Failed to save plot: {e}")
    else:
        print("\nDisplaying plot...")
        plt.show()

if __name__ == '__main__':
    # --- Set up command-line argument parser ---
    parser = argparse.ArgumentParser(
        description="Monitor and plot the convergence of ASE GA runs.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'ion_type', 
        type=str, 
        help='The chemical symbol of the intercalated ion (e.g., Li, Na, Mg).'
    )
    parser.add_argument(
        'lower_bound', 
        type=int, 
        help='The starting number of ions (inclusive).'
    )
    parser.add_argument(
        'upper_bound', 
        type=int, 
        help='The ending number of ions (inclusive).'
    )
    parser.add_argument(
        '--save_plot', 
        type=str, 
        default=None, 
        metavar='FILENAME',
        help='Optional: Save the plot to a file (e.g., "convergence.png") instead of displaying it.'
    )
    
    args = parser.parse_args()
    
    monitor_ga_runs(args.ion_type, args.lower_bound, args.upper_bound, args.save_plot)
