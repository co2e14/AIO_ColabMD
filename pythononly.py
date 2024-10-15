import openmm as mm
import openmm.app as app
from openmm import unit
from sys import stdout
import requests
import pdbfixer
from openmm.app import PDBFile

def download_pdb(pdb_id):
    """Downloads a PDB file from the RCSB Protein Data Bank"""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        with open(f"{pdb_id}.pdb", "w") as f:
            f.write(response.text)
        print(f"PDB file {pdb_id} downloaded successfully.")
        return f"{pdb_id}.pdb"
    else:
        raise Exception(f"Failed to download PDB file {pdb_id}.")

def setup_simulation(pdb_file):
    """Sets up the simulation system from a PDB file using OpenMM and PDBFixer"""
    # Use PDBFixer to handle missing atoms, residues, and terminal groups
    fixer = pdbfixer.PDBFixer(pdb_file)
    
    # Find missing residues, termini, and atoms
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens()
    
    # Use the fixer to create a fixed PDB file with necessary additions
    with open(f"fixed_{pdb_file}", 'w') as output:
        PDBFile.writeFile(fixer.topology, fixer.positions, output)
    
    # Load the fixed PDB file into OpenMM
    pdb = PDBFile(f"fixed_{pdb_file}")
    
    # Force field
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # Create the system using the fixed PDB file
    system = forcefield.createSystem(pdb.topology, 
                                     nonbondedMethod=app.PME, 
                                     nonbondedCutoff=1.0*unit.nanometers, 
                                     constraints=app.HBonds)
    
    # Set up integrator
    integrator = mm.LangevinIntegrator(
        300*unit.kelvin,       # Temperature
        1.0/unit.picoseconds,  # Friction coefficient
        0.002*unit.picoseconds # Time step
    )
    
    # Add barostat for constant pressure simulation
    system.addForce(mm.MonteCarloBarostat(1.0*unit.atmospheres, 300*unit.kelvin, 25))
    
    # Set up the simulation
    simulation = app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    
    return simulation

def minimize_energy(simulation):
    """Minimizes the energy of the system"""
    print("Minimizing energy...")
    simulation.minimizeEnergy()
    print("Energy minimization complete.")

def run_simulation(simulation, n_steps=10000):
    """Runs the molecular dynamics simulation"""
    # Output simulation to a DCD trajectory file
    simulation.reporters.append(app.DCDReporter('trajectory.dcd', 1000))
    simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True,
                                                      potentialEnergy=True, temperature=True))
    
    print(f"Running simulation for {n_steps} steps...")
    simulation.step(n_steps)
    print("Simulation complete.")

def main(pdb_id, n_steps=10000):
    # Step 1: Download the PDB file
    pdb_file = download_pdb(pdb_id)
    
    # Step 2: Set up the simulation using PDBFixer
    simulation = setup_simulation(pdb_file)
    
    # Step 3: Minimize energy
    minimize_energy(simulation)
    
    # Step 4: Run simulation
    run_simulation(simulation, n_steps)

if __name__ == "__main__":
    pdb_id = input("Enter the PDB ID: ").strip()  # Example: '1UBQ'
    n_steps = int(input("Enter the number of simulation steps: ").strip())
    
    main(pdb_id, n_steps)
    