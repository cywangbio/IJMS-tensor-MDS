from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import mdtraj as md
import os


def reline(inputpdb):
    pdb = PDBFile(f"{inputpdb}_fixed.pdb")
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.deleteWater()
    modeller.addExtraParticles(forcefield)
    modeller.addHydrogens(forcefield)
    modeller.addSolvent(forcefield, padding=1.0 * nanometer)

    system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1 * nanometer,
                                     constraints=HBonds)
    integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.004 * picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    simulation.reporters.append(PDBReporter(f'{inputpdb}_output_water.pdb', 800))
    simulation.reporters.append(StateDataReporter(stdout, 800, step=True, potentialEnergy=True, temperature=True))
    simulation.reporters.append(StateDataReporter(f"{inputpdb}_md_log_long.txt", 800, step=True, potentialEnergy=True, temperature=True))
    simulation.step(800000)
    os.system(f'grep -v HOH {inputpdb}_output_water.pdb > {inputpdb}_output.pdb')
    os.remove(f"{inputpdb}_output_water.pdb")
    os.system(f'grep -v CL {inputpdb}_output.pdb > {inputpdb}_nocl.pdb')
    os.remove(f"{inputpdb}_output.pdb")
    os.system(f'grep -v NA {inputpdb}_nocl.pdb > {inputpdb}_outputf.pdb')
    os.remove(f"{inputpdb}_nocl.pdb")
    traj = md.load(f"{inputpdb}_outputf.pdb", top=f"{inputpdb}_fixed.pdb")
    protein = [set(traj.topology.residue(0).atoms)]
    correct = traj.image_molecules(anchor_molecules=protein)
    correct.save_pdb(f"{inputpdb}_final.pdb")
    os.remove(f"{inputpdb}_outputf.pdb")


location = input('Give me conda absolute path (e.g.: /home/username/miniconda3)').strip()
line = input('Give me files: ').split()
line = line[0]
goto = "/".join(line.split("/")[:-1])
name = line.split("/")[-1][:-4]
os.chdir(goto)
os.system(f'{location}/bin/pdbfixer {name}.pdb --replace-nonstandard --add-residues --output={name}_fixed.pdb')
try:
    reline(name)
except Exception as error:
    print(error)
    os.remove(f"{name}_output_water.pdb")
    try:
        reline(name)
    except Exception as error:
        print(error)
        os.remove(f"{name}_output_water.pdb")
        try:
            reline(name)
        except Exception as error:
            print(error)
            os.remove(f"{name}_output_water.pdb")
            try:
                reline(name)
            except Exception as error:
                print(error)
                os.remove(f"{name}_output_water.pdb")
                try:
                    reline(name)
                except Exception as error:
                    print(error)
                    os.remove(f"{name}_output_water.pdb")
                    try:
                        reline(name)
                    except Exception as error:
                        print(error)
                        os.remove(f"{name}_output_water.pdb")
                        try:
                            reline(name)
                        except Exception as error:
                            print(error)
                            os.remove(f"{name}_output_water.pdb")
                            try:
                                reline(name)
                            except Exception as error:
                                print(error)
                                os.remove(f"{name}_output_water.pdb")
