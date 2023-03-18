#!/usr/bin/env python
import glob
import os
import shutil
import subprocess
from warnings import warn

from abfe.home import home
from typing import Union, Iterable
import tarfile

from toff import Parameterize
from abfe.scripts.preparation.topology_fixer import fix_topology
import BioSimSpace as bss
# from pdbfixer import PDBFixer
# from openmm.app import PDBFile
PathLike = Union[os.PathLike, str, bytes]

def get_gmx_ff(ff_code:str, out_dir:PathLike = '.') -> PathLike:
    """Get GROMACS Force Field


    Parameters
    ----------
    ff_code : PathLike
        The identification of the gromacs force field.
        For now only: Slipids_2020 and amber99sb-star-ildn are supported.
    out_dir : PathLike, optional
        Where the file will be decompress, by default '.'
    """
    out_dir = os.path.abspath(out_dir)
    supported_ff = [
        'Slipids_2020',
        'amber99sb-star-ildn',
    ]
    if ff_code not in supported_ff:
        raise ValueError(f"ff_code = {ff_code} is not valid. Chose between: {supported_ff}")
    else:
        fname = os.path.join(home(dataDir='gmx_ff'), f'{ff_code}.ff.tar.gz')
    tar = tarfile.open(fname, "r:gz")
    tar.extractall(out_dir)
    tar.close()
    return os.path.join(out_dir,  f'{ff_code}.ff')

def run(command:str, shell:bool = True, executable:str = '/bin/bash', Popen:bool = False) -> subprocess.CompletedProcess:
    """A simple wrapper around subprocess.Popen/subprocess.run

    Parameters
    ----------
    command : str
        The command line to be executed
    shell : bool, optional
        Create a shell section, by default True
    executable : str, optional
        what executable to use, pass `sys.executable` to check yours, by default '/bin/bash'
    Popen : bool, optional
        Use `Popen` (the PID could be access) instead of `run`, by default False

    Returns
    -------
    subprocess.CompletedProcess
        The process

    Raises
    ------
    RuntimeError
        In case that the command fails, the error is raised in a nice way
    """
    #Here I could make some modification in order that detect the operator system
    #NAd make the command compatible with the operator system
    #the function eval could be an option if some modification to the variable command
    #need to be done.... SOme flash ideas...

    if Popen:
        #In this case you could access the pid as: run.pid
        process = subprocess.Popen(command, shell = shell, executable = executable, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
    else:
        process = subprocess.run(command, shell = shell, executable = executable, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
        returncode = process.returncode
        if returncode != 0:
            print(f'Command {command} returned non-zero exit status {returncode}')
            raise RuntimeError(process.stderr)
    return process

def system_combiner(**md_elements):
    """This function simply sum up all the elements provided 
    as keyword arguments.

    Returns
    -------
    object
        any Python object with the method sum implemented. In case elements
        that evaluate as False in Python will not be taken into account:
        E.g. False, 0, '', None

    Raises
    ------
    RuntimeError
        In case all the elements evaluate as False
    """
    if any(md_elements.values()):
        # md_system = sum(element for element in md_elements.values() if element) # it does not work with sum
        for element in md_elements:
            if md_elements[element]:
                try:
                    md_system += md_elements[element]
                except NameError:
                    md_system = md_elements[element] 
    else:
        raise RuntimeError(f"system_combiner failed with the inputs: {md_elements}")
    print(f"The system was constructed as fallow: {' + '.join([key for key in md_elements if md_elements[key]])}")
    return md_system

# TODO, check what is the type of the bss_systems to add it as a HintType
def solvate(bss_system:object, out_dir:PathLike = '.', box:Iterable[float] = None, angles:Iterable[float] = None):
    """Solvate and add ions to the system, if box and angles are not provided,
    the system will be solvated as a truncated octahedron with a padding of 15 Angstroms.

    Parameters
    ----------
    bss_system : object
        The BSS system to solvate
    out_dir : PathLike, optional
        Where the files will be written: solvated.gro, solvated.top, by default '.'
    box : Iterable[float], optional
        This is the vectors of the bos in ANGSTROMS. It is important that the provided vector has the correct units, by default None
    angles : Iterable[float], optional
        This is the angles between the components of the vector in DEGREEs. It is important that the provided vector has the correct units, by default None

    Raises
    ------
    ValueError
        if box does not have three elements
    ValueError
        if angles does not have three elements
    """
    print("Solvating the system in: ", out_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    if box and angles:
        if len(box) != 3:
            raise ValueError(f"'box' must be a iterable of three values: [x,y,z] in Angstrom. Provided: box = {box}")
        elif len(angles) != 3:
            raise ValueError(f"'angles' must be a iterable of three values: [xy,yz,zx] in Degree. Provided: angles = {angles}")
        box = [bss.Types.Length(proj, 'angstrom') for proj in box]
        angles = [bss.Types.Angle(angle, 'degree') for angle in angles]
    else:
        box_min, box_max = bss_system.getAxisAlignedBoundingBox()
        box_size = [y - x for x, y in zip(box_min, box_max)]
        padding = 15 * bss.Units.Length.angstrom
        box_length = (max(box_size) + 1.5 * padding)
        box, angles = bss.Box.truncatedOctahedron(box_length.value() * bss.Units.Length.angstrom)
    
    solvated = bss.Solvent.tip3p(bss_system, box=box, angles=angles)
    
    cwd = os.getcwd()
    os.chdir(out_dir)
    bss.IO.saveMolecules('solvated', solvated, ["GroTop", "Gro87"])
    os.chdir(cwd)

def prepare_for_abfe(out_dir, ligand_dir, sys_dir):
    complex_out = out_dir + "/complex"
    ligand_out = out_dir + "/ligand"
    if (not os.path.exists(complex_out)): os.makedirs(complex_out)
    if (not os.path.exists(ligand_out)): os.makedirs(ligand_out)

    shutil.copyfile(src=sys_dir + "/solvated.gro", dst=complex_out + "/complex.gro")
    shutil.copyfile(src=sys_dir + "/solvated_fix.top", dst=complex_out + "/complex.top")

    for itp_file in glob.glob(sys_dir + "/*.itp"):
        shutil.copy(src=itp_file, dst=complex_out)

    shutil.copyfile(src=ligand_dir + "/solvated.gro", dst=ligand_out + "/ligand.gro")
    shutil.copyfile(src=ligand_dir + "/solvated_fix.top", dst=ligand_out + "/ligand.top")

    for itp_file in glob.glob(ligand_dir + "/*.itp"):
        shutil.copy(src=itp_file, dst=ligand_out)

# TODO create docs of the __init__ method
class MakeInputs:
    def __init__(self,
            protein_pdb:PathLike = None,
            membrane_pdb:PathLike = None,
            cofactor_mol:PathLike = None,
            hmr_factor:float = None,
            box:Iterable[float] = None,
            angles:Iterable[float] = None,
            keep_tmp_files_on:PathLike = None):

        self.protein_pdb = protein_pdb
        self.membrane_pdb = membrane_pdb
        self.cofactor_mol = cofactor_mol
        self.hmr_factor = hmr_factor
        self.box = box
        self.angles = angles
        self.keep_tmp_files_on = keep_tmp_files_on
        self.__self_was_called = False
        
        if self.membrane_pdb:
            if not self.box:
                warn(f"It looks that your are trying to set up a membrane system (membrane_pdb = {self.membrane_pdb}).\n"\
                    "However you have not provided a 'box' definition, the solvated system could not be right.\n"\
                    )
            if not self.angles:
                warn(f"It looks that your are trying to set up a membrane system (membrane_pdb = {self.membrane_pdb}).\n"\
                    "However you have not provided a 'angles' definition, the solvated system could not be right.\n"\
                    )
        # Working paths
        if keep_tmp_files_on:
            self.wd = keep_tmp_files_on
        else:
            self.wd = '.builder'
        if not os.path.exists(self.wd):
            os.makedirs(self.wd)
        self.wd = os.path.abspath(self.wd)

    def openff_process(self, mol_file:PathLike, name:str="MOL", safe_naming_prefix:str = None):
        """Get parameters for small molecules: ligands, cofactors, ...

        Parameters
        ----------
        mol_file : PathLike
            The path where the molecule is
        name : str, optional
            Name to give, by default "MOL"
        safe_naming_prefix : str, optional
            This is used to be sure that there will not happen any naming conflict in hte topologies, by default None

        Returns
        -------
        object
            The BioSimSpace system
        """

        if mol_file:
            print(f'Processing {mol_file}')
        else:
            return None
        
        parameterizer = Parameterize(
            force_field_code = 'openff_unconstrained-2.0.0.offxml',
            ext_types = ['top', 'gro'],
            hmr_factor = self.hmr_factor,
            overwrite = True,
            safe_naming_prefix = safe_naming_prefix,
            out_dir = self.wd,
        )
        # Actually you can pass to parameterize Chem.rdchem.Mol, *.inchi, *.smi, *.mol, *.mol2
        # An .sdf file is usually the same as .mol
        parameterizer(input_mol = mol_file,mol_resi_name = name)

        top_file = os.path.join(self.wd, f"{name}.top")
        gro_file = os.path.join(self.wd, f"{name}.gro")
        bss_system = bss.IO.readMolecules([top_file, gro_file])
        return bss_system

    def gmx_process(self, pdb_file:PathLike, pH:float = 7.0, is_membrane = False):
        """Used to process those biomolecules compatibles with amber99sb-ildn (protein, DNA, ..)
        and membrane compatibles with Slipid2020

        Parameters
        ----------
        pdb_file : PathLike
            Path to the file
        pH : float, optional
            pH for protonation (not working at the moment), by default 7.0
        is_membrane : bool, optional
            If True, Slipid2020 will be used instead of amber99sb-ildn, by default False

        Returns
        -------
        object
            The BioSimSpace system
        """

        if pdb_file:
            print(f'Processing {pdb_file}')
        else:
            return None
        
        name, _ = os.path.splitext(pdb_file)
        fixed_pdb = os.path.join(self.wd,f"{name}_fixed.pdb")

        # TODO, chake what is going on wrong, and use this kind od code, much better that call from the command line.
        # fixer = PDBFixer(filename=pdb_file)
        # if not is_membrane:
        #     fixer.findMissingResidues()
        #     fixer.findNonstandardResidues()
        #     fixer.replaceNonstandardResidues()
        #     fixer.removeHeterogens(True)
        # fixer.findMissingAtoms()
        # if not is_membrane:
        #     fixer.addMissingAtoms()
        # fixer.addMissingHydrogens(pH)
        # with open(fixed_pdb, 'w') as f:
        #     PDBFile.writeFile(fixer.topology, fixer.positions, f)

        env_prefix = os.environ["CONDA_PREFIX"]
        run(f"{env_prefix}/bin/pdbfixer {pdb_file} --output={fixed_pdb} --add-atoms=all --replace-nonstandard")

        gro_out = os.path.join(self.wd, f'{name}.gro')
        top_out = os.path.join(self.wd, f'{name}.top')
        posre_out = os.path.join(self.wd, f'{name}_posre.itp')

        if is_membrane:
            get_gmx_ff('Slipids_2020', out_dir=self.wd)
            run(f"gmx pdb2gmx -f {fixed_pdb} -ff Slipids_2020 -water none -o {top_out} -p {top_out} -i {posre_out}")
        else:
            run(f"gmx pdb2gmx -f {fixed_pdb} -merge all -ff amber99sb-ildn -water tip3p -o {gro_out} -p {top_out} -i {posre_out} -ignh")

        bss_system = bss.IO.readMolecules([gro_out,top_out])

        if self.hmr_factor:
            bss_system.repartitionHydrogenMass(factor=self.hmr_factor, water="no")

        cwd = os.getcwd()
        os.chdir(self.wd)
        bss.IO.saveMolecules(f"{name}_final", bss_system, ["GroTop"])
        bss_system = bss.IO.readMolecules([f'{name}.gro', f"{name}_final.top"])
        os.chdir(cwd)
        return bss_system

    def make_bss_system(self, ligand_mol:PathLike):
        """Create self.sys_ligand, self.sys_cofactor, self.sys_protein, self.sys_membrane
        and self.md_system (the combination of the available components). In case
        that the class was already called, it will be assumed that self.sys_cofactor, self.sys_protein, self.sys_membrane
        ere already calculated, only self.sys_ligand will be updated as well self.md_system

        Parameters
        ----------
        ligand_mol : PathLike
            The path of the where the ligand mol file
        """
        print("Processing components of the system...")
        self.sys_ligand = self.openff_process(
            mol_file = ligand_mol,
            name="LIG",
            safe_naming_prefix='w')
        
        # Only if the class has not yet called the full build will be carry out.
        if self.__self_was_called:
            print(f"Reusing components from cache")
        else:
            self.sys_cofactor = self.openff_process(
                mol_file = self.cofactor_mol,
                name="COF",
                safe_naming_prefix='z')
            self.sys_protein = self.gmx_process(pdb_file = self.protein_pdb)
            self.sys_membrane = self.gmx_process(pdb_file = self.membrane_pdb, is_membrane = True)
        print("Merging Components...")
        self.md_system = system_combiner(protein=self.sys_protein, membrane=self.sys_membrane, ligand=self.sys_ligand, cofactor=self.sys_cofactor)

    def clean(self):
        """Small cleaner, if keep_tmp_files_on was not provided,
        the intermediate steps saved on out_dir/.builder will be deleted
        """
        if not self.keep_tmp_files_on:
            shutil.rmtree(self.wd)

    def __call__(self, ligand_mol:PathLike, out_dir = 'abfe'):
        """The call implementation. It identify if it is needed to build
        all the components of the systems,
        In case that the class was already called, it will assume that all the components of the sytem,
        with the exception of the ligand, were already builded. This is usefull to call the class
        on several ligands that share the same components: protein, membrane and cofactor

        Parameters
        ----------
        ligand_mol : PathLike
            The path were the ligand is located
        out_dir : str, optional
            Were you would like to export the created files, by default 'abfe'
        """
        print(f"Processing ligand: {ligand_mol}")
        # Update (on multiple calls) or just create the out_dir (first call)
        self.out_dir = out_dir
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        
        # Construct MD system:
        self.make_bss_system(ligand_mol)
        system_dir = os.path.join(self.wd, 'system')
        ligand_dir = os.path.join(self.wd, 'ligand')

        # Solvate the system
        solvate(self.md_system, out_dir=system_dir, box = self.box, angles = self.box)
        solvate(self.sys_ligand, outdir=ligand_dir)

        # TODO Check what is done here and put it inside prepare_for_abfe
        fix_topology(input_topology_path=os.path.join(system_dir,'solvated.top'), out_topology_path=os.path.join(system_dir,'solvated_fix.top'))
        fix_topology(input_topology_path=os.path.join(ligand_dir,'solvated.top'), out_topology_path=os.path.join(ligand_dir,'solvated_fix.top'))
        
        # TODO Check what is done here
        # Construct ABFE system:
        prepare_for_abfe(out_dir=self.out_dir, ligand_dir=ligand_dir, sys_dir=system_dir)
        
        self.clean()

        # Change state
        self.__self_was_called = True


#############################################################################################

if __name__ == "__main__":...
    # box = [116.982, 116.982, 105.354]
    # angles = [90.0, 90.0, 60.0]
    # builder = MakeInputs(
    #         protein_pdb = 'protein.pdb', #None,#'protein.pdb',
    #         membrane_pdb = None,#'membrane.pdb',#'membrane.pdb',
    #         cofactor_mol = 'ligand.mol',# None
    #         hmr_factor = 3,
    #         keep_tmp_files_on = None)
    # builder(ligand_mol = 'ligand1.mol',out_dir='abfe/lig1')
    # builder(ligand_mol = 'ligand2.mol',out_dir='abfe/lig2')
    # builder(ligand_mol = 'ligand3.mol',out_dir='abfe/lig3')