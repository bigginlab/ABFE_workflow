import os, glob
import tqdm
import shutil
from time import sleep

orig_dir = os.getcwd()

#Paths! <-- Change!
in_root_dir = "/path/to/CyclophilinD"
out_root_dir = "/path/to/data/input/CyclophilinD"


print("move and prepare files:")
print("\tinput: ",in_root_dir)
print("\toutput: ",out_root_dir)

if(not os.path.isdir(out_root_dir) and os.path.isdir(os.path.dirname(out_root_dir))):
    print("\tgenerating out_root_dir")
    os.mkdir(out_root_dir)
print()

#copy files
for tar_file in tqdm.tqdm(glob.glob(in_root_dir+"/*tar.gz"), desc="copying", leave=True):
    out_name = out_root_dir+"/"+os.path.basename(tar_file)
    shutil.copy(tar_file, out_name)

#Untar
os.chdir(out_root_dir)
for tar_file in tqdm.tqdm(os.listdir(), desc="untar", leave=True):
    if("gz" in tar_file):
       os.system("tar -xzf "+tar_file)
    else:
       os.system("tar -xf "+tar_file)
        

#clean folder
iterate = list(os.walk(out_root_dir))
for path, dirs, files in tqdm.tqdm(iterate, desc="rename", leave=True):

    if("gromacs" in dirs):
        if(path.endswith("complex")):
            os.system("mv "+path+"/gromacs/complex*.gro "+path+"/complex.gro")
            os.system("mv "+path+"/gromacs/complex*.top "+path+"/complex.top")
            os.system("mv "+path+"/gromacs/*.itp "+path+"/")

        if(path.endswith("ligand")):
            os.system("mv "+path+"/gromacs/ligand*.gro "+path+"/ligand.gro")
            os.system("mv "+path+"/gromacs/ligand*.top "+path+"/ligand.top")
            os.system("mv "+path+"/gromacs/*.itp "+path+"/")

        for dir in dirs:        #cleanout all the rest:
            os.system("rm -r "+path+"/"+dir)

os.chdir(orig_dir)
print("\nDone")