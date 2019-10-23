# brewster
code for retrieval analysis of emission spectra from brown dwarfs and giant planets

## What's needed to run the code:
- Everything in this repository
- Atmospheric Models
- Line lists 
- Clouds (.mieff files) (This must be at the same level as the brewster directory)
- gFortran compiler (https://gcc.gnu.org/wiki/GFortranBinaries or https://hpc.sourceforge.net)
 
 The Atmospheric Models can live anywhere on your machine.

## Current Process for Installing Brewster
1. Get access to the Brewster Github. Either send Ben or Eileen your GitHub name or email associated with your GitHub account to get added to the Substellar team.
2. If you do not have a fortran compiler on your machine, get the proper version of the gfortran compiler from https://gcc.gnu.org/wiki/GFortranBinaries or https://hpc.sourceforge.net.
3. Git clone or fork the directory and place wherever your code lives. (If you plan to make any changes to the code you might like to incorporate into the master branch via a pull request, you should probably fork it.)
4. To keep various package versions clean and organized, you should create a Python 3 environment, ideally Python 3.7, as that is the version that will be needed to run the code on a cluster. (See list of packages needed in requirements.txt)
5. You will need to get access to the Linelist and Clouds files that are shared via Dropbox. Place the Clouds and Linelists folders separately at the same level as the Brewster code. i.e: 
```bash
|--Folder_Where_all_your_code_lives
   |-- Brewster
   |-- Clouds
   |-- Linelists
 ```
6. Build the code by running ``` ./build``` in the terminal in the Brewster directory. If you get an error or this is not a fresh build, run make clean before rebuilding the code. This will create ciamod and forwardmodel (which show up as missing imports in brewster.py and testkit.py), as well as various .o, .mod, and .so files.
7. To test that everything is properly installed, run the check_brewster.py file, by typing ```python check_brewster.py``` in the terminal making sure you are in the brewster working directory.
Common problems so far have included:
    1. Name of water linelist pickle should be renamed from h2o_xsecs_R10K.pic to  h2o_ucl2017_xsecs_R10K.pic
    2. Any others we ran into?
