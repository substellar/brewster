# brewster
code for retrieval analysis of emission spectra from brown dwarfs and giant planets

## What's needed to run the code:
- Everything in the Repo
- Atmospheric Models
- Line lists 
- Clouds (.mieff files) (This must be at the same level as the brewster directory)
- gFortran compiler (https://gcc.gnu.org/wiki/GFortranBinaries)
 
 The line list and Atmospheric Models can live anywhere on your machine.

## Getting Started
1. Run make clean in the terminal. This makes sure that nothing is left over from the previous build.
2. Build the code using ```./build```. This will create ciamod and forwardmodel (which show up as missing imports in 
brewster.py and testkit.py), as well as various .o, .mod, and .so files. 
3. Test that everything looks and is acting properly for your object by running the make_a_spectrum.ipynb
