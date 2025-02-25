import code_test as ct
import numpy as np

xpath = '../Linelists/'
xlist = 'data/gaslistR10K.dat'
results = [ct.NoCloud_Tdwarf(xpath,xlist),ct.MieClouds_Ldwarf(xpath,xlist)]

if np.all(results):
    print('     ')
    print('------------------------------------------------------------')
    print("Everything is as it was. No >1% differences with benchmark cases")

else:
    print('     ')
    print('------------------------------------------------------------')
    print("> 1% DIFFERENCES PRESENT BETWEEN NEW OUTPUTS AND BENCHMARKS.")
    print("IS THIS WHAT YOU WANTED??  ")
    print("PLEASE CHECK THAT YOUR CHANGES ARE HAVING THE DESIRED EFFECT.")
    print("IF SO, PLEASE PRODUCE NEW BENCHMARKS AND UPDATE THE TESTS ACCORDINGLY")

    
