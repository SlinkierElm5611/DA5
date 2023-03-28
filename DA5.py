# Imported libraries
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import scipy
FlexCode = """  """
#if polymer plate on bottom make this True:
poly_bot = False
# Path Variables
FlexFileName = "DA5" 
FlexLocation = "FlexPDE6s"
FlexVersion = 6
TextFile = (FlexVersion==7)*("DA5_output/")+'DA5_summary.txt'

def getSensitivity(Lbridge, hpercent, lpercent):
    disp = 0
    with open(FlexFileName + '.pde', 'w') as f: 
        #print(FlexCode %(Lbridge, Wbridge, hbridge,Lwindow), file = f)
        print(FlexCode %(Lbridge, hpercent, lpercent), file = f)
    # Runs FlexPDE
    subprocess.run([FlexLocation, "-S",FlexFileName],timeout=5)

    # Retrieves data
    with open(TextFile) as f:
        text = f.readlines()[-1].strip().split(" ")
        #print(text[1])
        disp = float(text[-1].strip())
        return disp
      
def grad(dimensions, delta):
    dx = [delta, 0]
    dy = [0,delta, 0]
    dz = [0,0,delta]
    
def gradientAscent(dimensions):
    #Similar to the gradient ascent algorithm but goes towards a maximum instead by moving towards a positive direction
    delta = 0.00001 
    Lbridge = dimensions[0][0] # initial guess
    hbridge = dimensions[0][1]*(Lbridge)
    Lwindow = dimensions[0][2]*(Lbridge)
    sensitivity = dimensions[1]
    Wbridge = Lbridge/4
    for i in range(20):
        g = grad(dimensions, delta) 
        thNext = th0 + delta*g
        #print("Grad:", g) #debugging
        print("The displacement is " + str(Disp0) + " at " + str(th0))
        DispNext = getSensitivity(th0)
        if (DispNext < Disp0):
            print("Passed maximum.")
            return th0, Disp0
        else:
            th0, Disp0 = thNext, DispNext
    
    return th0, Disp0

def gridSearch():
    Lbridge = np.linspace(5e-5,0.01,5)
    hper = np.linspace(0.05,0.2, 4)
    lper = np.linspace(0.05,0.4, 5)

    sensitivity = {}
    rel_sens = 0
    last_freq = 0
    for length in Lbridge:
        for h in hper:
            for l in lper:
                if h*(length) <= 0.001 or h*(length) >= 1e-5:
                    hbridge = h*length
                    Lwindow = l*(length)
                    Wbridge = length/4
                    #sensevity[getsensevity(Lbridge, Wbridge, hbridge,Lwindow)] = [Lbridge, Wbridge, hbridge,Lwindow]
                    freq = getSensitivity(length,h,l)
                    if last_freq == 0:
                        last_freq = freq
                    else:
                        rel_sens = ((last_freq - freq)/last_freq)*100 #Percent of change in relative frequency 
                    sensitivity[rel_sens] = [length, h, l]
    max_sens = max(sensitivity.keys())
    return [sensitivity[max_sens], max_sens]
def mat_thickness():
    return
def main():
    max_sens = []

    max_sens = gridSearch()

    next_max = []
    next_max = gradientAscent(max_sens)

main()