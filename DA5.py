# Imported libraries
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import differential_evolution
class CustomStopException(Exception):
    pass
FlexCode = """
TITLE 
'DA5'
SELECT
errlim=1e-3
ngrid=5
spectral_colors
modes = 1

COORDINATES
cartesian2

VARIABLES
u	!Displacement in x
v !Displacement in y

DEFINITIONS
mag=0.3*globalmax(magnitude(x,y))/globalmax(magnitude(u,v))

! Dimensions
Ly = %s ! length (0.5mm - 1 cm)
Lz = Ly/4 ! width

LxTotal = %s*Ly ! total thickness 
m = 1.163451761
LxTi = LxTotal/(m+1) ! thickness of titanium
LxQz = LxTotal - LxTi ! thickness of quartz

LyAu = %s*Ly ! window length 
LxAu = %s ! window thickness

! Material Properties
rho
E = 0
nu = 0
G=E/(2*(1+nu))

{ Stiffness matrix }
C11 =E*(1-nu)/(1+nu)/(1-2*nu)
C12 = E*nu/(1+nu)/(1-2*nu)	C13 = C12	C14=0	C15=0	C16=0
C21 = C12	C22 = C11	C23 = C12	C24=0 C25=0 C26=0
C31 = C12	C32 = C12 	C33 = C11	C34=0 C35=0 C36=0
C41=0 C42=0 C43=0 C44=G C45=0 C46=0
C51=0 C52=0 C53=0 C54=0 C55=G C56=0
C61=0 C62=0 C63=0 C64=0 C65=0 C66=G

{ Piezoelectric coupling coefficients }
d11 = 20e-12            d12 = -24e-12		d13 = 0             d14 = 0 d15 = 0 d16 = 0
d21 = 16e-12            d22 = -35e-12		d23 = 0             d24 = 0 d25 = 0 d26 = 0
d31 = 45e-12            d32 = 0				d33 = 70e-12    d34 = 0 d35 = 0 d36 = 0

{ Electric field }
voltage = 1e3
Efieldx = -voltage/LxQz

!! Strain
!Axial Strain
ex=dx(u)
ey=dy(v)
!Engineering Shear Strain
gxy=(dx(v)+dy(u))
!Mechanical strain
exm  = ex  - d11*Efieldx
eym  = ey  - d12*Efieldx
gxym = gxy - d16*Efieldx

!!Stress via Hooke's law
!Axial Stress
sx = C11*exm+C12*eym+C16*gxym
sy = C21*exm+C22*eym+C26*gxym

!Shear stress
sxy=C61*exm+C62*eym+C66*gxym

 f = sqrt(lambda)/(2*pi)

EQUATIONS
!FNet = 0
u:	dx(sx)+dy(sxy)=-rho*lambda^2*u 
v:	dx(sxy)+dy(sy)=-rho*lambda^2*v

BOUNDARIES
  	REGION 1 'Titanium'
	    rho = 4.506e3
		E = 106e9
        nu = 0.34
        d11 = 0           d12 = 0		d13 = 0    d14 = 0 d15 = 0 d16 = 0
		d21 = 0           d22 = 0		d23 = 0    d24 = 0 d25 = 0 d26 = 0
		d31 = 0           d32 = 0		d33 = 0	d34 = 0 d35 = 0 d36 = 0
    	START(0,0) !y=0 surface:
		load(u)=0
		value(v)=0
    	LINE TO (-LxTi,0) !x=Lx surface
		load(u)=0
		load(v)=0
		LINE TO (-LxTi,Ly/2) !y=Ly surface
		value(u)=0
		value(v)=0 
		LINE TO (0,Ly/2) !x=0 surface
		load(u)=0
		load(v)=0
		LINE TO CLOSE
        
    REGION 2 'Quartz'
    	rho = 2.684e3
		C11 = 8.67361769904436e10	C12 = 6.98526725701223e9	C13 = 1.19104335397808e10	C14 = 1.79081384131957e10	C15=0	C16=0
		C21 = 6.98526725701223e9	C22 = 8.67361769904436e10	C23 = 1.19104335397808e10	C24=-1.79081384131957e10	C25=0 C26=0
		C31 = 1.19104335397808e10	C32 = 1.19104335397808e10	C33 = 1.07193901858028e11	C34=0 C35=0 C36=0
		C41=1.79081384131957e10	C42 = -1.79081384131957e10	C43=-2.57028375615534e-7 C44=5.79427767324731e10 C45=0 C46=0
		C51=0 C52=0 C53=0 C54=0 C55=5.79491958802304e10 C56=1.79224317155352e10
		C61=0 C62=0 C63=0 C64=0 C65=1.79224317155352e10 C66=3.99072812865916e10
    	START(0,0) !y=0 surface:
		load(u)=0
		value(v)=0
    	LINE TO (LxQz,0) !x=Lx surface
		load(u)=0
		load(v)=0
		LINE TO (LxQz,Ly/2) !y=Ly surface
        value(u)=0
		value(v)=0 
		LINE TO (0,Ly/2) !x=0 surface
		load(u)=0
		load(v)=0
		LINE TO CLOSE
        
    REGION 3 'Deposit'
    	! Gold
    	rho = 19.3e3
        E = 79e9
		nu = 0.42
        d11 = 0           d12 = 0		d13 = 0    d14 = 0 d15 = 0 d16 = 0
		d21 = 0           d22 = 0		d23 = 0    d24 = 0 d25 = 0 d26 = 0
		d31 = 0           d32 = 0		d33 = 0	d34 = 0 d35 = 0 d36 = 0
    	START(LxQz,0) !y=0 surface:
		load(u)=0
		value(v)=0 
    	LINE TO (LxQz+LxAu,0) !x=Lx surface
		load(u)=0
		load(v)=0
		LINE TO (LxQz+LxAu,LyAu/2) !y=Ly surface
		value(u)=0
		value(v)=0 
		LINE TO (LxQz,LyAu/2) !x=0 surface
		load(u)=0
		load(v)=0
		LINE TO CLOSE


PLOTS
	grid(x+u*mag, y+v*mag)
 	!contour(u) painted
    elevation(u) from (0,-Ly/2) to (0, Ly/2)
	
	SUMMARY
        export file 'DA5_summary.txt'
		report f
end
"""
# Path Variables
FlexFileName = "DA5" 
FlexLocation = "FlexPDE6s"
FlexVersion = 6
TextFile = (FlexVersion==7)*("DA5_output/")+'DA5_summary_0.txt'
max_val = ""
def getFrequency(dimensions):
    freq = 0
    Lbridge = dimensions[0]
    hpercent = dimensions[1]
    lpercent = dimensions[2]
    mat_th = dimensions[3]
    with open(FlexFileName + '.pde', 'w') as f: 
        #print(FlexCode %(Lbridge, Wbridge, hbridge,Lwindow), file = f)
        print(FlexCode %(Lbridge, hpercent, lpercent, mat_th ), file = f)
    # Runs FlexPDE
    subprocess.run([FlexLocation, "-S",FlexFileName],timeout=5)

    # Retrieves data
    with open(TextFile) as f:
        text = f.readlines()[-1].strip().split(" ")
        #print(text[1])
        freq = float(text[-1].strip())
        return freq
    
def grad(dim, delta, freq0):
    dx = np.array([delta, 0, 0,0])
    dy = np.array([0, delta, 0,0])
    dz = np.array([0, 0, delta,0])
    dmat = np.array([0, 0, 0,0])
    print(dim)
    
    #rel_sens = ((last_freq - freq0)/last_freq)*100
    last_freq = getFrequency(dim)

    print(last_freq)
    print(getFrequency(dim + dx))
    gradient = np.array([
        (getFrequency(dim + dx) - last_freq) / delta,
        (getFrequency(dim + dy) - last_freq) / delta,
        (getFrequency(dim + dz) - last_freq) / delta,
        (getFrequency(dim + dmat)-last_freq)/delta
    ])
    print(gradient)
    return gradient
    
def gradientDescent(dimensions, sensGuess):
    #Similar to the gradient ascent algorithm but goes towards a maximum instead by moving towards a positive direction
    delta = 200000
    dim0 = np.array(dimensions)
    sens0 = sensGuess
    for i in range(5):
        last_freq = getFrequency(dim0)
        g = grad(dim0, delta, sens0) 
        dimNext = dim0 - delta*g
        print("Grad:", g) #debugging
        print(dimNext)
        freqNext = getFrequency(dimNext)
        sensNext = abs(((freqNext - last_freq )/last_freq))*100
        print("The sensitivity is " + str(sensNext) + " at " + str(dim0))
        print(freqNext)
        if (sensNext > sens0):
            print("Passed maximum.")
            return dim0, sens0
        else:
            dim0, sens0 = dimNext, sensNext
    return dim0, sens0

# Define a function to get the frequency based on the dimension and delta
def get_frequency(dim, delta):
    # Your implementation to calculate the frequency based on the dimension and delta
    return getFrequency(dim)
    

# Define the objective function to minimize the relative frequency

def objective_function(x):
    dim, delta = x[:-1], x[-1]  # Separate the dimensions and delta
    dx = delta+dim[0]
    dy = delta+dim[1]
    dz =delta+dim[2]
    dmat =0
    last_freq = getFrequency(dim)  # Get the frequency for the current dimensions
    new_dim = [dx, dy, dz, dmat]  # Update dimensions with delta
    new_freq = getFrequency(new_dim)  # Get the frequency for the updated dimensions
    threshold = 0.0002

    # Calculate the relative frequency
    relative_frequency = abs((new_freq - last_freq) / last_freq) * 100
    print(f'Dimensions: {dim},New Dimensions: {new_dim}, Delta: {delta}, Last Freq: {last_freq}, New Freq: {new_freq}, Relative Frequency: {relative_frequency}')
    
    return relative_frequency



def gridSearch():
    Lbridge = np.linspace(5e-4,0.01,5)
    hper = np.linspace(0.05,0.2, 4)
    lper = np.linspace(0.05,0.4, 5)
    last_freq = 0
    sensitivity = {}
    rel_sens = 0
    mat_th = 0 #Test value for now to find optimal dimensions
    
    for length in Lbridge:
        for h in hper:
            last_freq = 0
            for l in lper:
                if h*(length) <= 0.001 or h*(length) >= 1e-5:
                    hbridge = h*length
                    Lwindow = l*(length)
                    Wbridge = length/4
                    #sensivity[getsensevity(Lbridge, Wbridge, hbridge,Lwindow)] = [Lbridge, Wbridge, hbridge,Lwindow]
                    freq = getFrequency([length,h,l, mat_th])
                    if last_freq == 0:
                        last_freq = freq
                    else:
                        rel_sens = (abs(((last_freq - freq)/last_freq))*100) #Percent of change in relative frequency 
                        last_freq = freq
                    sensitivity[rel_sens] = [length, h, l, mat_th]
                    print("The sensitivity is "+ str(rel_sens) + " at " +  str(length) + " m length " + str(Wbridge)+ " m width "+ 
                          str(hbridge) + " m height " + str(Lwindow) + " m window \n")

    max_sens = min(sensitivity.keys())
    print("The maximum sensitivity is "+ str(max_sens)+ "at" + str(sensitivity[max_sens]))
    
    return [sensitivity[max_sens], max_sens]

def mat_thickness(dimensions, freq, properties):
    max_mat_thickness = np.linspace(0, 2e-5, 20) #guess for max material change
    Lbridge = dimensions[0]
    hpercent = dimensions[1]
    lpercent = dimensions[2]
    mat_th = dimensions[3]
    for th in max_mat_thickness:
        freq = getFrequency([Lbridge, hpercent, lpercent, th])
        if freq == 0:
            break
        mat_th = th
    return mat_th

def callback(xk, convergence):
    relative_frequency = objective_function(xk)
    if 0 < relative_frequency <= 0.0002:
        return True  # Stop the optimization
    return False  # Continue the optimization
def main():
    max_sens = []

    #Initial sweep of all values in ranges with large step sizes
    max_sens = gridSearch()

    #Further scope into a zoomed in section of the maximum value
    max_freq = 0
    print("Starting gradient descent with" + str(max_sens))
    opt_dim = []
    #freq = last_freq*(max_sens[1]/100)
    #opt_dim, max_freq = gradientDescent(max_sens[0], (max_sens[1] ))


    #Alternative method of optimization: Evolutionary Algorithm
    print("Starting differential evolution")
    bounds = [(5e-4, 0.01), (0.05, 0.2), (0.05, 0.4), (0, 0), (0, 0.00001)]

    # Use differential evolution to find the optimal dimensions and delta
    result = differential_evolution(objective_function, bounds, maxiter=1, callback=callback, disp=True)

    # Get the optimal dimensions, delta, and the corresponding relative frequency
    print(result)
    optimal_dimensions_and_delta = result.x
    optimal_relative_frequency = result.fun

    optimal_dimensions = optimal_dimensions_and_delta[:-1]
    optimal_delta = optimal_dimensions_and_delta[-1]

    print("Optimal dimensions:", optimal_dimensions)
    print("Optimal delta:", optimal_delta)
    print("Optimal relative frequency:", optimal_relative_frequency)

main()