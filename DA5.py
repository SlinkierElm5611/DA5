# Imported libraries
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import scipy
FlexCode = """TITLE 
'DA5'
SELECT
errlim=1e-3
ngrid=1
spectral_colors
modes = 1

COORDINATES
cartesian2

VARIABLES
u	!Displacement in x
v !Displacement in y
!w !Displacement in z

DEFINITIONS
mag=0.3*globalmax(magnitude(x,y))/globalmax(magnitude(u,v))
Ly = %s ! length
Lx = %s*(Ly)
LxQz = Lx/(0.85+1) ! thickness
LxTi = Lx-LxQz ! thickness h/(m+1) = a2
LyAu = %s*Ly ! window length
LxAu = %s ! window thickness
Lz = Ly/4 ! width

! Material Properties
rho
E = 0
nu = 0
G=E/(2*(1+nu))

C11 =E*(1-nu)/(1+nu)/(1-2*nu)
C12 = E*nu/(1+nu)/(1-2*nu)	C13 = C12	C14=0	C15=0	C16=0
C21 = C12	C22 = C11	C23 = C12	C24=0 C25=0 C26=0
C31 = C12	C32 = C12 	C33 = C11	C34=0 C35=0 C36=0
C41=0 C42=0 C43=0 C44=G C45=0 C46=0
C51=0 C52=0 C53=0 C54=0 C55=G C56=0
C61=0 C62=0 C63=0 C64=0 C65=0 C66=G

!! Strain
!Axial Strain
ex=dx(u) !partial derivatives
ey=dy(v)
!Engineering Shear Strain
gxy=(dx(v)+dy(u))

!!Stress via Hooke's law
!Axial Stress
sx = C11*ex+C12*ey+C16*gxy
sy = C21*ex+C22*ey+C26*gxy

!Shear stress
sxy=C61*ex+C62*ey+C66*gxy

omega = sqrt(lambda)
f = 2*pi*omega

EQUATIONS
!FNet = 0
u:	dx(sx)+dy(sxy)=-rho*lambda^2*u 
v:	dx(sxy)+dy(sy)=-rho*lambda^2*v

BOUNDARIES
  	REGION 1 'Titanium'
	    rho = 4.506e3
		E = 106e9
        nu = 0.34
    	START(0,-Ly/2) !y=0 surface:
		value(u)=0
		value(v)=0
    	LINE TO (LxTi,-Ly/2) !x=Lx surface
		load(u)=0
		load(v)=0
		LINE TO (LxTi,Ly/2) !y=Ly surface
		value(u)=0
		value(v)=0
		LINE TO (0,Ly/2) !x=0 surface
		load(u)=0
		load(v)=0
		LINE TO CLOSE
        
    REGION 2 'Quartz'
    	rho = 2.684
		C11 = 8.67361769904436e10	C12 = 6.98526725701223e9	C13 = 1.19104335397808e10	C14 = 1.79081384131957e10	C15=0	C16=0
		C21 = 6.98526725701223e9	C22 = 8.67361769904436e10	C23 = 1.19104335397808e10	C24=-1.79081384131957e10	C25=0 C26=0
		C31 = 1.19104335397808e10	C32 = 1.19104335397808e10	C33 = 1.07193901858028e11	C34=0 C35=0 C36=0
		C41=1.79081384131957e10	C42 = -1.79081384131957e10	C43=-2.57028375615534e-7 C44=5.79427767324731e10 C45=0 C46=0
		C51=0 C52=0 C53=0 C54=0 C55=5.79491958802304e10 C56=1.79224317155352e10
		C61=0 C62=0 C63=0 C64=0 C65=1.79224317155352e10 C66=3.99072812865916e10
    	START(LxTi,-Ly/2) !y=0 surface:
		value(u)=0
		value(v)=0
    	LINE TO (LxTi+LxQz,-Ly/2) !x=Lx surface
		load(u)=0
		load(v)=0
		LINE TO (LxTi+LxQz,Ly/2) !y=Ly surface
		value(u)=0
		value(v)=0
		LINE TO (LxTi,Ly/2) !x=0 surface
		load(u)=0
		load(v)=0
		LINE TO CLOSE
        
    REGION 3 'Gold'
    	rho = 19.3e3
        E = 79e9
		nu = 0.42
    	START(LxTi+LxQz,-LyAu/2) !y=0 surface:
		value(u)=0
		value(v)=0 
    	LINE TO (LxTi+LxQz+LxAu,-LyAu/2) !x=Lx surface
		load(u)=0
		load(v)=0
		LINE TO (LxTi+LxQz+LxAu,LyAu/2) !y=Ly surface
		value(u)=0
		value(v)=0
		LINE TO (LxTi+LxQz,LyAu/2) !x=0 surface
		load(u)=0
		load(v)=0
		LINE TO CLOSE


PLOTS
	grid(x+u*mag, y+v*mag)
 	!contour(u) on surface z=0
    !contour(sy) on surface y=0
	!elevation(sx,sy,sz) from (0,0,0) to (0,0,Lz)
    elevation(u) from (LxTi+LxQz/2,-Ly/2) to (LxTi+LxQz/2, Ly/2)
	
	SUMMARY
		export file 'DA4_summary.txt'
		report f
end"""
# Path Variables
FlexFileName = "DA5" 
FlexLocation = "FlexPDE6s"
FlexVersion = 6
TextFile = (FlexVersion==7)*("DA5_output/")+'DA5_summary.txt'

def getSensitivity(dimensions):
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
      
def grad(dim, delta):
    dx = np.array([delta, 0, 0])
    dy = np.array([0, delta, 0])
    dz = np.array([0, 0, delta])
    freq0 = getSensitivity(dim)
    grad = np.array([
        (getSensitivity(dim + dx) - freq0) / delta,
        (getSensitivity(dim + dy) - freq0) / delta,
        (getSensitivity(dim + dz) - freq0) / delta,
    ])
    return grad
    
def gradientAscent(dimensions, freqGuess):
    #Similar to the gradient ascent algorithm but goes towards a maximum instead by moving towards a positive direction
    delta = 0.00001 
    ''' Lbridge = dimensions[0][0] # initial guess
    hbridge = dimensions[0][1]*(Lbridge)
    Lwindow = dimensions[0][2]*(Lbridge)
    sensitivity = dimensions[1]
    Wbridge = Lbridge/4'''
    dim0 = np.array(dimensions)
    Freq0 = freqGuess
    for i in range(20):
        g = grad(dim0, delta) 
        dimNext = dim0 + delta*g
        #print("Grad:", g) #debugging
        print("The frequency is " + str(Freq0) + " at " + str(dim0))
        freqNext = getSensitivity(dim0)
        if (dimNext < dim0):
            print("Passed maximum.")
            return dim0, Freq0
        else:
            dim0, Freq0 = dimNext, freqNext
    return dim0, Freq0

def gridSearch():
    Lbridge = np.linspace(5e-5,0.01,5)
    hper = np.linspace(0.05,0.2, 4)
    lper = np.linspace(0.05,0.4, 5)

    sensitivity = {}
    rel_sens = 0
    last_freq = 0
    mat_th = 0 #Test value for now to find optimal dimensions
    for length in Lbridge:
        for h in hper:
            for l in lper:
                if h*(length) <= 0.001 or h*(length) >= 1e-5:
                    hbridge = h*length
                    Lwindow = l*(length)
                    Wbridge = length/4
                    #sensevity[getsensevity(Lbridge, Wbridge, hbridge,Lwindow)] = [Lbridge, Wbridge, hbridge,Lwindow]
                    freq = getSensitivity([length,h,l, mat_th])
                    if last_freq == 0:
                        last_freq = freq
                    else:
                        rel_sens = ((last_freq - freq)/last_freq)*100 #Percent of change in relative frequency 
                    sensitivity[rel_sens] = [length, h, l, mat_th]
                    print("The frequency and dimensions are as follows"+ str(sensitivity))
    max_sens = max(sensitivity.keys())
    return [sensitivity[max_sens], max_sens]

def mat_thickness(dimensions, freq, properties):
    max_mat_thickness = np.linspace(0, 2e-5, 20) #guess for max material change
    Lbridge = dimensions[0]
    hpercent = dimensions[1]
    lpercent = dimensions[2]
    mat_th = dimensions[3]

def main():
    max_sens = []

    #Initial sweep of all values in ranges with large step sizes
    max_sens = gridSearch()

    #Further scope into a zoomed in section of the maximum value
    max_freq = 0
    opt_dim = []
    opt_dim, max_freq = gradientAscent(max_sens)

    #Finding the longevity values for various material aka when it stops working
    gold = mat_thickness(opt_dim, max_freq, ) #add material value
    molybdenum = mat_thickness(opt_dim, max_freq, ) #add material value
    chromium = mat_thickness(opt_dim, max_freq, ) #add material value
    SiO2 =  mat_thickness(opt_dim, max_freq, ) #add material value
main()