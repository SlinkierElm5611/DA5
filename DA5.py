# Imported libraries
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import scipy
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
Ly = %s! length (0.5mm - 1 cm)
LxTotal = %s*Ly ! total thickness 
LyAu = %s*Ly ! window length 


Lz = Ly/4 ! width

m = 1.163451761 !0.8595113554
{LxQz = LxTotal/(m+1) ! thickness of quartz
LxTi = LxTotal - LxQz ! thickness of titanium}
LxTi = LxTotal/(m+1) ! thickness of quartz
LxQz = LxTotal - LxTi ! thickness of titanium
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
Efieldy = 0
Efieldz = 0

!! Strain
!Axial Strain
ex=dx(u)
ey=dy(v)
!Engineering Shear Strain
gxy=(dx(v)+dy(u))
!Mechanical strain
exm  = ex  - (d11*Efieldx+d21*Efieldy+d31*Efieldz)
eym  = ey  - (d12*Efieldx+d22*Efieldy+d32*Efieldz)
gxym = gxy - (d16*Efieldx+d26*Efieldy+d36*Efieldz)

!!Stress via Hooke's law
!Axial Stress
sx = C11*ex+C12*ey+C16*gxy
sy = C21*ex+C22*ey+C26*gxy

!Shear stress
sxy=C61*ex+C62*ey+C66*gxy

 f = sqrt(lambda)/2*pi

EQUATIONS
!FNet = 0
u:	dx(sx)+dy(sxy)=-rho*lambda^2*u 
v:	dx(sxy)+dy(sy)=-rho*lambda^2*v

BOUNDARIES
  	REGION 1 'Titanium'
	    rho = 4.506e3
		E = 106e9
        nu = 0.34
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
        
    REGION 3 'Gold'
    	rho = 19.3e3
        E = 79e9
		nu = 0.42
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
 	!contour(u) on surface z=0
    elevation(u) from (LxTi+LxQz/2,-Ly/2) to (LxTi+LxQz/2, Ly/2)
	
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
last_freq = 0
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
    dx = np.array([delta, 0, 0,0])
    dy = np.array([0, delta, 0,0])
    dz = np.array([0, 0, delta,0])
    dmat = np.array([0, 0, 0,0])
    print(dim)
    
    freq0 = getSensitivity(dim)
    #rel_sens = ((last_freq - freq0)/last_freq)*100
    last_freq = freq0

    print(freq0)
    
    print(dim + dx)
    print(dim + dy)
    print(dim+ dz)
    gradient = np.array([
        (getSensitivity(dim + dx) - freq0) / delta,
        (getSensitivity(dim + dy) - freq0) / delta,
        (getSensitivity(dim + dz) - freq0) / delta,
        (getSensitivity(dim + dmat)-freq0)/delta
    ])
    print(gradient)
    return gradient
    
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
        print(delta)
        dimNext = dim0 + delta*g
        #print("Grad:", g) #debugging
        print("The frequency is " + str(Freq0) + " at " + str(dim0))
        freqNext = getSensitivity(dim0)
        rel_sens = abs(((Freq0 - freqNext)/Freq0))*100
        if (freqNext < Freq0):
            print("Passed maximum.")
            return dim0, Freq0
        else:
            dim0, Freq0 = dimNext, freqNext
    return dim0, Freq0

def gridSearch():
    Lbridge = np.linspace(5e-5,0.01,10)
    hper = np.linspace(0.05,0.2, 4)
    lper = np.linspace(0.05,0.4, 5)
    last_freq = 0
    sensitivity = {}
    rel_sens = 0
    mat_th = 0 #Test value for now to find optimal dimensions
    
    for length in Lbridge:
        for h in hper:
            for l in lper:
                if h*(length) <= 0.001 or h*(length) >= 1e-5:
                    hbridge = h*length
                    Lwindow = l*(length)
                    Wbridge = length/4
                    #sensivity[getsensevity(Lbridge, Wbridge, hbridge,Lwindow)] = [Lbridge, Wbridge, hbridge,Lwindow]
                    freq = getSensitivity([length,h,l, mat_th])
                    if last_freq == 0:
                        last_freq = freq
                    else:
                        rel_sens = (abs(((last_freq - freq)/last_freq))*100) #Percent of change in relative frequency 
                    sensitivity[rel_sens] = [length, h, l, mat_th]
                    print("The freq is "+ str(rel_sens) + " at " +  str(length) + " m length " + str(Wbridge)+ " m width "+ 
                          str(hbridge) + " m height " + str(Lwindow) + " m window \n")

    max_sens = max(sensitivity.keys())

    print("The maximum sensitivity is "+ str(max_sens)+ "at" + str(sensitivity[max_sens]))
    
    return [sensitivity[max_sens], max_sens, max()]

def mat_thickness(dimensions, freq, properties):
    max_mat_thickness = np.linspace(0, 2e-5, 20) #guess for max material change
    Lbridge = dimensions[0]
    hpercent = dimensions[1]
    lpercent = dimensions[2]
    mat_th = dimensions[3]
    for th in max_mat_thickness:
        freq = getSensitivity([Lbridge, hpercent, lpercent, th])
        if freq == 0:
            break
        mat_th = th
    return mat_th

def main():
    max_sens = [[0.01, 0.05, 0.3125, 0],99.50004591047193]

    #Initial sweep of all values in ranges with large step sizes
    max_sens = gridSearch()

    #Further scope into a zoomed in section of the maximum value
    max_freq = 0
    print("Starting gradient descent with" + str(max_sens))
    opt_dim = []
    #opt_dim, max_freq = gradientAscent(max_sens[0], max_sens[1])

    #Finding the longevity values for various material aka when it stops working
   # gold = mat_thickness(opt_dim, max_freq, ) #add material value
   # molybdenum = mat_thickness(opt_dim, max_freq, ) #add material value
   # chromium = mat_thickness(opt_dim, max_freq, ) #add material value
   # SiO2 =  mat_thickness(opt_dim, max_freq, ) #add material value
main()