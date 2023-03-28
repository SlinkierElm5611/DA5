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
LxAu = 0 ! window thickness
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
TextFile = 'DA5_summary.txt'
FlexVersion = 6


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