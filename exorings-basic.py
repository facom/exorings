from exorings import *
"""
This script calculate exoring transit basic properties.

Using this properties it additionally provide the observed planetary
radius and observed stellar density.

Input parameters:

    Ms: Stellar mass.
    Rs: Stellar radius.
    Rp: Planetary radius.
    a: Semimajor axis.
    iorb: Orbital inclination.
    fi: Ring interior radius.
    fe: Ring exterior radius.
    tau: Ring normal opacity.
    theta: Projected tilt. 90 means that the image of the ring
           in perpendicular to orbit in the plane of the sky
    ir: Projected inclindation. 90 for an edge-on ring.

Usage:

    $ python exorings-basic.py [<par>=<value>;<par>=<value>;...]
"""

#############################################################
# INPUT PARAMETERS
#############################################################

#DEFAULT PARAMETER VALUES
default=dict(
    Rs=1.0,#Rsolar
    Ms=1.0,#Msolar
    a=1.0,#AU
    iorb=90.0,#degrees
    fi=1.5,#Rplanet
    fe=2.35,#Rplanet
    Rp=0.08,#Rstar
    tau=1.0,
    theta=30.0,#degrees
    ir=80.0,#degrees
    )

#GET NEW PARAMETER VALUES FROM COMMAND-LINE (IF PROVIDED)
par=getParameters(default)

#############################################################
# PHYSICAL PROPERTIES
#############################################################

#==============================
#ORBITAL PROPERTIES
#==============================
#Stellar density (kg/m^3)
rhotrue=(par.Ms*MSUN)/(4./3*np.pi*(par.Rs*RSUN)**3)

#Period (in seconds)
P=2*np.pi*np.sqrt((par.a*AU)**3/(GCONST*(par.Ms*MSUN)))

#Impact parameter (in stellar radii)
b=(par.a*AU)*np.cos(par.iorb*DEG)/(par.Rs*RSUN)

#Planetary scaled radius
p=par.Rp

#Scaled orbital semimajor axis
asc=(par.a*AU)/(par.Rs*RSUN)

#==============================
#EXTERNAL RING PROPERTIES
#==============================
#Projected semimajor axis
A=par.fe*par.Rp

#Projected semiminor axis
B=A*np.cos(par.ir*DEG)

#==============================
#CHECK TRANSIT CONDITION
#==============================
#Approximate semiheight of the box containing the ringed-planet
hp=max(par.Rp,A*np.sin(par.theta*DEG),B*np.cos(par.theta*DEG))

#Transit condition
if b>par.Rs-hp:
    print "No transit (or grazing transit) occurrs with this orbital inclination, iorb = %.2f deg."%par.iorb
    print "Actual impact parameter: b = %.3f"%b
    print "Maximum impact parameter: bmax ~ %.3f"%(par.Rs-hp)
    exit(1)

#############################################################
# BASIC TRANSIT PROPERTIES
#############################################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#TRANSIT DEPTH
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#==============================
#AUXILIAR VARIABLES
#==============================
cosir=np.cos(par.ir*DEG)
sinir=np.sin(par.ir*DEG)

#==============================
#ABSORTION/SCATTERING FACTOR
#==============================
beta=1-np.exp(-par.tau/cosir)

#==============================
#INTERNAL RING EFFECTIVE RADIUS
#==============================
if par.fi*cosir>1:
    ri2=par.fi**2*cosir-1
else:
    yi=np.sqrt(par.fi**2-1)/(par.fi*sinir)
    ri2=par.fi**2*cosir*2/np.pi*np.arcsin(yi)-\
        2/np.pi*np.arcsin(yi*par.fi*cosir)
ri2=beta*ri2

#==============================
#EXTERNAL RING EFFECTIVE RADIUS
#==============================
if par.fe*cosir>1:
    re2=par.fe**2*cosir-1
else:
    ye=np.sqrt(par.fe**2-1)/(par.fe*sinir)
    re2=par.fe**2*cosir*2/np.pi*np.arcsin(ye)-\
        2/np.pi*np.arcsin(ye*par.fe*cosir)
re2=beta*re2

#==============================
#RINGED-PLANET AREA
#==============================
ARp=np.pi*par.Rp**2+np.pi*(re2-ri2)*par.Rp**2

#==============================
#TRANSIT DEPTH (delta=ARp/As)
#==============================
delta=ARp/(np.pi*par.Rs**2)

#==============================
#OBSERVED RADIUS (p=Rp,obs/Rp)
#==============================
pobs=np.sqrt(delta)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#CONTACT POSITIONS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Contact positions for planetary disk
xp14=np.sqrt((1+p)**2-b**2)
xp1=-xp14;xp4=+xp14
xp23=np.sqrt((1-p)**2-b**2)
xp2=-xp23;xp3=+xp23

#Contact prosition for external ring (approximate solution)
xR13=1-A**2*(np.sin(par.theta*DEG)-b/A)**2*(1-B**2/A)
xR24=1-A**2*(np.sin(par.theta*DEG)+b/A)**2*(1-B**2/A)

xR1=-np.sqrt(xR13)-A*np.cos(par.theta*DEG)
xR2=-np.sqrt(xR24)+A*np.cos(par.theta*DEG)
xR3=+np.sqrt(xR13)-A*np.cos(par.theta*DEG)
xR4=+np.sqrt(xR24)+A*np.cos(par.theta*DEG)

#Determining which contact positions should be used
x1=min(xp1,xR1)
x2=max(xp2,xR2)

x3=min(xp3,xR3)
x4=max(xp4,xR4)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#TRANSIT TIMES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Non-ringed planet
T14p=P*np.arcsin((xp4-xp1)/(asc*np.sin(par.iorb*DEG)))/(2*np.pi)/HOUR
T23p=P*np.arcsin((xp3-xp2)/(asc*np.sin(par.iorb*DEG)))/(2*np.pi)/HOUR

#Ringed planet
T14=P*np.arcsin((x4-x1)/(asc*np.sin(par.iorb*DEG)))/(2*np.pi)/HOUR
T23=P*np.arcsin((x3-x2)/(asc*np.sin(par.iorb*DEG)))/(2*np.pi)/HOUR

#############################################################
# DERIVED TRANSIT PROPERTIES
#############################################################

#Observed scaled semimajor axis
aobs=2*(P/HOUR)/np.pi*delta**0.25/(T14**2-T23**2)**0.5

#Observed impact parameter
bobs=((T14**2*(1-np.sqrt(delta))-T23**2*(1+np.sqrt(delta)))/\
          (T14**2-T23**2))**0.5

#Observed stellar density
rhoobs=(3*np.pi/GCONST)*aobs**3/P**2

#Photo ring effect
PR=rhoobs/rhotrue
logPR=np.log10(PR)

#############################################################
# REPORT
#############################################################
L="*"*60+"\n";R="\n"+"*"*60
print L,"INPUT PARAMETERS",R
print "Stellar properties:"
print "\tRadius: Rs = %.3f Rsun."%par.Rs
print "\tMass: Ms = %.3f Msun."%par.Ms
print "\tDensity: rho_true = %.2f kg/m^3."%rhotrue
print "Planetary properties:"
print "\tRadius: Rp = %.3f Rsun."%par.Rp
print "\tScaled radius: p = %.4f."%p
print "Orbital properties:"
print "\tSemimajor axis: a = %.3f AU, a/Rs = %.2f."%(par.a,asc)
print "\tOrbital inclination: iorb = %.2f deg."%(par.iorb)
print "\tOrbital period: P = %.4f days."%(P/DAY)
print "Ring properties:"
print "\tInternal Ring Radius: Ri = %.3f Rp."%par.fi
print "\tExternal Ring Radius: Re = %.3f Rp."%par.fe
print "\tNormal opacity: tau = %.3f."%par.tau
print "\tProjected inclination (90 deg for edge on): ir = %.2f deg."%par.ir
print "\tProjected tilt: theta = %.2f deg."%par.theta
print "\tBlocking factor: beta = %e."%beta

print L,"DERIVATIVE PROPERTIES",R

print "Ring derivative properties:"
print "\tProjected exrenal ring axes: A = %.4f, B = %.4f Rs"%(A,B)
print "\tEffective ring interior radius: r_i = %e Rs"%np.sqrt(ri2)
print "\tEffective ring exterior radius: r_e = %e Rs"%np.sqrt(re2)
print "\tProjected ring Area: A_Rp = %e"%ARp

print L,"TRANSIT PROPERTIES",R
print "Transit Depth: delta = %.2f ppm"%(delta*1E6)
print "Observed Planetary Radius: pobs = %.6f"%(pobs)
print "Orbital period: P = %.2f days"%(P/DAY)
print "Impact parameter: b = %.4f Rs"%(b)
print "Scaled planetary radius: p = Rp/Rs = %.6f"%(p)
print "Scaled semimajor axis: (a/Rs) = %.2f"%(asc)
print "Planetary radius ratio: pobs/p = %.2f"%(pobs/p)
print "Planet contact positions:"
print "\tx1 = %.3f, x2 = %.3f"%(xp1,xp2)
print "\tx3 = %.3f, x4 = %.3f"%(xp3,xp4)
print "Ring contact positions:"
print "\tx1 = %.3f, x2 = %.3f"%(xR1,xR2)
print "\tx3 = %.3f, x4 = %.3f"%(xR3,xR4)
print "Final Contact positions:"
print "\tx1 = %.3f, x2 = %.3f"%(x1,x2)
print "\tx3 = %.3f, x4 = %.3f"%(x3,x4)
print "Transit duration (non-ringed): T14 = %.4f h, T23 = %.4f"%(T14p,
                                                                 T23p)
print "Transit duration (ringed): T14 = %.4f h, T23 = %.4f"%(T14,T23)
print "Observed scaled semimajor axis: (a/Rs)_obs = %.2f"%(aobs)
print "Observed impact parameter: b_obs = %.4f"%(bobs)
print "Observed stellar density: rho_obs = %.2f"%rhoobs
print "Photo-ring effect: PR = %.4f, log10(PR) = %.2f"%(PR,logPR)
