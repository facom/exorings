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
p=getParameters(default)

#############################################################
# BASIC TRANSIT PROPERTIES
#############################################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#TRANSIT DEPTH
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#==============================
#AUXILIAR VARIABLES
#==============================
cosir=np.cos(p.ir*DEG)
sinir=np.sin(p.ir*DEG)

#==============================
#ABSORTION/SCATTERING FACTOR
#==============================
beta=1-np.exp(-p.tau/cosir)

#==============================
#INTERNAL RING EFFECTIVE RADIUS
#==============================
if p.fi*cosir>1:
    ri2=p.fi**2*cosir-1
else:
    yi=np.sqrt(p.fi**2-1)/(p.fi*sinir)
    ri2=p.fi**2*cosir*2/np.pi*np.arcsin(yi)-\
        2/np.pi*np.arcsin(yi*p.fi*cosir)
ri2=beta*ri2

#==============================
#EXTERNAL RING EFFECTIVE RADIUS
#==============================
if p.fe*cosir>1:
    re2=p.fe**2*cosir-1
else:
    ye=np.sqrt(p.fe**2-1)/(p.fe*sinir)
    re2=p.fe**2*cosir*2/np.pi*np.arcsin(ye)-\
        2/np.pi*np.arcsin(ye*p.fe*cosir)
re2=beta*re2

#==============================
#RINGED-PLANET AREA
#==============================
ARp=np.pi*p.Rp**2+np.pi*(re2-ri2)*p.Rp**2

#==============================
#TRANSIT DEPTH (delta=ARp/As)
#==============================
delta=ARp/(np.pi*p.Rs**2)

#==============================
#OBSERVED RADIUS (p=Rp,obs/Rp)
#==============================
p=np.sqrt(delta*p.Rs**2/p.Rp**2)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#CONTACT TIMES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print "Blocking factor: beta = %e"%beta
print "Effective ring interior radius: r_i = %e"%np.sqrt(ri2)
print "Effective ring exterior radius: r_e = %e"%np.sqrt(re2)
print "Projected ring Area: A_Rp = %e"%ARp

print "*"*60,"\nTRANSIT PARAMETERS\n","*"*60
print "Transit Depth: delta = %.2f ppm"%(delta*1E6)
print "Observed Planetary Radius: p = %.3f"%(p)

print cosir
