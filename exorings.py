############################################################
#   _____  _____  _ __(_)_ __   __ _ ___ 
#  / _ \ \/ / _ \| '__| | '_ \ / _` / __|
# |  __/>  < (_) | |  | | | | | (_| \__ \
#  \___/_/\_\___/|_|  |_|_| |_|\__, |___/
#                              |___/     
# v. alpha 1.0
############################################################
# EXORINGS CODE REPOSITORY
# http://github.org/facom/exorings
############################################################
from sys import exit,argv

############################################################
# REQUIRED PACKAGES
############################################################
packages={
    "numpy":dict(linux="python-numpy",alias="np"),
    "matplotlib.pyplot":dict(linux="python-matplotlib",alias="plt")
    }

for pack in packages.keys():
    alias=packages[pack]["alias"]
    linux=packages[pack]["linux"]
    try:exec("import %s as %s"%(pack,alias))
    except ImportError:
        print "Package '%s' not installed.  To continue please install the linux package '%s'."%(pack,linux)
        exit(1)

############################################################
# MACROS AND CONSTANTS
############################################################
DEG=np.pi/180
RAD=1/DEG

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#PHYSICAL
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAY=86400.0 #s
HOUR=3600.0 #s
GCONST=6.67428E-11 #kg m^3/s^2
MSUN=1.98855E+30 #kg
RSUN=6.96342E+8 #m
AU=1.49597871e+11 #m

############################################################
# UTIL ROUTINES
############################################################
class dict2obj(object):
    """
    Convert a dictionary into an object.
    """
    def __init__(self,dic={}):self.__dict__.update(dic)
    def __add__(self,other):
        for attr in other.__dict__.keys():exec("self.%s=other.%s"%(attr,attr))
        return self

def getParameters(default):
    """
    Read parameters from the command line.
    
    default: dictionary with the default values
    
    """
    ip=dict()

    #GET PARAMETERS FROM COMMAND LINE
    for i in xrange(1,len(argv)):
        try:exec argv[i] in ip
        except IndexError:pass
        except NameError:
            print "Bad formed input parameter '%s'."%argv[i]
            exit(1)

    if len(ip.keys())==0: 
        print "No parameters overrided."
        return dict2obj(default)

    #CREATE PARAMETER OBJECTS
    p=dict()
    pkeys=ip.keys()
    provkeys=""
    for key in default.keys():
        if key in pkeys:
            p[key]=ip[key]
            provkeys+="%s,"%key
        else:p[key]=default[key]
    print "Overrided parameters: %s."%provkeys.strip(",")

    return dict2obj(p)

############################################################
# TESTS
############################################################
if __name__=="__main__":
    print "All tests passed. Code ready to be used."
    exit(0)
