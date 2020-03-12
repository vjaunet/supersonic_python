#
#     isentropic relations for python
#
#=============================================

from numpy import sqrt

def Mj_NPR(NPR):
    """ Mj_NPR(NPR)
    provides  fully expanded Mach Number from nozzle NPR
    """
    g=1.4
    return (2.0/(g-1)*(NPR**((g-1)/g)-1))**0.5

def Uj_Mj(Mj,T0=293):
    """ Uj_Mj(Mj,T0=293)
    provides fully expanded jet velocity
    from Mj and Total temperature
    """
    from numpy import sqrt
    g=1.4
    R=287.0
    return sqrt(g*R*T0*(1. + (g-1)/2.*Mj**2)**(-1))*Mj


# Get Dj from D,Mj,Mdesign
def Dj_Mj(D,Mj,Md):
    """ Dj_Mj(D,Mj,Md)
    provides fully expanded jet Diameter
    from D, Mj and design Mach number
    """
    g=1.4;
    return D*((1+0.2*Mj**2)/((1+0.2*Md**2)))**((g+1)/(4*(g-1)))*(Md/Mj)**(1/2);

# isentropic pressure ratio from Mach number
def tti_M(M):
    """ isentropic temperature ratio """
    g=1.4;
    return (1+(g-1)/2*M**2)**(-1)

def ppi_M(M):
    """
    Isentropic pressure ratio from Mach
            ppi_M(M)
    """
    g=1.4;
    return (1+(g-1)/2*M**2)**(-g/(g-1))

def M_ppi(ppi):
    """ Inverted isentropic pressure ratio """
    g=1.4;
    return sqrt(2/(g-1)*ppi**(g/(g-1) - 1 ))


def M_ppi(ppi):
    """ Inverted isentropic pressure ratio from M """
    g=1.4;
    return sqrt(2/(g-1)*ppi**(g/(g-1) - 1 ))

def aac_M(M):
    """ A/A* as a function of Mach
        returns: aac(M)
    """
    g=1.4;
    return 1.0/M*((1.0 + (g-1.0)/2.0*M**2)*2/(g+1))**((g+1)/2/(g-1))

def M_aac(Ac,A):
    """ M as function of A and Ac
        returns: M_aac(Ac,A)
    """
    from scipy.optimize import fsolve
    g=1.4;

    def func(M):
         return aac_M(M) - A/Ac

    M0_sub = fsolve(func,0.01,xtol=1e-10,maxfev=500)
    M0_sup = fsolve(func,1.02,xtol=1e-10,maxfev=500)

    return {'subsonic':M0_sub[0],'supersonic':M0_sup[0]}

def nu_M(M):
    """ Prandtl-Meyer relation """

    return 0
