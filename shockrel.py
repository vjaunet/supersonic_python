#
#     shock relations for python
#
#  author  : vincent jaunet,
#            vincent [dot] jaunet [at] ensma [dot] fr
#  date    : 03/2020
#  License : MIT
#
#=============================================

import numpy as np
import isentropic
from scipy.optimize import fsolve

def beta_theta_M(theta,M):
    """ beta_theta_M(theta,M)
    shock angle from deviation and Mach number
    """
    g=1.4

    theta = theta/180*np.pi

    def rel(beta):
        relation = np.tan(theta) - 2*np.cos(beta)/np.sin(beta)*(
            (M**2*np.sin(beta)**2-1)/(M**2*(g + np.cos(2*beta))+2))
        return relation

    beta_weak = fsolve(rel,0.01,xtol=1e-10,maxfev=500)
    beta_strong = fsolve(rel,np.pi/2+0.01,xtol=1e-10,maxfev=500)

    return {'weak':beta_weak[0]/np.pi*180,'strong':beta_strong[0]/np.pi*180}

def theta_beta_M(beta,M):
    """ beta_theta_M(theta,M)
    shock angle from deviation and Mach number
    """
    g=1.4

    beta = beta/180*np.pi

    def rel(theta):
        relation = np.tan(theta) - 2*np.cos(beta)/np.sin(beta)*(
            (M**2*np.sin(beta)**2-1)/(M**2*(g + np.cos(2*beta))+2))
        return relation

    theta = fsolve(rel,1,xtol=1e-10,maxfev=500)

    return abs(theta[0]/np.pi*180.0)


def M2_nshock(M1):
    """ Mach number after a normal shock"""
    g=1.4
    return np.sqrt((1.+(g-1.0)/2.0*M1**2)/(g*M1**2-(g-1.0)/2.0))


def t2t1_nshock(M1n):
    """ static pressure ratio across the shock """
    g=1.4
    return ((2*g*M1n**2-(g-1))*((g-1)*M1n**2+2)/((g+1)**2*M1n**2))

def p2p1_nshock(M1n):
    """ static pressure ratio across the shock """
    g=1.4
    return (1 + 2*g/(g+1.0)*(M1n**2 - 1.0))

def pt2pt1_nshock(M1n):
    """ static pressure ratio across the shock """
    g=1.4
    p2p1 = p2p1_nshock(M1n)
    pi_m1 = isentropic.ppi_M(M1n)
    pi_m2 = isentropic.ppi_M(M2_nshock(M1n))
    return p2p1/pi_m2*pi_m1
