#This code has for objective to calculate the steady state temperature of a wire of a certain diameter carrying a certain current.
#Author: Enrique Morell
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

#constant
#Stephan Boltzmann Constant
sigma = 5.67e-8


class Power:
    def __init__(self, name, symbol, equation, *arg, **kwargs):
        self.name = name
        self.symbol = symbol
        self.equation = equation

class Wire:
    def __init__(self, name, R, OD, emissivity, lunit = "imp"):
        #R is resistance PER UNIT LENGTH in (Ohm/km)
        #lunit = "imp" for imperial units (applies to OD) assumes OD is in inches, "m" for metric (it will be assumed that OD will be then given in centimeters)
        self.name = name
        self.resistance = R/1000 #transform to SI
        self.outdiameter = OD*2.54*1e-2 if lunit == "imp" else OD*1e-2
        print(f"Wire OD = {self.outdiameter}")
        self.emissivity = emissivity

class Ambient:
    def __init__(self, name, ambient_temp, tcond, windspeed, viscosity, wind_direction):
        self.name = name
        self.Ta = ambient_temp
        self.conductivity = tcond
        self.windspeed = windspeed
        self.viscosity = viscosity
        self.wind_direction = wind_direction*np.pi/180 #in rad

def _pr(epsilon, D, Tc, Ta = 25):
    return sigma*epsilon*np.pi*D*((Tc+273.15)**4-(Ta+273.15)**4)

def Tc(Pj, D, epsilon = 0.86, Ta =25):
    sqr = Pj/(sigma*epsilon*np.pi*D)+(Ta+273.15)**4
    return np.float_power(sqr, 0.25)

def Pj(wire : Wire, current):
    #returns power per unit length when this current goes through this Wire.
    #result is in Watts per meter
    return wire.resistance*current**2

def Pr(Tc, wire : Wire, ambient: Ambient):
    #TK = Tc+273.15#
    return sigma*wire.emissivity*np.pi*wire.outdiameter*((Tc+273.15)**4 - (ambient.Ta+273.15)**4)

def Pc(Tc, wire : Wire, ambient: Ambient):
    Re = wire.outdiameter*ambient.windspeed/ambient.viscosity
    Nu = 0.64*np.power(Re, 0.2)+0.2*np.power(Re, 0.61)
    Kwd = 1.194-np.sin(ambient.wind_direction)-0.194*np.cos(2*ambient.wind_direction)+0.364*np.sin(2*ambient.wind_direction)
    h = ambient.conductivity*Nu*Kwd
    return h*np.pi*(Tc-ambient.Ta)

wire = Wire("3/0", 0.202704, 0.675, 0.86)

ambient = Ambient("Lab air", 21.66, 0.0262, 0.1, 1.81e-5, 90)

current = 288/1 # current in amps per wire (Sr uses 288 A total)

def SteadyS_Equation(Tc, _wire = wire, _ambient = ambient, _current = current):
    value = Pr(Tc, _wire, _ambient) + Pc(Tc, _wire, _ambient) - Pj(_wire, _current)
    # print("Joule = {} W/m".format(Pj(_wire, _current)))
    # print("Radiation = {} W/m".format(Pr(Tc, _wire, _ambient)))
    # print("Convection = {} W/m".format(Pc(Tc, _wire, _ambient)))
    # print("Temperature = {} C".format(Tc))
    # print(value)
    return value

#find the zeros in function of Tc using numpy:
conv_temp_rootout = root_scalar(SteadyS_Equation, x0=20, x1=120, rtol=0.0001)
convection_temp = conv_temp_rootout.root
print(f"Temp of cable using convection (deg C)= {convection_temp}")

#Joule heating:
no_convection_temp = Tc(Pj(wire, current), 1.7145e-2)-273.15
print(f"Temp of cable, no convection (deg C) = {no_convection_temp}")

#epdm (insulation) has epsilon = 0.86

#print("at Ta = {}".format(SteadyS_Equation(ambient.Ta+6e-1)))

X = np.linspace(ambient.Ta, ambient.Ta + 80, 500)
Y = Pr(X, wire, ambient)

plt.plot(X, Y)
plt.show()