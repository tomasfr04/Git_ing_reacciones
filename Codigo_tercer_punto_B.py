# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 09:40:28 2024

@author: Tomas Fuentes R
"""

""" La super calculadora de concentraciones finales de un sistema de multiples reacciones """


from sympy import *  
from scipy.optimize import fsolve
import numpy as np

from scipy.optimize import fsolve
import numpy as np

# Definir las ecuaciones como funciones de Python
def Sistema_de_ecuaciones(incognitas):
    
    CA, CB, CC, CD, CE, CF = incognitas # Definimos las incognitas del sistema de ecuaciones
    
    Eq1 = 15 - (10*CA) - (12.5*CA*(CB**2)) - (15*CA*CD) # Ecuacion 1
    Eq2 = 20 - (10*CB) - (25*CA*CB**2) - (250*CB*CC**2) #Ecuacion 2
    Eq3 = (-10*CC) + (12.5*CA*CB**2) + (5*CA*CD) - (500*CB*CC**2) #Ecuacion 3
    Eq4 = (-10*CD) + (12.5*CA*CB**2) - (10*CA*CD) + (250*CB*CC**2) #Ecuacion 4
    Eq5 = -10*CE + 5*CA*CD # Ecuacion 5
    Eq6 = -10*CF + 250*CB*CC**2 # Ecuacion 6
    
    return [Eq1, Eq2, Eq3, Eq4, Eq5, Eq6] #Retornamos una lista con las concentraciones finales

# Estimación inicial
valores_iniciales = [1, 1, 1, 1, 1, 1] #Valores iniciales para cada especie del sistema

# Resolver el sistema de ecuaciones
solucion = fsolve(Sistema_de_ecuaciones, valores_iniciales) #LLamamos a la funcion

"""Formato de visualizacion bonito """

print("La concentración final para cada especie es:")
print()
print("Para A: {} mol/dm3".format(round(solucion[0],3)))
print("Para B: {} mol/dm3".format(round(solucion[1],3)))
print("Para C: {} mol/dm3".format(round(solucion[2],3)))
print("Para D: {} mol/dm3".format(round(solucion[3],3)))
print("Para E: {} mol/dm3".format(round(solucion[4],3)))
print("Para F: {} mol/dm3".format(round(solucion[5],3)))

 #Conversion   X = (Cantidad inicial de A - Cantidad que reaccino de A) / Cantidad inicial de A
FA0 = 15 #mol/min
vo = 10 #dm3/min
 
x = (FA0 - (solucion[0]*vo))/FA0
print()
print("La conversion final es de {}".format(round(x,3), )) 

 
 
 
 
 
 