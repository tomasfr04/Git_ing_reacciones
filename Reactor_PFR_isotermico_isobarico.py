# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 14:19:37 2024

@author: Tomas Fuentes R
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

"""La monstruosa calculadora de ecuaciones diferenciales para multiples reacciones"""

#Simulación de un reactor PFR isotérmico y isobárico con multiples reacciones

"Reaccoines"
# A + 2B -> C + D
# 2D + 3A -> C + E
# B + 2C -> D + F

# Definición de la función para resolver la ecuación diferencial
def ecuacion_diferencial(dep, v):
    
    # Variables dependientes (Concentración de cada especie)
    CA, CB, CC, CD, CE, CF = dep
   
    # Parámetros - constantes
    k1D = 0.25 # dm6/mol2*min #CAMBIAR 
    k2E = 0.1 # dm3/mol*min #CAMBIAR
    k3F = 5 # dm9/mol2*min #CAMBIAR
    vo = 10 #dm3/min #CAMBIAR
    
    # Ley de Velocidad usando relativas
    r1D = k1D * CA * CB**2
    r2E = k2E * CA * CD
    r3F = k3F * CB * CC**2
    
        
    # Sistema de Ecuaciones
    dAdV = (-r1D-(3*r2E)) / vo
    dBdV = (-2*r1D-r3F) / vo
    dCdV = (r1D+r2E-2*r3F) / vo 
    dDdV = (r1D-2*r2E+r3F) / vo
    dEdV = r2E / vo
    dFdV = r3F / vo
    
    return [dAdV, dBdV, dCdV,dDdV,dEdV,dFdV]

# Intervalo de integración
Vspam = np.linspace(0, 50, 100) # De 1 a 50 con 100 puntos

# Condiciones iniciales
Ci = [1.5, 2, 0, 0, 0, 0] #unicamente se alimenta A y B mol/dm3

# Resolución de la ecuación diferencial
solucion = odeint(ecuacion_diferencial, Ci, Vspam)

# Obtención de las soluciones
concentracion_A = solucion[:, 0]
concentracion_B = solucion[:, 1]
concentracion_C = solucion[:, 2]
concentracion_D = solucion[:, 3]
concentracion_E = solucion[:, 4]
concentracion_F = solucion[:, 5]

#print((concentracion_B))
#print((Vspam.tolist()))
# Graficar
fig, at = plt.subplots(figsize=(8, 6))  # Ajusta el tamaño de la figura
at.plot(Vspam, concentracion_A, label='Concentración de A') 
at.plot(Vspam, concentracion_B, label='Concentración de B')
at.plot(Vspam, concentracion_C, label='Concentración de C')
plt.plot(Vspam, concentracion_D, label='Concentración de D')
plt.plot(Vspam, concentracion_E, label='Concentración de E')
plt.plot(Vspam, concentracion_F, label='Concentración de F')
at.set_xlabel('Volumen del Reactor') #Titulo eje X
at.set_ylabel('Concentración') #Titulo eje Y
at.set_title('Perfil de Concentraciones') #Titulo gráfica
at.legend() #Leyenda
at.grid(True) #Cuadricula
plt.show()