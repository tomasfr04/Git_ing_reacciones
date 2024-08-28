# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 23:44:49 2024

@author: Tomas Fuentes R
"""

"Simulacion de un reactor PBR isotermico y no isobarico "


#Importaciones
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#Reaccion
#3A + 2B -> C + 2D

# Definición de la función para resolver la ecuación diferencial
def ecuacion_diferencial(dep, w):
    
    # Variables dependientes
    x , y = dep
    
    # Parámetros
    alpha = 0.0007 #CAMBIAR
    T0 = 150 + 273.15 # K #CAMBIAR
    P0 = 12 # atm #CAMBIAR
    CA0 = P0 / (0.08205 * T0) #mol/L
    FA0 = 6 # gmol/min #CAMBIAR
    k = 4.39 # dm^12/(mol^4 gcat min) #CAMBIAR
    
    # Estequiometría
    CA = (CA0 * (1 - x)/(1-(2/3)*x)) * y
    CB = (CA0 * (1 - x)/(1-(2/3)*x)) * y
    
    # Ley de Velocidad elemental
    rA = -k * ((CA ** 3) * (CB ** 2))
    
    # Sistema de Ecuaciones
    ecuacion = [ -rA / FA0, -alpha*(1-(1/3)*x) / (2 * y) if y > 0 else 0]
    
    return ecuacion

# Intervalo de integración
wspan = np.linspace(0,2000, 100)

# Condiciones iniciales
y0 = [0, 1]

# Resolución de la ecuación diferencial
solucion = odeint(ecuacion_diferencial, y0, wspan)

# Obtención de las soluciones
x = solucion[:, 0]
y = solucion[:, 1]
# Graficar
indice=0
indice_interseccion = np.argwhere(np.diff(np.sign(x - y))).flatten()[0]
indice_presion_cercana_a_0 = np.where(y < 0.01)[0][0]
plt.text(wspan[indice_presion_cercana_a_0], y[indice_presion_cercana_a_0], f'({wspan[indice_presion_cercana_a_0]:.2f}, {y[indice_presion_cercana_a_0]:.2f})')


plt.plot(wspan, y, label='Caida de presión')
plt.plot(wspan, x, label='Conversión')
plt.plot(wspan[indice_interseccion], y[indice_interseccion], 'ro')
plt.text(wspan[indice_interseccion], y[indice_interseccion], f'({y[indice_interseccion]:.2f})')
plt.plot(wspan[indice_presion_cercana_a_0], y[indice_presion_cercana_a_0], 'ro')
plt.text(wspan[indice_presion_cercana_a_0], y[indice_presion_cercana_a_0], f'({wspan[indice_presion_cercana_a_0]:.2f}')


plt.ylabel('Conversión y Caida de Presión')
plt.xlabel('Peso Catalizador (g)')
plt.legend()
plt.title('Punto 2')
plt.grid(True)
plt.minorticks_on()
plt.show()