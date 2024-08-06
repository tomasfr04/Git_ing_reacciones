# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 15:46:05 2024

@author: Tomas Fuentes R
"""
#Importaciones
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

"""La monstruosa calculadora de ecuaciones diferenciales para multiples reacciones"""

"Simulación de un reactor PBR isotérmico y no isobárico con multiples reacciones"

"Reacciones"
# 2A + 3B -> C
# 3C + A -> 2D
# D + 2C -> 2E

# Definición de la función para resolver la ecuación diferencial
def ecuacion_diferencial(dep, w):
    
    # Variables dependientes (Concentración de cada especie)
    FA, FB, FC, FD, FE, y = dep
    FTP = FA+FB+FC+FD+FE
    
    # Parámetros - constantes - etc
    k1B = 6
    k2D = 8
    k3E = 7
    alpha = 0.0006
    CA0=0.5
    vo = 165
    Po = CA0 * 0.082057 *(50+273)
    FA0 = vo * CA0
    P = (FA/vo) * 0.082057 *(50+273)
    # Concentraciones
    CA = 0.5 * (FA/82.5)*(P/Po)
    CC = 0.5 * (FC/82.5)*(P/Po)
    CD = 0.5 * (FC/82.5)*(P/Po)
    
    
    # Ley de Velocidad usando relativas
    r1B = k1B * CA **2
    r2D = k2D * CC * CA
    r3E = k3E * CD * CC
    
        
    # Sistema de Ecuaciones
    dFAdw = (-(2/3)*r1B)-((1/2)*r2D) 
    dFBdw = (r1B) 
    dFCdw = ((1/3)*r1B)-(-r3E)-((3/2)*r2D)  
    dFDdw = (r2D)-((1/2)*r3E) 
    dFEdw = r3E
    FT = dFAdw + dFBdw + dFCdw + dFDdw + dFEdw   
    if y > 0:
        dydw = (-alpha) / (2 * y) 
    else: 
        dydw = 0
    return [dFAdw, dFBdw, dFCdw,dFDdw,dFEdw,dydw]

# Intervalo de integración
Wspam = np.linspace(0, 2500, 2500) # De 1 a 50 con 100 puntos

# Condiciones iniciales
CA0=0.5 #Mol/L
vo = 165 #L/min
FA0 = vo * CA0 #Mol/min   

Ci = [FA0, 0, 0, 0, 0, 1] #unicamente se alimenta A y B

# Resolución de la ecuación diferencial
solucion = odeint(ecuacion_diferencial, Ci, Wspam)

# Obtención de las soluciones
Flujo_A = solucion[:, 0]
Flujo_B = solucion[:, 1]
Flujo_C = solucion[:, 2]
Flujo_D = solucion[:, 3]
Flujo_E = solucion[:, 4]
caida_de_presion = solucion[:,5]
# Graficar
fig, (at,ax) =plt.subplots(2, 1, figsize=(8, 10))  # Ajusta el tamaño de la figura

""""----------------------------------------------------------------------------"""

ax.plot(Wspam,caida_de_presion)
ax.set_xlabel('Peso del catalizador') # Título eje X
ax.set_ylabel('Caida de presión (y)') # Título eje Y
ax.set_title('Perfil de caida de presion') # Título gráfica
#ax.set_xlim(0, 50)
#ax.set_ylim(-1, 1)
ax.grid(True) # Cuadrícula
#ax.set_xticks(np.arange(0, 50, 5))  # De 0 a 2500, en intervalos de 250
#ax.set_yticks(np.arange(-1, 1, 0.5)) 

"""-----------------------------------------------------------------------------"""
at.plot(Wspam, Flujo_A, label='Flujo de A') 
at.plot(Wspam, Flujo_B, label='Flujo de B')
at.plot(Wspam, Flujo_C, label='Flujo de C')
at.plot(Wspam, Flujo_D, label='Flujo de D')
at.plot(Wspam, Flujo_E, label='Flujo de E')
at.set_xlabel('Peso del catalizador') # Título eje X
at.set_ylabel('Flujo por especie') # Título eje Y
at.set_title('Perfil de flujos') # Título gráfica
at.set_xlim(0, 2500)
at.set_ylim(0, 110)
at.legend() # Leyenda
at.grid(True) # Cuadrícula

at.set_xticks(np.arange(0, 2501, 250))  # De 0 a 2500, en intervalos de 250
at.set_yticks(np.arange(0, 111, 10)) # DE a 0 100 en intervalos de 10
plt.show()
"""-----------------------------------------------------------------------------------------------------------------------------"""

print("Podemos ver en la primera gráfica que la reacción se estabiliza mas o menos en 750 gramos de catalizador, por lo tanto es mas \
      conveniente para la reacción trabar con es cantidad. Aproximadamente en 1300 gramos de catalizador la caida de presión se vuelve 0 imposibilitando \
          la reaccion. Para concluir es nas efectivo traabajar el reactor entre 700 y 1500 gramos de catalizador")