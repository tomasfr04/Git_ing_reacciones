# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 15:08:23 2024

@author: Tomas Fuentes R
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Constantes - variables

k1Ao = 1 #dm^6/mol^2*s
k2Ao = 0.5 #s^-1    
k3Bo = 2 #dm^12/mol^4*s
delta_Hrx1o = -20000 #J/mol de A
delta_Hrx2o = 11500 #J/mol de B
delta_Hrx3o = -35000 #J/mol de C
E1 = 6000 #J/mol
E2 = 5000 #J/mol
E3 = 7000 #J/mol
CpA = 10 #J/mol*K
CpB = 20 #J/mol*K
CpC = 40 #J/mol*K
CpD = 20 #J/mol*K
CpE = 100 #J/mol*K
CpI = 30 #J/mol*K
Ua = 450 #J/dm^3*s*K
Vo = 10 #dm^3/s
To = 400 #K
Ta = 360 #K
R = 8.314 #J/mol*K

CA0 = 3 #mol/L
CB0 = 1.1 #mol/L
CC0 = 0 #mol/L
CD0 = 0 #mol/L
CE0 = 0 #mol/L
CI0 = 0.05 #mol/L

CT0 = CA0+CB0+CC0+CE0+CD0 #mol/L 
FT0 = CT0 * Vo #mol/s
FA0 = CA0 * Vo #mol/s
FB0 = CB0 * Vo #mol/s
FC0 = CC0 * Vo #mol/s
FD0 = CD0 * Vo #mol/s
FE0 = CE0 * Vo #mol/s
FI0 = CI0 * Vo #mol/s

delta_Cp1 = 0.5*CpC - 0.5*CpB - CpA
delta_Cp2 = CpD - CpA
delta_Cp3 = 0.5*CpE - CpC - 0.5*CpB- CpD


# Definición de la función para resolver la ecuación diferencial
def ecuacion_diferencial(dep, v):
    
    # Variables dependientes (Temperatura y flujo de cada especie)
    T, FA, FB, FC, FD, FE = dep
   
    #Concentraciones - estequiometría
    FT = FA+FB+FC+FD+FE+FI0
    CA = CT0*(FA/FT)*(To/T)
    CB = CT0*(FB/FT)*(To/T)
    CC = CT0*(FC/FT)*(To/T)
    CD = CT0*(FD/FT)*(To/T)
    
    #Corregimos los valores que dependen de la temperatura
    k1A = k1Ao * np.exp((E1/R)*((1/273)-(1/T))) #Siempre corregimos con respecto a la temperatura del ejercicio
    k2A = k2Ao * np.exp((E2/R)*((1/293)-(1/T))) #Siempre corregimos con respecto a la temperatura del ejercicio
    k3B = k3Bo * np.exp((E3/R)*((1/313)-(1/T))) #Siempre corregimos con respecto a la temperatura del ejercicio
    deltaHrx1 = delta_Hrx1o + (delta_Cp1*(T-400)) #Siempre corregimos con respecto a la temperatura del ejercicio
    deltaHrx2 = delta_Hrx2o + (delta_Cp2*(T-400)) #Siempre corregimos con respecto a la temperatura del ejercicio
    deltaHrx3 = delta_Hrx3o + (delta_Cp3*(T-400)) #Siempre corregimos con respecto a la temperatura del ejercicio

    
    deltaHrx = deltaHrx1 + deltaHrx2 + deltaHrx3 #J/mol
    
    #ley de velocidad 
    r1A = -k1A * ((CA**2)*CB) 
    r2A = -k2A * (CA) 
    r3B = -k3B * ((CC**2)*CB*(CD**2)) 
    r3C = -k3B * ((CC**2)*CB*(CD**2)) 

    #ley de velocidades netas
    rA = r1A + r2A
    rB = 0.5*r1A + r3B
    rC = -0.5*r1A + 2*r3B
    rD = -r2A + 2*r3B
    rE = -r3B

    # Sistema de ecuaciones diferenciales
    dTdV = (Ua*(Ta-T)+((r1A*deltaHrx1)+(r2A*deltaHrx2)+(2*r3B*deltaHrx3)))/((FA*CpA)+(FB*CpB)+(FC*CpC)+(FD*CpD)+(FE*CpE))
    dFAdV = rA #Balance de masa para A
    dFBdV = rB #Balance de masa para B
    dFCdV = rC #Balance de masa para C
    dFDdV = rD #Balance de masa para D
    dFEdV = rE #Balance de masa para E
    
    return [dTdV,dFAdV,dFBdV,dFCdV,dFDdV,dFEdV] #Devolvemos en forma de lista

# Intervalo de integración
Vspam = np.linspace(0, 100, 100) # De 0 a 100 con 100 puntos

# Condiciones iniciales
Ci = [400,FA0,FB0,FC0,CD0,FE0] #Temperatura y flujos iniciales

# Resolución de la ecuación diferencial
solucion = odeint(ecuacion_diferencial, Ci, Vspam) 

# Obtención de las soluciones
temperatura = solucion[:,0] #Obtenemos la lista de temperaturas
flujo_de_A = solucion[:,1] #Obtenemos la lista del flujo de A
flujo_de_B = solucion[:,2] #Obtenemos la lista del flujo de B
flujo_de_C = solucion[:,3] #Obtenemos la lista del flujo de C
flujo_de_D = solucion[:,4] #Obtenemos la lista del flujo de D
flujo_de_E = solucion[:,5] #Obtenemos la lista del flujo de E


fig, at = plt.subplots(figsize=(8, 6))  # Ajusta el tamaño de la figura
at.plot(Vspam,flujo_de_A, label='Flujo de A') #Definimos los valores a gráficar 
at.plot(Vspam,flujo_de_B, label='Flujo de B') #Definimos los valores a gráficar 
at.plot(Vspam,flujo_de_C, label='Flujo de C') #Definimos los valores a gráficar 
at.plot(Vspam,flujo_de_D, label='Flujo de D') #Definimos los valores a gráficar 
at.plot(Vspam,flujo_de_E, label='Flujo de E') #Definimos los valores a gráficar 
at.set_xlabel('Volumen del Reactor') #Titulo eje X
at.set_ylabel('Flujos (mol/s)') #Titulo eje Y
at.set_title('Perfil de flujos con Ta constante') #Titulo gráfica
at.legend() #Leyenda
at.grid(True) #Cuadricula


fig, ax = plt.subplots(figsize=(8, 6))  # Ajusta el tamaño de la figura
ax.plot(Vspam,temperatura, label='Temperatura a lo largo del reactor') #Definimos los valores a gráficar
ax.set_xlabel('Volumen del Reactor') #Titulo eje X
ax.set_ylabel('Temperatura (K)') #Titulo eje Y
ax.set_title('Perfil de temperatura con Ta constante') #Titulo gráfica
ax.legend() #Leyenda
ax.grid(True) #Cuadricula


plt.show()











