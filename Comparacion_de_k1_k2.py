# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 07:39:19 2024

@author: Tomas Fuentes R
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 14:19:37 2024

@author: Tomas Fuentes R
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

"""La monstruosa calculadora de ecuaciones diferenciales para multiples reacciones"""



# Definición de la función para resolver la ecuación diferencial
def ecuacion_diferencial(dep, t, k1, k2):
    
    # Variables dependientes (Concentración de cada especie)
    CA, CB, CC, = dep
   
    
    # Ley de Velocidad usando relativas
    
    r1 = k1 * CA
    r2 = k2 * CB
    
        
    # Sistema de Ecuaciones
    dCAdt = -r1
    dCBdt = r1 - r2
    dCCdt = r2 
    
    return [dCAdt, dCBdt, dCCdt] #Devolvemos la lista de valores de concentraciones dentro de una lista



#---------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------

    
    
def ejecutar_aumenta_k1_k2():
    """Función que muestra la gráfica de K1 y K2 aumentados en un valor específico"""
    # Intervalo de integración
    tspam = np.linspace(0, 10, 100) # De 1 a 10 con 100 puntos

    # Condiciones iniciales
    Ci = [1, 0, 0] #unicamente se alimenta A y B
    solucion = odeint(ecuacion_diferencial, Ci, tspam, args=(2, 1.5)) #Invocamos a la función que resuleve el sistema ODE
    
    #Calculamos las concentraciones
    concentracion_A = solucion[:, 0] #Obtenemos lista de valores para CA
    concentracion_B = solucion[:, 1] #Obtenemos lista de valores para CB
    concentracion_C = solucion[:, 2] #Obtenemos lista de valores para CC
    
    #Graficamos
    fig, at = plt.subplots(figsize=(10, 6))  # Ajusta el tamaño de la figura
    at.plot(tspam, concentracion_A, label='Concentración de A') #Definimos la primera curva 
    at.plot(tspam, concentracion_B, label='Concentración de B') #Definimos la segunda curva 
    at.plot(tspam, concentracion_C, label='Concentración de C') #Definimos la tercera curva
    at.set_xlabel('Tiempo de reacción (minutos)') #Titulo eje X
    at.set_ylabel('Concentración') #Titulo eje Y
    at.set_title('Perfil de Concentraciones en el tiempo (K aumentado)') #Titulo gráfica
    at.legend() #Leyenda
    at.grid(True) #Cuadricula
    plt.show() #Mostramos gráfica

def ejecutar_disminuir_k1_k2():
    """Función que muestra la gráfica de K1 y K2 disminuidos en un valor específico"""
    # Intervalo de integración
    tspam = np.linspace(0, 10, 100) # De 1 a 10 con 100 puntos

    # Condiciones iniciales
    Ci = [1, 0, 0] #unicamente se alimenta A y B
    solucion = odeint(ecuacion_diferencial, Ci, tspam, args=(0.5, 0.1)) #Invocamos la función que resuleve el sistema ODE
    
    #Calculamos las concentraciones
    concentracion_A = solucion[:, 0] #Obtenemos lista de valores para CA
    concentracion_B = solucion[:, 1] #Obtenemos lista de valores para CB
    concentracion_C = solucion[:, 2] #Obtenemos lista de valores para CC
    #Graficamos
    fig, at = plt.subplots(figsize=(10, 6))  # Ajusta el tamaño de la figura
    at.plot(tspam, concentracion_A, label='Concentración de A') #Definimos la primera curva 
    at.plot(tspam, concentracion_B, label='Concentración de B') #Definimos la segunda curva
    at.plot(tspam, concentracion_C, label='Concentración de C') #Definimos la tercera curva
    at.set_xlabel('Tiempo de reacción (minutos)') #Titulo eje X
    at.set_ylabel('Concentración') #Titulo eje Y
    at.set_title('Perfil de Concentraciones en el tiempo (K disminuido)') #Titulo gráfica
    at.legend() #Leyenda
    at.grid(True) #Cuadricula
    plt.show() #MOstramos gráfica

def ejecutar_normales_k1_k2():
    """Función que muestra la gráfica de K1 y K2 originales"""
    # Intervalo de integración
    tspam = np.linspace(0, 10, 100) # De 1 a 10 con 100 puntos

    # Condiciones iniciales
    Ci = [1, 0, 0] #unicamente se alimenta A y B
    solucion = odeint(ecuacion_diferencial, Ci, tspam, args=(1, 0.5)) #Invocamos la función que resuelve el sistema ODE
    
    #Calculamos las concentraciones
    concentracion_A = solucion[:, 0] #Obtenemos lista de valores para CA
    concentracion_B = solucion[:, 1] #Obtenemos lista de valores para CB
    concentracion_C = solucion[:, 2] #Obtenemos lista de valores para CC
    #Graficamos
    fig, at = plt.subplots(figsize=(10, 6))  # Ajusta el tamaño de la figura
    at.plot(tspam, concentracion_A, label='Concentración de A') #Definimos la primera curva
    at.plot(tspam, concentracion_B, label='Concentración de B') #Definimos la segunda curva 
    at.plot(tspam, concentracion_C, label='Concentración de C') #Definimos la tercera curva
    at.set_xlabel('Tiempo de reacción (minutos)') #Titulo eje X
    at.set_ylabel('Concentración') #Titulo eje Y
    at.set_title('Perfil de Concentraciones en el tiempo (K original)') #Titulo gráfica
    at.legend() #Leyenda
    at.grid(True) #Cuadricula
    plt.show() #Mostramos gráfica

def ejecutar_crear_tabla():
    """Función que genera 3 tablas de las concentraciones en el tiempo para CA, CB, CC.
    (Las graficas corresponden a K1 y K2 aumentados disminuidos y originales respectivamente)"""
    
    tspam = np.linspace(0, 10, 100)  # De 1 a 50 con 100 puntos
    Ci = [1, 0, 0]  # Únicamente se alimenta A y B
    solucion_aumentada = odeint(ecuacion_diferencial, Ci, tspam, args=(2, 1.5))  # Solución con K1 y K2 aumentados
    solucion_disminuida = odeint(ecuacion_diferencial, Ci, tspam, args=(0.5, 0.1))  # Solución con K1 y K2 disminuidos
    solucion_normal = odeint(ecuacion_diferencial, Ci, tspam, args=(1, 0.5))  # Solución con K1 y K2 originales
        
    
    print()
    print("Tabla de valores de K aumentados para CA, CB y CC")
    print("Tiempo (min)   |  Concentración de A    |  Concentración de B    |  Concentración de C")
    for i in range(0, len(tspam) - 1, 10):
        print(f"{tspam[i]:<14.0f} | {solucion_aumentada[i][0]:>20.4f} | {solucion_aumentada[i][1]:>20.4f} | {solucion_aumentada[i][2]:>20.4f}")

    print()
    print("Tabla de valores de K disminuidos para CA, CB y CC")
    print("Tiempo (min)   |  Concentración de A    |  Concentración de B    |  Concentración de C")
    for i in range(0, len(tspam) - 1, 10):
        print(f"{tspam[i]:<14.0f} | {solucion_disminuida[i][0]:>20.4f} | {solucion_disminuida[i][1]:>20.4f} | {solucion_disminuida[i][2]:>20.4f}")

    print()
    print("Tabla de valores de K originales para CA, CB y CC")
    print("Tiempo (min)   |  Concentración de A    |  Concentración de B    |  Concentración de C")
    for i in range(0, len(tspam) - 1, 10):
        print(f"{tspam[i]:<14.0f} | {solucion_normal[i][0]:>20.4f} | {solucion_normal[i][1]:>20.4f} | {solucion_normal[i][2]:>20.4f}")

    """Las anteriores lineas de código únicamente buscan generar un formato para los valores de las tres tablas usando metodos de strings"""

#---------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------

def mostrar_menu():
    """Funcióon que muestra el menú de opciones del programa"""
    print("Ingrese una opcion:")
    print("1. Ver gráfica de K1 y K2 aumentados (Aumento de 1 unidad)")
    print("2. Ver gráfica de K1 y K2 disminuidos (Disminución de 0.5 unidad para K1 y 0.4 para K2 )")
    print("3. ver tabla de valores para K")
    print("4. Ver gráfica con valores originales")
    print("5. salir del programa")


def iniciar_aplicacion():
    """ Funcion que permite elegir una de las opcines en bucle para comodida del usuario"""
    mostrar_menu()
    opcion = input("Ingrese una opcion: ")
    while opcion != "5":
         if opcion == "1":
             ejecutar_aumenta_k1_k2()
         elif opcion == "2":
             ejecutar_disminuir_k1_k2()
         elif opcion == "3":
             ejecutar_crear_tabla()
         elif opcion == "4":
             ejecutar_normales_k1_k2()

         mostrar_menu()     
         opcion = input("Ingrese una opcion nuevamente: ")
    
    
#Programa principal
iniciar_aplicacion() 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    