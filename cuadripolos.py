import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from sympy import lambdify
import pandas as pd
import seaborn as sns

sns.set()


class Cuadripolo(object):
    """Clase para definir cuadripolos como representación de sistemas de silenciadores
    :inputs
    RESPETAR LAS UNIDADES ESTABLECIDAS [m], [s], [kg]
        Existe la posibilidad de no definir ningun parametro y utilizar
        la clase simplemente como contenedor de un cuadripolo. Con el fin de poder multiplicar y
        guardar cuadripolos que reprensenten sistemas complejos
        s : float. seción interior.
        largo : Largo del tubo_recto de la camara, etc... Caso de ser helmholtz el largo es del tubo de entrada
        s1: float. sección de entrada
        s2: float. sección de salida
        c : float, 345 [m/s]
        rho_o : float, 1,225 [kg/m3]. Revisar si esta bien este kilogramo
        tipo : str. El tipo de elemento fundamental del silenciador
                    que va a reprenstar el cuadripolo:
                    camara, tubo_recto, helmholtz, extension_expansion, tubo_cerrado, Z_in.
                    En caso de que el cuadripolo sea producto de la multiplicación de varios cuadripolos el
                    tipo va a ser 'complejo' en relación a que representa una estructura de silenciador compleja y no fundamental
        Z_in (parametro): si en tipo se usa Z_in, Z_in debe ser definido como una función
                          sympy dependiente de la variable "f" (frecuencia). De otra manera no
                          va a ser posible evaluarla.
        vol: int [m3], volumen de resonador de helmholtz.
    """

    def __init__(self, s=None, largo=None, s1=None, s2=None, c=345, rho_o=1.225, tipo='complejo', vol=None, Z_in=None):
        # ninguna propiedad va a ser privada por lo cual no necesitamos uso de getters o setters
        self.c = c
        self.rho_o = rho_o
        self.s = s
        self.s1 = s1
        self.s2 = s2
        self.vol = vol
        self.tipo = tipo
        self.f = symbols('f')
        self.largo = largo
        self.Z_o = self.rho_o * c
        if tipo == 'helmholtz':
            # sumo el largo del tubo? esto es así?
            self.largo += np.square(self.s / np.pi) * 0.82
            #las secciones son todas las del tubo de entrada ?
            self.s1 = self.s
            self.s2 = None
            try:
                self.Z_in = self.rho_o * (((self.largo * 2 * np.pi * self.f) / self.s) - (
                        (self.c ** 2) / (self.vol * 2 * np.pi * self.f)))
            except Exception as e:
                raise print(e,
                            'el volumen y la sección del tubo_recto de ingreso deben estar definidas para el resonador de hemholtz')
        else:
            self.Z_in = Z_in
        if self.tipo == 'tubo_recto':
            self.s2 = self.s1 = self.s
        self.cuadri = [[0, 0], [0, 0]]
        self.tl = None
        self.tl_ = None

    def __getitem__(self, tup):
        i, j = tup
        return self.cuadri[i][j]

    def __setitem__(self, tup, value):
        # como setear cuadripolos
        i, j = tup
        self.cuadri[i][j] = value

    def __repr__(self):
        if self.tipo != 'complejo':
            return print('Cuadripolo {} con parametros sección, ' \
                         'sección de ingreso {}' \
                         'sección de salida {}' \
                         'sección media {} largo {}'.format(self.cuadri, self.s1, self.s2, self.s, self.largo))
        if self.tipo == 'complejo':
            return 'Cuadripolo de sistema multiple {}'.format(self.cuadri)

    def coeficientes(self):
        """función para obtener valores de cuadripolos para distintos tipos de silenciadores
        input:
            tpye= str. los tipos posibles: (tubo_recto: para un tubo_recto, camara: para camara de expansión, cambio: para cambio
            de seccion)"""
        # mando en función de f, cuando evaluo ? Me interesa el valor numerico
        # de los coeficientes
        # tenemos un problema con k, hay un k para cada caso.
        if self.tipo == 'helmholtz':
            k = (self.rho_o * (self.c ** 2) * (self.s ** 2)) / self.vol
        else:
            k = 2 * np.pi * self.f / self.c
        if (self.tipo == 'camara') | (self.tipo == 'tubo_recto'):
            try:
                A = cos(k * self.largo)
                B = (self.Z_o / self.s) * sin(k * self.largo)
                C = (self.s / self.Z_o) * sin(k * self.largo)
                D = cos(k * self.largo)
                self.cuadri = [[A, B, ], [C, D]]
            except Exception as e:
                print(e,
                      'se deben definir los valores adecuados para obtener los coeficientesn camara o tubo_recto, s, Z_o, largo, ')
        if (self.tipo == 'Z_in') | (self.tipo == 'helmholtz'):
            try:
                A = 1
                B = 0
                C = self.s / self.Z_in
                D = 1
                self.cuadri = [[A, B, ], [C, D]]
            except Exception as e:
                print(e, 'se deben definir los valores adecuados para helmholtz o algun Z_in')
        if self.tipo == 'extension_expansion':
            try:
                A = 1
                B = 0
                C = ((tan(k * self.largo) ** (-1)) * self.Z_o / (self.s1 - self.s)) ** (-1)
                D = 1
                self.cuadri = [[A, B, ], [C, D]]
            except Exception as e:
                print(e, 'deben estar definidos s, s1, Z_o y largo')
        if self.tipo == 'complejo':
            print('no se pueden calcular coeficientes de silenciadores multiples')
            print(self.cuadri)

    def obtencion_tl(self):
        '''Devuelve función de tl en frecuencia'''
        A = self.cuadri[0][0]
        B = self.cuadri[0][1]
        C = self.cuadri[1][0]
        D = self.cuadri[1][1]
        if self.s2:
            try:
                self.tl = 20 * log(abs(
                    0.5 * (A + (self.Z_o * C / self.s2) + (B * self.s1 / self.Z_o) + (
                            D * self.s1 / self.s2)))) + 10 * log(
                    self.s1 / self.s2)
                self.tl_ = lambdify(self.f, self.tl, ["numpy"])
            except ValueError:
                print('deben ser definidos los parametros correctos para obtener tl s1 y s2')
            except Exception as er:
                print(er)
        else:
            try:
                self.tl = 20 * log(abs(0.5 * (A + (self.Z_o * C / self.s1) + (B * self.s1 / self.Z_o) + D)))
                self.tl_ = lambdify(self.f, self.tl, ["numpy"])
            except Exception as er:
                print(er)

    def __mul__(self, other):
        # self por izquierda. Recordar que son todas expresiones sympy
        A = self.cuadri[0][0] * other.cuadri[0][0] + self.cuadri[0][1] * other.cuadri[1][0]
        B = self.cuadri[0][0] * other.cuadri[0][1] + self.cuadri[0][1] * other.cuadri[1][1]
        C = self.cuadri[1][0] * other.cuadri[0][0] + self.cuadri[1][1] * other.cuadri[1][0]
        D = self.cuadri[1][0] * other.cuadri[0][1] + self.cuadri[1][1] * other.cuadri[1][1]
        c = Cuadripolo()
        c.cuadri = [[A, B], [C, D]]
        c.tipo = 'complejo'
        return c

    def __rmul__(self, other):
        # self por derecha. Recordar que son todas expresiones sympy
        A = other.cuadri[0][0] * self.cuadri[0][0] + other.cuadri[0][1] * self.cuadri[1][0]
        B = other.cuadri[0][0] * self.cuadri[0][1] + other.cuadri[0][1] * self.cuadri[1][1]
        C = other.cuadri[1][0] * self.cuadri[0][0] + other.cuadri[1][1] * self.cuadri[1][0]
        D = other.cuadri[1][0] * self.cuadri[0][1] + other.cuadri[1][1] * self.cuadri[1][1]
        c = Cuadripolo()
        c.cuadri = [[A, B], [C, D]]
        c.tipo = 'complejo'
        return c

    def plot_tl(self, values=np.linspace(0, 10000, 10000)):
        if values is not None:
            data = pd.DataFrame({'TL [dbs]':self.tl_(values), 'Frecuencia [Hz]':values})
            sns.lineplot(x='Frecuencia [Hz]', y='TL [dbs]', data=data)
        else:
            plot_ = plotting.plot(self.tl, range=(self.f, 0, 150))
            plot_.show()
