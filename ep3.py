# -*- coding: latin-1 -*-#

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

class Magnetico:
    def __init__(self,a,b,c,d,g,h,delta,interno,externo,mi_r):
        #Dimensões
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.g = g
        self.h = h
        self.delta = delta

        #Potenciais e Constantes
        self.A_interno = interno
        self.A_externo = externo
        self.mi = mi_r * 1.2566 * 10**(-6)

        #Matrizes
        self.celulas = None
        self.nodos = None

        self.criarCelulasENodos()

    def criarCelulasENodos(self):
        #Colunas e linhas dos nós
        self.n_col = int(self.a/self.delta) + 1
        self.n_lin = int(self.b/self.delta) + 1

        #Cria as matrizes
        self.nodos = [[0.0 for _ in range(self.n_col)] for _ in range(self.n_lin)]

        #Inicializa valores dos potenciais
        self.criaBordaRetangular(self.nodos,self.A_externo, 0.0, 0.0, self.a, self.b)
        self.criaBordaRetangular(self.nodos,self.A_interno, self.g, self.h, self.g+self.c,
                            self.h+self.d)


    def criaBordaRetangular(self,matriz,potencial,x0,y0,x1,y1):
        ix0,ix1,iy0,iy1 = (int(x/self.delta) for x in (x0,x1,y0,y1))
        print(ix0,ix1,iy0,iy1)

        for i in range(ix0,ix1+1):
            matriz[iy0][i] = str(potencial)
            matriz[iy1][i] = str(potencial)
        for j in range(iy0,iy1+1):
            matriz[j][ix0] = str(potencial)
            matriz[j][ix1] = str(potencial)

    def getValorIdx(self, i, j):
        return float(self.nodos[j][i])

    def calcularPotencialIdx(self, i, j):
        A1 = self.getValorIdx(i-1, j)
        A2 = self.getValorIdx(i+1, j)
        A3 = self.getValorIdx(i, j+1)
        A4 = self.getValorIdx(i, j-1)
        return (A1+A2+A3+A4)/4

    def calcularHCelulaIdx(self, i, j):
        A1 = self.getValorIdx(i, j)
        A2 = self.getValorIdx(i+1, j)
        A3 = self.getValorIdx(i, j+1)
        A4 = self.getValorIdx(i+1, j+1)
        Hx = (1/(self.mi*1000000))*(A1+A3-A2-A4)/(2*self.delta)
        Hy = (1/(self.mi*1000000))*(A1+A2-A3-A4)/(2*self.delta)
        return (Hx, Hy)

    def getCampoVetorialH(self):
        celulas = [[None for _ in range(self.n_col-1)] for _ in range(self.n_lin-1)]
        for i in range(self.n_col-1):
            for j in range(self.n_lin-1):
                celulas[j][i] = self.calcularHCelulaIdx(i,j)
        return celulas

    def calcularIndutancia(self):
        #Calcula a Corrente
        campoH = self.getCampoVetorialH()

        #transforma coordenadas em índices matriciais
        ix0,ix1,iy0,iy1 = (
            int(x/self.delta) for x in (self.g,
                                        self.g+self.c+self.delta,
                                        self.h,
                                        self.h+self.d+self.delta)
        )

        #Integral de H*dl
        integral = 0
        for i in range(ix0,ix1+1):
            integral +=  abs(campoH[iy0][i][0]) #  |Ex| borda sup
            integral += -abs(campoH[iy1][i][0]) # -|Ex| borda inf
        for j in range(iy0,iy1+1):
            integral +=  abs(campoH[j][ix0][1]) #  |Ey| borda esq
            integral += -abs(campoH[j][ix1][1]) # -|Ey| borda dir

        corrente = integral
        print("corrente",corrente)
        fluxo = self.A_interno/1000000 #potencial estava em uWb/m
        indutancia = fluxo/corrente
        return indutancia

    def iterarPotenciais(self,n):
        """ Itera pelos nodos, sem passar pelas bordas externas ou internas
        de modo a evitar a condição de fronteira, para recalcular o potencial
        em cada nodo. """

        for num in range(1,n+1):
            if num%100 == 0:
                print("iteração #",num)
            for i in range(1,self.n_col-1):
                for j in range(1,self.n_lin-1):
                    if isinstance(self.nodos[j][i],float):
                        self.nodos[j][i] = self.calcularPotencialIdx(i,j)
                    else:
                        continue

#CALCULOS
problema = Magnetico(0.11, 0.06, 0.04, 0.03, 0.02, 0.015, 0.005, 100, 0, 1)
iteracoes = 2000
last_indut = 0
for _ in range(50):
    problema.iterarPotenciais(100)
    campoH = problema.getCampoVetorialH()
    indut = problema.calcularIndutancia()
    delta = abs(indut-last_indut)
    print("indutancia",indut,"Delta",delta)
    print("\n")
    last_indut = indut
    if( delta == 0 ):
        break

#EQUIPOTENCIAIS
fig, ax = plt.subplots()
ax.set_title(str(iteracoes)+" Iterações")
ax.set_xticks(np.arange(problema.n_col))
ax.set_yticks(np.arange(problema.n_lin))
ax.set_xticklabels(
    ( round(i,3) if k%4 == 0 else "" for k,i in
     enumerate(np.arange(0,problema.a+problema.delta,problema.delta)) )
)
ax.set_yticklabels(
    ( round(i,3) if k%4 == 0 else "" for k,i in
     enumerate(np.arange(0,problema.b+problema.delta,problema.delta)) )
)
fig.tight_layout()
fixed = [[float(x) for x in row] for row in problema.nodos]
cs = plt.contour(fixed,
                linewidths=0.7,
                levels=np.arange(problema.A_externo,problema.A_interno,11),
                colors="black")
#apenas para debugar, há o plot do campo
#sp = plt.streamplot(np.arange(problema.n_col-1),np.arange(problema.n_lin-1),
#                    np.array([[x[0] for x in row] for row in campoH]),
#                    np.array([[x[1] for x in row] for row in campoH])
#                     )
fig.set_dpi(96)
fig.set_size_inches(1280/96,720/96)
#fig.savefig("ep3_potencial",dpi='figure')
plt.show()

#INDUTANCIA






