import math as math
import cmath as cmath

#Dans un premier temps on souhaite tracer la courbe par interpolation de Lagrange, pour éviter les figures absurdes on souhaite que la dérivée de la courbe de tracée soit bornée entre des valeurs connu au préalable.

def epsi(a,obj,e):
    """prend en entrée une valeur, un objectif et un réel epsilon et renvoie si la valeur est epsilon proche de la valeur objectif"""
    if obj-e<=a and a<=obj+e:
        return True
    else : return False

def module(z):
    return math.sqrt(z.real**2 + z.imag**2)

def D(f,a):
    """calcul la dérivée de la fonction f en a"""
    h= 1*10**(-6)
    S = []
    flag = True
    while flag:
        S.append((f(a+h)-f(a-h))/(2*h))
        s = len(S)
        if s>=2:
            if epsi(module(S[s-1]),module(S[s-2]),10**(-9)):
                flag = False
        h = h/10
    return S[s-1]

#On représente le polynômes par la liste de leurs coefficients
def Prod_poly(P,Q):
    """fait le produit de deux polynômes P et Q de degré quelconques"""
    n = len(P) + len(Q) - 1
    result = [0] * n
    for i in range(len(P)):
        for j in range(len(Q)):
            result[i + j] += P[i] * Q[j]
    return result

def Poly_de_Lag(l):
    """définit les polynômes de lagrange d'une liste l de points"""
    n = len(l)
    L = []
    for i in range(n):
        P=[1]
        for M in l:
            if M[0]!=l[i][0]:
                P = Prod_poly(P,[(-M[0])/(l[i][0]-M[0]),1/(l[i][0]-M[0])])
        L.append(P)  
    return L

def Interpolation_Lagrange(l):
    """définit le polynôme interpolateur de lagrange d'une liste l de points"""
    PL = Poly_de_Lag(l)
    n = len(PL)
    IL = []
    for i in range(n):
        a=0
        for j in range(n):
            a += l[j][1]*PL[j][i]
        IL.append(a)
    return IL

def Poly_fun(P):
    """prend la liste des coefficients en entrée en renvoie la fonction polynomiale"""
    n =len(P)
    
    def Pol(a):
        S=0
        for i in range(n):
            S += P[i]*(a**i)
        return S
    return Pol
    
            
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


def Graph_C(f,view,Inter):
    #En premier lieu on construit un intervalle de points
    I= np.linspace(10*view,-10*view,num= 100*view)
    a = len(I)

    #On initialise I2 qui est I mais étiré en fonction de la dérivé de f, pour que les points ne soient pas trop espacé ou trop rapproché
    I2 = np.linspace(Inter,-Inter,num= 10*Inter*view)
    b=len(I2)
    #Ensuite on calcule les valeurs de la fonction au valeurs enregistrée dans I'
    fun = [f(i) for i in I2]
    im= [v.imag for v in fun ]
    re = [v.real for v in fun ]

    #On initialise un tableau de valeur dans lequel on lequel on place l'indice i de f(i) au coordonnées im[i] et re[i]
    Z = [[0 for _ in range(a)] for _ in range(a)]
    for i in range(b) :
        indice_im = 0
        indice_re = 0
        flag1 = False
        flag2 = False
        for j in range(a):
            if re[i]>=I[j] :
                indice_re = j
                flag1 = True
                break
        for l in range(a):
            if im[i]>=I[l] :
                indice_im = l
                flag2=True
                break
        if flag1 and flag2:
            Z[indice_im][indice_re] =(i/(10*view)) 

    #On créer une color map qui convient
    viridis = plt.get_cmap('viridis', 256)
    viridis_colors = viridis(np.linspace(0, 1, 256))
    custom_colors = np.vstack(([1, 1, 1, 1], viridis_colors[1:]))
    custom_cmap = ListedColormap(custom_colors)

    #On génére le graph
    plt.imshow(Z, extent=[-view,view,-view,view],origin="lower", cmap=custom_cmap, aspect="auto")
    plt.colorbar(label="sens de parcours")
    plt.xlabel("partie réel")
    plt.ylabel("partie imaginaire")
    plt.title("graphe de la fonction f dans le plan")
    plt.show ()


    
    