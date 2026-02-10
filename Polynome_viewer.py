import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

#Pour un polynôme P de la variable complexe on renvoie son graphe dans R2 sur [0,1] 

def Graph_C(f):
    I = np.linspace(0, 1, 1000)
    x = np.array([f(t).real for t in I])
    y = np.array([f(t).imag for t in I])

    plt.plot(x, y, color="black")
    plt.scatter(x[0], y[0], color="red", label="départ")
    plt.scatter(x[-1], y[-1], color="blue", label="arrivée")

    plt.axis("equal")
    plt.xlabel("partie réelle")
    plt.ylabel("partie imaginaire")
    plt.legend()
    plt.title("Trajectoire paramétrée")
    plt.show()

import matplotlib.pyplot as plt

def plot_values(values):
    # Préparer les données d'index
    x = [i for  i in range(len(values))]
    y = values

    # Création de la figure
    plt.figure(figsize=(8, 4))
    plt.plot(x, y, marker="o", linestyle="-")
    plt.scatter(x, y, s=40)  # met en évidence les points

    # Titres et grille
    plt.title("évolution de l'inversion du déterminant")
    plt.xlabel("itération")
    plt.ylabel("1/det(A)")
    plt.grid(True, linestyle="--", alpha=0.6)

    # Affichage des indices comme ticks entiers
    plt.xticks(x)
    plt.show()

def fun_p(P = np.array):
    """renvoie la fonction polynomiale du polynôme définie pas la liste de ses coefficients P"""
    n = len(P)
    def poly_fun(x):
        res = 0
        for k in range(n):
            res += (x**k)*P[k]
        return res
    return poly_fun


def tuple_to_complex(t):
    return complex(t[0], t[1])

def l_to_pol(l):
    """Transforme uns liste de liste de deux élèments en polynôme complexe"""
    P = [tuple_to_complex(z) for z in l]
    return P

def conv(P):
    """Convertit les polynômes sous forme liste de liste (comme dans les fichiers json) en polynômes"""
    return fun_p(l_to_pol(P))