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



#Meilleur Polynôme trouvé (degré 20) :
#P  =[[0.0,16],[0.0,-37],[-1111.256330664577,-2663.356239437416],[10446.48830826279,18936.326855807994],[-45421.92465239192,-67052.60038268479],[105173.00828515005,131554.83789995653],[-133624.26664757266,-144882.4240209144],[87794.15874118521,84072.40202928931],[-23276.43907705625,-19934.94057756686],[-8.880868670937435,0.2826354785487945],[39.05912488370479,6.700018259123785],[0.10961278411125934,0.5913301352990398],[-1.28431390287852,0.09640714378281563],[-2.3819361224242495,-0.2266614570853096],[-1.8415547681262163,-0.11557552639850564],[-1.625228806631387,-0.10919412133468381],[-0.53498079543762,0.010964288628949138],[-0.690750383795395,-0.012339143266521071],[-0.44815518309307883,-0.006514043912345639],[0.0019714373997672782,0.01065799812181316],[-0.6116666866382765,-0.01873499341067836],[0.0,0.0],[0.0,0.0] ]