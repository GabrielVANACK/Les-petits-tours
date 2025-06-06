import math as math
import cmath as cmath

##Fonctions usuels

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

def Poly_fun(P):
    """prend la liste des coefficients en entrée en renvoie la fonction polynomiale"""
    n =len(P)
    
    def Pol(a):
        S=0
        for i in range(n):
            S += P[i]*(a**i)
        return S
    return Pol

#On représente les matrices par un tableau de leur coefficients

def Prod_mat(A,B):
    """ prend les matrices A,B en entrée et renvoie le produit matrcielle AB, print un message d'erreur dans le terminal si toutefois les tailles de matrices ne sont pas compatibles """
    (n,p) = (len(A),len(A[0]))
    (q,m) = (len(B), len(B[0]))
    if p==q: 
        M = [[0 for j in range(m)] for i in range(n)]
        for i in range(n):
            for j in range(m):
                S=0
                for k in range(p):
                    S+= A[i][k]*B[k][j]
                M[i][j] = S
        return M
    else : raise ValueError("les tailles des matrices ne sont pas compatibles")

def Inversion_mat(A):
    #BY Chat GPT
    """Inverse une matrice carrée A en utilisant la méthode de Gauss-Jordan."""
    n = len(A)
    # Vérification si la matrice est carrée
    if any(len(row) != n for row in A):
        raise ValueError("La matrice doit être carrée pour être inversée.")
    
    # Création de la matrice augmentée [A | I]
    I = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    tab = [A[i] + I[i] for i in range(n)]
    
    # Réduction de la matrice augmentée
    for i in range(n):
        # Recherche du pivot
        pivot = tab[i][i]
        if pivot == 0:
            raise ValueError("La matrice ne peut pas être inversée.")
        
        # Normalisation de la ligne du pivot
        tab[i] = [x / pivot for x in tab[i]]
        
        # Élimination des autres lignes
        for j in range(n):
            if i != j:
                factor = tab[j][i]
                tab[j] = [tab[j][k] - factor * tab[i][k] for k in range(2 * n)]
    
    # Extraction de la matrice inverse (partie droite de la matrice augmentée)
    inverse = [row[n:] for row in tab]
    return inverse
    
##Interpolation de lagrange

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

##Interpolation par splines

#A tester
def Interpolation_splines(l):
    """prends une liste de points du plan complexe en entrée et renvoie sa fonction d'interpolation par des splines"""
    #Pour l'instant c'est du code retranscrit de https://fr.wikipedia.org/wiki/Spline
    n = len(l)
    h = [l[i+1][0]-l[i][0] for i in range(0,n-1)]
    F = [[0]] + [[(l[i+1][1]-l[i][1])/(h[i]) - (l[i][0]-l[i-1][0])/(h[i-1])] for i in range(1,n-1)]+[[0]]
    R = [[0 for i in range(n)] for j in range(n)]
    R[0][0]=1
    R[n-1][n-1]=1
    for i in range(1,n-1):
        R[i][i]=(h[i-1]+h[i])/3
        R[i][i+1]=h[i]/6
        R[i][i-1]=h[i-1]/6
    M = Prod_mat(Inversion_mat(R),F)
    C1 = [0 for _ in range(n-1)]
    C2 = [0 for _ in range(n-1)]
    for i in range(n-1):
        C1[i]= (l[i+1][1]-l[i][1])/h[i] -h[i]*(M[i+1][0]-M[i][0])/6
        C2[i]=l[i][1]-(M[i][0]*(h[i])**2)/2

    def fun(x):
        if l[0][0]<=x<=l[n-1][0]:
            for i in range(n-1):
                if l[i+1][0]>=x>=l[i][0]:
                    return (M[i][0]/(6*h[i]))*(l[i+1][0]-x)**3+(M[i+1][0]/(6*h[i]))*(x-l[i][0])**3 + C1[i]*(x-l[i][0]) +C2[i]

        else : return 0 #Dans le cadre du problème il me semble faire sens d'attribuer à une fonction la valeur lorsqu'elle n'est pas défini (pour garantir le bon fonctionnement de C_graph)
    return fun
    
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


def Graph_C(f,view,Inter):
    """prend en entrée un fonction f, view qui est la vue du plan complexe disponible sur le graph renvoyée et Inter un intervalle centrée en 0 sur lequel f est calculé; renvoie un graph. """
    #En premier lieu on construit un intervalle de points
    I= np.linspace(10*view,-10*view,num= 100*view)
    a = len(I)

    #On initialise I2 qui est Inter mais étiré, pour que les points ne soient pas trop espacé ou trop rapproché
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


    
    