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

def integral(f,a,b,n=1000):
    """réalise l'intégrale de f sur a b avec la méthode des rectangle."""
    S=0
    for i in range(n):
        S+= (1/n)*(f(i*(b-a)/n))
    return S

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
            S += (P[i])*(a**i)
        return S
    return Pol

def Deriv_poly(P):
    """prend un polynome en entrée et renvoie sa dérivé"""
    n= len(P)
    dP = [0]*n
    for i in range(1,n):
        dP[i-1]=i*P[i]
    return dP

#On représente les matrices par un tableau de leur coefficients
def Mat_add(A,B):
    """addtion de matrice A et B"""
    (n,m) = (len(A),len(A[0]))
    if (n,m) != (len(B),len(B[0])):
        raise ValueError ("Matrice de taille différente, impossible d'additionner")
    else : 
        res = [[0]*m]*n
        for i in range(n):
            for j in range(m):
                res[i][j]= A[i][j]+B[i][j]
    return res

def scal_mat(lam,A):
    """prend un scalaire et un matrice en entrée et réalise l'addition"""
    n,m = (len(A),len(A[0]))
    res = [[0]*m]*n
    for i in range(n):
        for j in range(m):
            res[i][j]= lam*A[i][j]
    return res

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
def Interpolation_splines(l):
    """prends une liste de points du plan complexe en entrée et renvoie sa fonction d'interpolation par des splines"""
    #Pour l'instant c'est du code retranscrit de https://fr.wikipedia.org/wiki/Spline
    #Ce premier algo fait donne fonctions qui ne sont absolument pas continue ce qui pose évidemment problème
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
                    return (M[i][0]/(6*h[i]))*(l[i+1][0]-x)**3 + (M[i+1][0]/(6*h[i]))*(x-l[i][0])**3 + C1[i]*(x-l[i][0]) +C2[i]

        else : return 0 #Dans le cadre du problème il me semble faire sens d'attribuer à une fonction la valeur lorsqu'elle n'est pas défini (pour garantir le bon fonctionnement de C_graph)
    return fun

##Méthode du Lagrangien (Dans l'espace des polynômes, pour résolution dans R^2)
# Sauf contre indication toutes les fonctions évoquée dans cette partie sont des fonctions de C_n[X]-> R, où n est un entier
deg = 6 #Degré maximal des polynômes


def Phi(P):
    return integral(module(Poly_fun(P)),0,1)

def dPhi(P,k):
    """renvoie la différentiel selon e_k de Phi de P"""
    
    p = Poly_fun(P) 
    if k<=deg and k>=0:
        def F(x):
            return (x**k)*((p(x)).real)/(module(p(x)))
    elif k<=2*deg+1 and k>=deg+1:
        def F(x):
            return (x**(k-deg-1))*((p(x)).imag)/(module(p(x)))
    else :
        raise ValueError("indice k trop grand ou trop petit")
    return integral(F,0,1)

def d2Phi(P,k,l):
    """renvoie la différentiel selon e_k et e_l de Phi de P"""

    p = Poly_fun(P) 
    
    if 0<=k<=deg and 0<=l<=deg: 
        def F(x):
            return -(x**(k+l))*((p(x).imag)**2)/((module(p(x)))**3)
    elif deg+1<=k<=2*deg+1 and deg+1<=l<=2*deg+1: 
        def F(x):
            return -(x**(k+l-2*deg-2))*((p(x).real)**2)/((module(p(x)))**3)
    elif 0<=k<=deg and deg+1<=l<=2*deg+1: 
        def F(x):
            return -(x**(k+l-deg-1))*((p(x).imag)*(p(x).real))/((module(p(x)))**3)
    else :
        raise ValueError("indice k ou l trop grand ou trop petit")
        
    return integral(F,0,1)

def Omega(P):
    return P(1).real-P(0).real + P(1).imag - P(0).imag

def dOmega(P,k):
    """renvoie la différentiel selon e_k de Omega de P"""
    return 1

def d2Omega(P,l,k):
    """renvoie la différentiel selon e_k et e_l de Omega de P"""
    return 0

def Psi(P):
    A = Poly_fun(Deriv_poly(P))
    return integral(module(Poly_fun(A)),0,1)

def dPsi(P,k):
    A = Poly_fun(Deriv_poly(P)) 
    if k==0 or k==deg+1:
        return 0
    elif k<=deg and k>=1:
        def F(x):
            return k*(x**(k-1))*((A(x)).real)/(module(A(x)))
    elif k<=2*deg+1 and k>=deg+2:
        def F(x):
            return (k-deg-1)*(x**(k-deg-2))*((A(x)).imag)/(module(A(x)))
    else :
        raise ValueError("indice k trop grand ou trop petit")
    return integral(F,0,1)

def d2Psi(P,l,k):
    A = Poly_fun(Deriv_poly(P))
    if k==1 or k==deg+1 or l==1 or l==deg+1:
        return 0
    elif 1<=k<=deg and 1<=l<=deg: 
        def F(x):
            return -(k*l)*(x**(k+l-2))*((A(x).imag)**2)/((module(A(x)))**3)
    elif deg+2<=k<=2*deg+1 and deg+2<=l<=2*deg+1: 
        def F(x):
            return -(k-deg-1)*(l-deg-1)*(x**(k+l-2*deg-4))*((A(x).real)**2)/((module(p(x)))**3)
    elif 1<=k<=deg and deg+2<=l<=2*deg+1: 
        def F(x):
            return -(l-deg-1)*k*(x**(k+l-deg-2))*((A(x).imag)*(A(x).real))/((module(A(x)))**3)
    else :
        raise ValueError("indice k ou l trop grand ou trop petit")
        
    return integral(F,0,1)

def Lagrangien(P,lambda1,lambda2):
    """appeler cette fonction est un abus de langage, elle représente la matrice dont on parle dans la section "problème à résoudre " des slides"""
    F = [[]*(deg+3)]
    for i in range(deg+1):
        F[0][i] = dPhi(P,i) - lambda1*dOmega(P,i)-lambda2*dPsi(P,i)
    F[deg+1]= Omega(P)
    F[deg+2] = Psi(P)

def Jacob(P,lambda1,lambda2):
    """ Cette fonction défini la matrice jacobienne du lagrangien."""
    F =[[0]*(deg+3)]*(deg+3)
    for i in range(deg+1):
        for j in range(deg+1):
            F[i][j]= d2Phi(P,i,j) - lambda1*d2Omega(P,i,j)-lambda2*d2Psi(P,i,j)
    
    for i in range(deg+3):
        if i!= deg+1:
            F[deg+1][i]= dOmega(P,i)
            F[i][deg+1]= dOmega(P,i)
        else : F[deg+1][deg+1]=0

    for i in range(deg+3):
        if i!= deg+2:
            F[deg+1][i]= dPsi(P,i)
            F[i][deg+1]= dPsi(P,i)
        else : F[deg+2][deg+2]=0

    return F 

def Raph_Newton(X0, round=100):
    """version crash test de la méthode de Raphson Newton, on demande à ce que X0 = (P,lambda1,lmabda2)"""
    a = len(X0)

    for i in range(round):
        X0 = Mat_add(X0,scal_mat((-1),Prod_mat((Inversion_mat(Jacob(X0[0],X0[1],X0[2]))),Lagrangien(Jacob(X0[0],X0[1],X0[2])))))



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

