import numpy as np
import scipy.integrate as integrate

###Généralités

def liste_shift(l,el):
    """prend un élèment el en entrée, l'ajoute au bout de la liste et supprime le premier élèment de la liste. Plante pour des listes de taille 1"""

    n = len(l)
    for i in range(n-1):
        l[i]=l[i+1]
    l[n-1]= el

def polprod(P = np.linspace,Q=np.array):
    """réalise le produit de cauchy des polynômes P et Q sous forme de liste de coefficients"""
    n=len(P)
    m=len(Q)
    PQ= np.ndarray.flatten(np.full((1,(n+m-1)), 0+0j)) 
    for i in range(n):
        for j in range(m):
            PQ[i+j] += P[i]*Q[j]

    return PQ

def poladd(P = np.array,Q=np.array):
    """réalise l'adition des polynômes P et Q représenté par des listes"""
    n=len(P)
    m=len(Q)
    p = max(n,m)
    l =min(n,m)
    PplusQ= np.ndarray.flatten(np.full((1,p), 0+0j)) 
    for i in range(l):
        PplusQ[i] = P[i] + Q[i]
    for i in range(l,p):
        if i >= n :
            PplusQ[i] = Q[i]
        else :
            PplusQ[i] = P[i]

    return PplusQ

def fun_p(P = np.array):
    """renvoie la fonction polynomiale du polynôme définie pas la liste de ses coefficients P"""
    n = len(P)
    def poly_fun(x):
        res = 0
        for k in range(n):
            res += (x**k)*P[k]
        return res
    return poly_fun

def dp(P = np.array):
    """prend un polynome en entrée et renvoie sa dérivé"""
    n= len(P)
    dP = np.ndarray.flatten(np.full((1,n), 0+0j))
    for k in range(1,n):
        dP[k-1]=k*P[k]
    return dP

def matprod(A=np.array,B=np.array):
    """renvoie le produit matricielle de A par B"""
    (n,m) = np.shape(A)
    (p,l) = np.shape(B)
    if p != m :
        raise ValueError ("Matrices de tailles non compatibles")
    else :
        return np.matmul(A,B)

def matadd(A=np.array,B=np.array):
    """réalise l'addtion de la matrice A et B"""
    if A.shape == B.shape :
        return A+B
    else : raise ValueError("les matrices ne sont pas de même taille, impossible de les addtionner")

def scal(mu=float,A=np.array):
    """multiplication de la matrice A par un scalaire"""
    return np.multiply(A,mu)

def matinv(A=np.array):
    print(np.linalg.det(A))
    return np.linalg.inv(A)

def module(z):
    """renvoie le module de z"""
    return np.abs(z)

def integ(f,a=float,b=float):
    """réalise l'intégrale de f sur a,b avec la méthode de simpson, par défaut it (="itération")=1000"""
    return integrate.quad(lambda x:f(x), a, b)[0]
    
### Fonctions du problème

deg = 8 # degré maximal des polynomes étudidés
longueur_du_parcours = 20 # Longueur du parcours étudié
reg = 3 #correction pour assurer que Phi soit coercif

## à Optimiser

def Phi(P = np.array):
    """fonctions rayon moyen"""
    def A(x):
        F=fun_p(P)
        d2F = fun_p(dp(dp(P)))
        return module(F(x)) + reg*(module(d2F(x)))**2
    
    return integ(A,0,1)

def dPhi(P = np.array, k=float):
    """renvoie la différentiel selon e_k de Phi de P"""
    p = fun_p(P)
    d2p =fun_p(dp(dp(P)))
    if k> 2*deg+1 or k<0:
        raise ValueError("indice k trop grand ou trop petit")
    
    elif k%2 == 0 :
        ex = k//2

        def F(x):
            return (x**ex)*((p(x)).real)/(module(p(x))) + reg*(ex-1)*ex*(x**(ex-2))*(d2p(x).real)
        
    elif k%2 == 1 :
        ex = k//2 

        def F(x):
            return (x**ex)*((p(x)).imag)/(module(p(x)))  + reg*(ex-1)*ex*(x**(ex-2))*(d2p(x).imag)
    return integ(F,0,1)

def d2Phi(P = np.array,k=float,l=float):
    """renvoie la différentiel selon e_k et e_l de Phi de P"""
    p = fun_p(P) 
    if (k> 2*deg +1 or k<0) or (l> 2*deg +1 or l<0) :
        raise ValueError("indice k ou l trop grand ou trop petit")
    
    elif k%2==0 and l%2 ==0: 
        ex = k//2 + l//2
        def F(x):
            return -(x**(ex))*((p(x).imag)**2)/((module(p(x)))**3) + reg*((ex-1)*ex*(x**(ex-2)))**2
        
    elif k%2 == 1 and l%2 == 1:
        ex = k//2 + l//2
        def F(x):
            return -(x**(ex))*((p(x).real)**2)/((module(p(x)))**3) + reg*((ex-1)*ex*(x**(ex-2)))**2

    elif (k%2==0 and l%2==1) or (k%2==1 and l%2==0): 
        ex = k//2 + l//2
        def F(x):
            return -(x**(ex))*((p(x).real)*(p(x).imag))/((module(p(x)))**3)  
    return integ(F,0,1)     

## Contraintes

#Contrainte de longueur

def Psi(P = np.array):
    def A(x):
        F = fun_p(dp(P))
        return module(F(x))
    return integ(A,0,1)-longueur_du_parcours

def dPsi(P = np.array,k=float):
    A = fun_p(dp(P)) 

    if k > 2*deg+1 or k<0 :
        raise ValueError("indice k trop grand ou trop petit")

    elif k==0 or k==1:
        return 0
    
    elif k%2==0 :
        ex = k//2
        def F(x):
            return ex*(x**(ex-1))*((A(x)).real)/(module(A(x)))
        
    elif k%2==1:
        ex = k//2
        def F(x):
            return (ex)*(x**(ex-1))*((A(x)).imag)/(module(A(x)))
        
    return integ(F,0,1)

def d2Psi(P = np.array,l=float,k=float):
    A = fun_p(dp(P))

    if k>2*deg+1 or k<0 or  l>2*deg+1 or l<0:
        raise ValueError("indice k ou l trop grand ou trop petit")
    
    elif k==0 or k==1 or l==0 or l==1:
        return 0
    
    elif k%2==0 and l%2==0: 
        ex = (k//2) + (l//2)
        facteur = (k//2) * (l//2)
        def F(x):
            return -(facteur)*(x**(ex-2))*((A(x).imag)**2)/((module(A(x)))**3)
        
    elif k%2==1 and l%2==1:
        ex = (k//2) + (l//2)
        facteur = (k//2) * (l//2 )
        def F(x):
            return -(facteur)*(x**(ex-2))*((A(x).real)**2)/((module(A(x)))**3)
        
    elif (k%2==0 and l%2==1) or (k%2==1 and l%2==0): 
        ex = (k//2) + (l//2)
        facteur = (k//2) * (l//2 )
        def F(x):
            return -(facteur)*(x**(ex-2))*((A(x).imag)*(A(x).real))/((module(A(x)))**3)
       
    return integ(F,0,1)

#Contrainte de cyclicité 

def Omega(P = np.array):
    p = fun_p(P)
    return (module(p(1)-p(0)))**2

def dOmega(P = np.array ,k = float):
    """renvoie la différentiel selon e_k de Omega de P"""
    p = fun_p(P)
    if k > 2*deg+1 or k<0 :
        raise ValueError("k est plus gand ou plus petit que le degré des polynôme (variable deg)")
    
    if k == 0 or k==1 :
        return 0
    
    elif k%2==0 :
        return 2*(p(1)-p(0)).real 
    
    elif k%2==1:
        return 2*(p(1)-p(0)).imag 

def d2Omega(P = np.array,l = float ,k = float):
    """renvoie la différentiel selon e_k et e_l de Omega de P"""
    if k>2*deg+1 or k<0 or  l>2*deg+1 or l<0:
        raise ValueError("indice k ou l trop grand ou trop petit")
    
    if k == 0 or k==1 or l ==0 or l==1:
        return 0
    elif (k%2 == 1 and l%2 ==0) or (k%2 == 0 and l%2 ==1):
        return 0
    else:
        return 2
    
## Matrices et fonction Lagrangienne

def Lagrangien(P = np.array ,lambda1 = float ,lambda2 = float):
    Lag = np.zeros((2*deg,1))

    for i in range(4,2*deg+2):
        Lag[i-4][0] = dPhi(P,i) - lambda1*dPsi(P,i) - lambda2*dOmega(P,i)

    Lag[2*deg-2][0] = Psi(P)
    Lag[2*deg -1][0] = Omega(P)

    return Lag

def Jacob(P = np.array,lambda1=float,lambda2=float):

    Jacobienne = np.zeros((2*deg,2*deg))

    for i in range(4,2*deg):
        for j in range(4,2*deg):
            Jacobienne[i-4][j-4] = d2Phi(P,i,j) - lambda1*d2Psi(P,i,j) - lambda2*d2Omega(P,i,j)
    
    for i in range(2*deg):
        Jacobienne[i][2*deg-2] = dPsi(P,i)
        Jacobienne[2*deg-2][i] = dPsi(P,i)
        Jacobienne[i][2*deg-1] = dOmega(P,i)
        Jacobienne[2*deg-1][i] = dOmega(P,i)
    
    #print(Jacobienne)
    return Jacobienne 

### Résolution

##Algorithme de Raphson Newton , p.24 (pdf)
mu =10**(-3)

def Raphson_Newton(P0 = np.array, lambda1 = float, lambda2 = float , it = 50 ):
    X = np.zeros((2*deg,1))
    for i in range(deg-1):
        X[2*i][0]= (P0[i+2]).real
        X[2*i+1][0]= (P0[i+2]).imag
    X[2*deg-2][0] = lambda1
    X[2*deg-1][0] = lambda2

    global mu

    for i in range(it):
        J =Jacob(P0,lambda1,lambda2)
        Jt = np.transpose(J)
        L = scal(-1,Lagrangien(P0,lambda1,lambda2))

        B = matprod(Jt,L)
        A = matprod(Jt,J) + mu*np.diag(np.diag(matprod(Jt, J)))
        deltaX = np.linalg.solve(A,B)

           
        if T(X+deltaX)<T(X):
            mu = mu/(10)
            X = X + deltaX

        elif mu> 10**(300):
            mu=10**(-100)
            print("mu devenu trop grand remise dans des limites raisonnable")
        else : mu = mu*(10) 
        #print("mu",mu)
        #print(deltaX)

        for i in range(deg-1):
            P0[i+2] = complex(X[2*i][0] , X[2*i+1][0])
        P0[0]= complex(0,0)
        P0[1]= complex(0,0)
        lambda1 = X[2*deg-2][0]
        lambda2 = X[2*deg-1][0]
        #print(P0)
        #print(satis_facteur(P0))
        #print("||Lag||",np.linalg.norm(Lagrangien(P0,lambda1,lambda2),2))
    return (P0,lambda1,lambda2)

def cv_reg(P0 = np.array, lambda1 = float, lambda2 = float):
    # Pour améliorer la convergence je vérifie que les itérations ne sont pas trop proche, bien sûr ça suppose qu'il n'y ait pas de "plateau" dans la convergence, mais bon entre ça et faire 20 itérations bêtement je sais pas quoi faire
    global reg
    i=0
    memoir_de_convergence = [0 if i == 0 else 1 for i in range(5)]
    while abs(memoir_de_convergence[4] - memoir_de_convergence[0]) > 10**(-5) or reg > 10**(-9):
        global mu 
        mu = 10**(-3)
        Lag_norm_mem = [0,1,float('inf')]
        trac = 0
        while abs(Lag_norm_mem[2] - Lag_norm_mem[0]) > 10**(-5):
            (P0,lambda1,lambda2) = Raphson_Newton(P0,lambda1,lambda2,35)
            trac +=35
            liste_shift(Lag_norm_mem, np.linalg.norm(Lagrangien(P0,lambda1,lambda2),2))

            print(f"Avec le facteur {np.round(reg,3)/(10**(int(np.log(reg)/(np.log(10)))))}x 10^{ (int(np.log(reg)/(np.log(10))))}, {trac} itérations :  ||Lag||",Lag_norm_mem[1])
        print(f" Fin de boucle {i}-eme) ||Lag|| ",Lag_norm_mem[1],"\n")

        liste_shift(memoir_de_convergence, np.linalg.norm(Lagrangien(P0,lambda1,lambda2),2))

        reg = reg/(10 + 1/(10**(-30)+abs(memoir_de_convergence[4] - memoir_de_convergence[3]))) #Diminue plus vite si la liste des lagrangiens sont très proches
        i+=1

    return (P0,lambda1,lambda2)

def cv_continue(P0 = np.array, lambda1 = float, lambda2 = float , deg_max = int):   
    global deg, reg
    original_deg = deg
    for i in range(deg, deg_max):
        deg = i
        print(f"\n \n Degré des polynômes étudiés : {deg}\n")
        reg = 10
        (P0, lambda1, lambda2) = cv_reg(P0, lambda1, lambda2)
        
        P0.append(complex(0,0))

    deg = original_deg
    return (P0, lambda1, lambda2)

#seed obtenue par interpolation  d'une solution discrète
seed=[complex(0.1,0.1),complex(35.1,140),complex(-1113,-2660),complex(10447,18904),complex(-45420,-67083),complex(105176,131527),complex(-133621,-144907),complex(87797,84052),complex(-23301,-19973)]
#fonctionne pour, deg =8
 
def satis_facteur(P):
    print(f"Cyclcité : \n distance entre P(1) et P(0) est {np.sqrt(Omega(P))}")
    print(f"Longueur : \n La longueur voulue est {longueur_du_parcours}, le parcours fait {Psi(P)+longueur_du_parcours} soit un différence de {Psi(P)}")
    print(f"Rayon moyen : \n le rayon moyen à l'origine du plan complexe de P sur [0,1] est {Phi(P)}")


## Autre méthode : page 114 du livre springer, page 127 pour le pdf

def T(X):
    (lambda2, lambda1) = (0,0)
    P = np.ndarray.flatten(np.full((1,deg+1), 0+0j)) 

    for i in range(deg-1):
        P[i+2] = complex(X[2*i][0] , X[2*i+1][0])
        lambda1 = X[2*deg-2][0]
        lambda2 = X[2*deg-1][0]

    L  = Lagrangien(P, lambda1 , lambda2)
    Lt = np.transpose(Lagrangien(P, lambda1 , lambda2))
    
    return scal(1/2,matprod(Lt,L)) 

