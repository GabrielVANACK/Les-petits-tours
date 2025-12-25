import numpy as np
import scipy.integrate as integrate

###Généralités

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

## à Optimiser

def Phi(P = np.array):
    """fonctions rayon moyen"""
    def A(x):
        F=fun_p(P)
        return module(F(x))
    
    return integ(A,0,1)

def dPhi(P = np.array, k=float):
    """renvoie la différentiel selon e_k de Phi de P"""
    p = fun_p(P) 
    if k> 2*deg+1 or k<0:
        raise ValueError("indice k trop grand ou trop petit")
    
    elif k%2 == 0 :
        ex = k//2

        def F(x):
            return (x**ex)*((p(x)).real)/(module(p(x)))
        
    elif k%2 == 1 :
        ex = k//2 

        def F(x):
            return (x**ex)*((p(x)).imag)/(module(p(x)))
    return integ(F,0,1)

def d2Phi(P = np.array,k=float,l=float):
    """renvoie la différentiel selon e_k et e_l de Phi de P"""

    p = fun_p(P) 
    if (k> 2*deg +1 or k<0) or (l> 2*deg +1 or l<0) :
        raise ValueError("indice k ou l trop grand ou trop petit")
    
    elif k%2==0 and l%2 ==0: 
        ex = k//2 + l//2
        def F(x):
            return -(x**(ex))*((p(x).imag)**2)/((module(p(x)))**3)
        
    elif k%2 == 1 and l%2 == 1:
        ex = k//2 + l//2
        def F(x):
            return -(x**(ex))*((p(x).real)**2)/((module(p(x)))**3)

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
    Lag = np.zeros((2*deg+4,1))

    for i in range(2*deg+2):
        Lag[i][0] = dPhi(P,i) - lambda1*dPsi(P,i) - lambda2*dOmega(P,i)

    Lag[2*deg + 2][0] = Psi(P)
    Lag[2*deg + 3][0] = Omega(P)

    return Lag

def Jacob(P = np.array,lambda1=float,lambda2=float):

    Jacobienne = np.zeros((2*deg+4,2*deg+4))

    for i in range(2*deg+2):
        for j in range(2*deg+2):
            Jacobienne[i][j] = d2Phi(P,i,j) - lambda1*d2Psi(P,i,j) - lambda2*d2Omega(P,i,j)
    
    for i in range(2*deg+2):
        Jacobienne[i][2*deg+2] = dPsi(P,i)
        Jacobienne[2*deg+2][i] = dPsi(P,i)
        Jacobienne[i][2*deg+3] = dOmega(P,i)
        Jacobienne[2*deg+3][i] = dOmega(P,i)
    
    #print(Jacobienne)
    return Jacobienne 

### Résolution

##Algorithme de Raphson Newton , p.24 (pdf)

def Raphson_Newton(P0 = np.array, lambda1 = float, lambda2 = float , it = 50 ):
    X = np.zeros((2*deg+4,1))
    for i in range(deg+1):
        X[2*i][0]= (P0[i]).real
        X[2*i+1][0]= (P0[i]).imag
    X[2*deg+2][0] = lambda1
    X[2*deg+3][0] = lambda2

    for i in range(it):
        J =Jacob(P0,lambda1,lambda2)
        deltaX = np.linalg.solve(matadd(J, scal((1/np.linalg.det(J))**6,np.identity(2*deg+4))),scal(-1,Lagrangien(P0,lambda1,lambda2)))
        X = X + deltaX

        for i in range(deg+1):
            P0[i] = complex(X[2*i][0] , X[2*i+1][0])
        lambda1 = X[2*deg+2][0]
        lambda2 = X[2*deg+3][0]
        print(satis_facteur(P0))
    return (P0,lambda1,lambda2)

#seed obtenue par interpolation  d'une solution discrète
#[complex(0.1,0.1),complex(35.1,140),complex(-1113,-2660),complex(10447,18904),complex(-45420,-67083),complex(105176,131527),complex(-133621,-144907),complex(87797,84052),complex(-23301,-19973)]
#fonctionne pour, deg =8
 
def satis_facteur(P):
    print(f"Cyclcité : \n distance entre P(1) et P(0) est {np.sqrt(Omega(P))}")
    print(f"Longueur : \n La longueur voulue est {longueur_du_parcours}, le parcours fait {Psi(P)+longueur_du_parcours} soit un différence de {Psi(P)}")
    print(f"Rayon moyen : \n le rayon moyen à l'origine du plan complexe de P sur [0,1] est {Phi(P)}")


## Autre méthode : page 114 du livre springer, page 127 pour le pdf

def T(X):
    (lambda2, lambda1) = (0,0)
    P = np.ndarray.flatten(np.full((1,deg+1), 0+0j)) 

    for i in range(deg+1):
        P[i] = complex(X[2*i][0] , X[2*i+1][0])
        lambda1 = X[2*deg+2][0]
        lambda2 = X[2*deg+3][0]
    
    return scal(1/2,matprod(np.transpose(Lagrangien(P, lambda1 , lambda2)), Lagrangien(P, lambda1 , lambda2)))

def find_steplength(x,dx ,s = float ,ls = list ,k=0.5 ):
    l=len(ls)
    m= T(x)

    if ls[l-2]<= ls[l-1]:
        s= min(max(ls), s/k)

    while T(matadd(x,scal(s,dx)))> m :
        s = k*s 
    ls.append(s)
    return s



def Steepest_descent(P = np.array,lambda1=float,lambda2=float, it = 50): 

    #Initialisation
    X = np.zeros((2*deg+4,1))
    for i in range(deg+1):
        X[2*i][0]= (P[i]).real
        X[2*i+1][0]= (P[i]).imag
    X[2*deg+2][0] = lambda1
    X[2*deg+3][0] = lambda2
    s=10
    ls = [s] #list des s

    #Itération
    for i in range(it):
        deltaX =scal(-1,matprod(np.transpose(Jacob(P, lambda1 , lambda2)), Lagrangien(P, lambda1 , lambda2)))
        s = find_steplength(X,deltaX,s,ls)
        print(s,ls)
        X = matadd(X,scal(s,deltaX))

        for i in range(deg+1):
            P[i] = complex(X[2*i][0] , X[2*i+1][0])
            lambda1 = X[2*deg+2][0]
            lambda2 = X[2*deg+3][0]
        print(satis_facteur(P))
    return (P,lambda1,lambda2)