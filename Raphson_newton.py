import numpy as np

###Généralités

def polprod(P=np.array,Q=np.array):
    """réalise le produit de cauchy des polynômes P et Q sous forme de liste de coefficients"""
    n=len(P)
    m=len(Q)
    PQ = np.zeros(1,(n+m-1))
    for i in range(n):
        for j in range(m):
            PQ[i+j] += P[i]*Q[j]

    return PQ

def poladd(P=np.array,Q=np.array):
    """réalise l'adition des polynômes P et Q représenté par des listes"""
    n=len(P)
    m=len(Q)
    p = max(n,m)
    PplusQ= np.zeros((1,p))
    for i in range(p):
        PplusQ[i] = P[i] + Q[i]
    return PplusQ

def fun_p(P=np.array):
    """renvoie la fonction polynomiale du polynôme définie pas la liste de ses coefficients P"""
    n = len(P)
    def poly_fun(x):
        res = 0
        for i in range(n):
            res = P[i]*x**i
        return res
    return poly_fun

def dp(P=np.array):
    """prend un polynome en entrée et renvoie sa dérivé"""
    n= len(P)
    dP = np.zeros((1,n))
    for i in range(1,n):
        dP[i-1]=i*P[i]
    return dP

def matprod(A=np.array,B=np.array):
    """renvoie le produit matricielle de A par B"""
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
    if np.linalg.det(A)==0 : 
        raise ValueError ("Determinant nul matrice non inversible")
    else :
        return np.linalg.inv(A)

def module(z):
    """renvoie le module de z"""
    return np.abs(z)

def integ(f,a=float,b=float,it=1000):
    """réalise l'intégrale de f sur a,b avec la méthode des rectangles, par défaut it (="itération")=1000"""
    S=0
    for i in range(it):
        S += (1/it)*(f(i*(a-b)/it))
    return S 
    
### Fonctions du problème

deg = 8 # degré maximal des polynomes étudidés
longueur_du_parcours = 20 # Longueur du parcours étudié

## à Optimiser

def Phi(P=np.array):
    """fonctions rayon moyen"""
    def A(x):
        F=fun_p(P)
        return module(F(x))
    
    return integ(A,0,1)

def dPhi(P=np.array, k=float):
    """renvoie la différentiel selon e_k de Phi de P"""
    p = fun_p(P) 
    if k> 2*deg+1 or k<0:
        raise ValueError("indice k trop grand ou trop petit")
    
    elif k%2 == 0 :
        ex = k//2

        def F(x):
            ep = 0
            if module(p(x)) == 0: 
                ep = 10**(-3)
            return (x**ex)*((p(x)).real)/(module(p(x+ep)))
        
    elif k%2 == 1 :
        ex = k//2 +1

        def F(x):
            ep = 0
            if module(p(x)) == 0: 
                ep = 10**(-3)
            return (x**ex)*((p(x)).imag)/(module(p(x+ep)))
        
    return integ(F,0,1)

def d2Phi(P=np.array,k=float,l=float):
    """renvoie la différentiel selon e_k et e_l de Phi de P"""

    p = fun_p(P) 
    if (k> 2*deg +1 or k<0) or (l> 2*deg +1 or l<0) :
        raise ValueError("indice k ou l trop grand ou trop petit")
    
    elif k%2==0 and l%2 ==1: 
        ex = k//2 + l//2
        def F(x):
            ep = 0
            if module(p(x)) == 0: 
                ep = 10**(-3)
            return -(x**(ex))*((p(x).imag)**2)/((module(p(x+ep)))**3)
        
    elif k%2 == 1 and l%2 == 1:
        ex = k//2 + l//2 + 2
        def F(x):
            ep = 0
            if module(p(x)) == 0: 
                ep = 10**(-3)
            return -(x**(ex))*((p(x).imag)*(p(x).real))/((module(p(x+ep)))**3)
        
    elif (k%2==0 and l%2==1) or (k%2==1 and l%2==0): 
        ex = k//2 + l//2 + 1
        def F(x):
            ep = 0
            if module(p(x)) == 0: 
                ep = 10**(-3)
            return -(x**(ex))*((p(x).real)**2)/((module(p(x+ep)))**3)     
        
    return integ(F,0,1)

## Contraintes

#Contrainte de longueur
def Psi(P=np.array):
    def A(x):
        F = fun_p(dp(P))
        return module(F(x))
    return integ(A,0,1)-longueur_du_parcours

def dPsi(P=np.array,k=float):
    A = fun_p(dp(P)) 

    if k > 2*deg+1 or k<0 :
        raise ValueError("indice k trop grand ou trop petit")

    elif k==0 or k==1:
        return 0
    
    elif k%2==0 :
        ex = k//2
        def F(x):
            ep = 0
            if module(A(x)) == 0: 
                ep = 10**(-3)
            return ex*(x**(ex-1))*((A(x)).real)/(module(A(x+ep)))
        
    elif k%2==1:
        ex = k//2 +1
        def F(x):
            ep = 0
            if module(A(x)) == 0: 
                ep = 10**(-3)
            return (ex)*(x**(ex-1))*((A(x)).imag)/(module(A(x+ep)))
        
    return integ(F,0,1)

def d2Psi(P=np.array,l=float,k=float):
    A = fun_p(dp(P))

    if k>2*deg+1 or k<0 or  l>2*deg+1 or l<0:
        raise ValueError("indice k ou l trop grand ou trop petit")
    
    elif k==0 or k==1 or l==0 or l==1:
        return 0
    
    elif k%2==0 or l%2==0: 
        ex = k//2 + l//2
        facteur = (k//2) * (l//2)
        def F(x):
            ep = 0
            if module(A(x)) == 0: 
                ep = 10**(-3)
            return -(facteur)*(x**(ex-2))*((A(x).imag)**2)/((module(A(x+ep)))**3)
        
    elif k%2==1 and l%2==1:
        ex = k//2 + l//2 +2
        facteur = ( 1 + k//2) * ( 1 + l//2 )
        def F(x):
            ep = 0
            if module(A(x)) == 0: 
                ep = 10**(-3)
            return -(facteur)*(x**(ex-2))*((A(x).real)**2)/((module(A(x+ep)))**3)
        
    elif (k%2==0 and l%2==1) or (k%2==1 and l%2==0): 
        ex = k//2 + l//2 +1
        facteur = ( k%2 + k//2) * ( l%2 + l//2 )

        def F(x):
            ep = 0
            if module(A(x)) == 0: 
                ep = 10**(-3)
            return -(facteur)*(x**(ex-2))*((A(x).imag)*(A(x).real))/((module(A(x+ep)))**3)
       
    return integ(F,0,1)

#Contrainte de cyclicité 

def Omega(P=np.array):
    p = fun_p(P)
    return (module(p(1)-p(0)))**2

def dOmega(P=np.array,k=float):
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

def d2Omega(P=np.array,l=float,k=float):
    """renvoie la différentiel selon e_k et e_l de Omega de P"""
    if k>2*deg+1 or k<0 or  l>2*deg+1 or l<0:
        raise ValueError("indice k ou l trop grand ou trop petit")
    
    if k == 0 or k==1 or l ==0 or l==1:
        return 0
    
    else:
        return 2
    
## Matrices et fonction Lagrangienne

def Lagrangien(P=np.array,lambda1=float,lambda2=float):
    Lag = np.zeros((1,2*deg+4))

    for i in range(2*deg+2):
        Lag[0][i] = dPhi(P,i) - lambda1*dPsi(P,i) - lambda2*dOmega(P,i)

    Lag[0][2*deg + 2] = Psi(P)
    Lag[0][2*deg + 3] = Omega(P)

    return Lag


def Jacob(P=np.array,lambda1=float,lambda2=float):

    Jacobienne = np.zeros((2*deg+4,2*deg+4))

    for i in range(2*deg+2):
        for j in range(2*deg+2):
            Jacobienne[i][j] =d2Phi(P,i,j) - lambda1*d2Psi(P,i,j) - lambda2*d2Omega(P,i,j)
    
    for i in range(2*deg+2):
        Jacobienne[i][2*deg+2] = dPsi(P,i)
        Jacobienne[2*deg+2][i] = dPsi(P,i)
        Jacobienne[i][2*deg+3] = dOmega(P,i)
        Jacobienne[2*deg+3][i] = dOmega(P,i)
    
    return Jacobienne 

### Résolution



##Algorithme de Raphson Newton

def Raphson_Newton(P0 = np.array, lambda1 =float, lambda2 = float , it = 50 ):
    X = np.zeros((1,2*deg+4))
    for i in range(deg+1):
        X[0][2*i]= (P0[i//2]).real
        X[0][2*i+1]= (P0[i//2]).imag
    X[0][2*deg+2] = lambda1
    X[0][2*deg+3] = lambda2

    for i in range(it):
        X =  X + scal(-1,matprod(matinv(Jacob(P0,lambda1,lambda2)),Lagrangien(P0,lambda1,lambda2)))

        for i in range(deg+1):
            P0[i] = complex(X[0][2*i] , X[0][2*i])
            lambda1 = X[0][2*deg+2]
            lambda2 = X[0][2*deg+3]

    return P0

#seed obtenue par interpolation  d'une solution discrète
#[complex(0.1,0.1),complex(35.1,140),complex(-1113,-2660),complex(10447,18904),complex(-45420,-67083),complex(105176,131527),complex(-133621,-144907),complex(87797,84052),complex(-23301,-19973)]
#fonctionne pour, deg =8