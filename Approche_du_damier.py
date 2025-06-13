#Essayons d'abord de trouver la solution sur un damier de nxn carrée (1x1 unités quelconques).
#On place le centre O au centre du damier (au centre de la case du centre si n impair et à l'intersection des cases centrales si n pair)
#On possède un "pion" qui peut se déplacer d'intersection en intersection (comme le roi aux échecs) et qui doit faire un tour de longueur m avec 0<m<4*n

#On représente le damier sous la forme d'une matrice contenant le sommet des carrées.

from math import *

##Fonctions usuels

def dist_O(n,M):
    """prend en entrée la taille du damier n un point M = (i,j) et renvoie la distance du point à l'origine""" 
    (i,j)=(M[0],M[1])
    return sqrt(((n/2)-i)**2+((n/2)-j)**2)

def longueur(M,P):
    """Prend en entrée deux point M =(i,j) et P=(p,q) et renvoie la distance qui les sépare"""
    (i,j)=(M[0],M[1])
    (p,q)=(P[0],P[1])
    return sqrt((i-q)**2+((j-p)**2))

def rectangle(f,a,b,n):
    """calcul d'intgrale par méthode des rectangle"""
    S=0
    for i in range(n+1):
        S += f(a + i*((b-a)/n))
    return S*((b-a)/n)

def moy_d_Or(n,M,P,precision=1000):
    """Prend en entrée la taille du damier n deux point M =(i,j) et M'=(p,q) et renvoie les distance moyenne à l'origine sur la ligne qui les sépare"""
    (i,j)=(M[0],M[1])
    (p,q)=(P[0],P[1])
    (x1,y1) = ((n/2)-i,(n/2)-j)
    (x2,y2) = ((n/2)-p,(n/2)-q)
    def f(x): return sqrt((x*x1+(1-x)*x2)**2 + (x*y1+(1-x)*y2)**2)
    return rectangle(f,0,1,precision)

#On note un parcours comme une liste de point par lequel on passe
# On dira qu'un parcours est cylcique si le point d'arrivée est également celui de départ.

def sans_repassage(l):
    """assure qu'un chemin ne repasse pas deux fois par le même point avec le même vecteur vitesse"""
    lngth = len(l)
    liste_doublons=[[] for i in l]
    for i in range(lngth):
        for j in range(lngth):
            if l[i]==l[j] and i!=j:
                (liste_doublons[i]).append(j)
    for a in range(lngth):
        if a!=0:
            v_a = (l[a][0]-l[a-1][0],l[a][1]-l[a-1][1]) 
            for k in liste_doublons[a]:
                if k!=0: 
                    v_k = (l[k][0]-l[k-1][0],l[k][1]-l[k-1][1]) 
                    if v_k == v_a:
                        return False
    #compléxité haute pour un algo simple comme ça je trouve

def est_cyclique(l):
    """prend un parcours en entrée renvoie true si il est cyclique et false sinon"""
    n = len(l)
    if l[0]== l[n-1]:
        return True
    else :
        return False
    
def longueur_parcours(l):  
    S=0
    m=len(l)
    for i in range(m-1):
        S+= longueur(l[i],l[i+1])
    return S

def moy_d_Or_parcours(l,n,precision=1000):
    """fait la distance moyenne à O du parcours"""
    S=0
    m = len(l)
    for i in range(m-1):
        S+=moy_d_Or(n,l[i],l[i+1],precision)
    return S/m

def min_parc(l,n):
    """prend une liste de parcours et n la taille du damier et renvoie le parcours avec la plus petite distance moyenne à l'origine"""
    print("start")
    Liste_moy_d_OR = [moy_d_Or_parcours(t,n,precision=100) for t in l]
    j=0
    n = len(l)
    for i in range(n):
        if Liste_moy_d_OR[j]>Liste_moy_d_OR[i]:
            j=i
            print(j)
    return l[j]

import matplotlib.pyplot as plt

def visualisation(p,n):
    """prend en entrée un parcours et la dimension du damier et renvoie le graph de ce parcours dans un damier de dimension n"""
    x_coords = [point[0] for point in p]
    y_coords = [point[1] for point in p]
    plt.figure(figsize=(8, 8))
    plt.plot(x_coords, y_coords, marker='o', label="Parcours")
    plt.plot(n/2,n/2, marker = 'o',label = "Origine")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xticks(range(n + 1))
    plt.yticks(range(n + 1))
    plt.xlim(-1, n+1)
    plt.ylim(-1, n+1)
    plt.title(f"Visualisation du Parcours (n={n})")
    plt.legend()
    plt.show()

def epsi(a,obj,e):
    """prend en entrée une valeur, un objectif et un réel epsilon et renvoie si la valeur est epsilon proche de la valeur objectif"""
    if obj-e<=a and a<=obj+e:
        return True
    else : return False

#L'objectif est maintenant de trouver un alogrithme qui génére des chemins d'une certaine longueur à +/- une erreur.
# Pour cela on commence par créer un lacet autour de l'origine que l'on agrandi jusqu'à ce qu'il soit de longueur l. Ce lacet est par nature cyclique

##Méthode du lacet

from random import randrange

def division(M,n):
    """prend un point et le divise en deux sous points"""
    (i,j)=(M[0],M[1])
    (d,c) = (randrange(-1,2,1),randrange(-1,2,1))
    (i1,j1,i2,j2)=(i+c,j+d,i-c,j-d,)
    for p in (i1,j2,i2,j2):
        if p>n:
            p=n
        elif p<0:
            p=0
    return [[i1,j1],[i2,j2]]
    
def homotéthie(M,n):
    O = (n/2,n/2)
    for i in range(2):
        if M[i]>O[i]:
            M[i]+=1
        elif M[i]<O[i]:
            M[i]-=1
    for p in M:
        if p>n:
            p=n
        elif p<0:
            p=0    

def parcours_lacet(l,e,n):
    """fait un parcour de longueur l (+/- e) par extension d'un lacet autour de l'origine dans un damier n*n"""
    parc = [(n/2,n/2)]
    max_opé = 100*l #pour que le programme ne tourne pas à l'infini
    opé=0
    while not epsi(longueur_parcours(parc),l,e) and opé<max_opé:    
        m=len(parc)
        if m==1: i=0
        else : i = randrange(0,m,1)
        coin = randrange(0,2,1)
        if coin==0:
            C = division(parc[i],n)
            if i==0: parc= C+parc[0:]+[C[0]]
            elif i==m-1: parc=  [C[1]] +parc[:m-2]+C
            else : parc = parc[:(i-1)] + C + parc[i:]
        else : 
            homotéthie(parc[i],n)
        opé+=1
    return parc

# Il est important de noté que pour l'instant cette algorithme comporte une grande portion d'aléatoire, il ne donc évidemment pas la solution optimale

##Méthodes naïves 

#A partir de maintenant et jusqu'à nouvel ordre on fera considère qu'on ne peut se déplacer que d'intersection adjacente en intersection adjacente.
#On se donne aussi une autre règle à partir de maintenant : interdiction de repasser par le même point avec le même vecteur vitesse (dans nos hypothèse cela revient à ne pas arriver à un point sans venir du même point deux fois)

#Dictionnaire contenant le meilleur chemin trouvé pour une clé [n,l] avec n la dimension et l la longueur du chemin
Opti_chemin = {}

#Initialisation du dictionnaire à partir d'un fichier texte.

f=open('Opti_dico.txt','r')
for ligne in f:
    c=ligne.split(';')
    simplify_char = "[()]"
    for i in simplify_char:
        c[2]= c[2].replace(i,"")
    coords = c[2].split(',')

    liste_points = []
    for i in range(len(coords)):
        if i%2==0:
            liste_points.append((int(coords[i]),int(coords[i+1])))
    
    Opti_chemin[(int(c[0]),float(c[1]))]=(liste_points,float(c[3]))

def save_dico():
    """enregistre le dictionnaire Opti_chemin mentioné dans le script, qui garde en mémoire les chemins les plus optimisé (d'un point de vue de la distance moyenne à l'origine) de longueur l dans un damier de taille n"""
    g=open('Opti_dico.txt','w')
    for x in Opti_chemin:
        g.write(f"{x[0]};{x[1]};{Opti_chemin[x][0]};{Opti_chemin[x][1]}"+"\n")



# Ce premier algorithme réalise un chemin aléatoirement, si il est cyclique et est plus optimal que celui du dictionnaire alors il l'enregistre dans celui ci 
def random_search(n,lngth):
    """enregistre dans un dictionnaire que tu peux acutalisé le meilleur parcours qu'il a trouvé pour une longueur n dans un damier de taille n, attention ici on ne s'autorise qu'à se déplacer vers des cases adjacentes"""
    #Ce première algorithme utilise encore de l'aléatoire ce qui ne garantit pas que la solution renvoyé soit optimal ou même qu'il renvoie une solution
    if (n,lngth) in Opti_chemin:
        opti = (Opti_chemin[(n,lngth)])[1]
    else :
        opti = 1000*n

    parc = [(randrange(0,n+1,1),randrange(0,n+1,1))]
    vitesses =[(0,0)]
    length = 1
    long=0
    while moy_d_Or_parcours(parc,n)<opti and long<lngth  and not epsi(long,lngth,2):
        vect = (randrange(-1,2,1),randrange(-1,2,1))
        M = (parc[length-1][0]+vect[0],parc[length-1][1]+vect[1])

        already_visited = False                     #faire un tableau de false pour savoir quels point j'ai déjà visité, mettre cette péripétie dans le cahier de TIPE
        for i in range(length):
            if M == parc[i] and vect == vitesses[i]:
                already_visited = True

        if not already_visited :
            parc.append(M)
            vitesses.append(vitesses)
            length+=1

        long=longueur_parcours(parc)

    moy = moy_d_Or_parcours(parc,n)
    if moy<opti and est_cyclique(parc): 
        Opti_chemin[(n,lngth)] = (parc,moy)
        return f"Nous avons trouvé un nouveau chemin {Opti_chemin[(n,lngth)]}"
    else : return "Pas de nouveau chemin trouvé"

#Ce deuxième algorithme test exhaustivement tout les chemins cyclique possibles et renvoie celui qui est le plus optimal

def brute_search_aux(n,parc,vitesse,lngth):
    m = len(parc)
    if epsi(lngth,0,sqrt(2)) :
        if est_cyclique(parc):
            return [parc]
        else :
            return ["fail"] 
            
    elif longueur(parc[0],parc[m-1])<=lngth :
        parc_liste = []
        vect_liste = [(1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1)]
        
        for v in vect_liste : 
            M = (parc[m-1][0] + v[0],parc[m-1][1] + v[1])

            already_visited = False
            for i in range(m):
                if M == parc[i] and v == vitesse[i]:
                    already_visited = True
            
            if not already_visited and (0<=M[0]<=n and 0<=M[1]<=n):
                rec_parc_liste =brute_search_aux(n,parc + [M],vitesse + [v],lngth-sqrt(v[0]**2+v[1]**2))
                parc_liste += [x for x in rec_parc_liste if x != "fail" ]
        return parc_liste
    else : 
        return ["fail"]

def brute_search(n,lngth):
    """enregistre dans un dictionnaire que tu peux acutalisé le meilleur parcours qu'il a trouvé pour une longueur n dans un damier de taille n, attention ici on ne s'autorise qu'à se déplacer vers des cases adjacentes """
    M_parc = []
    for x in range(int((n+1)/2)+1): #Ce problème est cconstant par symétrie centrale, trouver la solution optimal pour M (point de dépar) dans un quart de damier revient à la trouver quelque soit M
        for y in range(int((n+1)/2)+1):
            M =(x,y)  
            print(M) 
            result = brute_search_aux(n,[M],[(0,0)],lngth)
            if "fail" not in result:
                M_parc += result

    min = min_parc(M_parc,n)
    min_moy =moy_d_Or_parcours(min,n)
    if (n,lngth) in Opti_chemin:
        if min_moy<Opti_chemin[(n,lngth)][1]:
            print("succés")
            Opti_chemin[(n,lngth)]=(min,min_moy)
        else : print("Un chemin plus efficaces a déjà été trouvé")
    else : 
        print("succés")
        Opti_chemin[(n,lngth)]=(min,min_moy)

##DFS

#Par DFS, selon un article de Christine Solon on peut trouver les cycles d'un graphes à partir d'un DFS.
#WIP
#Marche pas pour l'instant
def cycle_BFS(n,source,l):
    """enregistre dans un dictionnaire que tu peux acutalisé le meilleur parcours qu'il a trouvé pour une longueur l dans un damier de taille n, attention ici on ne s'autorise qu'à se déplacer vers des cases adjacentes (diagonales exclus) """
    cycle =[]
    trajets = [(source,(0,0))]
    pile  =[trajets]
    while pile != []:
        t = pile.pop()
        lnght = len(t)
        M = t[lnght-1][0]
        V = [(M[0]+1,M[1]),(M[0]-1,M[1]),(M[0],M[1]-1),(M[0],M[1]+1)]
        if est_cyclique(t) and lnght==l:
            cycle.append(t)
        elif lnght < l :
            for v in V:
                if 0<=v[0]<=n and 0<=v[1]<=n:

                    new_traj = t + [(v,(v[0]-M[0],v[1]-M[1]))]
                    if new_traj not in trajets : 
                        trajets.append(new_traj)
                        pile.append(new_traj)
    return cycle

def opti_BFS(n,l):
    cycles = []
    for x in range(int((n+1)/2)+1): #Ce problème est cconstant par symétrie centrale, trouver la solution optimal pour M (point de dépar) dans un quart de damier revient à la trouver quelque soit M
        for y in range(int((n+1)/2)+1):
            cycles += cycle_BFS(n,(x,y),l)
    min = min_parc(cycles,n)
    min_moy =moy_d_Or_parcours(min,n)
    if (n,l) in Opti_chemin:
        if min_moy<Opti_chemin[(n,l)][1]:
            print("succés")
            Opti_chemin[(n,l)]=(min,min_moy)
        else : print("Un chemin plus efficaces a déjà été trouvé")
    else : 
        print("succés")
        Opti_chemin[(n,l)]=(min,min_moy)

#Après coup je ne suis pas sûr que cette technique nous fasse gagner grand chose face à Brute_search, à part de la simplicité d'écriture
#A retravailler     

##Récursion

#Pour trouver les chemins cyclique minimisant le rayon d'un point M vers M on s'intéresse à trouvers celui minimisant le rayon d'un point adjacent de M : P vers les autres points adjacents

#WIP, il marche vraiment pas du tout pour l'instant
def cycle_rec(n,l):
    #On initialise "naïvement" Min_mat
    
    def constructeur(i,j):
        if 0<abs(i%n-j%n) + abs(i//n - j//n)<=l:
            L = [(a,b) for a in range(min(i%n,j%n),max(i%n,j%n)+1) for b in range(min(i//n,j//n),max(i//n,j//n)+1)]#
            d_Or = moy_d_Or_parcours(L,n)
            return [L,d_Or,abs(i%n-j%n) + abs(i//n - j//n)]
        elif abs(i%n-j%n) + abs(i//n - j//n)==0:
            return [[(i%n,i//n)],pi*(n**2),1]
        else : return [[(-1,-1)],pi*(n**2),n*n] #valeurs absurde qui majorent ce qu'on veut trouver

    #Description Min_mat ->
    Min_mat = [[ constructeur(i,j) for i in range(n*n)] for j in range(n*n) ]

    #On applique Floyd Warshall
    for k in range(n*n):
        for i in range(n*n):
            for j in range(n*n):
                if epsi((Min_mat[i][k][2]+Min_mat[k][j][2]),l/(2**(n*n-k-1)),1): #je comprends pas ce qu'il ne va pas dans cette partie
                    moy_somme_parcours = (Min_mat[i][k][2]*Min_mat[i][k][1] + Min_mat[k][j][2]*Min_mat[k][j][1])/(Min_mat[k][j][2]+Min_mat[i][k][2]) #formule probablement fausse
                    if Min_mat[i][j][1]> moy_somme_parcours :
                        Min_mat[i][j] = [Min_mat[i][k][0]+Min_mat[k][j][0],moy_somme_parcours,Min_mat[i][k][2]+Min_mat[k][j][2]]

    #Recherche du minimum de min_Mat
    minimum = [[],pi*n*n,0]
    for i in range(n*n):
            print(Min_mat[i][i])
            if minimum[1]>Min_mat[i][i][1] and sans_repassage(Min_mat[i][i][0]):
                minimum=Min_mat[i][i]

    if (n,minimum[2]) in Opti_chemin:
        if minimum[1]<Opti_chemin[(n,minimum[2])][1]:
            print("succés")
            Opti_chemin[(n,minimum[2])]=(minimum[0],minimum[1])
        else : print("Un chemin plus efficaces a déjà été trouvé")
    else : 
        print("succés")
        Opti_chemin[(n,minimum[2])]=(minimum[0],minimum[1])

    return minimum

