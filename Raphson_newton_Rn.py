import numpy as np
import scipy.integrate as integrate

### Generalites

EPS_NORM = 1e-30
REG_MIN = 1e-9


def liste_shift(l, el):
    """Ajoute el en fin de liste et supprime le premier element."""

    n = len(l)
    for i in range(n - 1):
        l[i] = l[i + 1]
    l[n - 1] = el


def polprod(P=np.array, Q=np.array):
    """Produit de Cauchy de deux polynomes reels."""

    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    PQ = np.zeros(len(P) + len(Q) - 1)
    for i in range(len(P)):
        for j in range(len(Q)):
            PQ[i + j] += P[i] * Q[j]
    return PQ


def poladd(P=np.array, Q=np.array):
    """Addition de deux polynomes reels representes par leurs coefficients."""

    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    p = max(len(P), len(Q))
    out = np.zeros(p)
    out[: len(P)] += P
    out[: len(Q)] += Q
    return out


def polcomp(P, Q):
    """Retourne la composition P(Q(x)) pour deux polynomes reels."""

    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    result = np.zeros(1)
    current_power = np.array([1.0])

    for coeff in P:
        result = poladd(result, polprod([coeff], current_power))
        current_power = polprod(current_power, Q)

    return result


def fun_p(P=np.array):
    """Renvoie la fonction polynomiale associee a un polynome reel."""

    coeffs = np.asarray(P, dtype=float)
    powers = np.arange(len(coeffs))

    def poly_fun(x):
        return np.sum(coeffs * (x ** powers))

    return poly_fun


def dp(P=np.array):
    """Derivee d'un polynome reel."""

    coeffs = np.asarray(P, dtype=float)
    dP = np.zeros(len(coeffs))
    for k in range(1, len(coeffs)):
        dP[k - 1] = k * coeffs[k]
    return dP


def matprod(A=np.array, B=np.array):
    """Produit matriciel."""

    (n, m) = np.shape(A)
    (p, l) = np.shape(B)
    if p != m:
        raise ValueError("Matrices de tailles non compatibles")
    return np.matmul(A, B)


def matadd(A=np.array, B=np.array):
    """Addition matricielle."""

    if A.shape != B.shape:
        raise ValueError("les matrices ne sont pas de meme taille, impossible de les additionner")
    return A + B


def scal(mu=float, A=np.array):
    """Multiplication par un scalaire."""

    return np.multiply(A, mu)


def matinv(A=np.array):
    print(np.linalg.det(A))
    return np.linalg.inv(A)


def module(z):
    """Norme euclidienne."""

    if np.isscalar(z):
        return abs(z)
    return np.linalg.norm(np.asarray(z, dtype=float))


def safe_module(z):
    """Norme regularisee pour eviter les divisions par zero."""

    return max(module(z), EPS_NORM)


def normalize_curve(P):
    """Normalise l'entree en tableau reel de shape (m, deg + 1)."""

    arr = np.asarray(P)
    if np.iscomplexobj(arr):
        raise ValueError("P doit etre un tableau reel de shape (m, deg + 1).")
    arr = arr.astype(float)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.ndim != 2:
        raise ValueError("P doit etre un tableau numpy a deux dimensions.")
    return arr


# Seed issue de l'ancienne representation complexe, reecrite comme courbe de R^2.
seed = np.array(
    [
        [0.1, 35.1, -1113, 10447, -45420, 105176, -133621, 87797, -23301],
        [0.1, 140.0, -2660, 18904, -67083, 131527, -144907, 84052, -19973],
    ],
    dtype=float,
)


def regularization_for_degree(deg):
    """Calendrier de regularisation utilise pendant la continuation."""

    return 5 / (2 ** (deg // 10))


def seed_amplitude_for_degree(deg):
    """Amplitude de perturbation pour initialiser le degre suivant."""

    return regularization_for_degree(deg) / 1000


def default_seed_weights(m):
    """Poids deterministes pour perturber les composantes d'une seed."""

    weights = np.linspace(1.0, 0.5, m)
    weights[1::2] *= -1
    return weights


def curve_key(P):
    """Cle de memoisation pour une courbe polynomiale."""

    arr = normalize_curve(P)
    return tuple(tuple(float(value) for value in row) for row in arr)


def fun_curve(P=np.array):
    """Renvoie la fonction x -> P(x) dans R^m."""

    coeffs = normalize_curve(P)
    powers = np.arange(coeffs.shape[1])

    def poly_fun(x):
        return np.sum(coeffs * (x ** powers), axis=1)

    return poly_fun


def dp_curve(P=np.array):
    """Derivee composante par composante d'une courbe polynomiale."""

    coeffs = normalize_curve(P)
    m, n = coeffs.shape
    dP = np.zeros((m, n))
    for k in range(1, n):
        dP[:, k - 1] = k * coeffs[:, k]
    return dP


def coeff_count(deg):
    """Nombre de coefficients optimises par composante."""

    if deg < 1:
        raise ValueError("deg doit etre superieur ou egal a 1.")
    return deg - 1


def variable_count(m, deg):
    """Nombre total de coefficients optimises."""

    return m * coeff_count(deg)


def system_size(m, deg):
    """Taille du systeme de Newton avec les multiplicateurs."""

    return variable_count(m, deg) + 2


def decode_variable_index(k, deg, m):
    """Associe un indice de variable au couple (composante, degre)."""

    var_count = variable_count(m, deg)
    if k < 0 or k >= var_count:
        raise ValueError("indice k trop grand ou trop petit")
    per_component = coeff_count(deg)
    component = k // per_component
    degree = 2 + (k % per_component)
    return component, degree


def vector_to_curve(X, deg, m):
    """Reconstruit une courbe depuis les coefficients optimises."""

    coeffs = np.zeros((m, deg + 1))
    if deg >= 2:
        coeffs[:, 2:] = np.asarray(X, dtype=float).reshape(m, deg - 1)
    return coeffs


def curve_to_vector(P, deg):
    """Vectorise les coefficients de degre >= 2 composante par composante."""

    coeffs = normalize_curve(P)
    if coeffs.shape[1] != deg + 1:
        raise ValueError("Le degre ne correspond pas a la taille de P.")
    if deg < 2:
        return np.zeros((0, 1))
    return coeffs[:, 2:].reshape(-1, 1)


def make_seed(P=None, target_deg=None, amplitude=0.0, component_weights=None):
    """Construit une seed reelle de shape (m, deg + 1).

    L'entree peut etre :
    - un tableau reel de shape (m, deg + 1) ;
    - un vecteur complexe 1D, interprete comme une courbe de R^2 ;
    - None, auquel cas on repart de la seed par defaut du module.
    """

    if P is None:
        coeffs = seed.copy()
    else:
        raw = np.asarray(P)
        if np.iscomplexobj(raw):
            if raw.ndim != 1:
                raise ValueError("Une seed complexe doit etre un vecteur 1D.")
            coeffs = np.vstack((raw.real, raw.imag)).astype(float)
        else:
            coeffs = normalize_curve(raw).copy()

    current_deg = coeffs.shape[1] - 1
    if target_deg is None:
        target_deg = current_deg
    if target_deg < 1:
        raise ValueError("target_deg doit etre superieur ou egal a 1.")

    if target_deg > current_deg:
        coeffs = np.pad(coeffs, ((0, 0), (0, target_deg - current_deg)))
    elif target_deg < current_deg:
        coeffs = coeffs[:, : target_deg + 1].copy()

    coeffs[:, : min(2, target_deg + 1)] = 0.0

    if amplitude != 0.0 and target_deg >= 2:
        tcheb = polcomp(give_pol_T(target_deg), [-1.0, 2.0])
        if len(tcheb) < target_deg + 1:
            tcheb = np.pad(tcheb, (0, target_deg + 1 - len(tcheb)))
        if component_weights is None:
            component_weights = default_seed_weights(coeffs.shape[0])
        component_weights = np.asarray(component_weights, dtype=float)
        if component_weights.shape != (coeffs.shape[0],):
            raise ValueError("component_weights doit etre un vecteur de taille m.")
        for component, weight in enumerate(component_weights):
            coeffs[component] = poladd(coeffs[component], amplitude * weight * tcheb)
        coeffs[:, :2] = 0.0

    return coeffs


integ_mem = {}


def integ(f, a=float, b=float):
    """Realise l'integrale de f sur [a, b] avec scipy.integrate.quad."""

    if (f, a, b) in integ_mem:
        return integ_mem[f, a, b]
    integ_mem[f, a, b] = integrate.quad(lambda x: f(x), a, b)[0]
    return integ_mem[f, a, b]


### Fonctions du probleme

longueur_du_parcours = 20
reg = 3

Pol_fun_mem = {}


def get_curve_data(P):
    """Memoise P, P' et P'' comme fonctions vectorielles."""

    key = curve_key(P)
    if key not in Pol_fun_mem:
        coeffs = normalize_curve(P)
        dP = dp_curve(coeffs)
        d2P = dp_curve(dP)
        Pol_fun_mem[key] = (fun_curve(coeffs), fun_curve(dP), fun_curve(d2P))
    return Pol_fun_mem[key]


def Phi(P=np.array, deg=int):
    """Fonction rayon moyen generalisee a R^m."""

    coeffs = normalize_curve(P)

    def A(x):
        F, _, d2F = get_curve_data(coeffs)
        return module(F(x)) + reg * (module(d2F(x)) ** 2)

    return integ(A, 0, 1)


def dPhi(P=np.array, k=float, deg=int):
    """Differentielle de Phi selon un coefficient de P."""

    coeffs = normalize_curve(P)
    m = coeffs.shape[0]
    component, exponent = decode_variable_index(int(k), deg, m)
    p, _, d2p = get_curve_data(coeffs)

    def F(x):
        value = p(x)
        d2value = d2p(x)
        second_basis = exponent * (exponent - 1) * (x ** (exponent - 2))
        return (x ** exponent) * value[component] / safe_module(value) + reg * second_basis * d2value[component]

    return integ(F, 0, 1)


def d2Phi(P=np.array, k=float, l=float, deg=int):
    """Differentielle seconde de Phi."""

    coeffs = normalize_curve(P)
    m = coeffs.shape[0]
    component_k, exponent_k = decode_variable_index(int(k), deg, m)
    component_l, exponent_l = decode_variable_index(int(l), deg, m)
    p, _, _ = get_curve_data(coeffs)

    def F(x):
        value = p(x)
        norm_value = safe_module(value)
        total_power = x ** (exponent_k + exponent_l)
        if component_k == component_l:
            complement = (norm_value ** 2) - (value[component_k] ** 2)
            basis_k = exponent_k * (exponent_k - 1) * (x ** (exponent_k - 2))
            basis_l = exponent_l * (exponent_l - 1) * (x ** (exponent_l - 2))
            return -(total_power * complement) / (norm_value ** 3) + reg * basis_k * basis_l
        return -(total_power * value[component_k] * value[component_l]) / (norm_value ** 3)

    return integ(F, 0, 1)


### Contraintes


def Psi(P=np.array, deg=int):
    """Contrainte de longueur dans R^m."""

    coeffs = normalize_curve(P)

    def A(x):
        _, dF, _ = get_curve_data(coeffs)
        return module(dF(x))

    return integ(A, 0, 1) - longueur_du_parcours


def dPsi(P=np.array, k=float, deg=int):
    """Differentielle de la contrainte de longueur."""

    coeffs = normalize_curve(P)
    m = coeffs.shape[0]
    component, exponent = decode_variable_index(int(k), deg, m)
    _, A, _ = get_curve_data(coeffs)

    def F(x):
        value = A(x)
        return exponent * (x ** (exponent - 1)) * value[component] / safe_module(value)

    return integ(F, 0, 1)


def d2Psi(P=np.array, l=float, k=float, deg=int):
    """Differentielle seconde de la contrainte de longueur."""

    coeffs = normalize_curve(P)
    m = coeffs.shape[0]
    component_k, exponent_k = decode_variable_index(int(k), deg, m)
    component_l, exponent_l = decode_variable_index(int(l), deg, m)
    _, A, _ = get_curve_data(coeffs)

    def F(x):
        value = A(x)
        norm_value = safe_module(value)
        factor = exponent_k * exponent_l
        power = x ** (exponent_k + exponent_l - 2)
        if component_k == component_l:
            complement = (norm_value ** 2) - (value[component_k] ** 2)
            return -(factor * power * complement) / (norm_value ** 3)
        return -(factor * power * value[component_k] * value[component_l]) / (norm_value ** 3)

    return integ(F, 0, 1)


def Omega(P=np.array, deg=int):
    """Contrainte de cyclicite dans R^m."""

    coeffs = normalize_curve(P)
    p, _, _ = get_curve_data(coeffs)
    return module(p(1) - p(0)) ** 2


def dOmega(P=np.array, k=float, deg=int):
    """Differentielle de Omega selon un coefficient de P."""

    coeffs = normalize_curve(P)
    m = coeffs.shape[0]
    component, _ = decode_variable_index(int(k), deg, m)
    p, _, _ = get_curve_data(coeffs)
    return 2 * (p(1) - p(0))[component]


def d2Omega(P=np.array, l=float, k=float, deg=int):
    """Differentielle seconde de Omega."""

    coeffs = normalize_curve(P)
    m = coeffs.shape[0]
    component_k, _ = decode_variable_index(int(k), deg, m)
    component_l, _ = decode_variable_index(int(l), deg, m)
    if component_k == component_l:
        return 2
    return 0


### Polynomes de Tchebychev

Tcheb_pol = {}
dico_pass_mat = {}


def give_pol_T(n=int):
    """Renvoie le n-ieme polynome de Tchebychev."""

    if n == 0:
        return np.array([1.0])
    if n == 1:
        return np.array([0.0, 1.0])
    if n in Tcheb_pol:
        return Tcheb_pol[n]
    Tcheb_pol[n] = poladd(polprod([0.0, 2.0], give_pol_T(n - 1)), polprod([-1.0], give_pol_T(n - 2)))
    return Tcheb_pol[n]


def PassMmat_scalar(n=int):
    """Matrice de passage Tchebychev -> monomes pour un polynome reel."""

    if ("Tcheb to Mono scalar", n) in dico_pass_mat:
        return dico_pass_mat[("Tcheb to Mono scalar", n)]
    P = np.zeros((n + 1, n + 1))
    for j in range(n + 1):
        Tj = polcomp(give_pol_T(j), [-1.0, 2.0])
        P[: len(Tj), j] = Tj
    dico_pass_mat[("Tcheb to Mono scalar", n)] = P
    return P


def PassTmat_scalar(n=int):
    """Matrice de passage monomes -> Tchebychev pour un polynome reel."""

    if ("Mono to Tcheb scalar", n) in dico_pass_mat:
        return dico_pass_mat[("Mono to Tcheb scalar", n)]
    P = np.linalg.inv(PassMmat_scalar(n))
    dico_pass_mat[("Mono to Tcheb scalar", n)] = P
    return P


def PassMmat(n=int, m=int):
    """Matrice de passage bloc diagonale Tchebychev -> monomes dans R^m."""

    if ("Tcheb to Mono", n, m) in dico_pass_mat:
        return dico_pass_mat[("Tcheb to Mono", n, m)]
    P = np.kron(np.eye(m), PassMmat_scalar(n))
    dico_pass_mat[("Tcheb to Mono", n, m)] = P
    return P


def PassTmat(n=int, m=int):
    """Matrice de passage bloc diagonale monomes -> Tchebychev dans R^m."""

    if ("Mono to Tcheb", n, m) in dico_pass_mat:
        return dico_pass_mat[("Mono to Tcheb", n, m)]
    P = np.kron(np.eye(m), PassTmat_scalar(n))
    dico_pass_mat[("Mono to Tcheb", n, m)] = P
    return P


def PassMmat_reduced(n=int, m=int):
    """Passage Tchebychev -> monomes sur les degres 2..n pour chaque composante."""

    if ("Tcheb to Mono reduced", n, m) in dico_pass_mat:
        return dico_pass_mat[("Tcheb to Mono reduced", n, m)]
    reduced = np.kron(np.eye(m), PassMmat_scalar(n)[2:, 2:])
    dico_pass_mat[("Tcheb to Mono reduced", n, m)] = reduced
    return reduced


def PassTmat_reduced(n=int, m=int):
    """Passage monomes -> Tchebychev sur les degres 2..n pour chaque composante."""

    if ("Mono to Tcheb reduced", n, m) in dico_pass_mat:
        return dico_pass_mat[("Mono to Tcheb reduced", n, m)]
    reduced = np.kron(np.eye(m), PassTmat_scalar(n)[2:, 2:])
    dico_pass_mat[("Mono to Tcheb reduced", n, m)] = reduced
    return reduced


### Matrices et fonction lagrangienne


def Lagrangien(P=np.array, lambda1=float, lambda2=float, deg=int):
    coeffs = normalize_curve(P)
    m = coeffs.shape[0]
    var_count = variable_count(m, deg)
    Lag = np.zeros((var_count + 2, 1))

    for i in range(var_count):
        Lag[i][0] = dPhi(coeffs, i, deg) - lambda1 * dPsi(coeffs, i, deg) - lambda2 * dOmega(coeffs, i, deg)

    Lag[var_count][0] = Psi(coeffs, deg)
    Lag[var_count + 1][0] = Omega(coeffs, deg)
    return Lag


def Jacob(P=np.array, lambda1=float, lambda2=float, deg=int):
    coeffs = normalize_curve(P)
    m = coeffs.shape[0]
    var_count = variable_count(m, deg)
    Jacobienne = np.zeros((var_count + 2, var_count + 2))

    for i in range(var_count):
        for j in range(var_count):
            Jacobienne[i][j] = d2Phi(coeffs, i, j, deg) - lambda1 * d2Psi(coeffs, j, i, deg) - lambda2 * d2Omega(
                coeffs, j, i, deg
            )

    for i in range(var_count):
        Jacobienne[i][var_count] = dPsi(coeffs, i, deg)
        Jacobienne[var_count][i] = dPsi(coeffs, i, deg)
        Jacobienne[i][var_count + 1] = dOmega(coeffs, i, deg)
        Jacobienne[var_count + 1][i] = dOmega(coeffs, i, deg)

    return Jacobienne


### Resolution

mu = 10 ** (-3)


def lagrangian_norm(P, lambda1, lambda2, deg):
    """Norme euclidienne du lagrangien."""

    return float(np.linalg.norm(Lagrangien(P, lambda1, lambda2, deg), 2))


def Raphson_Newton(P0=np.array, lambda1=float, lambda2=float, it=50, deg=int):
    coeffs = make_seed(P0, target_deg=deg)
    m = coeffs.shape[0]
    var_count = variable_count(m, deg)

    X = np.zeros((var_count + 2, 1))
    X[:var_count] = curve_to_vector(coeffs, deg)
    X[var_count][0] = lambda1
    X[var_count + 1][0] = lambda2

    PassTM = PassMmat_reduced(deg, m)
    PassMT = PassTmat_reduced(deg, m)
    if var_count > 0:
        X[:var_count] = matprod(PassMT, X[:var_count])

    global mu

    for iteration in range(it):
        J = Jacob(coeffs, lambda1, lambda2, deg)
        if var_count > 0:
            J[:, :var_count] = matprod(J[:, :var_count], PassTM)
        L = scal(-1, Lagrangien(coeffs, lambda1, lambda2, deg))

        try:
            deltaX = np.linalg.solve(J, L)
        except np.linalg.LinAlgError as exc:
            raise np.linalg.LinAlgError(
                f"Jacobienne singuliere pendant Raphson_Newton (deg={deg}, iteration={iteration})."
            ) from exc

        X = X + deltaX

        if var_count > 0:
            a = matprod(PassTM, X[:var_count])
            coeffs = vector_to_curve(a[:, 0], deg, m)
        else:
            coeffs = np.zeros((m, deg + 1))
        coeffs[:, : min(2, deg + 1)] = 0.0
        lambda1 = X[var_count][0]
        lambda2 = X[var_count + 1][0]

    return coeffs, lambda1, lambda2


def cv_reg(
    P0=np.array,
    lambda1=float,
    lambda2=float,
    deg=int,
    newton_steps=35,
    inner_tol=10 ** (-5),
    outer_tol=10 ** (-5),
    max_inner_loops=12,
    max_outer_loops=12,
):
    """Boucle de continuation sur le parametre de regularisation."""

    coeffs = make_seed(P0, target_deg=deg)
    global reg
    best_coeffs = coeffs.copy()
    best_lambda1 = lambda1
    best_lambda2 = lambda2
    best_norm = float("inf")
    outer_history = []

    for outer_index in range(max_outer_loops):
        global mu
        mu = 10 ** (-3)
        previous_norm = None
        total_newton_iterations = 0

        for _ in range(max_inner_loops):
            coeffs, lambda1, lambda2 = Raphson_Newton(coeffs, lambda1, lambda2, newton_steps, deg)
            total_newton_iterations += newton_steps
            current_norm = lagrangian_norm(coeffs, lambda1, lambda2, deg)

            if not np.isfinite(current_norm):
                raise FloatingPointError("La norme du lagrangien n'est pas finie pendant cv_reg.")

            if current_norm < best_norm:
                best_norm = current_norm
                best_coeffs = coeffs.copy()
                best_lambda1 = lambda1
                best_lambda2 = lambda2

            print(
                f"Avec reg={reg:.3e}, {total_newton_iterations} iterations de Newton : ||Lag|| = {current_norm:.6e}"
            )

            if previous_norm is not None and abs(current_norm - previous_norm) <= inner_tol:
                break
            previous_norm = current_norm

        outer_history.append(current_norm)
        print(f"Fin de boucle externe {outer_index} : ||Lag|| = {current_norm:.6e}\n")

        if reg <= REG_MIN:
            break
        if len(outer_history) >= 2 and abs(outer_history[-1] - outer_history[-2]) <= outer_tol:
            break

        delta = 0.0 if len(outer_history) < 2 else abs(outer_history[-1] - outer_history[-2])
        next_reg = reg / (10 + 1 / (EPS_NORM + delta))
        reg = max(next_reg, REG_MIN)

    return best_coeffs, best_lambda1, best_lambda2


def cv_continue(P0=np.array, lambda1=float, lambda2=float, deg_max=int, deg=int):
    """Continuation en degre pour une courbe de R^m."""

    coeffs = make_seed(P0, target_deg=deg)
    global reg

    for current_deg in range(deg, deg_max + 1):
        print(f"\n \n Degre des polynomes etudies : {current_deg}\n")
        reg = regularization_for_degree(current_deg)
        coeffs, lambda1, lambda2 = cv_reg(coeffs, lambda1, lambda2, current_deg)
        if current_deg < deg_max:
            coeffs = make_seed(
                coeffs,
                target_deg=current_deg + 1,
                amplitude=seed_amplitude_for_degree(current_deg),
            )

    return coeffs, lambda1, lambda2


def satis_facteur(P):
    coeffs = normalize_curve(P)
    print(f"Cyclicite : \n distance entre P(1) et P(0) est {np.sqrt(Omega(coeffs))}")
    print(
        f"Longueur : \n La longueur voulue est {longueur_du_parcours}, le parcours fait {Psi(coeffs) + longueur_du_parcours} soit une difference de {Psi(coeffs)}"
    )
    print(f"Rayon moyen : \n le rayon moyen a l'origine de P sur [0,1] est {Phi(coeffs)}")


def T(X, deg=int, m=None):
    """Fonction auxiliaire de type moindres carres."""

    values = np.asarray(X, dtype=float)
    if values.ndim == 1:
        values = values.reshape(-1, 1)

    if m is None:
        if deg < 2:
            raise ValueError("m doit etre fourni si deg < 2.")
        var_count = values.shape[0] - 2
        per_component = coeff_count(deg)
        if per_component <= 0 or var_count % per_component != 0:
            raise ValueError("Impossible d'inferer m depuis X et deg.")
        m = var_count // per_component

    var_count = variable_count(m, deg)
    lambda1 = values[var_count][0]
    lambda2 = values[var_count + 1][0]

    P = np.zeros((m, deg + 1))
    if var_count > 0:
        a = matprod(PassMmat_reduced(deg, m), values[:var_count])
        P = vector_to_curve(a[:, 0], deg, m)

    L = Lagrangien(P, lambda1, lambda2, deg)
    J = Jacob(P, lambda1, lambda2, deg)
    if var_count > 0:
        J[:, :var_count] = matprod(J[:, :var_count], PassMmat_reduced(deg, m))

    M = matprod(np.transpose(J), L)
    Mt = np.transpose(M)
    return scal(1 / 2, matprod(Mt, M))
