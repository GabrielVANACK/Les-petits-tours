import json
import os

import numpy as np

from Raphson_newton_Rn import cv_reg
from Raphson_newton_Rn import make_seed
from Raphson_newton_Rn import regularization_for_degree
from Raphson_newton_Rn import seed as default_seed
from Raphson_newton_Rn import seed_amplitude_for_degree

STATE_FILE = "state_Rn.json"


def array_to_list(P):
    return np.asarray(P, dtype=float).tolist()


def list_to_array(values):
    return np.asarray(values, dtype=float)


def save_state(
    P0,
    lambda1,
    lambda2,
    current_deg,
    deg_max,
    last_completed_deg=None,
    last_solution=None,
    state_file=STATE_FILE,
):
    state = {
        "P0": array_to_list(P0),
        "lambda1": float(lambda1),
        "lambda2": float(lambda2),
        "current_deg": int(current_deg),
        "deg_max": int(deg_max),
        "last_completed_deg": None if last_completed_deg is None else int(last_completed_deg),
        "last_solution": None if last_solution is None else array_to_list(last_solution),
    }
    with open(state_file, "w", encoding="utf-8") as f:
        json.dump(state, f, indent=2)


def load_state(state_file=STATE_FILE):
    with open(state_file, "r", encoding="utf-8") as f:
        state = json.load(f)
    P0 = list_to_array(state["P0"])
    last_completed_deg = state.get("last_completed_deg")
    last_solution = state.get("last_solution")
    if last_solution is not None:
        last_solution = list_to_array(last_solution)
    return P0, state["lambda1"], state["lambda2"], state["current_deg"], state["deg_max"], last_completed_deg, last_solution


def initial_state(deg_max=20):
    P0 = make_seed(default_seed)
    current_deg = P0.shape[1] - 1
    return P0, 0.0, 0.0, current_deg, deg_max, current_deg - 1


def main():
    if os.path.exists(STATE_FILE):
        print("Reprise d'un calcul existant en dimension m.")
        P0, lambda1, lambda2, current_deg, deg_max, last_completed_deg, _ = load_state()
    else:
        print("Initialisation d'un nouveau calcul en dimension m.")
        P0, lambda1, lambda2, current_deg, deg_max, last_completed_deg = initial_state()
        save_state(P0, lambda1, lambda2, current_deg, deg_max, last_completed_deg)

    if current_deg > deg_max:
        print("Calcul deja termine.")
        return

    while current_deg <= deg_max:
        print(f"\n=== Calcul pour degre {current_deg} ===")
        print(f"Regularisation initiale : {regularization_for_degree(current_deg):.3e}")

        P0 = make_seed(P0, target_deg=current_deg)
        P0, lambda1, lambda2 = cv_reg(P0, lambda1, lambda2, current_deg)
        last_solution = P0.copy()
        last_completed_deg = current_deg

        if current_deg < deg_max:
            P0 = make_seed(
                P0,
                target_deg=current_deg + 1,
                amplitude=seed_amplitude_for_degree(current_deg),
            )

        current_deg += 1
        save_state(P0, lambda1, lambda2, current_deg, deg_max, last_completed_deg, last_solution=last_solution)
        print(f"Degre {last_completed_deg} sauvegarde dans {STATE_FILE}.")

    print("Calcul termine.")


if __name__ == "__main__":
    main()
