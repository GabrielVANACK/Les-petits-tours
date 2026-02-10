import json
import os
import time
from Raphson_newton import cv_continue

STATE_FILE = "state.json"

def complex_to_tuple(z):
    return [z.real, z.imag]

def tuple_to_complex(t):
    return complex(t[0], t[1])

def save_state(P0, lambda1, lambda2, current_deg, deg_max):
    state = {
        "P0": [complex_to_tuple(z) for z in P0],
        "lambda1": lambda1,
        "lambda2": lambda2,
        "current_deg": current_deg,
        "deg_max": deg_max
    }
    with open(STATE_FILE, "w") as f:
        json.dump(state, f, indent=2)

def load_state():
    with open(STATE_FILE, "r") as f:
        state = json.load(f)
    P0 = [tuple_to_complex(z) for z in state["P0"]]
    return P0, state["lambda1"], state["lambda2"], state["current_deg"], state["deg_max"]

def initial_state():
    # seed issue de ton fichier
    P0 = [
        complex(0.1,0.1), complex(35.1,140), complex(-1113,-2660),
        complex(10447,18904), complex(-45420,-67083),
        complex(105176,131527), complex(-133621,-144907),
        complex(87797,84052), complex(-23301,-19973)
    ]
    return P0, 0.0, 0.0, len(P0)-1, 20   # exemple deg_max = 20

def main():
    if os.path.exists(STATE_FILE):
        print("Reprise d’un calcul existant.")
        P0, lambda1, lambda2, current_deg, deg_max = load_state()
    else:
        print("Initialisation d’un nouveau calcul.")
        P0, lambda1, lambda2, current_deg, deg_max = initial_state()
        save_state(P0, lambda1, lambda2, current_deg, deg_max)

    for d in range(current_deg , deg_max + 1):
        print(f"\n=== Calcul pour degré {d} ===")
        P0, lambda1, lambda2 = cv_continue(P0, lambda1, lambda2, d+1,d )
        save_state(P0, lambda1, lambda2, d, deg_max)
        print(f"Degré {d} sauvegardé.")

    print("Calcul terminé.")

if __name__ == "__main__":
    main()
