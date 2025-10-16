from collections.abc import Iterable
from mxlpy import Model, Derived
from mxlpy.surrogates import qss
from mxlpy.fns import mass_action_1s
import numpy as np

__all__ = ["psii", "add_photosystem"]


def psii(
    pq: float,
    pqh_2: float,
    pfd: float,
    k_pqh2: float,
    k_eq_qapq: float,
    k_h: float,
    k_f: float,
    k_p: float,
    psii_tot: float,
    sigma_psii: float,
    b_active_ratio: float,
    conversion_factor: float,
) -> Iterable[float]:
    """Returns amounts of PSII states calculated using steady state assumption"""
    # light = np.ones(len(P)) * light
    b0 = pfd*sigma_psii*conversion_factor + k_pqh2 * pqh_2 / k_eq_qapq
    b1 = k_h + k_f
    b2 = k_h + k_f + k_p
    M = np.array(
        [
            [-b0, b1, k_pqh2 * pq, 0],  # B0
            [pfd*sigma_psii*conversion_factor, -b2, 0, 0],  # B1
            [0, 0, pfd*sigma_psii*conversion_factor, -b1],  # B3
            [1, 1, 1, 1], #PSIItot should PSII_tot be changed to B_active?
        ]
    )
    A = np.array([0, 0, 0, psii_tot])
    return np.linalg.solve(M, A) * b_active_ratio


def psii_repair(b: float, psii_tot: float, k: float) -> float:
    """Mass action based rate equation describing PSII repair"""
    return k * (1 - (b / psii_tot))


def mass_action_sum(sum1: float, sum2: float, k: float) -> float:
    """First order rate equation where the concentration is the sum of two comnpunds"""
    return k * (sum1 + sum2)


def add_photosystem(model: Model) ->Model:
    """adding all component to necessary for and depending on PSII"""
    return(
        model
        .add_variable("B_active", 2.5)
        .add_surrogate(
            "Bs",
            surrogate=qss.Surrogate(
                model=psii,
                args=[
                    "P_ox",
                    "P_red",
                    "PFD",
                    "k_P_red",
                    "K_eq_QAPQ",
                    "k_H",
                    "k_F",
                    "k_P",
                    "PSII_tot",
                    "sigma_PSII",
                    "B_active_ratio",
                    "conversion_factor",
                ],
                outputs=["B0", "B1", "B2", "B3"],
            ),
        )

        # degradation of photosystem II
        .add_reaction(
            "v_deg_PSII",
            fn=mass_action_sum,
            args=["B1", "B3", "k_deg"],
            stoichiometry={"B_active": -1},
        )

        # repair of photosystem II
        .add_reaction(
            "v_rep_PSII",
            fn=psii_repair,
            args=["B_active", "PSII_tot", "k_rep"],
            stoichiometry={"B_active": 1},
        )

        #activity of PSII
        .add_reaction(
            "v_PSII",
            fn=mass_action_1s,
            args=["B1", "k_P"],
            stoichiometry={"H_lumen": Derived(fn = "two_div", args= ["b_H"]), "P_red": 1},
        )

    )
