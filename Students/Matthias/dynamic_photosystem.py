from mxlpy import Model, Derived
from mxlpy.fns import mass_action_1s
from cfunctions import( 
    diff,
    ratio,
    sum_psii_states,
    two_div,
)


__all__ = [
            "psii_repair",
           "reversible_mass_action_2_2",
           "add_dynamic_psii",
]


def psii_repair(b: float, psii_tot: float, k: float) -> float:
    """Mass action based rate equation describing PSII repair"""
    return k * (1 - (b / psii_tot))


def reversible_mass_action_2_2(
    s1: float,
    s2: float,
    p1: float,
    p2: float,
    k_fwd: float,
    k_bwd: float,
) -> float:
    """Reversible mass action equation for two substrates and two products."""
    return (k_fwd * s1 * s2) - (k_bwd * p1 * p2)


def add_dynamic_psii(model: Model) -> Model:
    return(
        model

        #PSII added as single variables for dynamic simulation
        .add_variables(
            {
                "B0": 2.5,
                "B1": 0,
                "B2": 0,
                "B3": 0,
            }
        )

        #PSII states
        .add_derived("B_active",fn=sum_psii_states,args=["B0","B1","B2","B3"])
        .add_derived("B_inactive",fn=diff,args=["PSII_tot","B_active"])

        #backward rate constant for plastoquinone reduction
        .add_derived("k_P_red_bw",fn=ratio,args=["k_P_red","K_eq_QAPQ"])

        #degradation of B1
        .add_reaction(
              "v_deg_B1",
              fn = mass_action_1s,
              args = ["B1","k_deg"],
              stoichiometry = {"B1": -1}
        )

        #degradation of B3
        .add_reaction(
              "v_deg_B3",
              fn = mass_action_1s,
              args = ["B3","k_deg"],
              stoichiometry = {"B3": -1}
        )

        #repair of photosystem II
        .add_reaction(
                "v_rep_PSII",
                fn = psii_repair,
                args = ["B_active","PSII_tot","k_rep"],
                stoichiometry = {"B0": 1}
        )

        #heat dissipation of B1
        .add_reaction(
                "v_kH_B1",
                fn=mass_action_1s,
                args=["B1","k_H"],
                stoichiometry={"B0": 1, "B1": -1}
        )

        #fluorescence emission from B1
        .add_reaction(
            "v_kF_B1",
            fn=mass_action_1s,
            args=["B1","k_F"],
            stoichiometry={"B0": 1, "B1": -1}
        )

        #light activation of B0
        .add_reaction(
            "v_kLII_B1",
            fn=mass_action_1s,
            args=["B0","k_LII"],
            stoichiometry={"B0": -1, "B1": 1}
        )

        #Water splitting (O2 oxidation) by B1
        .add_reaction(
            "v_kP",
            fn=mass_action_1s,
            args=["B1","k_P"],
            stoichiometry={"B1": -1, "B2": 1, "H_lumen": Derived(fn = two_div, args = ["b_H"])}
        )

        #plastochinon reduction by B2
        #are actually 2 plastochinons reduced in every reaction? with only one excited state transition
        #TODO: check if stoichiometry is correct
        .add_reaction(
            "v_kPQred",
            fn=reversible_mass_action_2_2,
            args=["B2","P_ox","B0","P_red","k_P_red","k_P_red_bw"],
            stoichiometry={"B0": 1, "B2": -1,"P_red": 1}
        )

        #fluorescence emission by B3
        .add_reaction(
            "v_kF_B3",
            fn=mass_action_1s,
            args=["B3","k_F"],
            stoichiometry={"B2": 1, "B3": -1}
        )

        #heat dissipation by B3
        .add_reaction(
            "v_kH_B3",
            fn=mass_action_1s,
            args=["B3","k_H"],
            stoichiometry={"B2": 1, "B3": -1}
        )

        #light activation of B2 to B3
        .add_reaction(
            "v_kLII_B3",
            fn=mass_action_1s,
            args=["B2","k_LII"],
            stoichiometry={"B2": -1, "B3": 1}
        )
    )