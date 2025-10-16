import mxlpy as mxl
from mxlpy.fns import neg_div
from cfunctions import ratio, four_div, neg_one_div

__all__ = [
    "PSII_repair",
    "add_rates_to_model",
    "convenience_kinetics",
    "custom_ATPsynthase",
    "mass_action_3",
    "mass_action_sum",
    "rate_cytb6f",
    "reversible_mass_action_2_2",
    "reversible_mass_action_3_3_or_max",
    "state_transition_rate",
]


def convenience_kinetics(
    v_max: float, f_ox: float, f_red: float, n_ox: float, n_red: float, K_eq: float
) -> float:
        """Rate equation described by convenience kinetics"""
        return( (v_max *
        ((f_red ** 2) * n_ox - ((f_ox ** 2) * n_red) / K_eq) /
        ((1 + f_red+f_red ** 2) * (1 + n_ox) + (1 + f_ox+f_ox**2) * (1 + n_red) - 1)) )
        

def state_transition_rate(
    p_ox: float,
    p_tot: float,
    k_eq_m_st: float,
    n_st: float,
    p_antennae: float,
    k_stt7: float
) -> float:
    """Rate equation describing antennae phosphorylation"""
    a = ((p_ox / p_tot) / k_eq_m_st) ** n_st
    return (k_stt7 * (1 / (1 + a)) * (1 - p_antennae))


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


# rate of cytochrome b6f
def rate_cytb6f(p_red: float,
                p_ox: float, 
                pc_red: float, 
                pc_ox: float,
                k_fw: float,
                k_bw: float,
                vmin: float
                ) -> float: 
    """
    Reversible mass action equation for three substrates and three products
    Returns either calculated rate or minimal rate given v_min_cyt_b6f, whichever is larger
    Specific rate equation for cytochrome b6f
    """
    vb6f = k_fw * (p_red * pc_ox**2) - k_bw * (p_ox * pc_red**2)
    return max(vb6f, vmin)


def mass_action_3(s1: float, s2: float, s3: float, k_fwd: float) -> float:
    """Mass action equation for three substrates."""
    return k_fwd * s1 * s2 * s3


def custom_ATPsynthase(ATPase_active, AP_spent, ATP, k_atpase, K_eq):
    """Mass action based rate equation describing ATP synthase activity"""
    return ATPase_active * k_atpase * (AP_spent - (ATP / K_eq))


def add_rates_to_model(model: mxl.Model) -> mxl.Model:
    """Adds all reaction rates to the given model and returns the updated model. Except for rates from Photosystems II"""
    return (
        model
        #reduction of NADP and Fd oxidation
        .add_reaction(
            "v_FNR",
            fn=convenience_kinetics,
            args=[
                "v_FNR_max",
                "f_ox_conv",
                "f_red_conv",
                "n_ox_conv",
                "n_red_conv",
                "K_eq_FNR",
            ],
            stoichiometry={"NADP_red": 1, "Fd_red": -2},
        )
        # synthesis of ATP from spent AP pool by ATPsynthase
        .add_reaction(
            "v_ATPsynthase",
            fn=custom_ATPsynthase,
            args=[
                "ATPase_active",
                "AP_spent",
                "ATP",
                "k_ATPsynthase",
                "K_eq_ATPsynthase",
            ],
            stoichiometry={
                "ATP": 1,
                "H_lumen": "HPR_buffered_stoich",
    
            },  
        )

        # activation of ATPase by light 
        .add_reaction(
            "v_act_ATPase",
            fn=mxl.fns.mass_action_2s,
            args=["heaviside_PFD", "ATPase_inactive", "k_actATPase"],
            stoichiometry={"ATPase_active": 1},
        )

        # deactivation of ATPase by light
        .add_reaction(
            "v_deactATPase",
            fn=mxl.fns.mass_action_2s,
            args=["heaviside_inverse", "ATPase_active", "k_deactATPase"],
            stoichiometry={"ATPase_active": -1},
        )

        # consumption of ATP to AP
        .add_reaction(
            "v_ATPconsumption",
            fn=mxl.fns.mass_action_1s,
            args=["ATP", "k_ATPconsumption"],
            stoichiometry={"ATP": -1},
        )

        # oxidization of NADPH
        .add_reaction(
            "v_NADPHconsumption",
            fn=mxl.fns.mass_action_1s,
            args=["NADP_red", "k_NADPHconsumption"],
            stoichiometry={"NADP_red": -1},
        )

        # cyclic Plastoquinon reduction by Ferredoxin
        .add_reaction(
            "v_FQR",
            fn=mass_action_3,
            args=["Fd_red", "Fd_red", "P_ox", "k_FQR"],
            stoichiometry={"Fd_red": -2, "P_red": 1},
        )

        # reduction of plastoquinon pool by outside system sources
        .add_reaction(
            "v_NDH",
            fn=mxl.fns.mass_action_1s,
            args=["P_ox", "k_NDH"],
            stoichiometry={"P_red": 1},
        )

        #leak through membrane to stroma
        .add_reaction(
            "v_leak",
            fn=mxl.fns.mass_action_1s,
            args=["H_diff_absolute", "k_leak"],
            stoichiometry={"H_lumen":  mxl.Derived(fn = neg_one_div, args=["b_H"])},
        )

        # oxidation of Polyquinon pool by external oxygen
        .add_reaction(
            "v_PTOX",
            fn=mxl.fns.mass_action_2s,
            args=["O2_ex", "P_red", "k_PTOX"],
            stoichiometry={"P_red": -1},
        )
        
        # electron transport through cytochrom c
        .add_reaction(
            "v_cytb6f",
            fn=rate_cytb6f,
            args=[
                "P_red",
                "P_ox",
                "PC_red",
                "PC_ox",
                "k_cytb6f",
                "k_cytb6f_bw",
                "v_cytb6f_min",
            ],
            stoichiometry={"P_red": -1, "PC_red": 2, "H_lumen": mxl.Derived(fn = four_div, args=["b_H"])},
        )

        # moving antennae from PSI to PSII and dephosphyrlating them
        .add_reaction(
            "v_Pph1",
            fn=mxl.fns.mass_action_1s,
            args=["L", "k_Pph1"],
            stoichiometry={"L": -1},
        )

        # phosphorylating antaenna and moving from PSII to PSI
        .add_reaction(
            "v_Stt7",
            fn=state_transition_rate,
            args=["P_ox", "P_tot", "K_M_ST", "n_ST", "L","k_Stt7"],
            stoichiometry={"L": 1},
        )


        #activity of PSI
        .add_reaction(
            "v_PSI",
            fn=mxl.fns.mass_action_1s,
            args=["Y_0", "k_LI"],
            stoichiometry={"PC_red": -2, "Fd_red": 2},
        )
    )
