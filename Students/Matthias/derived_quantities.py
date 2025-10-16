from mxlpy import Model
from mxlpy.fns import neg

from cfunctions import (
    calc_ph,
    diff,
    heaviside,
    inverse,
    k_eq_atpsynth,
    k_eq_cytb6f,
    k_eq_fnr,
    k_eq_p700_c,
    k_eq_p700_fd,
    k_eq_qapq,
    k_li,
    k_lii,
    ratio,
    sigma_psi,
    sigma_psii,
    y_0,
    fluorescence,
    hpr_make_stoich,
    hpr_dy_stoich,
    psii_dy_stoich,
    b6f_dy_stoich,
    leak_dy_stoich,
    ph_diff,
    sum,
    charge_separations2,
    charge_separations
)
from delta_psi_variable import psi_from_ph

__all__ = ["add_derived_quantities"]


def add_derived_quantities(model: Model) -> Model:
    """Adds all derived quantities to the model necessary for base model without Photosystem description, which is seperately added"""
    return (
        model

        #two different versions of charge separations, one for PSI per PSI total and one for PSI+PSII per PSI total
        .add_derived(
            "charge_separations2",
            fn = charge_separations2,
            args = ["v_PSI", "v_kP", "PSI_tot"]
        )
        .add_derived(
            "charge_separations",
            fn = charge_separations,
            args = ["v_PSI", "PSI_tot"]
        )

        #excitation rates of PSII
        .add_derived(
            "excitations",
            fn=sum,
            args=["v_kLII_B1","v_kLII_B3"]
        )

        #quantum yield of PSII calculated by photchemistry rate per excitations
        .add_derived(
            "quantum_yield",
            fn=ratio,
            args=["v_kPQred","excitations"]
        )

        #fluorescence from PSII relaxations
        .add_derived(
            "PSII_fluorescence",
            fn=fluorescence,
            args=["k_F","k_H","k_P","B0","B2","sigma_PSII"]
        )
        
        # equilibrium constants calculated from Gibbs energy directly or redox potential
        .add_derived(
            "K_eq_ATPsynthase",
            fn=k_eq_atpsynth,
            args=["pH_lumen", "G0_ATP", "pH_stroma", "R", "T", "Pi_mol","HPR"],
        )

        #equilibrium constants for reactions
        #NOTE: k_fw/K_eq equals k_bw in this context
        .add_derived(
            "K_eq_FNR",
            fn=k_eq_fnr,
            args=["F", "R", "T", "E0_Fd_redox", "E0_NADP_redox", "pH_stroma"],
        )
        .add_derived(
            "K_eq_cytb6f",
            fn=k_eq_cytb6f,
            args=["pH_lumen", "F", "E0_PQ_redox", "R", "T", "E0_PC_redox", "pH_stroma"],
        )
        .add_derived(
            "K_eq_QAPQ",
            fn=k_eq_qapq,
            args=["F", "E0_QA_redox", "E0_PQ_redox", "pH_stroma", "R", "T"],
        )
        .add_derived("K_eq_p700_C",
                     fn=k_eq_p700_c,
                     args=["F", "R", "T", "E0_PC_redox","E0_P700_redox"]
        )
        .add_derived(
            "K_eq_p700_Fd",
            fn=k_eq_p700_fd,
            args=["F", "R", "T", "E0_Fd_redox","E0_FA_redox"]
        )

        #Derived values from heaviside step function used to determine wether light is on or off for ATPsynthase regulation
        #can be removed when different ATPsynthase regulation is used
        .add_derived("heaviside_PFD", fn=heaviside, args=["PFD"])
        .add_derived("heaviside_inverse", fn=inverse, args=["heaviside_PFD"])
        .add_derived("HPR_stoich", fn=neg, args=["HPR"])

        # important for parameter dependent stoichiometry, do not change over time
        .add_derived("HPR_buffered_stoich", fn=hpr_make_stoich, args=["HPR","b_H"])
        .add_derived("HPR_dy_stoich", fn=hpr_dy_stoich, args=["HPR", "F", "Capac"])
        .add_derived("PSII_dy_stoich", fn=psii_dy_stoich, args=["F", "Capac"])
        .add_derived("b6f_dy_stoich", fn=b6f_dy_stoich, args=["F", "Capac"])
        .add_derived("leak_dy_stoich", fn=leak_dy_stoich, args=["F", "Capac"])

        # backwards rate constants for reversible reactions calcualted with the equilibrium constant and forward rate constant
        .add_derived("k_cytb6f_bw", fn=ratio, args=["k_cytb6f", "K_eq_cytb6f"])

        # pools of compunds that have more than 1 state but total value is constant
        # for redox pairs this is oxidized + reduced
        # for ATP/ADP/Pi this is ATP + ADP 
        .add_derived("Fd_ox", fn=diff, args=["Fd_tot", "Fd_red"])
        .add_derived("NADP_ox", fn=diff, args=["NADP_tot", "NADP_red"])
        .add_derived("P_ox", fn=diff, args=["P_tot", "P_red"])
        .add_derived("AP_spent", fn=diff, args=["AP_tot", "ATP"])
        .add_derived("ATPase_inactive", fn=inverse, args=["ATPase_active"])
        .add_derived("PC_ox", fn=diff, args=["PC_tot", "PC_red"])

        # values for convenience kinetics
        .add_derived("f_ox_conv", fn=ratio, args=["Fd_ox", "K_M_F"])
        .add_derived("f_red_conv", fn=ratio, args=["Fd_red", "K_M_F"])
        .add_derived("n_ox_conv", fn=ratio, args=["NADP_ox", "K_M_N"])
        .add_derived("n_red_conv", fn=ratio, args=["NADP_red", "K_M_N"])

        # pH and H+ concentration conversions
        .add_derived("H_diff_absolute", fn=ph_diff, args=["pH_lumen", "pH_stroma"])
        .add_derived("pH_lumen", fn=calc_ph, args=["H_lumen"])

        # activation rates of photosystems
        .add_derived(
            "sigma_PSII",
            fn=sigma_psii,
            args=["sigma_0_PSI", "sigma_0_PSII", "L"],
        )
        .add_derived(
            "sigma_PSI",
            fn=sigma_psi,
            args=["sigma_0_PSI", "sigma_0_PSII", "prob_qt", "L"],
        )
        .add_derived("k_LI", fn=k_li, args=["conversion_factor", "sigma_PSI", "PFD"])
        .add_derived("k_LII", fn=k_lii, args=["conversion_factor", "sigma_PSII", "PFD"])

        # ratio of active/valid photosystems for easier plotting
        .add_derived("B_active_ratio", fn=ratio, args=["B_active", "PSII_tot"])
        .add_derived(
            "Y_0",
            fn=y_0,
            args=[
                "PSI_tot",
                "k_LI",
                "k_Fd_red",
                "Fd_ox",
                "Fd_red",
                "PC_red",
                "PC_ox",
                "K_eq_p700_C",
                "K_eq_p700_Fd",
                "k_PC_ox",
            ],
        )
    )
