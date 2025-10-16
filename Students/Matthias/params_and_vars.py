from mxlpy import Model

from derived_quantities import add_derived_quantities
from rates import add_rates_to_model
from photosystem import add_photosystem
from dynamic_photosystem import add_dynamic_psii

def get_model(ps_dynamic: bool = True) -> Model:
    """Returns the full model including all parameters, variables, derived quantities and rates."""
    model = (
        Model()

        # parameters of total pool sizes  unit = mmol/mol Chl
        .add_parameters(
            {
                "PSII_tot": 2.5,
                "PSI_tot": 2.5,
                "P_tot": 17.5,
                "PC_tot": 4,
                "Fd_tot": 5,
                "AP_tot": 60,
                "NADP_tot": 25,
                "O2_ex": 8,
                "Pi_mol": 0.01,
            }
        )

        # rate constants
        .add_parameters(
            {
                "k_P_red": 250,  # mmol^-1 (mol Chl) s^-1
                "k_PC_ox": 2500,  # mmol^-1 (mol Chl) s^-1
                "k_Fd_red": 2.5 * 10**5,  # mmol^-1 (mol Chl) s^-1
                "k_cytb6f": 2.5,  # mmol^-2 (mol Chl)^2s^-1
                "k_actATPase": 0.05,  # s^-1
                "k_deactATPase": 0.002,  # s^-1
                "k_ATPsynthase": 20,  # s^-1
                "k_ATPconsumption": 10,  # s^-1
                "k_NADPHconsumption": 15,
                "k_H": 5 * 10**9,  # s^-1
                "k_F": 6.25 * 10**8,  # s^-1
                "k_P": 5 * 10**9,  # s^-1
                "k_NDH": 0.004,  # s^-1
                "v_cytb6f_min": -2.5,  # mmol^-1 (mol Chl) s^-1
                "v_FNR_max": 1500,  # mmol^-1 (mol Chl) s^-1
                "k_FQR": 1,  # mmol^-2 (mol Chl)^2 s^-1
                "k_leak": 0.01,  # s^-1
                "b_H": 100,  # -
                "HPR": 14 / 3,  # -
                "k_PTOX": 0.01,  # mmol^-1(mol Chl)s^-1
                "k_outer_leak": 0.01,
                "conversion_factor": 4,  # should actually be fitted to data but I dont have values for now
            }
        )

        # physical constants
        .add_parameters(
            {
                "F": 96.485,  # C/mol leads to kj/mol when multiplied with E0[V]
                "R": 8.3e-3,  # kJ K^-1 mol^-1
                "T": 298,  # K
            }
        )

        # standard potentials
        # #V*F leads to j/mol, G0_ATP is kj/mol and RT
        .add_parameters(
            {
                "E0_QA_redox": -0.14,  # V
                "E0_PQ_redox": 0.354,  # V
                "E0_PC_redox": 0.38,  # V
                "E0_P700_redox": 0.48,  # V
                "E0_FA_redox": -0.55,  # V
                "E0_Fd_redox": -0.43,
                "E0_NADP_redox": -0.113,  # V
                "G0_ATP": 30.6,  # kj/mol*RT
                "Capac": 4336, #C/V membrane capacitance 
            }
        )

        # Michaelis constants
        .add_parameters(
            {
                "K_M_F": 1.56,  # mmol^-1 (mol Chl) s^-1
                "K_M_N": 0.22,  # mmol^-1 (mol Chl) s^-1
                "K_M_ST": 0.2,  # mmol^-1 (mol Chl) s^-1
                "K_M_fdST": 0.5,  # mmol^-1 (mol Chl) s^-1
            }
        )

        # parameters associated with state transitions
        .add_parameters(
            {
                "k_Stt7": 0.0035,  # s^-1
                "k_Pph1": 0.0013,  # s^-1
                "sigma_0_PSI": 0.37,  # -
                "sigma_0_PSII": 0.1,  # -
                "n_ST": 2,  # -
                "n_fdST": 2,  # -
                "prob_qt": 1,  # -
            }
        )

        # parameters associated with photoinhibition
        .add_parameters(
            {
                "k_deg": 100,  # s^-1
                "k_rep": 5.55 * 10**-4,
            }
        )

        # Xanophyll cycle parameters
        .add_parameters(
            {
                "k_kDeepoxV": 0.0024,  # s^-1
                "k_kEpoxZ": 0.00024,  # s-1
                "K_pHsat": 5.8,  # -
                "nHX": 5,  # -
                "KzSat": 0.12,  # -
            }
        )

        # other parameters
        .add_parameters(
            {
                "PFD": 100,  # umol photons m^-2 s^-1
                "pH_stroma": 7.2,
            }
        )

        .add_variables(
            {
                "P_red": 8.33e-1,  # reduced plastoquinon pool
                "PC_red": 3.99,  # reduced plastocyan pool
                "Fd_red": 3.32e-6,  # reduced Ferredoxin pool
                "ATPase_active": 1.7e-3,  # active ATPase pool
                "ATP": 4.45e-7,  # ATP pool normally 10**(-4)
                "NADP_red": 3.98e-10,  # reduced NADPH pool
                "H_lumen": 1.73e-4, # (10**(-6.5)) / (2.5 * 10**(-4)),  hydrogen pool
                "L": 0,  # pool of phosphorylated antennae
            }
        )
    )
    add_derived_quantities(model)
    add_rates_to_model(model)

    if ps_dynamic:
        add_dynamic_psii(model)
    else:
        add_photosystem(model)
        
    return model
