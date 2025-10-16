import math

__all__ = [
    "hpr_make_stoich",
    "k_eq_atpsynth",
    "k_eq_fnr",
    "k_eq_qapq",
    "k_eq_cytb6f",
    "k_eq_p700_c",
    "k_eq_p700_fd",
    "psii_make_stoich",
    "y_0",
    "b6f_make_stoich",
    "calc_ph",
    "diff",
    "heaviside",
    "inverse",
    "k_li",
    "k_lii",
    "leak_make_stoich",
    "negative",
    "ratio",
    "sigma_psi",
    "sigma_psii",
    "sum_psii_states",
    "pqred_make_stoich",
    "o2_ox_make_stoich",
    "ph_diff",
    "calc_ph_stroma",
    "sum",
    "charge_separations2",
    "charge_separations",
]


def charge_separations2(v_psi: float, v_psii: float, ps1_tot: float) -> float:
    """Returns the fraction of charge separations by PSII and PSI per total PSI"""
    #conversion factor converts mmol/mol CHl to mol/L
    conv_fact = 350 * 10**(-6) * 10**(-3)
    return ((v_psi + v_psii ) * conv_fact) / (ps1_tot * conv_fact)

def charge_separations(v_psi: float, ps1_tot: float) -> float:
    """Returns the fraction of charge separations by PSI per total PSI"""
    conv_fact = 350 * 10**(-6) * 10**(-3)
    return ((v_psi) * 350 * 10**(-6) * 10**(-3)) / (ps1_tot * conv_fact)

def sum(s1: float, s2: float) -> float:
    """Sums up two arguments"""
    return s1 + s2

def fluorescence(k_f: float,
                 k_h0: float,
                 k_p:float,
                 b0: float,
                 b2: float,
                 sigma_psii: float,
) -> float:
    """fluorescence of PSII"""
    a = k_f / (k_h0 + k_f + k_p)
    b = k_f / (k_h0 + k_f)
    return sigma_psii*(a * b0 + b * b2)

    
def heaviside(pfd: float) -> int:
    """Heaviside step function that returns 1 when arg is larger than 0 and otherwise 0"""
    if pfd > 0:
        return 1
    return 0

# factor of 2.5*10^-4 for calculating for concentration in lumen
def calc_ph(protons: float) -> float:
    """Returns pH in lumen corresponding to relative Proton concentration"""
    return -math.log10(protons * 2.5 * 10 ** (-4))


def ph_diff(ph_lumen: float, ph_stroma: float) -> float:
    """Returns difference in proton concentration between lumen and stroma as absolute values (not dependent on chlorophyll content)"""
    return 10**(-ph_lumen) - 10**(-ph_stroma)


def diff(total: float, subtractor: float) -> float:
    """Subtracts subtractor from total"""
    return total - subtractor

def ratio(numerator: float, denominator: float) -> float:
    """Returns ratio of the two arguments"""
    if denominator == 0:
        print("Zero division Problem encountered, displaying value as 0")
        return 0
    return numerator / denominator


def inverse(value: float) -> float:
    """Returns 1 minus the argument"""
    return 1 - value

def sigma_psi(
    sigma_i_0: float,
    sigma_ii_0: float,
    prob_qt: float,
    l: float,
) -> float:
    """Proportion of antennae on PSI"""
    return sigma_i_0 + prob_qt * (1 - sigma_i_0 - sigma_ii_0) * l


def sigma_psii(
    sigma_i_0: float,
    sigma_ii_0: float,
    l: float,
) -> float:
    """Proportion of antennae on PSII"""
    return sigma_ii_0 + (1 - sigma_i_0 - sigma_ii_0) * (1 - l)


def calc_ph_stroma(h_stroma: float) ->float:
    """Returns pH in stroma corresponding to relative proton concentration"""
    return -math.log10(h_stroma * 3.2 * 10 ** (-5))


def k_li(
    conversion_factor: float,
    sigma_psi: float,
    pfd: float,
) -> float:
    """Activation rate of PSI"""
    return conversion_factor * sigma_psi * pfd


def k_lii(
    conversion_factor: float,
    sigma_psii: float,
    pfd: float,
) -> float:
    """Activation rate of PSII"""
    return conversion_factor * sigma_psii * pfd


def k_eq_qapq(
    f: float,
    e0_qa: float,
    e0_pq: float,
    ph_st: float,
    r: float,
    t: float,
) -> float:
    "Equilibrium constant of quinone reduction"
    DG1 = -f * e0_qa
    DG2 = (-2 * f * e0_pq) + (2 * ph_st * math.log(10) * (r * t))
    DG0 = -2 * DG1 + DG2
    return math.exp(-DG0 / (r * t))


def k_eq_cytb6f(
    ph_lu: float,
    f: float,
    e0_pq: float,
    r: float,
    t: float,
    e0_pc: float,
    ph_st: float,
) -> float:
    """Equilibrium constant of Cytochrome b6f"""
    DG1 = -2 * f * e0_pq + 2 * r * t * math.log(10) * ph_lu
    DG2 = -f * e0_pc
    DG3 = r * t * math.log(10) * (ph_st - ph_lu)
    DG = -DG1 + 2 * DG2 + 2 * DG3
    return math.exp(-DG / (r * t))


def k_eq_atpsynth(
    ph_lu: float,
    g0_atp: float,
    ph_st: float,
    r: float,
    t: float,
    pi: float,
    hpr:float,
) -> float:
    """Equilibrium constant of ATP synthase."""
    RT = r * t
    G1 = g0_atp
    G2 = RT * math.log(10) * (ph_st - ph_lu)
    delta_G = G1 - (hpr) * G2
    return (pi * math.exp(-delta_G / RT))


def k_eq_fnr(
    f: float,
    r: float,
    t: float,
    e0_fd_redox: float,
    e0_nadp_redox: float,
    ph_stroma: float,
) -> float:
    "Equilibrium constant for Ferredoxin NADPH reductase"
    G0_1 = -f * e0_fd_redox
    G0_2 = -2 * f * e0_nadp_redox
    G0_3 = r * t * math.log(10) * (ph_stroma)
    delta_G = -2 * G0_1 + G0_2 + G0_3
    return math.exp(-delta_G / (r * t))


def k_eq_p700_fd(
    f: float,
    r: float,
    t: float,
    e0_fd_redox: float,
    e0_fa_redox: float
) -> float:
    """Equilibrium constant for ferredoxin reduction at PSI"""
    DG1 = -e0_fa_redox* f
    DG2 = -e0_fd_redox * f
    DG = -DG1+DG2
    return(math.exp(-DG/(r * t)))
    

def k_eq_p700_c(
    f: float,
    r: float,
    t: float,
    e0_pc_redox: float,
    e0_p700_redox: float
) -> float:
    """Equilibrium constant for Plastocyanin oxidation at PSI"""
    DG1 = -e0_pc_redox * f
    DG2 = -e0_p700_redox * f
    DG = -DG1+DG2
    K_eq = math.exp(-DG/(r * t)) 
    return(K_eq)


def y_0(
    psi_tot: float,
    k_li: float,
    k_fd_red: float,
    fd_ox: float,
    fd_red: float,
    pc_red: float,
    pc_ox: float,
    k_eq_p700_c: float,
    k_eq_p700_fd: float,
    k_pc_ox: float,
) -> float:
    """Returns total amount of PSI in ground state"""
    a1 = 1 + (k_li / (k_fd_red * fd_ox))
    a2 = 1 + (fd_ox / k_eq_p700_fd * fd_red)
    a3 = (pc_red / (k_eq_p700_c * pc_ox)) + (k_li / (k_pc_ox * pc_red))
    return psi_tot / (a1 + a2 * a3)


#only relevant for dynamic PSII
def sum_psii_states(b0: float,
                    b1: float,
                    b2: float,
                    b3: float,
                    ) -> float:
    """Sums up all four active states of PSII"""
    return(b0+b1+b2+b3)


#stiochiometric factors for rates concerning H values inclduing buffering
def hpr_make_stoich(hpr: float, bh: float) -> float:
    """stoichiometric factor with buffering taken into account for protons transported through ATP synthase"""
    return -hpr / bh

def two_div(denominator: float) -> float:
    """Returns 2 divided by the argument"""
    if denominator == 0:
        print("Zero division Problem encountered, displaying value as 0")
        return 0
    return 2 / denominator

def neg_two_div(denominator: float) -> float:
    """Returns -2 divided by the argument"""
    if denominator == 0:
        print("Zero division Problem encountered, displaying value as 0")
        return 0
    return -2 / denominator

def four_div(denominator: float) -> float:
    """Returns 4 divided by the argument"""
    if denominator == 0:
        print("Zero division Problem encountered, displaying value as 0")
        return 0
    return 4 / denominator

def neg_four_div(denominator: float) -> float:
    """Returns -4 divided by the argument"""
    if denominator == 0:
        print("Zero division Problem encountered, displaying value as 0")
        return 0
    return -4 / denominator

def neg_one_div(denominator: float) -> float:
    """Returns -1 divided by the argument"""
    if denominator == 0:
        print("Zero division Problem encountered, displaying value as 0")
        return 0
    return -1 / denominator


def hpr_dy_stoich(hpr: float, f: float, c:float) -> float:
    """stoichiometric factor for delta psi of HPR protons crossing membrane to stroma"""
    return -(f * hpr)/c

def psii_dy_stoich(f: float, c:float) -> float:
    """stoichiometric factor for delta psi of two protons crossing membrane to lumen"""
    return (f * 2)/c

def b6f_dy_stoich(f: float, c:float) -> float:
    """stoichiometric factor for delta psi of four protons crossing membrane to lumen"""
    return (f * 4)/c

def leak_dy_stoich(f: float, c:float) -> float:
    """stoichiometric factor for delta psi of one proton crossing membrane to stroma"""
    return -(f * 1)/c