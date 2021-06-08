# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import tellurium as te
import roadrunner
import numpy as np
r = te.loada("""
    //////////// Constants defining metabolite pools ////////////
    // Volume fractions and water space fractions
    V_c = 1.0          // buffer volume fraction        // L buffer (L cuvette)^(-1)
    V_m = 0.0005       // mitochondrial volume fraction // L mito (L cuvette)^(-1)
    V_m2c = V_m / V_c  // mito to cyto volume ratio     // L mito (L cuvette)^(-1)
    W_c = 1.0          // buffer water space            // L buffer water (L buffer)^(-1)
    W_m = 0.7238       // mitochondrial water space     // L mito water (L mito)^(-1)
    W_x = 0.9*W_m      // matrix water space            // L matrix water (L mito)^(-1)
    W_i = 0.1*W_m      // intermembrane water space     // L IM water (L mito)^(-1)

    // Define Compartments
    compartment matrix = W_x
    compartment ims = W_i
    compartment cyto = W_c / V_m2c
    compartment membrane = C_m

    // Define Species in compartments
    var sumATP_x in matrix
    var sumADP_x in matrix
    var sumPi_x in matrix
    var NADH_x  in matrix
    var NAD_x  in matrix
    var QH2_x  in matrix
    var Q_x  in matrix
    
    var cred_i  in ims
    var cox_i  in ims

    var sumATP_c  in cyto
    var sumADP_c  in cyto
    var sumPi_c  in cyto
    var DPsi in membrane

    // Biophysical constants
    R   = 8.314          // J (mol * K)^(-1)
    T   = 310.15         // K
    F   = 96485          // C mol^(-1)
    C_m = 3.1e-3         // mol (V * L mito)^(-1)

    // Total pool concentrations
    NAD_tot = 2.97e-3  // NAD+ and NADH conc            // mol (L matrix water)^(-1)
    Q_tot   = 1.35e-3  // Q and QH2 conc                // mol (L matrix water)^(-1)
    c_tot   = 2.7e-3   // cytochrome c ox and red conc  // mol (L IM water)^(-1)

    //////////// Set fixed pH, cation concentrations, and O2 partial pressure ////////////
    // pH
    pH_x = 7.40
    pH_c = 7.20
    // Hydrogen ion concentration
    H_x = 10^(-pH_x) // mol (L matrix water)^(-1)
    H_c = 10^(-pH_c) // mol (L cuvette water)^(-1)

    // K+ concentrations
    K_x  = 100e-3      // mol (L matrix water)^(-1)
    K_c  = 140e-3      // mol (L cyto water)^(-1)

    // Mg2+ concentrations
    Mg_x = 1.0e-3        // mol (L matrix water)^(-1)
    Mg_c = 1.0e-3        // mol (L cyto water)^(-1)

    // Oxygen partial pressure
    PO2 = 25 // mmHg
    // Oxygen concentration
    a_3  = 1.74e-6   // oxygen solubility in cuvette   // mol (L matrix water * mmHg)^(-1)
    O2_x = a_3*PO2   // mol (L matrix water)^(-1)

    //////////// Activity Parameter vector ////////////
    X_DH  = 0.1732
    X_C1  = 1.0e4
    X_C3  = 1.0e6
    X_C4  = 0.0125
    X_F   = 1.0e3
    E_ANT = 0.325
    E_PiC = 5.0e6
    X_H   = 1.0e3
    X_AtC = 3e-6 // Adjustable Parameter


    // Proton motive force parameters (dimensionless)
    n_F  = 8/3
    n_C1 = 4
    n_C3 = 2
    n_C4 = 4


    // Dissociation constants
    K_MgATP = 10^(-3.88)
    K_MgADP = 10^(-3.00)
    K_MgPi  = 10^(-1.66)
    K_HATP  = 10^(-6.33)
    K_HADP  = 10^(-6.26)
    K_HPi   = 10^(-6.62)
    K_KATP  = 10^(-1.02)
    K_KADP  = 10^(-0.89)
    K_KPi   = 10^(-0.42)


    //// Binding polynomials
    // Matrix species // mol (L mito water)^(-1)
    PATP_x := 1 + H_x/K_HATP + Mg_x/K_MgATP + K_x/K_KATP
    PADP_x := 1 + H_x/K_HADP + Mg_x/K_MgADP + K_x/K_KADP
    PPi_x  := 1 + H_x/K_HPi  + Mg_x/K_MgPi  + K_x/K_KPi

    // Cytosol species // mol (L cuvette water)^(-1)
    PATP_c := 1 + H_c/K_HATP + Mg_c/K_MgATP + K_c/K_KATP
    PADP_c := 1 + H_c/K_HADP + Mg_c/K_MgADP + K_c/K_KADP
    PPi_c  := 1 + H_c/K_HPi  + Mg_c/K_MgPi  + K_c/K_KPi

    //// Unbound species
    // Matrix species
    ATP_x := sumATP_x / PATP_x // [ATP4-]_x
    ADP_x := sumADP_x / PADP_x // [ADP3-]_x
    Pi_x  := sumPi_x  / PPi_x  // [HPO42-]_x

    // Cytosolic species
    ATP_c := sumATP_c / PATP_c // [ATP4-]_c
    ADP_c := sumADP_c / PADP_c // [ADP3-]_c
    Pi_c  := sumPi_c  / PPi_c  // [HPO42-]_c


    //////////// NADH Dehydrogenase ////////////
    // Constants
    r      = 6.8385
    k_Pi1  = 4.659e-4    // mol (L matrix water)^(-1)
    k_Pi2  = 6.578e-4    // mol (L matrix water)^(-1)
    NADH_Dehydronenase: NAD_x => NADH_x; X_DH * (r * NAD_x - NADH_x) * ((1 + sumPi_x / k_Pi1) / (1+sumPi_x / k_Pi2))

    //////////// Complex I ////////////
    DrGo_C1 = -109680
    DrGapp_C1 = DrGo_C1 - R * T * ln(H_x)
    Kapp_C1   := exp( -(DrGapp_C1 + n_C1 * F * DPsi) / (R * T)) * ((H_x / H_c)^n_C1)
    Complex_I: Q_x + NADH_x  => NAD_x + QH2_x + 4 DPsi; X_C1 * (Kapp_C1 * NADH_x * Q_x - NAD_x * QH2_x)
    

    //////////// Complex III ////////////
    // QH2_x + 2cuvetteC(ox)3+_i + 2H+_x <-> Q_x + 2cuvetteC(red)2+_i + 4H+_i + 2DPsi
    // Gibbs energy (J mol^(-1))
    DrGo_C3 = 46690
    DrGapp_C3 = DrGo_C3 + 2 * R * T * ln(H_c)
    // Apparent equilibrium constant
    Kapp_C3   := exp(-(DrGapp_C3 + n_C3 * F * DPsi) / (R * T)) * (H_x / H_c)^n_C3
    // Flux (mol (s * L mito)^(-1))
    Complex_III: 2 cox_i + QH2_x  => Q_x + 2 cred_i + 2 DPsi; X_C3 * (Kapp_C3 * cox_i^2 * QH2_x - cred_i^2 * Q_x)
    
    
    
    //////////// Complex IV ////////////
    // 2 cytoC(red)2+_i + 0.5O2_x + 4H+_x <-> cytoC(ox)3+_x + H2O_x + 2H+_i + 2DPsi
    // Constants
    k_O2 = 1.2e-4      // mol (L matrix water)^(-1)
    // Gibbs energy (J mol^(-1))
    DrGo_C4 = -202160  // J mol^(-1)
    DrGapp_C4 := DrGo_C4 - 2 * R * T * ln(H_c)
    // Apparent equilibrium constant
    Kapp_C4   := exp(-(DrGapp_C4 + n_C4 * F * DPsi) / (R * T)) * (H_x / H_c)^n_C4
    // Flux (mol (s * L mito)^(-1))
    Complex_IV: 2 cred_i  => 2 cox_i + 4 DPsi; X_C4 *(Kapp_C4^0.5 * cred_i * O2_x^0.25 - cox_i) * (1 / (1 + k_O2 / O2_x))
    // Convert Complex IV flux to oxygen consumption 
    J_O2 := Complex_IV / 2 * 60 * 1e9 * 0.0000012232

    //////////// F0F1-ATPase ////////////
    // ADP3-_x + HPO42-_x + H+_x + n_A*H+_i <-> ATP4- + H2O + n_A*H+_x
    // Gibbs energy (J mol^(-1))
    DrGo_F = 4990
    DrGapp_F := DrGo_F + R * T * ln( H_x * PATP_x / (PADP_x * PPi_x))
    // Apparent equilibrium constant
    Kapp_F   := exp( (DrGapp_F + n_F * F * DPsi ) / (R * T)) * (H_c / H_x)^n_F
    // Flux (mol (s * L mito)^(-1))
    F1F0: sumADP_x + sumPi_x + 2.6666667 DPsi => sumATP_x; X_F * (Kapp_F * sumADP_x * sumPi_x - sumATP_x)


    //////////// ANT ////////////
    // ATP4-_x + ADP3-_i <-> ATP4-_i + ADP3-_x
    // Constants
    del_D   = 0.0167
    del_T   = 0.0699
    k2o_ANT = 9.54/60      # s**(-1)
    k3o_ANT = 30.05/60     # s**(-1)
    K0o_D   = 38.89e-6     # mol (L cuvette water)**(-1)
    K0o_T   = 56.05e-6     # mol (L cuvette water)**(-1)
    A       = +0.2829
    B       = -0.2086
    C       = +0.2372
    phi := F * DPsi / (R * T)
    // Reaction rates
    k2_ANT := k2o_ANT * exp((A*(-3) + B*(-4) + C)*phi)
    k3_ANT := k3o_ANT * exp((A*(-4) + B*(-3) + C)*phi)
    // Dissociation constants
    K0_D := K0o_D * exp(3*del_D*phi)
    K0_T := K0o_T * exp(4*del_T*phi)
    q     := k3_ANT * K0_D * exp(phi) / (k2_ANT * K0_T)
    term1 := k2_ANT * ATP_x * ADP_c * q / K0_D
    term2 := k3_ANT * ADP_x * ATP_c / K0_T
    num   := term1 - term2
    den   := (1 + ATP_c/K0_T + ADP_c/K0_D) * (ADP_x + ATP_x * q)
    ANT: sumATP_x + sumADP_c + DPsi => sumATP_c + sumADP_x ; E_ANT * num / den


    //////////// H+-PI2 cotransporter ////////////
    // H2PO42-_x + H+_x = H2PO42-_c + H+_c
    // Constant
    k_PiC = 1.61e-3  // mol (L cuvette)^(-1)
    // H2P04- species
    HPi_c := Pi_c * (H_c / K_HPi)
    HPi_x := Pi_x * (H_x / K_HPi)
    // Flux (mol (s * L mito)^(-1))
    PiC: sumPi_c => sumPi_x; E_PiC * (H_c * HPi_c - H_x * HPi_x) / (k_PiC + HPi_c)

    //////////// H+ leak ////////////
    // Flux (mol (s * L mito)^(-1))
    Leak: DPsi => ; X_H * (H_c * exp(phi/2) - H_x * exp(-phi/2))

    //////////// ATPase ////////////
    // ATP4- + H2O = ADP3- + PI2- + H+
    //Flux (mol (s * L cyto)^(-1))
    sumATP_c => sumADP_c + sumPi_c ; X_AtC / V_m2c /V_c


    
    //////////// Initial Conditions /////////////
    // Membrane Potential 
    DPsi = 175e-3
    // Matrix species
    sumATP_x = 0.5e-3  // mol (L matrix water)**(-1)
    sumADP_x = 9.5e-3  // mol (L matrix water)**(-1)
    sumPi_x  = 0.3e-3  // mol (L matrix water)**(-1)
    NADH_x  = 0.1 * NAD_tot  // mol (L matrix water)**(-1)
    NAD_x = 0.9 * NAD_tot    // mol (L matrix water)**(-1)
    QH2_x  = 0.1 * Q_tot     // mol (L matrix water)**(-1)
    Q_x    = 0.9 * Q_tot     // mol (L matrix water)**(-1)
    
    // IMS species
    cred_i = 0.1 * c_tot // mol (L IMS water)**(-1)
    cox_i = 0.9 * c_tot // mol (L IMS water)**(-1)

    
    // Cytosoo=lic species
    sumATP_c = 5e-3       // mol (L cyto water)^(-1)
    sumADP_c = 0.       // mol (L cyto water)^(-1)
    sumPi_c  = 1.0e-3  // mol (L cyto water)^(-1)
""")


# Save as SMBL 
sbml_file = open('in_vitro.sbml','w')
sbml_file.write(r.getSBML())
sbml_file.close()



X_AtC = np.linspace(0e-6,6e-6, 60)           # Increase max hydrolysis to find apparent Km.
steady_state = np.zeros((len(X_AtC), 4))     # Save steady state concentrations 
JO2 = np.zeros(len(X_AtC))                   # Save JO2 steady state flux 
NADH_ss = np.zeros(len(X_AtC))
DPsi_ss = np.zeros(len(X_AtC))
cred_i_ss = np.zeros(len(X_AtC))
ADP_c_ss = np.zeros(len(X_AtC))


### LowPhosphate Experiment  
# Loop through range of ATPase rates
for i in range(len(X_AtC)):
    r.X_AtC = X_AtC[i]
    r.sumPi_c = 1e-3
    results = r.simulate(0,3000, 500)
    t, ATP_x, ADP_x,  Pi_x,   NADH_x,    NAD_x,     QH2_x,      Q_x,    cred_i,    cox_i, ATP_c,  ADP_c,  Pi_c,   DPsi = results.T
    NADH_ss[i] = NADH_x[-1]
    DPsi_ss[i] = DPsi[-1]
    cred_i_ss[i] =  cred_i[-1]
    ADP_c_ss[i] = ADP_c[-1]
    JO2[i] = r.J_O2
    r.reset()



fig, ax = plt.subplots(2,2, figsize = (8,6))
# Low Pi Plotting  Plotting
ax[0,0].plot(JO2, NADH_ss/r.NAD_tot, 'C0')   # NADH
ax[0,1].plot(JO2, DPsi_ss * 1000, 'C0', label = '1 mM Pi Model')      # Membrane Potential
ax[1,0].plot(JO2, cred_i_ss/r.c_tot,'C0')      # Cytochrome C
ax[1,1].plot(JO2, ADP_c_ss * 1000, 'C0')     # ADP


### High Phosphate Experiment 
# Loop through range of ATPase rates
for i in range(len(X_AtC)):
    r.X_AtC = X_AtC[i]
    r.sumPi_c = 5e-3
    results = r.simulate(0,3000, 500)
    t, ATP_x, ADP_x,  Pi_x,   NADH_x,    NAD_x,     QH2_x,      Q_x,    cred_i,    cox_i, ATP_c,  ADP_c,  Pi_c,   DPsi = results.T
    NADH_ss[i] = NADH_x[-1]
    DPsi_ss[i] = DPsi[-1]
    cred_i_ss[i] =  cred_i[-1]
    ADP_c_ss[i] = ADP_c[-1]
    JO2[i] = r.J_O2
    r.reset()



# High Pi Plotting  Plotting
ax[0,0].plot(JO2, NADH_ss/r.NAD_tot, 'red')   # NADH
ax[0,1].plot(JO2, DPsi_ss * 1000, 'red', label = '5 mM Pi Model')      # Membrane Potential
ax[1,0].plot(JO2, cred_i_ss/r.c_tot,'red')      # Cytochrome C
ax[1,1].plot(JO2, ADP_c_ss * 1000, 'red')     # ADP

ax[1,1].set_ylim([-0.01,0.8])

ax[0,0].set_ylabel('NADH (normalized)')
ax[0,1].set_ylabel('Membrane Potential $\Psi$ (mV)')
ax[1,0].set_xlabel('OCR (nmol O$_2$ min$^{-1}$ U CS$^{-1}$)')
ax[1,0].set_ylabel('Cyt C $^{2+}$ (Normalized)')
ax[1, 1].set_xlabel('OCR (nmol O$_2$ min$^{-1}$ U CS$^{-1}$)')
ax[1, 1].set_ylabel('Buffer ADP (Normalized)')
plt.tight_layout()
ax[0,1].legend()
plt.show()


