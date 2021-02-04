#!/usr/bin/env python
# coding: utf-8

# # Basic chemical, electrical, and thermodynamic principles
# 
# To develop a quantitative understanding of how these processes work, we start with a set of definitions of the basic quantities and concepts with which we are concerned.
# 
# ```{figure} Figure1.png 
# ------
# name: mitofig
# ------
# Diagram of a mitochondrion with the cytosol, intermembrane space (IMS), and matrix indicated. *Inset from left to right:* Protein channels and complexes associated with oxidative phosphorylation in the cristae of the mitochondrion. Complex I (C1) catalyzes the oxidation of NADH$^{2-}$ to NAD$^{-}$ and reduction of ubiquinone (Q) to QH$_2$. Complex II (C2) catalyzes the oxidation of FADH$_2$ to FAD coupled to the reduction of Q. Complex III (C3) catalyzes the oxidation of QH$_2$ coupled to the reduction of cytochrome c (Cyt c). Complex IV (C4) catalyzes the oxidation of Cyt c coupled to the reduction of oxygen to water. These redox transfers drive pumping of H$^+$ ions out of the matrix, establishing the proton motive force across the inner mitochondrial membrane (IMM) that drives ATP synthesis at complex V, or the F$_0$F$_1$-ATPase (F$_0$F$_1$). The adenine nucleotide translocase (ANT) exchanges matrix ATP for IMS ADP. The inorganic phosphate cotransporter (PiC) brings protons and Pi from the IMS to the matrix. Lastly, there is a passive H$^{+}$ leak across the IMM. (Figure created with Biorender.com.)
# ```

# ## Mitochondrial anatomy
# 
# The mitochondrion is a membrane-bound, rod-shaped organelle that is responsible for generating most of the chemical energy needed to power the cell's biochemical reactions by respiration {cite}`Nicholls2013`. Mitochondria are comprised of an outer and inner membrane that are separated by the intermembrane space (IMS) ({numref}`mitofig`). The outer mitochondrial membrane is freely permeable to small molecules and ions. The IMM folds inward to make cristae that extend into the matrix. Transmembrane channels called porins and the respiratory complexes involved in oxidative phosphorylation and ATP synthesis allow for more selective IMM permeability. The IMM encloses the mitochondrial matrix, which contains mitochondrial deoxyribonucleic acid (DNA), the majority of mitochondrial proteins, soluble metabolic intermediates including ATP, ADP, and Pi, and the enzymes catalyzing the tricarboxylic acid (TCA) cycle and $\beta$-oxidation. 

# ## IMM capacitance
# 
# The IMM acts as an electrical capacitor to store energy in an electrostatic potential difference between the milieu on each side. Electrical capacitance of a membrane ($C_m$) is the proportionality between the rate of charge transport across the membrane, i.e. current ($I$), to the rate of membrane potential ($\Delta \Psi$) change, that is, 
# ```{math}
#     C_m \dfrac{ {\rm d} {\Delta\Psi}}{{\rm d} t} = I.
# ```
# In the model and associated calculations presented below, we express fluxes in units of moles per unit time per unit volume of mitochondria. Thus, it is convenient to obtain an estimate of $C_m$ in units of mole per volt per volume of mitochondria. Approximating a mitochondrion as a sphere with radius $r = 1 \ \mu\text{m}$, we obtain a surface area-to-volume ratio of $3 \ \mu\text{m}^{-1}$. Furthermore, we estimate that the IMM has ten-fold greater surface area than the outer membrane, yielding a surface area to volume ratio of $30 \ \mu\text{m}^{-1}$ for the IMM. Since the capacitance density of biological membranes ranges from $0.5\text{-}1.0 \mu\text{F cm}^{-2}$, or $0.5 \text{-} 1.0 \times \ 10^{-8} \ \mu\text{F} \ \mu\text{m}^{-2}$ {cite}`Nicholls2013`, $C_m$ is approximately $3 \times 10^{-8} \ \mu\text{F} \ \mu\text{m}^{-3} = 300 \ \text{F (L mito)}^{-1}$. To convert to the units used in the calculations below, we have 
# ```{math}
#   C_m = 300 \ \frac{\rm F}{\rm L \ mito} = 300 \ \frac{\rm C}{\rm V \cdot L \, mito}\cdot
#   \frac{1}{F}\, \frac{\rm mol}{\rm C} =
#   3.1 \times 10^{-3} \, 
#   \frac{\rm mol}{\rm V \cdot  L \, mito}, \,
# ```
# where $F = 96,485 \ \text{C mol}^{-1}$ is Faraday's constant. 

# ## Gibbs free energy
# 
# A *free energy* is a thermodynamic quantity that relates a change in the thermodynamic state of a system to an associated change in total entropy of the system plus its environment. Chemical reaction processes necessarily proceed in the direction associated with a reduction in free energy {cite}`Nicholls2013`. When free energy of a system is reduced, total entropy (of the universe) is increased. The form of free energy that is operative in constant-temperature and constant-pressure systems (most relevant for biochemistry) is the Gibbs free energy, or simply the *Gibbs energy*. 
# 
# For a chemical reaction of reactants $A_i$ and products $B_j$, 
# ```{math}
#     \sum_{i = 1}^M m_i A_i \rightleftharpoons \sum_{j = 1}^N n_j B_j
# ```
# where $M$ and $N$ are the total number of reactants and products, respectively, and $m_i$ and $n_j$ are the coefficients of reactant $i$ and product $j$, respectively, the Gibbs energy can be expressed as
# ```{math}
# :label: Delta_rG     
#     \Delta_r G = \Delta_r G^\circ + R{\rm T} \ln \left( \dfrac{ \prod_{i = 1}^{N} [\text{B}_j]^{n_i}}{ \prod_{i = 1}^{M} [\text{A}_i]^{m_i}} \right), 
# ```   
# where $\Delta_r G^\circ$ is the reference Gibbs energy for the reaction (a constant at given constant chemical conditions of temperature, pressure, ionic conditions, etc.), $R = 8.314 \ \text{J mol}^{-1} \ \text{K}^{-1}$ is the gas constant, and $\text{T} = 310.15 \ \text{K}$ is the temperature. The second term on the right hand side of Equation {eq}`Delta_rG` governs how changes in concentrations of species affects $\Delta_r G$.
# 
# A system is in chemical equilibrium when there is no thermodynamic driving force, that is, $\Delta_r G = 0$. Thus, for this chemical reaction the reference Gibbs energy is related to the equilibrium constant as
# ```{math} 
# :label: eq:equilibrium
#     K_{eq} = \left( \frac{\prod_{i = 1}^{N} [\text{B}_j]^{n_i}}{\prod_{i = 1}^{M} [\text{A}_i]^{m_i}} \right)_{eq}
#            = \exp\left\{ -\frac{\Delta_r G^\circ}{R{\rm T}} \right\} .
# ```

# ## Membrane potential and proton motive force
# 
# Free energy associated with the oxidation of primary fuels is transduced to generate the chemical potential across the IMM known as the *proton motive force*, which is used to synthesize ATP in the matrix and transport ATP out of the matrix to the cytosol {cite}`Nicholls2013`. The thermodynamic driving force for translocation of hydrogen ions ($\text{H}^{+}$) across the IMM has two components: the difference in electrostatic potential across the membrane, $\Delta\Psi$ (V), and the difference in $\text{H}^{+}$ concentration (or activity) between the media on either side of the membrane, $\Delta\text{pH}$, that is 
# ```{math}
# :label: DG_H
#     \Delta G_{\rm H} &=& -F\Delta\Psi + R{\rm T}\ln\left( [{\rm H}^+]_x/[{\rm H}^+]_c \right)  \nonumber \\
#     &=&  -F\Delta\Psi - 2.3 R{\rm T} \, \Delta{\rm pH},
# ```
# where the subscripts $x$ and $c$ indicate matrix and external (cytosol) spaces. $\Delta\Psi$ is defined as the cytosolic potential minus matrix potential, yielding a negative change in free energy for a positive potential. Membrane potential in respiring mitochondria is approximately $150 \text{-} 200 \ \text{mV}$, yielding a contribution to $\Delta G_{\rm H}$ on the order of $15 \text{-} 20 \ \text{kJ mol}^{-1}$ {cite}`Bazil2016`. Under in vitro conditions, $\Delta\text{pH}$ between the matrix and external buffer is on the order of $0.1 \ \text{pH}$ units {cite}`Bazil2016`. Thus, the contribution to proton motive force from a pH difference is less than $1 \ \text{kJ mol}^{-1}$ and substantially smaller than that from $\Delta\Psi$.

# ## Thermodynamics of ATP synthesis/hydrolysis
# 
# Under physiological conditions the ATP hydrolysis reaction
# ```{math} 
# :label: ATP1
# 	\text{ATP}^{4-} + \text{H}_2\text{O} \rightleftharpoons 
# 	    \text{ADP}^{3-} + \text{HPO}_4^{2-} + \text{H}^{+} 
# ```
# is thermodynamically favored to proceed from the left-to-right direction. The Gibbs energy associated with turnover of this reaction is 
# ```{math}
# :label: DrG_ATP
# 	\Delta_r G_{\rm ATP} = \Delta_r G^o_\text{ATP} + R{\rm T} \ln 
# 	\left( \frac{ [\text{ADP}^{3-}] [\text{HPO}_4^{2-}] [{\rm H}^{+}] }
# 	{ [\text{ATP}^{4-}] }\right),
# ```
# where the Gibbs energy for ATP hydrolysis under physiological conditions is approximately $\Delta_r G^o_\text{ATP} = 4.99 \ \text{kJ mol}^{-1}$ {cite}`Li2011`.

# ### Calculation of the ATP hydrolysis potential
# 
# Equation {eq}`DrG_ATP` expresses the Gibbs energy of chemical Equation {eq}`ATP1` in terms of its *chemical species*. In practice, biochemistry typically deals with biochemical *reactants*, which are comprised of sums of rapidly interconverting chemical species. We calculate the total ATP concentration, $[\Sigma \text{ATP}]$, in terms of its bound and unbound species, that is, 
# ```{math}  
# :label: sumATP
# 	[\Sigma \text{ATP}] &=& [\text{ATP}^{4-}] + [\text{MgATP}^{2-}] + [\text{HATP}^{3-}] + [\text{KATP}^{3-}] \nonumber\\
# 	&=& [\text{ATP}^{4-}] + \frac{[\text{Mg}^{2+}] [\text{ATP}^{4-}]}{K_{\text{MgATP}}} + \frac{ [\text{H}^{+}] [\text{ATP}^{4-}]}{K_{\text{HATP}}} + \frac{ [\text{K}^{+}] [\text{ATP}^{4-}]}{K_{\text{KATP}}} \nonumber \\
# 	&=& [\text{ATP}^{4-}] \left( 1 + \frac{[\text{Mg}^{2+}]}{K_{\text{MgATP}}} + \frac{ [\text{H}^{+}]}{K_{\text{HATP}}} + \frac{ [\text{K}^{+}]}{K_{\text{KATP}}} \right) \nonumber \\
# 	&=& [\text{ATP}^{4-}] P_{\text{ATP}},
# ```
# where $P_{\text{ATP}}$ is a *binding polynomial*. Here, we we account for only the single cation-bound species. (Free $\text{H}^+$ in solution associates with water to form $\text{H}_3\text{O}^+$. Here we use [$\text{H}^+$] to indicate hydrogen ion activity, which is equal to $10^{-\text{pH}}$.)  {numref}`table-dissociationconstants` lists the dissociation constants used in this study from {cite}`Li2011`. Similarly, total ADP, [$\Sigma \text{ADP}$], and inorganic phosphate, [$\Sigma \text{Pi}$], concentrations are 
# ```{math} 
# :label: sumADP
#     [\Sigma {\rm ADP} ] &=& [{\rm ADP}^{3-}]\left( 1 + \frac{[{\rm Mg}^{2+}]}{K_{\rm MgADP}} + \frac{ [{\rm H}^{+}]}{K_{\rm HADP}} + \frac{ [{\rm K}^{+}]}{K_{\rm KADP}} \right) \nonumber \\
#     &=& [{\rm ADP}^{3-}]P_{\rm ADP} 
# ```
# and
# ```{math}
# :label: sumPi
#     [\Sigma {\rm Pi} ] &=& [{\rm HPO}_4^{2-}] \left( 1 + \frac{[{\rm Mg}^{2+}]}{K_{\rm MgPi}} + \frac{ [{\rm H}^{+}]}{K_{\rm HPi}} + \frac{ [{\rm K}^{+}]}{K_{\rm KPi}} \right) \nonumber \\
#     &=& [{\rm HPO}_4^{2-}] P_{\rm Pi},
# ```
# for binding polynomials $P_{\text{ADP}}$ and $P_{\text{Pi}}$. 
# 
# Expressing the Gibbs energy of ATP hydrolysis in Equation {eq}`ATP1` in terms of biochemical reactant concentrations, we obtain
# ```{math} 
# :label: ATP2
#     \Delta_r G_{\rm ATP} &=& \Delta_r G^o_\text{ATP} + R{\rm T} \ln \left(
#     \frac{[\Sigma{\rm ADP}][\Sigma{\rm Pi}]}
#     {[\Sigma{\rm ATP}]}\cdot\frac{[{\rm H}^+]P_{\rm ATP}}{P_{\rm ADP}P_{\rm Pi}}
#     \right) \nonumber \\ 
#     &=& \Delta_r G^o_\text{ATP}
#     + R{\rm T} \ln \left(\frac{[{\rm H}^+]P_{\rm ATP}}{P_{\rm ADP}P_{\rm Pi}} \right)
#     + R{\rm T} \ln \left(\frac{[\Sigma{\rm ADP}][\Sigma{\rm Pi}]}
#     {[\Sigma{\rm ATP}]}\right) \nonumber \\ 
#     &=& \Delta_r G'^o_\text{ATP}
#     + R{\rm T} \ln \left(\frac{[\Sigma{\rm ADP}][\Sigma{\rm Pi}]}
#     {[\Sigma{\rm ATP}]}\right)
# ```
# where $\Delta_r G'^o_\text{ATP}$ is a transformed, or *apparent*, reference Gibbs energy for the reaction. 

# ```{list-table} Dissociation constants given as 10$^{-\text{p}K_a}$.
# :header-rows: 2
# :name: table-dissociationconstants
# 
# * - 
#   - 
#   - Ligand ($L$)
#   - 
# * -
#   - Mg$^{2+}$ 
#   - H$^{+}$ 
#   -	K$^{+}$	
# * - $K_{L-\text{ATP}}$ 
#   - $10^{-3.88}$ 
#   - $10^{-6.33}$ 
#   - $10^{-1.02}$ 
# * - $K_{L-\text{ADP}}$ 
#   - $10^{-3.00}$ 
#   - $10^{-6.26}$ 
#   - $10^{-0.89}$  
# * - $K_{L-\text{Pi}}$ 
#   - $10^{-1.66}$ 
#   - $10^{-6.62}$ 
#   - $10^{-0.42}$  
# ```

# The following code computes the apparent Gibbs energy with $\text{pH} = 7$, $[\text{K}^{+}] = 150 \ \text{mM}$, and $[\text{Mg}^{2+}] = 1 \ \text{mM}$. Biochemical reactant concentrations are set such that the total adenine nucleotide (TAN) pool inside the mitochondrion is $10 \ \text{mM}$, $[\Sigma \text{ATP}] = 0.5 \ \text{mM}$, $[\Sigma \text{ADP}] = 9.5 \ \text{mM}$, and $[\Sigma \text{Pi}] = 1 \ \text{mM}$. Here, we obtain a value of approximately $\text{-}45 \ \text{kJ mol}^{-1}$.  

# In[1]:


# Import numpy package for calculations 
import numpy as np

# Dissociation constants
K_MgATP = 10**(-3.88)
K_MgADP = 10**(-3.00)
K_MgPi  = 10**(-1.66)
K_HATP  = 10**(-6.33)
K_HADP  = 10**(-6.26)
K_HPi   = 10**(-6.62)
K_KATP  = 10**(-1.02)
K_KADP  = 10**(-0.89)
K_KPi   = 10**(-0.42)

#  Gibbs energy under physiological conditions(J mol^(-1))
DrGo_ATP = 4990

# Thermochemical constants
R = 8.314           # J (mol * K)**(-1)
T = 310.15          # K
F = 96485           # C mol**(-1)

# Environment concentrations 
pH = 7
H  = 10**(-pH)      # Molar 
K  = 150e-3         # Molar 
Mg = 1e-3           # Molar 

# Binding polynomials
P_ATP = 1 + H/K_HATP + K/K_KATP + Mg/K_MgATP # equation 6
P_ADP = 1 + H/K_HADP + K/K_KADP + Mg/K_MgADP # equation 7 
P_Pi  = 1 + H/K_HPi  + K/K_KPi  + Mg/K_MgPi  # equation 8 

# Total measureable concentrations 
sumATP = 0.5e-3         # Molar
sumADP = 9.5e-3         # Molar
sumPi  = 1.0e-3         # Molar

# Reaction:
# ATP4− + H2O ⇌ ADP3− + HPO2−4 + H+

# Use equation 9 to calcuate apparent Gibbs energy 
DrG_ATP_apparent = DrGo_ATP + R * T * np.log(H * P_ATP / (P_ADP * P_Pi))

# Use equation 9 to calculate actual Gibbs energy 
DrG_ATP = DrG_ATP_apparent + R * T * np.log((sumADP * sumPi / sumATP))

print('Gibbs energy of ATP hydrolysis (kJ mol^(-1))')
print(DrG_ATP / 1000) 


# The reactant concentrations used in the above calculation represent reasonable values for concentrations in the mitochondrial matrix. In the cytosol, the ATP/ADP ratio is on the order of 100:1, yielding a $\Delta_r G_\text{ATP}$ of approximately $\text{-}45 \ \text{kJ mol}^{-1}$. 

# ### ATP synthesis in the mitochondrial matrix
# 
# The F$_0$F$_1$ ATP synthase catalyzes the synthesis of ATP from ADP and Pi by coupling to the translocation of  $n_{\text{F}} = 8/3$ protons from the cytosol to the matrix via the combined reaction
# ```{math} 
# :label: ATP3
#     ({\rm ADP}^{3-})_x + ({\rm HPO_4}^{2-})_x + ({\rm H}^+)_x + n_{\text{F}} (\text{H}^{+})_c 
# 	\rightleftharpoons 
# 	({\rm ATP})^{4-}_x + {\rm H_2O} +  n_{\text{F}} (\text{H}^{+})_x \, .
# ``` 
# Using the Gibbs energy of the reaction of Equation {eq}`ATP2` and the proton motive force in Equation {eq}`DG_H`, the overall Gibbs energy for the coupled process of ATP synthesis and proton transport via the F$_0$F$_1$ ATP synthase is 
# ```{math} 
# :label: DG_F
# 	\Delta G_{\text{F}} &=& -\Delta_r G_{\rm ATP} + n_\text{F} \Delta G_{\rm H} \nonumber \\
# 	&=& -\Delta_r G'^o_\text{ATP} - R{\rm T} \ln \left(\frac{[\Sigma{\rm ADP}]_x[\Sigma{\rm Pi}]_x}
#     {[\Sigma{\rm ATP}]_x}\right) - n_\text{F} F \Delta \Psi + R{\rm T} \ln \left( 
#     \frac{ [{\rm H}^{+}]_x }{ [{\rm H}^{+}]_c } \right)^{n_{\rm F}} . 
# ```
# Note that the negative before $\Delta_r G_\text{ATP}$ indicates that the reaction of Equation {eq}`ATP1` is reversed in Equation {eq}`ATP3`. The equilibrium concentration ratio occurs when $\Delta G_{\text{F}} = 0$. Solving for the second term in Equation {eq}`DG_F`, we calculate the apparent equilibrium constant for ATP synthesis as 
# ```{math} 
# 	K_{eq,\text{F}}^\prime = 
# 	    \left( \frac{[\Sigma{\rm ATP}]_x}{[\Sigma{\rm ADP}]_x[\Sigma{\rm Pi}]_x} \right)_{eq}
# 	    = \exp\left\{\frac{ \Delta_rG'^o_{\rm ATP}  + n_{\rm F} F \Delta\Psi}{R{\rm T}}\right\}
#         \left( \frac{[{\rm H^+}]_c}{[{\rm H^+}]_x} \right)^{n_{\rm F}}. 
# ```

# (modelATPsynthesis)= 
# ### Mathematical modeling ATP synthesis
# 
# A simple model of ATP synthesis kinetics can be constructed using the apparent equilibrium constant and mass-action kinetics in the form
# ```{math} 
#     J_{\text{F}} = X_{\text{F}} (K_{eq,\text{F}}^\prime [\Sigma \text{ADP}]_x [\Sigma \text{Pi}]_x - [\Sigma \text{ATP}]_x), 
# ``` 
# where $X_{\text{F}} = 1000 \ \text{mol s}^{-1} \ \text{(L mito)}^{-1}$ is a rate constant set to an arbitrarily high value that maintains the reaction in equilibrium in model simulations. To simulate ATP synthesis at a given membrane potential, matrix pH, cytosolic pH, and cation concentrations, we have
# ```{math} 
# :label: system-ATPase
#     \left\{       
#         \renewcommand{\arraystretch}{2}
#         \begin{array}{rl}
#             \dfrac{ {\rm d} [\Sigma \text{ATP}]_x }{{\rm d} t} &= J_\text{F} / W_x  \\
#             \dfrac{ {\rm d} [\Sigma \text{ADP}]_x }{{\rm d} t} &= -J_\text{F} / W_x  \\ 
#             \dfrac{ {\rm d} [\Sigma \text{Pi}]_x }{{\rm d} t}  &= -J_\text{F} / W_x, 
#         \end{array} 
#         \renewcommand{\arraystretch}{1}
#     \right. 
# ```
# where $W_x \ \text{((L matrix water) (L mito)}^{-1}$) is the fraction of water volume in the mitochondrial matrix to total volume of the mitochondrion. Dissociation constants are listed in {numref}`table-dissociationconstants` and all other parameters are listed in {numref}`table-biophysicalconstants`. 

# ```{list-table} Parameters for ATP synthesis in vitro.
# :header-rows: 1
# :name: table-biophysicalconstants
# 
# * - Symbol 
#   - Units 
#   - Description 
#   - Value 
#   - Source
# * - F$_0$F$_1$ ATP synthase constants 
#   -
#   -
#   -
#   -
# * - $n_{\text{F}}$ 
#   - 
#   - Protons translocated  
#   - $8/3 $
#   - {cite}`Nicholls2013` 
# * - $X_\text{F}$ 
#   - mol s$^{-1}$ (L mito)$^{-1}$ 
#   - Rate constant 
#   - $1000 $
#   - 
# * - $\Delta_r G_\text{ATP}^\circ$ 
#   - kJ mol$^{-1}$ 
#   - Reference Gibbs energy  
#   - $4.99 $
#   - {cite}`Li2011` 
# * - Biophysical constants 
#   -
#   -
#   -
#   -
# * - $R$ 
#   - J mol$^{-1}$ K$^{-1}$ 
#   - Gas constant 
#   - $8.314 $
#   - 
# * - $T$ 
#   - K 
#   - Temperature 
#   - $310.15 $
#   - 
# * - $F$ 
#   - C mol$^{-1}$ 
#   - Faraday's constant 
#   - $96485$ 
#   - 
# * - $C_m$ 
#   - mol V$^{-1}$ (L mito)$^{-1}$ 
#   - IMM capacitance 
#   - $3.1\text{e-}3$ 
#   - {cite}`Beard2005` 
# * - Volume ratios 
#   -
#   -
#   -
#   -
# * - $V_c$ 
#   - (L cyto) (L cell)$^{-1}$ 
#   - Cyto to cell ratio 
#   - $0.6601$
#   - {cite}`Bazil2016`
# * - $V_m$ 
#   - (L mito) (L cell)$^{-1}$ 
#   - Mito to cell ratio 
#   - $0.2882$ 
#   - {cite}`Bazil2016`
# * - $V_{m2c}$ 
#   - (L mito) (L cyto)$^{-1}$ 
#   - Mito to cyto ratio 
#   - $V_m / V_c$ 
#   - 
# * - $W_c$ 
#   - (L cyto water) (L cyto)$^{-1}$ 
#   - Cyto water space ratio 
#   - $0.8425$ 
#   - {cite}`Bazil2016` 
# * - $W_m$ 
#   - (L mito water) (L mito)$^{-1}$ 
#   - Mito water space ratio 
#   - $0.7238 $
#   - {cite}`Bazil2016` 
# * -	$W_x$ 
#   - (L matrix water) (L mito)$^{-1}$ 
#   - Mito matrix water space ratio 
#   - $0.9$ $W_m$
#   - {cite}`Bazil2016` 
# * - $W_i$ 
#   - (L IM water) (L mito)$^{-1}$ 
#   - IMS water space ratio 
#   - $0.1$ $W_m$ 
#   - {cite}`Bazil2016` 
# ```

# The following code simulates steady state ATP, ADP, and Pi concentrations for $\Delta \Psi = 175 \ \text{mV}$. Here, a pH gradient is fixed across the IMM such that the pH in the matrix is slightly more basic than the cytosol, $\text{pH}_x = 7.4$ and $\text{pH}_c = 7.2$. All other conditions remain unchanged. 

# In[2]:


from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np

# Define system of ordinary differential equations from equation (12)
def dXdt(t, X, DPsi, pH_c):
    # Unpack X state variable
    sumATP, sumADP, sumPi = X
    
    # Biophysical constants 
    R   = 8.314          # J (mol * K)**(-1)
    T   = 310.15         # K
    F   = 96485          # C mol**(-1)
    C_m = 3.1e-3         # mol (V * L mito)**(-1)
    
    # F0F1 constants 
    n_F = 8/3
    X_F = 1000            # mol (s * L mito)**(-1)
    
    # Dissociation constants
    K_MgATP = 10**(-3.88)
    K_MgADP = 10**(-3.00)
    K_MgPi  = 10**(-1.66)
    K_HATP  = 10**(-6.33)
    K_HADP  = 10**(-6.26)
    K_HPi   = 10**(-6.62)
    K_KATP  = 10**(-1.02)
    K_KADP  = 10**(-0.89)
    K_KPi   = 10**(-0.42)

    # Environment concentrations 
    pH_x = 7.4          # pH in matrix
    H_x  = 10**(-pH_x)  # M 
    H_c  = 10**(-pH_c)  # M 
    K_x  = 150e-3       # M 
    Mg_x = 1e-3         # M 

    # Binding polynomials
    P_ATP = 1 + H_x/K_HATP + K_x/K_KATP + Mg_x/K_MgATP # equation 6
    P_ADP = 1 + H_x/K_HADP + K_x/K_KADP + Mg_x/K_MgADP # equation 7 
    P_Pi  = 1 + H_x/K_HPi  + K_x/K_KPi  + Mg_x/K_MgPi  # equation 8 

    # Volume ratios
    W_m = 0.7238         # (L mito water) (L mito)**(-1)
    W_x = 0.9 * W_m      # (L matrix water) (L mito)**(-1)
    
    # Gibbs energy (equation 9)
    DrGo_F   = 4990      # (J mol**(-1))
    DrGapp_F = DrGo_F + R * T * np.log(H_x * P_ATP / (P_ADP * P_Pi))
    
    # Apparent equilibrium constant 
    Kapp_F = np.exp((DrGapp_F + n_F * F * DPsi)/ (R * T)) * (H_c / H_x) ** n_F
    
    # Flux (mol (s * L mito)**(-1))  
    J_F = X_F * (Kapp_F * sumADP * sumPi - sumATP)
       
    ###### Differential equations (equation 12) ######
    dATP = J_F / W_x
    dADP = -J_F / W_x
    dPi  = -J_F / W_x
    
    dX = (dATP, dADP, dPi)
    return dX


# Simple steady state simulation at 175 mV membrane potential 

# Initial conditions (M)
sumATP_0 = 0.5e-3
sumADP_0 = 9.5e-3
sumPi_0  = 1e-3

X_0 = np.array([sumATP_0, sumADP_0, sumPi_0])

# Inputs  
DPsi = 175e-3 # Constant membrane potential (V)
pH_c = 7.2    # IMS/buffer pH 

solutions = solve_ivp(dXdt, [0, 1], X_0, method = 'Radau', args = (DPsi,pH_c))
t = solutions.t
results = solutions.y 
results = results * 1000

# Plot figure 
plt.figure()
plt.plot(t, results[0,:], label = '[$\Sigma$ATP]$_x$')
plt.plot(t, results[1,:], label = '[$\Sigma$ADP]$_x$')
plt.plot(t, results[2,:], label = '[$\Sigma$Pi]$_x$')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Concentration (mM)')
plt.ylim(0, 10)
plt.show()


# The above simulation shows that under the clamped pH and $\Delta\Psi$ conditions simulated here, the model is nearly in equilibrium at the initial conditions of ATP, ADP, and Pi concentrations. Most of the adenine nucleotide remains in the form of ADP and the ATP/ADP ratio in the matrix is approximately $1$:$20$, with the inorganic phosphate concentration of approximately $1 \ \text{mM}$.
# 
# To explore how the equilibrium changes with membrane potential, the following code computes the predicted equilibrium steady-state over a ranges of $\Delta\Psi$ from $100$ to $250 \ \text{mV}$ and over a range of matrix pH from $6.5$ to $8$. 

# In[3]:


### Simulate over a range of Membrane potential from 100 mV to 250 mV ###

# Define array to iterate over
membrane_potential = np.linspace(100,250)    # mV

# Constant external pH
pH_c = 7.2 # IMS/buffer pH

# Define arrays to store steady state results 
ATP_steady_DPsi = np.zeros(len(membrane_potential))
ADP_steady_DPsi = np.zeros(len(membrane_potential))
Pi_steady_DPsi  = np.zeros(len(membrane_potential))

# Iterate through range of membrane potentials 
for i in range(len(membrane_potential)):
    DPsi = membrane_potential[i] / 1000      # convert to V
    temp_results = solve_ivp(dXdt, [0, 5], X_0, method = 'Radau', args = (DPsi, pH_c,)).y*1000  # Concentration in mM
    ATP_steady_DPsi[i] = temp_results[0,-1] 
    ADP_steady_DPsi[i] = temp_results[1,-1] 
    Pi_steady_DPsi[i] = temp_results[2,-1] 
    
# Concentration vs DPsi
plt.figure()
plt.plot(membrane_potential, ATP_steady_DPsi, label = '[$\Sigma$ATP]$_x$')
plt.plot(membrane_potential, ADP_steady_DPsi, label = '[$\Sigma$ADP]$_x$')
plt.plot(membrane_potential, Pi_steady_DPsi, label = '[$\Sigma$Pi]$_x$')
plt.legend()
plt.xlabel('Membrane potential (mV)')
plt.ylabel('Concentration (mM)')
plt.xlim([100, 250])
plt.show()    
    


# The above simulations show that under physiological levels of $\Delta$pH, matrix ATP concentrations become essentially zero for values of the membrane potential less than approximately $150 \ \text{mV}$. At high levels of $\Delta\Psi$ all of the available phosphate is used to phosphorylate ADP to ATP. Since the initial $[\text{Pi}]$ and $[\text{ATP}]$ are $1 \ \text{mM}$ and $0.5 \ \text{mM}$, respectively, the maximum ATP obtained at the maximal $\Delta\Psi$ is $1.5 \ \text{mM}$.

# In[ ]:




