import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# ----------------------------
# Parameter settings
# ----------------------------
T = 25  # Temperature °C
KH = 3.3e-2  # CO2 Gas-liquid partition coefficient mol/(L·atm)
K1 = 4.3e-7  # H2CO3 ⇌ H+ + HCO3-
K2 = 4.8e-11 # HCO3- ⇌ H+ + CO3^2-

# Sample Alkalinity (mol/L)
alk_DI = 0       # DIW
alk_tap = 2.0e-3 # Tap water
alk_waste = 1e-2 # Wastewater

# CO2 partial pressuer (atm)
P_CO2 = 1.0  # under pure CO2 

# ----------------------------
# Charge balance equation of the carbonic acid system
# ----------------------------
def carbon_system(pH, KH, P_CO2, K1, K2, Alk):
    H = 10**(-pH)
    CO2_aq = KH * P_CO2
    HCO3 = K1 * CO2_aq / H
    CO3 = K2 * HCO3 / H
    OH = 1e-14 / H
    # Charge balance：H+ + Other cations = OH- + HCO3- + 2CO3^2- + Alk
    charge_balance = H - (HCO3 + 2*CO3 + OH) + Alk
    return charge_balance

# ----------------------------
# Solve for pH using different initial values ​​to ensure convergence.
# ----------------------------
def solve_pH(Alk):
    if Alk < 1e-4:
        pH_guess = 4.0  # DIW
    elif Alk < 4e-3:
        pH_guess = 5.0  # Tap water
    else:
        pH_guess = 7.0  # Wastewater
    pH_solution = fsolve(carbon_system, pH_guess, args=(KH, P_CO2, K1, K2, Alk))
    return pH_solution[0]

# ----------------------------
# Calculate different water pH
# ----------------------------
pH_DI = solve_pH(alk_DI)
pH_tap = solve_pH(alk_tap)
pH_waste = solve_pH(alk_waste)

# ----------------------------
# print results
# ----------------------------
print(f"DIW balance pH: {pH_DI:.2f}")
print(f"Tap water balance pH: {pH_tap:.2f}")
print(f"Wastewater balance pH : {pH_waste:.2f}")

# ----------------------------
# Visual bar chart
# ----------------------------
waters = ['DIW', 'Tap water', 'Wastewater']
pH_values = [pH_DI, pH_tap, pH_waste]

plt.bar(waters, pH_values, color=['skyblue', 'lightgreen', 'salmon'])
plt.ylabel('Balance pH')
plt.title('pH equilibrium in different water bodies after CO₂ dissolution')
plt.ylim(0, 14)
plt.show()


# ----------------------------
# Coding by Shuwen Zhang : 999023104 / Zijie Zhou : 999023187
# ----------------------------
