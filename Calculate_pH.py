import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# ----------------------------
# 参数设置
# ----------------------------
T = 25  # 温度 °C
KH = 3.3e-2  # CO2 气-液分配系数 mol/(L·atm)
K1 = 4.3e-7  # H2CO3 ⇌ H+ + HCO3-
K2 = 4.8e-11 # HCO3- ⇌ H+ + CO3^2-

# 水体碱度 (mol/L)
alk_DI = 0       # DI 水
alk_tap = 2.0e-3 # 自来水
alk_waste = 1e-2 # 废水 (假设值)

# CO2 分压 (atm)
P_CO2 = 1.0  # 纯 CO2 下

# ----------------------------
# 碳酸体系电荷平衡方程
# ----------------------------
def carbon_system(pH, KH, P_CO2, K1, K2, Alk):
    H = 10**(-pH)
    CO2_aq = KH * P_CO2
    HCO3 = K1 * CO2_aq / H
    CO3 = K2 * HCO3 / H
    OH = 1e-14 / H
    # 电荷平衡：H+ + 其他阳离子 = OH- + HCO3- + 2CO3^2- + Alk
    charge_balance = H - (HCO3 + 2*CO3 + OH) + Alk
    return charge_balance

# ----------------------------
# 求解 pH，使用不同初值确保收敛
# ----------------------------
def solve_pH(Alk):
    if Alk < 1e-4:
        pH_guess = 4.0  # DI 水
    elif Alk < 4e-3:
        pH_guess = 5.0  # 自来水
    else:
        pH_guess = 7.0  # 废水
    pH_solution = fsolve(carbon_system, pH_guess, args=(KH, P_CO2, K1, K2, Alk))
    return pH_solution[0]

# ----------------------------
# 计算三种水 pH
# ----------------------------
pH_DI = solve_pH(alk_DI)
pH_tap = solve_pH(alk_tap)
pH_waste = solve_pH(alk_waste)

# ----------------------------
# 输出结果
# ----------------------------
print(f"DI 水平衡 pH: {pH_DI:.2f}")
print(f"自来水平衡 pH: {pH_tap:.2f}")
print(f"废水平衡 pH (假设碱度): {pH_waste:.2f}")

# ----------------------------
# 可视化柱状图
# ----------------------------
waters = ['DI 水', '自来水', '废水']
pH_values = [pH_DI, pH_tap, pH_waste]

plt.bar(waters, pH_values, color=['skyblue', 'lightgreen', 'salmon'])
plt.ylabel('平衡 pH')
plt.title('CO₂ 溶解后不同水体平衡 pH')
plt.ylim(0, 14)
plt.show()