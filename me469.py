import numpy as np
import pandas as pd
import math
from sympy import solve, Function, symbols, cos, sin, diff, simplify, solveset, Eq, S

# Defining time variable
t = symbols('t')

# Defining parameter symbols

# limb lengths
l_1, l_2, l_3, l_4 = symbols('l_1 l_2 l_3 l_4')

# limb masses
m_1, m_2, m_3 = symbols('m_1, m_2, m_3')

# joint stiffness's
k_B, k_C = symbols('k_B, k_C')

# joint damping
b_B, b_C = symbols('b_B, b_C')

# joint neutral angles
psi_B, psi_C = symbols('psi_B, psi_C')

# foot angle
phi_D = symbols('phi_E')

# Walking speed
v = symbols('v')

# Defining coordinate symbols

# limb orientations
theta_1 = Function('theta_1')(t)
theta_2 = Function('theta_2')(t)
theta_3 = Function('theta_3')(t)

# hip height
y = Function('y')(t)

# Calculate positions

# rates of change
omega_1 = diff(theta_1)
omega_2 = diff(theta_2)
omega_3 = diff(theta_3)
v = diff(y)

# heel angle
theta_4 = theta_3 - phi_D
omega_4 = diff(theta_4, t)

# joint positions
x_A = 0
y_A = y
x_B = x_A + l_1 * sin(theta_1)
y_B = y_A - l_1 * cos(theta_1)
x_C = x_B + l_2 * sin(theta_2)
y_C = y_B - l_2 * cos(theta_2)
x_D = x_C + l_3 * sin(theta_3)
y_D = y_C - l_3 * cos(theta_3)
x_E = x_C + l_4 * sin(theta_4)
y_E = y_C - l_4 * cos(theta_4)

# joint velocities
u_A = diff(x_A, t)
v_A = diff(y_A, t)
u_B = diff(x_B, t)
v_B = diff(y_B, t)
u_C = diff(x_C, t)
v_C = diff(y_C, t)
u_D = diff(x_D, t)
v_D = diff(y_D, t)
u_E = diff(x_E, t)
v_E = diff(y_E, t)

# Centers of masses
x_1 = x_A + l_1 * sin(theta_1) / 2
y_1 = y_A - l_1 * cos(theta_1) / 2
x_2 = x_B + l_2 * sin(theta_2) / 2
y_2 = y_B - l_2 * cos(theta_2) / 2
x_3 = x_C + l_3 * sin(theta_3) / 2
y_3 = y_C - l_3 * cos(theta_3) / 2

# COM velocities
u_1 = diff(x_1, t)
v_1 = diff(y_1, t)
u_2 = diff(x_2, t)
v_2 = diff(y_2, t)
u_3 = diff(x_3, t)
v_3 = diff(y_3, t)

# Moments of inertia
I_1 = m_1 * l_1 ** 2 / 12
I_2 = m_2 * l_2 ** 2 / 12
I_3 = m_3 * l_3 ** 2 / 12

# read the expression for x_E
x_E

# Kinetic energy of the system
T = m_1 * (u_1 ** 2 + v_1 ** 2) / 2 + m_2 * (u_2 ** 2 + v_2 ** 2) / 2 + m_3 * (
        u_3 ** 2 + v_3 ** 2) / 2 + I_1 * omega_1 ** 2 / 2 + I_2 * omega_2 ** 2 / 2 + I_3 * omega_3 ** 2 / 2

# Potential energy of the system
V = k_B * (theta_1 - theta_2 - psi_B) ** 2 / 2 + k_C * (theta_2 - theta_3 - psi_C) ** 2 / 2

# Content energy of the system
D = b_B * (omega_1 - omega_2) ** 2 / 2 + b_C * (omega_2 - omega_3) ** 2 / 2

# Constraints

# Phase A - Swing
y_A_SW = y
const_FP_y_A = y_A - y_A_SW

# Phase B - Heel Strike
const_HS_y_E = y_E
const_HS_u_E = u_E + v

# Phase E - Foot Planted
const_FP_y_E = y_E
const_FP_y_D = y_D
const_FP_u_E = u_E + v
const_FP_u_D = u_D + v

# Phase D - Toe Push-Off
const_TLO_y_D = y_D
const_TLO_u_D = u_D + v

'''
# Equations of Motion

# y
EOM_1 = simplify(diff(diff(T, v_A), t) - diff(T, y) + diff(D, v_A) + diff(V, y))
print(EOM_1)

# theta_1
EOM_2 = simplify(diff(diff(T, omega_1), t) - diff(T, theta_1) + diff(D, omega_1) + diff(V, theta_1))
print(EOM_2)

# theta_2
EOM_3 = simplify(diff(diff(T, omega_2), t) - diff(T, theta_2) + diff(D, omega_2) + diff(V, theta_2))
print(EOM_3)

# theta_3
EOM_4 = simplify(diff(diff(T, omega_3), t) - diff(T, theta_3) + diff(D, omega_3) + diff(V, theta_3))
print(EOM_4)
'''

# Defining Lagrange's Equation

L1_t1 = solveset(Eq(const_HS_y_E, y_E), y)
print('Lagrange expression in terms of y')
print(L1_t1)

L1_t2 = solveset(Eq(const_HS_y_E, y_E), theta_2)
print('Lagrange expression in terms of theta_2')
print(L1_t2)

L1_t3 = solveset(Eq(const_HS_y_E, y_E), theta_3)
print('Lagrange expression in terms of theta_3')
print(L1_t3)
