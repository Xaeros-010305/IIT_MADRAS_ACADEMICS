# Setting up imports 
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
#=============================================================================================================================
# Question 1
# Setting up parameters
k1 = 1.0
k2= 0.5
k3 = 0.2

# Inital conditions [G0, Py0, L0, W0]
y0 = [1, 0, 0, 0]

def derivatives_1(t, y):
    G, Py, L, W = y

    dG = -((k1*G) + (k3*G**2))
    dPy = (k1*G) - (k2*Py)
    dL = (k2*Py)
    dW = (k3*G**2)

    return [dG, dPy, dL, dW]

sol_1 = solve_ivp(derivatives_1, [0, 10], y0)
G_1, Py_1, L_1, W_1 = sol_1.y
#=============================================================================================================================
plt.figure(figsize=(8,5))

plt.plot(sol_1.t, G_1, linewidth=2, label='Glucose (G)')

plt.xlabel('Time (min)', fontsize=12)
plt.ylabel('Concentration (mol/L)', fontsize=12)
plt.title('Glucose Concentration Profile over time', fontsize=14)

plt.grid(True, linestyle='--', alpha=0.6)

plt.legend(shadow=True, fontsize=10)

plt.tight_layout()
plt.show()
#=============================================================================================================================
plt.figure(figsize=(8,5))

plt.plot(sol_1.t, Py_1, linewidth=2, label='Pyruvate (Py)')

plt.xlabel('Time (min)', fontsize=12)
plt.ylabel('Concentration (mol/L)', fontsize=12)
plt.title('Pyruvate Concentration Profile', fontsize=14)

plt.grid(True, linestyle='--', alpha=0.6)

plt.legend(shadow=True, fontsize=10)

plt.tight_layout()
plt.show()
#=============================================================================================================================
plt.figure(figsize=(8,5))

plt.plot(sol_1.t, L_1, linewidth=2, label='Lactate (L)')

plt.xlabel('Time (min)', fontsize=12)
plt.ylabel('Concentration (mol/L)', fontsize=12)
plt.title('Lactate Concentration Profile', fontsize=14)

plt.grid(True, linestyle='--', alpha=0.6)

plt.legend(shadow=True, fontsize=10)

plt.tight_layout()
plt.show()
#=============================================================================================================================
plt.figure(figsize=(8,5))

plt.plot(sol_1.t, W_1, linewidth=2, label='Waste-product (W)')

plt.xlabel('Time (min)', fontsize=12)
plt.ylabel('Concentration (mol/L)', fontsize=12)
plt.title('Waste-product Concentration Profile', fontsize=14)

plt.grid(True, linestyle='--', alpha=0.6)

plt.legend(shadow=True, fontsize=10)

plt.tight_layout()
plt.show()
#=============================================================================================================================
# Part d
# Setting up parameters
k2= 0.5
k3 = 0.2

Vm, Km = 1.2, 0.5

# Inital conditions [G0, Py0, L0, W0]
y0 = [1, 0, 0, 0]

def derivatives_2(t, y):
    G, Py, L, W = y

    dG = -(((Vm*G)/(Km + G)) + (k3*G**2))
    dPy = ((Vm*G)/(Km + G)) - (k2*Py)
    dL = (k2*Py)
    dW = (k3*G**2)

    return [dG, dPy, dL, dW]

sol_2 = solve_ivp(derivatives_2, [0, 10], y0)
G_2, Py_2, L_2, W_2 = sol_2.y
#=============================================================================================================================
plt.figure(figsize=(8,5))

plt.plot(sol_1.t, G_1, linestyle='--', linewidth=2, label='Glucose (G, First-order)')
plt.plot(sol_2.t, G_2, linestyle='-', linewidth=2, label='Glucose (G, Michaelis-Menten)')

plt.xlabel('Time (min)', fontsize=12)
plt.ylabel('Concentration (mol/L)', fontsize=12)
plt.title('Glucose Concentration Comparison', fontsize=14)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(shadow=True)

plt.tight_layout()
plt.show()
#=============================================================================================================================
plt.figure(figsize=(8,5))

plt.plot(sol_1.t, Py_1, linestyle='--', linewidth=2, label='Pyruvate (Py, First-order)')
plt.plot(sol_2.t, Py_2, linestyle='-', linewidth=2, label='Pyruvate (Py, Michaelis-Menten)')

plt.xlabel('Time (min)', fontsize=12)
plt.ylabel('Concentration (mol/L)', fontsize=12)
plt.title('Pyruvate Concentration Comparison', fontsize=14)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(shadow=True)

plt.tight_layout()
plt.show()
#=============================================================================================================================
plt.figure(figsize=(8,5))

plt.plot(sol_1.t, L_1, linestyle='--', linewidth=2, label='Lactate (L, First-order)')
plt.plot(sol_2.t, L_2, linestyle='-', linewidth=2, label='Lactate (L, Michaelis-Menten)')

plt.xlabel('Time (min)', fontsize=12)
plt.ylabel('Concentration (mol/L)', fontsize=12)
plt.title('Lactate Concentration Comparison', fontsize=14)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(shadow=True)

plt.tight_layout()
plt.show()
#=============================================================================================================================
plt.figure(figsize=(8,5))

plt.plot(sol_1.t, W_1, linestyle='--', linewidth=2, label='Waste-product (W, First-order)')
plt.plot(sol_2.t, W_2, linestyle='-', linewidth=2, label='Waste-product (W, Michaelis-Menten)')

plt.xlabel('Time (min)', fontsize=12)
plt.ylabel('Concentration (mol/L)', fontsize=12)
plt.title('Waste-product Concentration Comparison', fontsize=14)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(shadow=True)

plt.tight_layout()
plt.show()
#=============================================================================================================================
# Part e
# Setting up parameters
k2= 0.5
k3 = 0.2

Vm = 1.2
Kms = [0.2, 1.0, 2.0]

# Inital conditions [G0, Py0, L0, W0]
y0 = [1, 0, 0, 0]

def derivatives_2(t, y):
    G, Py, L, W = y

    dG = -(((Vm*G)/(Km + G)) + (k3*G**2))
    dPy = ((Vm*G)/(Km + G)) - (k2*Py)
    dL = (k2*Py)
    dW = (k3*G**2)

    return [dG, dPy, dL, dW]

G_dict, Py_dict, L_dict, W_dict = {}, {}, {}, {}
Km = 0.5

sol = solve_ivp(derivatives_2, [0, 10], y0)
G, Py, L, W = sol.y

G_dict[Km] = (sol.t, G)
Py_dict[Km] = (sol.t, Py)
L_dict[Km] = (sol.t, L)
W_dict[Km] = (sol.t, W)

for Km in Kms:
    sol = solve_ivp(derivatives_2, [0, 10], y0)
    G, Py, L, W = sol.y
    
    G_dict[Km] = (sol.t, G)
    Py_dict[Km] = (sol.t, Py)
    L_dict[Km] = (sol.t, L)
    W_dict[Km] = (sol.t, W)
#=============================================================================================================================
plt.figure(figsize=(10,6))

for i in G_dict.keys():
    t, G = G_dict[i]
    plt.plot(t, G, linewidth=2, label=f'Km = {i}')

plt.xlabel('Time (min)', fontsize=12)
plt.ylabel('Concentration (mol/L)', fontsize=12)
plt.title('Glucose Concentration Profile', fontsize=14)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(shadow=True, fontsize=10)

plt.tight_layout()
plt.show()
# ---------------------------------------------------
plt.figure(figsize=(10,6))

for i in Py_dict.keys():
    t, Py = Py_dict[i]
    plt.plot(t, Py, linewidth=2, label=f'Km = {i}')

plt.xlabel('Time (min)', fontsize=12)
plt.ylabel('Concentration (mol/L)', fontsize=12)
plt.title('Pyruvate Concentration Profile', fontsize=14)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(shadow=True, fontsize=10)

plt.tight_layout()
plt.show()
# ---------------------------------------------------
plt.figure(figsize=(10,6))

for i in L_dict.keys():
    t, L = L_dict[i]
    plt.plot(t, L, linewidth=2, label=f'Km = {i}')

plt.xlabel('Time (min)', fontsize=12)
plt.ylabel('Concentration (mol/L)', fontsize=12)
plt.title('Lactate Concentration Profile', fontsize=14)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(shadow=True, fontsize=10)

plt.tight_layout()
plt.show()
# ---------------------------------------------------
plt.figure(figsize=(10,6))

for i in W_dict.keys():
    t, W = W_dict[i]
    plt.plot(t, W, linewidth=2, label=f'Km = {i}')

plt.xlabel('Time (min)', fontsize=12)
plt.ylabel('Concentration (mol/L)', fontsize=12)
plt.title('Waste-product Concentration Profile', fontsize=14)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(shadow=True, fontsize=10)

plt.tight_layout()
plt.show()
#=============================================================================================================================
# Question 2
# We use a 4-gene system:
# CycD, Rb, E2F, and CycE.

def update_CycD(CycD, Rb, E2F, CycE):
    return (CycD and not(CycE))

def update_Rb(CycD, Rb, E2F, CycE):
    return (not CycD)

def update_E2F(CycD, Rb, E2F, CycE):
    return (not Rb)  

def update_CycE(CycD, Rb, E2F, CycE):
    return (E2F)

def step(state):
    CycD, Rb, E2F, CycE = state
    
    new_CycD = update_CycD(CycD, Rb, E2F, CycE)
    new_Rb = update_Rb(CycD, Rb, E2F, CycE)
    new_E2F = update_E2F(CycD, Rb, E2F, CycE)
    new_CycE = update_CycE(CycD, Rb, E2F, CycE)
        
    int(new_CycD)
    
    return (int(new_CycD), int(new_Rb), int(new_E2F), int(new_CycE))

def simulate(initial_state, steps=10):
    state = initial_state
    trajectory = [state]
    
    for _ in range(steps):
        state = step(state)
        trajectory.append(state)
    
    return trajectory

traj = simulate((1, 0, 1, 0), steps=15)

for t, state in enumerate(traj):
    print(f"t={t}: {state}")

def plot_trajectory(traj):
    
    plt.figure(figsize=(10, 6))
    
    CycD = [s[0] for s in traj]
    Rb = [s[1] for s in traj]
    E2F = [s[2] for s in traj]
    CycE = [s[3] for s in traj]

    plt.plot(CycD, label="Cyclin D")
    plt.plot(Rb, label="Retinoblastoma protein")
    plt.plot(E2F, label="E2F transcription factor")
    plt.plot(CycE, label="Cyclin E")
    
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("State (0/1)")
    plt.title("Boolean Network Dynamics")
    plt.show()

plot_trajectory(traj)
#=============================================================================================================================
Initial_states = {}

Initial_states[1] = (0, 0, 0, 0)
Initial_states[2] = (0, 0, 0, 1)
Initial_states[3] = (0, 0, 1, 0)
Initial_states[4] = (0, 0, 1, 1)
Initial_states[5] = (0, 1, 0, 0)
Initial_states[6] = (0, 1, 0, 1)
Initial_states[7] = (0, 1, 1, 0)
Initial_states[8] = (0, 1, 1, 1)

Initial_states[9] = (1, 0, 0, 0)
Initial_states[10] = (1, 0, 0, 1)
Initial_states[11] = (1, 0, 1, 0)
Initial_states[12] = (1, 0, 1, 1)
Initial_states[13] = (1, 1, 0, 0)
Initial_states[14] = (1, 1, 0, 1)
Initial_states[15] = (1, 1, 1, 0)
Initial_states[16] = (1, 1, 1, 1)
#=============================================================================================================================
for key in Initial_states:
    
    traj = simulate(Initial_states[key], steps=10)
    
    print(f"Initial State: {Initial_states[key]}")
    for t, state in enumerate(traj):
        print(f"t={t}: {state}")
    
    plt.figure(figsize=(10, 6))
    
    CycD = [s[0] for s in traj]
    Rb = [s[1] for s in traj]
    E2F = [s[2] for s in traj]
    CycE = [s[3] for s in traj]

    plt.plot(CycD, label="Cyclin D")
    plt.plot(Rb, label="Retinoblastoma protein")
    plt.plot(E2F, label="E2F transcription factor")
    plt.plot(CycE, label="Cyclin E")
    
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("State (0/1)")
    plt.title(f"Boolean Network Dynamics (Initial [CycD, Rb, E2F, CycE] = {Initial_states[key]})")
    plt.show()
    print("="*125)
#=============================================================================================================================