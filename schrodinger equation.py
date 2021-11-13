import numpy as np
import matplotlib.pyplot as plt

a = 1
V_infinity = 10  # Functional value for infinite potential

# Defining x values
delta_x = 0.01
x_min = -5
x_max = 3
x_values = np.arange(x_min, x_max+delta_x, delta_x)

# Indices where potential becomes infinite, finite, and zero
potential_wall_index = np.where(x_values > 0)[0][0]
potential_zero_index = np.where(x_values > a)[0][0]
# Building potential form
V = np.zeros(len(x_values)) - 1
V[:potential_wall_index] = V_infinity
V[potential_zero_index:] = 0

E_min = -5
E_max = 0
N = 1000
delta_E = (abs(E_max)+abs(E_min))/N
E_values = np.arange(E_min, E_max+delta_E, delta_E)

psi_array_of_lists = np.ndarray(N, dtype=object)

# Initial values for psi and dpsi/dx (u)
E = 1
psi_init = 0
psi = psi_init
u = 1

# Lists for plotting later
psi_list = [psi]
x = x_values[potential_wall_index:]

for j in range(0, N):
    E = E_values[j]
    for i in range(len(x)):

        psih = psi + (delta_x/2) * u
        uh = u + (delta_x/2) * -((E-V[i])/(a**2)) * psi
        # New values
        psi = psi + delta_x * uh
        u = u + delta_x * -((E-V[i])/(a**2)) * psih
        
        psi_list.append(psi)
    # Delting final value   
    del psi_list[-1]
    
    psi_list = np.array(psi_list)
    psi_squared = psi_list**2
    normalisation = np.trapz(psi_squared, x)
    psi_list = psi_list/np.sqrt(normalisation)
    
    psi_array_of_lists[j] = psi_list
    psi_list = [psi_init]

fig_wavefunction = plt.figure(figsize=(10,8))
ax_wavefunction = fig_wavefunction.gca()
ax_wavefunction.plot(x_values, V, 'k')
ax_wavefunction.set_xlim(-1, 10)
ax_wavefunction.set_ylim(-2, 2)

ax_wavefunction.plot(x, psi_array_of_lists[900])

plt.show()

#%%
fig_wavefunction = plt.figure(figsize=(10,8))
ax_wavefunction = fig_wavefunction.gca()
ax_wavefunction.plot(x, psi_list, 'r--')
ax_wavefunction.plot(x_values, V, 'k')
ax_wavefunction.set_xlim(-1, x_max)
ax_wavefunction.set_ylim(-2, 2)
plt.show()