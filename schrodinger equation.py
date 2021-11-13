import numpy as np
import matplotlib.pyplot as plt

a = 1
V_infinity = 10  # Functional value for infinite potential

# Defining x values
delta_x = 0.01
x_min = -5
x_max = 15
x_values = np.arange(x_min, x_max+delta_x, delta_x)

# Indices where potential becomes infinite, finite, and zero
potential_wall_index = np.where(x_values > 0)[0][0]
potential_zero_index = np.where(x_values > a)[0][0]
# Building potential form
V = np.zeros(len(x_values)) - 1
V[:potential_wall_index] = V_infinity
V[potential_zero_index:] = 0

E_min = -10
E_max = 10
delta_E = 0.01
E_values = np.arange(E_min, E_max+delta_E, delta_E)

# Initial values for psi and dpsi/dx (u)
E = 5
psi = 0
u = 1

# Lists for plotting later
psi_list = [psi]
x = x_values[potential_wall_index:]

for i in range(len(x)):
    # Derivatives
    psi_prime = u
    u_prime = ((V[i]-E)/(a**2)) * psi

    # Halfway values
    psih = psi + (delta_x/2) * u
    uh = u + (delta_x/2) * u_prime

    # Re-evaluating u_prime at psih
    u_prime = ((V[i]-E)/(a**2)) * psih

    # New values
    psi = psi + delta_x * uh
    u = u + delta_x * u_prime
    
    psi_list.append(psi)
# Delting final value
del psi_list[-1]

constant = 

fig_wavefunction = plt.figure(figsize=(10,8))
ax_wavefunction = fig_wavefunction.gca()
ax_wavefunction.plot(x, psi_list, 'r--')
ax_wavefunction.plot(x_values, V, 'k')
ax_wavefunction.set_xlim(-1, 10)
ax_wavefunction.set_ylim(-2, 2)
plt.show()