
import numpy as np
import matplotlib.pyplot as plt

def schrodinger_solver(x_vals, delta_x, psi, u):
    for i in range(len(x_vals)):
        # Halfway values
        psi_h = psi[i] + delta_x/2 * u[i]
        u_h = u[i] + delta_x/2 * -(E-V[i])*psi[i]

        # New values at x + delta x
        new_psi = psi[i] + delta_x*u_h
        new_u = u[i] + delta_x* -(E-V[i])*psi_h

        psi.append(new_psi)
        u.append(new_u)

    # Normalising the wave function
    psi = np.array(psi)
    psi_squared = psi**2 # Psi is real here
    constant = np.trapz(psi_squared[:-1], x_vals)
    psi = psi/np.sqrt(constant)

    return psi

# Setting up x values from 0 to some maximum
x_min = 0
x_max = 4
delta_x = 0.001
x_vals = np.arange(x_min, x_max+delta_x, delta_x)

# Establishing potential
V_0 = 10
V = np.zeros(len(x_vals)) - V_0
potential_step_index = np.where(x_vals > 1)[0][0]
V[potential_step_index:] = 0
V[0] = 100000

# Establishing energies to check
E_min = -5
E_max = -3
E_step = 0.0001
Energies = np.arange(E_min, E_max+E_step, E_step)

# Initial conditions of psi and u (u=dpsi/dx)
u_init = 1
psi_init = 0
psi = [psi_init]
u = [u_init]

# Tolerance within which to check if psi has settled on 0
tolerance = 0.01
found = False # Whether the ground state was found

# Creating figure
fig = plt.figure(figsize=(10, 8))
ax = fig.gca()
ax.set_ylim(-11, 1.5)
ax.set_xlabel("$x$", fontsize=16)
ax.set_ylabel("$\psi$ and $V$ (J)", fontsize=16)
ax.set_title("Ground State Wave Function for a Half-Finite Potential Well", fontsize=20)
ax.plot(x_vals, V, 'k', label="Potential $V$")

# Checks each energy and calculates wave function until ground state energy is found
for E in Energies:
    psi = schrodinger_solver(x_vals, delta_x, psi, u)
    
    # Checking if psi has tended towards 0. Stops further calculations if so
    if abs(psi[-1]) < tolerance:
        # Modulus squared of psi plotted at an arbitrary point within the potential.
        psi_squared = psi**2
        ax.plot(x_vals, V, 'k')
        ax.plot(x_vals, psi_squared[:-1]-6, 'r', label="Ground state $\psi$")
        ax.set_xlabel("$x$", fontsize=16)
        ax.set_ylabel("$V$ (J)", fontsize=16)
        ax.set_title("Ground State Wave Function for a Finite Potential Well", fontsize=20)
        ax.legend(fontsize=16)
        print("Energy value was found to be {:.4f}".format(E))
        found = True
        break
    
    # Resets lists for use in next loop
    psi = [psi_init]
    u = [u_init]

plt.show()

# Plots potential and the wave function found for ground state energy
if found == False:
    print("No energy eigenvalue was found within the tolerance.")
