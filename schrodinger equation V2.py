import numpy as np
import matplotlib.pyplot as plt

xmin = 0
xmax = 3
delta_x = 0.001

V_0 = 10

x_vals = np.arange(0, xmax+delta_x, delta_x)

V = np.zeros(len(x_vals)) - V_0
potential_step_index = np.where(x_vals > 1)[0][0]
V[potential_step_index:] = 0

E = -3
Energies = np.arange(-4.7, -4.6, 0.00001)
u_init = 1
psi_init = 0
psi = [psi_init]
u = [u_init]

fig = plt.figure()
ax = fig.gca()

for E in Energies:
    for i in range(len(x_vals)):
        psi_h = psi[i] + delta_x/2 * u[i]
        u_h = u[i] + delta_x/2 * -(E-V[i])*psi[i]

        new_psi = psi[i] + delta_x*u_h
        new_u = u[i] + delta_x* -(E-V[i])*psi_h

        psi.append(new_psi)
        u.append(new_u)

    psi = np.array(psi)
    psi_squared = psi**2
    constant = np.trapz(psi_squared[:-1], x_vals)
    psi = psi/np.sqrt(constant)
    
    if abs(psi[-1]) < 0.001:
        ax.plot(x_vals, psi[:-1])
        print("Energy value was found to be {:.4f}".format(E))

    psi = [psi_init]
    u = [u_init]

ax.plot(x_vals, V)
ax.set_ylim(-11, 1.5)
plt.show()
