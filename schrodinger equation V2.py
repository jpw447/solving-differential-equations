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

E_upper = -4
E_lower = -6
u_init = 1
psi_init = 0
psi = [psi_init]
u = [u_init]

searching = True
below = True

fig = plt.figure()
ax = fig.gca()

ax.plot(x_vals, V)
counter = 0

while searching == True:
    # Checks whether previous iteration was below or above 0
    if below == True:
        E_upper = np.mean([E_upper, E_lower])
        E = E_lower
        print("Investigating energy {:.3f}".format(E))
    else:
        E_lower = np.mean([E_upper, E_lower])
        E = E_upper
        print("Investigating energy {:.3f}".format(E))

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

    # Checks whether or not psi is above or below
    if psi[-1] < 0: 
        below = True    
    elif psi[-1] > 0:
        below = False
    else: 
        # In case the final value of psi is 0
        searching = False
    
    # Resets lists for another loop
    ax.plot(x_vals, psi[:-1])
    psi = [psi_init]
    u = [u_init]
    counter += 1
    if counter == 10:
        searching = False
    plt.pause(0.5)



E = 4.5
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

ax.plot(x_vals, psi[:-1], 'k')
print("Energy value was found to be {:.4f}".format(E))



ax.plot(x_vals, V)
ax.set_ylim(-11, 1.5)
plt.show()
