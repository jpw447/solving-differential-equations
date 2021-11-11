import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Initial values
J0 = 1
v0 = 0

# Values for v and J: v = dJ/dz
v_values = [v0]
J_values = [J0]

# Setting up values for z
z_max = 100
delta_z = 0.01
z_values = np.arange(0.0001, z_max+delta_z, delta_z)

for i in range(len(z_values)):
    # J and v evaluated at z
    z = z_values[i]
    v = v_values[i]
    J = J_values[i]

    # Coupled DEs to get derivatives from analytical equations
    J_prime = v
    v_prime = -(1/z)*v-J

    # Values at half step
    Jh = J + (delta_z/2)*J_prime # halfway J
    vh = v + (delta_z/2)*v_prime

    # Using halfway values to get new values at full step
    new_J = J + delta_z*vh
    new_v = v + delta_z*v_prime

    J_values.append(new_J)
    v_values.append(new_v)

del J_values[-1] # Last value corresponds to a z value not in the z array, so delete it
fig_numerical = plt.figure(figsize=(8,10))
ax_numerical = fig_numerical.gca()
ax_numerical.plot(z_values, J_values)

# Function used by odeint to find the derivative, given an initial conditon
def dU_dz(U, z):
    '''
    Function of interest is given as U, which contains U[0]=J, U[1]=J'.
    Function then returns J' and J'' for use in odeint.
    '''
    v = U[1]                    # J'
    v_prime = -(1/z)*U[1]-U[0]  # J''
    return v, v_prime
    
U = [J0, v0]                   # Initial values with which dU_dz calculates the derivative
sol = odeint(dU_dz, U, z_values) # odeint now uses dU_dz to calculate the derivatives for every step of the z array
J = sol[:, 0]
ax_numerical.plot(z_values, J, 'r--')

J_error = J - J_values          # Error in J compared to real values

# Figure for axis containing errors in numerical solution
fig_error = plt.figure(figsize=(8, 10))
ax_error = fig_error.gca()
ax_error.plot(z_values, J_error)
ax_error.hlines(0, 0, z_max, linestyle='dashed', color='red')

plt.show()