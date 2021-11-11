import numpy as np
import matplotlib.pyplot as plt

# Initial values
J0 = 1
v0 = 0

# Values for v and J: v = dJ/dz
v_values = [v0]
J_values = [J0]

# Setting up values for z
z_max = 10
delta_z = 0.01
z_values = np.arange(0.0001, z_max+delta_z, delta_z)

for i in range(len(z_values)):
    # J and v evaluated at z
    z = z_values[i]
    v = v_values[i]
    J = J_values[i]

    # Coupled DEs to get derivatives
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

del J_values[-1] # Last value corresponds to a z value not in the z array
plt.plot(z_values, J_values)
plt.show()