import numpy as np
import matplotlib.pyplot as plt

# Initial conditions and target. Variables are updated in loop
V_init = 200
V = V_init
Vx = V_init*np.cos(theta_init)
Vy = V_init.np.sin(theta_init)

x = 0
y = 0
x_list = []
y_list = []

N = 20                  # Number of iterations to perform
theta_step = np.pi/N
theta_init = 0
theta_max = np.pi/2

target = 3000   # Target distance

theta = np.arange(theta_init, theta_max + theta_step, theta_step)
# Time step and start. Create array after loop
delta_t = 0.01
t = 0

for i in range(len(theta)):
    while (y > 0):
        V = np.sqrt(Vx**2 + Vy**2)

        # Calculations for Vx and x values
        x_prime = Vx
        Vx_prime = -k*V*Vx

        xh = x + (delta_t/2)*Vx
        Vxh = Vx + (delta_t/2)*Vx_prime

        x = x + delta_t*Vxh
        Vx = Vx + delta_t*Vx_prime
        x_list.append(x)
        
        y_prime = Vy
        Vy_prime = g - k*V*Vy

        yh = y + (delta_t/2)*Vy
        Vyh = Vy + (delta_t/2)*Vy_prime
        
        y = y + delta_t*Vyh
        Vy = Vy + delta_t*Vy_prime
        y_list.append(y)

        t += delta_t
