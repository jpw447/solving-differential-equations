import numpy as np 
import matplotlib.pyplot as plt
import time

# Defining angles to test
N = 10000                  # Number of iterations to perform
theta_init = 0
theta_max = np.pi/8
theta_step = theta_max/N
theta_values = np.arange(theta_init, theta_max + theta_step, theta_step)

# Target, tolerance and constants
target = 3000 
tolerance = 0.5
k = 10e-4
g = -9.8

# Initial conditions and target
V = 2000

x_list = [0]
y_list = [0]

x_array_of_lists = np.ndarray(N, dtype=object)
y_array_of_lists = np.ndarray(N, dtype=object)
# Used to check if trajectory hit the target
check_list = np.ndarray(N, dtype=bool)*False

# Time step and start. Create array after loop
delta_t = 0.001
t = 0
start = time.time()
for i in range(0, N):
    theta = theta_values[i]
    # Updating components of V
    Vx = V*np.cos(theta)
    Vy = V*np.sin(theta)

    while (y_list[-1] >= 0):
        V_mag = np.sqrt(Vx**2 + Vy**2)
        x = x_list[-1]
        y = y_list[-1]

        # Calculations for Vx and x values
        x_prime = Vx
        Vx_prime = -k*V_mag*Vx

        xh = x + (delta_t/2)*Vx
        Vxh = Vx + (delta_t/2)*Vx_prime

        x = x + delta_t*Vxh
        Vx = Vx + delta_t*Vx_prime
        x_list.append(x)
        
        y_prime = Vy
        Vy_prime = g - k*V_mag*Vy

        yh = y + (delta_t/2)*Vy
        Vyh = Vy + (delta_t/2)*Vy_prime
        
        y = y + delta_t*Vyh
        Vy = Vy + delta_t*Vy_prime
        y_list.append(y)

        t += delta_t
    
    x_array_of_lists[i] = x_list
    y_array_of_lists[i] = y_list

    if (x_list[-1] < target+tolerance) and (x_list[-1] > target-tolerance):
        print("Final x value is "+str(x_list[-1]))
        check_list[i] = True
        break

    x_list = [0]
    y_list = [0]
    
end = time.time()
duration = np.round(end-start, 3)
print("That took "+str(duration)+" seconds")
fig_displacement = plt.figure(figsize=(10, 8))
ax_displacement = fig_displacement.gca()
ax_displacement.set_xlabel("$x$ (m)", fontsize=14)
ax_displacement.set_ylabel("$y$ (m)", fontsize=14)
ax_displacement.set_title("$x$ versus $y$ positions over time", fontsize=18)

ax_displacement.hlines(0, target-tolerance, target+tolerance, 'k')

# Checks where the loop stopped
j = np.where(check_list==True)[0][0] + 1

for i in range(j):
    x = x_array_of_lists[i]
    y = y_array_of_lists[i]
    if check_list[i] == True:
        ax_displacement.plot(x, y, 'g')
    else:
        ax_displacement.plot(x, y, 'r--')
plt.show()
print(check_list)