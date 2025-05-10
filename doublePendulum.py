import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_ivp
import os

# Create directory on desktop if it doesn't exist
desktop_path = os.path.join(os.path.expanduser('~'), 'Desktop', 'doublePendulum')
os.makedirs(desktop_path, exist_ok=True)

# Pendulum parameters
L1 = 1.0      # Length of first rod (m)
L2 = 0.8      # Length of second rod (m)
m1 = 1.0      # Mass of first bob (kg)
m2 = 0.8      # Mass of second bob (kg)
g = 9.81      # Acceleration due to gravity (m/s²)

# Initial conditions [theta1, omega1, theta2, omega2]
theta1_0 = np.pi/2    # Initial angle (rad)
theta2_0 = np.pi/2    # Initial angle (rad)
omega1_0 = 0.0        # Initial angular velocity (rad/s)
omega2_0 = 0.0        # Initial angular velocity (rad/s)

# Time parameters
duration = 10.0       # Duration of simulation (s)
fps = 30              # Frames per second
dt = 1.0/fps          # Time step
t_eval = np.arange(0, duration, dt)

def derivatives(t, y):
    """Calculate derivatives for the double pendulum system"""
    theta1, omega1, theta2, omega2 = y
    
    delta = theta2 - theta1
    den1 = (m1 + m2)*L1 - m2*L1*np.cos(delta)*np.cos(delta)
    den2 = (L2/L1)*den1
    
    dtheta1 = omega1
    domega1 = (m2*L1*omega1*omega1*np.sin(delta)*np.cos(delta) +
               m2*g*np.sin(theta2)*np.cos(delta) +
               m2*L2*omega2*omega2*np.sin(delta) -
               (m1 + m2)*g*np.sin(theta1)) / den1
    
    dtheta2 = omega2
    domega2 = (-m2*L2*omega2*omega2*np.sin(delta)*np.cos(delta) +
               (m1 + m2)*g*np.sin(theta1)*np.cos(delta) -
               (m1 + m2)*L1*omega1*omega1*np.sin(delta) -
               (m1 + m2)*g*np.sin(theta2)) / den2
    
    return [dtheta1, domega1, dtheta2, domega2]

# Solve the differential equations
sol = solve_ivp(derivatives, (0, duration),
                [theta1_0, omega1_0, theta2_0, omega2_0],
                t_eval=t_eval, rtol=1e-6, atol=1e-6)

# Extract the solution
theta1 = sol.y[0]
theta2 = sol.y[2]

# Convert angles to Cartesian coordinates
x1 = L1 * np.sin(theta1)
y1 = -L1 * np.cos(theta1)
x2 = x1 + L2 * np.sin(theta2)
y2 = y1 - L2 * np.cos(theta2)

# Set up the figure and axis
fig, ax = plt.subplots(figsize=(10, 8))
ax.set_xlim(-(L1 + L2 + 0.1), L1 + L2 + 0.1)
ax.set_ylim(-(L1 + L2 + 0.1), L1 + L2 + 0.1)
ax.set_aspect('equal')
ax.grid()

# Create pendulum elements
line1, = ax.plot([], [], 'b-', lw=2)
line2, = ax.plot([], [], 'r-', lw=2)
bob1, = ax.plot([], [], 'bo', ms=10*m1)
bob2, = ax.plot([], [], 'ro', ms=10*m2)
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

def init():
    """Initialize the animation"""
    line1.set_data([], [])
    line2.set_data([], [])
    bob1.set_data([], [])
    bob2.set_data([], [])
    time_text.set_text('')
    return line1, line2, bob1, bob2, time_text

def animate(i):
    """Update the animation frame"""
    # Update pendulum rods
    line1.set_data([0, x1[i]], [0, y1[i]])
    line2.set_data([x1[i], x2[i]], [y1[i], y2[i]])
    
    # Update pendulum bobs (as arrays)
    bob1.set_data([x1[i]], [y1[i]])
    bob2.set_data([x2[i]], [y2[i]])
    
    # Update time display
    time_text.set_text(f'Time = {t_eval[i]:.2f}s')
    
    return line1, line2, bob1, bob2, time_text

# Create the animation
ani = FuncAnimation(fig, animate, frames=len(t_eval),
                    init_func=init, blit=True, interval=1000/fps)

plt.title('Double Pendulum Animation')
plt.xlabel('x position (m)')
plt.ylabel('y position (m)')

# Save the animation as GIF
output_path = os.path.join(desktop_path, 'double_pendulum.gif')
ani.save(output_path, writer='pillow', fps=fps, dpi=100)

print(f"Animation saved to: {output_path}")
plt.close()
