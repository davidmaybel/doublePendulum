mport numpy as np
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
