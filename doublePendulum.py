import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button, TextBox
from scipy.integrate import solve_ivp

# --- Pendulum Parameters ---
L1, L2 = 1.0, 0.8
m1, m2 = 1.0, 0.8
g = 9.81

# --- Initial Conditions ---
theta1_0, theta2_0 = np.pi/2, np.pi/2
omega1_0, omega2_0 = 0.0, 0.0

# --- Time Parameters ---
duration, fps = 10, 30
dt = 1.0 / fps
t_eval = np.arange(0, duration, dt)

# --- Solving Differential Equations ---
def get_solution(g_val):
    def derivatives(t, y):
        theta1, omega1, theta2, omega2 = y
        delta = theta2 - theta1
        den1 = (m1 + m2)*L1 - m2*L1*np.cos(delta)**2
        den2 = (L2 / L1) * den1

        dtheta1 = omega1
        domega1 = (m2*L1*omega1**2*np.sin(delta)*np.cos(delta) +
                   m2*g_val*np.sin(theta2)*np.cos(delta) +
                   m2*L2*omega2**2*np.sin(delta) -
                   (m1 + m2)*g_val*np.sin(theta1)) / den1

        dtheta2 = omega2
        domega2 = (-m2*L2*omega2**2*np.sin(delta)*np.cos(delta) +
                   (m1 + m2)*g_val*np.sin(theta1)*np.cos(delta) -
                   (m1 + m2)*L1*omega1**2*np.sin(delta) -
                   (m1 + m2)*g_val*np.sin(theta2)) / den2

        return [dtheta1, domega1, dtheta2, domega2]

    sol = solve_ivp(derivatives, (0, duration),
                    [theta1_0, omega1_0, theta2_0, omega2_0],
                    t_eval=t_eval, rtol=1e-6, atol=1e-6)
    return sol.y[0], sol.y[2]

# Initial solution
theta1, theta2 = get_solution(g)
x1 = L1 * np.sin(theta1)
y1 = -L1 * np.cos(theta1)
x2 = x1 + L2 * np.sin(theta2)
y2 = y1 - L2 * np.cos(theta2)

# --- Plot Setup ---
fig, ax = plt.subplots(figsize=(8, 6))
plt.subplots_adjust(bottom=0.25)
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_aspect('equal')
line1, = ax.plot([], [], 'b-', lw=2)
line2, = ax.plot([], [], 'r-', lw=2)
bob1, = ax.plot([], [], 'bo', ms=10)
bob2, = ax.plot([], [], 'ro', ms=10)
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

# --- Animation Controls ---
is_running = [True]  # use list for mutability in closure

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    bob1.set_data([], [])
    bob2.set_data([], [])
    time_text.set_text('')
    return line1, line2, bob1, bob2, time_text

def animate(i):
    if is_running[0]:
        line1.set_data([0, x1[i]], [0, y1[i]])
        line2.set_data([x1[i], x2[i]], [y1[i], y2[i]])
        bob1.set_data([x1[i]], [y1[i]])
        bob2.set_data([x2[i]], [y2[i]])
        time_text.set_text(f'Time = {t_eval[i]:.2f}s')
    return line1, line2, bob1, bob2, time_text

ani = FuncAnimation(fig, animate, frames=len(t_eval),
                    init_func=init, blit=True, interval=1000/fps)

# --- Pause Button ---
pause_ax = plt.axes([0.7, 0.05, 0.1, 0.075])
pause_button = Button(pause_ax, 'Pause/Start')

def toggle_run(event):
    is_running[0] = not is_running[0]

pause_button.on_clicked(toggle_run)

# --- Gravity TextBox ---
text_ax = plt.axes([0.1, 0.05, 0.15, 0.075])
text_box = TextBox(text_ax, 'Gravity g:', initial=str(g))

def update_gravity(text):
    global theta1, theta2, x1, x2, y1, y2
    try:
        g_val = float(text)
        theta1, theta2 = get_solution(g_val)
        x1 = L1 * np.sin(theta1)
        y1 = -L1 * np.cos(theta1)
        x2 = x1 + L2 * np.sin(theta2)
        y2 = y1 - L2 * np.cos(theta2)
    except ValueError:
        print("Invalid gravity value.")

text_box.on_submit(update_gravity)

plt.title("Interactive Double Pendulum")
plt.show()
