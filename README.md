# doublePendulum
This project simulates and visualizes the chaotic motion of a double pendulum system using numerical methods. It demonstrates key characteristics of nonlinear dynamic systems such as sensitivity to initial conditions and energy conservation. The entire solution is written in Python using open-source scientific libraries.

---

## Abstract

This project numerically simulates and visualizes the chaotic motion of a double pendulum in a computational physics setting. The system's nonlinear equations of motion are solved using a high-precision Runge-Kutta method, and the resulting motion is animated to illustrate the unpredictable behavior of chaotic systems.

---

## Methodology

### Physics Modeling
- Derived equations of motion using Lagrangian mechanics
- Modeled as coupled nonlinear ODEs with:
  - Customizable rod lengths (`L1`, `L2`)
  - Masses (`m1`, `m2`)
  - Gravity constant `g`

### Numerical Solution
- Used `scipy.integrate.solve_ivp` with `RK45` solver
- Fine-tuned accuracy with:
  - `rtol = 1e-6`
  - `atol = 1e-6`
- Simulated 10 seconds at 30 FPS

### Visualization
- Real-time animation via Matplotlib
- Converted polar coordinates to Cartesian
- Added time display and scaled bobs by mass

### Technical Features
- Modular, readable code structure
- Parameters easily adjustable (mass, gravity, angles)
- Output: `.gif` saved directly to desktop

---

## Key Results
- Chaotic behavior clearly demonstrated from simple initial conditions
- Verified energy conservation within numerical precision
- Generated professional-quality animation
- System exhibits hallmark of chaos: sensitive dependence on initial conditions

---

