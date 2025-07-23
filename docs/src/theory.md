```@meta
CurrentModule = PSDTFE
```

# Theory

Inline equation: $E = m c^2$

Display equation: $$v = H d$$


$$\begin{equation}
e^{i \pi} = - 1
\end{equation}$$

$$a_0 x^n + a_1 x^{n-1} + \dots + a_{n-1} x + a_n = 0,$$

```@example
using Plots

x = range(0, 10, length=100)
y = sin.(x)
plot(x, y)
```

Here's an equation:

```math
\frac{n!}{k!(n - k)!} = \binom{n}{k}
```

This is the binomial coefficient.

---

To write a system of equations, use the `aligned` environment:

```math
\begin{aligned}
\nabla\cdot\mathbf{E}  &= 4 \pi \rho \\
\nabla\cdot\mathbf{B}  &= 0 \\
\nabla\times\mathbf{E} &= - \frac{1}{c} \frac{\partial\mathbf{B}}{\partial t} \\
\nabla\times\mathbf{B} &= - \frac{1}{c} \left(4 \pi \mathbf{J} + \frac{\partial\mathbf{E}}{\partial t} \right)
\end{aligned}
```

These are Maxwell's equations.

## Lagrangian fluid dynamics

The evolution of our universe can be described in Lagrangian fluid dynamics in terms of the Lagrangian map $x_t(q) = q + s_t(q)$ mapping a point in the space of initial conditions (Lagrangian space) to a point in the current universe (Eulerian space). The displacement field $s_t(q)$ captures the displacement of a particle starting at $q$ in time $t$. Given the Lagrangian map $x_t$, the density field is given by 

```math
\rho(x') = \sum_{q \in x_t^{-1}(x')} \frac{\rho_u}{\| \nabla x_t(q)\|}
```

The Phase-Space DTFE method implements this density field, and the associated velocity fields for cosmological $N$-body simulations using both Delaunay tessellations and phase-space methods.

### The Delaunay Tessellation Field Estimator (DTFE)
The DTFE method converts the current positions of the $N$-body particles and their velocities in to an estimate of the density and velocity fields in Eulerian space.

### The Phase-Space Delaunay Tessellation Field Estimator (PS-DTFE)
The PS-DTFE method extends the DTFE method to phase-space.