```@meta
CurrentModule = PSDTFE
```

# Theory

The evolution of our universe can be described in Lagrangian fluid dynamics in terms of the Lagrangian map $x_t(q) = q + s_t(q)$ mapping a point in the space of initial conditions (Lagrangian space) to a point in the current universe (Eulerian space). The displacement field $s_t(q)$ captures the displacement of a particle starting at $q$ in time $t$. Given the Lagrangian map $x_t$, the density field is given by 

```math
\rho(x') = \sum_{q \in x_t^{-1}(x')} \frac{\rho_u}{\| \nabla x_t(q)\|}
```

The Phase-Space DTFE method implements this density field, and the associated velocity fields for cosmological $N$-body simulations using both Delaunay tessellations and phase-space methods.

## The Delaunay Tessellation Field Estimator (DTFE)
The DTFE method converts the current positions of the $N$-body particles and their velocities in to an estimate of the density and velocity fields in Eulerian space.

## The Phase-Space Delaunay Tessellation Field Estimator (PS-DTFE)
The PS-DTFE method extends the DTFE method to phase-space.


```@example
using Plots

x = range(0, 10, length=100)
y = sin.(x)
plot(x, y)
```

```@example
A = 5
```

```@example
B = A + 1
```