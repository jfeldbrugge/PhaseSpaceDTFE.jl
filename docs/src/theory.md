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