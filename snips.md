

$$
\begin{align}
  \min \quad f(x) \newline
  \textrm{s.t.} \quad g_{i}(x, \omega), \quad x \in \mathbb{R}, \quad \omega \in \Omega, \quad |\Omega| = \infty
\end{align}
$$


$$
\begin{align}
&\min_{w,b,\xi} \quad \frac{1}{2}w^{t}w+C\sum_{i=1}^{N}{\xi_{i}}\\\\
&\textrm{subject to} \quad y_{i}(w\phi(x_{i}+b))+\xi_{i}-1\\\\
\end{align}
$$


El problema _minimax_ en $ \eqref{eq:chebyshevproblem} $ puede verse de la siguiente manera. La minimización se encuentra en la reducción del valor de $ f(x) $. La maximización se observa al momento de encontrar una diferencia $ |d(w) - F(\tilde{x}, w)| $ lo suficientemente amplia para sobrepasar $ z $ y satisfacer la restricción $ g(x, w) $.

[//]
for $ p_j \in \bb{N} $ and $ d \in \bb{N} \cup \lbrace 0 \rbrace $. $$\tag*{$\Box$}$$

Note that at this point Lagrange Multipliers are only used when your constrains are equalities. But if you make the equality constant a parameter that belongs to an continuous real interval you can get a SIP. For instance, for problem $ \min f(x, y)\; \textrm{s.t.}\; x^2 + y^2 = 3 $ you can get a constraint area by making $ x^2 + y^2 = c $ where $ c \in [0, 3] $.