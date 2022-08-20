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


When $a \ne 0$, there are two solutions to $ax^2 + bx + c = 0$ and they are
$$x = {-b \pm \sqrt{b^2-4ac} \over 2a}.$$


El problema _minimax_ en $ \eqref{eq:chebyshevproblem} $ puede verse de la siguiente manera. La minimización se encuentra en la reducción del valor de $ f(x) $. La maximización se observa al momento de encontrar una diferencia $ |d(w) - F(\tilde{x}, w)| $ lo suficientemente amplia para sobrepasar $ z $ y satisfacer la restricción $ g(x, w) $.