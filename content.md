# 1
Un problema de programación semi-infinita (SIP) tiene la siguiente 
descripción:

[//]: # (https://www.overleaf.com/learn/latex/Spacing_in_math_mode)
[//]: # (https://tex.stackexchange.com/questions/235382/optimization-formulas-in-latex)
[//]: # (https://math.meta.stackexchange.com/questions/11720/new-line-within-mathjax)
[//]: # (https://en.m.wikibooks.org/wiki/LaTeX/Mathematics  )


$$
\begin{equation}
\displaylines{
  \min \quad f(x) \newline
  \textrm{s.t.} \quad g_{i}(x, \omega), \quad x \in \mathbb{R}, \quad \omega \in \Omega, \quad |\Omega| = \infty 
}
\label{eq:problemdef}
\end{equation}
$$

Donde para $\eqref{eq:problemdef}$
  $ f: K \mapsto \mathbb{R} $ y
  $ g: K \times \Omega \mapsto \mathbb{R} $ son funciones continuas,
  $ \Omega \subseteq \mathbb{R} $, 
  $ K \subseteq \mathbb{R} $.

También


## 1.1
### 1.1.2

$$
\begin{align}
&\min_{w,b,\xi} \quad \frac{1}{2}w^{t}w+C\sum_{i=1}^{N}{\xi_{i}}\\\\
&\textrm{subject to} \quad y_{i}(w\phi(x_{i}+b))+\xi_{i}-1\\\\
\end{align}
$$

When $a \ne 0$, there are two solutions to $ax^2 + bx + c = 0$ and they are
$$x = {-b \pm \sqrt{b^2-4ac} \over 2a}.$$