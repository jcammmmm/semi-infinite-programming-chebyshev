[//]: # (https://www.overleaf.com/learn/latex/Spacing_in_math_mode)
[//]: # (https://tex.stackexchange.com/questions/235382/optimization-formulas-in-latex)
[//]: # (https://math.meta.stackexchange.com/questions/11720/new-line-within-mathjax)
[//]: # (https://en.m.wikibooks.org/wiki/LaTeX/Mathematics  )

## 1
Un problema de programación semi-infinita (SIP) tiene la siguiente 
descripción:

$$
\begin{equation}
\displaylines{
  \min \quad f(x) \newline
  \textrm{s.a.} \quad g_{i}(x, w) \leq 0, \quad x \in \mathbb{R}, \quad w \in \Omega, \quad |\Omega| = \infty 
}
\label{eq:problemdef}
\end{equation}
$$

Donde 
  $ f: K \mapsto \mathbb{R} $ y
  $ g: K \times \Omega \mapsto \mathbb{R} $ son funciones continuas,
  $ \Omega \subseteq \mathbb{R}^{n-1} $, 
  $ K \subseteq \mathbb{R}^{m} $.

También el problema $\eqref{eq:problemdef}$ puede escribirse de la siguiente manera:
$$
\begin{equation}
\displaylines{
  \min \quad z \newline
  \textrm{s.a.} \quad g_{i}(x, w) &\leq 0, \newline
  \qquad \quad f(x) &\leq z 
}
\label{eq:problemdef2}
\end{equation}
$$
donde $ g, f, x, w, K, \Omega $ tienen la misma forma que en $\eqref{eq:problemdef}$ 
y $ z \in K $. Esta nueva representación es equivalente porque una vez se cumplen las
restricciones de la función $ g(x, w) $, ahora el menor valor lo impone $ f(x) $,
como en el problema original. 

## 2
El problema de aproximación de Chebyshev puede ser formulado de maneras varias, una 
de ellas es la siguiente se describe con el siguiente problema _minimax_:
$$
\begin{equation}
  \min \max_{w \in \Omega} |d(w) - F(x, w)|
  \quad \textrm{con} \quad x \in K_{n-1}
  \label{eq:chebyshevproblem}
\end{equation}
$$
donde 
  $ K_{n-1} \subseteq \mathbb{R}^{n-1} $ y 
  $ \Omega \subseteq \mathbb{R}^{m} $ son no vacíos y compactos, 
  $ d: \Omega \mapsto \mathbb{R} $ y 
  $ F: K_{n-1}\times\Omega \mapsto \mathbb{R} $ son funciones continuas y dadas
  como entrada al problema.

## 3
El problema $ \eqref{eq:chebyshevproblem} $ puede expresarse de la siguiente manera:
$$
\begin{equation}
\displaylines{
  \min \quad f(x) := z, \; x := (\tilde{x}, z),\; \tilde{x} \in \mathbb{R}^{n-1},\; z \in \mathbb{R} \newline
  \textrm{s.a.} \quad g(x, w) := |d(w) - F(\tilde(x), w)| - z \leq 0,\; w \in \Omega \newline
}
\end{equation}
$$
El problema _minimax_ en $ \eqref{eq:chebyshevproblem} $ puede verse de la siguiente manera. 
La minimizaciónse encuentra en la reducción del valor de $ f(x) $. La maximización se observa 
al momento de encontrar una diferencia $ |d(w) - F(\tilde{x}, w)| $ lo suficientemente 
amplia para sobrepasar $ z $ y satisfacer la restricción $ g(x, w) $.