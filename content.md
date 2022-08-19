[//]: # (https://www.overleaf.com/learn/latex/Spacing_in_math_mode)
[//]: # (https://tex.stackexchange.com/questions/235382/optimization-formulas-in-latex)
[//]: # (https://math.meta.stackexchange.com/questions/11720/new-line-within-mathjax)
[//]: # (https://en.m.wikibooks.org/wiki/LaTeX/Mathematics  )

$$
\newcommand{\bb}[1]{\mathbb{#1}}
\newcommand{\cc}[1]{\mathcal{#1}}
\newcommand{\txt}[1]{\;\textrm{#1}\;}
\newcommand{\problemUnconstr}{\min f(x)\; \txt{s.t.}\; x \in \bb{R}^n}
\newcommand{\problemMinimizeSingle}[2]{\min #1\; \txt{s.t.}\; #2 }
\newcommand{\problemMinimizeMulti}[2]{
  \displaylines{
    \min #1 \newline
    \txt{s.t.} #2
  }
}
$$

# Numerical solution of function approximation problems as semi-infinite programming problems

>> ABSTRACT: In this document a Chebyshev's approximation to a real valued function is performed through a semi-infinite programming problem. This restatement uses the tools available for optimization problems to compute the approximation. In particular, the computer program employed to execute the optimization task relies heavily in Sequential Quadratic Programming (SQP) method. In order to made this document self-contained, the definitions and techniques that composes the SQP method are described. In the first section, the problem restatement into semi-infinite programming terms is detailed, and some problem examples are portrayed. The following section describes the SQP techniques and core concepts that makes the method. In the final section, the sample problems shown are computed.  
> KEYWORDS: chebyshev's aproximation, semi-infinite programming, sequential quadratic programming, constrained optimization.


## 1. Problem's definition
In this section we will provide the necessesary definitions in order to put the _Chebyshev's approximation problem_ (CAP) in terms of a _semi-infinite programming problem_ (SIP).

### 1.1 Semi-infinite programming problem _(SIP)_
In general terms, a SIP is an optimization problem described as follows:

$$
\begin{equation}
\displaylines{
  \min \; f(x) \newline
  \textrm{s.t.} \; g(x, w) \leq 0, \; x \in \bb{R}, \; w \in \Omega, \; |\Omega| = \infty 
}
\label{eq:problemdef}
\end{equation}
$$

where 
  $ f: K \mapsto \bb{R} $ y
  $ g: K \times \Omega \mapsto \bb{R} $ are smooth functions, 
  with
  $ \Omega \subseteq \bb{R}^{n-1} $, 
  $ K \subseteq \bb{R}^{m} $.

In other words, a SIP is just a minimization problem where one of the constraints is parametrized with a variable that belongs to an infinite set, leaving virtually an infinite number of constrains (one for each possible parameter value).

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
donde $ g, f, x, w, K, \Omega $ tienen la misma forma que en $\eqref{eq:problemdef}$ y $ z \in K $. Esta nueva representación es equivalente porque una vez se cumplen las restricciones de la función $ g(x, w) $, ahora el menor valor lo impone $ f(x) $, como en el problema original. 

### 1.2 Chebyshev approximation problem _(CAP)_
El problema de aproximación de Chebyshev puede ser formulado de maneras varias, una de ellas es la siguiente se describe con el siguiente problema _minimax_:
$$
\begin{equation}
  \min \max_{w \in \Omega} |d(w) - F(x, w)|
  \;\textrm{with}\; x \in K_{n-1}
  \label{eq:chebyshevproblem}
\end{equation}
$$
donde 
  $ K_{n-1} \subseteq \mathbb{R}^{n-1} $ y 
  $ \Omega \subseteq \mathbb{R}^{m} $ son no vacíos y compactos, 
  $ d: \Omega \mapsto \mathbb{R} $ y 
  $ F: K_{n-1}\times\Omega \mapsto \mathbb{R} $ son funciones continuas y dadas como entrada al problema.

### 1.3 CAP in terms of SIP
El problema $ \eqref{eq:chebyshevproblem} $ puede expresarse de la siguiente manera:
$$
\begin{equation}
\displaylines{
  \min \quad f(x) := z, \; x := (\tilde{x}, z),\; \tilde{x} \in \mathbb{R}^{n-1},\; z \in \mathbb{R} \newline
  \textrm{s.a.} \quad g(x, w) := |d(w) - F(\tilde{x}, w)| - z \leq 0,\; w \in \Omega \newline
}
\end{equation}
$$

El problema _minimax_ en $ \eqref{eq:chebyshevproblem} $ puede verse de la siguiente manera. La minimización se encuentra en la reducción del valor de $ f(x) $. La maximización se observa al momento de encontrar una diferencia $ |d(w) - F(\tilde{x}, w)| $ lo suficientemente amplia para sobrepasar $ z $ y satisfacer la restricción $ g(x, w) $.

Another example that can help to understand this _reformulation technique_ [6:307] can be seen when the kinks (no differentiable points) are removed from the following optimization problem:
$$
\problemMinimizeSingle{\max(x^2, x)}{x \in \bb{R}}
$$
in order to have a smooth one, by adding the artificial variable $ t $:
$$
\problemMinimizeSingle{t}{t \geq x,\; t \geq x^2,\; t \in \bb{R}}
$$


### 1.4 CAPs to solve
After defining the approximation problem in terms of SIP, only left to pour the functions $ F(x, w) $ and $ d(w) $ to $ \eqref{eq:chebyshevproblem} $. In order to compare this document's results to other autors, $ d(w) $ is defined as the following polynomial:
$$
\begin{equation}
x_1^{i_1} x_2^{i_2} = \Sigma
\end{equation}
$$

## 2 SQP Method
CAP will be computed with an open source software implementation of Sequential Programming Method (SQP). In particular the implementation provided by MATLAB will be used, that can be found in _fseminf_ routine that belongs to the Optimization Package Extension. 

Current section's aim is to explain the core gears of SQP method, that relies on the concepts of _Newton's Method_ for polynomial root approximation, _Langrage multipliers_ for local constraint optimization, _Kuhn-Karush-Tucker_ First Order conditions and _Quadratic Programming_. Those methods and techniques requires that the objective functions and constraints must be smooth; this ensures a predictable algorithm behaviour because they are designed on top of the essence of _Calculus Theory_.

### 2.1 Newton's Method

The Newton's method is a numerical method that approximates the roots of a smooth function (i.e. $ x $ where $ f(x) $ vanishes). Since this is a necessary condition for a maximizer point, this method is employed intensively within the optimization theory.

Having the nonlinear uncostrained minimization problem $ \problemUnconstr $ one necessary condition for the optimal point $ x^\star $ is $ g(x^\star) = \nabla{f(x^\star)} = 0 $. This means that we have a system of $ n $ non-linear equations that must be equal to $\min f(x)$ s.t. $g(x) = c$.

The idea of the method is to employ a linear approximation for the point that should be the function's root. This is written as usual:
$$
\begin{equation}
g(x + h) = g(x) + g'(x)h + R(x, h)
\end{equation}
$$

Where $ h $ is the extent from the point $ x $ that we want to approximate and $ R(x, h) $ is the remainder of this linear approximation. Then, with this approximation the root can be found as:
$$
\begin{align}
\displaylines{
  g(x) + g'(x)h &= 0 \newline
  g'(x)h &= -g(x) \newline
  h &= -\frac{g(x)}{g'(x)h} \newline
}
\end{align}
$$

Since this result is an approximation to the true $ x $ that vanishes $ g(x) $, the value $ x^N = x + h $ can be employed as starting point for the next iteration. The method finishes when a precision threshold is met.

As pointed out by several authors **[3]** this is one of the most important techniques in numerical optimization because of its fast rate of convergence. In fact, some optimization books (**[2]**, **[3]**) at least one chapter is devoted to develop better convergence and to lease undesired behaviors of this technique.

As stated before, and for the following subsections, here we will only provide an overview of the core method. Further details can be found in the addressed references.

### 2.2 Lagrange Multipliers
This method appears naturally when the unconstrained problem $ \min f(x) $ $ \txt{s.t.} $ $ x \in \bb{R} $ becomes constrained by equalities, i.e: $ \min f(x) $ $ \txt{s.t.} $ $ g(x) = c $, $ c \in \bb{R} $. 

At its core the method imposes a necessary condition to any critical point $ x^\star $ as follows. Let $ f: K^n \mapsto \bb{R} $ and $ g: K^n \mapsto \bb{R} $ be $ C^1 $ real functions, $ K^n \subseteq \bb{R}^n $, $ x^\star \in K $, $ g(x^\star) = c $, $ S = \lbrace x \in \bb{R} \;|\; g(x) = c \rbrace $ (i.e. the level set) and $ \nabla g(x^\star) ≠ 0 $. If $ f|S $ ($ f $ restricted to $ S $) has an optimal value at $ x^\star $, then there is a real number $ \lambda $ such that
$$
\begin{equation}
  \nabla f(x^\star) = \lambda \nabla g(x^\star).
  \label{lagrangeoptcriteria}
\end{equation}
$$

Seemingly, if there are several constraints, for the minimization problem:
$$
\begin{equation}
  \problemMinimizeMulti
  {f(x),\; x \in \bb{R}^n}
  {g_i(x) = c_i, \; c_i \in \bb{R}^n, \; i \in \cc{I} \subset \bb{N}}
\end{equation}
$$
the optimal criteria $\eqref{lagrangeoptcriteria} $ becomes:
$$
\nabla f(x^\star) = \sum_{i \in \cc{I}} \lambda_{i} \nabla g_i(x^\star)
\begin{equation}

\end{equation}
$$

Note that at this point Lagrange Multipliers are only used when your constrains are equalities. But if you make the equality constant a parameter that belongs to an continuous real interval you can get a SIP. For instance, for problem $ \min f(x, y)\; \textrm{s.t.}\; x^2 + y^2 = 3 $ you can get a constraint area by making $ x^2 + y^2 = c $ where $ c \in [0, 3] $.

### 2.3 Condiciones KKT
Knonw also as _First-Order Necessary Conditions_
Las condiciones KKT son condiciones suficientes que debe satisfacer un punto para que sea considerado como óptimo

### 2.4 Programación Quadrática
Es un método de aproximación local

### 2.5 Método SQP
Existen varios métodos SQP, el IQP y el EQP. Actualmente la librería emplea un método...

## 3 Code reading overview
In this section we will examine the source code of _fseminf_ routine. Having the 
knoledge of SQP Methods, we will point out what ideas are applied and in which 
parts.

## TODO
- Change 'continuous' by 'smooth', or belongs to C^i, because 'continuos' does not  mean 'diffrenciable'.
- Expand the idea of necessary and sufficient conditions, to expand the Lagrange multipliers method to region constraints.
- Add the examples to the end of the first section.




## Referencias
**[1]** MathWorks - https://www.mathworks.com/help/optim/ug/fseminf.html    
**[2]** Numerical Optimization Jorge Nocedal, Stephen Wright - (2006)    
**[3]** Numerical optimization theoretical and practical - J. Bonnans, J. Gilbert, C. Lemarechal, C. Sagastizábal - (2006)   
**[4]** Reemtsen R., Discretizations Methods for the Solutions of Semi-
In nite Programming Problems, J. Optim. Theory Appl, 71 (1991),
pp. 85-103.    
**[6]** Nocedal