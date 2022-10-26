[//]: # (https://www.overleaf.com/learn/latex/Spacing_in_math_mode)
[//]: # (https://tex.stackexchange.com/questions/235382/optimization-formulas-in-latex)
[//]: # (https://math.meta.stackexchange.com/questions/11720/new-line-within-mathjax)
[//]: # (https://en.m.wikibooks.org/wiki/LaTeX/Mathematics  )

$$
\newcommand{\bb}[1]{\mathbb{#1}}
\newcommand{\cc}[1]{\mathcal{#1}}
\newcommand{\txt}[1]{\;\;\textrm{#1}\;\;}
\newcommand{\problemUnconstr}{\min f(x)\; \txt{s.t.}\; x \in \bb{R}^n}
\newcommand{\problemMinimizeSingle}[2]{\min #1\; \txt{s.t.}\; #2 }
\newcommand{\problemMinimizeMulti}[2]{
  \displaylines{
    \min #1 \newline
    \txt{s.t.} #2
  }
}
\newcommand{\sctnot}[2]{
  #1\mathrm{e}{#2} 
}
$$

# Numerical computation of function approximations with semi-infinite programming

>> ABSTRACT: In this document a Chebyshev's approximation to a real valued function is performed through a semi-infinite programming problem. This reformulation uses the software tools available for optimization problems to compute the approximation. In particular, the computer program employed to execute the optimization task relies heavily in Sequential Quadratic Programming (SQP) method. In order to made this document self-contained, the definitions and techniques that composes the SQP method are described. In the first section, the problem restatement into semi-infinite programming terms is detailed, and some problem examples are portrayed. The following section describes the SQP techniques and core concepts that makes the method. In the final section, the sample problems shown are computed.  
> KEYWORDS: chebyshev's aproximation, semi-infinite programming, sequential quadratic programming, constrained optimization.

## Table of Contents
[TOC]

1. Problem's definition
--------------------------------------------------------------------------------------
In this section we will provide the necessesary definitions in order to put the _Chebyshev's approximation problem_ (CAP) in terms of a _semi-infinite programming problem_ (SIP).

### 1.1 Chebyshev approximation problem _(CAP)_
The _Chebyshev approximation problem_ _(CAP)_ can be formulated in serveral ways, one of them is described as the following _minimax_ problem:
$$
\begin{equation}
  \min_{x \in K_{n-1}} \max_{w \in \Omega} |d(w) - F(x, w)|
  \label{chebyshevproblem}
\end{equation}
$$
where 
  $ K_{n-1} \subseteq \mathbb{R}^{n-1} $ and
  $ \Omega \subseteq \mathbb{R}^{m} $ are non-empty and compact sets, 
  $ d: \Omega \mapsto \mathbb{R} $ and
  $ F: K_{n-1}\times\Omega \mapsto \mathbb{R} $ are smooth functions given as input to the problem. 

Here, $ d(w) $ and $ F(x, w) $ represents the function to aproximate and the approximation function respectively, where $ x $ is the vector of parameters or coefficients that we want to optimize and $ w $ is the function input or variables. 
  
Note that the approximation error given by $ |d(w) - F(x, w)| $ is not squared and is computed linearly, and it is performed in every point of the approximation interval $ \Omega $. In fact, each _semi-infinite constraint_ represent each possible function value of $ d $ and $ F $ over the approximation interval $ \Omega $.

As a side note, most of the numerical methods shown in this document works thanks to the approximation domain set convexity, for that reason authors commonly refer and name that set with a symbol that also has a convex shape such as $ \Omega $.

### 1.2 Semi-infinite programming problem _(SIP)_
In general terms, a SIP is an optimization problem described as follows:

$$
\begin{equation}
\displaylines{
  \min \; f(x) \newline
  \textrm{s.t.} \; g(x, w) \leq 0, \; x \in K, \; w \in \Omega, \; |\Omega| = \infty 
}
\label{def-sip}
\end{equation}
$$

where 
  $ \Omega \subseteq \bb{R}^{n-1} $, 
  $ K \subseteq \bb{R}^{m} $, and
  $ f: K \mapsto \bb{R} $ and
  $ g: K \times \Omega \mapsto \bb{R} $ are smooth functions ($ g $ will be referred as semi-infinite constraint).
In general, this problem definition can have other kind constraints, i.e. equality, inequality and several semi-infinite constraints, but at least must be one semi-infinite constraint to have a SIP.

In other words, a SIP is just a minimization problem where one of the constraints is parametrized with a variable that belongs to an infinite set, leaving virtually an infinite number of constrains (one for each possible parameter value).

Also, the problem $ \eqref{def-sip} $ can be reformulated as follows:
$$
  \begin{equation}
  \displaylines{
    \min \;\; h(y) := y \newline
    \textrm{s.t.} \quad g(x, w) \leq 0, \newline
    \qquad \; f(x) \leq y
  }
  \label{def-sip2}
  \end{equation}
$$
where $ g, f, x, w, K, \Omega $ have the same types as in $\eqref{def-sip}$ and $ y \in K $. This new reformulation is equivalent since once the constraint $ g(x, w) $ is satisfied, now the least value is imposed by $ f(x) $, as in the original problem.

### 1.3 CAP in terms of SIP
The following example can help to understand the _reformulation technique_ [6:307] employed at the end of this subsection that will allow to pose CAP in terms of SIP. In the following example, this reformulation removes the non-differentiable points at $ x = 0 $ and $ x = 1 $:
$$
\begin{equation}
\problemMinimizeSingle{\max(x^2, x)}{x \in \bb{R}}
\label{capsipex1}
\end{equation}
$$
The reformulation consist only in adding the artificial variable $ t $:
$$
\problemMinimizeSingle{t \in \bb{R}}{t \geq x,\; t \geq x^2,\; x \in \bb{R}}
$$
The trick happens when minimization of $ t $ is done from above for $ x $ and $ x^2 $ i.e. $ t \geq x $ and $ t \geq x^2 $. If $ t_i^\star $ denotes a feasible $ t $, then $ t_i^\star $ must satisfy all the constraints, which means that each $ t_i^\star $ necessarily forms a set $ F $ of feasible upper bounds for the constraints set $ C =  \lbrace x \in \bb{R} | \; t \geq x,\; t \geq x^2,\; t \in \bb{R} \rbrace $. Since $ F $ is a partially ordered set, the query for the least element (minimum) can be made. Moreover, by the _supremum property_ $ F $ has an infimum (minimum) value because $ F $ is bounded below by $ C $.

The same idea can be applied easily for several constraint rules, for instance:
$$
\begin{equation}
\displaylines {
  \min \max (       x^7 - x^6, \newline
    \qquad \qquad   x^5 - x^4, \newline
    \qquad \qquad   x^3 - x^2, \newline
    \qquad \qquad   x \;\, - 1 ) \newline
    \txt{with}      x \in \bb{R}
}
\label{capsipex2}
\end{equation}
$$
can be reformulated into the equivalent one:
$$
\displaylines{ 
  \min {t \in \bb{R}} \newline
  \txt{s.t.} \;\; x^7 - x^6 & \leq t, \newline
  \qquad \quad x^5 - x^4 & \leq t, \newline
  \qquad \quad x^3 - x^2 & \leq t, \newline
  \qquad \quad x \;\, - 1 & \leq t
}
$$

And naturally, the same technique used in $ \eqref{capsipex1} $ and $ \eqref{capsipex2} $ can be applied for the semi-infinite programming problem in $ \eqref{capsipex3} $:
$$
\begin{equation}
\min_{x \in \bb{R}} \max_{i \in \bb{N}} g_i(x)
\label{capsipex3}
\end{equation}
$$
that is equivalent to the following problem:
$$
\min t \in \bb{R} \txt{s.t.} g_i(x) \leq t,\; i \in \bb{N}.
$$

Now, problem $ \eqref{chebyshevproblem} $ can be expressed without loss of generality **[7:200]** in the following terms 
$$
\displaylines{
  \min \quad f(x) := t \txt{with} \; x := (\tilde{x}, t)\;  \newline
  \txt{s.t.} \quad \tilde{x} \in K^{n-1},\; t \in \mathbb{R}, \newline
  \qquad \quad \; g(x, w) := |d(w) - F(\tilde{x}, w)| \leq t ,\; w \in \Omega. \newline
}
$$
Note that the objective function and also the restriction $ g(\tilde{x}, w) $ are both linear. Finally, the previous problem can be reformulated again by splitting the absolute value restriction $ g(\tilde{x}, w) $ to get the kinks removed as in $ \eqref{capsipex1} $:

$$
\begin{equation}
\displaylines{
  \min \quad f(x) := t \txt{with} \; x := (\tilde{x}, t)\;  \newline
  \txt{s.t.} \quad \tilde{x} \in K_{n-1},\; t \in \mathbb{R} \newline
  \qquad \quad g_1(x, w) := d(w) - F(\tilde{x}, w) \leq t  \newline
  \qquad \quad g_2(x, w) := F(\tilde{x}, w) - d(w) \leq t  \newline
  \qquad \quad w \in \Omega.
}
\label{capsip}
\end{equation}
$$

This reformulation will allow to compute Chebyshev's approximations with computer software and frameworks available to solve semi-infinite programming problems. In next section the problem model components are defined.

### 1.4 Approximation function definition
The problem shown in $ \eqref{capsip} $, requires that both _approximation function_ $ F(x, w) $ and _target function_ $ d(w) $ must be provided as problem's input. From now, the next [definition](#def-approximation-function) will provide the function employed in numerical examples to approximate the example targets.

**Definition 1.4.1** _(Multivariate approximation)_{: #def-approximation-function}     
For given $ d \in \bb{N} \cup \lbrace 0 \rbrace $, let $ F $ be the _multivariate polynomial_:
$$ 
\begin{equation}
\displaylines {
  F(x, w; d) &= x^Tz(w) \newline
            &= (x_1, .., x_{n-1})^T(z_1(w), .., z_{n-1}(w))
}
\label{multivariate-approx}
\end{equation}
$$
where $ x \in K_{n-1} $ is the polynomial coefficients tuple, $ w \in \Omega $, $ z: \Omega \mapsto K_{n-1} $ a vector function with $ \Omega \subseteq \mathbb{R}^{m} $, $ K_{n-1} \subseteq \mathbb{R}^{n-1} $, and each $ z_i $ the following monomial:
$$ 
z_i = w_1^{p_1} w_2^{p_2} ... w_{m}^{p_{m}},\; \sum_{j=0}^{m} p_j \leq d.
$$
with $ i = 1, ..., n-1 $ and $ p_j \in \bb{N} $. $ \quad \Box $

For instance $ F(x, w; 2) $ with $ x \in K_{n-1} \subseteq \bb{R}^6 $ and $ w \in \Omega \subseteq \bb{R}^2 $:
$$
\displaylines {
  F(x, w; 2) & = x_1w_1^0w_2^0 + x_2w_1^0w_2^1 + x_3w_1^0w_2^2 \newline
    & +x_4w_1^1w_2^0 + x_5w_1^1w_2^1 + x_6w_1^2w_2^0.
}
$$

Note that the number of monomials $ k $ that will add up the polynomial is:
$$
\begin{equation}
k = \binom{m + d}{d}.
\label{monomnum}
\end{equation}
$$
Since each monomial has a coefficient and the factors involved in the dot product for $ \eqref{multivariate-approx} $ must be dimensional consistent, the coefficients tuple dimension $ dim(x) = k $, that is $ n - 1 = k $.

### 1.5 Problem model specification
Aiming to employ the SIP reformulation $ \eqref{capsip} $ within the SQP framework, here we define the problem model (_objective function_, _decision variables_ and _constraints_) to feed the SQP method exposed in [section 2](#2-sqp-method). 

**Definition 1.5.1** _(Objective function)_{: #objective-function}   
Recalling the SIP minimization problem $ \eqref{capsip} $, the objective function for CAP is
$$
f(\tilde{x}, t) = t
$$ 
where $ f: K_{n-1} \times \bb{R} \mapsto \bb{R}. \quad \Box $

**Definition 1.5.2** _(Decision variables)_{: #decision-variables}   
Recalling the SIP minimization problem $ \eqref{capsip} $, the decision variables for CAP are $ \tilde{x} $ and $ w $. $ \Box $

**Definition 1.5.3** _(Constraints)_{: #constraints}   
Recalling the SIP minimization problem $ \eqref{capsip} $, the constraints for CAP are:
$$
  \displaylines{
    \tilde{x} \in K_{n-1}, \newline
    t \in \mathbb{R}, \newline
    w \in \Omega, \newline
    d(w) - F(\tilde{x}, w) \leq t  \newline
    F(\tilde{x}, w) - d(w) \leq t.
  }
$$
where the target function $ d(w) $ is a given preset example (see next [subsection](#16-problem-instances)) and approximation function $ F $ is defined as in $ \eqref{multivariate-approx} $: $ F(\tilde{x}, w) := F(\tilde{x}, w; d) $. Also, note that the _semi-infinite_ parameter is $ w $; the approximation function must be near to the target function for each posible value of $ \Omega $.  $ \Box $

### 1.6 Problem instances
In order to compare the results of this work with [[4]](#ref4) and [[8]](#ref8), refer to CAP-SIP $ \eqref{capsip} $ and let $ F(x, w) $ as in $ \eqref{multivariate-approx} $, $ k $ as in $ \eqref{monomnum} $ and let $ d(w) $ with $ w \in \Omega \subseteq \bb{R}^m $ be any of the following approximation targets:

**1.6.1. Example 1**{: #example1}   
$$
\displaylines{
  d(w) := \log(w_1 + w_2)\sin(w_1) \newline
  w \in \Omega = [0, 1] \times [1, 2.5] \newline
  F(x, w; i) = \prod_{j=1}^{k} x_j w_1^{p1} w_2^{p2}; \newline
  i \in \lbrace 2, 3, 4, 5, 6 \rbrace
} 
$$

**1.6.2 Example 2**{: #example2}    
$$
\displaylines{
  d(w) := (1 + w_1)^{w_2} \qquad \qquad \newline
  w \in \Omega = [0, 1] \times [1, 2.5] \newline
  F(x, w; i) = \prod_{j=1}^{k} x_j w_1^{p1} w_2^{p2}; \newline
  i \in \lbrace 2, 3, 4, 5, 6 \rbrace
}
$$

**1.6.3 Example 3**{: #example3}   
$$
\displaylines{
  d(w) := e^{(w_1^2 + w_1w_2)} \newline
  w \in \Omega = [-1, 1] \times [-1, 1] \newline
  F(x, w; 9) = \prod_{j=1}^{k} x_j w_1^{p1} w_2^{p2}; \newline
}
$$

**1.6.4 Example 4**{: #example4}  
$$
\displaylines{
  d(w) = \cos(w_3 (1 + w_1))^{w_2} \newline
  w \in \Omega = [0, 1] \times [1, 2] \times [0, 1] \newline
  F(x, w; i) = \prod_{j=1}^{k} x_j w_1^{p1} w_2^{p2} w_3^{p3} \newline
  i \in \lbrace 2, 3, 4, 5 \rbrace
}
$$

**1.6.5 Example 5**{: #example5}   
$$
\displaylines{
  d(w) = \Big| log\frac{w_1 w_2 + 1}{x_1 + 0.5} \Big| x_2^{\frac{x_3 + 1}{2}} \newline
  w \in \Omega = [0, 1] \times [0, 1] \times [0, 1] \newline
  F(x, w; i) = \prod_{j=1}^{k} x_j w_1^{p1} w_2^{p2} w_3^{p3} \newline
  i \in \lbrace 2, 3, 4, 5 \rbrace
}
$$

2. Numerical Optimization Background
--------------------------------------------------------------------------------------
CAP will be computed with an open source software implementation of Sequential Programming Method (SQP). In particular the implementation provided by MATLAB will be used, that can be found in _fseminf_ routine that belongs to the Optimization Package Extension. 

Current section's aim is to explain the core gears of SQP method, that relies on the concepts of _Newton's Method_ for polynomial root approximation, _Langrage multipliers_ for local constraint optimization, _Kuhn-Karush-Tucker_ optimality conditions and _Quadratic Programming_. Those methods and techniques requires that the objective functions and constraints must be smooth; this ensures a predictable algorithm behaviour because they are designed on top of the essence of _Differential Calculus_.

### 2.1 Newton's Method

The Newton's method is a numerical method that approximates the roots of a smooth function (i.e. $ x $ where $ f(x) $ vanishes). Since this is a necessary condition for a maximizer point, this method is employed intensively within the optimization theory.

Having the nonlinear uncostrained minimization problem $ \problemUnconstr $ one necessary condition for the optimal point $ x^\star $ is $ g(x^\star) = \nabla{f(x^\star)} = 0 $. That means we have a system of $ n $ non-linear equations that must be equal to $ 0 $ in order to minimize $ f(x) $.

The idea of the method is to employ a linear approximation for the point that should be the function's root. This is written as usual:
$$
\begin{equation}
g(x + h) = g(x) + g'(x)h + R(x, h)
\end{equation}
$$

Where $ h $ is a step vector from the point $ x $ that we want to approximate and $ R(x, h) $ is the remainder of this linear approximation. Then, with this approximation the root can be found as:
$$
\begin{equation}
\displaylines{
  g(x) + g'(x)h &= 0 \newline
  g'(x)h &= -g(x) \newline
  h &= -\frac{g(x)}{g'(x)h} \newline
}
\end{equation}
$$

Since this result is an approximation to the true $ x $ that vanishes $ g(x) $, the value $ x^N = x + h $ can be employed as starting point for the next iteration. The method finishes when a precision threshold is met.

As pointed out by several authors **[3]** this is one of the most important techniques in numerical optimization because of its fast rate of convergence. In fact, some optimization books (**[2]**, **[3]**) devote at least one chapter to develop better convergence rates and to lease techinique's undesired behaviors.

As stated before, and for the following subsections, here we will only provide an overview of the core method. Further details can be found in the addressed references.

### 2.2 Quasi Newton's methods (BFGS)
These methods are an alternative to compute the zeroes of functions. They appear as a novel technique when _Jacobian_ or _Hessian_ matrices are hard to compute. The method that is currently employed within the optimization software tool employs BFGS as method to find zeros. BFGS has the main advantage that it does not use second derivates to perform the optimization, for that reason it is an efficient algorithm to compute the CAP-SIP.

At its core, the BFGS method tries to find the zeros iteratively of the following quadratic model at $ f_k = f(x_k) $:

$$
q(p) = f_{k} + \nabla {f_{k}^T}{p} + \frac{1}{2} p B_{k} p^T
$$

by updating the matrix $ B_k $ in a way that uses the gradient information obtained at each iteration and merging it to the matrix $ B_{k + 1} $ and avoid to recompute the $ n^2 $ derivatives ($ n $ is the number of variables of $ f $) that composes the $ B $ matrix. The minimizer's direction can be found at $ p_k $:
 
$$
  p_k = -B_{k}^{-1} \nabla f_{k} 
$$ 

where $ p_k $ is used as search direction to update the model to this next domain point:

$$
  x_{k + 1} = x_{k} + \alpha_{k} p_{k}
$$

where the step length $ \alpha_{k} $ is chose to met the Wolf condtions and $ k + 1 $ refers to model's new components that we want to compute. 

The actual method follows the W. C. Davidon's **[2]** idea of preserving the previous iteration gradient in order to incorporate curvature information to the $ B_k $ matrix:

$$ 
\displaylines{
  \nabla q_{k + 1} (-\alpha p_k) & = \nabla f_k \newline
    & = \nabla f_{k + 1} - \alpha B_{k + 1} p_k
}
$$

By re-arranging we get:

$$
  B_{k+1} \alpha_{k} p_{k} = \nabla f_{k + 1} - \nabla f_{k}
$$

if we define $ s_k = x_{k + 1} - x_{k} = \alpha_{k} p_{k} $ and $ y_k = \nabla f_{k+1} - \nabla f$, the previous expression simplifies to the _secant equation_:

$$
  B_{k+1} s_{k} = y_{k}
$$

The next important step is made thanks to the _Sherman-Morrison-Woodbury_ formula, that states the following theorem.

**Theorem 2.2.1** _(Sherman-Morrison-Woodbury formula)_{:#thrm-smw-formula}    
Let $ U $ and $ V $ be matrices in $ \bb{R}^{n \times p} $ for some $ p $ between $ 1 $ and $ n $. If we define 
$$
  \hat{A} = A + U V^{T}
$$
then $ \hat{A} $ is nonsingular if and only if $ (I + V^{T} A^{-1} U) $ is nonsingular, and in this case we have
$$
\begin{equation}
  \hat{A^{-1}} = A^{-1} - A^{-1} U (I + V^{T} A^{-1} U)^{-1} V^{T} A^{-1} .
  \label{eq-smw-form} 
  \quad \Box
\end{equation}
$$

If $ H_{k} = B_{k}^{-1} $ and following equation $ \eqref{eq-smw-form} $ we get that:
$$
  H_{k+1} = H_{k} - \frac{H_k y_k y_k^T H_k}{y_k^T H_k y_k} - \frac{s_k s_k^T}{y_k^T y_k}. 
$$
the inverse matrix of $ B_k $ can be computed in a cheap way only throught matrix multiplications.

The BFGS also imposes that requirement that at each step $ H_{k+1} $ must be symmetric and positive definite, and must satisfy the _secant equation_:
$$  
  H_{k+1} y_k = s_k
$$

Also, it is desired that at each step the nearest $ H_{k+1} $ is computed. This definition of nearest $ H_{k + 1} $ to the previous $ H_{k} $ is stated with the following minimization problem:
$$
  \min_{H} || H - H_{k} || \txt{s.t.} H = H^T, H s_k = y_k
$$

whose unique solution is given by:
$$
\begin{equation}
  H_{k + 1} = (I - \rho_{k} s_{k} y_{k}^T) H_k (I - \rho_{k} y_{k} s_{k}^T) + \rho_{k} s_{k} s_{k}^T
  \label{eq-hkk}
\end{equation}
$$

Now, here is the algorithm BFGS Quasi-Newton method:

```{.python}
def BFGS(x0, eps, H0):
  """
  x0 : starting point
  eps: convergence tolerance
  H0 : inverse Hessian approximation 
  """
  
  k = 0
  while abs(gradient(fk) > eps):
    # compute search direction.
    pk = matmult(-Hk, gradient(fk))

    # ak is computed from a line search procedure to satisfy 
    # Wolfe condtions.
    xkk = xk + dot(ak, pk)
    sk = xkk - xk
    yk = grad(fkk) - grad(fk)

    # Hkk is computed with equation (13)
    Hkk = computeHkk(Hk)
    k = k + 1
  return Hkk
```

### 2.2 Lagrange Multipliers
This is one of the fundamental techniques of constrained optimization. Here the technique is reviewed incrementally, starting from one equality constraint, then with several equality constraints and finally with unequality constraints.

This method appears naturally when the unconstrained problem $ \min_{x \in \bb{R}} f(x) $ becomes constrained by one equality: 
$$
\displaylines {
  \min_{x \in \bb{R}^n} f(x) \txt{s.t.} g(x) = 0
}
$$ 
where $ f $ and $ g $ are real-valued smooth functions on a subset of $ \bb{R}^n $. At its core the method imposes a necessary condition to any critical point $ x^\star $, as stated in the following theorem. 

**Theorem 2.2.1** _(Lagrange Multipliers)_{: #thrm-lagrange-mult-uni}    
Let $ f: K^n \mapsto \bb{R} $ and $ g: K^n \mapsto \bb{R} $ be $ C^1 $ real functions, $ K^n \subseteq \bb{R}^n $, $ x^\star \in K $, $ g(x^\star) = c $, $ S = \lbrace x \in \bb{R} \;|\; g(x) = c \rbrace $ (i.e. the level set) and $ \nabla g(x^\star) ≠ 0 $. If $ f|S $ ($ f $ restricted to $ S $) has an optimal value at $ x^\star $, then there is a real number $ \lambda $ such that
$$
  \nabla f(x^\star) = \lambda \nabla g(x^\star). \quad \Box
  \label{lagrangeoptcriteria}
$$

Seemingly, if there are several equality constraints, for the problem:
$$
  \min_{x \in \bb{R}^n} f(x) \txt{s.t} g_i(x) = 0, \; i \in \cc{E}
$$
where $ f $ and $ g $ are real-valued smooth functions on a subset of $ \bb{R}^n $ and  $ \cc{E} $ is a finite set indices, the previous theorem can be extended as follows.

**Theorem 2.2.2** _(Lagrange Multipliers)_{: #thrm-lagrange-mult}    
If $ f $ has a maximum or minimum at $ x^\star $ on $ S = \lbrace x \in \bb{R} \;|\; g_i(x) = 0, i \in \cc{E} \rbrace $ and the vectors $ \nabla g_i $ for $ i \in \cc{E} $ are linearly independent, then must exist constants $ \lambda_i $  such that:
$$
\nabla f(x^\star) = \sum_{i \in \cc{E}} \lambda_{i} \nabla g_i(x^\star). \quad \Box
$$

Naturally, Lagrange multipliers can be employed to solve optimization problems that also involves unequality constraints, such as this general formulation:
$$
\begin{equation}
\min_{x \in \bb{R}}f(x) \txt{s.t.}
  \begin{cases}
    \; g_i(x) =    0, \; i \in \cc{E} \newline
    \; g_i(x) \leq 0, \; i \in \cc{I} 
  \end{cases}
  \label{eq-gralproblem}
\end{equation}
$$
where $ f $ and $ g_i $ are real-valued smooth functions on a subset of $ \bb{R}^n $ for $ i \in \cc{E} \cup \cc{I}$, and $ \cc{E}$ and $ \cc{I} $ are finite set of indices for equality and inequality constraints respectively.

<!-- TODO: Write down the strategy to optimize regions using lagrange multipliers -->

The equality criteria in previous theorem can be stated more compactly by introducing the following function called _Lagrangian function_:
$$
\begin{equation}
  \cc{L}(x, \lambda) = f(x) - \sum_{i \in \cc{E}} \lambda_{i} g_i(x)
  \label{eq-lagrangian}
\end{equation}
$$

And then taking the gradient respect to $ x $:
$$
\nabla \cc{L_x} (x, \lambda) = \nabla f(x) - \sum_{i \in \cc{E}} \lambda_{i} \nabla g_i(x)
$$

Finally, if $ x^{\star} $ is an optimal point, then:
$$
\begin{equation}
  \nabla \cc{L_x} (x^{\star}, \lambda) = 0
\end{equation}
$$



### 2.3 KKT optimality conditions
Known as _First-Order Necessary Conditions_, are conditions concerned to the gradients of a local solution $ x^\star $. As such, they are statements that must be true for a given $ x^\star $. At its core, this conditions are based on a step-wise strategy when the function value is optimized.

This strategy states that if you cannot find a step vector $ s $ such that $ \nabla f(x) s < 0 $ and also $ \nabla g_i(x) s \geq 0 $ for the general problem $ \eqref{eq-gralproblem} $, it means that you have reached an optimal value bacause the function cannot get a lower value while satisfying the constraints. In fact, geometrically this point is met when both $ \nabla g_i(x) $ and $ \nabla f(x) $ point to the same direction: $ \nabla f(x) = \lambda \nabla g_i(x) $ for $ \lambda \in \bb{R} $. That is $ \nabla_{x} \cc{L}(x, \lambda) $, the _Lagrangian_ gradient respect $ x $ of $ \eqref{eq-lagrangian} $.

As this method is based on _first-order approximations_ i.e. _gradients_, a qualification must be satisfied in order to capture correctly the geometric features of the neighborhood of $ x $. This qualification is provided in the following definition [[2]](#ref2).

**Definition 2.3.1** _(LICQ)_{: #def-licq}   
Given the point $ x $ and the active set $ \cc{A}(x) $ $ = $ $ \cc{E} $ $ \cup $ $ \lbrace i \in \cc{I} $ $ | $ $ c_i(x) = 0 \rbrace $, we say that the _linear independence constraint qualification_ (LICQ) holds if the set of active constraint gradients $ \lbrace \nabla c_i(x), \; i \in \cc{A}(x) \rbrace $ is linearly independent. $ \Box $

Finally, the _first-order necessary conditions_ are given in the following theorem. This conditions are often known as _Kuhn-Karush-Tucker Condtions_ or _KKT Conditions_ for short.

**Theorem 2.3.2**{: #thrm-kkt } _(First-Order Necessary Conditions)_  
Suppose that $ x^\star $ is a local solution of $ \eqref{eq-gralproblem} $, that the functions $ f $ and $ g_i $ are continuously differentiable, and that the LICQ holds at $ x^\star $. Then there is a Lagrange multiplier vector $ \lambda^\star $, with components $ \lambda_{i}^{\star} $,  $ i \in \cc{E} \cup \cc{I} $, such that the following conditions are satisfied at $ (x^\star, \lambda^\star) $:

$$
\begin{equation}
\displaylines {
  \nabla_x \cc{L}(x^\star, \lambda^\star) & = 0, \newline
  \qquad \;                 c_i(x^\star)  & = 0, \quad \forall i \in \cc{E}, \newline 
  \qquad \;               c_i(x^\star) & \geq 0, \quad \forall i \in \cc{I}, \newline
  \qquad \qquad          \lambda^\star & \geq 0, \quad \forall i \in \cc{I}, \newline
  \quad \lambda^\star        c_i(x^\star) & = 0, \quad \forall i \in \cc{E} \cup \cc{I}. \quad \Box
}
\label{kkt-conditions}
\end{equation}
$$

3. SQP Method and Software Implementation
--------------------------------------------------------------------------------------
_Sequential Quadratic Programming_ is an iterative numerical optimization method employed to compute constrained optimization problems. At each step it solves a _Quadratic Programming (QP)_ suproblem throught its the KKT conditions $ \eqref{kkt-conditions} $, that as commented before gives optimality necessary conditions for constrained optimization problems. Additionally, also applies a _quasi-Newton BFGS_ method that as commented before ensures a super-linear convergence by accumulating function's second-order behavior information.

In this section the _SQP_ will be presented and how the software tool implements the technique to solve this kind of problems. The method's implementation and the mathematical details presented were taken from official documentation [[1b]](#ref1b). 

The method at each iteration it is similar to the Newthon's method presented before, in the way that each time the optimal value is searched with the help of an approximation made to the objective and contrain functions. This approximation is made to the Hessian of the Lagrangian function through a quasi-Newton method. Then, with this information a QP problem is formulated and solved in order to obtain a search direction as described in Newthon's method. Finally, this process is repeated until a threshold criteria is met.

This process at each iteration has the following stages, presented in the next subsections as are executed in the software implementation:

1. Hessian Matrix Update
2. Quadratic Programming solution
3. Initialization
4. Line Search and Merit Function

### 3.1 Hessian Matrix Update
In each loop unit an approximation to the Hessian's Lagrangian function is made with the the following BFGS expression $ \eqref{eq-hkk} $ presented in previous chapters:

$$
\begin{equation}
  H_{k + 1} = (I - \rho_{k} s_{k} y_{k}^T) H_k (I - \rho_{k} y_{k} s_{k}^T) + \rho_{k} s_{k} s_{k}^T
\end{equation}
$$

This update is made while trying to keep the matrix $ H_{k + 1} $ positive definite in order to follow the convergence recomendations proposed by method authors. This positive- definiteness is achieved when $ y_k^T s_k > 0 $. If this condition is not met, the most negative element of $ y_k $ is halved iteratively until the dot product be positive. If the value is not positive after this procedure, $ y_k $ is updated with the following expression:
$$
  y_k = y_k + w v_i
$$
where
$$
\color{red}{
\displaylines{
  v_i = \nabla g_i(x_{k+1}) g_i(x_{k+1}) − \nabla g_i(x_k) g_i(x_k) \newline
  \txt (q_k)_i w < 0 and (q_k)_i (s_k)_i < 0, i=1,...,mvi=0
}
}
$$

### 3.2 Quadratic Programming solution
As commented at the beginning of this chapter, at each iteration a _Quadratic Programming (QP)_ problem is formulated and solved. QP is a topic within numerical optimization that studies techniques to optimize linear-constrainded cuadratic polynomials:

$$
\begin{equation}
\displaylines {
  \min_{x} \;    q(x) &    = \frac{1}{2} x^T G x + x^T c \newline
  \txt{s.t} a_{i}^T x &    = b_i, \quad i \in \cc{E},    \newline
  \quad \quad  a_{i}^T x & \geq b_i, \quad i \in \cc{I}.
}
\label{def-qp}
\end{equation}
$$

Within the literature there are many techniques to solve this special kind of approximation problems that includes _trust region_, _linear search_, _interior point_ and _active set_ methods. The technique employed in the software tool is based on the _active set_ method. Active set methods identifies which inequality constraints are active at each iteration, that means to identify the constraint set. This reduces the complexity of the problem since active constraints are treated as equality constraints that are more easy to handle.

The solution of this problem involves two stages. The first stage consists in to stablish an starting feasible point. The second stage consists in iteratively update this point in order to reach the optimality conditions. In this process an active set $ \overline{A}_k $ of constraints are mantained.

At each stage a search direction $ \hat{d}_k $ is mantained that depends on the set $ \overline{A}_k $ that also conform a space basis. Since this is an active set method, $ \hat{d}_k $ is computed always over the constraint set boundaries (active constraints). Then, a feasible subspace $ Z_k $ for $ \hat{d}_k $ is built with a set of vectors that are orthongonal to the active set $ \overline{A}_k $, i.e. $ \overline{A}_k Z_k = 0 $. When $ \hat{d}_k $ is computed with $ Z_k $ is warranteed that the next $ \hat{d}_k $ will lie on the constraint boundaries.

The matrix $ Z_k $ is built from the last $ m - l $ columns of the matrix $ \overline{A}^T_k $ QR decomposition, where $ l $ is the number of active constraints:
$$
  Z_k = Q[:, l + 1 : n]
$$
where
$$
  Q^T A^{T}_k = 
 \begin{bmatrix}
  R \\\
  0 \\\
 \end{bmatrix}
$$

When the subspace $ Z_k $ is found, a new search direction $ \hat{d}_k $ is coumputed such that minimizes the current quadratic problem at iteration $ k $. This means that $ \hat{d}_k \in gen(Z_k) $, where $ gen(Z) $ denotes the subspace generated by $ Z $.

Knowing that $ \hat{d}_k = Z_k p $ is dimentionally consistent for the QP, the objective function in $ \eqref{def-qp} $ can be written in terms of $ Z_k $ and $ p $ as:
$$
\begin{equation}
  q(p) = \frac{1}{2} p^T Z^T_k H Z_k p + c^T Z_k p.
  \label{sqp-quad}
\end{equation}
$$

Then, the derivative respect to $ p $ gives:
$$
 \nabla q(p) = Z^T_k H Z_k p + Z^T_k c.
$$

Given that the matrix $ H $ is postitive definite (as will always happen for the presented SQP method), the minimum of $ \eqref{sqp-quad} $ will be found at $ \nabla q(p) = 0 $, that is:
$$
  Z^T_k H Z_k p = - Z^T_k c.
$$

This will let the following iteration step to be taken:
$$
  x_{k+1} = x_k + \alpha \hat{d}_k,  \txt{where} d_k = Z_k p.
$$

Because the objective function has quadratic form, the step factor $ \alpha $ can be chosen from two options. The first involves $ | \alpha | = 1 $ in the same direction of $ \hat{d}_k $. If such $ \alpha $ exists, that is the solution for the QP. If not, $ | \alpha < 1 $ must be selected and a new active constraint is added to the active set $ \overline{A} $.

Finally, if $ n $ constraints satisfy the _LICQ_ definition and the solution point is far from a solution, Lagrange multipliers from $ \eqref{eq-lagrangian} $ are computed from the following nonsingular linear system of equations:
$$
  A^T_k \lambda_k = c + H x_k
$$

If $ \lambda_k > 0 $ for every $ \lambda_i $ that solves the previous system of equations, $ x_k $ is the solution of the QP $ \eqref{def-qp} $; else, that $ \lambda_i $ is discarded, the set $ \overline{A} $ is modified and a new iteration is started.


### 3.3 Initialization
### 3.4 Line Search and Merit Function



4. Optimization Software Walkthrough
--------------------------------------------------------------------------------------
As aforementioned, the core gear is based on the SQP method. In this section, the software implementation considerations and the previously introduced mathematical concepts are discused. 

It is well known that _chevishev approximation problem (CAP)_ is one of the first examples shown to understand the gears behind the _semi-infinite programming (SIP)_ problems. Looking at the equation $ \eqref{chebyshevproblem} $, it is seen that its nuts and bolts are keep together in two stages: a maximization stage that looks for the biggest difference among all the approximation interval, and a second stage that takes the lowest difference while modifying the approximation function parameters. In fact, this mechanism is employed in the software optimization tool implementation [[1]](#ref1).

### 4.1 First stage: maximization

The software tool solves iteratively semi-infinite programming problems such as the described in $ \eqref{def-sip} $, after doing two major problem reformulations in the first stage at each iteration. 

First, perform a piecewise quadratic or cubic approximation of $ g $ over a discretizion of the domain set $ \Omega $. Then, reformulates each semi-infinite constraint $ g $ with a fixed $ x $ into a maximization problem, such as the core idea in _CAP_ $ \eqref{chebyshevproblem} $. 

In mathematical terms the semi-infinite constraint $ g $ in problem $ \eqref{def-sip} $ should look like as this:

$$
  \max_{w \in \hat{\Omega}} \; \hat{g}(x, w) \leq 0, \; x = c, \; |\hat{\Omega}| \in \bb{N} 
$$

where $ \hat{\Omega} \subset \Omega $ is an user provided discretization of $ \Omega $ (e.g. an user defined grid), $ \hat{g} $ is a piecewise quadratic and cubic approximation of the original $ g $ and $ c \in K $, with $ \Omega $, $ g $ and $ K $ defined as in $ \eqref{def-sip} $.

By applying this reformulations, the original problem with infinitely constraints was translated into a problem that has finite and piecewise polinomial constraints. Once the maximization is computed, the result is passed to a constrained non-linear solver that will perform the minimization stage.

### 4.2 Second stage: minimization
The minimization stage is perfomed in a constrained non-linear solver that executes an SQP optimization algorithm, whose details were discussed in previous sections. This minimization problem is composed of the original objective function $ f $ in $ \eqref{def-sip} $ and a new constraint set made of the results computed in the maximization stage.

### 4.3 Software implementation summary
This two stages are condensed in the following algorithm that composes the actual implementation of the software tool. 

**Algorithm 3.3.1** _fseminf implementation overview_

1. piece wise approximation, then maximizes
2. compute the minimization with SQP $ \min f(x) \txt{s.t.} c(x) $ where c(x) is expanded with the results of step 1.
3. check for any stoping criteria, if not step 4.
4. update constraints and Langrange multipliers.


4. Numerical Experiments
--------------------------------------------------------------------------------------
The numerical experiments were performed to approximate 2-dimensional and 3-dimensional functions. For each example a two tile figure is provided. The left figure is an error plot that shows the difference between the CAP-SIP approximation and the target function. 

The 2-dimensional plot is rotated around the _z-axis_ to appreciate the approximation from every side. On the other hand, as 3-dimensional plots belongs to 4th dimension, the 3-dimensional levels are plotted and this is the actual animation of those pictures.

### 3.1 Example 1 _(results)_    
The results for this [example](#example1) where $ d(w_1, w_2) = \log(w_1 + w_2)\sin(w_1) $ over $ [0, 1] $ $ \times $ $ [1, 2.5] $ with a second degree polynomial are shown in the following picture:

![gif image](results/ex1.gif "2d degree approx.")
  >> **Figure 1:**{# #fig1} Left image shows the absolute difference between __target function__ and __approximation__ values. Right animation shows the target function __surface__ (in color) and the computed approximation __grid__ (in white), seen from different angles. 

This example is also developed in [[4]](#ref4) and [[8]](#ref8). In [4] the author approximates incrementally the target function with a 7th degree polynomial on $ 12 $ iterations with a minimum error of $ \sctnot{1.00478}{-5} $. Also, in [8] the approximation is computed with a lowest error of $ \sctnot{2.8063}{-2} $ in $ 95 $ iterations. The pictured approximation was obtained with this document's method in $ 93 $ iterations with a minimal error of $ \sctnot{1.5}{-6} $: 

![gif image](results/ex1a.gif "2d degree approx.")
>> **Figure 2:** 7th degree polynomial approximation for example 1.6.1 with SQP method. 

Table1 summarizes the findings. $ n $, $ d $, $ e $ stands for _number of iterations_, _polynomial degree_ and _error_ respectively; _DM_  stands for _Discretization Methods_ the method employed in [[4]](#ref4); and _EA_ stands for _External approximations_ the method employed in [[8]](#ref8):

d     |SQP n | DM n | EA n  | SQP e   | DM e    | EA e    |
------|------|------|-------|---------|---------|---------|
2     | 22   |  4   | -     | 2.81e-2 | 2.80e-2 | -       |
3     | 35   |  6   | -     | 3.50e-5 | 3.47e-3 | -       |
4     | 161  |  4   | -     | 9.00e-4 | 6.96e-4 | -       |
5     | 91   | 11   | -     | 2.00e-3 | 1.62e-4 | -       |
6     | 92   |  8   | -     | 2.50e-3 | 3.96e-5 | -       |
7     | 93   | 13   | 95    | 1.50e-3 | 1.00e-5 | 2.80e-1 |

>> **Table1**. The results of the method exposed in this document compared to other authors.

Now, some clarifications about the meaning of $ e $. The technique employed to compute the results is based in the _minimax_ problem $ \eqref{chebyshevproblem} $. Then the maximization part of the problem was refomulated into the artificial variable $ t $ in $ \eqref{capsip} $, being now $ t $ the objective function. The value of $ e $ is the minimal value that $ f(x, t) = t $ holds after the minimization process, this means that $ e $ is the lower difference among the maximum differences that were obtained while perfoming the optimization process. This value can be seen in the brightests zones within the 2-dimensional error plots shown in the left side of each figure. Moreover, the value that appears in the table is the maxiumum value that the colorbar displays.

### 3.2 Example 2 _(results)_   
For this problem instance the function $ d(w_1, w_2) = (1 + w_1)^{w_2} $ is approximated within the interval $ [0, 1] \times [1, 2.5] $. The **Figure2** shows the approximation achieved with a 2nd-degree polynomial in $ 27 $ iterations and with a minimum error $ e = 0.1776 $: 

![example2](./results/ex2.gif)
> **Figure3:** 2nd-degree polynomial approximation for $ d(w_1, w_2) = (1 + w_1)^{w_2} $

As with the previous problem, next figure shows the approximation with a 7th-degree polynomial. Note that the patterns formed in error plot are similar to those displayed in previous example, due to the polynomial's degree. Additionaly, check that the error also improves.

![example2a](./results/ex2a.gif)
> **Figure4:** 7th-degree polynomial approximation for $ d(w_1, w_2) = (1 + w_1)^{w_2} $

The following table puts side-to-side the results computed by other authors:

d     |SQP n | DM n | EA n  | SQP e   | DM e    | EA e    |
------|------|------|-------|---------|---------|---------|
2     | 27   | 4    | -     | 1.77e-1 | 1.77e-1 | -       |
3     | 50   | 5    | -     | 3.66e-2 | 3.65e-2 | -       |
4     | 82   | 5    | -     | 7.40e-3 | 4.67e-3 | -       |
5     | 88   | 7    | -     | 2.00e-3 | 7.38e-4 | -       |
6     | 93   | 6    | -     | 2.48e-3 | 7.67e-5 | -       |
7     | 95   | 8    | 48    | 9.80e-3 | 8.80e-6 | 2.80e-1 |

Starting from 2nd-degree polynomial approximations in some examples, the polynomial approximation plot is near enough to target surface  that causes [_z-fighting_](https://en.wikipedia.org/wiki/Z-fighting) within the plot. This validates the approximation quality, because this phenomena occurs when two graphical primitives are extremely close. This picture highlights the _z-fighting_ for the current example:

![z-fighting](./results/z-fight.sm.png)
> **Figure5:** Z-fighting on 5th-degree polynomial approximation. The white grid is the approximation. The color surface is the target function. 

### 3.3 Example 3 _(results)_    
For this approximation example a 9th-degree polynomial was chosen. The results compared to [[4]](#ref4) are discussed after following figure:

![example3](./results/ex3.gif)
> **Figure6:** A 9th-degree polynomial approximation for $ d(w) = e^{(w_1^2 + w_1w_2)} $.

This approximation was achieved after $ 96 $ iterations and has a maximum error of $ \sctnot{8.0}{-4} $, on the other hand [[4]](#ref4) has an error of $ \sctnot{7.35}{-1} $ in $ 10 $ iterations. The aim of this example is to show how good can be this method to approximate a surface with pronounced peaks with a high degree polynomial. Remark, the more black the error plot is, more uniform and smooth the approximation is.

### 3.4 Example 4 _(results)_   
In this example a 3-dimensional function is approximated. Certainly, the computation time is increased due to more dense matrices, but the number of iterations remains similar to the 2-dimensional examples.

Since 3-dimensional functions plots belongs to 4th dimension, here the 3-dimensional _level_ or _contour_ plots are shown. The third axis i.e. $ w_3 $ will give de level for the level plots. In particular, for this example the level values will be in $ [0, 1] $. 

As in previous examples, the surface-grid plot is rotated to look the approximation form several angles. Since, the 4-dimensinal plot is composed of several 3-dimensional plots, each time that the figure rotates, the level plotted is incremented by $ 0.1 $. At the same time, the error plot on the left is updated each time the level changes.

![example4](./results/ex8.gif)
> **Figure7:** 5th-degree approximation of $ d(w) = \cos(w_3 (1 + w_1))^{w_2} $.

This approximation took $ 90 $ iterations and had a minimum error of $ \sctnot{6.50}{-3} $ while the approximation made in [[4]](ref#4) had a minimum error of $ \sctnot{7.08}{-4} $ in only $ 11 $ iterations. The following table shows the number of iterations $ n $ required to approximate $ d(w) $ with a polynomial of degree $ d $ and its minimum error $ e $:

d     |SQP n | DM n | SQP e   | DM e    |
------|------|------|---------|---------|
2     | 76   |  4   | 5.12e-0 | 1.52e-1 |
3     | 56   |  7   | 2.05e-0 | 3.11e-2 |
4     | 85   | 12   | 1.95e-1 | 4.87e-3 |
5     | 90   | 11   | 6.50e-3 | 7.08e-4 |

> **Table3:** Results comparative with [[4]](#ref4).

### 3.5 Example 5 _(results)_
For this target, the author in [[4]](#ref4) stated that their method had stability ussues, and for that reason it was only approximated up to the 3rd-degree polynomial.
The following figure shows our approximation for the 5th-degree polynomial:

![example9](./results/ex9.gif)

Finally, this table shows the computation results compared to author in [[4]](#ref4):

d     |SQP n | DM n | SQP e   | DM e    |
------|------|------|---------|---------|
2     | 60   |  5   | 1.08e-1 | 8.89e-2 |
3     | 79   | 12   | 9.71e-2 | 4.81-e2 |
4     | 89   | -    | 7.49e-2 | -       |
5     | 87   | -    | 7.81e-2 | -       |

> **Table4:** Results comparative.

Only left to review and compare the iteration complexity in each one of the methods _SQP_ used in the current work, _DM_ used in [[4]](#ref4) and _AE_ in [[8]](#ref8). This comparation will be the object of a future article.


5. Source code revision
--------------------------------------------------------------------------------------
In this section we will examine the source code of _fseminf_ routine. Having the 
knoledge of SQP Methods, we will point out what ideas are applied and in which 
parts.
After defining the approximation problem in terms of SIP, only left to pour the functions $ F(x, w) $ and $ d(w) $ to $ \eqref{chebyshevproblem} $. 

6. Conclusions
--------------------------------------------------------------------------------------
The software implementation employed to compute the numerical experiments actually uses a reformulation that has the same idea that _CAP_, that is a two-stage _maximization-then-minimization_ problem. For that reason, the reformulation made in [[7]](#ref7) can be ommitted and instead use a computer software solution that handles the problem directly such as in the actual software tool implementation of _fseminf_. Avoiding this double reformulation reduces the reformulation overhead that results in better computing times.

A good technique to solve _SIP_ it is related to a good problem's reformulation. . Providing an elaborate and adapatative discretization in order to obtain a better piecewise cubic and quadratic approximation, could improve the results.

## References
**[1]**{: #ref1} MathWorks - [fseminf - Algorithms (SQP)](https://www.mathworks.com/help/optim/ug/fseminf.html)
**[1b]**{: #ref1b} Mathworks - [Constrained Nonlinear Optimization Algorithms](https://www.mathworks.com/help/optim/ug/constrained-nonlinear-optimization-algorithms.html)
**[1c]**{: #ref1c} MathWorks - []
**[2]**{: #ref2} Numerical Optimization Jorge Nocedal, Stephen Wright - (2006)    
**[3]** Numerical optimization theoretical and practical - J. Bonnans, J. Gilbert, C. Lemarechal, C. Sagastizábal - (2006)   
**[4]**{: #ref4} Reemtsen R., Discretizations Methods for the Solutions of Semi-
In nite Programming Problems, J. Optim. Theory Appl, 71 (1991),
pp. 85-103.    
**[6]** Nocedal   
**[7]**{: #ref7} Moradito
**[8]**{: #ref8} Carlos Gamboa