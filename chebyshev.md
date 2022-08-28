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
$$

# Numerical solution of function approximations with semi-infinite programming

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
  \min_{ x \in K_{n-1}} \max_{w \in \Omega} |d(w) - F(x, w)|
  \label{chebyshevproblem}
\end{equation}
$$
where 
  $ K_{n-1} \subseteq \mathbb{R}^{n-1} $ and
  $ \Omega \subseteq \mathbb{R}^{m} $ are non-empty and compact sets, 
  $ d: \Omega \mapsto \mathbb{R} $ and
  $ F: K_{n-1}\times\Omega \mapsto \mathbb{R} $ are smooth functions given as input to the problem. Here, $ d(w) $ and $ F(x, w) $ represents the function to aproximate and the approximation function respectively, where $ x $ is the vector of parameters or coefficients that we want to optimize. Note that the approximation error given by $ |d(w) - F(x, w)| $ is not squared and computed linearly.

### 1.2 Semi-infinite programming problem _(SIP)_
In general terms, a SIP is an optimization problem described as follows:

$$
\begin{equation}
\displaylines{
  \min \; f(x) \newline
  \textrm{s.t.} \; g(x, w) \leq 0, \; x \in \bb{R}, \; w \in \Omega, \; |\Omega| = \infty 
}
\label{def-sip}
\end{equation}
$$

where 
  $ f: K \mapsto \bb{R} $ and
  $ g: K \times \Omega \mapsto \bb{R} $ are smooth functions ($ g $ will be referred as semi-infinite constraint), with
  $ \Omega \subseteq \bb{R}^{n-1} $, 
  $ K \subseteq \bb{R}^{m} $. In general, this problem definition can have other kind constraints, i.e. equality, inequality and several semi-infinite constraints, but at least must be one semi-infinite constraint to have a SIP.

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
\problemMinimizeSingle{\max(x^7 - x^6, x^5 - x^4, x^3 - x^2, x - 1)}{x \in \bb{R}}
\label{capsipex2}
\end{equation}
$$
can be reformulated into the equivalent one:
$$
\problemMinimizeSingle{t \in \bb{R}}{x^7 - x^6 \leq t,\; x^5 - x^4 \leq t,\; x^3 - x^2 \leq t,\; x - 1 \leq t,\; x \in \bb{R}}.
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
  \txt{s.t.} \quad \tilde{x} \in K^{n-1},\; t \in \mathbb{R},\; g(x, w) := |d(w) - F(\tilde{x}, w)| \leq t ,\; w \in \Omega. \newline
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
F(x, w; d) = x^Tz(w) = (x_1, ..., x_{n-1})^T(z_1(w), ..., z_{n-1}(w))
\label{multivariate-approx}
\end{equation}
$$
where $ x \in K_{n-1} $ is the polynomial coefficients tuple, $ w \in \Omega $, $ z: \Omega \mapsto K_{n-1} $ a vector function with $ \Omega \subseteq \mathbb{R}^{m} $, $ K_{n-1} \subseteq \mathbb{R}^{n-1} $, and each $ z_i $ the following monomial:
$$ 
z_i = w_1^{p_1} w_2^{p_2} ... w_{m}^{p_{m}},\; \sum_{j=0}^{m} p_j \leq d.
$$
with $ i = 1, ..., n-1 $ and $ p_j \in \bb{N} $. $ \quad \Box $

Note that the number of monomials $ k $ that will add up the polynomial is:
$$
k = \binom{m + d}{d}.
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
In order to compare the results of this work with [[4]](#ref4) and [[8]](#ref8), let $ d(w) $ be any of the following approximation targets:

**Example 1.6.1**{: #example1}   
Refer to CAP-SIP $ \eqref{capsip} $ and let $ d(w) := d_1(w) $, $ w \in \Omega $ and $ F(x, w) $ as in $ \eqref{multivariate-approx} $:
$$
\displaylines{
  d_1(w) := d(w_1, w_2) = \log(w_1 + w_2)\sin(w_1) \newline
  w \in \Omega = [0, 1] \times [1, 2.5] \newline
  F(x, w; 2) = 
     x_1w_1^0w_2^0
    +x_2w_1^0w_2^1
    +x_3w_1^0w_2^2
    +x_4w_1^1w_2^0
    +x_5w_1^1w_2^1
    +x_6w_1^2w_2^0
}
$$

**Example 1.6.2**{: #example2}    
Refer to CAP-SIP $ \eqref{capsip} $ and let $ d(w) := d_2(w) $, $ w \in \Omega $ and $ F(x, w) $ as in $ \eqref{multivariate-approx} $:
$$
\displaylines{
  d_2(w) := d(w_1, w_2) = (1 + w_1)^{w_2} \newline
  w \in \Omega = [0, 1] \times [1, 2.5] \newline
  F(x, w; 2) = 
     x_1w_1^0w_2^0
    +x_2w_1^0w_2^1
    +x_3w_1^0w_2^2
    +x_4w_1^1w_2^0
    +x_5w_1^1w_2^1
    +x_6w_1^2w_2^0
}
$$


2. SQP Method
--------------------------------------------------------------------------------------
CAP will be computed with an open source software implementation of Sequential Programming Method (SQP). In particular the implementation provided by MATLAB will be used, that can be found in _fseminf_ routine that belongs to the Optimization Package Extension. 

Current section's aim is to explain the core gears of SQP method, that relies on the concepts of _Newton's Method_ for polynomial root approximation, _Langrage multipliers_ for local constraint optimization, _Kuhn-Karush-Tucker_ optimality conditions and _Quadratic Programming_. Those methods and techniques requires that the objective functions and constraints must be smooth; this ensures a predictable algorithm behaviour because they are designed on top of the essence of _Calculus Theory_.

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
\min_{x \in \bb{R}}f(x) \txt{s.t.}
\left\lbrace\displaylines{
  \quad g_i(x) =    0, \; i \in \cc{E} \quad \newline
  \quad g_i(x) \leq 0, \; i \in \cc{I} \quad
}\right\rbrace
$$
where $ f $ and $ g_i $ are real-valued smooth functions on a subset of $ \bb{R}^n $ for $ i \in \cc{E} \cup \cc{I}$, and $ \cc{E}$ and $ \cc{I} $ are finite set of indices for equality and inequality constraints respectively.

<!-- TODO: Write down the strategy to optimize regions using lagrange multipliers -->

The equality criteria stated in previous theorem is referred as _Lagrangian function_:
$$
\begin{equation}
  \cc{L}(x, \lambda) = f(x) - \sum_{i \in \cc{E}} \lambda_{i} \nabla g_i(x)
\end{equation}
$$

### 2.3 KKT optimality conditions
Known as _First-Order Necessary Conditions_, are conditions concerned to the gradients of a local solution $ x^\star $.

**Theorem 2.3.1**{: #theorem-kkt }
_Let $ F(x, w) $ an approximation function and let x be an arbirary number in $ \Omega $. Then there exists one element in $ F_i $ such that x < 3_

### 2.4 Quadratic programming
Es un método de aproximación local

### 2.5 Sequential Quadratic Programming method SQP
Existen varios métodos SQP, el IQP y el EQP. Actualmente la librería emplea un método...

3. Numerical Experiments
--------------------------------------------------------------------------------------
At degree 5, the approximation is near enough that causes _z-fighting_ within the plot.
This let show how the approximation is good, because this phenomenon occurs when two
or more primitives are near to the render camera, which means that are almost at the
same distance. This was for $ d = (1 + W1).^W2 $.

4. Source code revision
--------------------------------------------------------------------------------------
In this section we will examine the source code of _fseminf_ routine. Having the 
knoledge of SQP Methods, we will point out what ideas are applied and in which 
parts.
After defining the approximation problem in terms of SIP, only left to pour the functions $ F(x, w) $ and $ d(w) $ to $ \eqref{chebyshevproblem} $. 

## TODO
- Change 'continuous' by 'smooth', or belongs to C^i, because 'continuos' does not  mean 'diffrenciable'.
- Expand the idea of necessary and sufficient conditions, to expand the Lagrange multipliers method to region constraints.
- Add the examples to the end of the first section.
- Review merit function at fsemiinf.m 321




## Referencias
**[1]** MathWorks - https://www.mathworks.com/help/optim/ug/fseminf.html    
**[2]** Numerical Optimization Jorge Nocedal, Stephen Wright - (2006)    
**[3]** Numerical optimization theoretical and practical - J. Bonnans, J. Gilbert, C. Lemarechal, C. Sagastizábal - (2006)   
**[4]**{: #ref4} Reemtsen R., Discretizations Methods for the Solutions of Semi-
In nite Programming Problems, J. Optim. Theory Appl, 71 (1991),
pp. 85-103.    
**[6]** Nocedal   
**[7]** Moradito
**[8]**{: ref8} Carlos Gamboa