# Optimization Techniques

This repository contains MATLAB implementations of various optimization techniques from scratch. These methods can be applied to multi-dimensional functions and are designed to solve different types of optimization problems.

## Files and Optimization Methods

### 1. LagrangeNewton.m

**Method:**  
The Lagrange-Newton method is used for constrained optimization problems, where the goal is to optimize an objective function subject to equality constraints. It combines the Newton-Raphson method with Lagrange multipliers, which adds a constraint to the function being optimized by introducing new variables (Lagrange multipliers) to the system of equations.

**Advantages:**
- Effective for constrained optimization problems with equality constraints.
- Utilizes second-order information (Hessian), leading to faster convergence compared to first-order methods like gradient descent.

**Disadvantages:**
- Requires the computation of the Hessian matrix, which can be computationally expensive for large problems.
- Not suitable for problems with inequality constraints.

---

### 2. Optimisation_Secant.m

**Method:**  
The Secant method is an iterative root-finding algorithm that approximates the root of a function by using a sequence of secant lines. It is an extension of the Newton-Raphson method but does not require the computation of derivatives.

**Advantages:**
- Does not require derivatives, making it useful for functions that are difficult to differentiate.
- Typically faster than bisection or other derivative-free methods.

**Disadvantages:**
- The method can be slow to converge, especially for poorly conditioned problems.
- Requires two initial approximations, and poor choices for initial guesses may lead to divergence.

---

### 3. Optimisation_newton.m

**Method:**  
The Newton-Raphson method is a widely used optimization method that iteratively solves for the stationary points of a function using its first and second derivatives (gradient and Hessian). It is a second-order method that uses curvature information to find the minimum more efficiently than first-order methods.

**Advantages:**
- Converges faster than gradient-based methods, especially near the optimum.
- Suitable for smooth, differentiable objective functions.

**Disadvantages:**
- Requires the calculation of the Hessian matrix, which can be computationally expensive for high-dimensional problems.
- May not converge for poorly conditioned functions or functions with non-smooth behavior.

---

### 4. Projected_gradient_descent.m

**Method:**  
The Projected Gradient Descent (PGD) method is an extension of gradient descent for constrained optimization problems. It projects the gradients onto the feasible set defined by the constraints to ensure that iterates remain within the feasible region.

**Advantages:**
- Efficient for constrained optimization problems, especially when dealing with box constraints.
- Can handle a wide range of constraint types (e.g., linear, convex).

**Disadvantages:**
- The performance heavily depends on the constraints' complexity.
- May struggle with non-convex problems where the projection is computationally expensive.

---

### 5. Simplex_Algorithm.m

**Method:**  
The Simplex method is a well-known algorithm for solving linear programming problems. It iteratively moves along the edges of the feasible region to find the optimal vertex, solving the problem of maximizing or minimizing a linear objective function subject to linear constraints.

**Advantages:**
- Widely used and very effective for linear programming problems.
- Can handle large linear optimization problems efficiently.

**Disadvantages:**
- Can be slow in practice for some large-scale problems, especially if the problem is degenerate.
- Not suitable for non-linear optimization problems.

---

### 6. broyden.m

**Method:**  
Broyden's method is a quasi-Newton method for solving nonlinear equations. It approximates the Jacobian matrix iteratively, providing a good balance between convergence speed and computational cost, especially for large problems where computing the full Jacobian is infeasible.

**Advantages:**
- Does not require the full Jacobian matrix, making it more computationally efficient for large-scale problems.
- Can handle non-linear optimization problems.

**Disadvantages:**
- The convergence rate may not be as fast as full Newton's method.
- May struggle with problems where the Jacobian approximation is inaccurate.

---

### 7. grad_descent_with_step_size_search.m

**Method:**  
This script implements the Gradient Descent method with step size search. It is an iterative optimization algorithm that finds the local minimum of a function. The step size is dynamically adjusted during the optimization process to improve convergence.

**Advantages:**
- Simple to implement and works well for differentiable convex functions.
- Step size search helps in dynamically adjusting the learning rate for better convergence.

**Disadvantages:**
- Can converge slowly for large or complex problems.
- Requires careful tuning of hyperparameters, such as the maximum number of iterations and stopping criteria.

---

### 8. subgrad_descent.m

**Method:**  
The Subgradient Descent method is a variant of gradient descent designed for non-differentiable functions. It generalizes the gradient descent method by using subgradients, which can be applied to convex but non-differentiable functions.

**Advantages:**
- Useful for convex optimization problems where the objective function is not differentiable (e.g., L1 regularization).
- Can handle a broader range of optimization problems compared to standard gradient descent.

**Disadvantages:**
- Convergence is slower than standard gradient descent, especially for non-smooth functions.
- Not as efficient for highly smooth problems as gradient descent or Newton's method.

---

## Usage
Each script is a standalone implementation of the respective optimization technique. You can use these functions by passing your objective function and initial guesses as arguments. They are designed to handle multi-dimensional functions and can be customized further based on specific requirements.

