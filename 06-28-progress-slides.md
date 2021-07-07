% Progress Report
% Drew Allen McNeely
% June 28, 2021

# Reproducing the Automatica Paper


Finite-horizon covariance control for discrete-time stochastic linear systems subject to input constraints

## Work so far

- Read the paper thoroughly and rederived chunks of it by hand for better understanding
- Attempting to code the general vector state version of the controller using CVXPY as the convex optimizer
- Running into a snag with Disciplined Convex Programming

-------------------

## Output from CVX

~~~
cvxpy.error.DCPError:
Problem does not follow DCP rules. Specifically:
The objective is not DCP.
Its following subexpressions are not:
...
~~~

## Testing in pdb

![Output from pdb](/home/drew/research/summer-2021/bakolas-automatica-paper/is_dcp.png)

# Source of issues

The issue is in the convexity of the composition of the components of the cost function.
As in the paper, $\mathrm{trace}(\mathcal{A}\mathcal{A}^T)$ is convex, but $\mathcal{A}\mathcal{A}^T$ is not necessarily convex in $\mathcal{A}$.

I am currently trying to rederive the performance index to be able to utilize one of the built-in CVX atoms. Another avenue I am exploring is in looking at the code that Boyd himself uses for his paper.

-------------------
# Future Plans

- Translate the Boyd code into Python and see if I can adapt it to the Automatica paper while getting it to run correctly
- Expand the codebase to include the separation based paper
- Run the numbers from the illustrative example in the Automatica paper through the code
