# Report

In this miniature report, I will walk through my code along with a bit of math from the Automatica paper to display where my issue is.

Before I go into that, I'll briefly go over the overarching problem I'm having with Diciplined Convex Programming and the curvature/sign of the sub-expressions within the performance index.

The matrices that are formed via affine functions of the desision variable, ie. $\mathcal{X}_0(\Psi)$, $\mathcal{X}_w(\Psi)$, $\mathcal{U}_0(\Psi)$, $\mathcal{U}_w(\Psi)$ get multiplied by their respective transposes, whose trace is then taken in the cost function.
While $\mathrm{tr}(AA^T)$ is a convex function in $A$, it does not follow the rules of Disciplined Convex Programming as laid out by CVXPY because the sign of $AA^T$ is unknown.
From what I've read online, a way to get around this is to reformulate the cost in terms of a Kronecker product of two matrices. I am currently working on getting this to work and using existing documentation as a guide. I was hoping to get feedback on my thought process and formulation so that I may move forward.
I know that the Kronecker product is an alternative way to express the tensor product as a matrix instead of a rank 4 covariant tensor, but I fail to see how I can turn my reformulation into a Kronecker product. I have my doubts that this process would actually change my results. The matrix $\mathrm{tr](AA^T)$ itself has positive sign and curvature, but CVXPY has no way of knowing what its sign is. Currently, I am using the Boyd paper as a refrence point in order to get my code functioning and eliminate the error that I'm having with sign/curvature.

Below, I'm going to document what I've done and the issues that I am currently having. Please let me know of any comments or suggestions you have of how to make this work.

We have our standard import statements
~~~{.python}
import numpy as np
from functools import reduce
import operator
import cvxpy as cp
import scipy as sp
~~~

And we build up a generic 

N = 2
n_x = 1
n_u = 1

As = [np.array([[2]]) for k in range(N)]
Bs = [np.array([[1]]) for k in range(N)]
Cs = [np.identity(n_x) for k in range(N)]

Sigma_0 = np.identity(n_x)
Sigma_f = 0.25 * np.identity(n_x)

Qs = [0.5*np.identity(n_x) for k in range(N)]
Rs = [0.5*np.identity(n_u) for k in range(N)]
Rcss = [[0.1*np.identity(n_u) for k in range(N)]]
cbars = [5]

Q = sp.linalg.block_diag(*Qs, np.zeros((n_x,n_x)))
R = sp.linalg.block_diag(*Rs)
Rcs = [sp.linalg.block_diag(*mats) for mats in Rcss]

def transition(matrices,t,tau):
    dimension = matrices[0].shape[0]
    return reduce(operator.matmul, reversed(matrices[tau:t]), np.identity(dimension))

def fill_upper_zeros(rows):
    z = np.zeros(rows[-1][0].shape)
    N = len(rows[-1])

    blocks = [row + [z]*(N - len(row)) for row in rows]
    return np.block(blocks)

def history(Ns, As):
    '''
    Build a history matrix.
    As is the sequence of dynamics operators.
    Ns is the sequence of control or noise matrices.
    '''

    N = len(As)
    histlist = [[transition(As,k,j+1) @ Ns[j] for j in range(k)]
            for k in range(N+1)]
    return fill_upper_zeros(histlist)


H = history(Bs, As)
G = history(Cs, As)

Gamma = np.vstack([transition(As, t, 0) for t in range(len(As)+1)])


psis = [[cp.Variable((n_u, n_x)) for col in range(row+1)] +
        [np.zeros((n_u, n_x))]*(N-row) for row in range(N)]
psi = cp.bmat(psis)

X0 = np.identity(H.shape[0]) + H@psi
Xw = X0 @ G

U0 = psi
Uw = psi @ G

def quad(Mw, M0):
    return Mw@Mw.T + M0@Gamma@Sigma_0@cp.transpose(M0@Gamma) 

xterm = cp.trace(quad(Xw, X0)@Q)
uterm = cp.trace(quad(Uw, U0)@R)
import pdb; pdb.set_trace()
performance_index = cp.Minimize(xterm + uterm)

input_constraints = [(cp.trace(quad(Uw, U0) @ Rc) <= cbar)
        for (Rc, cbar) in zip(Rcs, cbars)]

tcweight = G@G.T + Gamma@Sigma_0@Gamma.T
pns = [np.zeros((n_x,n_x))]*N + [np.identity(n_x)]
PN = np.hstack(pns)
tcquad = PN @ X0
terminal_constraint = cp.PSD(Sigma_f - tcquad@tcweight@cp.transpose(tcquad))

constraints = input_constraints + [terminal_constraint]
problem = cp.Problem(performance_index, constraints)
problem.solve()

F = np.inverse(X0.value) @ psi.value
~~~
