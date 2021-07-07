import numpy as np
from functools import reduce
import operator
import cvxpy as cp
import scipy as sp

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
