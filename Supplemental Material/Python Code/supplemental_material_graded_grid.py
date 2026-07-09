#graded grid regular shooting

import numpy as np
import matplotlib.pyplot as plt

#the grid
def graded_grid(n, ratio):
    h = ratio**np.arange(n)
    h = h/h.sum()*np.pi
    return np.concatenate([[0.0], np.cumsum(h)])

#the problem
def shoot_on_grid(grid, lam):
    n = len(grid) - 1
    ts = [grid[0]]
    ys = [0.0]
    y, v = 0.0, 1.0
    for i in range(n):
        d = grid[i+1] - grid[i]

        k1y = v
        k1v = -lam * y

        k2y = v + 0.5*d*k1v
        k2v = -lam * (y + 0.5*d*k1y)

        k3y = v + 0.5*d*k2v
        k3v = -lam * (y + 0.5*d*k2y)

        k4y = v + d*k3v
        k4v = -lam * (y + d*k3y)

        y = y + (d/6) * (k1y + 2*k2y + 2*k3y + k4y)
        v = v + (d/6) * (k1v + 2*k2v + 2*k3v + k4v)
        ts.append(grid[i+1])
        ys.append(y)
    return ts, ys

#bisection
def bisect(f, a, b, xtol=1e-11, maxiter=200):
    fa, fb = f(a), f(b)
    if fa == 0.0: return a
    if fb == 0.0: return b
    if fa * fb > 0:
        raise ValueError("f(a) and f(b) must straddle zero")
    for _ in range(maxiter):
        c = 0.5 * (a + b)
        fc = f(c)
        if fc == 0.0 or 0.5 * (b - a) < xtol:
            return c
        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return 0.5 * (a + b)

#parameters
ratio = 1.3
numeigs  = 10
n = 100 #number of points on grids
grid = graded_grid(n, ratio)

#residual
def residual(guess_lam):
    _, ys = shoot_on_grid(grid, guess_lam)
    return ys[-1]

#Sweep
lam_max = (numeigs + 1)**2
guess_lams = np.linspace(0.1, lam_max, 1000)
F = np.array([residual(L) for L in guess_lams])

#Bracket
brackets = []
for i in range(len(guess_lams) - 1):
    if F[i] * F[i+1] < 0:
        brackets.append((guess_lams[i], guess_lams[i+1]))
        if len(brackets) == numeigs:
            break

#Bisect
eigs = []
for a, b in brackets:
    eigs.append(bisect(residual, a, b, xtol=1e-10))
eigs = np.array(eigs)

exact = np.array([k*k for k in range(1, numeigs+1)])
print(f"\nGraded grid: n={n}, ratio={ratio}")
print(f"  k   shooting lambda_k     exact k^2        error")
for k, (a, e) in enumerate(zip(eigs, exact), start=1):
    print(f"  {k}   {a:18.10f}   {e:14.10f}   {abs(a-e):.2e}")

#residual plot
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(guess_lams, F, lw=0.9, label=r"$F(\lambda)=y(\pi;\lambda)$")
ax.axhline(0, color="k", lw=0.5)
for k, e in enumerate(eigs, start=1):
    ax.plot(e, 0, "ro", ms=6)
    ax.annotate(fr"$\lambda_{k}$", xy=(e, 0), xytext=(e, 0.05),
                ha="center", color="r", fontsize=10)
ax.set_xlabel(r"$\lambda$ (guess)"); ax.set_ylabel(r"residual $y(\pi;\lambda)$")
ax.set_title(f"Residual sweep on graded grid (n={n}, ratio={ratio}) for {numeigs} eigenvalues")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f'Graded_Residual_n{n}_{numeigs}eigs.pdf', bbox_inches='tight')
plt.show()

#eigen function plot
print(f"\n{'k':>3} {'λ_k':>14} {'y(π)':>14} {'y(π)/‖y‖∞':>14}")
fig, ax_curves = plt.subplots(figsize=(30, 10))
for k, lam in enumerate(eigs, start=1):
    ts, ys = shoot_on_grid(grid, lam)
    peak = max(abs(y) for y in ys)
    rel  = ys[-1] / peak
    print(f"{k:>3} {lam:>14.6f} {ys[-1]:>14.4e} {rel:>14.4e}")
    ys_norm = [y / peak for y in ys]
    ax_curves.plot(ts, ys_norm, marker='.', ms=4, lw=1,
                    label=fr"$k={k}$, $\lambda\approx{lam:.3f}$")

ax_curves.axhline(0, color='k', lw=0.5)
ax_curves.plot(np.pi, 0, 'ro', ms=7, label="target (π, 0)")
ax_curves.set_xlabel('t'); ax_curves.set_ylabel('y (normalized)')
plt.title(f"Eigenfunctions on graded grid (ratio={ratio}, n={n}) for {numeigs} eigenvalues")
ax_curves.legend(loc='lower left', fontsize=9, ncol=4)
ax_curves.grid(alpha=0.3)
ax_curves.set_xticks(grid)
plt.tight_layout()
plt.savefig(f'Graded_Eigenfunctions_n{n}_{numeigs}eigs.pdf', bbox_inches='tight')
plt.show()

from google.colab import files
files.download(f'Graded_Residual_n{n}_{numeigs}eigs.pdf')
files.download(f'Graded_Eigenfunctions_n{n}_{numeigs}eigs.pdf')

#graded grid prufer shooting

import numpy as np
import matplotlib.pyplot as plt

#the grid
def graded_grid(n, ratio):
    h = ratio**np.arange(n)
    h = h/h.sum()*np.pi
    return np.concatenate([[0.0], np.cumsum(h)])

#scaled Prüfer
def prufer_trace(grid, lam, q=lambda s: 0.0):
    sl = np.sqrt(lam)
    f = lambda s, th: sl - (q(s) / sl) * np.sin(th)**2
    ts  = [grid[0]]
    ths = [0.0]
    th = 0.0
    for i in range(len(grid) - 1):
        d = grid[i+1] - grid[i]
        s = grid[i]
        k1 = f(s,         th)
        k2 = f(s + 0.5*d, th + 0.5*d*k1)
        k3 = f(s + 0.5*d, th + 0.5*d*k2)
        k4 = f(s + d,     th + d*k3)
        th = th + (d/6) * (k1 + 2*k2 + 2*k3 + k4)
        ts.append(grid[i+1])
        ths.append(th)
    return ts, ths

#bisection
def bisect(f, a, b, xtol=1e-10, maxiter=200):
    fa, fb = f(a), f(b)
    if fa == 0.0: return a
    if fb == 0.0: return b
    if fa * fb > 0:
        raise ValueError("f(a) and f(b) must straddle zero")
    for _ in range(maxiter):
        c = 0.5 * (a + b)
        fc = f(c)
        if fc == 0.0 or 0.5 * (b - a) < xtol:
            return c
        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return 0.5 * (a + b)

#parameters
ratio = 1.3
numeigs = 10
n = 100 #number of points on grids
grid = graded_grid(n, ratio)

#residual
def residual(guess_lam):
    _, ths = prufer_trace(grid, guess_lam)
    return ths[-1] / np.pi

#Sweep
lam_max = (numeigs + 1)**2
guess_lams = np.linspace(0.1, lam_max, 100000)
F = np.array([residual(L) for L in guess_lams])

#Bracket
brackets = []
for m in range(1, numeigs + 1):
    for i in range(len(guess_lams) - 1):
        if (F[i] - m) * (F[i+1] - m) < 0:
            brackets.append((guess_lams[i], guess_lams[i+1], m))
            break

#Bisect
eigs = []
for a, b, m in brackets:
    res_m = lambda lam, m=m: residual(lam) - m
    eigs.append(bisect(res_m, a, b, xtol=1e-10))
eigs = np.array(eigs)

exact = np.array([k*k for k in range(1, numeigs+1)])
print(f"\nGraded grid (Prüfer): n={n}, ratio={ratio}")
print(f"  k   Prufer lambda_k       exact k^2        error")
for k, (a, e) in enumerate(zip(eigs, exact), start=1):
    print(f"  {k}   {a:18.10f}   {e:14.10f}   {abs(a-e):.2e}")

#residual plot
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(guess_lams, F, lw=0.9, label=r"$T(\lambda)=\theta(\pi;\lambda)/\pi$")
for m in range(1, numeigs + 1):
    ax.axhline(m, color="k", lw=0.3, alpha=0.3)
for k, e in enumerate(eigs, start=1):
    ax.plot(e, k, "ro", ms=6)
    ax.annotate(fr"$\lambda_{k}$", xy=(e, k), xytext=(e, k + 0.2),
                ha="center", color="r", fontsize=10)
ax.set_xlabel(r"$\lambda$ (guess)"); ax.set_ylabel(r"Prüfer phase $\theta(\pi;\lambda)/\pi$")
ax.set_title(f"Prüfer phase sweep on graded grid (n={n}, ratio={ratio}) for {numeigs} eigenvalues")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f'Graded_Prufer_Residual_n{n}_{numeigs}eigs.pdf', bbox_inches='tight')
plt.show()

#eigen function plot
print(f"\n{'k':>3} {'λ_k':>14} {'θ(π)/π':>14} {'θ(π)/π − k':>14}")
fig, ax_curves = plt.subplots(figsize=(30, 10))
for k, lam in enumerate(eigs, start=1):
    ts, ths = prufer_trace(grid, lam)
    ys = [np.sin(th) for th in ths]
    peak = max(abs(y) for y in ys)
    rel  = ths[-1] / np.pi - k
    print(f"{k:>3} {lam:>14.6f} {ths[-1]/np.pi:>14.6f} {rel:>14.4e}")
    ys_norm = [y / peak for y in ys]
    ax_curves.plot(ts, ys_norm, marker='.', ms=4, lw=1,
                    label=fr"$k={k}$, $\lambda\approx{lam:.3f}$")

ax_curves.axhline(0, color='k', lw=0.5)
ax_curves.plot(np.pi, 0, 'ro', ms=7, label="target (π, 0)")
ax_curves.set_xlabel('t'); ax_curves.set_ylabel('y (normalized)')
plt.title(f"Eigenfunctions via Prüfer on graded grid (ratio={ratio}, n={n}) for {numeigs} eigenvalues")
ax_curves.legend(loc='lower left', fontsize=9, ncol=4)
ax_curves.grid(alpha=0.3)
plt.tight_layout()
ax_curves.set_xticks(grid)
plt.savefig(f'Graded_Prufer_Eigenfunctions_n{n}_{numeigs}eigs.pdf', bbox_inches='tight')
plt.show()

from google.colab import files
files.download(f'Graded_Prufer_Residual_n{n}_{numeigs}eigs.pdf')
files.download(f'Graded_Prufer_Eigenfunctions_n{n}_{numeigs}eigs.pdf')