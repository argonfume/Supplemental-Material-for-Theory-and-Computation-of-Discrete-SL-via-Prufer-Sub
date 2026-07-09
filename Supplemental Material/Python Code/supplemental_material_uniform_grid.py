#uniform grid regular shooting

import numpy as np
import matplotlib.pyplot as plt

#the grid
def uniform_grid(n):
    return np.linspace(0, np.pi, n + 1)

#the problem: y'' + lambda y = 0, y(0) = y'(0) = 0
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


n = 100 #number of points on grids
numeig = 10 #number of desired eigen values

grid = uniform_grid(n)

#residual
def residual(lam):
    _, ys = shoot_on_grid(grid, lam)
    return ys[-1]


#Sweep
lam_max = (numeig + 1)**2
spacing = 1000
guess_lams = np.linspace(0.1, lam_max, spacing)
F = np.array([residual(L) for L in guess_lams])

#Bracket
brackets = []
for i in range(len(guess_lams) - 1):
    if F[i] * F[i+1] < 0:
        brackets.append((guess_lams[i], guess_lams[i+1]))
        if len(brackets) == numeig:
            break

#Bisect
eigs = []
for a, b in brackets:
    eigs.append(bisect(residual, a, b, xtol=1e-11))
eigs = np.array(eigs)

#compare against known exact eigenvalues
exact = np.array([k*k for k in range(1, numeig+1)])
print(f"\nUniform grid: n={n}")
print(f"  k   shooting lambda_k     exact k^2        error")
for k, (a, e) in enumerate(zip(eigs, exact), start=1):
    print(f"  {k}   {a:18.10f}   {e:14.10f}   {abs(a-e):.2e}")

#residual plot
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(guess_lams, F, lw=0.9, label=r"$y(\pi;\lambda)$")
ax.axhline(0, color="k", lw=0.5)
for k, e in enumerate(eigs, start=1):
    ax.plot(e, 0, "ro", ms=6)
    ax.annotate(fr"$\lambda_{k}$", xy=(e, 0), xytext=(e, 0.05),
                ha="center", color="r", fontsize=10)
ax.set_xlabel(r"Eigenvalue $\lambda$"); ax.set_ylabel(r"Residual $y(\pi;\lambda)$")
ax.set_title(f"Residual sweep on uniform grid (n={n}) for {numeig} eigenvalues")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f'Uniform_Residual_n{n}_eig{numeig}.pdf', bbox_inches='tight')
plt.show()

#eigen function plot
print(f"\n{'k':>3} {'λ_k':>14} {'y(π)':>14} {'y(π)/‖y‖∞':>14}")
fig, ax_curves = plt.subplots(figsize=(12, 5))
for k, lam in enumerate(eigs, start=1):
    ts, ys = shoot_on_grid(grid, lam)
    peak = max(abs(y) for y in ys)
    rel  = ys[-1] / peak
    print(f"{k:>3} {lam:>14.6f} {ys[-1]:>14.4e} {rel:>14.4e}")
    ys_norm = [y / peak for y in ys]
    ax_curves.plot(ts, ys_norm, marker='.', ms=3, lw=1,
                   label=fr"$k={k}$, $\lambda\approx{lam:.3f}$")

ax_curves.axhline(0, color='k', lw=0.5)
ax_curves.plot(np.pi, 0, 'ro', ms=7, label="target (π, 0)")
ax_curves.set_xlabel('t'); ax_curves.set_ylabel('y (normalized)')
ax_curves.set_title(f"Eigenfunctions on uniform grid (n={n}) for {numeig} eigenvalues")
ax_curves.legend(loc='lower left', fontsize=9, ncol=4)
ax_curves.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f'Uniform_Eigenfunctions_n{n}_eig{numeig}.pdf', bbox_inches='tight')
plt.show()

from google.colab import files
files.download(f'Uniform_Residual_n{n}_eig{numeig}.pdf')
files.download(f'Uniform_Eigenfunctions_n{n}_eig{numeig}.pdf')

#uniform grid prufer shooting

import numpy as np
import matplotlib.pyplot as plt

#the grid
def uniform_grid(n):
    return np.linspace(0, np.pi, n + 1)

#scaled Prüfer: theta' = sqrt(lam) - (q(s)/sqrt(lam)) sin^2(theta)
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

n = 100 #number of points on grids
num = 10 #number of desired eigen values

grid = uniform_grid(n)

#residual
def residual(lam):
    _, ths = prufer_trace(grid, lam)
    return ths[-1] / np.pi

#Sweep
lam_max = (num + 1)**2
guess_lams = np.linspace(0.1, lam_max, 100000)
F = np.array([residual(L) for L in guess_lams])

#Bracket
brackets = []
for m in range(1, num + 1):
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

#compare known exact eigenvalues
exact = np.array([k*k for k in range(1, num+1)])
print(f"\nUniform grid (Prüfer): n={n}")
print(f"  k   Prufer lambda_k       exact k^2        error")
for k, (a, e) in enumerate(zip(eigs, exact), start=1):
    print(f"  {k}   {a:18.10f}   {e:14.10f}   {abs(a-e):.2e}")

#residual plot
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(guess_lams, F, lw=0.9, label=r"$\theta(\pi;\lambda)/\pi$")
for m in range(1, num + 1):
    ax.axhline(m, color="k", lw=0.3, alpha=0.3)
for k, e in enumerate(eigs, start=1):
    ax.plot(e, k, "ro", ms=6)
    ax.annotate(fr"$\lambda_{k}$", xy=(e, k), xytext=(e, k + 0.2),
                ha="center", color="r", fontsize=10)
ax.set_xlabel(r"Eigenvalue $\lambda$"); ax.set_ylabel(r"Prüfer phase $\theta(\pi;\lambda)/\pi$ ")
ax.set_title(f"Prüfer phase sweep on uniform grid (n={n}) for {num} eigenvalues")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f'Uniform_Prufer_Residual_n{n}_eig{num}.pdf', bbox_inches='tight')
plt.show()

#eigen function plot
print(f"\n{'k':>3} {'λ_k':>14} {'θ(π)/π':>14} {'θ(π)/π − k':>14}")
fig, ax_curves = plt.subplots(figsize=(12, 5))
for k, lam in enumerate(eigs, start=1):
    ts, ths = prufer_trace(grid, lam)
    ys = [np.sin(th) for th in ths]
    peak = max(abs(y) for y in ys)
    rel  = ths[-1] / np.pi - k
    print(f"{k:>3} {lam:>14.6f} {ths[-1]/np.pi:>14.6f} {rel:>14.4e}")
    ys_norm = [y / peak for y in ys]
    ax_curves.plot(ts, ys_norm, marker='.', ms=3, lw=1,
                   label=fr"$k={k}$, $\lambda\approx{lam:.3f}$")

ax_curves.axhline(0, color='k', lw=0.5)
ax_curves.plot(np.pi, 0, 'ro', ms=7, label="target (π, 0)")
ax_curves.set_xlabel('t'); ax_curves.set_ylabel('y (normalized)')
ax_curves.set_title(f"Eigenfunctions via Prüfer on uniform grid (n={n}) for {num} eigenvalues")
ax_curves.legend(loc='lower left', fontsize=9, ncol=4)
ax_curves.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f'Uniform_Prufer_Eigenfunctions_n{n}_eig{num}.pdf', bbox_inches='tight')
plt.show()

from google.colab import files
files.download(f'Uniform_Prufer_Residual_n{n}_eig{num}.pdf')
files.download(f'Uniform_Prufer_Eigenfunctions_n{n}_eig{num}.pdf')