---
title: 'Mathematics and anthropogenic extinction II: PDEs'
date: 2024-09-26
permalink: /posts/2024/09/mammoth_extinction_II/
tags:
  - Partial differential equations
  - Mathematical Biology
  - Spectral methods
---

Modelling extinction events
======

Starting from the ODE model presented by Frank et. al. in ["Investigating Anthropogenic Mammoth Extinction with Mathematical Models"](https://ir.library.illinoisstate.edu/spora/vol1/iss1/3/), and discussed in previous [posts](https://javierchico.github.io/posts/2024/09/mammoth_extinction/), we include spatial variation by extending the model into a continuum. We remark that the equations have already been non-dimensionalised.

Now, \\( H(x,t)\\) and \\(M(x,t)\\) denote human and mammoth densities, and not the total number of individuals. By using a conservation of mass argument and integrating over an aribitrary domain and using the divergence theorem, we can derive a PDE system

$$\frac{\partial H}{\partial t}=A_H H(1-H) + B_H \frac{M^2H}{M^2+r^2}+\boldsymbol{\nabla}\cdot \boldsymbol{J}_H$$

$$\frac{\partial M}{\partial t}=A_M M(1-M)(\alpha M-1) - B_M \frac{M^2H}{M^2+r^2} +\boldsymbol{\nabla}\cdot \boldsymbol{J}_M$$

Here, the fluxes \\(\boldsymbol{J}\\) represent the density of each species leaving a point of the domain. We must relate the flux to the densities and their spatial derivatives. Assuming mammoths try to distribute themsleves uniformly to maximise grazing grounds and minimise intraspecies conflict, we can postulate a simple consitute law akin to Fourier's law for heat flux, 

$$ \boldsymbol{J}_M=- D_M \boldsymbol{\nabla}M$$ 

so that mammoths migrate from areas of high density to areas of low density (and their movements are not directly influenced by what humans are doing). Although we could assume a similar law for humans, we will build a "tracking" term into the flux. In particular, we want humans to go from areas of low mammoth density to areas of high mammoth density, at a rate proportional to the number of humans in that location, so that a simple flux for our movement can be given by

$$\boldsymbol{J}_H=H\boldsymbol{\nabla}M$$

This choice for the flux makes the flux nonlinear, and the fact that we do not care about the gradient of human population density makes this term slighlty atypical. As the PDEs where already nonlinear this is not that much of a sacrifise. 

We now need to decice on the boundary conditions. [Periodic boundary conditions usually allow for nice spectral methods](https://javierchico.github.io/posts/2012/08/vorticity_eq/), but continents are not periodic, so in the sake of realism we go for Neumann (no flux) boundary conditions on the edges of the domain. This means \\(\boldsymbol{n}\cdot\boldsymbol{\nabla}H,M=0\\) with n the unit outward normal to the domain. If we had a powerful enough computer, we could discretise the continental United States, compute the unit normal at each element and solve the equations in that domain. Unfortunately we do not so we just stick to the 1D case for now. 

We now must choose how to discretise the domain. The method of lines (discretising space and leaving a system of many time-dependent ODE problems) is already a good start, but the nonlienar growth, coupling and fluxes means we need to be really careful with the next step. Finite differences is likely to struggle, and finite elements and volumes are too engineery for the writer's taste, so we try an spectral method. The FFT is not rescuing us here because of our boundary conditions, but the [Discrete Cosine Transform might](https://en.wikipedia.org/wiki/Discrete_cosine_transform). 

The DCT world is quite complicated (SciPy has like 4 different types, each with its own sign convention), so my recommendation is just pickking one, type II in our case, and just sticking with it. 

Unlike in the ODE case, our investigation here is not so focused on wether mammoths go extinct or not, but we just want to see some cool solutions. Cool solutions in this field means a travelling wave: a front of humans advancing onto a mammoth population. 

![ODE_example!](/images/mammoth_figures/ODE_example_1.jpeg)

```
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def RHS(t,vector,A_H, B_H, r,alpha, A_M,B_M):
    H = vector[0];M = vector[-1]
    dHdt = A_H*H*(1-H) + B_H*M**2*H/(M**2+r**2)
    dMdt = A_M*M*(1-M)*(alpha*M-1) - B_M*M**2*H/(M**2+r**2)
    return np.array([dHdt,dMdt])

H0=0.1;M0=1
y0 = np.array([H0,M0])
params = {
    'A_H': 1.0, 'B_H': .1, 'A_M': .01, 'B_M': 0.1,'r':.1,'alpha':10
}
T = 20;Nt=300
t_span = (0, T)  # From t=0 to t=10
t_eval = np.linspace(0, T, Nt)
# Solve the PDE system
sol = solve_ivp(
    lambda t, vector: RHS(t, vector, **params),
    [0, T],
    y0,
    t_eval = t_eval,
    method='BDF'
)
H = sol.y[0,:];M = sol.y[-1,:]

```

By sweeping through values of the initial conditions and integrating both forwards and backwards in time we can construct phase diagrams, plots in \\( (H,M)\\) space where we can clearly see the stability of the fixed points. 

![ODE_s](/images/mammoth_figures/extinction_phase_plots_ODEs.jpeg "test")

Here, different values of \\(A_M \\) and \\(B_M \\) are represented in the grid. 