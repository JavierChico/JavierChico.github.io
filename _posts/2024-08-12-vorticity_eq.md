---
title: 'Two dimensional fluid mechanics in your laptop!'
date: 2024-08-10
permalink: /posts/2012/08/vorticity_eq/
tags:
  - CFD
  - fluid dynamics
  - spectral methods
---

Solving Navier-Stokes via the vorticity equation
======

In this short post I provide some tricks and code to simulate high accuracy two dimensional flows using a laptop. Starting from the Navier-Stokes equations for a Newtonian incompressible fluid:

$$\frac{\partial \boldsymbol{u}}{\partial t}+\boldsymbol{u}\cdot\nabla\boldsymbol{u}=-\nabla p +\frac{1}{\text{Re}}\nabla^2\boldsymbol{u}\quad \boldsymbol{\nabla}\cdot\boldsymbol{u}=0$$

The pressure gradient can be problematic in CFD algorithms, as there are no explicit conditions for the pressure at boundaries. Hence, a commonly used trick is to take the curl of the momentum equations, leading to the [vorticity equation](https://en.wikipedia.org/wiki/Vorticity_equation):

$$\frac{D\boldsymbol{\omega}}{Dt} = \boldsymbol{\omega}\cdot\nabla \boldsymbol{u}+\frac{1}{Re}\nabla^2\boldsymbol{\omega}$$

In 2 dimensions, the vorticity can be described by a single scalar, \\(\boldsymbol{\omega} = \nabla\times \boldsymbol{u}=\omega \boldsymbol{k}\\). It is easy to verify that the vortex stretching term ( \\(\boldsymbol{\omega}\cdot\nabla\boldsymbol{u}\\) ) is identically zero, so that we have an advection diffusion equation for the (scalar) vorticity \\(\omega\\):

$$\frac{D\omega}{Dt} = \frac{1}{Re}\nabla^2\omega$$

Finally, it is well known that in two dimensions a stream function \\(\psi\\) is guaranteed to exist and is related to the vorticity by

$$\omega = -\nabla^2 \psi$$

$$u = -\frac{\partial \psi}{\partial y}\quad v = \frac{\partial \psi}{\partial x}$$

So that the continuity equation is satisfied identically and we need not worry about it. The numerical algorithm is now clear, at each time-step, given a vorticity field \\(\omega\\), we can solve for the streamfunction using a Poisson solver, and given the streamfunction we can compute the velocities required to update the vorticity. 

The FFT and the pseudospectral method
======
At the heart of the solution strategy we find the [Fast Fourier Transform](https://en.wikipedia.org/wiki/Fast_Fourier_transform), developed by Cooley and Tukey. Without going into too much detail, the FFT allows us to compute Fourier transforms and their inverses in \\(\mathcal{O}(N\log N)\\) time (fast). This is useful as evaluating derivatives in Fourier space is easy (multiplication by \\(i^k\\), where \\(k\\) is the derivative order. Thus, we can evaluate a derivative by moving into Fourier space, multiplication by \\(i^k\\) and inverting back to real space. 

The pseudospectral method makes full use of this, by evaluating derivatives in Fourier space and nonlinearities in real space. Furthermore, for the Poisson solver, it is well known that the FFT diagonalizes the Laplacian, so that we can really easily solve the Poisson equation.

The main drawback of our method (and spectral methods more generally) is that enforcing boundary conditions can be very hard, and generally involves methods [beyond the scope](https://ntrs.nasa.gov/api/citations/19960029104/downloads/19960029104.pdf) of this post. For this reason we choose to work with the natural conditions for spectral methods: periodic boundary conditions. 

Given a spectral discretization in space, the system can be sovled in time using a built-in integrator, like scipy's solve_ivp. A good method to prescribe is [RK45](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods), although in some instances a stiff solver like the Backward Differentiation Formula, [BDF](https://en.wikipedia.org/wiki/Backward_differentiation_formula) might be better suited.

Some videos of Solutions
------
We simulate the system for different Reynolds numbers and different spectral accuracies. My laptop can sovle the equations for N = 64 in less than a minute, but struggles a bit more for higher N. A collection of movies is available in [this Youtube playlist](https://www.youtube.com/playlist?list=PLLYwsGNINCFRoQVrj4AtngQL8uvaAS3RL).

Here for Re = 1000, we don't observe any fancy formations and the solutions seem to converge to a steady pattern. 

<video width="800" height="800" controls>
  <source src="/videos/spectral_vortcity/vorticity_evolution_Re_1000_N_64_HD_test.mp4" type="video/mp4">
Your browser does not support the video tag.
</video>

For Re = 1 million we observe far more interesting unsteady phenomena.

<video width="800" height="800" controls>
  <source src="/videos/spectral_vortcity/vorticity_evolution_Re_1000000_N_256.mp4" type="video/mp4">
Your browser does not support the video tag.
</video>

Using a cluster, we can obtain solutions with N=512

<video width="800" height="800" controls>
  <source src="/videos/spectral_vortcity/vorticity_evolution_Re_10000000_N_512_HD.mp4" type="video/mp4">
Your browser does not support the video tag.
</video>

The code required to solve the problem, including a smoothed random initial condition is available here 

```
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft2, ifft2, fftfreq
from scipy.integrate import solve_ivp
import matplotlib.animation as animation
import time

print('Parameters')
# Parameters
N = 64  # Grid size 
L = 1.0  # Domain size
dx = L / N
dy = L / N
Re = 1e5  # Reynolds number
T = 100  # Total time
time_steps = 101  # Number of time steps
t = np.linspace(0,T,time_steps); t_span = (0,T)
dt = t[1]-t[0]
# Grid points
x = np.linspace(0, L, N, endpoint=False)
y = np.linspace(0, L, N, endpoint=False)
dx=x[1]-x[0];dy=y[1]-y[0]
X, Y = np.meshgrid(x, y)

# Frequency components
kx = fftfreq(N, d=dx) * 2 * np.pi
ky = fftfreq(N, d=dy) * 2 * np.pi
KX, KY = np.meshgrid(kx, ky)
K2 = KX**2 + KY**2
K2[0, 0] = 1  # To avoid division by zero later

def RHS(t, vorticity_vector):
    w = vorticity_vector.reshape((N,N), order='C')
    w_hat = fft2(w)
    
    psi_hat = -w_hat / K2
    psi = np.real(ifft2(psi_hat))

    # Compute velocity: u = curl(psi k)
    u = np.real(ifft2(1j * KY * psi_hat))  # u = d(psi)/dy
    v = -np.real(ifft2(1j * KX * psi_hat))  # v = -d(psi)/dx

    # Compute nonlinear term: (u.grad)w
    w_x = np.real(ifft2(1j * KX * w_hat))
    w_y = np.real(ifft2(1j * KY * w_hat))
    nonlinear_term = u * w_x + v * w_y
    w_x = np.real(ifft2(1j * KX * w_hat))  # ∂w/∂x
    w_y = np.real(ifft2(1j * KY * w_hat))  # ∂w/∂y
    nonlinear_term = u * w_x + v * w_y

    # Compute the diffusion term
    diffusion_term = np.real(ifft2(-K2 * w_hat)) / Re
    dwdt = diffusion_term-nonlinear_term

    # Explicit update of vorticity
    return dwdt.flatten(order='C')

# Initial conditions for vorticity w (random smooth initial condition)
random_coefficients = (np.random.normal(size=(N, N)) + 1j * np.random.normal(size=(N, N)))
gaussian_filter = np.exp(-0.01 * (KX**2 + KY**2))
w_hat = random_coefficients * gaussian_filter
w_hat = (w_hat + np.conj(np.flipud(np.fliplr(w_hat)))) / 2
w_hat[0, 0] = 0
w_ic = np.real(ifft2(w_hat))+np.sin(10*X)*np.sin(4*Y)
w = w_ic/np.max(abs(w_ic))
vorticity_vector_ic = w.flatten(order='C')
w_sol_vec = np.zeros((N**2, len(t)))
w_sol_vec[:,0] = vorticity_vector_ic
#solve
start = time.time();print('Start')
for j in range(len(t)-1):
    t_span_j = (t[j],t[j+1])
    solution = solve_ivp(RHS, t_span_j, w_sol_vec[:,j],method='RK45', atol=1e-5)
    w_sol_vec[:,j+1] = solution.y[:,-1]
    print(f"{round(100*(1+j)/len(t),3)} %")
w_sol = w_sol_vec.reshape((N,N,len(t)),order='C')
stop = time.time()
print(f'Equations solved in {round(stop-start,1)} seconds')
```
