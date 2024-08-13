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

Some videos of Solutions
------
We simulate the system for different Reynolds numbers and different spectral accuracies. My laptop can sovle the equations for N = 64 in less than a minute, but struggles a bit more for higher N. A collection of movies is available in [this Youtube playlist](https://www.youtube.com/playlist?list=PLLYwsGNINCFRoQVrj4AtngQL8uvaAS3RL).

Here for Re = 1000, we don't observe any fancy formations and the solutions seem to converge to a steady pattern. 

For Re = 1 million we observe far more interesting unsteady phenomena.

<video width="800" height="800" controls>
  <source src="/videos/spectral_vortcity/vorticity_evolution_Re_1000000_N_256.mp4" type="video/mp4">
Your browser does not support the video tag.
</video>
