---
title: 'Two dimensional fluid mechanics in your laptop!'
date: 2024-08-12
permalink: /posts/2012/08/vorticity_eq/
tags:
  - cool posts
  - category1
  - category2
---

Theory
======

In this short post I provide some tricks and code to simulate high accuracy two dimensional flows using a laptop. Starting from the Navier-Stokes equations:
$$\frac{\partial \boldsymbol{u}}{\partial t}+\boldsymbol{u}\cdot\nabla\boldsymbol{u}=-\nabla p +\frac{1}{\text{Re}}\nabla^2\boldsymbol{u}$$
The pressure gradient can be problematic in CFD algorithms, as there are no explicit conditions for the pressure at boundaries. Hence, a commonly used trick is to take the curl of the momentum equations, leading to the vorticity equation:
$$\frac{D\omega}{Dt} = \omega\cdot\nabla \boldsymbol{u}+\frac{1}{Re}\nabla^2\omega$$
In 2 dimensions, the vortex stretching term is identically zero, so that we have an advection diffusion equation for the (scalar) vortcity $\omega$:
$$\frac{D\omega}{Dt} = \frac{1}{Re}\nabla^2\omega$$
Finally, it is well known that in two dimensions a stream function $\psi$ is guaranteed to exist and is related to the vorticity by
$$\omega = -\nabla^2 \psi$$
$$u = -\frac{\partial \psi}{\partial y}\quad v = \frac{\partial \psi}{\partial x}$$

The FFT and the pseudospectral method
======
At the heart of the solution strategy we find the Fast Fourier Transform
Aren't headings cool?
------