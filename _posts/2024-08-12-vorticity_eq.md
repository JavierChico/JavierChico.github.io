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
$$\frac{\partial \bold{u}}P{\partial t}$$
The pressure gradient can be problematic in CFD algorithms, as there are no explicit conditions for the pressure at boundaries. Does a commonly used trick is to take the curl of the momentum equations, leading to the vorticity equation:

In 2 dimensions, the vortex stretching term is identically zero, so that we have an advection diffusion equation for the (scalar) vortcity $\omega$. Finally, it is well known that in two dimensions a stream function $\psi$ is guaranteed to exist and is related to the vorticity by
$$\omega = -\nabla^2 \psi$$


You can have many headings
======

Aren't headings cool?
------
