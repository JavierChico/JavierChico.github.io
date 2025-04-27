---
title: 'Stochastic Simulations for Dummies: Why √(Δt) is the Key'
date: 2025-04-29
permalink: /posts/2025/04/stochastic_simulation/
tags:
  - Stochastic differential equations
  - Mathematical Biology
  - Brownian Motion
---

Deterministic motion
======

# Stochastic Simulations for Dummies: Why Random Motion Needs √(Δt)

Let's talk about randomness. Not just any randomness—the kind that makes tiny pollen grains jiggle in water (Brownian motion) or stock prices wiggle unpredictably.  

## 1. Deterministic vs. Random Motion

### The Predictable World

Suppose you are moving around with velocity $v$. How much will you have traveled after a small time $dt$? The answer is clearly $dx= v dt$, provided that the time increment is small enough and the velocity does not change much (is continuous). This makes sense: double the time, double the distance. Today, we will see how randomness breaks this very intuitive scaling - and of course how to remedy this.

### The Random World

Now suppose you take random steps—like a drunk person staggering left and right. Each step is independent, and on average, you go nowhere. But you don’t stay perfectly still—you spread out over time. By this I mean if you do this with 100 friends, you expect people to be more spread out over time.

**Key Idea:**  
- The **average position** stays zero  
- The **spread** (variance) grows with time




