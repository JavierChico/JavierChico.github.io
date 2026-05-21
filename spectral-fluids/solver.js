/**
 * solver.js - Pseudospectral 2D Navier-Stokes solver with periodic boundary conditions.
 * Mirrors the exact physics and equations from spectral_vorticity_integrator.py.
 */

import { FFT2D } from './fft.js';

export class FluidSolver {
  constructor(N = 128, Re = 1e4) {
    this.N = N;
    this.Re = Re;
    this.L = 1.0;
    this.initBuffers();
  }

  /**
   * Allocates all grid-specific buffers and precomputes spectral wavenumbers.
   */
  initBuffers() {
    const N = this.N;
    const size = N * N;

    this.fft = new FFT2D(N);

    // 1. Physical fields
    this.w = new Float64Array(size);         // Vorticity field
    this.w_temp = new Float64Array(size);    // Intermediate RK2 vorticity
    this.u = new Float64Array(size);         // Velocity X
    this.v = new Float64Array(size);         // Velocity Y
    this.rhs = new Float64Array(size);       // Primary RK2 RHS
    this.rhs_temp = new Float64Array(size);  // Secondary RK2 RHS
    this.psi = new Float64Array(size);       // Streamfunction (physical)
    this.velMag = new Float64Array(size);    // Preallocated velocity magnitude

    // 2. Fourier space coefficients (real/imaginary pairs)
    this.w_re = new Float64Array(size);
    this.w_im = new Float64Array(size);
    
    this.psi_re = new Float64Array(size);
    this.psi_im = new Float64Array(size);

    this.u_re = new Float64Array(size);
    this.u_im = new Float64Array(size);

    this.v_re = new Float64Array(size);
    this.v_im = new Float64Array(size);

    this.wx_re = new Float64Array(size);
    this.wx_im = new Float64Array(size);

    this.wy_re = new Float64Array(size);
    this.wy_im = new Float64Array(size);

    this.diff_re = new Float64Array(size);
    this.diff_im = new Float64Array(size);

    // 3. Precompute Wavenumbers (kx, ky) and Laplacian eigenvalues (K2)
    // dx = dy = L / N
    const dx = this.L / N;
    this.kx = new Float64Array(size);
    this.ky = new Float64Array(size);
    this.K2 = new Float64Array(size);

    for (let r = 0; r < N; r++) {
      // Ky wavenumber (columns)
      const ky_val = 2 * Math.PI * (r < N / 2 ? r : r - N);
      for (let c = 0; c < N; c++) {
        // Kx wavenumber (rows)
        const kx_val = 2 * Math.PI * (c < N / 2 ? c : c - N);
        const idx = r * N + c;
        this.kx[idx] = kx_val;
        this.ky[idx] = ky_val;
        
        let k2 = kx_val * kx_val + ky_val * ky_val;
        if (r === 0 && c === 0) k2 = 1.0; // Avoid division by zero
        this.K2[idx] = k2;
      }
    }
  }

  /**
   * Resizes the grid to a new resolution N, bilinearly interpolating the current state.
   * @param {number} newN - New grid size (must be power of 2).
   */
  resize(newN) {
    if (newN === this.N) return;
    
    console.log(`Resizing solver grid from ${this.N} to ${newN}...`);
    const oldW = this.w;
    const oldN = this.N;

    this.N = newN;
    this.initBuffers();

    // Bilinearly interpolate the existing vorticity field to the new grid size
    const scale = oldN / newN;
    for (let cy = 0; cy < newN; cy++) {
      const oy = cy * scale;
      const y0 = Math.floor(oy);
      const y1 = (y0 + 1) % oldN;
      const dy = oy - y0;

      for (let cx = 0; cx < newN; cx++) {
        const ox = cx * scale;
        const x0 = Math.floor(ox);
        const x1 = (x0 + 1) % oldN;
        const dx = ox - x0;

        const w00 = oldW[y0 * oldN + x0];
        const w10 = oldW[y0 * oldN + x1];
        const w01 = oldW[y1 * oldN + x0];
        const w11 = oldW[y1 * oldN + x1];

        const interpolatedVorticity = 
          (1 - dx) * (1 - dy) * w00 +
          dx * (1 - dy) * w10 +
          (1 - dx) * dy * w01 +
          dx * dy * w11;

        this.w[cy * newN + cx] = interpolatedVorticity;
      }
    }
  }

  /**
   * Sets the Reynolds number.
   * @param {number} Re - Reynolds number.
   */
  setReynolds(Re) {
    this.Re = Re;
  }

  /**
   * Computes the Right-Hand Side of the vorticity equation:
   * RHS = D(w) - [u*w_x + v*w_y]
   * @param {Float64Array} w_in - Input physical vorticity field.
   * @param {Float64Array} rhs_out - Output array to write the physical RHS into.
   */
  computeRHS(w_in, rhs_out) {
    const N = this.N;
    const size = N * N;

    // 1. Copy input physical vorticity into our Fourier real buffer, zeroing imaginary part
    for (let i = 0; i < size; i++) {
      this.w_re[i] = w_in[i];
      this.w_im[i] = 0.0;
    }

    // 2. Transform vorticity to spectral space
    this.fft.forward(this.w_re, this.w_im);

    // 3. Solve Poisson for streamfunction: psi_hat = -w_hat / K^2
    for (let i = 0; i < size; i++) {
      this.psi_re[i] = -this.w_re[i] / this.K2[i];
      this.psi_im[i] = -this.w_im[i] / this.K2[i];
    }
    // Zero out the DC component (mean flow streamfunction)
    this.psi_re[0] = 0.0;
    this.psi_im[0] = 0.0;

    // 4. Compute velocities in Fourier space:
    // u_hat = i * ky * psi_hat  =>  u_re = -ky * psi_im,  u_im = ky * psi_re
    // v_hat = -i * kx * psi_hat =>  v_re = kx * psi_im,   v_im = -kx * psi_re
    for (let i = 0; i < size; i++) {
      const kx = this.kx[i];
      const ky = this.ky[i];
      
      this.u_re[i] = -ky * this.psi_im[i];
      this.u_im[i] = ky * this.psi_re[i];
      
      this.v_re[i] = kx * this.psi_im[i];
      this.v_im[i] = -kx * this.psi_re[i];
      
      // 5. Compute vorticity gradients in Fourier space:
      // wx_hat = i * kx * w_hat => wx_re = -kx * w_im, wx_im = kx * w_re
      // wy_hat = i * ky * w_hat => wy_re = -ky * w_im, wy_im = ky * w_re
      this.wx_re[i] = -kx * this.w_im[i];
      this.wx_im[i] = kx * this.w_re[i];
      
      this.wy_re[i] = -ky * this.w_im[i];
      this.wy_im[i] = ky * this.w_re[i];

      // 6. Compute Diffusion in Fourier space:
      // diff_hat = -K^2 * w_hat / Re
      this.diff_re[i] = -this.K2[i] * this.w_re[i] / this.Re;
      this.diff_im[i] = -this.K2[i] * this.w_im[i] / this.Re;
    }

    // 7. Transform velocities, gradients, and diffusion back to physical space
    // NOTE: Streamfunction (psi_re) IFFT is skipped here and handled on-demand to save huge processing time!
    this.fft.inverse(this.u_re, this.u_im);
    this.fft.inverse(this.v_re, this.v_im);
    this.fft.inverse(this.wx_re, this.wx_im);
    this.fft.inverse(this.wy_re, this.wy_im);
    this.fft.inverse(this.diff_re, this.diff_im);

    // 8. Compute total RHS in physical space: RHS = Diffusion - Advection
    // Advection = u * w_x + v * w_y
    for (let i = 0; i < size; i++) {
      this.u[i] = this.u_re[i];     // Store velocities
      this.v[i] = this.v_re[i];

      const advection = this.u_re[i] * this.wx_re[i] + this.v_re[i] * this.wy_re[i];
      const diffusion = this.diff_re[i];
      rhs_out[i] = diffusion - advection;
    }
  }

  /**
   * Steps the vorticity field forward in time using Runge-Kutta 2 (Heun's Method).
   * @param {number} dt - Timestep.
   */
  step(dt) {
    const size = this.N * this.N;

    // 1. Compute RHS at current state: rhs1 = RHS(w)
    this.computeRHS(this.w, this.rhs);

    // 2. Predict intermediate state: w_temp = w + dt * rhs1
    for (let i = 0; i < size; i++) {
      this.w_temp[i] = this.w[i] + dt * this.rhs[i];
    }

    // 3. Compute RHS at intermediate state: rhs2 = RHS(w_temp)
    this.computeRHS(this.w_temp, this.rhs_temp);

    // 4. Correct final state: w = w + 0.5 * dt * (rhs1 + rhs2)
    for (let i = 0; i < size; i++) {
      this.w[i] = this.w[i] + 0.5 * dt * (this.rhs[i] + this.rhs_temp[i]);
    }
  }

  /**
   * Steps the vorticity field forward in time using a high-performance 1st-order Euler step.
   * Only requires 1 RHS calculation per step, doubling physics performance!
   * @param {number} dt - Timestep.
   */
  stepEuler(dt) {
    const size = this.N * this.N;

    // 1. Compute RHS at current state: rhs = RHS(w)
    this.computeRHS(this.w, this.rhs);

    // 2. Advance state: w = w + dt * rhs
    for (let i = 0; i < size; i++) {
      this.w[i] = this.w[i] + dt * this.rhs[i];
    }
  }

  /**
   * Computes the physical streamfunction field psi on-demand from the current physical vorticity.
   * Only needs to be executed once per frame, and only when streamfunction contours or heightmaps are rendered.
   */
  updatePhysicalPsi() {
    const size = this.N * this.N;

    // 1. Copy current vorticity to Fourier buffers
    for (let i = 0; i < size; i++) {
      this.w_re[i] = this.w[i];
      this.w_im[i] = 0.0;
    }

    // 2. Transform to spectral space
    this.fft.forward(this.w_re, this.w_im);

    // 3. Solve Poisson for streamfunction: psi_hat = -w_hat / K^2
    for (let i = 0; i < size; i++) {
      this.psi_re[i] = -this.w_re[i] / this.K2[i];
      this.psi_im[i] = -this.w_im[i] / this.K2[i];
    }
    this.psi_re[0] = 0.0;
    this.psi_im[0] = 0.0;

    // 4. Transform back to physical space
    this.fft.inverse(this.psi_re, this.psi_im);

    // 5. Store in the physical streamfunction array
    for (let i = 0; i < size; i++) {
      this.psi[i] = this.psi_re[i];
    }
  }

  /**
   * Computes the velocity magnitude in-place into the preallocated velMag buffer.
   * Eliminates the need to allocate Float64Arrays dynamically during rendering.
   */
  updateVelocityMagnitude() {
    const size = this.N * this.N;
    for (let i = 0; i < size; i++) {
      this.velMag[i] = Math.sqrt(this.u[i] * this.u[i] + this.v[i] * this.v[i]);
    }
  }

  /**
   * Automatically calculates a stable timestep based on the CFL condition.
   * dt = C * dx / max(|u|, |v|)
   * @param {number} maxCFL - Safety Courant number (typically 0.1 - 0.3 for explicit schemes).
   * @returns {number} Stable timestep.
   */
  calculateStableTimestep(maxCFL = 0.2) {
    const N = this.N;
    const dx = this.L / N;
    let maxVelSq = 1e-8; // Avoid division by zero

    for (let i = 0; i < N * N; i++) {
      const u = this.u[i];
      const v = this.v[i];
      const velSq = u * u + v * v;
      if (velSq > maxVelSq) {
        maxVelSq = velSq;
      }
    }

    const maxVel = Math.sqrt(maxVelSq);
    // Limit dt to be between 1e-4 and 5e-3 to keep simulation running nicely
    const dt = (maxCFL * dx) / maxVel;
    return Math.max(1e-4, Math.min(5e-3, dt));
  }

  /**
   * Resets the vorticity field to zero (clear canvas).
   */
  clear() {
    const size = this.N * this.N;
    this.w.fill(0.0);
    this.u.fill(0.0);
    this.v.fill(0.0);
    this.psi.fill(0.0);
    this.velMag.fill(0.0);
  }

  /**
   * Injects a localized Gaussian packet of vorticity.
   * Handles periodic boundary conditions naturally.
   * @param {number} cx - Normalized coordinate X [0, 1].
   * @param {number} cy - Normalized coordinate Y [0, 1].
   * @param {number} radius - Vortex core radius (normalized scale).
   * @param {number} strength - Vorticity amplitude.
   */
  injectVortex(cx, cy, radius, strength) {
    const N = this.N;
    const size = N * N;

    for (let r = 0; r < N; r++) {
      let dy = (r / N) - cy;
      // Periodic wraps in y
      if (dy > 0.5) dy -= 1.0;
      if (dy < -0.5) dy += 1.0;

      for (let c = 0; c < N; c++) {
        let dx = (c / N) - cx;
        // Periodic wraps in x
        if (dx > 0.5) dx -= 1.0;
        if (dx < -0.5) dx += 1.0;

        const r2 = dx * dx + dy * dy;
        const rLimit = radius * radius;
        
        // Compute Gaussian factor up to 3 sigma
        if (r2 < rLimit * 9) {
          const factor = Math.exp(-r2 / (2 * rLimit));
          const idx = r * N + c;
          this.w[idx] += strength * factor;
        }
      }
    }
  }

  /**
   * Injects an elegant vortex dipole perpendicular to the drag direction.
   * Simulates fluid stirring by creating self-propelling vortex pairs.
   * @param {number} cx - Normalized start coordinate X.
   * @param {number} cy - Normalized start coordinate Y.
   * @param {number} dx - Drag delta X.
   * @param {number} dy - Drag delta Y.
   * @param {number} radius - Dipole core size.
   * @param {number} strength - Stirring multiplier.
   */
  injectDipole(cx, cy, dx, dy, radius, strength) {
    const len = Math.sqrt(dx * dx + dy * dy);
    if (len < 1e-6) return;

    // Perpendicular vector for the dipole separation
    const px = -dy / len;
    const py = dx / len;

    // Spacing between the vortex centers
    const spacing = radius * 0.7;

    const posCx = cx + px * spacing;
    const posCy = cy + py * spacing;
    const negCx = cx - px * spacing;
    const negCy = cy - py * spacing;

    // Injected vorticity is proportional to drag speed
    const amplitude = strength * len * 15.0;

    // Inject positive and negative pair
    this.injectVortex(posCx, posCy, radius, amplitude);
    this.injectVortex(negCx, negCy, radius, -amplitude);
  }

  /**
   * Injects random small-scale vorticity noise to kickstart turbulence.
   */
  injectNoise(density = 0.05) {
    const N = this.N;
    for (let i = 0; i < N * N; i++) {
      if (Math.random() < density) {
        this.w[i] += (Math.random() - 0.5) * 5.0;
      }
    }
  }

  /**
   * Calculates enstrophy (vorticity variance: \Omega = 0.5 * \int w^2 dx dy)
   * Normalized by grid cell count.
   */
  getEnstrophy() {
    const size = this.N * this.N;
    let sum = 0.0;
    for (let i = 0; i < size; i++) {
      sum += this.w[i] * this.w[i];
    }
    return 0.5 * sum / size;
  }

  /**
   * Calculates average kinetic energy (E = 0.5 * \int (u^2 + v^2) dx dy)
   * Normalized by grid cell count.
   */
  getKineticEnergy() {
    const size = this.N * this.N;
    let sum = 0.0;
    for (let i = 0; i < size; i++) {
      sum += this.u[i] * this.u[i] + this.v[i] * this.v[i];
    }
    return 0.5 * sum / size;
  }
}
