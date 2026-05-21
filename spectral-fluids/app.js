/**
 * app.js - Main Application Controller for SpectralFluids2D.
 * Connects the UI controls, manages the simulation loop, handles mouse interaction,
 * and renders the fluid vorticity and particle trails.
 */

import { FluidSolver } from './solver.js';

// --- APPLICATION STATE ---
let solver = null;
let isPlaying = true;
let currentN = 128;
let currentPalette = 'coolwarm';
let currentRenderMode = 'vorticity';
let enableStreamlines = true;
let enableContours = false;
let stepRequested = false;

// Optimization and Performance States
let currentIntegrationScheme = 'rk2';
let currentCfl = 0.20;

// Brush settings
let brushType = 'dipole';
let brushRadius = 0.04;
let brushStrength = 1.0;

// Performance diagnostics
let lastFrameTime = performance.now();
let fpsInterval = 1000;
let lastFpsUpdate = performance.now();
let framesThisSecond = 0;
let averageSolverTime = 0.0;
let solverTimeAccumulator = 0;
let solverStepsCount = 0;

// Canvas details
const canvas = document.getElementById('fluid-canvas');
const ctx = canvas.getContext('2d');
let canvasWidth = canvas.width;
let canvasHeight = canvas.height;

// Helper canvas for hardware-accelerated bilinear upscaling
const helperCanvas = document.createElement('canvas');
const helperCtx = helperCanvas.getContext('2d');

// Mouse dragging state
let isDragging = false;
let lastMouseX = 0;
let lastMouseY = 0;
let mouseStartInterval = 0;

// Particle System for Streamlines
const MAX_PARTICLE_COUNT = 2000;
let currentParticleCount = 1200;
const TRAIL_LENGTH = 6; // Historical positions to keep for trailing
class FlowParticle {
  constructor() {
    this.reset();
  }

  reset() {
    this.x = Math.random();
    this.y = Math.random();
    this.history = [];
    for (let i = 0; i < TRAIL_LENGTH; i++) {
      this.history.push({ x: this.x, y: this.y });
    }
    this.life = Math.random() * 150 + 50; // Random lifetime
    this.speedFactor = Math.random() * 0.4 + 0.8;
  }
}
const particles = Array.from({ length: MAX_PARTICLE_COUNT }, () => new FlowParticle());

// --- INITIALIZATION ---
function init() {
  // 1. Create the physics solver
  solver = new FluidSolver(currentN, getReynoldsFromSlider());
  
  // Set up initial conditions with some random nice vortices
  resetToNoise();

  // 2. Setup DOM Event Listeners
  setupUIEventListeners();

  // 3. Initialize the particle system
  resizeHelperCanvas();

  // 4. Start the simulation loop
  requestAnimationFrame(loop);
}

// --- CORE SIMULATION & RENDERING LOOP ---
function loop(timestamp) {
  // A. Calculate Frame FPS
  const delta = timestamp - lastFrameTime;
  lastFrameTime = timestamp;
  
  framesThisSecond++;
  if (timestamp - lastFpsUpdate >= fpsInterval) {
    const fps = Math.round((framesThisSecond * 1000) / (timestamp - lastFpsUpdate));
    document.getElementById('stat-fps').textContent = fps;
    framesThisSecond = 0;
    lastFpsUpdate = timestamp;
    
    // Update average solver step time
    if (solverStepsCount > 0) {
      averageSolverTime = solverTimeAccumulator / solverStepsCount;
      document.getElementById('stat-time').textContent = `${averageSolverTime.toFixed(1)} ms`;
      solverTimeAccumulator = 0;
      solverStepsCount = 0;
    } else {
      document.getElementById('stat-time').textContent = `0.0 ms`;
    }
  }

  // B. Physics Step (if playing)
  let currentDt = 0.002;
  
  if (isPlaying) {
    const t0 = performance.now();
    
    // Calculate adaptive time-step based on CFL condition for stability
    currentDt = solver.calculateStableTimestep(currentCfl);
    
    // Step equations forward using selected scheme
    if (currentIntegrationScheme === 'euler') {
      solver.stepEuler(currentDt);
    } else {
      solver.step(currentDt);
    }
    
    const t1 = performance.now();
    solverTimeAccumulator += (t1 - t0);
    solverStepsCount++;

    // C. Update diagnostics in the UI
    updateDiagnostics();
  }

  // D. Flow Particles Update (Active only when playing or stepping)
  if (enableStreamlines && (isPlaying || stepRequested)) {
    if (stepRequested) {
      currentDt = solver.calculateStableTimestep(currentCfl);
    }
    updateParticles(currentDt);
    if (!isPlaying) {
      stepRequested = false; // Reset stepping flag after one update
    }
  }

  // E. Update physical Streamfunction only on-demand when actively rendering
  if (currentRenderMode === 'streamfunction' || enableContours) {
    solver.updatePhysicalPsi();
  }

  // F. Render Fluid Field
  renderFluid();

  // G. Render Streamline Contours and Particles on top
  if (enableStreamlines) {
    renderParticleTrails();
  }

  // Queue next frame
  requestAnimationFrame(loop);
}

// --- PHYSICS INJECTION & INTERACTIONS ---
function handleCanvasInteraction(clientX, clientY, isMoving) {
  const rect = canvas.getBoundingClientRect();
  const cx = (clientX - rect.left) / rect.width;
  const cy = (clientY - rect.top) / rect.height;

  // Clamp within bounds [0, 1]
  const clampedX = Math.max(0, Math.min(1, cx));
  const clampedY = Math.max(0, Math.min(1, cy));

  if (!isMoving) {
    // Initial click/touch - inject a single vortex
    const strength = brushStrength * (brushType === 'negative' ? -3.0 : 3.0);
    if (brushType === 'dipole') {
      // Small default dipole on tap
      solver.injectDipole(clampedX, clampedY, 0.01, 0.01, brushRadius, brushStrength);
    } else if (brushType === 'noise') {
      solver.injectNoise(0.1);
    } else {
      solver.injectVortex(clampedX, clampedY, brushRadius, strength);
    }
  } else {
    // Dragging - calculate velocity
    const dx = clampedX - lastMouseX;
    const dy = clampedY - lastMouseY;

    if (brushType === 'dipole') {
      solver.injectDipole(clampedX, clampedY, dx, dy, brushRadius, brushStrength);
    } else if (brushType === 'positive') {
      solver.injectVortex(clampedX, clampedY, brushRadius, brushStrength * 0.4);
    } else if (brushType === 'negative') {
      solver.injectVortex(clampedX, clampedY, brushRadius, -brushStrength * 0.4);
    } else if (brushType === 'noise') {
      solver.injectNoise(0.01 * brushStrength);
    }
  }

  lastMouseX = clampedX;
  lastMouseY = clampedY;
}

// --- PARTICLE DYNAMIC VECTOR SYSTEM ---
function updateParticles(dt) {
  const N = solver.N;
  const u = solver.u;
  const v = solver.v;
  
  // Visual tuning multiplier to make trails look elegant and active
  const visualMultiplier = 2.0;

  for (let i = 0; i < currentParticleCount; i++) {
    const p = particles[i];
    p.life--;

    if (p.life <= 0) {
      p.reset();
      continue;
    }

    // 1. Bilinearly interpolate fluid velocity at the particle's fractional position (p.x, p.y)
    const px = p.x * N;
    const py = p.y * N;
    
    const x0 = Math.floor(px) % N;
    const x1 = (x0 + 1) % N;
    const y0 = Math.floor(py) % N;
    const y1 = (y0 + 1) % N;
    
    const tx = px - Math.floor(px);
    const ty = py - Math.floor(py);

    // Fetch velocity grid values
    const u00 = u[y0 * N + x0];
    const u10 = u[y0 * N + x1];
    const u01 = u[y1 * N + x0];
    const u11 = u[y1 * N + x1];

    const v00 = v[y0 * N + x0];
    const v10 = v[y0 * N + x1];
    const v01 = v[y1 * N + x0];
    const v11 = v[y1 * N + x1];

    // Bilinear velocity values
    const interpU = (1 - tx) * (1 - ty) * u00 + tx * (1 - ty) * u10 + (1 - tx) * ty * u01 + tx * ty * u11;
    const interpV = (1 - tx) * (1 - ty) * v00 + tx * (1 - ty) * v10 + (1 - tx) * ty * v01 + tx * ty * v11;

    // 2. Shift history
    for (let h = TRAIL_LENGTH - 1; h > 0; h--) {
      p.history[h].x = p.history[h - 1].x;
      p.history[h].y = p.history[h - 1].y;
    }
    p.history[0].x = p.x;
    p.history[0].y = p.y;

    // 3. Move particle with periodic wrapping, scaled by the solver's physical timestep
    p.x = (p.x + interpU * dt * visualMultiplier * p.speedFactor + 1.0) % 1.0;
    p.y = (p.y + interpV * dt * visualMultiplier * p.speedFactor + 1.0) % 1.0;
  }
}

// --- OPTIMIZED IMAGE RENDERING PIPELINE ---
function resizeHelperCanvas() {
  helperCanvas.width = currentN;
  helperCanvas.height = currentN;
}

function renderFluid() {
  const N = solver.N;
  const size = N * N;
  
  // A. Create/Get high-performance ImageData array on helper canvas size
  const imgData = helperCtx.createImageData(N, N);
  const data = imgData.data;

  // Select target field to render
  let renderField = solver.w; // default vorticity
  if (currentRenderMode === 'streamfunction') {
    renderField = solver.psi;
  } else if (currentRenderMode === 'velocity') {
    solver.updateVelocityMagnitude();
    renderField = solver.velMag;
  }

  // Find max amplitude for normalization to map colors cleanly
  let maxAbs = 1e-4;
  for (let i = 0; i < size; i++) {
    const abs = Math.abs(renderField[i]);
    if (abs > maxAbs) maxAbs = abs;
  }

  // Pre-calculate streamfunction properties for contours
  let maxPsi = 1e-4;
  if (enableContours) {
    for (let i = 0; i < size; i++) {
      const abs = Math.abs(solver.psi[i]);
      if (abs > maxPsi) maxPsi = abs;
    }
  }

  // Pre-calculate inverse factors to avoid divisions in pixel loop (speeds up CPU significantly)
  const invMaxAbs = 1.0 / maxAbs;
  const invMaxPsi = enableContours ? 1.0 / maxPsi : 1.0;

  // B. Write pixel color mappings in-place
  for (let y = 0; y < N; y++) {
    for (let x = 0; x < N; x++) {
      const idx = y * N + x;
      const pixelOffset = idx * 4;

      const val = renderField[idx];
      const norm = val * invMaxAbs; // scaled [-1.0, 1.0] or [0.0, 1.0]

      // Colors
      let r = 0, g = 0, b = 0;

      // Apply selected Color Palette
      if (currentPalette === 'neon') {
        // Neon Vortex: positive (magenta/orange), negative (cyan/blue), zero (charcoal dark)
        if (norm >= 0) {
          // positive vorticity - Magenta
          r = Math.round(14 + norm * 241); // 14 -> 255
          g = Math.round(14 + norm * 0);   // 14 -> 14
          b = Math.round(20 + norm * 160); // 20 -> 180
        } else {
          // negative vorticity - Cyan
          const mag = -norm;
          r = Math.round(14 + mag * 0);
          g = Math.round(14 + mag * 210);
          b = Math.round(20 + mag * 235);
        }
      } else if (currentPalette === 'coolwarm') {
        // Traditional scientific blue-white-red
        if (norm >= 0) {
          // positive - red
          r = 255;
          g = Math.round(250 * (1 - norm));
          b = Math.round(250 * (1 - norm));
        } else {
          // negative - blue
          const mag = -norm;
          r = Math.round(250 * (1 - mag));
          g = Math.round(250 * (1 - mag));
          b = 255;
        }
      } else if (currentPalette === 'magma') {
        // Volcanic Magma: black -> purple -> orange -> bright yellow
        const mag = Math.abs(norm); // use magnitude
        r = Math.round(Math.pow(mag, 1.5) * 255);
        g = Math.round(Math.pow(mag, 2.5) * 200);
        b = Math.round(Math.pow(mag, 4.0) * 100 + mag * 40);
      } else {
        // Monochrome slate gray
        const mag = Math.round((Math.abs(norm) * 200) + 15);
        r = mag; g = mag; b = mag;
      }

      // C. Live Streamfunction Contour Line Overlay Hack
      // Draws mathematical streamlines instantly at zero extra compute cost!
      if (enableContours) {
        const psiVal = solver.psi[idx] * invMaxPsi;
        // Multiplication-free rounding modulo check instead of Math.sin (hundreds of times faster!)
        const scaled = (psiVal + 1.0) * 16.0;
        const frac = scaled - Math.round(scaled);
        if (Math.abs(frac) < 0.12) {
          // Overlay bright contour line
          r = Math.round(r * 0.5 + 120);
          g = Math.round(g * 0.5 + 120);
          b = Math.round(b * 0.5 + 120);
        }
      }

      data[pixelOffset] = r;     // Red
      data[pixelOffset + 1] = g; // Green
      data[pixelOffset + 2] = b; // Blue
      data[pixelOffset + 3] = 255; // Alpha (fully opaque)
    }
  }

  // C. Put ImageData on the N x N helper canvas
  helperCtx.putImageData(imgData, 0, 0);

  // D. Upscale helper canvas to main display canvas with hardware bilinear filtering
  ctx.imageSmoothingEnabled = true;
  ctx.drawImage(helperCanvas, 0, 0, canvasWidth, canvasHeight);
}

// Render moving flow particles on top
function renderParticleTrails() {
  ctx.lineWidth = 1.8;
  ctx.lineCap = 'round';

  for (let i = 0; i < currentParticleCount; i++) {
    const p = particles[i];
    
    // Draw each point in the trail history
    for (let h = 0; h < TRAIL_LENGTH - 1; h++) {
      const ptCurrent = p.history[h];
      const ptNext = p.history[h + 1];

      // Skip drawing if the particle wrapped around the periodic boundary
      if (Math.abs(ptCurrent.x - ptNext.x) > 0.5 || Math.abs(ptCurrent.y - ptNext.y) > 0.5) {
        continue;
      }

      const x1 = ptCurrent.x * canvasWidth;
      const y1 = ptCurrent.y * canvasHeight;
      const x2 = ptNext.x * canvasWidth;
      const y2 = ptNext.y * canvasHeight;

      // Taper opacity and size based on trail age (glowing trail effect)
      const opacity = (1.0 - (h / TRAIL_LENGTH)) * 0.45;
      
      if (currentPalette === 'neon') {
        ctx.strokeStyle = `hsla(185, 95%, 75%, ${opacity})`;
      } else if (currentPalette === 'magma') {
        ctx.strokeStyle = `hsla(45, 95%, 65%, ${opacity})`;
      } else if (currentPalette === 'coolwarm') {
        ctx.strokeStyle = `rgba(30, 30, 45, ${opacity * 1.5})`;
      } else {
        ctx.strokeStyle = `rgba(255, 255, 255, ${opacity * 0.8})`;
      }

      ctx.beginPath();
      ctx.moveTo(x1, y1);
      ctx.lineTo(x2, y2);
      ctx.stroke();
    }
  }
}

// --- DIAGNOSTICS & UI UPDATES ---
function updateDiagnostics() {
  const energy = solver.getKineticEnergy();
  const enstrophy = solver.getEnstrophy();

  document.getElementById('stat-energy').textContent = energy.toFixed(4);
  document.getElementById('stat-enstrophy').textContent = enstrophy.toFixed(2);
}

function resetToNoise() {
  solver.clear();
  // Inject some random initial large and mid-scale vortices
  for (let i = 0; i < 8; i++) {
    const cx = Math.random();
    const cy = Math.random();
    const radius = Math.random() * 0.05 + 0.03;
    const strength = (Math.random() - 0.5) * 8.0;
    solver.injectVortex(cx, cy, radius, strength);
  }
  // And some small-scale noise
  solver.injectNoise(0.08);
}

function getReynoldsFromSlider() {
  const sliderVal = parseFloat(document.getElementById('slider-reynolds').value);
  // Logarithmic scaling: 10^sliderVal
  return Math.pow(10, sliderVal);
}

// --- UI INTERACTION EVENT HANDLERS ---
function setupUIEventListeners() {
  // A. Simulation Actions
  const btnPlayPause = document.getElementById('btn-play-pause');
  const txtPlayPause = document.getElementById('txt-play-pause');
  const iconPlay = btnPlayPause.querySelector('.icon-play');
  const iconPause = btnPlayPause.querySelector('.icon-pause');

  btnPlayPause.addEventListener('click', () => {
    isPlaying = !isPlaying;
    if (isPlaying) {
      txtPlayPause.textContent = 'Pause';
      iconPause.classList.remove('hidden');
      iconPlay.classList.add('hidden');
    } else {
      txtPlayPause.textContent = 'Play';
      iconPause.classList.add('hidden');
      iconPlay.classList.remove('hidden');
    }
  });

  document.getElementById('btn-step').addEventListener('click', () => {
    if (!isPlaying) {
      const dt = solver.calculateStableTimestep(currentCfl);
      if (currentIntegrationScheme === 'euler') {
        solver.stepEuler(dt);
      } else {
        solver.step(dt);
      }
      updateDiagnostics();
      stepRequested = true; // Trigger exactly one-frame particle update
    }
  });

  document.getElementById('btn-clear').addEventListener('click', () => {
    solver.clear();
    updateDiagnostics();
  });

  document.getElementById('btn-reset').addEventListener('click', () => {
    resetToNoise();
    updateDiagnostics();
  });

  // B. Physical Controls
  const sliderRe = document.getElementById('slider-reynolds');
  const valRe = document.getElementById('val-reynolds');
  sliderRe.addEventListener('input', () => {
    const Re = getReynoldsFromSlider();
    solver.setReynolds(Re);
    
    // Pretty formatting
    if (Re >= 1e6) {
      valRe.textContent = '1,000,000 (Turbulent)';
    } else if (Re >= 1e3) {
      valRe.textContent = Re.toLocaleString(undefined, { maximumFractionDigits: 0 });
    } else {
      valRe.textContent = Re.toFixed(0);
    }
  });

  const resRadios = document.getElementsByName('resolution');
  resRadios.forEach(radio => {
    radio.addEventListener('change', (e) => {
      const newN = parseInt(e.target.value);
      currentN = newN;
      solver.resize(newN);
      resizeHelperCanvas();
      
      // Reset all drifting particles to ensure they stay on coordinate grid
      particles.forEach(p => p.reset());
    });
  });

  // Integration Scheme Selector Listener
  const integrationRadios = document.getElementsByName('integration');
  integrationRadios.forEach(radio => {
    radio.addEventListener('change', (e) => {
      currentIntegrationScheme = e.target.value;
    });
  });

  // CFL Timestep Selector Listener
  const sliderCfl = document.getElementById('slider-cfl');
  const valCfl = document.getElementById('val-cfl');
  sliderCfl.addEventListener('input', (e) => {
    const val = parseFloat(e.target.value);
    currentCfl = val;
    valCfl.textContent = val.toFixed(2);
  });

  // C. Brush controls
  const selectBrush = document.getElementById('select-brush');
  selectBrush.addEventListener('change', (e) => {
    brushType = e.target.value;
  });

  const sliderBrushSize = document.getElementById('slider-brush-size');
  const valBrushSize = document.getElementById('val-brush-size');
  sliderBrushSize.addEventListener('input', (e) => {
    const val = parseFloat(e.target.value);
    brushRadius = val / 100.0; // convert percentage to normalized scale [0, 1]
    valBrushSize.textContent = `${val}%`;
  });

  const sliderBrushStrength = document.getElementById('slider-brush-strength');
  const valBrushStrength = document.getElementById('val-brush-strength');
  sliderBrushStrength.addEventListener('input', (e) => {
    const val = parseFloat(e.target.value);
    brushStrength = val;
    valBrushStrength.textContent = val.toFixed(1);
  });

  // D. Rendering & Visualization Controls
  const toggleStreamlines = document.getElementById('toggle-streamlines');
  toggleStreamlines.addEventListener('change', (e) => {
    enableStreamlines = e.target.checked;
    if (!enableStreamlines) {
      particles.forEach(p => p.reset());
    }
  });

  // Flow Particle Count Listener
  const sliderParticles = document.getElementById('slider-particles');
  const valParticles = document.getElementById('val-particles');
  sliderParticles.addEventListener('input', (e) => {
    const val = parseInt(e.target.value);
    currentParticleCount = val;
    valParticles.textContent = val.toLocaleString();
  });

  const toggleContours = document.getElementById('toggle-contours');
  toggleContours.addEventListener('change', (e) => {
    enableContours = e.target.checked;
  });

  const selectPalette = document.getElementById('select-palette');
  selectPalette.addEventListener('change', (e) => {
    currentPalette = e.target.value;
  });

  const selectRenderMode = document.getElementById('select-render-mode');
  selectRenderMode.addEventListener('change', (e) => {
    currentRenderMode = e.target.value;
  });

  // E. Canvas Interactions
  canvas.addEventListener('mousedown', (e) => {
    isDragging = true;
    handleCanvasInteraction(e.clientX, e.clientY, false);
  });

  window.addEventListener('mousemove', (e) => {
    if (isDragging) {
      handleCanvasInteraction(e.clientX, e.clientY, true);
    }
  });

  window.addEventListener('mouseup', () => {
    isDragging = false;
  });

  // Mobile Touch Support
  canvas.addEventListener('touchstart', (e) => {
    isDragging = true;
    const touch = e.touches[0];
    handleCanvasInteraction(touch.clientX, touch.clientY, false);
    e.preventDefault(); // Stop mobile scrolling
  });

  window.addEventListener('touchmove', (e) => {
    if (isDragging && e.touches.length > 0) {
      const touch = e.touches[0];
      handleCanvasInteraction(touch.clientX, touch.clientY, true);
    }
  }, { passive: false });

  window.addEventListener('touchend', () => {
    isDragging = false;
  });
}

// Run app init on load
window.addEventListener('DOMContentLoaded', init);
