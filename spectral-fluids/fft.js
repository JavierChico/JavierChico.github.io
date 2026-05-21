/**
 * fft.js - Highly optimized, zero-allocation, strided 2D Fast Fourier Transform.
 * Specialized for row (stride=1) and column (stride=N) operations.
 * Bypasses all multiplications/divisions in index calculations within hot loops.
 */

export class FFT1D {
  constructor(N) {
    this.N = N;
    this.levels = Math.log2(N);
    if (!Number.isInteger(this.levels)) {
      throw new Error(`FFT size N must be a power of 2, got ${N}`);
    }

    // Precompute bit-reversal lookup table
    this.rev = new Int32Array(N);
    for (let i = 0; i < N; i++) {
      let r = 0;
      for (let j = 0; j < this.levels; j++) {
        if ((i & (1 << j)) !== 0) {
          r |= (1 << (this.levels - 1 - j));
        }
      }
      this.rev[i] = r;
    }

    // Precompute twiddle factors (cos/sin tables) to avoid math overhead in hot loops
    this.cosTable = new Float64Array(N * this.levels);
    this.sinTable = new Float64Array(N * this.levels);
    for (let level = 1; level <= this.levels; level++) {
      const m = 1 << level;
      const theta = -2 * Math.PI / m;
      const wOffset = (level - 1) * N;
      for (let j = 0; j < m / 2; j++) {
        this.cosTable[wOffset + j] = Math.cos(j * theta);
        this.sinTable[wOffset + j] = Math.sin(j * theta);
      }
    }
  }

  /**
   * Highly optimized Row FFT/IFFT (stride = 1).
   * Sequential memory accesses, zero index multiplications in hot loops.
   */
  fftRow(re, im, offset, isInverse = false) {
    const N = this.N;
    const rev = this.rev;

    // 1. Bit-reversal permutation (stride = 1)
    for (let i = 0; i < N; i++) {
      const j = rev[i];
      if (i < j) {
        const iIdx = offset + i;
        const jIdx = offset + j;
        
        const tempRe = re[iIdx];
        re[iIdx] = re[jIdx];
        re[jIdx] = tempRe;
        
        const tempIm = im[iIdx];
        im[iIdx] = im[jIdx];
        im[jIdx] = tempIm;
      }
    }

    // 2. Cooley-Tukey Decimation-in-Time Butterfly (stride = 1)
    const levels = this.levels;
    for (let level = 1; level <= levels; level++) {
      const m = 1 << level;
      const mHalf = m >> 1;
      const wOffset = (level - 1) * N;

      for (let k = 0; k < N; k += m) {
        const uStart = offset + k;
        const tStart = uStart + mHalf;
        for (let j = 0; j < mHalf; j++) {
          const wr = this.cosTable[wOffset + j];
          let wi = this.sinTable[wOffset + j];
          if (isInverse) wi = -wi;

          const uIdx = uStart + j;
          const tIdx = tStart + j;

          const tRe = re[tIdx];
          const tIm = im[tIdx];

          // Complex multiplication: t * w
          const tr = tRe * wr - tIm * wi;
          const ti = tRe * wi + tIm * wr;

          re[tIdx] = re[uIdx] - tr;
          im[tIdx] = im[uIdx] - ti;

          re[uIdx] = re[uIdx] + tr;
          im[uIdx] = im[uIdx] + ti;
        }
      }
    }

    // 3. Scaling for Inverse FFT
    if (isInverse) {
      for (let i = 0; i < N; i++) {
        const idx = offset + i;
        re[idx] /= N;
        im[idx] /= N;
      }
    }
  }

  /**
   * Highly optimized Column FFT/IFFT (stride = N).
   * Pointer increments for index updates, zero index multiplications in hot loops.
   */
  fftCol(re, im, offset, isInverse = false) {
    const N = this.N;
    const rev = this.rev;

    // 1. Bit-reversal permutation (stride = N)
    for (let i = 0; i < N; i++) {
      const j = rev[i];
      if (i < j) {
        const iIdx = offset + i * N;
        const jIdx = offset + j * N;
        
        const tempRe = re[iIdx];
        re[iIdx] = re[jIdx];
        re[jIdx] = tempRe;
        
        const tempIm = im[iIdx];
        im[iIdx] = im[jIdx];
        im[jIdx] = tempIm;
      }
    }

    // 2. Cooley-Tukey Butterfly (stride = N)
    const levels = this.levels;
    for (let level = 1; level <= levels; level++) {
      const m = 1 << level;
      const mHalf = m >> 1;
      const wOffset = (level - 1) * N;

      const mHalfN = mHalf * N;

      for (let k = 0; k < N; k += m) {
        const uStart = offset + k * N;
        const tStart = uStart + mHalfN;
        
        let uIdx = uStart;
        let tIdx = tStart;
        for (let j = 0; j < mHalf; j++) {
          const wr = this.cosTable[wOffset + j];
          let wi = this.sinTable[wOffset + j];
          if (isInverse) wi = -wi;

          const tRe = re[tIdx];
          const tIm = im[tIdx];

          // Complex multiplication: t * w
          const tr = tRe * wr - tIm * wi;
          const ti = tRe * wi + tIm * wr;

          re[tIdx] = re[uIdx] - tr;
          im[tIdx] = im[uIdx] - ti;

          re[uIdx] = re[uIdx] + tr;
          im[uIdx] = im[uIdx] + ti;

          uIdx += N;
          tIdx += N;
        }
      }
    }

    // 3. Scaling for Inverse FFT
    if (isInverse) {
      for (let i = 0; i < N; i++) {
        const idx = offset + i * N;
        re[idx] /= N;
        im[idx] /= N;
      }
    }
  }
}

export class FFT2D {
  constructor(N) {
    this.N = N;
    this.fft1d = new FFT1D(N);
  }

  /**
   * Performs an in-place 2D Forward FFT using row-column decomposition.
   */
  forward(re, im) {
    const N = this.N;
    // 1. FFT along each row (stride = 1)
    for (let r = 0; r < N; r++) {
      this.fft1d.fftRow(re, im, r * N, false);
    }
    // 2. FFT along each column (stride = N)
    for (let c = 0; c < N; c++) {
      this.fft1d.fftCol(re, im, c, false);
    }
  }

  /**
   * Performs an in-place 2D Inverse FFT using row-column decomposition.
   */
  inverse(re, im) {
    const N = this.N;
    // 1. IFFT along each row (stride = 1)
    for (let r = 0; r < N; r++) {
      this.fft1d.fftRow(re, im, r * N, true);
    }
    // 2. IFFT along each column (stride = N)
    for (let c = 0; c < N; c++) {
      this.fft1d.fftCol(re, im, c, true);
    }
  }
}
