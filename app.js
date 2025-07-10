// Ultimate Complex Audio Transformer - Always Real-time
let ctx;
let audioBuffer;
let processedBuffer;
let audioSource;
let realtimeProcessor;
let currentSource;
let isPlaying = false;
let isRealTimeMode = false;
let selectedFunction = 'power';
let analyser;
let canvas;
let canvasCtx;
let animationId;
let scriptProcessor;

// Real-time processing parameters
let realtimeParams = {
  fftSize: 2048,
  amount: 1.0,
  smoothing: 0.1,
  mix: 1.0,
  gain: 1.0,
  frequencyRange: 0, // 0: Normal, 1: Extended, 2: Unlimited
  antiAliasingMode: 0, // 0: Filter, 1: Wrap, 2: Off
  frequencyMapping: 0, // 0: Normal, 1: Human Range Scaled
  modEnabled: false,   // Whether amount modulation is enabled
  modPeriod: 5.0,      // Modulation period in seconds
  modStartTime: 0,     // Time when modulation started
  modMinAmount: 0.0,   // Minimum modulation amount
  modMaxAmount: 1.0    // Maximum modulation amount
};

// Memory management pool
let memoryPool = {
  realBuffers: [],
  imagBuffers: [],
  windowBuffers: new Map(),
  fftBuffers: []
};

// Complex function definitions with complete implementations
const complexFunctions = {
  // === BASIC MATHEMATICAL OPERATIONS ===
  'power': { name: 'Power (f^n)', desc: 'Frequency stretching/compression' },
  'sqrt': { name: '√z', desc: 'Square root expansion' },
  'cbrt': { name: '∛z', desc: 'Cube root expansion' },
  'nroot': { name: 'ⁿ√z', desc: 'N-th root (power-controlled)' },
  'square': { name: 'z²', desc: 'Square compression' },
  'cube': { name: 'z³', desc: 'Cube transformation' },
  'quartic': { name: 'z⁴', desc: 'Fourth power' },
  'quintic': { name: 'z⁵', desc: 'Fifth power' },
  'inverse': { name: '1/z', desc: 'Reciprocal inversion' },

  // === LOGARITHMIC & EXPONENTIAL ===
  'log': { name: 'ln(z)', desc: 'Natural logarithm' },
  'log10': { name: 'log₁₀(z)', desc: 'Base-10 logarithm' },
  'log2': { name: 'log₂(z)', desc: 'Base-2 logarithm' },
  'exp': { name: 'e^z', desc: 'Natural exponential' },
  'exp2': { name: '2^z', desc: 'Base-2 exponential' },
  'exp10': { name: '10^z', desc: 'Base-10 exponential' },

  // === TRIGONOMETRIC FUNCTIONS ===
  'sin': { name: 'sin(z)', desc: 'Sine wave warping' },
  'cos': { name: 'cos(z)', desc: 'Cosine wave warping' },
  'tan': { name: 'tan(z)', desc: 'Tangent sharp bending' },
  'csc': { name: 'csc(z)', desc: 'Cosecant (1/sin)' },
  'sec': { name: 'sec(z)', desc: 'Secant (1/cos)' },
  'cot': { name: 'cot(z)', desc: 'Cotangent (1/tan)' },

  // === INVERSE TRIGONOMETRIC ===
  'asin': { name: 'arcsin(z)', desc: 'Inverse sine mapping' },
  'acos': { name: 'arccos(z)', desc: 'Inverse cosine mapping' },
  'atan': { name: 'arctan(z)', desc: 'S-curve frequency mapping' },
  'atan2': { name: 'atan2(y,x)', desc: 'Two-argument arctangent' },

  // === HYPERBOLIC FUNCTIONS (FIXED) ===
  'sinh': { name: 'sinh(z)', desc: 'Hyperbolic sine growth' },
  'cosh': { name: 'cosh(z)', desc: 'Hyperbolic cosine catenary' },
  'tanh': { name: 'tanh(z)', desc: 'Hyperbolic tangent saturation' },
  'csch': { name: 'csch(z)', desc: 'Hyperbolic cosecant' },
  'sech': { name: 'sech(z)', desc: 'Hyperbolic secant' },
  'coth': { name: 'coth(z)', desc: 'Hyperbolic cotangent' },

  // === INVERSE HYPERBOLIC ===
  'asinh': { name: 'asinh(z)', desc: 'Inverse hyperbolic sine' },
  'acosh': { name: 'acosh(z)', desc: 'Inverse hyperbolic cosine' },
  'atanh': { name: 'atanh(z)', desc: 'Inverse hyperbolic tangent' },

  // === COMPLEX NUMBER OPERATIONS ===
  'complex_magnitude': { name: '|z|', desc: 'Distance from origin' },
  'complex_phase': { name: 'arg(z)', desc: 'Angle transformation' },
  'conjugate': { name: 'z*', desc: 'Complex conjugate mirror' },
  'reciprocal': { name: '1/z', desc: 'Complex reciprocal' },
  'abs_square': { name: '|z|²', desc: 'Squared magnitude' },

  // === FAMOUS PHYSICS EQUATIONS ===
  'schrodinger': { name: 'Schrödinger Ψ', desc: 'Quantum wave function' },
  'heisenberg': { name: 'Heisenberg Δx·Δp', desc: 'Uncertainty principle' },
  'einstein_mass_energy': { name: 'E=mc²', desc: 'Mass-energy equivalence' },
  'planck_radiation': { name: 'Planck B(ν,T)', desc: 'Blackbody radiation' },
  'maxwell_boltzmann': { name: 'Maxwell-Boltzmann', desc: 'Velocity distribution' },
  'wave_equation': { name: 'Wave ∂²u/∂t²', desc: 'String vibrations' },
  'heat_equation': { name: 'Heat ∂u/∂t', desc: 'Thermal diffusion' },
  'lorentz_factor': { name: 'Lorentz γ', desc: 'Relativistic factor' },
  'schwarzschild': { name: 'Schwarzschild r_s', desc: 'Black hole radius' },
  'compton_scattering': { name: 'Compton Δλ', desc: 'Photon scattering' },

  // === SPECIAL FUNCTIONS ===
  'gamma': { name: 'Γ(z)', desc: 'Gamma function' },
  'factorial': { name: 'z!', desc: 'Factorial function' },
  'zeta': { name: 'ζ(z)', desc: 'Riemann zeta function' },
  'eta': { name: 'η(z)', desc: 'Dirichlet eta function' },
  'bessel_j0': { name: 'J₀(z)', desc: 'Bessel function of 1st kind' },
  'bessel_y0': { name: 'Y₀(z)', desc: 'Bessel function of 2nd kind' },
  'bessel_i0': { name: 'I₀(z)', desc: 'Modified Bessel function' },
  'bessel_k0': { name: 'K₀(z)', desc: 'Modified Bessel function 2nd' },
  'airy_ai': { name: 'Ai(z)', desc: 'Airy function Ai' },
  'airy_bi': { name: 'Bi(z)', desc: 'Airy function Bi' },
  'erf': { name: 'erf(z)', desc: 'Error function' },
  'erfc': { name: 'erfc(z)', desc: 'Complementary error function' },

  // === GEOMETRIC TRANSFORMATIONS ===
  'mobius': { name: 'Möbius (az+b)/(cz+d)', desc: 'Conformal mapping' },
  'joukowsky': { name: 'Joukowsky z+1/z', desc: 'Airfoil transformation' },
  'schwarz_christoffel': { name: 'Schwarz-Christoffel', desc: 'Polygon mapping' },
  'cayley': { name: 'Cayley (z-i)/(z+i)', desc: 'Disk to half-plane' },

  // === SPIRAL TRANSFORMATIONS ===
  'spiral': { name: 'Log Spiral', desc: 'Logarithmic spiral' },
  'fibonacci_spiral': { name: 'Fibonacci Spiral', desc: 'Golden ratio spiral' },
  'archimedean_spiral': { name: 'Archimedean rθ', desc: 'Linear spiral' },
  'fermat_spiral': { name: 'Fermat r²=aθ', desc: 'Parabolic spiral' },
  'hyperbolic_spiral': { name: 'Hyperbolic r=a/θ', desc: 'Reciprocal spiral' },
  'ulam_spiral': { name: 'Ulam Spiral', desc: 'Prime number spiral' },

  // === FRACTAL FUNCTIONS ===
  'mandelbrot': { name: 'Mandelbrot z²+c', desc: 'Classic fractal' },
  'julia': { name: 'Julia Set', desc: 'Julia set fractal' },
  'burning_ship': { name: 'Burning Ship |z|²+c', desc: 'Absolute value fractal' },
  'tricorn': { name: 'Tricorn z̄²+c', desc: 'Conjugate fractal' },
  'multibrot': { name: 'Multibrot z^n+c', desc: 'Generalized Mandelbrot' },
  'newton_raphson': { name: 'Newton z³-1', desc: 'Newton fractal' },
  'nova': { name: 'Nova Fractal', desc: 'Newton-style fractal' },
  'phoenix': { name: 'Phoenix z²+c+b', desc: 'Phoenix fractal' },
  'magnet': { name: 'Magnet (z²+c-1)²', desc: 'Magnet fractal' },
  'lambda': { name: 'Lambda λz(1-z)', desc: 'Lambda fractal' },
  'spider': { name: 'Spider z²+c/z', desc: 'Spider fractal' },
  'henon_map': { name: 'Hénon Map', desc: '2D chaotic attractor' },
  'lorenz_attractor': { name: 'Lorenz Attractor', desc: 'Butterfly effect chaos' },
  'logistic_map': { name: 'Logistic rx(1-x)', desc: 'Period-doubling chaos' },
  'tent_map': { name: 'Tent Map', desc: 'Piecewise linear chaos' },

  // === WAVE FUNCTIONS ===
  'wave': { name: 'Wave Interference', desc: 'Constructive/destructive' },
  'ripple': { name: 'Ripple Effect', desc: 'Concentric waves' },
  'standing_wave': { name: 'Standing Wave', desc: 'Node/antinode pattern' },
  'beat_frequency': { name: 'Beat Frequency', desc: 'Amplitude modulation' },
  'doppler_effect': { name: 'Doppler f(1±v/c)', desc: 'Moving source shift' },
  'gaussian_wave': { name: 'Gaussian Wave', desc: 'Bell curve packet' },
  'sinc_function': { name: 'Sinc sin(x)/x', desc: 'Cardinal sine' },
  'chirp': { name: 'Chirp f(t)', desc: 'Frequency sweep' },
  'morlet_wavelet': { name: 'Morlet Wavelet', desc: 'Gaussian windowed sine' },
  'mexican_hat': { name: 'Mexican Hat', desc: 'Ricker wavelet' },

  // === SPECIAL PATTERNS ===
  'dirac_delta': { name: 'Dirac δ(x)', desc: 'Impulse function' },
  'heaviside_step': { name: 'Heaviside H(x)', desc: 'Unit step function' },
  'sawtooth': { name: 'Sawtooth Wave', desc: 'Linear ramp' },
  'square_wave': { name: 'Square Wave', desc: 'Binary switching' },
  'triangle_wave': { name: 'Triangle Wave', desc: 'Triangular oscillation' },
  'pulse_train': { name: 'Pulse Train', desc: 'Periodic impulses' },

  // === STEREO PANNING EFFECTS ===
  'stereo_real_imag': { name: '🎧 Real/Imaginary Split', desc: 'Real→Left, Imag→Right' },
  'stereo_mag_phase': { name: '🎧 Magnitude/Phase Split', desc: 'Mag→Left, Phase→Right' },
  'stereo_rotation': { name: '🎧 Complex Rotation', desc: 'Rotating panning field' },
  'stereo_spiral_pan': { name: '🎧 Spiral Panning', desc: 'Frequency-dependent pan' },
  'stereo_doppler_pan': { name: '🎧 Doppler Panning', desc: 'Moving source simulation' },
  'stereo_binaural': { name: '🎧 Binaural Beats', desc: 'Frequency difference beats' },

  // === QUANTUM & EXOTIC ===
  'quantum_harmonic': { name: 'Quantum Harmonic', desc: 'QHO wave functions' },
  'hydrogen_orbital': { name: 'Hydrogen Orbital', desc: 'Atomic wave functions' },
  'tunneling_effect': { name: 'Quantum Tunneling', desc: 'Barrier penetration' },
  'klein_gordon': { name: 'Klein-Gordon', desc: 'Relativistic wave equation' },
  'dirac_equation': { name: 'Dirac Equation', desc: 'Relativistic electron' },
  'feynman_diagram': { name: 'Feynman Propagator', desc: 'Particle interactions' }
};

// Memory management functions
function getBuffer(size, type = 'real') {
  // Ensure memoryPool exists
  if (!memoryPool) {
    memoryPool = {
      realBuffers: [],
      imagBuffers: [],
      windowBuffers: new Map(),
      fftBuffers: []
    };
  }

  const pool = type === 'real' ? memoryPool.realBuffers : memoryPool.imagBuffers;
  if (!pool) {
    console.warn(`Memory pool for type '${type}' is undefined, creating new array`);
    if (type === 'real') {
      memoryPool.realBuffers = [];
      return new Float32Array(size);
    } else {
      memoryPool.imagBuffers = [];
      return new Float32Array(size);
    }
  }

  let buffer = pool.find(b => b && b.length >= size);
  if (!buffer) {
    buffer = new Float32Array(size);
    pool.push(buffer);
  }
  return buffer.subarray(0, size);
}

function getWindow(size) {
  // Ensure memoryPool exists
  if (!memoryPool || !memoryPool.windowBuffers) {
    if (!memoryPool) {
      memoryPool = {
        realBuffers: [],
        imagBuffers: [],
        windowBuffers: new Map(),
        fftBuffers: []
      };
    } else if (!memoryPool.windowBuffers) {
      memoryPool.windowBuffers = new Map();
    }
  }

  if (!memoryPool.windowBuffers.has(size)) {
    const window = new Float32Array(size);
    // Hann window
    for (let i = 0; i < size; i++) {
      window[i] = 0.5 * (1 - Math.cos(2 * Math.PI * i / (size - 1)));
    }
    memoryPool.windowBuffers.set(size, window);
  }
  return memoryPool.windowBuffers.get(size);
}

function updateMemoryInfo() {
  try {
    // Ensure memoryPool exists and has required properties
    if (!memoryPool) {
      console.warn('memoryPool is undefined, reinitializing...');
      memoryPool = {
        realBuffers: [],
        imagBuffers: [],
        windowBuffers: new Map(),
        fftBuffers: []
      };
    }

    // Ensure each property exists
    if (!memoryPool.realBuffers) memoryPool.realBuffers = [];
    if (!memoryPool.imagBuffers) memoryPool.imagBuffers = [];
    if (!memoryPool.windowBuffers) memoryPool.windowBuffers = new Map();
    if (!memoryPool.fftBuffers) memoryPool.fftBuffers = [];

    const realCount = memoryPool.realBuffers.length;
    const imagCount = memoryPool.imagBuffers.length;
    const windowCount = memoryPool.windowBuffers.size;
    const fftCount = memoryPool.fftBuffers.length;

    const memoryInfoElement = document.getElementById('memoryInfo');
    if (memoryInfoElement) {
      memoryInfoElement.textContent =
        `Memory Pool: ${realCount} real, ${imagCount} imag, ${windowCount} window, ${fftCount} FFT buffers`;
    }
  } catch (error) {
    console.error('Error in updateMemoryInfo:', error);
    // Don't throw the error, just log it to prevent breaking the audio processing
  }
}

// Audio context initialization
async function initAudioContext() {
  try {
    if (!ctx) {
      ctx = new (window.AudioContext || window.webkitAudioContext)();
      console.log('Audio context created');
    }

    if (ctx.state === 'suspended') {
      console.log('Audio context suspended, attempting to resume...');
      await ctx.resume();
      console.log('Audio context resumed successfully');
    }

    console.log(`Audio context state: ${ctx.state}, sample rate: ${ctx.sampleRate}`);
    return ctx;

  } catch (error) {
    console.error('Failed to initialize audio context:', error);
    showStatus(`Audio context error: ${error.message}. Try clicking a button first.`, 'error');
    throw error;
  }
}

// Status display
function showStatus(message, type = 'info') {
  const statusDiv = document.getElementById('status');
  statusDiv.textContent = message;
  statusDiv.className = type;
  console.log(`[${type}] ${message}`);
}

// Real-time Audio Processor with large buffer for smooth processing
class RealtimeAudioProcessor {
  constructor(audioContext, bufferDuration = 3.0) {
    this.audioContext = audioContext;
    this.sampleRate = audioContext.sampleRate;
    this.bufferDuration = bufferDuration;

    // Use larger buffer size for smoother processing
    this.bufferSize = 8192; // Larger ScriptProcessorNode buffer

    // Large circular buffer for smooth processing (3+ seconds)
    this.circularBufferSize = Math.ceil(this.sampleRate * bufferDuration);
    this.inputCircularBuffer = new Float32Array(this.circularBufferSize);
    this.outputCircularBuffer = new Float32Array(this.circularBufferSize);

    // Buffer pointers
    this.writePtr = 0;
    this.readPtr = 0;
    this.bufferFilled = false;

    // Initialize with 3 second latency for smooth processing
    this.latencyBufferSize = Math.floor(this.sampleRate * 1.5); // 1.5 second initial latency
    this.readPtr = this.circularBufferSize - this.latencyBufferSize;

    console.log(`Initialized large buffer: ${this.circularBufferSize} samples (${(this.circularBufferSize / this.sampleRate).toFixed(2)}s) at ${this.sampleRate}Hz`);

    // Processing parameters
    this.params = {
      function: 'power',
      amount: 1.0,
      smoothing: 0.1,
      mix: 1.0,
      gain: 1.0
    };

    // FFT setup for real-time processing
    this.fftSize = 2048;
    this.windowFunc = this.createWindow(this.fftSize);

    // Processing state
    this.lastSample = 0;
    this.processCounter = 0;

    // Create ScriptProcessorNode with larger buffer
    this.processor = audioContext.createScriptProcessor(this.bufferSize, 1, 1);
    this.processor.onaudioprocess = this.processAudio.bind(this);

    // Gain node for output control
    this.gainNode = audioContext.createGain();
    this.processor.connect(this.gainNode);

    console.log(`RealtimeAudioProcessor initialized with ${this.bufferSize} sample buffer and ${(this.circularBufferSize / this.sampleRate).toFixed(2)}s circular buffer`);
  }

  createWindow(size) {
    const window = new Float32Array(size);
    for (let i = 0; i < size; i++) {
      window[i] = 0.5 * (1 - Math.cos(2 * Math.PI * i / (size - 1)));
    }
    return window;
  }

  updateParams(newParams) {
    // Update parameters safely
    Object.assign(this.params, newParams);
    this.gainNode.gain.value = this.params.gain;

    // Update FFT size if changed
    if (newParams.fftSize && newParams.fftSize !== this.fftSize) {
      this.fftSize = newParams.fftSize;
      this.windowFunc = this.createWindow(this.fftSize);
      this.lastSample = 0;
      console.log(`FFT size updated to ${this.fftSize}`);
    }

    // Log parameter updates for debugging
    if (this.processCounter % 50 === 0) {
      console.log('Updated params:', {
        function: this.params.function,
        amount: this.params.amount,
        mix: this.params.mix,
        fftSize: this.fftSize
      });
    }
  }

  processAudio(event) {
    const inputData = event.inputBuffer.getChannelData(0);
    const outputData = event.outputBuffer.getChannelData(0);
    const bufferLength = inputData.length;

    // Write input data to circular buffer
    for (let i = 0; i < bufferLength; i++) {
      this.inputCircularBuffer[this.writePtr] = inputData[i];
      this.writePtr = (this.writePtr + 1) % this.circularBufferSize;
    }

    // Simple direct processing to avoid circular buffer complexity
    // Process the current input buffer directly with proper windowing

    if (bufferLength >= this.fftSize) {
      // Process in overlapping windows
      const hopSize = this.fftSize / 4; // 75% overlap for very smooth processing
      const tempOutput = new Float32Array(bufferLength);
      tempOutput.fill(0); // Initialize to zero for overlap-add

      for (let windowStart = 0; windowStart < bufferLength - this.fftSize + hopSize; windowStart += hopSize) {
        const windowEnd = Math.min(windowStart + this.fftSize, bufferLength);
        const actualWindowSize = windowEnd - windowStart;

        if (actualWindowSize < this.fftSize / 2) continue; // Skip tiny windows

        // Extract window
        const windowBuffer = new Float32Array(this.fftSize);
        windowBuffer.fill(0);

        for (let i = 0; i < actualWindowSize; i++) {
          windowBuffer[i] = inputData[windowStart + i];
        }

        // Process window
        const processedWindow = this.processChunk(windowBuffer);

        // Overlap-add to output with proper scaling
        const overlapScale = 1.0 / 3.0; // Scale for 75% overlap
        for (let i = 0; i < actualWindowSize; i++) {
          if (windowStart + i < bufferLength) {
            tempOutput[windowStart + i] += processedWindow[i] * overlapScale;
          }
        }
      }

      // Copy to final output with mixing and gain
      for (let i = 0; i < bufferLength; i++) {
        let processedSample = tempOutput[i];

        // Apply smoothing to reduce artifacts
        if (this.lastSample !== undefined) {
          const smoothingAmount = Math.min(0.1, this.params.smoothing);
          processedSample = smoothingAmount * this.lastSample + (1 - smoothingAmount) * processedSample;
        }
        this.lastSample = processedSample;

        // Apply dry/wet mix
        const originalSample = inputData[i];
        let mixedSample = this.params.mix * processedSample + (1 - this.params.mix) * originalSample;

        // Apply gain with soft limiting
        mixedSample *= this.params.gain;
        if (Math.abs(mixedSample) > 0.95) {
          mixedSample = Math.sign(mixedSample) * (0.95 + 0.05 * Math.tanh((Math.abs(mixedSample) - 0.95) * 10));
        }

        outputData[i] = mixedSample;
      }

    } else {
      // Buffer smaller than FFT size - simple processing
      const paddedBuffer = new Float32Array(this.fftSize);
      paddedBuffer.fill(0);

      for (let i = 0; i < bufferLength; i++) {
        paddedBuffer[i] = inputData[i];
      }

      const processedBuffer = this.processChunk(paddedBuffer);

      for (let i = 0; i < bufferLength; i++) {
        let processedSample = processedBuffer[i];

        // Apply smoothing
        if (this.lastSample !== undefined) {
          const smoothingAmount = Math.min(0.1, this.params.smoothing);
          processedSample = smoothingAmount * this.lastSample + (1 - smoothingAmount) * processedSample;
        }
        this.lastSample = processedSample;

        // Apply mix and gain
        const originalSample = inputData[i];
        let mixedSample = this.params.mix * processedSample + (1 - this.params.mix) * originalSample;

        mixedSample *= this.params.gain;
        if (Math.abs(mixedSample) > 0.95) {
          mixedSample = Math.sign(mixedSample) * (0.95 + 0.05 * Math.tanh((Math.abs(mixedSample) - 0.95) * 10));
        }

        outputData[i] = mixedSample;
      }
    }

    // Update processing counter for diagnostics
    this.processCounter++;
    if (this.processCounter % 200 === 0) {
      const bufferFillPercent = (this.getAvailableData() / this.circularBufferSize * 100).toFixed(1);
      console.log(`Processing stable - Buffer: ${bufferFillPercent}%`);
    }
  }

  // Helper method to get available data in circular buffer
  getAvailableData() {
    if (this.writePtr >= this.readPtr) {
      return this.writePtr - this.readPtr;
    } else {
      return this.circularBufferSize - this.readPtr + this.writePtr;
    }
  }

  // Process a single chunk through FFT
  processChunk(inputChunk) {
    const real = new Float32Array(this.fftSize);
    const imag = new Float32Array(this.fftSize);

    // Copy input with windowing
    for (let i = 0; i < Math.min(inputChunk.length, this.fftSize); i++) {
      real[i] = inputChunk[i] * this.windowFunc[i];
      imag[i] = 0;
    }

    // Forward FFT
    fft(real, imag);

    // Apply frequency domain transformation
    try {
      applyComplexFunction(real, imag, this.fftSize, this.sampleRate, this.params.function, this.params.amount);
    } catch (error) {
      console.warn('Frequency transformation error:', error);
    }

    // Inverse FFT
    ifft(real, imag);

    // Apply windowing again and return
    const output = new Float32Array(this.fftSize);
    for (let i = 0; i < this.fftSize; i++) {
      output[i] = real[i] * this.windowFunc[i];
    }

    return output;
  }

  connect(destination) {
    this.gainNode.connect(destination);
    return this.processor;
  }

  disconnect() {
    if (this.processor) {
      this.processor.disconnect();
    }
    if (this.gainNode) {
      this.gainNode.disconnect();
    }
  }
}

// Complex function applications with proper frequency transformations
function applyComplexFunction(real, imag, fftSize, sampleRate, funcName, amount = 1.0) {
  try {
    const nyquist = sampleRate / 2;
    const binSize = sampleRate / fftSize;

    // Frequency range limits based on user settings
    let minFreq, maxFreq;
    switch (realtimeParams.frequencyRange) {
      case 0: // Normal (audible range)
        minFreq = 20;
        maxFreq = Math.min(20000, nyquist);
        break;
      case 1: // Extended
        minFreq = 5;
        maxFreq = Math.min(40000, nyquist);
        break;
      case 2: // Unlimited
        minFreq = 0;
        maxFreq = nyquist;
        break;
    }

    // Create temporary arrays for frequency remapping
    const tempReal = new Float32Array(fftSize);
    const tempImag = new Float32Array(fftSize);

    // Copy original spectrum
    for (let i = 0; i < fftSize; i++) {
      tempReal[i] = real[i];
      tempImag[i] = imag[i];
    }

    // Clear output arrays
    real.fill(0);
    imag.fill(0);

    // Process only positive frequencies (first half of spectrum)
    for (let i = 1; i < fftSize / 2; i++) {
      const frequency = i * binSize;

      // Skip frequencies outside range based on current settings
      if (frequency < minFreq || frequency > maxFreq) continue;

      // Normalize frequency to 0-1 within selected range
      const normalizedFreq = (frequency - minFreq) / (maxFreq - minFreq);

      // Get magnitude and phase from original spectrum
      const magnitude = Math.sqrt(tempReal[i] * tempReal[i] + tempImag[i] * tempImag[i]);
      const phase = Math.atan2(tempImag[i], tempReal[i]);

      if (magnitude < 1e-12) continue; // Skip silent bins

      let mappedFreq = normalizedFreq; // Default: no change
      let magnitudeMultiplier = 1.0;
      let phaseShift = 0;

      // Apply dramatic frequency transformations based on function type
      // All functions use 'amount' parameter to control intensity (0-10 typical range)
      switch (funcName) {
        // === BASIC OPERATIONS ===
        case 'power':
          // Frequency stretching/compression: amount controls exponent (0.1 to 5.0)
          const exponent = 0.1 + (amount * 0.5);
          mappedFreq = Math.pow(normalizedFreq, exponent);
          break;

        case 'sqrt':
          // Square root expansion: spreads low frequencies
          mappedFreq = Math.pow(normalizedFreq, 0.5 - amount * 0.03);
          break;

        case 'cbrt':
          // Cube root: even more expansion
          mappedFreq = Math.pow(normalizedFreq, 0.33 - amount * 0.02);
          break;

        case 'nroot':
          // N-th root: amount controls root (1-20)
          const rootPower = 1 + amount;
          mappedFreq = Math.pow(normalizedFreq, 1 / rootPower);
          break;

        case 'square':
          // Square compression: compresses to higher frequencies
          mappedFreq = Math.pow(normalizedFreq, 2 + amount * 0.5);
          break;

        case 'cube':
          // Cube: severe compression
          mappedFreq = Math.pow(normalizedFreq, 3 + amount);
          break;

        case 'inverse':
          // Frequency inversion with scaling
          mappedFreq = 1 - normalizedFreq * (1 - amount * 0.1);
          magnitudeMultiplier = 1 + amount * 0.3;
          break;

        // === LOGARITHMIC & EXPONENTIAL ===
        case 'log':
          // Logarithmic stretching: expands low frequencies dramatically
          mappedFreq = Math.log(1 + normalizedFreq * (Math.E - 1) * amount) / amount;
          break;

        case 'log10':
          // Base-10 log: different curve shape
          mappedFreq = Math.log10(1 + normalizedFreq * 9 * amount) / Math.log10(1 + 9 * amount);
          break;

        case 'log2':
          // Base-2 log: computer-friendly stretching
          mappedFreq = Math.log2(1 + normalizedFreq * 3 * amount) / Math.log2(1 + 3 * amount);
          break;

        case 'exp':
          // Exponential compression: compresses to high frequencies
          mappedFreq = (Math.exp(normalizedFreq * amount) - 1) / (Math.exp(amount) - 1);
          break;

        case 'exp2':
          // Base-2 exponential: sharp compression
          mappedFreq = (Math.pow(2, normalizedFreq * amount) - 1) / (Math.pow(2, amount) - 1);
          break;

        case 'exp10':
          // Base-10 exponential: extreme compression
          mappedFreq = (Math.pow(10, normalizedFreq * amount * 0.5) - 1) / (Math.pow(10, amount * 0.5) - 1);
          break;

        // === TRIGONOMETRIC FUNCTIONS ===
        case 'sin':
          // Sine wave warping: creates multiple bands
          mappedFreq = (1 + Math.sin((normalizedFreq * Math.PI * 2 * amount) - Math.PI / 2)) / 2;
          magnitudeMultiplier = 1 + 0.3 * Math.sin(normalizedFreq * Math.PI * amount);
          phaseShift = Math.sin(normalizedFreq * Math.PI * amount) * 0.5;
          break;

        case 'cos':
          // Cosine wave warping: different phase than sine
          mappedFreq = (1 + Math.cos(normalizedFreq * Math.PI * amount)) / 2;
          magnitudeMultiplier = 1 + 0.3 * Math.cos(normalizedFreq * Math.PI * amount);
          break;

        case 'tan':
          // Tangent: creates sharp frequency bands
          const tanInput = (normalizedFreq - 0.5) * Math.PI * 0.8 * amount;
          mappedFreq = (Math.atan(Math.tan(tanInput)) / (Math.PI * 0.8 * amount)) + 0.5;
          mappedFreq = Math.max(0, Math.min(1, mappedFreq));
          break;

        case 'asin':
          // Arcsine: S-curve stretching
          const asinInput = (normalizedFreq * 2 - 1) * Math.min(0.99, amount * 0.1);
          mappedFreq = (Math.asin(asinInput) / Math.asin(Math.min(0.99, amount * 0.1)) + 1) / 2;
          break;

        case 'acos':
          // Arccosine: inverse S-curve
          const acosInput = Math.min(0.99, amount * 0.1);
          mappedFreq = Math.acos(normalizedFreq * acosInput) / Math.acos(acosInput);
          break;

        case 'atan':
          // Arctangent: gentle S-curve with amount scaling
          const atanInput = (normalizedFreq - 0.5) * amount * 2;
          mappedFreq = (Math.atan(atanInput) / Math.atan(amount * 2)) * 0.5 + 0.5;
          break;

        // === HYPERBOLIC FUNCTIONS ===
        case 'sinh':
          // Hyperbolic sine: exponential-like growth
          const sinhInput = (normalizedFreq - 0.5) * amount;
          mappedFreq = Math.tanh(Math.sinh(sinhInput) * 0.5) * 0.5 + 0.5;
          magnitudeMultiplier = 1 + 0.2 * Math.abs(Math.sinh(sinhInput));
          break;

        case 'cosh':
          // Hyperbolic cosine: U-shaped transformation
          const coshInput = (normalizedFreq - 0.5) * amount;
          mappedFreq = (Math.cosh(coshInput) - 1) / (Math.cosh(amount) - 1);
          break;

        case 'tanh':
          // Hyperbolic tangent: saturation function
          const tanhInput = (normalizedFreq - 0.5) * amount * 3;
          mappedFreq = (Math.tanh(tanhInput) + 1) / 2;
          magnitudeMultiplier = 1 + 0.3 * Math.abs(Math.tanh(tanhInput));
          break;

        // === COMPLEX OPERATIONS ===
        case 'complex_magnitude':
          // Distance from origin with phase modulation
          const compReal = normalizedFreq * Math.cos(amount);
          const compImag = normalizedFreq * Math.sin(amount);
          mappedFreq = Math.sqrt(compReal * compReal + compImag * compImag);
          break;

        case 'complex_phase':
          // Phase-based frequency mapping
          const phaseAngle = normalizedFreq * Math.PI * amount;
          mappedFreq = (Math.atan2(Math.sin(phaseAngle), Math.cos(phaseAngle)) + Math.PI) / (2 * Math.PI);
          break;

        case 'conjugate':
          // Complex conjugate: frequency mirroring
          mappedFreq = 1 - normalizedFreq * (amount * 0.1);
          phaseShift = -phase * amount * 0.1;
          break;

        // === PHYSICS EQUATIONS ===
        case 'schrodinger':
          // Quantum harmonic oscillator
          const x_pos = (normalizedFreq - 0.5) * 6;
          const n_level = Math.floor(amount) % 10;
          const hermite = Math.exp(-x_pos * x_pos / 2) * Math.pow(x_pos, n_level);
          mappedFreq = normalizedFreq * (1 + 0.5 * hermite * (amount * 0.1));
          magnitudeMultiplier = 1 + 0.3 * Math.abs(hermite);
          break;

        case 'heisenberg':
          // Uncertainty principle: Δx·Δp ≥ ħ/2
          const uncertainty = amount / (normalizedFreq * 10 + 1);
          mappedFreq = normalizedFreq * (1 + uncertainty * 0.3);
          magnitudeMultiplier = 1 + uncertainty * 0.2;
          break;

        case 'einstein_mass_energy':
          // E = mc²: relativistic effects
          const velocity = normalizedFreq * 0.99;
          const gamma = 1 / Math.sqrt(1 - velocity * velocity);
          mappedFreq = normalizedFreq * (1 + (gamma - 1) * amount * 0.1);
          magnitudeMultiplier = Math.sqrt(gamma) * (1 + amount * 0.1);
          break;

        case 'planck_radiation':
          // Blackbody radiation curve
          const temp = 1000 + amount * 500; // Temperature
          const planckCurve = Math.pow(normalizedFreq, 3) / (Math.exp(normalizedFreq * 20 / temp) - 1);
          mappedFreq = normalizedFreq * (1 + planckCurve * amount * 0.1);
          magnitudeMultiplier = 1 + planckCurve * amount * 0.2;
          break;

        case 'wave_equation':
          // Standing wave patterns
          const modeNum = Math.floor(amount) % 20 + 1;
          const standingWave = Math.sin(normalizedFreq * Math.PI * modeNum);
          mappedFreq = normalizedFreq * (1 + standingWave * amount * 0.1);
          magnitudeMultiplier = 1 + Math.abs(standingWave) * amount * 0.3;
          break;

        // === FRACTALS ===
        case 'mandelbrot':
          // Mandelbrot set escape time
          let zr = (normalizedFreq - 0.5) * 3;
          let zi = 0;
          const cr = -0.5 + amount * 0.1;
          const ci = 0.6 * amount * 0.1;
          let iter = 0;
          const maxIter = 50;

          while (iter < maxIter && zr * zr + zi * zi < 4) {
            const temp = zr * zr - zi * zi + cr;
            zi = 2 * zr * zi + ci;
            zr = temp;
            iter++;
          }

          const escape = iter / maxIter;
          mappedFreq = normalizedFreq * (0.5 + escape * 0.5 * amount * 0.1);
          magnitudeMultiplier = 0.5 + escape * amount * 0.1;
          break;

        case 'julia':
          // Julia set
          const juliaC_r = -0.7 + amount * 0.05;
          const juliaC_i = 0.27 + amount * 0.05;
          let jr = (normalizedFreq - 0.5) * 2;
          let ji = 0.1 * amount;
          let juliaIter = 0;

          while (juliaIter < 30 && jr * jr + ji * ji < 4) {
            const temp = jr * jr - ji * ji + juliaC_r;
            ji = 2 * jr * ji + juliaC_i;
            jr = temp;
            juliaIter++;
          }

          const juliaEscape = juliaIter / 30;
          mappedFreq = normalizedFreq * (0.3 + juliaEscape * 0.7 * amount * 0.1);
          break;

        // === GEOMETRIC TRANSFORMATIONS ===
        case 'mobius':
          // Möbius transformation: z' = (az + b) / (cz + d)
          const a = 1 + amount * 0.1;
          const b = amount * 0.05;
          const c = amount * 0.03;
          const d = 1;
          const mobZ = normalizedFreq + 0.1;
          mappedFreq = Math.abs((a * mobZ + b) / (c * mobZ + d)) % 1;
          break;

        case 'joukowsky':
          // Joukowsky airfoil transformation
          const jouZ = normalizedFreq * 2 - 1; // -1 to 1
          const jouResult = 0.5 * (jouZ + amount * 0.2 / (jouZ + 0.01));
          mappedFreq = (jouResult + 1) / 2;
          mappedFreq = Math.max(0, Math.min(1, mappedFreq));
          break;

        case 'spiral':
          // Logarithmic spiral
          const spiralR = normalizedFreq;
          const spiralTheta = normalizedFreq * Math.PI * 2 * amount;
          const spiralX = spiralR * Math.cos(spiralTheta);
          mappedFreq = (spiralX + 1) / 2;
          break;

        // Add more functions...
        default:
          // For unimplemented functions, apply a default transformation
          mappedFreq = Math.pow(normalizedFreq, 1 + amount * 0.1);
          break;
      }

      // Clamp mapped frequency to valid range based on anti-aliasing mode
      let newBin;

      switch (realtimeParams.antiAliasingMode) {
        case 0: // Filter - clamp to 0-1 range
          mappedFreq = Math.max(0, Math.min(1, mappedFreq));
          // Convert back to actual frequency in the selected range
          const newFrequency = minFreq + mappedFreq * (maxFreq - minFreq);
          newBin = Math.round(newFrequency / binSize);
          break;

        case 1: // Wrap - wrap around frequencies (modulo)
          // Wrap mapped frequency to 0-1 range
          mappedFreq = ((mappedFreq % 1) + 1) % 1;
          // Convert back to actual frequency in the selected range
          const wrappedFrequency = minFreq + mappedFreq * (maxFreq - minFreq);
          newBin = Math.round(wrappedFrequency / binSize);
          break;

        case 2: // Off - allow any frequency (may cause aliasing)
          // Convert back to actual frequency, not limiting to max
          const unlimitedFrequency = minFreq + mappedFreq * (maxFreq - minFreq);
          newBin = Math.round(unlimitedFrequency / binSize);
          // Handle extreme frequencies with human range mapping if enabled
          if (realtimeParams.frequencyMapping === 1) {
            // Scale frequencies beyond hearing range back to audible range
            if (newBin >= fftSize / 2 || newBin < 0) {
              // Map large frequencies logarithmically back to 20Hz-20kHz
              const extremeFreq = unlimitedFrequency;
              const scaledFreq = 20 * Math.pow(1000, (Math.log(Math.abs(extremeFreq) / 20) / Math.log(1000)) % 1);
              newBin = Math.round(scaledFreq / binSize);
            }
          }
          break;
      }

      // Only place frequency content in valid bins
      if (newBin > 0 && newBin < fftSize / 2) {
        // Apply magnitude and phase modifications
        const finalMagnitude = magnitude * magnitudeMultiplier;
        const finalPhase = phase + phaseShift;

        // Convert back to real/imaginary
        real[newBin] += finalMagnitude * Math.cos(finalPhase);
        imag[newBin] += finalMagnitude * Math.sin(finalPhase);

        // Mirror for negative frequencies (maintain symmetry)
        const negBin = fftSize - newBin;
        if (negBin < fftSize) {
          real[negBin] = real[newBin];
          imag[negBin] = -imag[newBin];
        }
      }
    }

    // Preserve DC component
    real[0] = tempReal[0];
    imag[0] = 0;

  } catch (error) {
    console.error(`Error in applyComplexFunction (${funcName}):`, error);
    // On error, just copy original data
    for (let i = 0; i < fftSize; i++) {
      real[i] = tempReal[i];
      imag[i] = tempImag[i];
    }
  }
}

// STFT Implementation
function stft(audioData, fftSize, hopSize, windowFunc) {
  const numFrames = Math.floor((audioData.length - fftSize) / hopSize) + 1;
  const frames = [];

  for (let i = 0; i < numFrames; i++) {
    const start = i * hopSize;
    const frame = audioData.slice(start, start + fftSize);

    // Apply window
    for (let j = 0; j < fftSize; j++) {
      frame[j] *= windowFunc[j];
    }

    frames.push(frame);
  }

  return frames;
}

// Custom FFT Implementation
function fft(real, imag) {
  const n = real.length;
  if (n <= 1) return;

  // Ensure n is a power of 2
  if ((n & (n - 1)) !== 0) {
    console.error('FFT size must be a power of 2');
    return;
  }

  // Bit-reversal permutation
  for (let i = 0; i < n; i++) {
    let j = 0;
    let temp = i;
    for (let k = 0; k < Math.log2(n); k++) {
      j = (j << 1) | (temp & 1);
      temp >>= 1;
    }
    if (j > i) {
      [real[i], real[j]] = [real[j], real[i]];
      [imag[i], imag[j]] = [imag[j], imag[i]];
    }
  }

  // Cooley-Tukey FFT
  for (let len = 2; len <= n; len *= 2) {
    const angle = -2 * Math.PI / len;
    const wlen = { re: Math.cos(angle), im: Math.sin(angle) };

    for (let i = 0; i < n; i += len) {
      let w = { re: 1, im: 0 };

      for (let j = 0; j < len / 2; j++) {
        const u = { re: real[i + j], im: imag[i + j] };
        const v = {
          re: real[i + j + len / 2] * w.re - imag[i + j + len / 2] * w.im,
          im: real[i + j + len / 2] * w.im + imag[i + j + len / 2] * w.re
        };

        real[i + j] = u.re + v.re;
        imag[i + j] = u.im + v.im;
        real[i + j + len / 2] = u.re - v.re;
        imag[i + j + len / 2] = u.im - v.im;

        const newW = {
          re: w.re * wlen.re - w.im * wlen.im,
          im: w.re * wlen.im + w.im * wlen.re
        };
        w = newW;
      }
    }
  }
}

// Inverse FFT
function ifft(real, imag) {
  const n = real.length;

  // Conjugate
  for (let i = 0; i < n; i++) {
    imag[i] = -imag[i];
  }

  // Forward FFT
  fft(real, imag);

  // Conjugate and scale
  for (let i = 0; i < n; i++) {
    real[i] /= n;
    imag[i] = -imag[i] / n;
  }
}

// Main audio processing function with proper frequency transformations
async function processAudio() {
  if (!audioBuffer) {
    showStatus('No audio loaded', 'error');
    return null;
  }

  showStatus('Processing audio...', 'info');

  try {
    const fftSize = Math.pow(2, parseInt(document.getElementById('fftSizeSlider').value));
    const overlapFactor = Math.pow(2, parseInt(document.getElementById('overlapSlider').value));
    const hopSize = Math.floor(fftSize / overlapFactor);
    const amount = parseFloat(document.getElementById('powerSlider').value);
    const smoothing = parseFloat(document.getElementById('smoothingSlider').value);
    const mixValue = parseFloat(document.getElementById('mixSlider').value) / 100;
    const gain = parseFloat(document.getElementById('gainSlider').value);

    console.log(`Processing with FFT size: ${fftSize}, hop size: ${hopSize}, amount: ${amount}`);

    const audioData = audioBuffer.getChannelData(0);
    const sampleRate = audioBuffer.sampleRate;
    const windowFunc = getWindow(fftSize);

    // Perform STFT
    const frames = stft(audioData, fftSize, hopSize, windowFunc);
    const processedFrames = [];

    console.log(`Processing ${frames.length} frames`);

    for (let frameIdx = 0; frameIdx < frames.length; frameIdx++) {
      const frame = frames[frameIdx];
      const real = getBuffer(fftSize, 'real');
      const imag = getBuffer(fftSize, 'imag');

      // Copy frame data and zero-pad if necessary
      for (let i = 0; i < fftSize; i++) {
        real[i] = i < frame.length ? frame[i] : 0;
        imag[i] = 0;
      }

      // Forward FFT
      fft(real, imag);

      // Apply complex frequency transformations
      applyComplexFunction(real, imag, fftSize, sampleRate, selectedFunction, amount);

      // Apply smoothing in frequency domain
      if (smoothing > 0) {
        for (let i = 1; i < fftSize / 2 - 1; i++) {
          const prevReal = real[i - 1];
          const prevImag = imag[i - 1];
          const nextReal = real[i + 1];
          const nextImag = imag[i + 1];

          real[i] = real[i] * (1 - smoothing) + (prevReal + nextReal) * smoothing / 2;
          imag[i] = imag[i] * (1 - smoothing) + (prevImag + nextImag) * smoothing / 2;

          // Update mirrored bins
          real[fftSize - i] = real[i];
          imag[fftSize - i] = -imag[i];
        }
      }

      // Inverse FFT
      ifft(real, imag);

      // Store raw IFFT output (no additional windowing)
      processedFrames.push(real.slice(0, fftSize));
    }

    // Overlap-add reconstruction with proper normalization
    const outputLength = (frames.length - 1) * hopSize + fftSize;
    const outputBuffer = new Float32Array(outputLength);
    const normBuffer = new Float32Array(outputLength);

    // Calculate window normalization factor for overlap-add
    for (let i = 0; i < frames.length; i++) {
      const start = i * hopSize;
      for (let j = 0; j < fftSize; j++) {
        if (start + j < outputLength) {
          normBuffer[start + j] += windowFunc[j] * windowFunc[j];
        }
      }
    }

    // Reconstruct with proper windowing and normalization
    for (let i = 0; i < processedFrames.length; i++) {
      const start = i * hopSize;
      for (let j = 0; j < fftSize; j++) {
        if (start + j < outputLength) {
          // Apply synthesis window and accumulate
          outputBuffer[start + j] += processedFrames[i][j] * windowFunc[j];
        }
      }
    }

    // Normalize to prevent amplitude modulation
    for (let i = 0; i < outputLength; i++) {
      if (normBuffer[i] > 1e-12) {
        outputBuffer[i] /= normBuffer[i];
      }
    }

    // Mix with original signal
    const mixedBuffer = new Float32Array(Math.min(outputBuffer.length, audioData.length));
    for (let i = 0; i < mixedBuffer.length; i++) {
      const dry = audioData[i];
      const wet = outputBuffer[i] * gain;
      mixedBuffer[i] = dry * (1 - mixValue) + wet * mixValue;
    }

    // Create new audio buffer
    const newBuffer = ctx.createBuffer(1, mixedBuffer.length, sampleRate);
    newBuffer.copyToChannel(mixedBuffer, 0);

    showStatus('Audio processing complete!', 'success');
    updateMemoryInfo();

    return newBuffer;

  } catch (error) {
    showStatus(`Error processing audio: ${error.message}`, 'error');
    console.error('Audio processing error:', error);
    return null;
  }
}

// Visualization
function setupVisualization() {
  canvas = document.getElementById('visualizerCanvas');
  canvasCtx = canvas.getContext('2d');
  canvas.width = canvas.offsetWidth;
  canvas.height = canvas.offsetHeight;
}

function drawVisualization() {
  if (!analyser) return;

  const bufferLength = analyser.frequencyBinCount;
  const dataArray = new Uint8Array(bufferLength);
  const timeArray = new Uint8Array(analyser.fftSize);

  analyser.getByteFrequencyData(dataArray);
  analyser.getByteTimeDomainData(timeArray);

  // Clear canvas with fade effect
  canvasCtx.fillStyle = 'rgba(10, 10, 10, 0.1)';
  canvasCtx.fillRect(0, 0, canvas.width, canvas.height);

  // Draw frequency spectrum (top half)
  const freqHeight = canvas.height / 2;
  const barWidth = canvas.width / bufferLength;
  let x = 0;

  for (let i = 0; i < bufferLength; i++) {
    const barHeight = (dataArray[i] / 255) * freqHeight;

    // Color based on frequency and transformation
    const hue = (i / bufferLength) * 360;
    const saturation = 70 + (dataArray[i] / 255) * 30;
    const lightness = 40 + (dataArray[i] / 255) * 20;

    canvasCtx.fillStyle = `hsl(${hue}, ${saturation}%, ${lightness}%)`;
    canvasCtx.fillRect(x, freqHeight - barHeight, barWidth, barHeight);

    x += barWidth;
  }

  // Draw waveform (bottom half)
  canvasCtx.strokeStyle = '#4ecdc4';
  canvasCtx.lineWidth = 2;
  canvasCtx.beginPath();

  const sliceWidth = canvas.width / timeArray.length;
  x = 0;

  for (let i = 0; i < timeArray.length; i++) {
    const v = (timeArray[i] - 128) / 128;
    const y = freqHeight + (v * freqHeight / 2) + freqHeight / 2;

    if (i === 0) {
      canvasCtx.moveTo(x, y);
    } else {
      canvasCtx.lineTo(x, y);
    }

    x += sliceWidth;
  }

  canvasCtx.stroke();

  // Draw function name and parameters
  canvasCtx.fillStyle = 'rgba(255, 255, 255, 0.8)';
  canvasCtx.font = '16px Segoe UI';
  canvasCtx.fillText(`Function: ${complexFunctions[selectedFunction].name}`, 10, 25);

  if (isPlaying) {
    canvasCtx.fillStyle = 'rgba(78, 205, 196, 0.8)';
    canvasCtx.fillText('● REAL-TIME', 10, 50);
  }

  if (isPlaying) {
    animationId = requestAnimationFrame(drawVisualization);
  }
}

// Robust real-time audio processor with proper overlap-add
function setupRealtimeAudio() {
  if (!ctx || !audioBuffer) return;

  try {
    // Create script processor for real-time effects
    scriptProcessor = ctx.createScriptProcessor(4096, 1, 1);
    analyser = ctx.createAnalyser();
    analyser.fftSize = 2048;
    analyser.smoothingTimeConstant = 0.3;

    // FFT processing setup
    const fftSize = 2048;
    const hopSize = fftSize / 4; // 75% overlap
    const windowFunc = createOptimalHannWindow(fftSize, hopSize);

    // Delay line for overlap-add
    let inputDelay = new Float32Array(fftSize);
    let outputDelay = new Float32Array(fftSize);
    let frameCounter = 0;

    // Initialize processor parameters
    scriptProcessor.params = { ...realtimeParams };
    // Set global reference to processor
    realtimeProcessor = scriptProcessor;

    scriptProcessor.onaudioprocess = function (audioProcessingEvent) {
      const input = audioProcessingEvent.inputBuffer.getChannelData(0);
      const output = audioProcessingEvent.outputBuffer.getChannelData(0);

      // Get current parameters from the processor object
      let amount = scriptProcessor.params.amount;
      const mixValue = scriptProcessor.params.mix;
      const gain = scriptProcessor.params.gain;

      // Apply modulation to amount if enabled
      if (scriptProcessor.params.modEnabled) {
        const minAmount = scriptProcessor.params.modMinAmount || 0.0;
        const maxAmount = scriptProcessor.params.modMaxAmount || amount;
        const currentTime = ctx.currentTime;
        const startTime = scriptProcessor.params.modStartTime;
        const modPeriod = scriptProcessor.params.modPeriod;

        // Calculate elapsed time and position in modulation cycle
        const elapsedTime = currentTime - startTime;
        const normalizedPosition = (elapsedTime % modPeriod) / modPeriod;

        // Sinusoidal modulation between min and max bounds
        const modulationFactor = Math.sin(normalizedPosition * 2 * Math.PI);
        amount = minAmount + (maxAmount - minAmount) * (modulationFactor * 0.5 + 0.5);
      }

      for (let n = 0; n < input.length; n++) {
        // Shift delay lines
        for (let i = 0; i < fftSize - 1; i++) {
          inputDelay[i] = inputDelay[i + 1];
          outputDelay[i] = outputDelay[i + 1];
        }
        inputDelay[fftSize - 1] = input[n];
        outputDelay[fftSize - 1] = 0;

        // Process frame every hopSize samples
        if (frameCounter % hopSize === 0) {
          // Prepare FFT buffers
          const real = new Float32Array(fftSize);
          const imag = new Float32Array(fftSize);

          // Apply window to input
          for (let i = 0; i < fftSize; i++) {
            real[i] = inputDelay[i] * windowFunc[i];
            imag[i] = 0;
          }

          // Forward FFT
          fft(real, imag);

          // Apply frequency transformation
          applyComplexFunction(real, imag, fftSize, ctx.sampleRate, selectedFunction, amount);

          // Inverse FFT
          ifft(real, imag);

          // Apply synthesis window and overlap-add
          for (let i = 0; i < fftSize; i++) {
            outputDelay[i] += real[i] * windowFunc[i];
          }
        }

        // Output with proper mix
        const wet = outputDelay[0] * gain;
        output[n] = input[n] * (1 - mixValue) + wet * mixValue;

        frameCounter++;
      }
    };

    return scriptProcessor;

  } catch (error) {
    console.error('Error setting up real-time audio:', error);
    return null;
  }
}

function createOptimalHannWindow(size, hopSize) {
  const window = new Float32Array(size);

  // Create basic Hann window
  for (let i = 0; i < size; i++) {
    window[i] = 0.5 * (1 - Math.cos(2 * Math.PI * i / (size - 1)));
  }

  // Calculate the optimal normalization for constant amplitude with overlap-add
  // For Hann window with 75% overlap (hopSize = size/4), we need to compensate
  // for the overlapping window energy

  // Test the overlap sum at the center of the window
  let overlapSum = 0;
  const center = Math.floor(size / 2);

  // Sum all overlapping windows at the center point
  for (let offset = -size; offset <= size; offset += hopSize) {
    const index = center + offset;
    if (index >= 0 && index < size) {
      overlapSum += window[index] * window[index];
    }
  }

  // Normalize so that overlapping windows sum to constant amplitude
  const normFactor = 1.0 / Math.sqrt(overlapSum);
  for (let i = 0; i < size; i++) {
    window[i] *= normFactor;
  }

  return window;
}

// Hann window creation for analysis
function createHannWindow(size) {
  const window = new Float32Array(size);
  for (let i = 0; i < size; i++) {
    window[i] = 0.5 * (1 - Math.cos(2 * Math.PI * i / (size - 1)));
  }

  // Normalize for overlap-add with 75% overlap (hopSize = size/4)
  // This ensures constant amplitude reconstruction
  let sum = 0;
  const hopSize = size / 4;
  for (let i = 0; i < size; i++) {
    let overlap = 0;
    for (let j = -3; j <= 3; j++) { // Check overlapping windows
      const pos = i + j * hopSize;
      if (pos >= 0 && pos < size) {
        overlap += window[pos] * window[pos];
      }
    }
    if (overlap > sum) sum = overlap;
  }

  // Apply normalization factor
  const normFactor = 1.0 / Math.sqrt(sum);
  for (let i = 0; i < size; i++) {
    window[i] *= normFactor;
  }

  return window;
}

function createSimpleHannWindow(size) {
  const window = new Float32Array(size);
  for (let i = 0; i < size; i++) {
    window[i] = 0.5 * (1 - Math.cos(2 * Math.PI * i / (size - 1)));
  }

  // Simple normalization for 75% overlap
  // For Hann window with 75% overlap, the normalization factor is approximately 2/3
  const normFactor = 2.0 / 3.0;
  for (let i = 0; i < size; i++) {
    window[i] *= normFactor;
  }

  return window;
}

// Update real-time amount display with modulated value
function updateModulatedDisplay() {
  if (!isPlaying || !realtimeProcessor || !realtimeProcessor.params.modEnabled) {
    // If modulation is off, restore original amount display
    if (realtimeProcessor && realtimeProcessor.params) {
      const powerValue = document.getElementById('powerValue');
      if (powerValue) {
        powerValue.textContent = realtimeProcessor.params.amount.toFixed(1);
      }
    }
    return;
  }

  const powerValue = document.getElementById('powerValue');
  if (!powerValue) return;

  const amount = realtimeProcessor.params.amount;
  const minAmount = realtimeProcessor.params.modMinAmount || 0.0;
  const maxAmount = realtimeProcessor.params.modMaxAmount || amount;
  const currentTime = ctx.currentTime;
  const startTime = realtimeProcessor.params.modStartTime;
  const modPeriod = realtimeProcessor.params.modPeriod;

  // Calculate elapsed time and position in modulation cycle
  const elapsedTime = currentTime - startTime;
  const normalizedPosition = (elapsedTime % modPeriod) / modPeriod;

  // Sinusoidal modulation between min and max bounds
  const modulationFactor = Math.sin(normalizedPosition * 2 * Math.PI);
  const modulatedAmount = minAmount + (maxAmount - minAmount) * (modulationFactor * 0.5 + 0.5);

  // Update display with the current modulated value and bounds
  powerValue.textContent = `${modulatedAmount.toFixed(3)} (${minAmount.toFixed(3)}→${maxAmount.toFixed(3)})`;

  // Request next animation frame
  requestAnimationFrame(updateModulatedDisplay);
}

// Setup real-time parameter updates
function setupRealtimeParameterUpdates() {
  // This function can be empty for bypass mode
  // In full mode, it would update FFT processing parameters in real-time
}

// Add toggle button handlers
function setupToggleButtons() {
  // Anti-aliasing mode toggles
  const filterToggle = document.getElementById('filterToggle');
  const wrapToggle = document.getElementById('wrapToggle');
  const noFilterToggle = document.getElementById('noFilterToggle');

  const antiAliasingBtns = [filterToggle, wrapToggle, noFilterToggle];

  antiAliasingBtns.forEach((btn, index) => {
    btn.addEventListener('click', () => {
      antiAliasingBtns.forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      realtimeParams.antiAliasingMode = index;
      console.log(`Anti-aliasing mode set to: ${index}`);
    });
  });

  // Frequency mapping toggles
  const normalFreqBtn = document.getElementById('normalFreqBtn');
  const humanRangeBtn = document.getElementById('humanRangeBtn');

  const freqMapBtns = [normalFreqBtn, humanRangeBtn];

  freqMapBtns.forEach((btn, index) => {
    btn.addEventListener('click', () => {
      freqMapBtns.forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      realtimeParams.frequencyMapping = index;
      console.log(`Frequency mapping set to: ${index}`);
    });
  });
}

// Event listeners and initialization
document.addEventListener('DOMContentLoaded', function () {
  // Initialize function grid
  function initFunctionSelect() {
    const select = document.getElementById('functionSelect');

    // Set initial value
    select.value = selectedFunction;

    // Add event listener
    select.addEventListener('change', function () {
      selectedFunction = this.value;
      showStatus(`Function changed to: ${complexFunctions[selectedFunction].name}`, 'info');
    });
  }

  function initSliders() {
    const sliders = {
      rangeSlider: (val) => {
        const values = ['Normal', 'Extended', 'Unlimited'];
        document.getElementById('rangeValue').textContent = values[val];
        realtimeParams.frequencyRange = parseInt(val);
        if (realtimeProcessor && realtimeProcessor.params) {
          realtimeProcessor.params.frequencyRange = parseInt(val);
        }
      },
      powerSlider: (val) => {
        document.getElementById('powerValue').textContent = parseFloat(val).toFixed(3);
        realtimeParams.amount = parseFloat(val);
        if (realtimeProcessor && realtimeProcessor.params) {
          realtimeProcessor.params.amount = parseFloat(val);
        }
      },
      modPeriodSlider: (val) => {
        const value = parseFloat(val);
        if (value < 0) {
          // Modulation off
          document.getElementById('modStatus').textContent = 'Off';
          document.getElementById('modPeriodValue').textContent = 'Off';
          realtimeParams.modEnabled = false;
          if (realtimeProcessor && realtimeProcessor.params) {
            realtimeProcessor.params.modEnabled = false;
          }

          // Reset amount display to show normal value
          if (realtimeProcessor && realtimeProcessor.params) {
            const powerValue = document.getElementById('powerValue');
            if (powerValue) {
              powerValue.textContent = realtimeProcessor.params.amount.toFixed(3);
            }
          }
        } else {
          // Convert logarithmic slider value to period (0.1s to 10m)
          // value range is 0-3, mapping to log10(0.1) to log10(600)
          const minLog = Math.log10(0.1);
          const maxLog = Math.log10(600); // 10 minutes in seconds
          const logPeriod = minLog + (value * (maxLog - minLog) / 3);
          const period = Math.pow(10, logPeriod);

          // Format display based on duration
          let displayValue;
          if (period < 60) {
            displayValue = `${period.toFixed(1)}s`;
          } else {
            const minutes = Math.floor(period / 60);
            const seconds = Math.floor(period % 60);
            displayValue = `${minutes}m ${seconds}s`;
          }

          document.getElementById('modStatus').textContent = 'On';
          document.getElementById('modPeriodValue').textContent = displayValue;

          realtimeParams.modEnabled = true;
          realtimeParams.modPeriod = period;
          realtimeParams.modStartTime = ctx ? ctx.currentTime : 0;

          if (realtimeProcessor && realtimeProcessor.params) {
            realtimeProcessor.params.modEnabled = true;
            realtimeProcessor.params.modPeriod = period;
            realtimeProcessor.params.modStartTime = ctx ? ctx.currentTime : 0;
          }

          // Start modulation display updates if playing
          if (isPlaying) {
            requestAnimationFrame(updateModulatedDisplay);
          }
        }
      },
      mixSlider: (val) => {
        document.getElementById('mixValue').textContent = `${val}%`;
        realtimeParams.mix = parseInt(val) / 100;
        if (realtimeProcessor && realtimeProcessor.params) {
          realtimeProcessor.params.mix = parseInt(val) / 100;
        }
      },
      gainSlider: (val) => {
        document.getElementById('gainValue').textContent = parseFloat(val).toFixed(1);
        realtimeParams.gain = parseFloat(val);
        if (realtimeProcessor && realtimeProcessor.params) {
          realtimeProcessor.params.gain = parseFloat(val);
        }
      }
    };

    Object.entries(sliders).forEach(([id, callback]) => {
      const slider = document.getElementById(id);
      if (slider) {
        slider.addEventListener('input', (e) => callback(e.target.value));
        // Initialize display without triggering processor updates
        const currentValue = slider.value;
        if (id === 'rangeSlider') {
          const values = ['Normal', 'Extended', 'Unlimited'];
          document.getElementById('rangeValue').textContent = values[currentValue];
          realtimeParams.frequencyRange = parseInt(currentValue);
        } else if (id === 'powerSlider') {
          document.getElementById('powerValue').textContent = parseFloat(currentValue).toFixed(3);
          realtimeParams.amount = parseFloat(currentValue);
        } else if (id === 'modPeriodSlider') {
          // Initialize modulation slider display
          const value = parseFloat(currentValue);
          if (value < 0) {
            document.getElementById('modStatus').textContent = 'Off';
            document.getElementById('modPeriodValue').textContent = 'Off';
            realtimeParams.modEnabled = false;
          } else {
            const minLog = Math.log10(0.1);
            const maxLog = Math.log10(600);
            const logPeriod = minLog + (value * (maxLog - minLog) / 3);
            const period = Math.pow(10, logPeriod);

            let displayValue;
            if (period < 60) {
              displayValue = `${period.toFixed(1)}s`;
            } else {
              const minutes = Math.floor(period / 60);
              const seconds = Math.floor(period % 60);
              displayValue = `${minutes}m ${seconds}s`;
            }

            document.getElementById('modStatus').textContent = 'On';
            document.getElementById('modPeriodValue').textContent = displayValue;
            realtimeParams.modEnabled = true;
            realtimeParams.modPeriod = period;
          }
        } else if (id === 'mixSlider') {
          document.getElementById('mixValue').textContent = `${currentValue}%`;
          realtimeParams.mix = parseInt(currentValue) / 100;
        } else if (id === 'gainSlider') {
          document.getElementById('gainValue').textContent = parseFloat(currentValue).toFixed(1);
          realtimeParams.gain = parseFloat(currentValue);
        }
      }
    });

    // Initialize dual-handle modulation bounds sliders
    const modMinSlider = document.getElementById('modMinSlider');
    const modMaxSlider = document.getElementById('modMaxSlider');
    const modRangeValue = document.getElementById('modRangeValue');
    const dualSliderContainer = document.querySelector('.dual-slider-container');

    if (modMinSlider && modMaxSlider && modRangeValue) {
      function updateModulationBounds() {
        const minVal = parseFloat(modMinSlider.value);
        const maxVal = parseFloat(modMaxSlider.value);

        // Ensure min is always less than max
        if (minVal >= maxVal) {
          if (this === modMinSlider) {
            modMaxSlider.value = (minVal + 0.001).toFixed(3);
          } else {
            modMinSlider.value = (maxVal - 0.001).toFixed(3);
          }
        }

        const finalMin = parseFloat(modMinSlider.value);
        const finalMax = parseFloat(modMaxSlider.value);

        // Update display
        modRangeValue.textContent = `${finalMin.toFixed(3)} → ${finalMax.toFixed(3)}`;

        // Update parameters
        realtimeParams.modMinAmount = finalMin;
        realtimeParams.modMaxAmount = finalMax;
        if (realtimeProcessor && realtimeProcessor.params) {
          realtimeProcessor.params.modMinAmount = finalMin;
          realtimeProcessor.params.modMaxAmount = finalMax;
        }
      }

      modMinSlider.addEventListener('input', updateModulationBounds);
      modMaxSlider.addEventListener('input', updateModulationBounds);

      // Initialize display
      updateModulationBounds();

      // Toggle visibility based on modulation state
      function updateModulationVisibility() {
        if (realtimeParams.modEnabled) {
          dualSliderContainer.classList.remove('modulation-disabled');
        } else {
          dualSliderContainer.classList.add('modulation-disabled');
        }
      }

      // Initial visibility update
      updateModulationVisibility();

      // Update visibility when modulation period changes
      const modPeriodSlider = document.getElementById('modPeriodSlider');
      if (modPeriodSlider) {
        modPeriodSlider.addEventListener('input', () => {
          setTimeout(updateModulationVisibility, 0); // Update after the main callback
        });
      }
    }
  }

  // Enhanced file input with drag-and-drop support
  function setupFileInput() {
    const fileInput = document.getElementById('fileInput');
    const dropZone = document.getElementById('fileDropZone');

    // Handle file selection (both click and drag-drop)
    async function handleFile(file) {
      if (!file) return;

      showStatus(`📁 Loading: ${file.name} (${(file.size / 1024 / 1024).toFixed(1)}MB)...`, 'info');
      console.log(`Loading file: ${file.name}, size: ${file.size} bytes, type: ${file.type}`);

      try {
        // Enhanced file type checking - check both MIME type and extension
        const supportedExtensions = ['wav', 'mp3', 'ogg', 'webm', 'm4a', 'aac', 'flac', 'opus'];
        const supportedMimeTypes = [
          'audio/wav', 'audio/wave', 'audio/x-wav',
          'audio/mp3', 'audio/mpeg', 'audio/mp4',
          'audio/ogg', 'audio/opus',
          'audio/webm', 'audio/m4a', 'audio/aac',
          'audio/flac', 'audio/x-flac'
        ];

        const fileExtension = file.name.toLowerCase().split('.').pop();
        const mimeType = file.type.toLowerCase();

        const isValidExtension = supportedExtensions.includes(fileExtension);
        const isValidMimeType = supportedMimeTypes.some(type => mimeType.includes(type.replace('audio/', '')));

        if (!isValidExtension && !isValidMimeType) {
          throw new Error(`Unsupported file format: ${fileExtension || mimeType}.\nSupported formats: ${supportedExtensions.join(', ').toUpperCase()}`);
        }

        // Initialize audio context first
        await initAudioContext();
        console.log('Audio context initialized, reading file...');
        showStatus(`🔄 Reading file data...`, 'info');

        // Read file as array buffer with progress
        const arrayBuffer = await file.arrayBuffer();
        console.log(`File read complete, buffer size: ${arrayBuffer.byteLength} bytes`);
        showStatus(`🎵 Decoding audio (${fileExtension?.toUpperCase() || 'unknown format'})...`, 'info');

        // Decode audio data with comprehensive error handling
        try {
          audioBuffer = await ctx.decodeAudioData(arrayBuffer);

          const duration = audioBuffer.duration;
          const sampleRate = audioBuffer.sampleRate;
          const channels = audioBuffer.numberOfChannels;
          const samples = audioBuffer.length;

          console.log(`Audio decoded successfully:
            - File: ${file.name}
            - Duration: ${duration.toFixed(2)}s
            - Sample Rate: ${sampleRate}Hz
            - Channels: ${channels} (${channels === 1 ? 'Mono' : channels === 2 ? 'Stereo' : channels + '-channel'})
            - Samples: ${samples.toLocaleString()}
            - Size: ${(file.size / 1024 / 1024).toFixed(1)}MB`);

          showStatus(`✅ Ready: ${file.name} | ${duration.toFixed(1)}s | ${sampleRate}Hz | ${channels}ch`, 'success');
          document.getElementById('playBtn').disabled = false;

          // Update drop zone text to show loaded file
          const dropText = dropZone.querySelector('.file-drop-text');
          dropText.textContent = `✅ ${file.name} (${duration.toFixed(1)}s)`;
          dropZone.style.borderColor = 'rgba(78, 205, 196, 0.7)';
          dropZone.style.background = 'rgba(78, 205, 196, 0.1)';

          // Show audio info in status for user
          setTimeout(() => {
            showStatus(`🎵 ${file.name} loaded successfully (${duration.toFixed(1)}s, ${channels === 1 ? 'Mono' : 'Stereo'})`, 'success');
          }, 2000);

        } catch (decodeError) {
          console.error('Audio decoding failed:', decodeError);

          // Provide specific error messages for common issues
          let errorMessage = 'Failed to decode audio file.';
          if (decodeError.message.includes('Unable to decode')) {
            errorMessage = `Cannot decode ${fileExtension?.toUpperCase() || 'this'} file. It may be corrupted, DRM-protected, or in an unsupported codec.`;
          } else if (decodeError.message.includes('InvalidStateError')) {
            errorMessage = 'Audio context error. Try refreshing the page and loading the file again.';
          } else if (decodeError.message.includes('DataError')) {
            errorMessage = 'Invalid audio data. The file may be corrupted or not a valid audio file.';
          }

          throw new Error(`${errorMessage}\n\nTechnical details: ${decodeError.message}`);
        }

      } catch (error) {
        const errorMsg = error.message.replace('Error: ', '');
        showStatus(`❌ ${errorMsg}`, 'error');
        console.error('Audio loading error:', error);
        console.error('File details:', {
          name: file?.name,
          size: file?.size,
          type: file?.type,
          lastModified: file?.lastModified ? new Date(file.lastModified).toISOString() : 'unknown'
        });

        // Reset file input and disable play button
        fileInput.value = '';
        document.getElementById('playBtn').disabled = true;
        audioBuffer = null;

        // Reset drop zone appearance
        const dropText = dropZone.querySelector('.file-drop-text');
        dropText.textContent = '🎵 Drop audio file here or click to browse';
        dropZone.style.borderColor = 'rgba(255, 255, 255, 0.3)';
        dropZone.style.background = 'rgba(255, 255, 255, 0.05)';
      }
    }

    // File input change event
    fileInput.addEventListener('change', async function (e) {
      const file = e.target.files[0];
      await handleFile(file);
    });

    // Drag and drop events
    dropZone.addEventListener('dragover', (e) => {
      e.preventDefault();
      dropZone.classList.add('dragover');
    });

    dropZone.addEventListener('dragleave', (e) => {
      e.preventDefault();
      dropZone.classList.remove('dragover');
    });

    dropZone.addEventListener('drop', async (e) => {
      e.preventDefault();
      dropZone.classList.remove('dragover');

      const files = e.dataTransfer.files;
      if (files.length > 0) {
        const file = files[0];
        // Also update the file input to reflect the dropped file
        const dt = new DataTransfer();
        dt.items.add(file);
        fileInput.files = dt.files;
        await handleFile(file);
      }
    });

    // Click on drop zone to trigger file input
    dropZone.addEventListener('click', () => {
      fileInput.click();
    });
  }

  // File input handler (legacy - kept for compatibility)
  document.getElementById('fileInput').addEventListener('change', async function (e) {
    // This is now handled by setupFileInput() but kept for safety
    // The new handler will take precedence
  });

  // Play/Stop button (single button that toggles with real-time processing)
  document.getElementById('playBtn').addEventListener('click', async function () {
    const playBtn = document.getElementById('playBtn');

    // If currently playing, stop
    if (isPlaying) {
      // Stop real-time processing
      if (currentSource) {
        currentSource.stop();
        currentSource = null;
      }
      if (realtimeProcessor) {
        realtimeProcessor.disconnect();
        realtimeProcessor = null;
      }
      if (analyser) {
        analyser.disconnect();
      }

      isPlaying = false;
      isRealTimeMode = false;
      playBtn.textContent = '▶️ Play';
      cancelAnimationFrame(animationId);
      showStatus('Playback stopped', 'info');
      return;
    }

    // If not playing, start real-time processing
    if (!audioBuffer) {
      showStatus('Please load an audio file first', 'error');
      return;
    }

    try {
      await initAudioContext();

      playBtn.textContent = '⏸️ Starting...';

      // Create real-time processor
      realtimeProcessor = setupRealtimeAudio();
      if (!realtimeProcessor) {
        throw new Error('Failed to create real-time processor');
      }

      // Connect audio chain for real-time processing
      const source = ctx.createBufferSource();
      source.buffer = audioBuffer;
      source.loop = true;

      // Connect: source -> realtimeProcessor -> analyser -> destination
      source.connect(realtimeProcessor);
      realtimeProcessor.connect(analyser);
      analyser.connect(ctx.destination);

      currentSource = source;
      isPlaying = true;
      isRealTimeMode = true;
      playBtn.textContent = '⏹️ Stop';

      // Start visualization and modulation display updates
      animationId = requestAnimationFrame(drawVisualization);
      if (realtimeParams.modEnabled) {
        requestAnimationFrame(updateModulatedDisplay);
      }

      source.onended = () => {
        if (isPlaying) { // Only if not manually stopped
          isPlaying = false;
          isRealTimeMode = false;
          playBtn.textContent = '▶️ Play';
          currentSource = null;
          if (realtimeProcessor) {
            realtimeProcessor.disconnect();
            realtimeProcessor = null;
          }
          if (analyser) {
            analyser.disconnect();
          }
          cancelAnimationFrame(animationId);
          showStatus('Playback finished', 'info');
        }
      };

      source.start();

      // Setup real-time parameter updates
      setupRealtimeParameterUpdates();

      // Start visualization
      setupVisualization();
      drawVisualization();

      document.getElementById('downloadBtn').disabled = false;
      showStatus(`🎵 ULTIMATE FFT PROCESSING: ${Object.keys(complexFunctions).length} mathematical functions available!`, 'success');

      // Start updating the modulated display
      updateModulatedDisplay();

    } catch (error) {
      isPlaying = false;
      isRealTimeMode = false;
      playBtn.textContent = '▶️ Play';
      showStatus(`Error during playback: ${error.message}`, 'error');
      console.error('Playback error:', error);
    }
  });

  // Download button
  document.getElementById('downloadBtn').addEventListener('click', function () {
    if (!processedBuffer) {
      showStatus('No processed audio to download', 'error');
      return;
    }

    try {
      const channelData = processedBuffer.getChannelData(0);
      const length = channelData.length;
      const buffer = new ArrayBuffer(44 + length * 2);
      const view = new DataView(buffer);

      // WAV header
      const writeString = (offset, string) => {
        for (let i = 0; i < string.length; i++) {
          view.setUint8(offset + i, string.charCodeAt(i));
        }
      };

      writeString(0, 'RIFF');
      view.setUint32(4, 36 + length * 2, true);
      writeString(8, 'WAVE');
      writeString(12, 'fmt ');
      view.setUint32(16, 16, true);
      view.setUint16(20, 1, true);
      view.setUint16(22, 1, true);
      view.setUint32(24, processedBuffer.sampleRate, true);
      view.setUint32(28, processedBuffer.sampleRate * 2, true);
      view.setUint16(32, 2, true);
      view.setUint16(34, 16, true);
      writeString(36, 'data');
      view.setUint32(40, length * 2, true);

      // Convert float to 16-bit PCM
      let offset = 44;
      for (let i = 0; i < length; i++) {
        const sample = Math.max(-1, Math.min(1, channelData[i]));
        view.setInt16(offset, sample < 0 ? sample * 0x8000 : sample * 0x7FFF, true);
        offset += 2;
      }

      // Create download
      const blob = new Blob([buffer], { type: 'audio/wav' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `transformed_${selectedFunction}_${Date.now()}.wav`;
      a.click();
      URL.revokeObjectURL(url);

      showStatus('Audio downloaded successfully!', 'success');

    } catch (error) {
      showStatus(`Error downloading audio: ${error.message}`, 'error');
      console.error('Download error:', error);
    }
  });

  // Test audio button
  document.getElementById('testAudioBtn').addEventListener('click', async function () {
    try {
      await initAudioContext();

      showStatus('Generating test audio...', 'info');
      console.log('Generating test audio with enhanced harmonics...');

      const duration = 10; // 10 seconds
      const sampleRate = ctx.sampleRate;
      const length = duration * sampleRate;

      console.log(`Creating buffer: ${length} samples at ${sampleRate}Hz`);

      audioBuffer = ctx.createBuffer(1, length, sampleRate);
      const channelData = audioBuffer.getChannelData(0);

      // Create a rich harmonic test signal
      const fundamentalFreq = 220; // A3
      const amplitude = 0.3;

      for (let i = 0; i < length; i++) {
        const t = i / sampleRate;
        let sample = 0;

        // Add multiple harmonics for a rich sound
        sample += Math.sin(2 * Math.PI * fundamentalFreq * t) * amplitude; // Fundamental
        sample += Math.sin(2 * Math.PI * fundamentalFreq * 2 * t) * amplitude * 0.5; // 2nd harmonic
        sample += Math.sin(2 * Math.PI * fundamentalFreq * 3 * t) * amplitude * 0.3; // 3rd harmonic
        sample += Math.sin(2 * Math.PI * fundamentalFreq * 4 * t) * amplitude * 0.2; // 4th harmonic
        sample += Math.sin(2 * Math.PI * fundamentalFreq * 5 * t) * amplitude * 0.1; // 5th harmonic

        // Add some gentle frequency modulation
        const fmAmount = 5; // Hz modulation depth
        const fmRate = 2; // Hz modulation rate
        const fm = Math.sin(2 * Math.PI * fmRate * t) * fmAmount;
        sample += Math.sin(2 * Math.PI * (fundamentalFreq + fm) * t) * amplitude * 0.1;

        channelData[i] = sample;
      }

      showStatus('Test audio generated! Click Play to start real-time processing.', 'success');
      document.getElementById('playBtn').disabled = false;
      console.log('Test audio generation complete');

    } catch (error) {
      showStatus(`Error generating test audio: ${error.message}`, 'error');
      console.error('Test audio generation error:', error);
      console.error('Error details:', {
        name: error.name,
        message: error.message,
        stack: error.stack
      });
    }
  });

  // Global error handler for all audio operations
  window.addEventListener('error', function (event) {
    if (event.error && event.error.message &&
      (event.error.message.includes('audio') ||
        event.error.message.includes('AudioContext') ||
        event.error.message.includes('TypeError'))) {
      console.error('Audio error caught by global handler:', event.error);
      showStatus(`Error: ${event.error.message}. Try clicking the play button first.`, 'error');
      event.preventDefault();
    }
  });

  // Initialize everything
  initFunctionSelect();
  initSliders();
  setupFileInput();
  setupVisualization();
  setupToggleButtons();

  showStatus('Ultimate Complex Audio Transformer ready! Load an audio file or generate test audio.', 'info');
  updateMemoryInfo();
});
