# 🎧 Ultimate Complex Audio Transformer

A sophisticated web-based audio processing application that applies complex mathematical transformations to audio signals in real-time. Transform your audio through fractal mathematics, transcendental functions, and geometric transformations.

## ✨ Features

### 🎵 Audio Processing
- **Real-time audio transformation** using Web Audio API
- **STFT (Short-Time Fourier Transform)** for frequency domain analysis
- **Custom FFT implementation** optimized for audio processing
- **Overlap-add reconstruction** for seamless audio output
- **Memory-efficient buffer management** with pooling system

### 🔢 Mathematical Functions
Over 50 complex mathematical transformations including:

#### Basic Operations
- Power functions (f^n)
- Square root, cube root
- Logarithmic and exponential functions
- Trigonometric functions (sin, cos, tan, sinh, cosh, tanh)

#### Fractal Functions
- **Mandelbrot Set** - Classic fractal iteration
- **Julia Sets** - Dynamic fractal transformations
- **Burning Ship** - Unique fractal variant
- **Phoenix, Tricorn, Multibrot** - Advanced fractal types
- **Newton, Nova, Magnet** - Mathematical convergence patterns

#### Geometric Transformations
- **Möbius transformations** - Conformal mappings
- **Joukowsky transformation** - Airfoil mathematics
- **Spiral, Wave, Ripple** - Dynamic distortions
- **Swirl, Pinch, Bulge** - Geometric effects
- **Kaleidoscope** - Symmetrical patterns

#### Advanced Functions
- **Bessel functions** - Cylindrical wave solutions
- **Gamma function** - Generalized factorial
- **Riemann Zeta** - Number theory applications
- **Hypergeometric functions** - Special mathematical functions

### 🎛️ Controls
- **FFT Size**: 512 to 16384 samples (adjustable resolution)
- **Overlap Factor**: 2x to 8x (seamless processing)
- **Frequency Power**: f^0.1 to f^5.0 (harmonic shaping)
- **Smoothing**: 0% to 100% (artifact reduction)
- **Dry/Wet Mix**: Original to fully processed
- **Output Gain**: -20dB to +10dB (volume control)

### 📊 Visualization
- **Real-time frequency spectrum** display
- **Color-coded frequency bands** for visual feedback
- **Smooth animation** with requestAnimationFrame
- **Responsive design** for all screen sizes

### 💾 Export Features
- **WAV file download** of processed audio
- **16-bit PCM encoding** for compatibility
- **Automatic filename generation** with function name

## 🚀 Getting Started

### Prerequisites
- Modern web browser with Web Audio API support
- Local web server (for file loading security)

### Installation
1. Clone or download the project files
2. Serve the files using a local web server:
   ```bash
   # Python 3
   python -m http.server 8000
   
   # Node.js (http-server)
   npx http-server
   
   # VS Code Live Server extension
   # Right-click index.html → "Open with Live Server"
   ```
3. Open `http://localhost:8000` in your browser

### Usage
1. **Load Audio**: Click "Select Audio File" and choose an audio file
2. **Select Function**: Choose from 50+ mathematical transformations
3. **Adjust Parameters**: Fine-tune FFT size, overlap, power, and other settings
4. **Process**: Click "Play Transformed Audio" to hear the result
5. **Download**: Save the processed audio as a WAV file

## 🔧 Technical Details

### Audio Processing Pipeline
```
Audio File → PCM Data → Windowing → FFT → Complex Transform → IFFT → Overlap-Add → Output
```

### Key Technologies
- **Web Audio API** - Modern browser audio processing
- **Custom FFT** - Optimized Cooley-Tukey implementation
- **Complex Mathematics** - Using math.js library
- **Canvas API** - Real-time visualization
- **Typed Arrays** - Efficient memory management

### Performance Optimizations
- **Memory Pooling** - Reuse audio buffers to prevent garbage collection
- **Windowing Cache** - Pre-computed Hann windows for different sizes
- **Efficient FFT** - Bit-reversal and in-place operations
- **Real-time Processing** - Optimized for smooth audio playback

## 🎨 Supported Audio Formats
- **WAV** - Uncompressed audio
- **MP3** - Compressed audio
- **OGG** - Open-source format
- **M4A** - AAC compressed audio
- **FLAC** - Lossless compression

## 🔍 Mathematical Background

### Complex Numbers
All transformations work with complex numbers z = a + bi, where:
- **Real part (a)**: Magnitude information
- **Imaginary part (b)**: Phase information
- **Magnitude**: |z| = √(a² + b²)
- **Phase**: arg(z) = atan2(b, a)

### FFT Processing
The application uses the **Cooley-Tukey FFT algorithm** with:
- **Bit-reversal permutation** for efficient computation
- **In-place operations** to minimize memory usage
- **Windowing** to reduce spectral leakage
- **Overlap-add** for seamless reconstruction

### STFT Analysis
Short-Time Fourier Transform provides:
- **Time-frequency representation** of audio signals
- **Adjustable window sizes** for different resolutions
- **Overlap processing** for smooth transitions
- **Phase coherence** preservation

## 🎯 Use Cases

### Music Production
- **Harmonic enhancement** with power functions
- **Spectral filtering** using complex transformations
- **Creative effects** with fractal mathematics
- **Frequency domain processing** for unique sounds

### Audio Research
- **Mathematical modeling** of audio phenomena
- **Complex signal analysis** in frequency domain
- **Algorithm development** for audio processing
- **Educational demonstrations** of mathematical concepts

### Sound Design
- **Unique textures** through fractal transformations
- **Geometric audio effects** for spatial processing
- **Mathematical synthesis** of complex timbres
- **Experimental audio processing** techniques

## 🛠️ Development

### Project Structure
```
├── index.html          # Main HTML interface
├── app.js             # Core application logic
├── .github/
│   └── copilot-instructions.md  # AI coding guidelines
└── README.md          # Documentation
```

### Memory Management
The application implements efficient memory management:
```javascript
// Buffer pooling system
const memoryPool = {
  realBuffers: [],      // Real-valued audio buffers
  imagBuffers: [],      // Imaginary-valued buffers
  windowBuffers: new Map(),  // Cached window functions
  fftBuffers: []        // FFT working buffers
};
```

### Adding New Functions
To add a new mathematical transformation:
1. Add function definition to `complexFunctions` object
2. Implement transformation in `applyComplexFunction`
3. Test with various audio inputs
4. Document mathematical properties

## 🎪 Advanced Features

### Fractal Audio Processing
Fractal mathematics creates unique audio textures:
- **Self-similar patterns** across frequency bands
- **Infinite detail** at different scales
- **Chaotic dynamics** for unpredictable effects
- **Mathematical beauty** in audio form

### Real-time Visualization
Advanced visualization features:
- **Frequency spectrum** with color mapping
- **Smooth animations** at 60 FPS
- **Responsive design** for mobile devices
- **GPU acceleration** where available

### Professional Audio Quality
- **32-bit float** internal processing
- **16-bit PCM** export for compatibility
- **Anti-aliasing** through proper windowing
- **Phase coherence** preservation

## 🎨 Customization

### Styling
The interface uses modern CSS with:
- **Glassmorphism** design elements
- **Gradient backgrounds** for visual appeal
- **Smooth transitions** and hover effects
- **Responsive grid** layouts

### Mathematical Extensions
Easy to extend with:
- **Custom complex functions** 
- **New fractal algorithms**
- **Advanced windowing functions**
- **Alternative FFT implementations**

## 🔧 Browser Support

### Minimum Requirements
- **Chrome 66+** (full Web Audio API support)
- **Firefox 60+** (Web Audio API compatibility)
- **Safari 12+** (modern JavaScript features)
- **Edge 79+** (Chromium-based)

### Optimal Performance
- **Desktop browsers** for complex processing
- **8GB+ RAM** for large audio files
- **Modern CPU** for real-time processing
- **Dedicated graphics** for smooth visualization

## 🎵 Audio Examples

### Recommended Test Files
- **Sine waves** - Clean mathematical analysis
- **Musical instruments** - Harmonic content
- **Vocals** - Complex spectral structure
- **Percussion** - Transient analysis
- **White noise** - Frequency distribution

### Processing Tips
- **Start with power=1.0** for baseline comparison
- **Use moderate overlap** (4x) for most applications
- **Adjust FFT size** based on frequency resolution needs
- **Monitor output levels** to prevent clipping

## 🎯 Performance Tips

### For Large Files
- **Use smaller FFT sizes** (1024-2048) for faster processing
- **Reduce overlap factor** for less CPU usage
- **Close other applications** for maximum performance
- **Use Chrome** for best Web Audio API performance

### For Real-time Use
- **Enable hardware acceleration** in browser settings
- **Use dedicated audio interface** for low latency
- **Minimize other browser tabs** to free resources
- **Consider desktop app** for professional use

## 🔮 Future Enhancements

### Planned Features
- **Multi-channel processing** for stereo/surround
- **MIDI control** for real-time parameter adjustment
- **Preset system** for saving favorite configurations
- **Batch processing** for multiple files
- **Plugin system** for custom transformations

### Advanced Mathematics
- **Quantum field theory** applications
- **Differential geometry** transformations
- **Topology** preserving operations
- **Machine learning** integration

## 🤝 Contributing

### Getting Started
1. **Fork the repository**
2. **Create feature branch**
3. **Implement changes**
4. **Add tests** for new functions
5. **Submit pull request**

### Guidelines
- **Follow ES6+ standards**
- **Maintain performance** for real-time use
- **Document mathematical** properties
- **Test with various** audio formats
- **Preserve backward** compatibility

## 📚 Mathematical References

### Complex Analysis
- **Ahlfors, L.** - Complex Analysis
- **Conway, J.** - Functions of One Complex Variable
- **Rudin, W.** - Real and Complex Analysis

### Digital Signal Processing
- **Oppenheim, A.** - Discrete-Time Signal Processing
- **Proakis, J.** - Digital Signal Processing
- **Smith, S.** - The Scientist and Engineer's Guide to DSP

### Fractal Mathematics
- **Mandelbrot, B.** - The Fractal Geometry of Nature
- **Falconer, K.** - Fractal Geometry
- **Edgar, G.** - Measure, Topology, and Fractal Geometry

## 📄 License

This project is released under the **MIT License**. See LICENSE file for details.

## 🙏 Acknowledgments

- **Math.js** library for complex number operations
- **Web Audio API** specification contributors
- **Fractal mathematics** research community
- **Open source** audio processing projects

---

**Transform your audio with the power of complex mathematics!** 🎵✨
