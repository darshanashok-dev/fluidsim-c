# FluidSim 

A high-performance, real-time 2D fluid simulation built with C, [raylib](https://www.raylib.com/), and [OpenMP](https://www.openmp.org/).

Based on Jos Stam's "Stable Fluids" method, this project focuses on providing a smooth, interactive experience with vibrant visuals and multi-threaded performance.

<img width="804" height="799" alt="image" src="https://github.com/user-attachments/assets/12c8da88-56c2-46db-acd1-9d4ccde9d845" />

<img width="804" height="799" alt="1778916684266253451" src="https://github.com/user-attachments/assets/4545f2f1-1db6-4576-8475-666d63f52c04" />



##  Features

- **Real-time Interaction**: Add fluid density and velocity dynamically using your mouse.
- **Multi-threaded Solver**: Powered by OpenMP for high-resolution simulations on multi-core CPUs.
- **Dynamic Resolution**: Hot-swap between 64x64, 128x128, and 256x256 grids on the fly.
- **Physics-based Simulation**: Includes diffusion, viscosity, and obstacle interaction.
- **Vibrant Aesthetics**: Smooth color gradients and dissipation effects for a "glowy" look.

##  Getting Started

### Prerequisites

You'll need the following installed on your system:

- **C Compiler** (GCC or Clang recommended)
- **CMake** (3.14 or higher)
- **OpenMP** (Available in most modern compilers)

### Building from Source

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/FluidSim.git
   cd FluidSim
   ```

2. Create a build directory and run CMake:
   ```bash
   mkdir build
   cd build
   cmake ..
   ```

3. Compile the project:
   ```bash
   make
   ```

4. Run the executable:
   ```bash
   ./FluidSim
   ```

##  Controls

| Key / Mouse | Action |
|-------------|--------|
| **Left Mouse** | Add fluid density and velocity |
| **Right Mouse** | Add temporary obstacles |
| **Key 1** | Switch to 64x64 resolution |
| **Key 2** | Switch to 128x128 resolution |
| **Key 3** | Switch to 256x256 resolution |
| **Key R** | Reset the simulation |
| **Esc** | Exit |

##  Performance Tips

- The simulation uses `-O3 -ffast-math -march=native` flags for maximum performance.
- Ensure your environment allows OpenMP to use all available cores for the best experience at 256x256 resolution.

##  License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

##  Acknowledgments

- Jos Stam for the groundbreaking "Stable Fluids" algorithm.
- [raylib](https://github.com/raysan5/raylib) for the excellent graphics library.
