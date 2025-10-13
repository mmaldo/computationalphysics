computationalphysics
====================

This repository contains numerical solvers and computational physics code for various classical physics problems. All code has been comprehensively tested and verified for physical accuracy and numerical correctness.

## Repository Status

✅ **All physics code verified and tested** - See [COMPREHENSIVE_TEST_REPORT.md](COMPREHENSIVE_TEST_REPORT.md)

## Physics Codes

### 1. **FinalProject/** - 2D Laplace Equation with Thermal Effects
Solves driven classical diffusion with spatially correlated disorder using thermal Gauss-Seidel/SOR method.

**Physics:**
- 2D Laplace equation with Boltzmann-weighted thermal effects
- External driving field parameter (F)
- Temperature parameter (T) for thermal fluctuations
- Spatially correlated disorder via Gaussian convolution

**Status:** ✅ Fully verified - See [FinalProject/VERIFICATION_REPORT.md](FinalProject/VERIFICATION_REPORT.md)

### 2. **logmap/** - Logistic Map Bifurcation Diagram
Demonstrates route to chaos through period-doubling bifurcations in the logistic map.

**Physics:**
- Logistic map: x_{n+1} = r * x_n * (1 - x_n)
- Bifurcation diagram for r ∈ [2.5, 4.0]
- Period-doubling cascade
- Chaotic regime analysis

**Status:** ✅ Verified and tested

### 3. **newtonsMethod/** - Newton's Method for Root Finding
Fast root-finding algorithm demonstrating quadratic convergence.

**Physics/Math:**
- Solves x² - 5 = 0 (finds sqrt(5))
- Newton's method with quadratic convergence
- Achieves machine precision in ~5 iterations
- Babylonian square root method

**Status:** ✅ Verified and tested

### 4. **fixedPointSolve/** - Fixed Point Iteration Method
Fundamental iterative method demonstrating linear convergence.

**Physics/Math:**
- Solves x² - 5 = 0 using fixed point iteration
- Linear convergence (slower than Newton's)
- Demonstrates convergence criterion: 0 < A < 2/sqrt(5)
- Shows trade-off between simplicity and speed

**Status:** ✅ Verified and tested

### 5. **pendulum/** - Driven Damped Pendulum with RK4
Classical mechanics example of deterministic chaos.

**Physics:**
- Driven damped nonlinear pendulum
- RK4 integration (4th-order accuracy)
- Chaos and sensitivity to initial conditions
- Lyapunov exponent calculation framework

**Status:** ✅ Verified and tested

## Running Tests

### Run all tests:
```bash
./run_all_tests.sh
```

### Run individual tests:
```bash
# FinalProject
cd FinalProject && ./test_physics_validation.sh

# logmap
cd logmap && ./test_logmap.sh

# newtonsMethod
cd newtonsMethod && ./test_newtonsMethod.sh

# fixedPointSolve
cd fixedPointSolve && ./test_fixedPointSolve.sh

# pendulum
cd pendulum && ./test_pendulum.sh
```

## Numerical Methods Implemented

- **Gauss-Seidel / SOR:** Relaxation method for elliptic PDEs
- **Newton's Method:** Quadratic convergence root finding
- **Fixed Point Iteration:** Linear convergence iterative method
- **RK4 Integration:** 4th-order Runge-Kutta for ODEs
- **Gaussian Convolution:** Spatial correlation generation
- **Logistic Map:** Discrete dynamical system iteration

## Documentation

- [COMPREHENSIVE_TEST_REPORT.md](COMPREHENSIVE_TEST_REPORT.md) - Complete verification results for all physics codes
- [FinalProject/VERIFICATION_REPORT.md](FinalProject/VERIFICATION_REPORT.md) - Detailed thermal physics verification
- [FinalProject/VERIFICATION_SUMMARY.md](FinalProject/VERIFICATION_SUMMARY.md) - Executive summary
- [FinalProject/TESTING.md](FinalProject/TESTING.md) - Testing guide for FinalProject

## Requirements

- GCC compiler
- GNU Scientific Library (GSL) - for FinalProject only
- Standard math library (libm)

Install GSL:
```bash
# Ubuntu/Debian
sudo apt-get install libgsl-dev

# macOS
brew install gsl
```

## Compilation Examples

```bash
# FinalProject (requires GSL)
cd FinalProject
gcc -o finalproject finalproject.c -lgsl -lgslcblas -lm -Wall

# Other programs (standard libraries only)
cd logmap && gcc -o logmap logmap.c -Wall
cd newtonsMethod && gcc -o newtonsMethod newtonsMethod.c -lm -Wall
cd fixedPointSolve && gcc -o fixedPointSolve fixedPointSolve.c -lm -Wall
cd pendulum && gcc -o test test.c -lm -Wall
```

## Verification Status

All code has been comprehensively tested for:
- ✅ Physical accuracy
- ✅ Numerical correctness
- ✅ Proper implementation
- ✅ Expected behavior

**Status:** APPROVED ⭐⭐⭐⭐⭐ for educational and research use

## Author

Matthew Maldonado

## License

This code was written during undergraduate studies and is provided for educational purposes.
