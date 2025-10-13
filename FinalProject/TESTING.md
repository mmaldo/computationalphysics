# Testing Guide for finalproject.c

This document describes how to test the thermal physics implementation in `finalproject.c`.

## Quick Test

To quickly verify that the program works:

```bash
# Compile
gcc -o finalproject finalproject.c -lgsl -lgslcblas -lm -Wall

# Run with default parameters
./finalproject
```

Expected output:
- Convergence message: `# Converged after N iterations (error < eps)`
- 3-column data format (x, y, potential)
- Charge output: `# charges: value1 value2`

## Comprehensive Test Suite

Run the comprehensive test script:

```bash
./test_physics_validation.sh
```

This tests:
1. Default parameters (T=2.0, F=2.0)
2. Different grid spacings (hx=0.05, 0.2)
3. SOR parameter variations (alpha=1.5)
4. Combined parameter changes

## Manual Testing

### Test 1: Verify Compilation
```bash
gcc -o finalproject finalproject.c -lgsl -lgslcblas -lm -Wall
```
Expected: No errors or warnings

### Test 2: Basic Execution
```bash
./finalproject > output.txt 2>&1
```
Expected: 
- Program exits successfully
- Contains convergence message
- Contains charge output

### Test 3: Parameter Variations

Test different grid spacing:
```bash
./finalproject hx=0.2
```

Test different SOR parameter:
```bash
./finalproject alpha=1.5
```

Test combined parameters:
```bash
./finalproject hx=0.15 alpha=1.2
```

### Test 4: Verify Output Format

Check that output is suitable for gnuplot:
```bash
./finalproject > data.txt 2> error.txt
head -20 data.txt
```

Expected format:
```
x    y    potential
0    0    247.945
0    0.1  247.397
...
```

### Test 5: Check Convergence

Extract convergence information:
```bash
./finalproject 2>&1 > /dev/null | grep Converged
```

Expected output:
```
# Converged after 11 iterations (error=1.08943e-11 < eps=1e-10)
```

### Test 6: Verify Charge Calculation

Extract charge values:
```bash
./finalproject 2>&1 | grep "# charges:"
```

Expected output:
```
# charges:  -2258.36 -2954.68
```

## Physics Validation

### Thermal Physics Formula

The code implements:
```
A(i,j) = 0.25 * Σ A(neighbor) * exp(-(V(center)+V(neighbor))/(2*T))
```

where:
- `V(i,j) = A(i,j) + F*i` (total potential including external field)
- `T = 2.0` (temperature parameter)
- `F = 2.0` (external field strength)

**Verification:** Check lines 382-392 in `finalproject.c`

### Correlated Disorder

The code uses Gaussian convolution to create spatially correlated initial conditions:

1. **Random field generation:** Lines 233-237
   - Uniform random values in [0,1)
   
2. **Gaussian kernel:** Lines 240-249
   - Correlation length: L = hx*10
   - Form: exp(-(x²+y²)/(2*L²))

3. **Convolution:** Lines 278-294
   - Creates A = B ⊗ C
   - Smooth, correlated disorder

**Verification:** Run program and observe that initial conditions are smoothly varying

### Boundary Conditions

- **Outer electrode (domain boundary):** V = 1.0
- **Inner circular electrode:** V = 2.0 (radius R=1.0, center at (1.5, 1.5))
- **Interior points:** Updated via thermal Gauss-Seidel

**Verification:** Check output - boundary values should remain fixed

## Expected Results

### Default Parameters
- Convergence: ~11 iterations
- Final error: ~1e-11 (below eps=1e-10)
- Charges: approximately -2258 and -2955 (varies due to random initial conditions)
- Execution time: ~0.01 seconds

### Grid Spacing Effects
- Finer grid (smaller hx): More iterations, more accurate
- Coarser grid (larger hx): Fewer iterations, less accurate

### SOR Parameter Effects
- alpha = 1.0: Gauss-Seidel (standard)
- alpha > 1.0: Over-relaxation (may converge faster)
- alpha must be < 2.0 for stability

## Troubleshooting

### Compilation Errors

**Error:** `gsl/gsl_rng.h: No such file or directory`
- **Solution:** Install GSL library: `sudo apt-get install libgsl-dev`

**Error:** `undefined reference to 'gsl_rng_alloc'`
- **Solution:** Link with GSL: `-lgsl -lgslcblas -lm`

### Runtime Errors

**Error:** Program doesn't converge (max iterations reached)
- **Cause:** Grid too fine or tolerance too strict
- **Solution:** Increase Nt parameter or relax eps tolerance

**Error:** Segmentation fault
- **Cause:** Grid size exceeds array allocation
- **Solution:** Use hx ≥ 0.1 (grid size limited by GRID_SIZE_A=61)

### Output Issues

**Issue:** No charge output
- **Check:** Ensure final output includes `# charges:` line
- **Verify:** Lines 436-437 in main() function

**Issue:** Unexpected convergence behavior
- **Check:** Parameter values (T, F, hx, alpha)
- **Verify:** Physical reasonableness of setup

## Performance Notes

- Default configuration converges quickly (~11 iterations, 0.01 seconds)
- Memory usage scales with grid size: O((Nx+1)*(Ny+1))
- Computation scales as: O(Nt * Nx * Ny)
- For larger domains, increase Nt to allow convergence

## Regression Testing

To ensure no functionality is broken:

1. Run test suite: `./test_physics_validation.sh`
2. Verify all tests pass (except fine grid which may exceed memory)
3. Check VERIFICATION_REPORT.md for detailed results
4. Compare with reference implementation (finalproj.c) if needed

## See Also

- `VERIFICATION_REPORT.md` - Comprehensive verification results
- `README.md` - Project description
- `finalproj.c` - Reference implementation for comparison
