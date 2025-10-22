# AxiSymFEA

**Axisymmetric Finite Element Analysis Library for Swift**

A high-performance finite element analysis library specialized for axisymmetric structures using the Transfer Matrix Method.

## Overview

AxiSymFEA provides a complete FEA framework for analyzing axisymmetric structures (disks, shells, and shafts) commonly found in rotating machinery, pressure vessels, and mechanical components.

### Key Features

- **Transfer Matrix Method (TMM)**: Efficient semi-analytical approach for axisymmetric structures
- **Direct Stiffness Method**: Classical FEA for beam and shaft elements
- **Fourier Mode Decomposition**: Handles non-axisymmetric loading through modal analysis
- **LAPACK Integration**: High-performance linear algebra via Apple's Accelerate framework
- **Multiple Element Types**: Disk, Shell, and Shaft elements with various formulations

## Element Types

### DiskElement
Radial disk bending analysis using Transfer Matrix Bending (TMB) equations:
- Linear and power-law thickness variation
- Analytical integration for uniform thickness
- Numerical integration for variable thickness
- Gravity loading support

### ShellElement
Cylindrical shell analysis with two theoretical models:
- **Ventsel-Krauthammer**: General shell theory
- **TWK (Timoshenko-Woinowsky-Krieger)**: Sparse formulation optimized for cylindrical shells

### ShaftElement
Combined beam and bar element for shaft analysis:
- **Euler-Bernoulli**: Classical beam theory (neglects shear)
- **Timoshenko**: Includes transverse shear deformation
- Handles bending, axial, and torsional loads
- Direct 8×8 stiffness matrix formulation

## Theoretical Background

### Transfer Matrix Method

The Transfer Matrix Method relates the state vectors at two points along the structure:

```
q₁ = T · q₀ + L
```

where:
- `q₀`, `q₁`: State vectors (displacements, rotations, forces, moments)
- `T`: Transfer matrix
- `L`: Load vector

For axisymmetric structures with mode `n`, the transfer matrix is computed by integrating:

```
dq/dξ = Hₘ(ξ) · q
```

using matrix exponential: `T = exp(∫ Hₘ dξ)`

### Fourier Mode Decomposition

Non-axisymmetric loads are decomposed into Fourier modes:

```
F(θ) = F₀ + Σ [Fₙcos(nθ) + Gₙsin(nθ)]
```

Each mode is solved independently, then superposed.

## Usage

### Basic Example: Cantilever Beam

```swift
import AxiSymFEA

// Create assembly
let assembly = FEAssembly()

// Create shaft element
let shaft = ShaftElement(
    diameter: 50.0,           // mm
    axialPositionStart: 0.0,
    axialPositionEnd: 1000.0,
    youngsModulus: 210000.0,  // MPa
    poissonsRatio: 0.3,
    mode: 0,
    model: .timoshenko
)

assembly.addElement(shaft)
assembly.assemble(mode: 0)

// Apply boundary conditions (fixed at left)
let node0DOFs = assembly.getDOFIndices(forNode: shaft.nodes[0])
assembly.fixDOF(node0DOFs[0])  // w = 0
assembly.fixDOF(node0DOFs[1])  // gamma = 0
assembly.fixDOF(node0DOFs[2])  // u = 0
assembly.fixDOF(node0DOFs[3])  // beta = 0

// Apply load at free end
let node1DOFs = assembly.getDOFIndices(forNode: shaft.nodes[1])
assembly.applyForce(-1000.0, atDOF: node1DOFs[0])

// Solve
if assembly.solve() {
    if let deflection = assembly.getDisplacement(atDOF: node1DOFs[0]) {
        print("End deflection: \(deflection) mm")
    }
}
```

### Disk Analysis

```swift
let disk = DiskElement(
    innerRadius: 100.0,
    outerRadius: 500.0,
    thicknessBegin: 20.0,
    thicknessEnd: 10.0,
    youngsModulus: 210000.0,
    poissonsRatio: 0.3,
    useNumericalIntegration: true,
    mode: 0
)

// Add gravity loading
disk.addGravity(9.81, density: 7850.0)  // m/s², kg/m³
```

## Architecture

### Core Components

- **Element Protocol**: Base interface for all element types
- **Node**: Manages degrees of freedom and connectivity
- **FEAssembly**: Global assembly and solver
- **Matrix/Vector**: Linear algebra types using Accelerate
- **SparseMatrix**: Efficient sparse matrix storage

### Solver

The linear system `Kx = F` is solved using LAPACK's `dgesv`:
- LU factorization with partial pivoting
- Handles general dense systems
- Column-major format (LAPACK standard)

## Performance

- **Accelerate Framework**: BLAS operations for matrix multiplication
- **Sparse Assembly**: Only non-zero entries stored during assembly
- **Matrix Exponential**: Binary scaling (2¹² subdivisions) for stability
- **Logarithmic Spacing**: Optimal integration point distribution

## Installation

### Swift Package Manager

Add to your `Package.swift`:

```swift
dependencies: [
    .package(url: "https://github.com/andrewhustrulid/AxiSymFEA.git", from: "1.0.0")
]
```

## Requirements

- Swift 5.9+
- macOS 13+ or iOS 16+
- Accelerate framework (included with Apple platforms)

## Testing

```bash
swift test
```

The test suite includes:
- Matrix and vector operations (10 tests)
- Node index management (3 tests)
- Element initialization and computation (11 tests)
- Transfer matrix validation (8 tests)
- Assembly and solver (5 tests)

## Theory References

1. **Transfer Matrix Method**:
   - Pestel & Leckie, "Matrix Methods in Elastomechanics"

2. **Shell Theory**:
   - Ventsel & Krauthammer, "Thin Plates and Shells"
   - Timoshenko & Woinowsky-Krieger, "Theory of Plates and Shells"

3. **Beam Theory**:
   - Timoshenko, "Strength of Materials"

## Applications

AxiSymFEA is ideal for:
- Rotating machinery (pulleys, flywheels, turbine disks)
- Pressure vessels and piping
- Shaft stress analysis
- Conveyor system design
- Axisymmetric structural components

## License

MIT License - see LICENSE file

## Contributing

Contributions welcome! Please open an issue or pull request.

## Authors

Originally developed as part of the PulleySwift project, extracted as a standalone general-purpose FEA library.
