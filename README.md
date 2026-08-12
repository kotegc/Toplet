# Toplet

Voxel-based topology optimization toy solver for Rhino.

Toplet is a toy solver intended for experimentation and learning. It is not
numerically robust or validated for engineering use. It implements a
simplified, voxel-based solver inspired by the classic Top99 algorithm (Ole
Sigmund), adapted for interactive use directly inside the Rhino viewport.

<img width="300" height="200" alt="Screenshot 2026-04-04 031519" src="https://github.com/user-attachments/assets/da3d6674-33fc-4081-b6e3-f1ee97e7c883" />
<img width="300" height="200" alt="Screenshot 2026-04-04 031536" src="https://github.com/user-attachments/assets/024178cf-ff3c-476d-9c82-07973b53f8b9" />

The goal is not production-grade analysis, but a clear and flexible
environment for exploring topology optimization concepts within a CAD
workflow.

**Status:** active WIP. The Rhino/Grasshopper-side workflow — voxelizing a
Brep into a design domain, setting point/face loads and supports, running a
solve, viewing the result — works end to end. The solver itself doesn't
produce valid results yet; I'm still tuning methods and parameters. Part of
that: 2D and 3D currently use two different linear-solve strategies as an
active trade study —

- **2D** (`TopletSolverNative/src/solve_2d.cpp`) uses Eigen's Conjugate
  Gradient with an incomplete-Cholesky preconditioner.
- **3D** (`TopletSolverNative/src/solve_3d.cpp`) uses a hand-rolled
  geometric-multigrid-preconditioned CG, since plain CG didn't scale to 3D
  grid sizes.
- **AMGCL** (vendored under `TopletSolverNative/include/amgcl/`) is a third
  candidate — algebraic multigrid, no hand-tuned grid hierarchy needed — but
  isn't wired into either solve path yet.

None of the three has been declared a winner. The main open problem is
getting a solve to converge to a *valid* result fast enough to be usable
interactively, not just fast in isolation.

## How it works

- **[`TopOpt2D/`](TopOpt2D)** / **[`TopOpt3D/`](TopOpt3D)** — Rhino command
  implementations (`Toplet2DCommand.cs`, `Toplet3DCommand.cs`) and their
  supporting types: domain/problem/result POCOs, voxelization
  (`VoxelDomainBuilder3D`), support-connectivity filtering, viewport display.
- **[`Interop/`](Interop)** — the P/Invoke boundary to the native solver:
  marshals mask/force/constraint arrays across, runs the solve on a
  background task while a progress dialog owns the UI thread, and unpacks
  the result back into managed arrays.
- **[`TopletSolverNative/`](TopletSolverNative)** — the native C++ solver
  core (SIMP + optimality-criteria update + sensitivity filter), built as a
  separate DLL and loaded by the C# plugin at runtime. Eigen and AMGCL are
  vendored in `include/` rather than pulled via a package manager, since
  there's no native C++ package manager as standard in the Rhino plugin
  ecosystem as NuGet is for .NET.

## Building

- **C# plugin**: .NET Framework 4.8, [RhinoCommon](https://www.nuget.org/packages/RhinoCommon)
  7.0.20314.3001 (pinned in `Toplet.csproj`).
- **Native solver**: Visual Studio 2022 (`v143` toolset), C++17, **x64
  only**. Build this project *first* — the C# plugin P/Invokes the resulting
  `TopletSolverNative.dll` at runtime and won't find it otherwise.
- From the command line, build both with an **explicit platform**:
  ```
  MSBuild Toplet.sln /t:Restore,Build /p:Configuration=Debug /p:Platform=x64
  ```
  `Toplet.sln` only defines a `Build.0` mapping for `TopletSolverNative`
  under `x64` — omitting `/p:Platform=x64` silently builds just the C#
  project and skips the native one with no error.
- `Properties/launchSettings.json` assumes Rhino is installed at its default
  location (`C:\Program Files\Rhino 7` or `Rhino 8`) — edit the profile's
  `executablePath` if yours lives elsewhere.

## Features

- Voxelized optimization domain from closed Breps
- Point and face-based constraints
- Point loads/supports for localized control
- Face loads/supports with distributed forces along surface normals
- Compliance-driven material distribution (SIMP-style)
- Direct Rhino interaction for geometry and constraint setup

## Roadmap

1. **Solver validity** *(current)* — get 2D and 3D both converging to a
   physically sensible result; resolve the CG/GMG/AMGCL trade study.
2. **Performance** — once results are valid, get the solve loop fast enough
   for interactive use at useful grid resolutions.
3. **Everything else** — UI polish, more constraint types, is on hold until
   the solver itself is trustworthy.

## License

MIT — see [LICENSE](LICENSE). The vendored [Eigen](https://eigen.tuxfamily.org)
(MPL2) and [amgcl](https://github.com/ddemidov/amgcl) (MIT) libraries under
`TopletSolverNative/include/` retain their own upstream licenses.

## Attribution

Based on the Top99 topology optimization method:
Ole Sigmund, *A 99 line topology optimization code written in MATLAB*.
