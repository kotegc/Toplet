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

**Status:** active WIP, and honestly a bit messy under the hood — the solver
internals have been rewritten several times chasing a setup that actually
works, and that hasn't landed yet. The Rhino/Grasshopper-side workflow works
end to end; the solvers themselves don't produce valid results yet — still
tuning methods and parameters. See **Commands** below for what the plugin
does, and **Solver backends** for what's still unresolved under the hood.

## Commands

Toplet adds two Rhino commands:

- **`Toplet2D`** (`Toplet2DCommand.cs`) — takes a closed curve and voxelizes
  the region it encloses into a 2D grid.
- **`Toplet3D`** (`Toplet3DCommand.cs`) — takes a closed Brep and voxelizes
  its interior into a 3D grid.

Both then follow the same workflow: set point/face loads and supports on the
voxelized domain, run a solve, and view the resulting material-density
result in the Rhino viewport.

## Solver backends

Each solve ultimately comes down to solving a large sparse linear system,
and there are a few different ways to do that:

- **Conjugate Gradient (CG) with an incomplete-Cholesky preconditioner** —
  Eigen's built-in implementation. Straightforward, no setup cost, but
  doesn't scale well as grid size grows.
- **Geometric multigrid (GMG)**, hand-rolled for this project, used as a CG
  preconditioner — builds an explicit coarse-grid hierarchy over the voxel
  domain to accelerate convergence at larger grid sizes.
- **AMGCL** — a vendored algebraic multigrid library
  (`TopletSolverNative/include/amgcl/`). Algebraic multigrid builds its
  hierarchy from the matrix itself rather than an explicit geometric grid,
  so it doesn't need the hand-tuned coarsening logic GMG does here.

Right now, `solve_2d.cpp` uses the first and `solve_3d.cpp` uses the second;
AMGCL is vendored but not wired into either. That pairing isn't settled —
which backend ends up serving which command is still an open question, not
a deliberate or final choice. The one *active* comparison is GMG vs. AMGCL
for the 3D backend specifically; neither has been declared a winner. The
main open problem is getting a solve to converge to a *valid* result fast
enough to be usable interactively, not just fast in isolation.

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

- Voxelized optimization domain from a closed curve (2D) or closed Brep (3D)
- Point and face-based constraints
- Point loads/supports for localized control
- Face loads/supports with distributed forces along surface normals
- Compliance-driven material distribution (SIMP-style)
- Direct Rhino interaction for geometry and constraint setup

## Roadmap

1. **Solver validity** *(current)* — get 2D and 3D both converging to a
   physically sensible result; resolve the GMG-vs-AMGCL trade study for the
   3D backend. Whether 2D should end up on the same backend as 3D is an open
   question, not yet decided.
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
