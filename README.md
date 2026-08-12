<h1>Toplet</h1>
Voxel-based topology optimization toy solver for Rhino

Toplet is a toy solver intended for experimentation and learning. It is not numerically robust or validated for engineering use.
It implements a simplified, voxel-based solver inspired by the classic Top99 algorithm (Ole Sigmund), adapted for interactive use directly inside the Rhino viewport.

<img width="300" height="200" alt="Screenshot 2026-04-04 031519" src="https://github.com/user-attachments/assets/da3d6674-33fc-4081-b6e3-f1ee97e7c883" />
<img width="300" height="200" alt="Screenshot 2026-04-04 031536" src="https://github.com/user-attachments/assets/024178cf-ff3c-476d-9c82-07973b53f8b9" />

The goal is not production-grade analysis, but a clear and flexible environment for exploring topology optimization concepts within a CAD workflow.

<h3>Features</h3>

Voxelized optimization domain from closed Breps
Point and face-based constraints
Point loads/supports for localized control
Face loads/supports with distributed forces along surface normals
Compliance-driven material distribution (SIMP-style)
Direct Rhino interaction for geometry and constraint setup

<h3>Attribution</h3>

Based on the Top99 topology optimization method:
Ole Sigmund, A 99 line topology optimization code written in MATLAB