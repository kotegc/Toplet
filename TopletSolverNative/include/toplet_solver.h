#pragma once

#ifdef _WIN32
  #define TOPLET_API __declspec(dllexport)
#else
  #define TOPLET_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

// -----------------------------------------------------------------------
// 3D hexahedral topology optimization (SIMP + OC + sensitivity filter)
//
// Inputs
//   mask_flat   : nelx*nely*nelz bytes, index = ex*(nely*nelz)+ey*nelz+ez
//                 1 = active element (inside design domain), 0 = void
//   forces      : dofCount doubles  (dofCount = (nelx+1)*(nely+1)*(nelz+1)*3)
//   fixed_dofs  : dofCount bytes    (1 = constrained DOF, 0 = free)
//   nelx/y/z    : grid element counts
//   vol_frac    : target volume fraction (e.g. 0.5)
//   penal       : SIMP penalty exponent (e.g. 3.0)
//   filter_rad  : sensitivity filter radius in elements (e.g. 1.5)
//   max_iter    : maximum OC iterations
//   E0          : Young's modulus of solid material (typically 1.0)
//   Emin        : Young's modulus of void material (typically 1e-9)
//   nu          : Poisson's ratio (typically 0.3)
//
// Outputs (caller must pre-allocate to correct size)
//   density_out : nelx*nely*nelz doubles, same index layout as mask_flat
//   compliance  : pointer to one double (final compliance value)
//   iterations  : pointer to one int   (number of iterations performed)
//
// progress_cb : optional callback fired at the end of each iteration.
//               Signature: void cb(int iter, int max_iter, double compliance)
//               Pass NULL to disable.
// Returns: 0 = success, 1 = error
// -----------------------------------------------------------------------
typedef void (*progress_callback_t)(int iter, int max_iter, double compliance);

TOPLET_API int solve_3d(
    const unsigned char* mask_flat,
    const double*        forces,
    const unsigned char* fixed_dofs,
    int   nelx, int nely, int nelz,
    double vol_frac,
    double penal,
    double filter_rad,
    int    max_iter,
    double E0,
    double Emin,
    double nu,
    double* density_out,
    double* compliance,
    int*    iterations,
    progress_callback_t progress_cb
);

// -----------------------------------------------------------------------
// 2D quad topology optimization (plane stress, SIMP + OC + sensitivity filter)
//
//   mask_flat   : nelx*nely bytes, index = ex*nely+ey
//   forces      : (nelx+1)*(nely+1)*2 doubles
//   fixed_dofs  : same length as forces, bytes
//   density_out : nelx*nely doubles
// -----------------------------------------------------------------------
TOPLET_API int solve_2d(
    const unsigned char* mask_flat,
    const double*        forces,
    const unsigned char* fixed_dofs,
    int   nelx, int nely,
    double vol_frac,
    double penal,
    double filter_rad,
    int    max_iter,
    double E0,
    double Emin,
    double nu,
    double* density_out,
    double* compliance,
    int*    iterations,
    progress_callback_t progress_cb
);

#ifdef __cplusplus
}
#endif
