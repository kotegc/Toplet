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
//   density_out    : nelx*nely*nelz doubles, same index layout as mask_flat
//   compliance_out : pointer to one double (final compliance value)
//   iterations_out : pointer to one int   (number of iterations performed).
//                     Note: hitting max_iter without satisfying the
//                     convergence threshold still returns 0 (success), with
//                     iterations_out == max_iter — this field alone doesn't
//                     distinguish "converged on the last iteration" from
//                     "never converged," so treat iterations_out == max_iter
//                     as a signal to double-check the result.
//
// progress_cb : optional callback fired at the end of each iteration.
//               Signature: void cb(int iter, int max_iter, double compliance)
//               Pass NULL to disable.
//
// diag_buf_out      : optional caller-provided buffer for solver diagnostics
//                      (GMG hierarchy sizes, V-cycle quality check, PCG
//                      convergence per SIMP iteration). Null-terminated,
//                      truncated to fit. Pass NULL + 0 to disable — the
//                      caller decides how/whether to surface this text (e.g.
//                      printed to Rhino's command line); this function never
//                      shows UI of its own.
// diag_buf_capacity : size of diag_buf_out in bytes (recommend >= 8192).
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
    double* compliance_out,
    int*    iterations_out,
    progress_callback_t progress_cb,
    char*   diag_buf_out,
    int     diag_buf_capacity
);

// -----------------------------------------------------------------------
// 2D quad topology optimization (plane stress, SIMP + OC + sensitivity filter)
//
//   mask_flat   : nelx*nely bytes, index = ex*nely+ey
//   forces      : (nelx+1)*(nely+1)*2 doubles
//   fixed_dofs  : same length as forces, bytes
//   density_out : nelx*nely doubles
//
// Same iterations_out / progress_cb semantics as solve_3d above. No
// diagnostics buffer here — this path has no equivalent to GMG to report on.
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
    double* compliance_out,
    int*    iterations_out,
    progress_callback_t progress_cb
);

#ifdef __cplusplus
}
#endif
