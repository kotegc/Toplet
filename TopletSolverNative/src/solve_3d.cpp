#include <vector>
#include <cmath>
#include <algorithm>

#include "../include/Eigen/Sparse"
#include "../include/toplet_solver.h"

#define AMGCL_NO_BOOST
#include <amgcl/make_solver.hpp>
#include <amgcl/backend/eigen.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/cg.hpp>

// ---------------------------------------------------------------------------
// Index helpers — match the C# node indexing formulas exactly
// ---------------------------------------------------------------------------
static inline int node_index(int ix, int iy, int iz, int nny, int nnz)
{
    return ix * (nny * nnz) + iy * nnz + iz;
}

static void element_dofs(int ex, int ey, int ez, int nny, int nnz, int dofs[24])
{
    int n[8];
    n[0] = node_index(ex,   ey,   ez,   nny, nnz);
    n[1] = node_index(ex+1, ey,   ez,   nny, nnz);
    n[2] = node_index(ex+1, ey+1, ez,   nny, nnz);
    n[3] = node_index(ex,   ey+1, ez,   nny, nnz);
    n[4] = node_index(ex,   ey,   ez+1, nny, nnz);
    n[5] = node_index(ex+1, ey,   ez+1, nny, nnz);
    n[6] = node_index(ex+1, ey+1, ez+1, nny, nnz);
    n[7] = node_index(ex,   ey+1, ez+1, nny, nnz);
    for (int i = 0; i < 8; i++) {
        dofs[3*i]   = 3*n[i];
        dofs[3*i+1] = 3*n[i]+1;
        dofs[3*i+2] = 3*n[i]+2;
    }
}

// ---------------------------------------------------------------------------
// 24x24 element stiffness matrix for a unit-cube hex element.
// Uses 2x2x2 Gauss quadrature. E=1 (SIMP scaling applied during assembly).
// ---------------------------------------------------------------------------
static void build_ke(double nu, double KE[24][24])
{
    for (int i = 0; i < 24; i++)
        for (int j = 0; j < 24; j++)
            KE[i][j] = 0.0;

    // 3D isotropic elasticity matrix C (6x6), E=1
    double fac = 1.0 / ((1.0 + nu) * (1.0 - 2.0*nu));
    double a = (1.0 - nu) * fac;
    double b = nu * fac;
    double c = 0.5 * (1.0 - 2.0*nu) * fac;
    double C[6][6] = {};
    C[0][0]=a; C[0][1]=b; C[0][2]=b;
    C[1][0]=b; C[1][1]=a; C[1][2]=b;
    C[2][0]=b; C[2][1]=b; C[2][2]=a;
    C[3][3]=c; C[4][4]=c; C[5][5]=c;

    // Node positions of unit cube — order matches element_dofs above
    static const double xc[8] = {0,1,1,0, 0,1,1,0};
    static const double yc[8] = {0,0,1,1, 0,0,1,1};
    static const double zc[8] = {0,0,0,0, 1,1,1,1};

    // Sign arrays drive the trilinear shape function partial derivatives
    static const double sx[8] = {-1, 1, 1,-1, -1, 1, 1,-1};
    static const double sy[8] = {-1,-1, 1, 1, -1,-1, 1, 1};
    static const double sz[8] = {-1,-1,-1,-1,  1, 1, 1, 1};

    const double g = 1.0 / std::sqrt(3.0);
    const double gauss[2] = {-g, g};

    for (int ix = 0; ix < 2; ix++)
    for (int iy = 0; iy < 2; iy++)
    for (int iz = 0; iz < 2; iz++) {
        double xi   = gauss[ix];
        double eta  = gauss[iy];
        double zeta = gauss[iz];

        // Shape function derivatives in natural coordinates (3 x 8)
        double dN[3][8];
        for (int i = 0; i < 8; i++) {
            dN[0][i] = 0.125*sx[i]*(1.0+sy[i]*eta) *(1.0+sz[i]*zeta);
            dN[1][i] = 0.125*sy[i]*(1.0+sx[i]*xi)  *(1.0+sz[i]*zeta);
            dN[2][i] = 0.125*sz[i]*(1.0+sx[i]*xi)  *(1.0+sy[i]*eta);
        }

        // Jacobian J[3][3]
        double J[3][3] = {};
        for (int i = 0; i < 8; i++) {
            J[0][0]+=dN[0][i]*xc[i]; J[0][1]+=dN[0][i]*yc[i]; J[0][2]+=dN[0][i]*zc[i];
            J[1][0]+=dN[1][i]*xc[i]; J[1][1]+=dN[1][i]*yc[i]; J[1][2]+=dN[1][i]*zc[i];
            J[2][0]+=dN[2][i]*xc[i]; J[2][1]+=dN[2][i]*yc[i]; J[2][2]+=dN[2][i]*zc[i];
        }

        double detJ =
            J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])
           -J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])
           +J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);

        double d = 1.0 / detJ;
        double invJ[3][3];
        invJ[0][0]= d*(J[1][1]*J[2][2]-J[1][2]*J[2][1]);
        invJ[0][1]=-d*(J[0][1]*J[2][2]-J[0][2]*J[2][1]);
        invJ[0][2]= d*(J[0][1]*J[1][2]-J[0][2]*J[1][1]);
        invJ[1][0]=-d*(J[1][0]*J[2][2]-J[1][2]*J[2][0]);
        invJ[1][1]= d*(J[0][0]*J[2][2]-J[0][2]*J[2][0]);
        invJ[1][2]=-d*(J[0][0]*J[1][2]-J[0][2]*J[1][0]);
        invJ[2][0]= d*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
        invJ[2][1]=-d*(J[0][0]*J[2][1]-J[0][1]*J[2][0]);
        invJ[2][2]= d*(J[0][0]*J[1][1]-J[0][1]*J[1][0]);

        // Physical derivatives (3 x 8)
        double dNxyz[3][8] = {};
        for (int r = 0; r < 3; r++)
        for (int n = 0; n < 8; n++)
        for (int k = 0; k < 3; k++)
            dNxyz[r][n] += invJ[r][k] * dN[k][n];

        // Strain-displacement matrix B (6 x 24)
        double B[6][24] = {};
        for (int n = 0; n < 8; n++) {
            int col = 3*n;
            double Nx=dNxyz[0][n], Ny=dNxyz[1][n], Nz=dNxyz[2][n];
            B[0][col]   = Nx;
            B[1][col+1] = Ny;
            B[2][col+2] = Nz;
            B[3][col]   = Ny; B[3][col+1] = Nx;
            B[4][col+1] = Nz; B[4][col+2] = Ny;
            B[5][col]   = Nz; B[5][col+2] = Nx;
        }

        // KE += B^T C B detJ   (Gauss weight = 1 for 2-point rule)
        double CB[6][24] = {};
        for (int r = 0; r < 6; r++)
        for (int s = 0; s < 24; s++)
        for (int k = 0; k < 6; k++)
            CB[r][s] += C[r][k] * B[k][s];

        for (int i = 0; i < 24; i++)
        for (int j = 0; j < 24; j++)
        for (int k = 0; k < 6; k++)
            KE[i][j] += B[k][i] * CB[k][j] * detJ;
    }
}

// ---------------------------------------------------------------------------
// Main exported function
// ---------------------------------------------------------------------------
int solve_3d(
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
    progress_callback_t progress_cb)
{
    try {
        const int nny = nely + 1;
        const int nnz = nelz + 1;
        const int nnx = nelx + 1;
        const int node_count = nnx * nny * nnz;
        const int dof_count  = node_count * 3;
        const int elem_count = nelx * nely * nelz;
        const int nyz = nely * nelz;

        // Precompute element stiffness (same for all elements)
        double KE[24][24];
        build_ke(nu, KE);

        // --- Identify which nodes are touched by at least one active element ---
        std::vector<bool> active_node(node_count, false);
        for (int ex = 0; ex < nelx; ex++)
        for (int ey = 0; ey < nely; ey++)
        for (int ez = 0; ez < nelz; ez++) {
            if (!mask_flat[ex*nyz + ey*nelz + ez]) continue;
            int dofs[24];
            element_dofs(ex, ey, ez, nny, nnz, dofs);
            for (int i = 0; i < 24; i++) active_node[dofs[i]/3] = true;
        }

        // --- Build free DOF list (active node AND not fixed) ---
        std::vector<int> free_dofs;
        free_dofs.reserve(dof_count);
        for (int node = 0; node < node_count; node++) {
            if (!active_node[node]) continue;
            for (int d = 0; d < 3; d++) {
                int dof = 3*node + d;
                if (!fixed_dofs[dof]) free_dofs.push_back(dof);
            }
        }
        const int nfree = (int)free_dofs.size();

        // Map global DOF index → free DOF index (-1 if not free)
        std::vector<int> dof_to_free(dof_count, -1);
        for (int i = 0; i < nfree; i++) dof_to_free[free_dofs[i]] = i;

        // --- Initialize element densities ---
        std::vector<double> x(elem_count, 0.0);
        for (int e = 0; e < elem_count; e++)
            if (mask_flat[e]) x[e] = vol_frac;

        // --- Precompute sensitivity filter stencil ---
        const int rfloor = (int)std::floor(filter_rad);
        struct FilterEntry { int dx, dy, dz; double weight; };
        std::vector<FilterEntry> stencil;
        for (int ddx = -rfloor; ddx <= rfloor; ddx++)
        for (int ddy = -rfloor; ddy <= rfloor; ddy++)
        for (int ddz = -rfloor; ddz <= rfloor; ddz++) {
            double dist = std::sqrt((double)(ddx*ddx + ddy*ddy + ddz*ddz));
            double w = filter_rad - dist;
            if (w > 0.0) stencil.push_back({ddx, ddy, ddz, w});
        }

        std::vector<double> dc(elem_count, 0.0);
        std::vector<double> dc_filt(elem_count, 0.0);
        std::vector<double> x_new(elem_count, 0.0);
        double compliance = 0.0;

        // Count active elements once (mask never changes)
        int active_count = 0;
        for (int e = 0; e < elem_count; e++) if (mask_flat[e]) active_count++;

        // ===================================================================
        // Pre-build the sparsity pattern of Kff ONCE.
        // The active set is fixed, so (row,col) pairs never change — only
        // the values change as densities evolve.  We cache the position of
        // each element (i,j) contribution in the compressed value array so
        // each subsequent iteration can fill K with a single O(N) pass and
        // zero memory allocation.
        // ===================================================================
        Eigen::SparseMatrix<double, Eigen::RowMajor> Kff(nfree, nfree);
        {
            typedef Eigen::Triplet<double> T;
            std::vector<T> pat;
            pat.reserve((size_t)active_count * 576);
            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++)
            for (int ez = 0; ez < nelz; ez++) {
                const int eidx = ex*nyz + ey*nelz + ez;
                if (!mask_flat[eidx]) continue;
                int edofs[24];
                element_dofs(ex, ey, ez, nny, nnz, edofs);
                for (int i = 0; i < 24; i++) {
                    const int fi = dof_to_free[edofs[i]]; if (fi < 0) continue;
                    for (int j = 0; j < 24; j++) {
                        const int fj = dof_to_free[edofs[j]]; if (fj < 0) continue;
                        pat.emplace_back(fi, fj, 0.0);
                    }
                }
            }
            Kff.setFromTriplets(pat.begin(), pat.end());
        }
        Kff.makeCompressed();

        // Cache: for each active element, the list of (value-array-index, KE entry).
        // Binary search on the compressed column to map (fi,fj) → valuePtr position.
        struct KEntry { int vidx; double ke; };
        std::vector<std::vector<KEntry>> elem_k(elem_count);
        {
            const int* outer = Kff.outerIndexPtr();
            const int* inner = Kff.innerIndexPtr();
            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++)
            for (int ez = 0; ez < nelz; ez++) {
                const int eidx = ex*nyz + ey*nelz + ez;
                if (!mask_flat[eidx]) continue;
                int edofs[24];
                element_dofs(ex, ey, ez, nny, nnz, edofs);
                auto& ev = elem_k[eidx];
                ev.reserve(576);
                for (int i = 0; i < 24; i++) {
                    const int fi = dof_to_free[edofs[i]]; if (fi < 0) continue;
                    for (int j = 0; j < 24; j++) {
                        const int fj = dof_to_free[edofs[j]]; if (fj < 0) continue;
                        int lo = outer[fi], hi = outer[fi + 1];
                        int pos = (int)(std::lower_bound(inner + lo, inner + hi, fj) - inner);
                        ev.push_back({pos, KE[i][j]});
                    }
                }
            }
        }

        // RHS is constant (forces don't change between iterations)
        Eigen::VectorXd Ff(nfree);
        for (int i = 0; i < nfree; i++) Ff(i) = forces[free_dofs[i]];

        using AmgBackend = amgcl::backend::eigen<double>;
        using AmgSolver = amgcl::make_solver<
            amgcl::amg<
                AmgBackend,
                amgcl::coarsening::smoothed_aggregation,
                amgcl::relaxation::spai0
            >,
            amgcl::solver::cg<AmgBackend>
        >;
        AmgSolver::params amg_prm;
        amg_prm.solver.tol     = 1e-3;
        amg_prm.solver.maxiter = 100;

        Eigen::VectorXd Uf = Eigen::VectorXd::Zero(nfree);

        // ===================================================================
        // Main iteration loop
        // ===================================================================
        for (int iter = 0; iter < max_iter; iter++) {

            // --- Fill Kff values directly (zero then accumulate) ---
            double* kvals = Kff.valuePtr();
            std::fill(kvals, kvals + Kff.nonZeros(), 0.0);

            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++)
            for (int ez = 0; ez < nelz; ez++) {
                const int eidx = ex*nyz + ey*nelz + ez;
                if (!mask_flat[eidx]) continue;
                const double Ee = Emin + std::pow(x[eidx], penal) * (E0 - Emin);
                for (const auto& e : elem_k[eidx])
                    kvals[e.vidx] += Ee * e.ke;
            }

            // --- Solve Kff * Uf = Ff (warm-started from previous step) ---
            AmgSolver amg(Kff, amg_prm);
            size_t amg_iters; double amg_error;
            std::tie(amg_iters, amg_error) = amg(Ff, Uf);

            // Expand to full displacement vector
            Eigen::VectorXd U = Eigen::VectorXd::Zero(dof_count);
            for (int i = 0; i < nfree; i++) U(free_dofs[i]) = Uf(i);

            // --- Compute compliance and sensitivities ---
            compliance = 0.0;
            std::fill(dc.begin(), dc.end(), 0.0);

            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++)
            for (int ez = 0; ez < nelz; ez++) {
                const int eidx = ex*nyz + ey*nelz + ez;
                if (!mask_flat[eidx]) continue;

                int edofs[24];
                element_dofs(ex, ey, ez, nny, nnz, edofs);

                double ue[24];
                for (int i = 0; i < 24; i++) ue[i] = U(edofs[i]);

                // ce = ue^T * KE * ue
                double Keue[24] = {};
                for (int i = 0; i < 24; i++)
                for (int j = 0; j < 24; j++)
                    Keue[i] += KE[i][j] * ue[j];
                double ce = 0.0;
                for (int i = 0; i < 24; i++) ce += ue[i] * Keue[i];

                const double rho = x[eidx];
                const double Ee  = Emin + std::pow(rho, penal) * (E0 - Emin);
                compliance += Ee * ce;
                dc[eidx] = -penal * std::pow(rho, penal - 1.0) * (E0 - Emin) * ce;
            }

            // --- Sensitivity filter ---
            std::fill(dc_filt.begin(), dc_filt.end(), 0.0);
            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++)
            for (int ez = 0; ez < nelz; ez++) {
                const int eidx = ex*nyz + ey*nelz + ez;
                if (!mask_flat[eidx]) continue;

                double sum_w = 0.0, val = 0.0;
                for (const auto& s : stencil) {
                    int i2 = ex + s.dx, j2 = ey + s.dy, k2 = ez + s.dz;
                    if (i2 < 0 || i2 >= nelx || j2 < 0 || j2 >= nely || k2 < 0 || k2 >= nelz) continue;
                    const int nidx = i2*nyz + j2*nelz + k2;
                    if (!mask_flat[nidx]) continue;
                    sum_w += s.weight;
                    val   += s.weight * x[nidx] * dc[nidx];
                }
                const double denom = x[eidx] * sum_w;
                dc_filt[eidx] = (denom > 1e-12) ? val / denom : dc[eidx];
            }

            // --- OC density update with bisection for volume constraint ---
            const double target_vol = vol_frac * active_count;
            const double move = 0.2, xmin = 0.001;
            double l1 = 0.0, l2 = 1e9;

            while ((l2 - l1) / (l1 + l2 + 1e-12) > 1e-4) {
                const double lmid = 0.5 * (l1 + l2);
                double total = 0.0;

                for (int ex = 0; ex < nelx; ex++)
                for (int ey = 0; ey < nely; ey++)
                for (int ez = 0; ez < nelz; ez++) {
                    const int eidx = ex*nyz + ey*nelz + ez;
                    if (!mask_flat[eidx]) { x_new[eidx] = 0.0; continue; }

                    const double rho = x[eidx];
                    double b = -dc_filt[eidx] / lmid;
                    if (b < 1e-12) b = 1e-12;
                    double cand = rho * std::sqrt(b);
                    const double lo = std::max(xmin, rho - move);
                    const double hi = std::min(1.0,  rho + move);
                    cand = std::max(lo, std::min(hi, cand));
                    x_new[eidx] = cand;
                    total += cand;
                }

                if (total > target_vol) l1 = lmid; else l2 = lmid;
            }

            // --- Convergence check, then update densities ---
            double change = 0.0;
            for (int e = 0; e < elem_count; e++) {
                if (!mask_flat[e]) { x[e] = 0.0; continue; }
                const double diff = std::abs(x_new[e] - x[e]);
                if (diff > change) change = diff;
                x[e] = x_new[e];
            }

            if (progress_cb) progress_cb(iter + 1, max_iter, compliance);

            if (change < 0.01 && iter > 10) {
                std::copy(x.begin(), x.end(), density_out);
                *compliance_out  = compliance;
                *iterations_out  = iter + 1;
                return 0;
            }
        }

        std::copy(x.begin(), x.end(), density_out);
        *compliance_out  = compliance;
        *iterations_out  = max_iter;
        return 0;
    }
    catch (...) {
        return 1;
    }
}
