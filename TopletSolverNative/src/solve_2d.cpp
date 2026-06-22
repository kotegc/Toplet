#include <vector>
#include <cmath>
#include <algorithm>

#include "../include/Eigen/Sparse"
#include "../include/Eigen/IterativeLinearSolvers"
#include "../include/toplet_solver.h"

// ---------------------------------------------------------------------------
// Index helpers — match the C# 2D node indexing formulas exactly
// ---------------------------------------------------------------------------
static inline int node_index_2d(int ix, int iy, int nny)
{
    return ix * nny + iy;
}

static void element_dofs_2d(int ex, int ey, int nny, int dofs[8])
{
    int n[4];
    n[0] = node_index_2d(ex,   ey,   nny);
    n[1] = node_index_2d(ex+1, ey,   nny);
    n[2] = node_index_2d(ex+1, ey+1, nny);
    n[3] = node_index_2d(ex,   ey+1, nny);
    for (int i = 0; i < 4; i++) {
        dofs[2*i]   = 2*n[i];
        dofs[2*i+1] = 2*n[i]+1;
    }
}

// ---------------------------------------------------------------------------
// 8x8 element stiffness matrix for a unit-square quad element.
// Plane stress, E=1, Poisson ratio nu. Uses 2x2 Gauss quadrature.
// ---------------------------------------------------------------------------
static void build_ke_2d(double nu, double KE[8][8])
{
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            KE[i][j] = 0.0;

    // Plane-stress elasticity matrix C (3x3), E=1
    const double c1 = 1.0 / (1.0 - nu*nu);
    const double c2 = nu / (1.0 - nu*nu);
    const double c3 = 0.5 * (1.0 - nu) / (1.0 - nu*nu);
    double C[3][3] = {};
    C[0][0]=c1; C[0][1]=c2;
    C[1][0]=c2; C[1][1]=c1;
    C[2][2]=c3;

    // Node positions of unit square — order matches element_dofs_2d above
    static const double xc[4] = {0, 1, 1, 0};
    static const double yc[4] = {0, 0, 1, 1};

    // Sign arrays for bilinear shape function derivatives
    static const double sx[4] = {-1, 1, 1, -1};
    static const double sy[4] = {-1,-1, 1,  1};

    const double g = 1.0 / std::sqrt(3.0);
    const double gauss[2] = {-g, g};

    for (int ix = 0; ix < 2; ix++)
    for (int iy = 0; iy < 2; iy++) {
        const double xi  = gauss[ix];
        const double eta = gauss[iy];

        // Shape function natural-coordinate derivatives (2 x 4)
        double dN[2][4];
        for (int i = 0; i < 4; i++) {
            dN[0][i] = 0.25 * sx[i] * (1.0 + sy[i]*eta);
            dN[1][i] = 0.25 * sy[i] * (1.0 + sx[i]*xi);
        }

        // Jacobian J[2][2]
        double J[2][2] = {};
        for (int i = 0; i < 4; i++) {
            J[0][0] += dN[0][i]*xc[i]; J[0][1] += dN[0][i]*yc[i];
            J[1][0] += dN[1][i]*xc[i]; J[1][1] += dN[1][i]*yc[i];
        }

        const double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
        const double d = 1.0 / detJ;
        double invJ[2][2];
        invJ[0][0] =  d*J[1][1]; invJ[0][1] = -d*J[0][1];
        invJ[1][0] = -d*J[1][0]; invJ[1][1] =  d*J[0][0];

        // Physical derivatives (2 x 4)
        double dNxy[2][4] = {};
        for (int r = 0; r < 2; r++)
        for (int n = 0; n < 4; n++)
        for (int k = 0; k < 2; k++)
            dNxy[r][n] += invJ[r][k] * dN[k][n];

        // Strain-displacement matrix B (3 x 8) for plane stress
        double B[3][8] = {};
        for (int n = 0; n < 4; n++) {
            int col = 2*n;
            const double Nx = dNxy[0][n], Ny = dNxy[1][n];
            B[0][col]   = Nx;
            B[1][col+1] = Ny;
            B[2][col]   = Ny; B[2][col+1] = Nx;
        }

        // KE += B^T C B detJ   (Gauss weight = 1)
        double CB[3][8] = {};
        for (int r = 0; r < 3; r++)
        for (int s = 0; s < 8; s++)
        for (int k = 0; k < 3; k++)
            CB[r][s] += C[r][k] * B[k][s];

        for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
        for (int k = 0; k < 3; k++)
            KE[i][j] += B[k][i] * CB[k][j] * detJ;
    }
}

// ---------------------------------------------------------------------------
// Main exported function
// ---------------------------------------------------------------------------
int solve_2d(
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
    progress_callback_t progress_cb)
{
    try {
        const int nny = nely + 1;
        const int nnx = nelx + 1;
        const int node_count = nnx * nny;
        const int dof_count  = node_count * 2;
        const int elem_count = nelx * nely;

        // Precompute element stiffness (same for all elements)
        double KE[8][8];
        build_ke_2d(nu, KE);

        // --- Identify which nodes are touched by at least one active element ---
        std::vector<bool> active_node(node_count, false);
        for (int ex = 0; ex < nelx; ex++)
        for (int ey = 0; ey < nely; ey++) {
            if (!mask_flat[ex*nely + ey]) continue;
            int dofs[8];
            element_dofs_2d(ex, ey, nny, dofs);
            for (int i = 0; i < 8; i++) active_node[dofs[i]/2] = true;
        }

        // --- Build free DOF list ---
        std::vector<int> free_dofs;
        free_dofs.reserve(dof_count);
        for (int node = 0; node < node_count; node++) {
            if (!active_node[node]) continue;
            for (int d = 0; d < 2; d++) {
                int dof = 2*node + d;
                if (!fixed_dofs[dof]) free_dofs.push_back(dof);
            }
        }
        const int nfree = (int)free_dofs.size();

        std::vector<int> dof_to_free(dof_count, -1);
        for (int i = 0; i < nfree; i++) dof_to_free[free_dofs[i]] = i;

        // --- Initialize element densities ---
        std::vector<double> x(elem_count, 0.0);
        for (int e = 0; e < elem_count; e++)
            if (mask_flat[e]) x[e] = vol_frac;

        // --- Precompute sensitivity filter stencil ---
        const int rfloor = (int)std::floor(filter_rad);
        struct FilterEntry { int dx, dy; double weight; };
        std::vector<FilterEntry> stencil;
        for (int ddx = -rfloor; ddx <= rfloor; ddx++)
        for (int ddy = -rfloor; ddy <= rfloor; ddy++) {
            double dist = std::sqrt((double)(ddx*ddx + ddy*ddy));
            double w = filter_rad - dist;
            if (w > 0.0) stencil.push_back({ddx, ddy, w});
        }

        std::vector<double> dc(elem_count, 0.0);
        std::vector<double> dc_filt(elem_count, 0.0);
        std::vector<double> x_new(elem_count, 0.0);
        double compliance = 0.0;

        // ===================================================================
        // Main iteration loop
        // ===================================================================
        for (int iter = 0; iter < max_iter; iter++) {

            // --- Assemble sparse Kff ---
            typedef Eigen::Triplet<double> T;
            std::vector<T> triplets;
            triplets.reserve(elem_count * 8);

            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++) {
                const int eidx = ex*nely + ey;
                if (!mask_flat[eidx]) continue;

                const double rho = x[eidx];
                const double Ee  = Emin + std::pow(rho, penal) * (E0 - Emin);

                int edofs[8];
                element_dofs_2d(ex, ey, nny, edofs);

                for (int i = 0; i < 8; i++) {
                    const int fi = dof_to_free[edofs[i]];
                    if (fi < 0) continue;
                    for (int j = 0; j < 8; j++) {
                        const int fj = dof_to_free[edofs[j]];
                        if (fj < 0) continue;
                        triplets.emplace_back(fi, fj, Ee * KE[i][j]);
                    }
                }
            }

            Eigen::SparseMatrix<double> Kff(nfree, nfree);
            Kff.setFromTriplets(triplets.begin(), triplets.end());

            Eigen::VectorXd Ff(nfree);
            for (int i = 0; i < nfree; i++) Ff(i) = forces[free_dofs[i]];

            Eigen::ConjugateGradient<
                Eigen::SparseMatrix<double>,
                Eigen::Lower | Eigen::Upper
            > cg;
            cg.setMaxIterations(std::min(5000, nfree / 2 + 1000));
            cg.setTolerance(1e-6);
            cg.compute(Kff);
            const Eigen::VectorXd Uf = cg.solve(Ff);

            Eigen::VectorXd U = Eigen::VectorXd::Zero(dof_count);
            for (int i = 0; i < nfree; i++) U(free_dofs[i]) = Uf(i);

            // --- Compliance and sensitivities ---
            compliance = 0.0;
            std::fill(dc.begin(), dc.end(), 0.0);

            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++) {
                const int eidx = ex*nely + ey;
                if (!mask_flat[eidx]) continue;

                int edofs[8];
                element_dofs_2d(ex, ey, nny, edofs);

                double ue[8];
                for (int i = 0; i < 8; i++) ue[i] = U(edofs[i]);

                double Keue[8] = {};
                for (int i = 0; i < 8; i++)
                for (int j = 0; j < 8; j++)
                    Keue[i] += KE[i][j] * ue[j];
                double ce = 0.0;
                for (int i = 0; i < 8; i++) ce += ue[i] * Keue[i];

                const double rho = x[eidx];
                const double Ee  = Emin + std::pow(rho, penal) * (E0 - Emin);
                compliance += Ee * ce;
                dc[eidx] = -penal * std::pow(rho, penal - 1.0) * (E0 - Emin) * ce;
            }

            // --- Sensitivity filter ---
            std::fill(dc_filt.begin(), dc_filt.end(), 0.0);
            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++) {
                const int eidx = ex*nely + ey;
                if (!mask_flat[eidx]) continue;

                double sum_w = 0.0, val = 0.0;
                for (const auto& s : stencil) {
                    int i2 = ex + s.dx, j2 = ey + s.dy;
                    if (i2 < 0 || i2 >= nelx || j2 < 0 || j2 >= nely) continue;
                    const int nidx = i2*nely + j2;
                    if (!mask_flat[nidx]) continue;
                    sum_w += s.weight;
                    val   += s.weight * x[nidx] * dc[nidx];
                }
                const double denom = x[eidx] * sum_w;
                dc_filt[eidx] = (denom > 1e-12) ? val / denom : dc[eidx];
            }

            // --- OC update ---
            int active_count = 0;
            for (int e = 0; e < elem_count; e++) if (mask_flat[e]) active_count++;

            const double target_vol = vol_frac * active_count;
            const double move = 0.2, xmin = 0.001;
            double l1 = 0.0, l2 = 1e9;

            while ((l2 - l1) / (l1 + l2 + 1e-12) > 1e-4) {
                const double lmid = 0.5 * (l1 + l2);
                double total = 0.0;

                for (int ex = 0; ex < nelx; ex++)
                for (int ey = 0; ey < nely; ey++) {
                    const int eidx = ex*nely + ey;
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

            // --- Convergence check ---
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
                *compliance_out = compliance;
                *iterations_out = iter + 1;
                return 0;
            }
        }

        std::copy(x.begin(), x.end(), density_out);
        *compliance_out = compliance;
        *iterations_out = max_iter;
        return 0;
    }
    catch (...) {
        return 1;
    }
}
