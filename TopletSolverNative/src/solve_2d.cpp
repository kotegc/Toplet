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
// Index helpers
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

    const double c1 = 1.0 / (1.0 - nu*nu);
    const double c2 = nu / (1.0 - nu*nu);
    const double c3 = 0.5 * (1.0 - nu) / (1.0 - nu*nu);
    double C[3][3] = {};
    C[0][0]=c1; C[0][1]=c2;
    C[1][0]=c2; C[1][1]=c1;
    C[2][2]=c3;

    static const double xc[4] = {0, 1, 1, 0};
    static const double yc[4] = {0, 0, 1, 1};
    static const double sx[4] = {-1, 1, 1, -1};
    static const double sy[4] = {-1,-1, 1,  1};

    const double g = 1.0 / std::sqrt(3.0);
    const double gauss[2] = {-g, g};

    for (int ix = 0; ix < 2; ix++)
    for (int iy = 0; iy < 2; iy++) {
        const double xi  = gauss[ix];
        const double eta = gauss[iy];

        double dN[2][4];
        for (int i = 0; i < 4; i++) {
            dN[0][i] = 0.25 * sx[i] * (1.0 + sy[i]*eta);
            dN[1][i] = 0.25 * sy[i] * (1.0 + sx[i]*xi);
        }

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

        double dNxy[2][4] = {};
        for (int r = 0; r < 2; r++)
        for (int n = 0; n < 4; n++)
        for (int k = 0; k < 2; k++)
            dNxy[r][n] += invJ[r][k] * dN[k][n];

        double B[3][8] = {};
        for (int n = 0; n < 4; n++) {
            int col = 2*n;
            const double Nx = dNxy[0][n], Ny = dNxy[1][n];
            B[0][col]   = Nx;
            B[1][col+1] = Ny;
            B[2][col]   = Ny; B[2][col+1] = Nx;
        }

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

        double KE[8][8];
        build_ke_2d(nu, KE);

        std::vector<bool> active_node(node_count, false);
        for (int ex = 0; ex < nelx; ex++)
        for (int ey = 0; ey < nely; ey++) {
            if (!mask_flat[ex*nely + ey]) continue;
            int dofs[8];
            element_dofs_2d(ex, ey, nny, dofs);
            for (int i = 0; i < 8; i++) active_node[dofs[i]/2] = true;
        }

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

        std::vector<double> x(elem_count, 0.0);
        for (int e = 0; e < elem_count; e++)
            if (mask_flat[e]) x[e] = vol_frac;

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

        int active_count = 0;
        for (int e = 0; e < elem_count; e++) if (mask_flat[e]) active_count++;

        // ===================================================================
        // Pre-build sparsity pattern once
        // ===================================================================
        Eigen::SparseMatrix<double, Eigen::RowMajor> Kff(nfree, nfree);
        {
            typedef Eigen::Triplet<double> T;
            std::vector<T> pat;
            pat.reserve((size_t)active_count * 64);
            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++) {
                const int eidx = ex*nely + ey;
                if (!mask_flat[eidx]) continue;
                int edofs[8];
                element_dofs_2d(ex, ey, nny, edofs);
                for (int i = 0; i < 8; i++) {
                    const int fi = dof_to_free[edofs[i]]; if (fi < 0) continue;
                    for (int j = 0; j < 8; j++) {
                        const int fj = dof_to_free[edofs[j]]; if (fj < 0) continue;
                        pat.emplace_back(fi, fj, 0.0);
                    }
                }
            }
            Kff.setFromTriplets(pat.begin(), pat.end());
        }
        Kff.makeCompressed();

        struct KEntry { int vidx; double ke; };
        std::vector<std::vector<KEntry>> elem_k(elem_count);
        {
            const int* outer = Kff.outerIndexPtr();
            const int* inner = Kff.innerIndexPtr();
            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++) {
                const int eidx = ex*nely + ey;
                if (!mask_flat[eidx]) continue;
                int edofs[8];
                element_dofs_2d(ex, ey, nny, edofs);
                auto& ev = elem_k[eidx];
                ev.reserve(64);
                for (int i = 0; i < 8; i++) {
                    const int fi = dof_to_free[edofs[i]]; if (fi < 0) continue;
                    for (int j = 0; j < 8; j++) {
                        const int fj = dof_to_free[edofs[j]]; if (fj < 0) continue;
                        int lo = outer[fi], hi = outer[fi + 1];
                        int pos = (int)(std::lower_bound(inner + lo, inner + hi, fj) - inner);
                        ev.push_back({pos, KE[i][j]});
                    }
                }
            }
        }

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

            double* kvals = Kff.valuePtr();
            std::fill(kvals, kvals + Kff.nonZeros(), 0.0);

            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++) {
                const int eidx = ex*nely + ey;
                if (!mask_flat[eidx]) continue;
                const double Ee = Emin + std::pow(x[eidx], penal) * (E0 - Emin);
                for (const auto& e : elem_k[eidx])
                    kvals[e.vidx] += Ee * e.ke;
            }

            AmgSolver amg(Kff, amg_prm);
            size_t amg_iters; double amg_error;
            std::tie(amg_iters, amg_error) = amg(Ff, Uf);

            Eigen::VectorXd U = Eigen::VectorXd::Zero(dof_count);
            for (int i = 0; i < nfree; i++) U(free_dofs[i]) = Uf(i);

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
