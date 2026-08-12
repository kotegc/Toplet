#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdarg>
#include <cstring>
#include <string>
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

#include "../include/Eigen/Sparse"
#include "../include/Eigen/SparseCholesky"
#include "../include/toplet_solver.h"

// ---------------------------------------------------------------------------
// Throwaway diagnostics: buffer collected during solve, shown in a MessageBox
// (auto-copied to clipboard) when solve finishes. No file artifact.
// ---------------------------------------------------------------------------
static thread_local std::string gmg_diag_buf;

static void gmg_log(const char* fmt, ...) {
    char buf[1024]; buf[sizeof(buf)-1] = 0;
    va_list args; va_start(args, fmt);
    vsnprintf(buf, sizeof(buf)-1, fmt, args);
    va_end(args);
    gmg_diag_buf += buf;
}
static void gmg_show_diag() {
    if (gmg_diag_buf.empty()) return;
    // Auto-copy to clipboard so user can paste it
    if (OpenClipboard(NULL)) {
        EmptyClipboard();
        size_t len = gmg_diag_buf.size() + 1;
        HGLOBAL hMem = GlobalAlloc(GMEM_MOVEABLE, len);
        if (hMem) {
            void* ptr = GlobalLock(hMem);
            if (ptr) { std::memcpy(ptr, gmg_diag_buf.c_str(), len); GlobalUnlock(hMem); }
            SetClipboardData(CF_TEXT, hMem);
        }
        CloseClipboard();
    }
    std::string msg = gmg_diag_buf + "\n(text copied to clipboard)";
    MessageBoxA(NULL, msg.c_str(), "GMG Diagnostics", MB_OK | MB_SETFOREGROUND);
    gmg_diag_buf.clear();
}

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

    double fac = 1.0 / ((1.0 + nu) * (1.0 - 2.0*nu));
    double a = (1.0 - nu) * fac;
    double b = nu * fac;
    double c = 0.5 * (1.0 - 2.0*nu) * fac;
    double C[6][6] = {};
    C[0][0]=a; C[0][1]=b; C[0][2]=b;
    C[1][0]=b; C[1][1]=a; C[1][2]=b;
    C[2][0]=b; C[2][1]=b; C[2][2]=a;
    C[3][3]=c; C[4][4]=c; C[5][5]=c;

    static const double xc[8] = {0,1,1,0, 0,1,1,0};
    static const double yc[8] = {0,0,1,1, 0,0,1,1};
    static const double zc[8] = {0,0,0,0, 1,1,1,1};
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

        double dN[3][8];
        for (int i = 0; i < 8; i++) {
            dN[0][i] = 0.125*sx[i]*(1.0+sy[i]*eta) *(1.0+sz[i]*zeta);
            dN[1][i] = 0.125*sy[i]*(1.0+sx[i]*xi)  *(1.0+sz[i]*zeta);
            dN[2][i] = 0.125*sz[i]*(1.0+sx[i]*xi)  *(1.0+sy[i]*eta);
        }

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

        double dNxyz[3][8] = {};
        for (int r = 0; r < 3; r++)
        for (int n = 0; n < 8; n++)
        for (int k = 0; k < 3; k++)
            dNxyz[r][n] += invJ[r][k] * dN[k][n];

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
// Geometric Multigrid (GMG) data structures and helpers
// ---------------------------------------------------------------------------
struct KEntry { int vidx; double ke; };

struct GmgLevel {
    int nelx, nely, nelz;
    int nfree;
    std::vector<int>    free_dofs;    // free_dofs[i]    = local DOF index
    std::vector<int>    dof_to_free;  // local DOF index → free DOF index (-1 if not free)
    std::vector<bool>   mask;         // active element flags
    std::vector<double> density;      // element densities

    Eigen::SparseMatrix<double, Eigen::RowMajor> K;     // stiffness (refilled each SIMP iter)
    Eigen::VectorXd                               diag;  // diagonal of K (for Jacobi smoother)
    Eigen::SparseMatrix<double, Eigen::RowMajor> P;     // prolongation: coarser → this level
    Eigen::SparseMatrix<double, Eigen::RowMajor> R;     // restriction = P^T

    std::vector<std::vector<KEntry>> elem_k;  // precomputed K value-array positions
};

// Coarse element mask: active if any of the 2x2x2 fine children is active.
static std::vector<bool> coarsen_mask(
    const std::vector<bool>& fm, int fnx, int fny, int fnz)
{
    int cnx=fnx/2, cny=fny/2, cnz=fnz/2;
    int cnyz=cny*cnz, fnyz=fny*fnz;
    std::vector<bool> cm((size_t)cnx*cny*cnz, false);
    for (int cx=0; cx<cnx; cx++)
    for (int cy=0; cy<cny; cy++)
    for (int cz=0; cz<cnz; cz++) {
        int ci=cx*cnyz+cy*cnz+cz;
        for (int dx=0; dx<2&&!cm[ci]; dx++)
        for (int dy=0; dy<2; dy++)
        for (int dz=0; dz<2; dz++)
            if (fm[(2*cx+dx)*fnyz+(2*cy+dy)*fnz+(2*cz+dz)]) cm[ci]=true;
    }
    return cm;
}

// Build free-DOF maps for one GMG level.
// scale = 2^l: a coarse node at (cx,cy,cz) corresponds to fine (level-0) node
// at (cx*scale, cy*scale, cz*scale).  fixed_dofs_fine is the original fine-grid
// fixed-DOF array; we inherit fixity from the matching fine node.
static void build_level_free_dofs(
    GmgLevel& L,
    const unsigned char* fixed_dofs_fine,
    int fine_nny, int fine_nnz,
    int scale)
{
    const int nnx=L.nelx+1, nny=L.nely+1, nnz=L.nelz+1;
    const int node_count=nnx*nny*nnz;
    const int nyz=L.nely*L.nelz;
    const int local_dof_count=node_count*3;

    std::vector<bool> active_node(node_count, false);
    for (int ex=0; ex<L.nelx; ex++)
    for (int ey=0; ey<L.nely; ey++)
    for (int ez=0; ez<L.nelz; ez++) {
        if (!L.mask[ex*nyz+ey*L.nelz+ez]) continue;
        int dofs[24]; element_dofs(ex,ey,ez,nny,nnz,dofs);
        for (int i=0; i<24; i++) active_node[dofs[i]/3]=true;
    }

    L.dof_to_free.assign(local_dof_count, -1);
    L.free_dofs.clear();
    for (int node=0; node<node_count; node++) {
        if (!active_node[node]) continue;
        int cix = node/(nny*nnz);
        int ciy = (node/nnz)%nny;
        int ciz = node%nnz;
        int fine_node = cix*scale*(fine_nny*fine_nnz) + ciy*scale*fine_nnz + ciz*scale;
        for (int d=0; d<3; d++) {
            int ldof = 3*node+d;
            int fdof = 3*fine_node+d;
            if (!fixed_dofs_fine[fdof]) {
                L.dof_to_free[ldof] = (int)L.free_dofs.size();
                L.free_dofs.push_back(ldof);
            }
        }
    }
    L.nfree = (int)L.free_dofs.size();
}

// Build the K sparsity pattern and the elem_k position cache for a level.
static void build_elem_k_cache(GmgLevel& L, const double KE[24][24])
{
    const int nny=L.nely+1, nnz=L.nelz+1, nyz=L.nely*L.nelz;
    int active_count=0;
    for (auto b : L.mask) if (b) active_count++;

    typedef Eigen::Triplet<double> T;
    std::vector<T> pat;
    pat.reserve((size_t)active_count*576);
    for (int ex=0; ex<L.nelx; ex++)
    for (int ey=0; ey<L.nely; ey++)
    for (int ez=0; ez<L.nelz; ez++) {
        int eidx=ex*nyz+ey*L.nelz+ez;
        if (!L.mask[eidx]) continue;
        int edofs[24]; element_dofs(ex,ey,ez,nny,nnz,edofs);
        for (int i=0; i<24; i++) { int fi=L.dof_to_free[edofs[i]]; if(fi<0) continue;
        for (int j=0; j<24; j++) { int fj=L.dof_to_free[edofs[j]]; if(fj<0) continue;
            pat.emplace_back(fi,fj,0.0); }}
    }
    L.K.resize(L.nfree, L.nfree);
    L.K.setFromTriplets(pat.begin(), pat.end());
    L.K.makeCompressed();

    const int* outer=L.K.outerIndexPtr();
    const int* inner=L.K.innerIndexPtr();
    L.elem_k.resize(L.mask.size());
    for (int ex=0; ex<L.nelx; ex++)
    for (int ey=0; ey<L.nely; ey++)
    for (int ez=0; ez<L.nelz; ez++) {
        int eidx=ex*nyz+ey*L.nelz+ez;
        if (!L.mask[eidx]) continue;
        int edofs[24]; element_dofs(ex,ey,ez,nny,nnz,edofs);
        auto& ev=L.elem_k[eidx]; ev.reserve(576);
        for (int i=0; i<24; i++) { int fi=L.dof_to_free[edofs[i]]; if(fi<0) continue;
        for (int j=0; j<24; j++) { int fj=L.dof_to_free[edofs[j]]; if(fj<0) continue;
            int lo=outer[fi],hi=outer[fi+1];
            int pos=(int)(std::lower_bound(inner+lo,inner+hi,fj)-inner);
            ev.push_back({pos,KE[i][j]}); }}
    }
}

// Trilinear prolongation P (nfree_fine x nfree_coarse): coarse → fine.
// Fine node at integer position (ix,iy,iz) interpolates from up to 8 surrounding
// coarse nodes.  Odd coordinates split weight 0.5/0.5; even coordinates weight 1.
static void build_prolongation(GmgLevel& fine, const GmgLevel& coarse)
{
    const int cnny=coarse.nely+1, cnnz=coarse.nelz+1;
    const int fnny=fine.nely+1,   fnnz=fine.nelz+1;

    typedef Eigen::Triplet<double> T;
    std::vector<T> trips;
    trips.reserve(fine.nfree*4);

    for (int fi=0; fi<fine.nfree; fi++) {
        int ldof=fine.free_dofs[fi];
        int d   =ldof%3;
        int node=ldof/3;
        int ix  =node/(fnny*fnnz);
        int iy  =(node/fnnz)%fnny;
        int iz  =node%fnnz;

        int    cx[2]={ix/2,(ix%2)?ix/2+1:ix/2};
        int    cy[2]={iy/2,(iy%2)?iy/2+1:iy/2};
        int    cz[2]={iz/2,(iz%2)?iz/2+1:iz/2};
        double wx[2]={(ix%2)?0.5:1.0,(ix%2)?0.5:0.0};
        double wy[2]={(iy%2)?0.5:1.0,(iy%2)?0.5:0.0};
        double wz[2]={(iz%2)?0.5:1.0,(iz%2)?0.5:0.0};

        for (int a=0;a<2;a++) for (int b=0;b<2;b++) for (int cc=0;cc<2;cc++) {
            double w=wx[a]*wy[b]*wz[cc];
            if (w<1e-12) continue;
            if (cx[a]<0||cx[a]>=coarse.nelx+1) continue;
            if (cy[b]<0||cy[b]>=coarse.nely+1) continue;
            if (cz[cc]<0||cz[cc]>=coarse.nelz+1) continue;
            int cnode=cx[a]*(cnny*cnnz)+cy[b]*cnnz+cz[cc];
            int cdof =3*cnode+d;
            if (cdof>=(int)coarse.dof_to_free.size()) continue;
            int cfi=coarse.dof_to_free[cdof];
            if (cfi<0) continue;
            trips.emplace_back(fi,cfi,w);
        }
    }
    fine.P.resize(fine.nfree, coarse.nfree);
    fine.P.setFromTriplets(trips.begin(), trips.end());
    fine.R=fine.P.transpose();
}

// Fill K values from current element densities and update the Jacobi diagonal.
static void assemble_level_K(GmgLevel& L, double E0, double Emin, double penal)
{
    double* kvals=L.K.valuePtr();
    std::fill(kvals, kvals+L.K.nonZeros(), 0.0);
    const int nyz=L.nely*L.nelz;
    for (int ex=0; ex<L.nelx; ex++)
    for (int ey=0; ey<L.nely; ey++)
    for (int ez=0; ez<L.nelz; ez++) {
        int eidx=ex*nyz+ey*L.nelz+ez;
        if (!L.mask[eidx]||L.elem_k[eidx].empty()) continue;
        double Ee=Emin+std::pow(L.density[eidx],penal)*(E0-Emin);
        for (const auto& e : L.elem_k[eidx]) kvals[e.vidx]+=Ee*e.ke;
    }
    L.diag.resize(L.nfree);
    for (int i=0; i<L.nfree; i++) {
        double dv=L.K.coeff(i,i);
        L.diag(i)=(dv>1e-15)?dv:1e-15;
    }
}

// Coarsen densities: each coarse element gets the average of its active fine children.
static void coarsen_densities(GmgLevel& coarse, const GmgLevel& fine)
{
    const int cnyz=coarse.nely*coarse.nelz;
    const int fnyz=fine.nely*fine.nelz;
    for (int cx=0; cx<coarse.nelx; cx++)
    for (int cy=0; cy<coarse.nely; cy++)
    for (int cz=0; cz<coarse.nelz; cz++) {
        int cidx=cx*cnyz+cy*coarse.nelz+cz;
        if (!coarse.mask[cidx]) continue;
        double sum=0.0; int cnt=0;
        for (int dx=0;dx<2;dx++) for (int dy=0;dy<2;dy++) for (int dz=0;dz<2;dz++) {
            int fidx=(2*cx+dx)*fnyz+(2*cy+dy)*fine.nelz+(2*cz+dz);
            if (fine.mask[fidx]) { sum+=fine.density[fidx]; cnt++; }
        }
        coarse.density[cidx]=(cnt>0)?sum/cnt:0.0;
    }
}

// V-cycle: pre-smooth (damped Jacobi) → restrict residual → coarse correction
//          (direct solve at coarsest) → prolongate → post-smooth.
static void vcycle(
    std::vector<GmgLevel>& levels,
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>& coarse_ldlt,
    int l,
    const Eigen::VectorXd& rhs,
    Eigen::VectorXd& x)
{
    const double omega=2.0/3.0;
    if (l==(int)levels.size()-1) {
        x=coarse_ldlt.solve(rhs);
        return;
    }
    GmgLevel& L=levels[l];

    for (int s=0;s<2;s++) { Eigen::VectorXd res=rhs-L.K*x; x+=omega*res.cwiseQuotient(L.diag); }

    Eigen::VectorXd rf=rhs-L.K*x;
    Eigen::VectorXd rc=L.R*rf;
    Eigen::VectorXd ec=Eigen::VectorXd::Zero(levels[l+1].nfree);
    vcycle(levels,coarse_ldlt,l+1,rc,ec);
    x+=L.P*ec;

    for (int s=0;s<2;s++) { Eigen::VectorXd res=rhs-L.K*x; x+=omega*res.cwiseQuotient(L.diag); }
}

// PCG with one GMG V-cycle as the preconditioner per CG iteration.
static void pcg_gmg(
    std::vector<GmgLevel>& levels,
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>& coarse_ldlt,
    const Eigen::VectorXd& rhs,
    Eigen::VectorXd& x,
    double tol, int maxiter,
    bool verbose_first = false)
{
    const GmgLevel& L0=levels[0];
    Eigen::VectorXd r=rhs-L0.K*x;
    Eigen::VectorXd z=Eigen::VectorXd::Zero(L0.nfree);
    vcycle(levels,coarse_ldlt,0,r,z);
    Eigen::VectorXd p=z;
    double rho=r.dot(z);
    double rhs_norm=rhs.norm();
    bool converged=false;
    int final_iter=maxiter;
    if (verbose_first) gmg_log("[PCG1] rho0=%.3e\n", rho);
    for (int i=0;i<maxiter;i++) {
        Eigen::VectorXd q=L0.K*p;
        double pq=p.dot(q);
        if (verbose_first && i<5) gmg_log("[PCG1] i=%d rho=%.3e pq=%.3e alpha=%.3e ||r||=%.3e\n",
                i, rho, pq, (std::abs(pq)>1e-30)?rho/pq:0.0, r.norm());
        if (std::abs(pq)<1e-30) { final_iter=i; break; }
        double alpha=rho/pq;
        x+=alpha*p;
        r-=alpha*q;
        if (r.norm()<tol*rhs_norm) { converged=true; final_iter=i+1; break; }
        z.setZero();
        vcycle(levels,coarse_ldlt,0,r,z);
        double rho_new=r.dot(z);
        if (std::abs(rho)<1e-30) { final_iter=i+1; break; }
        double beta=rho_new/rho;
        p=z+beta*p;
        rho=rho_new;
    }
    double final_resid = (rhs_norm>0) ? r.norm()/rhs_norm : 0.0;
    gmg_log("[PCG] %s: %d iters, rel_resid=%.3e, rhs_norm=%.3e\n",
            converged ? "OK" : "FAILED", final_iter, final_resid, rhs_norm);
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
    gmg_diag_buf.clear();

    try {
        const int nny = nely + 1;
        const int nnz = nelz + 1;
        const int nnx = nelx + 1;
        const int node_count = nnx * nny * nnz;
        const int dof_count  = node_count * 3;
        const int elem_count = nelx * nely * nelz;
        const int nyz = nely * nelz;

        double KE[24][24];
        build_ke(nu, KE);

        // --- Initialize element densities ---
        std::vector<double> x(elem_count, 0.0);
        for (int e = 0; e < elem_count; e++)
            if (mask_flat[e]) x[e] = vol_frac;

        // --- Sensitivity filter stencil ---
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

        int active_count = 0;
        for (int e = 0; e < elem_count; e++) if (mask_flat[e]) active_count++;

        // ===================================================================
        // GMG hierarchy setup (once — geometry is fixed across SIMP iterations)
        // 4 levels: L0=fine, L1, L2, L3=coarsest (direct solve via LDLT)
        // Dimensions halved each level; requires nelx/y/z divisible by 8.
        // ===================================================================
        const int num_levels = 4;
        std::vector<GmgLevel> gmg(num_levels);

        // Level 0: fine grid
        gmg[0].nelx=nelx; gmg[0].nely=nely; gmg[0].nelz=nelz;
        gmg[0].mask.assign(mask_flat, mask_flat+elem_count);
        gmg[0].density.assign(elem_count, vol_frac);
        build_level_free_dofs(gmg[0], fixed_dofs, nny, nnz, 1);
        build_elem_k_cache(gmg[0], KE);

        // Levels 1, 2, 3
        for (int l=1; l<num_levels; l++) {
            GmgLevel& prev=gmg[l-1];
            GmgLevel& cur =gmg[l];
            cur.nelx=prev.nelx/2; cur.nely=prev.nely/2; cur.nelz=prev.nelz/2;
            cur.mask=coarsen_mask(prev.mask, prev.nelx, prev.nely, prev.nelz);
            build_level_free_dofs(cur, fixed_dofs, nny, nnz, 1<<l);
            // No elem_k cache: coarse K computed by Galerkin (R * K_fine * P)
        }

        // Prolongation and restriction for each level pair (fixed geometry → build once)
        for (int l=0; l<num_levels-1; l++)
            build_prolongation(gmg[l], gmg[l+1]);

        // --- Initial fine K assembly (needed to bootstrap Galerkin coarsenings) ---
        assemble_level_K(gmg[0], E0, Emin, penal);

        // --- Galerkin coarsening: K_l = R_{l-1} * K_{l-1} * P_{l-1} ---
        // Guaranteed SPD because K_0 is SPD and P has full column rank.
        // Eliminates rigid-body modes that plague rediscretization on thin/coarse grids.
        for (int l=1; l<num_levels; l++) {
            Eigen::SparseMatrix<double> Kf_cm(gmg[l-1].K);   // ColMajor copy
            Eigen::SparseMatrix<double> P_cm(gmg[l-1].P);    // ColMajor copy
            Eigen::SparseMatrix<double> R_cm(gmg[l-1].R);    // ColMajor copy
            Eigen::SparseMatrix<double> Kl_cm = R_cm * (Kf_cm * P_cm);
            gmg[l].K = Kl_cm;   // convert result to RowMajor
            gmg[l].diag.resize(gmg[l].nfree);
            for (int i=0; i<gmg[l].nfree; i++) {
                double dv = gmg[l].K.coeff(i,i);
                gmg[l].diag(i) = std::max(dv, 1e-15);
            }
        }

        // Symbolic factorization of coarsest Galerkin K (sparsity fixed: depends only on R,P pattern)
        Eigen::SparseMatrix<double> Kc_col;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> coarse_ldlt;
        Kc_col = gmg[num_levels-1].K;
        coarse_ldlt.analyzePattern(Kc_col);

        // Diagnostic: report hierarchy sizes so we can see if coarsest level is degenerate
        for (int l=0; l<num_levels; l++)
            gmg_log("[GMG] Level %d: %dx%dx%d, nfree=%d\n",
                    l, gmg[l].nelx, gmg[l].nely, gmg[l].nelz, gmg[l].nfree);

        const int nfree = gmg[0].nfree;
        const std::vector<int>& free_dofs = gmg[0].free_dofs;

        // Constant RHS (forces don't change between iterations)
        Eigen::VectorXd Ff(nfree);
        for (int i = 0; i < nfree; i++) Ff(i) = forces[free_dofs[i]];

        Eigen::VectorXd Uf = Eigen::VectorXd::Zero(nfree);

        // --- One-time V-cycle quality test at uniform initial densities ---
        // Fine K and all Galerkin coarse K already assembled in setup above; just factorize.
        {
            Kc_col = gmg[num_levels-1].K;
            coarse_ldlt.factorize(Kc_col);

            gmg_log("[VCTEST] fine K diag min=%.3e max=%.3e\n",
                    gmg[0].diag.minCoeff(), gmg[0].diag.maxCoeff());
            gmg_log("[VCTEST] coarsest LDLT min D = %.3e (negative=indefinite)\n",
                    coarse_ldlt.vectorD().minCoeff());

            Eigen::VectorXd z_vc = Eigen::VectorXd::Zero(nfree);
            vcycle(gmg, coarse_ldlt, 0, Ff, z_vc);
            double rho_test    = Ff.dot(z_vc);
            double z_norm      = z_vc.norm();
            double ff_norm     = Ff.norm() > 0 ? Ff.norm() : 1.0;
            double resid_after = (Ff - gmg[0].K * z_vc).norm() / ff_norm;
            gmg_log("[VCTEST] z_norm=%.3e, rho=Ff.z=%.3e, resid_after_vc=%.3e\n",
                    z_norm, rho_test, resid_after);
            gmg_log("[VCTEST] (resid_after > 1 means V-cycle actively makes things worse)\n");
        }

        // ===================================================================
        // Main SIMP iteration loop
        // ===================================================================
        for (int iter = 0; iter < max_iter; iter++) {

            // --- Assemble fine K then Galerkin-coarsen all other levels ---
            gmg[0].density = x;
            assemble_level_K(gmg[0], E0, Emin, penal);
            for (int l=1; l<num_levels; l++) {
                Eigen::SparseMatrix<double> Kf_cm(gmg[l-1].K);
                Eigen::SparseMatrix<double> P_cm(gmg[l-1].P);
                Eigen::SparseMatrix<double> R_cm(gmg[l-1].R);
                Eigen::SparseMatrix<double> Kl_cm = R_cm * (Kf_cm * P_cm);
                gmg[l].K = Kl_cm;
                for (int i=0; i<gmg[l].nfree; i++) {
                    double dv = gmg[l].K.coeff(i,i);
                    gmg[l].diag(i) = std::max(dv, 1e-15);
                }
            }

            // --- Refactorize coarsest Galerkin K (guaranteed SPD) ---
            Kc_col = gmg[num_levels-1].K;
            coarse_ldlt.factorize(Kc_col);
            if (coarse_ldlt.info() != Eigen::Success)
                gmg_log("[GMG] WARNING: coarsest Galerkin LDLT failed at SIMP iter %d\n", iter+1);

            // Start from zero to avoid contamination from a previous bad Uf
            Uf.setZero();
            // --- PCG with GMG V-cycle preconditioner ---
            pcg_gmg(gmg, coarse_ldlt, Ff, Uf, 1e-3, 300, iter == 0);

            // --- Expand to full displacement vector ---
            Eigen::VectorXd U = Eigen::VectorXd::Zero(dof_count);
            for (int i = 0; i < nfree; i++) U(free_dofs[i]) = Uf(i);

            // --- Compliance and sensitivities ---
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
                    if (i2<0||i2>=nelx||j2<0||j2>=nely||k2<0||k2>=nelz) continue;
                    const int nidx = i2*nyz + j2*nelz + k2;
                    if (!mask_flat[nidx]) continue;
                    sum_w += s.weight;
                    val   += s.weight * x[nidx] * dc[nidx];
                }
                const double denom = x[eidx] * sum_w;
                dc_filt[eidx] = (denom > 1e-12) ? val / denom : dc[eidx];
            }

            // --- OC density update with bisection ---
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
                    double bv = -dc_filt[eidx] / lmid;
                    if (bv < 1e-12) bv = 1e-12;
                    double cand = rho * std::sqrt(bv);
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
                *compliance_out = compliance;
                *iterations_out = iter + 1;
                gmg_show_diag();
                return 0;
            }
        }

        std::copy(x.begin(), x.end(), density_out);
        *compliance_out = compliance;
        *iterations_out = max_iter;
        gmg_show_diag();
        return 0;
    }
    catch (...) {
        gmg_show_diag();
        return 1;
    }
}
