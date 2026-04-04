using System;
using System.Collections.Generic;

namespace Toplet_v0_Alpha.TopOpt3D
{
    public class TopOptSolver3D
    {
        public TopOptResult3D Solve(TopOptProblem3D problem, TopOptDomain3D domain)
        {
            if (problem == null)
                throw new ArgumentNullException(nameof(problem));
            if (domain == null)
                throw new ArgumentNullException(nameof(domain));

            int nelx = problem.NelX;
            int nely = problem.NelY;
            int nelz = problem.NelZ;

            if (nelx < 1 || nely < 1 || nelz < 1)
                throw new ArgumentException("NelX, NelY, and NelZ must all be >= 1.");

            double volfrac = problem.VolumeFraction;
            double penal = problem.Penal;
            double rmin = problem.FilterRadius;
            int maxIter = problem.MaxIterations;

            double E0 = problem.YoungsModulusSolid;
            double Emin = problem.YoungsModulusMin;
            double nu = problem.PoissonRatio;

            int nnx = nelx + 1;
            int nny = nely + 1;
            int nnz = nelz + 1;
            int nodeCount = nnx * nny * nnz;
            int dofCount = nodeCount * 3;

            ValidateInputs(problem, domain, dofCount);

            bool[,,] mask = domain.DesignMask;
            double[] forces = domain.Forces;
            bool[] isFixed = domain.FixedDofs;

            int activeElementCount = CountActiveElements(mask, nelx, nely, nelz);
            if (activeElementCount == 0)
                throw new InvalidOperationException("Design mask contains no active elements.");

            bool[] activeNodes = BuildActiveNodes(mask, nelx, nely, nelz, nny, nnz);

            if (!HasLoadOnActiveNode(forces, activeNodes))
                throw new InvalidOperationException("No load is applied to any active node.");

            if (!HasFixedDofOnActiveNode(isFixed, activeNodes))
                throw new InvalidOperationException("No fixed DOFs were found on active nodes.");

            List<int> freeDofs = BuildFreeDofs(isFixed, activeNodes);
            if (freeDofs.Count == 0)
                throw new InvalidOperationException("No free DOFs found on active connected nodes.");

            double[,,] x = new double[nelx, nely, nelz];
            double[,,] xNew = new double[nelx, nely, nelz];
            double[,,] dc = new double[nelx, nely, nelz];
            double[,,] dcFiltered = new double[nelx, nely, nelz];

            InitializeDensity(x, mask, nelx, nely, nelz, volfrac);

            double[,] KE = BuildElementStiffness(nu);

            double compliance = 0.0;

            for (int iter = 0; iter < maxIter; iter++)
            {
                double[] U = SolveDisplacements(
                    nelx,
                    nely,
                    nelz,
                    nny,
                    nnz,
                    x,
                    mask,
                    penal,
                    E0,
                    Emin,
                    KE,
                    forces,
                    freeDofs);

                compliance = 0.0;
                ClearArray(dc, nelx, nely, nelz);
                ClearArray(dcFiltered, nelx, nely, nelz);

                for (int ex = 0; ex < nelx; ex++)
                {
                    for (int ey = 0; ey < nely; ey++)
                    {
                        for (int ez = 0; ez < nelz; ez++)
                        {
                            if (!mask[ex, ey, ez])
                                continue;

                            int[] edofs = ElementDofs(ex, ey, ez, nny, nnz);

                            double[] ue = new double[24];
                            for (int i = 0; i < 24; i++)
                                ue[i] = U[edofs[i]];

                            double ce = Dot(ue, Multiply(KE, ue));

                            double rho = x[ex, ey, ez];
                            double Ee = Emin + Math.Pow(rho, penal) * (E0 - Emin);

                            compliance += Ee * ce;
                            dc[ex, ey, ez] = -penal * Math.Pow(rho, penal - 1.0) * (E0 - Emin) * ce;
                        }
                    }
                }

                ApplySensitivityFilter(nelx, nely, nelz, rmin, mask, x, dc, dcFiltered);
                UpdateDensitiesOC(nelx, nely, nelz, volfrac, mask, x, xNew, dcFiltered);

                double change = 0.0;

                for (int ex = 0; ex < nelx; ex++)
                {
                    for (int ey = 0; ey < nely; ey++)
                    {
                        for (int ez = 0; ez < nelz; ez++)
                        {
                            if (!mask[ex, ey, ez])
                            {
                                x[ex, ey, ez] = 0.0;
                                xNew[ex, ey, ez] = 0.0;
                                continue;
                            }

                            double diff = Math.Abs(xNew[ex, ey, ez] - x[ex, ey, ez]);
                            if (diff > change)
                                change = diff;

                            x[ex, ey, ez] = xNew[ex, ey, ez];
                        }
                    }
                }

                if (change < 0.01 && iter > 10)
                {
                    return new TopOptResult3D
                    {
                        NelX = nelx,
                        NelY = nely,
                        NelZ = nelz,
                        Density = x,
                        Compliance = compliance,
                        Iterations = iter + 1
                    };
                }
            }

            return new TopOptResult3D
            {
                NelX = nelx,
                NelY = nely,
                NelZ = nelz,
                Density = x,
                Compliance = compliance,
                Iterations = maxIter
            };
        }

        private static void ValidateInputs(TopOptProblem3D problem, TopOptDomain3D domain, int dofCount)
        {
            if (domain.DesignMask == null)
                throw new ArgumentException("DesignMask cannot be null.");
            if (domain.FixedDofs == null)
                throw new ArgumentException("FixedDofs cannot be null.");
            if (domain.Forces == null)
                throw new ArgumentException("Forces cannot be null.");

            if (domain.DesignMask.GetLength(0) != problem.NelX ||
                domain.DesignMask.GetLength(1) != problem.NelY ||
                domain.DesignMask.GetLength(2) != problem.NelZ)
            {
                throw new ArgumentException("DesignMask size must match NelX x NelY x NelZ.");
            }

            if (domain.FixedDofs.Length != dofCount)
                throw new ArgumentException("FixedDofs length must match total DOF count.");

            if (domain.Forces.Length != dofCount)
                throw new ArgumentException("Forces length must match total DOF count.");
        }

        private static int CountActiveElements(bool[,,] mask, int nelx, int nely, int nelz)
        {
            int count = 0;

            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    for (int ez = 0; ez < nelz; ez++)
                    {
                        if (mask[ex, ey, ez])
                            count++;
                    }
                }
            }

            return count;
        }

        private static void InitializeDensity(double[,,] x, bool[,,] mask, int nelx, int nely, int nelz, double volfrac)
        {
            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    for (int ez = 0; ez < nelz; ez++)
                    {
                        x[ex, ey, ez] = mask[ex, ey, ez] ? volfrac : 0.0;
                    }
                }
            }
        }

        private static void ClearArray(double[,,] a, int nelx, int nely, int nelz)
        {
            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    for (int ez = 0; ez < nelz; ez++)
                    {
                        a[ex, ey, ez] = 0.0;
                    }
                }
            }
        }

        private static bool[] BuildActiveNodes(bool[,,] mask, int nelx, int nely, int nelz, int nny, int nnz)
        {
            int nodeCount = (nelx + 1) * (nely + 1) * (nelz + 1);
            bool[] activeNodes = new bool[nodeCount];

            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    for (int ez = 0; ez < nelz; ez++)
                    {
                        if (!mask[ex, ey, ez])
                            continue;

                        int n1 = NodeIndex(ex, ey, ez, nny, nnz);
                        int n2 = NodeIndex(ex + 1, ey, ez, nny, nnz);
                        int n3 = NodeIndex(ex + 1, ey + 1, ez, nny, nnz);
                        int n4 = NodeIndex(ex, ey + 1, ez, nny, nnz);
                        int n5 = NodeIndex(ex, ey, ez + 1, nny, nnz);
                        int n6 = NodeIndex(ex + 1, ey, ez + 1, nny, nnz);
                        int n7 = NodeIndex(ex + 1, ey + 1, ez + 1, nny, nnz);
                        int n8 = NodeIndex(ex, ey + 1, ez + 1, nny, nnz);

                        activeNodes[n1] = true;
                        activeNodes[n2] = true;
                        activeNodes[n3] = true;
                        activeNodes[n4] = true;
                        activeNodes[n5] = true;
                        activeNodes[n6] = true;
                        activeNodes[n7] = true;
                        activeNodes[n8] = true;
                    }
                }
            }

            return activeNodes;
        }

        private static bool HasLoadOnActiveNode(double[] forces, bool[] activeNodes)
        {
            int nodeCount = activeNodes.Length;

            for (int node = 0; node < nodeCount; node++)
            {
                if (!activeNodes[node])
                    continue;

                int dofX = 3 * node;
                int dofY = 3 * node + 1;
                int dofZ = 3 * node + 2;

                if (Math.Abs(forces[dofX]) > 1e-12 ||
                    Math.Abs(forces[dofY]) > 1e-12 ||
                    Math.Abs(forces[dofZ]) > 1e-12)
                {
                    return true;
                }
            }

            return false;
        }

        private static bool HasFixedDofOnActiveNode(bool[] isFixed, bool[] activeNodes)
        {
            int nodeCount = activeNodes.Length;

            for (int node = 0; node < nodeCount; node++)
            {
                if (!activeNodes[node])
                    continue;

                int dofX = 3 * node;
                int dofY = 3 * node + 1;
                int dofZ = 3 * node + 2;

                if (isFixed[dofX] || isFixed[dofY] || isFixed[dofZ])
                    return true;
            }

            return false;
        }

        private static List<int> BuildFreeDofs(bool[] isFixed, bool[] activeNodes)
        {
            List<int> freeDofs = new List<int>();

            int nodeCount = activeNodes.Length;

            for (int node = 0; node < nodeCount; node++)
            {
                if (!activeNodes[node])
                    continue;

                int dofX = 3 * node;
                int dofY = 3 * node + 1;
                int dofZ = 3 * node + 2;

                if (!isFixed[dofX]) freeDofs.Add(dofX);
                if (!isFixed[dofY]) freeDofs.Add(dofY);
                if (!isFixed[dofZ]) freeDofs.Add(dofZ);
            }

            return freeDofs;
        }

        private double[] SolveDisplacements(
            int nelx,
            int nely,
            int nelz,
            int nny,
            int nnz,
            double[,,] x,
            bool[,,] mask,
            double penal,
            double E0,
            double Emin,
            double[,] KE,
            double[] F,
            List<int> freeDofs)
        {
            int nodeCount = (nelx + 1) * (nely + 1) * (nelz + 1);
            int dofCount = nodeCount * 3;

            double[,] K = new double[dofCount, dofCount];

            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    for (int ez = 0; ez < nelz; ez++)
                    {
                        if (!mask[ex, ey, ez])
                            continue;

                        double rho = x[ex, ey, ez];
                        double Ee = Emin + Math.Pow(rho, penal) * (E0 - Emin);

                        int[] edofs = ElementDofs(ex, ey, ez, nny, nnz);

                        for (int i = 0; i < 24; i++)
                        {
                            int ii = edofs[i];

                            for (int j = 0; j < 24; j++)
                            {
                                int jj = edofs[j];
                                K[ii, jj] += Ee * KE[i, j];
                            }
                        }
                    }
                }
            }

            int nfree = freeDofs.Count;
            double[,] Kff = new double[nfree, nfree];
            double[] Ff = new double[nfree];

            for (int i = 0; i < nfree; i++)
            {
                int gi = freeDofs[i];
                Ff[i] = F[gi];

                for (int j = 0; j < nfree; j++)
                {
                    int gj = freeDofs[j];
                    Kff[i, j] = K[gi, gj];
                }
            }

            double[] Uf = SolveLinearSystem(Kff, Ff);

            double[] U = new double[dofCount];
            for (int i = 0; i < nfree; i++)
                U[freeDofs[i]] = Uf[i];

            return U;
        }

        private static int NodeIndex(int ix, int iy, int iz, int nny, int nnz)
        {
            return ix * (nny * nnz) + iy * nnz + iz;
        }

        private static int[] ElementDofs(int ex, int ey, int ez, int nny, int nnz)
        {
            int n1 = NodeIndex(ex, ey, ez, nny, nnz);
            int n2 = NodeIndex(ex + 1, ey, ez, nny, nnz);
            int n3 = NodeIndex(ex + 1, ey + 1, ez, nny, nnz);
            int n4 = NodeIndex(ex, ey + 1, ez, nny, nnz);
            int n5 = NodeIndex(ex, ey, ez + 1, nny, nnz);
            int n6 = NodeIndex(ex + 1, ey, ez + 1, nny, nnz);
            int n7 = NodeIndex(ex + 1, ey + 1, ez + 1, nny, nnz);
            int n8 = NodeIndex(ex, ey + 1, ez + 1, nny, nnz);

            return new int[]
            {
                3 * n1, 3 * n1 + 1, 3 * n1 + 2,
                3 * n2, 3 * n2 + 1, 3 * n2 + 2,
                3 * n3, 3 * n3 + 1, 3 * n3 + 2,
                3 * n4, 3 * n4 + 1, 3 * n4 + 2,
                3 * n5, 3 * n5 + 1, 3 * n5 + 2,
                3 * n6, 3 * n6 + 1, 3 * n6 + 2,
                3 * n7, 3 * n7 + 1, 3 * n7 + 2,
                3 * n8, 3 * n8 + 1, 3 * n8 + 2
            };
        }

        private static double[,] BuildElementStiffness(double nu)
        {
            double[,] ke = new double[24, 24];

            double[,] C = BuildElasticityMatrix(nu);

            double g = 1.0 / Math.Sqrt(3.0);
            double[] gauss = new double[] { -g, g };

            for (int ix = 0; ix < 2; ix++)
            {
                for (int iy = 0; iy < 2; iy++)
                {
                    for (int iz = 0; iz < 2; iz++)
                    {
                        double xi = gauss[ix];
                        double eta = gauss[iy];
                        double zeta = gauss[iz];

                        double[,] dN_dNatural = ShapeFunctionDerivativesNatural(xi, eta, zeta);
                        double[,] J = BuildJacobian(dN_dNatural);
                        double detJ = Determinant3x3(J);

                        if (detJ <= 0.0)
                            throw new InvalidOperationException("Invalid hexahedral element Jacobian.");

                        double[,] invJ = Inverse3x3(J);
                        double[,] dN_dXYZ = Multiply3x3_3x8(invJ, dN_dNatural);
                        double[,] B = BuildBMatrix(dN_dXYZ);

                        double[,] Bt = Transpose(B);
                        double[,] BtC = Multiply(Bt, C);
                        double[,] BtCB = Multiply(BtC, B);

                        for (int i = 0; i < 24; i++)
                        {
                            for (int j = 0; j < 24; j++)
                            {
                                ke[i, j] += BtCB[i, j] * detJ;
                            }
                        }
                    }
                }
            }

            return ke;
        }

        private static double[,] BuildElasticityMatrix(double nu)
        {
            double factor = 1.0 / ((1.0 + nu) * (1.0 - 2.0 * nu));
            double a = 1.0 - nu;
            double b = nu;
            double c = 0.5 * (1.0 - 2.0 * nu);

            return new double[,]
            {
                { a, b, b, 0, 0, 0 },
                { b, a, b, 0, 0, 0 },
                { b, b, a, 0, 0, 0 },
                { 0, 0, 0, c, 0, 0 },
                { 0, 0, 0, 0, c, 0 },
                { 0, 0, 0, 0, 0, c }
            }.Scale(factor);
        }

        private static double[,] ShapeFunctionDerivativesNatural(double xi, double eta, double zeta)
        {
            // rows: d/dxi, d/deta, d/dzeta
            // cols: node 1..8
            double[,] dN = new double[3, 8];

            double[] sx = { -1, 1, 1, -1, -1, 1, 1, -1 };
            double[] sy = { -1, -1, 1, 1, -1, -1, 1, 1 };
            double[] sz = { -1, -1, -1, -1, 1, 1, 1, 1 };

            for (int i = 0; i < 8; i++)
            {
                dN[0, i] = 0.125 * sx[i] * (1.0 + sy[i] * eta) * (1.0 + sz[i] * zeta);
                dN[1, i] = 0.125 * sy[i] * (1.0 + sx[i] * xi) * (1.0 + sz[i] * zeta);
                dN[2, i] = 0.125 * sz[i] * (1.0 + sx[i] * xi) * (1.0 + sy[i] * eta);
            }

            return dN;
        }

        private static double[,] BuildJacobian(double[,] dN_dNatural)
        {
            // Unit cube element with node coords:
            // 1:(0,0,0) 2:(1,0,0) 3:(1,1,0) 4:(0,1,0)
            // 5:(0,0,1) 6:(1,0,1) 7:(1,1,1) 8:(0,1,1)
            double[] x = { 0, 1, 1, 0, 0, 1, 1, 0 };
            double[] y = { 0, 0, 1, 1, 0, 0, 1, 1 };
            double[] z = { 0, 0, 0, 0, 1, 1, 1, 1 };

            double[,] J = new double[3, 3];

            for (int i = 0; i < 8; i++)
            {
                J[0, 0] += dN_dNatural[0, i] * x[i];
                J[0, 1] += dN_dNatural[0, i] * y[i];
                J[0, 2] += dN_dNatural[0, i] * z[i];

                J[1, 0] += dN_dNatural[1, i] * x[i];
                J[1, 1] += dN_dNatural[1, i] * y[i];
                J[1, 2] += dN_dNatural[1, i] * z[i];

                J[2, 0] += dN_dNatural[2, i] * x[i];
                J[2, 1] += dN_dNatural[2, i] * y[i];
                J[2, 2] += dN_dNatural[2, i] * z[i];
            }

            return J;
        }

        private static double Determinant3x3(double[,] A)
        {
            return
                A[0, 0] * (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]) -
                A[0, 1] * (A[1, 0] * A[2, 2] - A[1, 2] * A[2, 0]) +
                A[0, 2] * (A[1, 0] * A[2, 1] - A[1, 1] * A[2, 0]);
        }

        private static double[,] Inverse3x3(double[,] A)
        {
            double det = Determinant3x3(A);
            if (Math.Abs(det) < 1e-14)
                throw new InvalidOperationException("3x3 matrix is singular.");

            double invDet = 1.0 / det;
            double[,] inv = new double[3, 3];

            inv[0, 0] = (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]) * invDet;
            inv[0, 1] = -(A[0, 1] * A[2, 2] - A[0, 2] * A[2, 1]) * invDet;
            inv[0, 2] = (A[0, 1] * A[1, 2] - A[0, 2] * A[1, 1]) * invDet;

            inv[1, 0] = -(A[1, 0] * A[2, 2] - A[1, 2] * A[2, 0]) * invDet;
            inv[1, 1] = (A[0, 0] * A[2, 2] - A[0, 2] * A[2, 0]) * invDet;
            inv[1, 2] = -(A[0, 0] * A[1, 2] - A[0, 2] * A[1, 0]) * invDet;

            inv[2, 0] = (A[1, 0] * A[2, 1] - A[1, 1] * A[2, 0]) * invDet;
            inv[2, 1] = -(A[0, 0] * A[2, 1] - A[0, 1] * A[2, 0]) * invDet;
            inv[2, 2] = (A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]) * invDet;

            return inv;
        }

        private static double[,] Multiply3x3_3x8(double[,] A, double[,] B)
        {
            double[,] C = new double[3, 8];

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    double sum = 0.0;
                    for (int k = 0; k < 3; k++)
                        sum += A[i, k] * B[k, j];
                    C[i, j] = sum;
                }
            }

            return C;
        }

        private static double[,] BuildBMatrix(double[,] dN)
        {
            // dN rows: dN/dx, dN/dy, dN/dz
            double[,] B = new double[6, 24];

            for (int i = 0; i < 8; i++)
            {
                int c = 3 * i;

                double dNx = dN[0, i];
                double dNy = dN[1, i];
                double dNz = dN[2, i];

                B[0, c] = dNx;
                B[1, c + 1] = dNy;
                B[2, c + 2] = dNz;

                B[3, c] = dNy;
                B[3, c + 1] = dNx;

                B[4, c + 1] = dNz;
                B[4, c + 2] = dNy;

                B[5, c] = dNz;
                B[5, c + 2] = dNx;
            }

            return B;
        }

        private static void ApplySensitivityFilter(
            int nelx,
            int nely,
            int nelz,
            double rmin,
            bool[,,] mask,
            double[,,] x,
            double[,,] dc,
            double[,,] dcFiltered)
        {
            int r = (int)Math.Floor(rmin);

            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    for (int ez = 0; ez < nelz; ez++)
                    {
                        if (!mask[ex, ey, ez])
                        {
                            dcFiltered[ex, ey, ez] = 0.0;
                            continue;
                        }

                        double sum = 0.0;
                        double val = 0.0;

                        int iMin = Math.Max(ex - r, 0);
                        int iMax = Math.Min(ex + r, nelx - 1);
                        int jMin = Math.Max(ey - r, 0);
                        int jMax = Math.Min(ey + r, nely - 1);
                        int kMin = Math.Max(ez - r, 0);
                        int kMax = Math.Min(ez + r, nelz - 1);

                        for (int i = iMin; i <= iMax; i++)
                        {
                            for (int j = jMin; j <= jMax; j++)
                            {
                                for (int k = kMin; k <= kMax; k++)
                                {
                                    if (!mask[i, j, k])
                                        continue;

                                    double dx = ex - i;
                                    double dy = ey - j;
                                    double dz = ez - k;
                                    double dist = Math.Sqrt(dx * dx + dy * dy + dz * dz);
                                    double weight = rmin - dist;

                                    if (weight > 0.0)
                                    {
                                        sum += weight;
                                        val += weight * x[i, j, k] * dc[i, j, k];
                                    }
                                }
                            }
                        }

                        double denom = x[ex, ey, ez] * sum;
                        if (denom > 1e-12)
                            dcFiltered[ex, ey, ez] = val / denom;
                        else
                            dcFiltered[ex, ey, ez] = dc[ex, ey, ez];
                    }
                }
            }
        }

        private static void UpdateDensitiesOC(
            int nelx,
            int nely,
            int nelz,
            double volfrac,
            bool[,,] mask,
            double[,,] x,
            double[,,] xNew,
            double[,,] dc)
        {
            int activeCount = CountActiveElements(mask, nelx, nely, nelz);
            if (activeCount == 0)
                throw new InvalidOperationException("No active elements available for OC update.");

            double move = 0.2;
            double xmin = 0.001;
            double l1 = 0.0;
            double l2 = 1e9;
            double targetVolume = volfrac * activeCount;

            while ((l2 - l1) / (l1 + l2 + 1e-12) > 1e-4)
            {
                double lmid = 0.5 * (l1 + l2);
                double total = 0.0;

                for (int ex = 0; ex < nelx; ex++)
                {
                    for (int ey = 0; ey < nely; ey++)
                    {
                        for (int ez = 0; ez < nelz; ez++)
                        {
                            if (!mask[ex, ey, ez])
                            {
                                xNew[ex, ey, ez] = 0.0;
                                continue;
                            }

                            double rho = x[ex, ey, ez];
                            double b = -dc[ex, ey, ez] / lmid;
                            if (b < 1e-12)
                                b = 1e-12;

                            double candidate = rho * Math.Sqrt(b);

                            double lower = Math.Max(xmin, rho - move);
                            double upper = Math.Min(1.0, rho + move);

                            if (candidate < lower) candidate = lower;
                            if (candidate > upper) candidate = upper;

                            xNew[ex, ey, ez] = candidate;
                            total += candidate;
                        }
                    }
                }

                if (total > targetVolume)
                    l1 = lmid;
                else
                    l2 = lmid;
            }
        }

        private static double[] Multiply(double[,] A, double[] x)
        {
            int rows = A.GetLength(0);
            int cols = A.GetLength(1);

            if (x.Length != cols)
                throw new ArgumentException("Matrix/vector size mismatch.");

            double[] y = new double[rows];

            for (int i = 0; i < rows; i++)
            {
                double sum = 0.0;
                for (int j = 0; j < cols; j++)
                    sum += A[i, j] * x[j];
                y[i] = sum;
            }

            return y;
        }

        private static double Dot(double[] a, double[] b)
        {
            if (a.Length != b.Length)
                throw new ArgumentException("Vector size mismatch.");

            double sum = 0.0;
            for (int i = 0; i < a.Length; i++)
                sum += a[i] * b[i];
            return sum;
        }

        private static double[] SolveLinearSystem(double[,] A, double[] b)
        {
            int n = b.Length;

            if (A.GetLength(0) != n || A.GetLength(1) != n)
                throw new ArgumentException("Matrix must be square and match rhs length.");

            double[,] M = new double[n, n];
            double[] rhs = new double[n];

            for (int i = 0; i < n; i++)
            {
                rhs[i] = b[i];
                for (int j = 0; j < n; j++)
                    M[i, j] = A[i, j];
            }

            for (int k = 0; k < n; k++)
            {
                int pivotRow = k;
                double pivotMax = Math.Abs(M[k, k]);

                for (int i = k + 1; i < n; i++)
                {
                    double val = Math.Abs(M[i, k]);
                    if (val > pivotMax)
                    {
                        pivotMax = val;
                        pivotRow = i;
                    }
                }

                if (pivotMax < 1e-14)
                    throw new InvalidOperationException("Linear solve failed: singular or near-singular matrix.");

                if (pivotRow != k)
                {
                    for (int j = k; j < n; j++)
                    {
                        double tmp = M[k, j];
                        M[k, j] = M[pivotRow, j];
                        M[pivotRow, j] = tmp;
                    }

                    double tmpRhs = rhs[k];
                    rhs[k] = rhs[pivotRow];
                    rhs[pivotRow] = tmpRhs;
                }

                double pivot = M[k, k];

                for (int i = k + 1; i < n; i++)
                {
                    double factor = M[i, k] / pivot;
                    if (Math.Abs(factor) < 1e-20)
                        continue;

                    M[i, k] = 0.0;

                    for (int j = k + 1; j < n; j++)
                        M[i, j] -= factor * M[k, j];

                    rhs[i] -= factor * rhs[k];
                }
            }

            double[] x = new double[n];
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = rhs[i];
                for (int j = i + 1; j < n; j++)
                    sum -= M[i, j] * x[j];

                x[i] = sum / M[i, i];
            }

            return x;
        }

        private static double[,] Transpose(double[,] A)
        {
            int rows = A.GetLength(0);
            int cols = A.GetLength(1);
            double[,] T = new double[cols, rows];

            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    T[j, i] = A[i, j];

            return T;
        }

        private static double[,] Multiply(double[,] A, double[,] B)
        {
            int aRows = A.GetLength(0);
            int aCols = A.GetLength(1);
            int bRows = B.GetLength(0);
            int bCols = B.GetLength(1);

            if (aCols != bRows)
                throw new ArgumentException("Matrix size mismatch.");

            double[,] C = new double[aRows, bCols];

            for (int i = 0; i < aRows; i++)
            {
                for (int j = 0; j < bCols; j++)
                {
                    double sum = 0.0;
                    for (int k = 0; k < aCols; k++)
                        sum += A[i, k] * B[k, j];
                    C[i, j] = sum;
                }
            }

            return C;
        }
    }

    internal static class MatrixExtensions3D
    {
        public static double[,] Scale(this double[,] A, double factor)
        {
            int rows = A.GetLength(0);
            int cols = A.GetLength(1);
            double[,] B = new double[rows, cols];

            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    B[i, j] = A[i, j] * factor;

            return B;
        }
    }
}