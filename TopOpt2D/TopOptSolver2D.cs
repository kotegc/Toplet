using System;
using System.Collections.Generic;

namespace Toplet_v0_Alpha.TopOpt2D
{
    public class TopOptSolver2D
    {
        public TopOptResult2D Solve(TopOptProblem2D problem, TopOptDomain2D domain)
        {
            if (problem == null)
                throw new ArgumentNullException(nameof(problem));
            if (domain == null)
                throw new ArgumentNullException(nameof(domain));

            int nelx = problem.NelX;
            int nely = problem.NelY;

            if (nelx < 1 || nely < 1)
                throw new ArgumentException("NelX and NelY must both be >= 1.");

            double volfrac = problem.VolumeFraction;
            double penal = problem.Penal;
            double rmin = problem.FilterRadius;
            int maxIter = problem.MaxIterations;

            double E0 = problem.YoungsModulusSolid;
            double Emin = problem.YoungsModulusMin;
            double nu = problem.PoissonRatio;

            int nnx = nelx + 1;
            int nny = nely + 1;
            int nodeCount = nnx * nny;
            int dofCount = nodeCount * 2;

            ValidateInputs(problem, domain, dofCount);

            bool[,] mask = domain.DesignMask;
            double[] forces = domain.Forces;
            bool[] isFixed = domain.FixedDofs;

            int activeElementCount = CountActiveElements(mask, nelx, nely);
            if (activeElementCount == 0)
                throw new InvalidOperationException("Design mask contains no active elements.");

            bool[] activeNodes = BuildActiveNodes(mask, nelx, nely, nny);

            if (!HasLoadOnActiveNode(forces, activeNodes))
                throw new InvalidOperationException("No load is applied to any active node.");

            bool hasAnyFixedOnActiveNode = HasFixedDofOnActiveNode(isFixed, activeNodes);
            if (!hasAnyFixedOnActiveNode)
                throw new InvalidOperationException("No fixed DOFs were found on active nodes.");

            List<int> freeDofs = BuildFreeDofs(isFixed, activeNodes);
            if (freeDofs.Count == 0)
                throw new InvalidOperationException("No free DOFs found on active connected nodes.");

            double[,] x = new double[nelx, nely];
            double[,] xNew = new double[nelx, nely];
            double[,] dc = new double[nelx, nely];
            double[,] dcFiltered = new double[nelx, nely];

            InitializeDensity(x, mask, nelx, nely, volfrac);

            double[,] KE = BuildElementStiffness(nu);

            double compliance = 0.0;

            for (int iter = 0; iter < maxIter; iter++)
            {
                double[] U = SolveDisplacements(
                    nelx,
                    nely,
                    nny,
                    x,
                    mask,
                    penal,
                    E0,
                    Emin,
                    KE,
                    forces,
                    freeDofs);

                compliance = 0.0;
                ClearArray(dc, nelx, nely);
                ClearArray(dcFiltered, nelx, nely);

                for (int ex = 0; ex < nelx; ex++)
                {
                    for (int ey = 0; ey < nely; ey++)
                    {
                        if (!mask[ex, ey])
                            continue;

                        int[] edofs = ElementDofs(ex, ey, nny);

                        double[] ue = new double[8];
                        for (int i = 0; i < 8; i++)
                            ue[i] = U[edofs[i]];

                        double ce = Dot(ue, Multiply(KE, ue));

                        double rho = x[ex, ey];
                        double Ee = Emin + Math.Pow(rho, penal) * (E0 - Emin);

                        compliance += Ee * ce;
                        dc[ex, ey] = -penal * Math.Pow(rho, penal - 1.0) * (E0 - Emin) * ce;
                    }
                }

                ApplySensitivityFilter(nelx, nely, rmin, mask, x, dc, dcFiltered);
                UpdateDensitiesOC(nelx, nely, volfrac, mask, x, xNew, dcFiltered);

                double change = 0.0;
                for (int ex = 0; ex < nelx; ex++)
                {
                    for (int ey = 0; ey < nely; ey++)
                    {
                        if (!mask[ex, ey])
                        {
                            x[ex, ey] = 0.0;
                            xNew[ex, ey] = 0.0;
                            continue;
                        }

                        double diff = Math.Abs(xNew[ex, ey] - x[ex, ey]);
                        if (diff > change)
                            change = diff;

                        x[ex, ey] = xNew[ex, ey];
                    }
                }

                if (change < 0.01 && iter > 10)
                {
                    return new TopOptResult2D
                    {
                        NelX = nelx,
                        NelY = nely,
                        Density = x,
                        Compliance = compliance,
                        Iterations = iter + 1
                    };
                }
            }

            return new TopOptResult2D
            {
                NelX = nelx,
                NelY = nely,
                Density = x,
                Compliance = compliance,
                Iterations = maxIter
            };
        }

        private static void ValidateInputs(TopOptProblem2D problem, TopOptDomain2D domain, int dofCount)
        {
            if (domain.DesignMask == null)
                throw new ArgumentException("DesignMask cannot be null.");
            if (domain.FixedDofs == null)
                throw new ArgumentException("FixedDofs cannot be null.");
            if (domain.Forces == null)
                throw new ArgumentException("Forces cannot be null.");

            if (domain.DesignMask.GetLength(0) != problem.NelX || domain.DesignMask.GetLength(1) != problem.NelY)
                throw new ArgumentException("DesignMask size must match NelX x NelY.");

            if (domain.FixedDofs.Length != dofCount)
                throw new ArgumentException("FixedDofs length must match total DOF count.");

            if (domain.Forces.Length != dofCount)
                throw new ArgumentException("Forces length must match total DOF count.");
        }

        private static int CountActiveElements(bool[,] mask, int nelx, int nely)
        {
            int count = 0;
            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    if (mask[ex, ey])
                        count++;
                }
            }
            return count;
        }

        private static void InitializeDensity(double[,] x, bool[,] mask, int nelx, int nely, double volfrac)
        {
            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    x[ex, ey] = mask[ex, ey] ? volfrac : 0.0;
                }
            }
        }

        private static void ClearArray(double[,] a, int nelx, int nely)
        {
            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    a[ex, ey] = 0.0;
                }
            }
        }

        private static bool[] BuildActiveNodes(bool[,] mask, int nelx, int nely, int nny)
        {
            int nodeCount = (nelx + 1) * (nely + 1);
            bool[] activeNodes = new bool[nodeCount];

            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    if (!mask[ex, ey])
                        continue;

                    int n1 = NodeIndex(ex, ey, nny);
                    int n2 = NodeIndex(ex + 1, ey, nny);
                    int n3 = NodeIndex(ex + 1, ey + 1, nny);
                    int n4 = NodeIndex(ex, ey + 1, nny);

                    activeNodes[n1] = true;
                    activeNodes[n2] = true;
                    activeNodes[n3] = true;
                    activeNodes[n4] = true;
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

                int dofX = 2 * node;
                int dofY = 2 * node + 1;

                if (Math.Abs(forces[dofX]) > 1e-12 || Math.Abs(forces[dofY]) > 1e-12)
                    return true;
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

                int dofX = 2 * node;
                int dofY = 2 * node + 1;

                if (isFixed[dofX] || isFixed[dofY])
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

                int dofX = 2 * node;
                int dofY = 2 * node + 1;

                if (!isFixed[dofX])
                    freeDofs.Add(dofX);

                if (!isFixed[dofY])
                    freeDofs.Add(dofY);
            }

            return freeDofs;
        }

        private double[] SolveDisplacements(
            int nelx,
            int nely,
            int nny,
            double[,] x,
            bool[,] mask,
            double penal,
            double E0,
            double Emin,
            double[,] KE,
            double[] F,
            List<int> freeDofs)
        {
            int nodeCount = (nelx + 1) * (nely + 1);
            int dofCount = nodeCount * 2;

            double[,] K = new double[dofCount, dofCount];

            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    if (!mask[ex, ey])
                        continue;

                    double rho = x[ex, ey];
                    double Ee = Emin + Math.Pow(rho, penal) * (E0 - Emin);

                    int[] edofs = ElementDofs(ex, ey, nny);

                    for (int i = 0; i < 8; i++)
                    {
                        int ii = edofs[i];
                        for (int j = 0; j < 8; j++)
                        {
                            int jj = edofs[j];
                            K[ii, jj] += Ee * KE[i, j];
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
            {
                U[freeDofs[i]] = Uf[i];
            }

            return U;
        }

        private static int NodeIndex(int ix, int iy, int nny)
        {
            return ix * nny + iy;
        }

        private static int[] ElementDofs(int ex, int ey, int nny)
        {
            int n1 = NodeIndex(ex, ey, nny);
            int n2 = NodeIndex(ex + 1, ey, nny);
            int n3 = NodeIndex(ex + 1, ey + 1, nny);
            int n4 = NodeIndex(ex, ey + 1, nny);

            return new int[]
            {
                2 * n1, 2 * n1 + 1,
                2 * n2, 2 * n2 + 1,
                2 * n3, 2 * n3 + 1,
                2 * n4, 2 * n4 + 1
            };
        }

        private static double[,] BuildElementStiffness(double nu)
        {
            double[] k = new double[]
            {
                0.5 - nu / 6.0,
                0.125 + nu / 8.0,
                -0.25 - nu / 12.0,
                -0.125 + 3.0 * nu / 8.0,
                -0.25 + nu / 12.0,
                -0.125 - nu / 8.0,
                nu / 6.0,
                0.125 - 3.0 * nu / 8.0
            };

            double factor = 1.0 / (1.0 - nu * nu);

            double[,] KE = new double[8, 8]
            {
                { k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7] },
                { k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2] },
                { k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1] },
                { k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4] },
                { k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3] },
                { k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6] },
                { k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5] },
                { k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0] }
            };

            for (int i = 0; i < 8; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    KE[i, j] *= factor;
                }
            }

            return KE;
        }

        private static void ApplySensitivityFilter(
            int nelx,
            int nely,
            double rmin,
            bool[,] mask,
            double[,] x,
            double[,] dc,
            double[,] dcFiltered)
        {
            int r = (int)Math.Floor(rmin);

            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    if (!mask[ex, ey])
                    {
                        dcFiltered[ex, ey] = 0.0;
                        continue;
                    }

                    double sum = 0.0;
                    double val = 0.0;

                    int iMin = Math.Max(ex - r, 0);
                    int iMax = Math.Min(ex + r, nelx - 1);
                    int jMin = Math.Max(ey - r, 0);
                    int jMax = Math.Min(ey + r, nely - 1);

                    for (int i = iMin; i <= iMax; i++)
                    {
                        for (int j = jMin; j <= jMax; j++)
                        {
                            if (!mask[i, j])
                                continue;

                            double dist = Math.Sqrt((ex - i) * (ex - i) + (ey - j) * (ey - j));
                            double weight = rmin - dist;

                            if (weight > 0.0)
                            {
                                sum += weight;
                                val += weight * x[i, j] * dc[i, j];
                            }
                        }
                    }

                    double denom = x[ex, ey] * sum;
                    if (denom > 1e-12)
                        dcFiltered[ex, ey] = val / denom;
                    else
                        dcFiltered[ex, ey] = dc[ex, ey];
                }
            }
        }

        private static void UpdateDensitiesOC(
            int nelx,
            int nely,
            double volfrac,
            bool[,] mask,
            double[,] x,
            double[,] xNew,
            double[,] dc)
        {
            int activeCount = CountActiveElements(mask, nelx, nely);
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
                        if (!mask[ex, ey])
                        {
                            xNew[ex, ey] = 0.0;
                            continue;
                        }

                        double rho = x[ex, ey];
                        double b = -dc[ex, ey] / lmid;
                        if (b < 1e-12)
                            b = 1e-12;

                        double candidate = rho * Math.Sqrt(b);

                        double lower = Math.Max(xmin, rho - move);
                        double upper = Math.Min(1.0, rho + move);

                        if (candidate < lower) candidate = lower;
                        if (candidate > upper) candidate = upper;

                        xNew[ex, ey] = candidate;
                        total += candidate;
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
    }
}