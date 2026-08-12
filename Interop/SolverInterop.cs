using System;
using System.Diagnostics;
using System.IO;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Threading.Tasks;
using System.Windows.Forms;
using Toplet_v0_Alpha.TopOpt2D;
using Toplet_v0_Alpha.TopOpt3D;

namespace Toplet_v0_Alpha.Interop
{
    internal static class NativeMethods
    {
        private const string DllName = "TopletSolverNative";

        // Pre-load the DLL using its full path before P/Invoke tries to resolve
        // it by name. This is necessary inside Rhino because the plugin directory
        // is not automatically on the Win32 DLL search path.
        static NativeMethods()
        {
            string dir = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location);
            string dllPath = Path.Combine(dir, DllName + ".dll");
            if (File.Exists(dllPath))
                LoadLibraryW(dllPath);
        }

        [DllImport("kernel32.dll", CharSet = CharSet.Unicode)]
        private static extern IntPtr LoadLibraryW(string lpFileName);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate void ProgressCallback(int iter, int maxIter, double compliance);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        internal static extern int solve_3d(
            [In]  byte[]   mask_flat,
            [In]  double[] forces,
            [In]  byte[]   fixed_dofs,
            int nelx, int nely, int nelz,
            double vol_frac, double penal, double filter_rad, int max_iter,
            double E0, double Emin, double nu,
            [Out] double[] density_out,
            out   double   compliance,
            out   int      iterations,
            ProgressCallback progressCb);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        internal static extern int solve_2d(
            [In]  byte[]   mask_flat,
            [In]  double[] forces,
            [In]  byte[]   fixed_dofs,
            int nelx, int nely,
            double vol_frac, double penal, double filter_rad, int max_iter,
            double E0, double Emin, double nu,
            [Out] double[] density_out,
            out   double   compliance,
            out   int      iterations,
            ProgressCallback progressCb);
    }

    public static class NativeSolver3D
    {
        public static TopOptResult3D Solve(TopOptProblem3D problem, TopOptDomain3D domain)
        {
            int nelx = problem.NelX, nely = problem.NelY, nelz = problem.NelZ;
            int elemCount = nelx * nely * nelz;

            byte[] maskFlat = new byte[elemCount];
            bool[,,] mask = domain.DesignMask;
            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++)
            for (int ez = 0; ez < nelz; ez++)
                maskFlat[ex*(nely*nelz) + ey*nelz + ez] = mask[ex, ey, ez] ? (byte)1 : (byte)0;

            bool[] fixedBool = domain.FixedDofs;
            byte[] fixedFlat = new byte[fixedBool.Length];
            for (int i = 0; i < fixedBool.Length; i++)
                fixedFlat[i] = fixedBool[i] ? (byte)1 : (byte)0;

            double[] densityFlat = new double[elemCount];
            double compliance = 0;
            int iterations = 0;
            int rc = 0;
            Exception thrownEx = null;
            var sw = new Stopwatch();
            bool wasCancelled3d = false;

            using (var progressForm = new SolverProgressForm(problem.MaxIterations))
            {
                NativeMethods.ProgressCallback callback = (iter, maxIter, comp) =>
                {
                    if (progressForm.IsHandleCreated)
                        progressForm.BeginInvoke(new Action(() =>
                            progressForm.UpdateProgress(iter, maxIter, comp)));
                };

                progressForm.Shown += (s, e) =>
                {
                    sw.Start();
                    Task.Run(() =>
                    {
                        try
                        {
                            rc = NativeMethods.solve_3d(
                                maskFlat, domain.Forces, fixedFlat,
                                nelx, nely, nelz,
                                problem.VolumeFraction, problem.Penal, problem.FilterRadius, problem.MaxIterations,
                                problem.YoungsModulusSolid, problem.YoungsModulusMin, problem.PoissonRatio,
                                densityFlat, out compliance, out iterations, callback);
                        }
                        catch (Exception ex) { thrownEx = ex; }
                        finally
                        {
                            sw.Stop();
                            if (progressForm.IsHandleCreated)
                                progressForm.BeginInvoke(new Action(() => progressForm.Close()));
                        }
                    });
                };

                progressForm.ShowDialog();
                GC.KeepAlive(callback);
                wasCancelled3d = progressForm.WasCancelled;
            }

            if (wasCancelled3d) return null;
            if (thrownEx != null) throw thrownEx;
            if (rc != 0) throw new InvalidOperationException($"Native solve_3d returned error code {rc}.");

            bool converged = iterations < problem.MaxIterations;
            using (var doneForm = new SolverCompletedForm(iterations, problem.MaxIterations, compliance, sw.Elapsed, converged))
                doneForm.ShowDialog();

            double[,,] density = new double[nelx, nely, nelz];
            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++)
            for (int ez = 0; ez < nelz; ez++)
                density[ex, ey, ez] = densityFlat[ex*(nely*nelz) + ey*nelz + ez];

            return new TopOptResult3D {
                NelX = nelx, NelY = nely, NelZ = nelz,
                Density = density, Compliance = compliance, Iterations = iterations
            };
        }
    }

    public static class NativeSolver2D
    {
        public static TopOptResult2D Solve(TopOptProblem2D problem, TopOptDomain2D domain)
        {
            int nelx = problem.NelX, nely = problem.NelY;
            int elemCount = nelx * nely;

            byte[] maskFlat = new byte[elemCount];
            bool[,] mask = domain.DesignMask;
            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++)
                maskFlat[ex*nely + ey] = mask[ex, ey] ? (byte)1 : (byte)0;

            bool[] fixedBool = domain.FixedDofs;
            byte[] fixedFlat = new byte[fixedBool.Length];
            for (int i = 0; i < fixedBool.Length; i++)
                fixedFlat[i] = fixedBool[i] ? (byte)1 : (byte)0;

            double[] densityFlat = new double[elemCount];
            double compliance = 0;
            int iterations = 0;
            int rc = 0;
            Exception thrownEx = null;
            var sw = new Stopwatch();
            bool wasCancelled2d = false;

            using (var progressForm = new SolverProgressForm(problem.MaxIterations))
            {
                NativeMethods.ProgressCallback callback = (iter, maxIter, comp) =>
                {
                    if (progressForm.IsHandleCreated)
                        progressForm.BeginInvoke(new Action(() =>
                            progressForm.UpdateProgress(iter, maxIter, comp)));
                };

                progressForm.Shown += (s, e) =>
                {
                    sw.Start();
                    Task.Run(() =>
                    {
                        try
                        {
                            rc = NativeMethods.solve_2d(
                                maskFlat, domain.Forces, fixedFlat,
                                nelx, nely,
                                problem.VolumeFraction, problem.Penal, problem.FilterRadius, problem.MaxIterations,
                                problem.YoungsModulusSolid, problem.YoungsModulusMin, problem.PoissonRatio,
                                densityFlat, out compliance, out iterations, callback);
                        }
                        catch (Exception ex) { thrownEx = ex; }
                        finally
                        {
                            sw.Stop();
                            if (progressForm.IsHandleCreated)
                                progressForm.BeginInvoke(new Action(() => progressForm.Close()));
                        }
                    });
                };

                progressForm.ShowDialog();
                GC.KeepAlive(callback);
                wasCancelled2d = progressForm.WasCancelled;
            }

            if (wasCancelled2d) return null;
            if (thrownEx != null) throw thrownEx;
            if (rc != 0) throw new InvalidOperationException($"Native solve_2d returned error code {rc}.");

            bool converged = iterations < problem.MaxIterations;
            using (var doneForm = new SolverCompletedForm(iterations, problem.MaxIterations, compliance, sw.Elapsed, converged))
                doneForm.ShowDialog();

            double[,] density = new double[nelx, nely];
            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++)
                density[ex, ey] = densityFlat[ex*nely + ey];

            return new TopOptResult2D {
                NelX = nelx, NelY = nely,
                Density = density, Compliance = compliance, Iterations = iterations
            };
        }
    }
}
