using System;
using System.Diagnostics;
using System.IO;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using Rhino;
using Toplet.TopOpt2D;
using Toplet.TopOpt3D;

namespace Toplet.Interop
{
    internal static class NativeMethods
    {
        private const string DllName = "TopletSolverNative";

        // Must match toplet_solver.h's recommended diag_buf_capacity.
        internal const int DiagBufCapacity = 8192;

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

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
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
            ProgressCallback progressCb,
            [Out] StringBuilder diagBufOut,
            int diagBufCapacity);

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

    // Shared by NativeSolver2D/NativeSolver3D: both drive a native solve call
    // from a background task while a modal progress form owns the UI thread,
    // then surface cancellation/exceptions/error codes the same way. Only the
    // native call itself (and the mask flattening around it) differs by
    // dimensionality, so that part stays in each Solve() method below.
    internal readonly struct SolveOutcome
    {
        public readonly bool WasCancelled;
        public readonly int Rc;
        public readonly double Compliance;
        public readonly int Iterations;
        public readonly TimeSpan Elapsed;

        public SolveOutcome(bool wasCancelled, int rc, double compliance, int iterations, TimeSpan elapsed)
        {
            WasCancelled = wasCancelled;
            Rc = rc;
            Compliance = compliance;
            Iterations = iterations;
            Elapsed = elapsed;
        }
    }

    internal static class NativeSolverRunner
    {
        public static SolveOutcome Run(int maxIterations, Func<NativeMethods.ProgressCallback, (int rc, double compliance, int iterations)> nativeCall)
        {
            int rc = 0;
            double compliance = 0;
            int iterations = 0;
            Exception thrownEx = null;
            var sw = new Stopwatch();
            bool wasCancelled;

            using (var progressForm = new SolverProgressForm(maxIterations))
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
                            (rc, compliance, iterations) = nativeCall(callback);
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
                wasCancelled = progressForm.WasCancelled;
            }

            if (thrownEx != null) throw thrownEx;
            return new SolveOutcome(wasCancelled, rc, compliance, iterations, sw.Elapsed);
        }
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
            var diagBuf = new StringBuilder(NativeMethods.DiagBufCapacity);

            SolveOutcome outcome = NativeSolverRunner.Run(problem.MaxIterations, callback =>
            {
                int rc = NativeMethods.solve_3d(
                    maskFlat, domain.Forces, fixedFlat,
                    nelx, nely, nelz,
                    problem.VolumeFraction, problem.Penal, problem.FilterRadius, problem.MaxIterations,
                    problem.YoungsModulusSolid, problem.YoungsModulusMin, problem.PoissonRatio,
                    densityFlat, out double compliance, out int iterations, callback,
                    diagBuf, diagBuf.Capacity);
                return (rc, compliance, iterations);
            });

            // GMG solver diagnostics: printed to Rhino's command line rather than
            // a blocking dialog, so a solve never stalls waiting on user input.
            if (diagBuf.Length > 0)
                RhinoApp.WriteLine(diagBuf.ToString());

            if (outcome.WasCancelled) return null;
            if (outcome.Rc != 0) throw new InvalidOperationException($"Native solve_3d returned error code {outcome.Rc}.");

            bool converged = outcome.Iterations < problem.MaxIterations;
            using (var doneForm = new SolverCompletedForm(outcome.Iterations, problem.MaxIterations, outcome.Compliance, outcome.Elapsed, converged))
                doneForm.ShowDialog();

            double[,,] density = new double[nelx, nely, nelz];
            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++)
            for (int ez = 0; ez < nelz; ez++)
                density[ex, ey, ez] = densityFlat[ex*(nely*nelz) + ey*nelz + ez];

            return new TopOptResult3D {
                NelX = nelx, NelY = nely, NelZ = nelz,
                Density = density, Compliance = outcome.Compliance, Iterations = outcome.Iterations
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

            SolveOutcome outcome = NativeSolverRunner.Run(problem.MaxIterations, callback =>
            {
                int rc = NativeMethods.solve_2d(
                    maskFlat, domain.Forces, fixedFlat,
                    nelx, nely,
                    problem.VolumeFraction, problem.Penal, problem.FilterRadius, problem.MaxIterations,
                    problem.YoungsModulusSolid, problem.YoungsModulusMin, problem.PoissonRatio,
                    densityFlat, out double compliance, out int iterations, callback);
                return (rc, compliance, iterations);
            });

            if (outcome.WasCancelled) return null;
            if (outcome.Rc != 0) throw new InvalidOperationException($"Native solve_2d returned error code {outcome.Rc}.");

            bool converged = outcome.Iterations < problem.MaxIterations;
            using (var doneForm = new SolverCompletedForm(outcome.Iterations, problem.MaxIterations, outcome.Compliance, outcome.Elapsed, converged))
                doneForm.ShowDialog();

            double[,] density = new double[nelx, nely];
            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++)
                density[ex, ey] = densityFlat[ex*nely + ey];

            return new TopOptResult2D {
                NelX = nelx, NelY = nely,
                Density = density, Compliance = outcome.Compliance, Iterations = outcome.Iterations
            };
        }
    }
}
