using Toplet_v0_Alpha.TopOpt2D;
using Eto.Forms;
using Rhino;
using Rhino.Commands;
using Rhino.DocObjects;
using Rhino.Geometry;
using Rhino.Input;
using Rhino.Input.Custom;
using System;
using System.Drawing;
using System.Windows.Forms;

namespace Toplet_v0_Alpha
{
    public class Toplet2DCommand : Rhino.Commands.Command
    {
        public Toplet2DCommand()
        {
            Instance = this;
        }

        public static Toplet2DCommand Instance { get; private set; }

        public override string EnglishName => "Toplet2D";

        protected override Result RunCommand(RhinoDoc doc, RunMode mode)
        {
            // 0. Ask user for a closed curve
            ObjRef curveRef;
            var rc = RhinoGet.GetOneObject(
                "Select closed curve for design domain",
                false,
                ObjectType.Curve,
                out curveRef);

            if (rc != Result.Success || curveRef == null)
                return rc;

            Curve curve = curveRef.Curve();
            if (curve == null || !curve.IsClosed)
            {
                RhinoApp.WriteLine("Curve must be closed.");
                return Result.Failure;
            }

            // 1. Get curve bounding box
            BoundingBox bbox = curve.GetBoundingBox(true);

            double minX = bbox.Min.X;
            double minY = bbox.Min.Y;
            double maxX = bbox.Max.X;
            double maxY = bbox.Max.Y;

            double width = maxX - minX;
            double height = maxY - minY;

            if (width <= 0.0 || height <= 0.0)
            {
                RhinoApp.WriteLine("Curve bounding box is invalid.");
                return Result.Failure;
            }

            // 2. Define square cell size in model units
            double cellSize = 1.0;

            // Safety cap for this dense toy solver
            int maxElementsX = 60;
            int maxElementsY = 40;

            // If the curve is too large, automatically increase cell size
            if (width / cellSize > maxElementsX || height / cellSize > maxElementsY)
            {
                double cellSizeX = width / maxElementsX;
                double cellSizeY = height / maxElementsY;
                cellSize = Math.Max(cellSizeX, cellSizeY);
            }

            int nelx = Math.Max(1, (int)Math.Ceiling(width / cellSize));
            int nely = Math.Max(1, (int)Math.Ceiling(height / cellSize));

            // Square cells
            double dx = cellSize;
            double dy = cellSize;

            // 3. Build optimization problem using the computed grid size
            var problem = new TopOptProblem2D
            {
                NelX = nelx,
                NelY = nely,
                VolumeFraction = 0.5,
                Penal = 3.0,
                FilterRadius = 1.5,
                MaxIterations = 30
            };

            int nnx = nelx + 1;
            int nny = nely + 1;
            int nodeCount = nnx * nny;
            int dofCount = nodeCount * 2;

            if (dofCount > 12000)
            {
                RhinoApp.WriteLine("Grid is too large for the current dense solver.");
                RhinoApp.WriteLine("Try a larger cell size or a smaller curve.");
                RhinoApp.WriteLine("Current DOF count: {0}", dofCount);
                return Result.Failure;
            }
            // 4. Build raw design mask from curve containment
            bool[,] rawMask = new bool[nelx, nely];

            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    double cx = minX + (ex + 0.5) * dx;
                    double cy = minY + (ey + 0.5) * dy;

                    Point3d pt = new Point3d(cx, cy, 0.0);
                    PointContainment containment = curve.Contains(pt, Plane.WorldXY, doc.ModelAbsoluteTolerance);

                    rawMask[ex, ey] =
                        containment == PointContainment.Inside ||
                        containment == PointContainment.Coincident;
                }
            }

            // Keep only the region that is edge-connected to the leftmost active support side
            bool[,] designMask = ExtractSupportConnectedMask(rawMask, nelx, nely);

            // 5. Build supports and load from the active mask
            double[] forces = new double[dofCount];
            bool[] fixedDofs = new bool[dofCount];

            // Find leftmost active column
            int leftActiveColumn = -1;
            for (int ex = 0; ex < nelx; ex++)
            {
                bool hasActive = false;
                for (int ey = 0; ey < nely; ey++)
                {
                    if (designMask[ex, ey])
                    {
                        hasActive = true;
                        break;
                    }
                }

                if (hasActive)
                {
                    leftActiveColumn = ex;
                    break;
                }
            }

            if (leftActiveColumn < 0)
            {
                RhinoApp.WriteLine("No active cells found inside the selected curve.");
                return Result.Failure;
            }

            // Fix all nodes on the left boundary of the leftmost active column
            for (int ey = 0; ey <= nely; ey++)
            {
                int node = leftActiveColumn * nny + ey;
                fixedDofs[2 * node] = true;
                fixedDofs[2 * node + 1] = true;
            }

            // Find rightmost active column
            int rightActiveColumn = -1;
            for (int ex = nelx - 1; ex >= 0; ex--)
            {
                bool hasActive = false;
                for (int ey = 0; ey < nely; ey++)
                {
                    if (designMask[ex, ey])
                    {
                        hasActive = true;
                        break;
                    }
                }

                if (hasActive)
                {
                    rightActiveColumn = ex;
                    break;
                }
            }

            if (rightActiveColumn < 0)
            {
                RhinoApp.WriteLine("Could not find a valid load column.");
                return Result.Failure;
            }

            // Find row nearest the middle in the rightmost active column
            int targetRow = -1;
            int midRow = nely / 2;
            int bestDist = int.MaxValue;

            for (int ey = 0; ey < nely; ey++)
            {
                if (!designMask[rightActiveColumn, ey])
                    continue;

                int dist = Math.Abs(ey - midRow);
                if (dist < bestDist)
                {
                    bestDist = dist;
                    targetRow = ey;
                }
            }

            if (targetRow < 0)
            {
                RhinoApp.WriteLine("Could not find a valid load location.");
                return Result.Failure;
            }

            // Apply downward load at the right-side node of that active cell
            int loadNode = (rightActiveColumn + 1) * nny + targetRow;
            int loadDofY = 2 * loadNode + 1;
            forces[loadDofY] = -1.0;

            // 6. Debug info
            int rawActiveCount = 0;
            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    if (rawMask[ex, ey])
                        rawActiveCount++;
                }
            }

            int activeCount = 0;
            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    if (designMask[ex, ey])
                        activeCount++;
                }
            }

            RhinoApp.WriteLine("Raw active cells: {0}", rawActiveCount);
            RhinoApp.WriteLine("Connected active cells: {0}", activeCount);
            RhinoApp.WriteLine("Grid: {0} x {1}", nelx, nely);
            RhinoApp.WriteLine("Cell size used: {0}", cellSize);
            RhinoApp.WriteLine("Support column: {0}", leftActiveColumn);
            RhinoApp.WriteLine("Load column: {0}, load row: {1}", rightActiveColumn, targetRow);

            // 7. Solve
            TopOptDomain2D domain = new TopOptDomain2D
            {
                DesignMask = designMask,
                Forces = forces,
                FixedDofs = fixedDofs
            };

            TopOptSolver2D solver = new TopOptSolver2D();
            TopOptResult2D result = solver.Solve(problem, domain);

            if (result == null || result.Density == null)
            {
                RhinoApp.WriteLine("Solver returned no result.");
                return Result.Failure;
            }

            // 8. Draw result
            ShowDensityBitmap(result, 4, 0.3);

            RhinoApp.WriteLine("TopOpt completed.");
            RhinoApp.WriteLine("Iterations: {0}", result.Iterations);
            RhinoApp.WriteLine("Final compliance: {0:F4}", result.Compliance);

            doc.Views.Redraw();
            return Result.Success;
        }
        private bool[,] ExtractSupportConnectedMask(bool[,] rawMask, int nelx, int nely)
        {
            bool[,] connected = new bool[nelx, nely];

            // Find leftmost active column
            int leftActiveColumn = -1;
            for (int ex = 0; ex < nelx; ex++)
            {
                bool hasActive = false;
                for (int ey = 0; ey < nely; ey++)
                {
                    if (rawMask[ex, ey])
                    {
                        hasActive = true;
                        break;
                    }
                }

                if (hasActive)
                {
                    leftActiveColumn = ex;
                    break;
                }
            }

            if (leftActiveColumn < 0)
                return connected;

            // Flood fill from every active cell in the leftmost active column
            System.Collections.Generic.Queue<(int ex, int ey)> queue =
                new System.Collections.Generic.Queue<(int ex, int ey)>();

            for (int ey = 0; ey < nely; ey++)
            {
                if (rawMask[leftActiveColumn, ey])
                {
                    connected[leftActiveColumn, ey] = true;
                    queue.Enqueue((leftActiveColumn, ey));
                }
            }

            int[] nx = new int[] { -1, 1, 0, 0 };
            int[] ny = new int[] { 0, 0, -1, 1 };

            while (queue.Count > 0)
            {
                var cell = queue.Dequeue();

                for (int k = 0; k < 4; k++)
                {
                    int ex2 = cell.ex + nx[k];
                    int ey2 = cell.ey + ny[k];

                    if (ex2 < 0 || ex2 >= nelx || ey2 < 0 || ey2 >= nely)
                        continue;

                    if (!rawMask[ex2, ey2])
                        continue;

                    if (connected[ex2, ey2])
                        continue;

                    connected[ex2, ey2] = true;
                    queue.Enqueue((ex2, ey2));
                }
            }

            return connected;
        }
        private void ShowDensityBitmap(TopOptResult2D result, int scale, double threshold)
        {
            int width = result.NelX;
            int height = result.NelY;

            int imgW = width * scale;
            int imgH = height * scale;

            Bitmap bmp = new Bitmap(imgW, imgH);

            for (int ex = 0; ex < width; ex++)
            {
                for (int ey = 0; ey < height; ey++)
                {
                    double rho = result.Density[ex, ey];
                    Color color = (rho >= threshold) ? Color.Black : Color.White;

                    int imgY = height - 1 - ey;

                    for (int sx = 0; sx < scale; sx++)
                    {
                        for (int sy = 0; sy < scale; sy++)
                        {
                            int px = ex * scale + sx;
                            int py = imgY * scale + sy;
                            bmp.SetPixel(px, py, color);
                        }
                    }
                }
            }

            System.Windows.Forms.Form form = new System.Windows.Forms.Form();
            form.Text = "Topology Result";
            form.ClientSize = new Size(imgW, imgH);

            System.Windows.Forms.PictureBox pb = new System.Windows.Forms.PictureBox();
            pb.Dock = System.Windows.Forms.DockStyle.Fill;
            pb.Image = bmp;
            pb.SizeMode = System.Windows.Forms.PictureBoxSizeMode.Zoom;

            form.Controls.Add(pb);
            form.Show();
        }
    }
}