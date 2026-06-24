using Toplet_v0_Alpha.TopOpt3D;
using Toplet_v0_Alpha.Interop;
using Rhino;
using Rhino.Commands;
using Rhino.DocObjects;
using Rhino.Geometry;
using Rhino.Input;
using Rhino.Input.Custom;
using System;
using System.Collections.Generic;

namespace Toplet_v0_Alpha
{
    public class Toplet3DCommand : Command
    {
        private enum LoadInputMode
        {
            Point,
            Face
        }

        private enum SupportInputMode
        {
            Points,
            Faces
        }

        private TopletDisplayConduit _activeConduit;
        private ResultViewerForm     _activeViewer;

        public Toplet3DCommand()
        {
            Instance = this;
        }

        public static Toplet3DCommand Instance { get; private set; }

        public override string EnglishName => "Toplet3D";

        protected override Result RunCommand(RhinoDoc doc, RunMode mode)
        {
            // Close any viewer left open from a previous run
            if (_activeViewer != null && !_activeViewer.IsDisposed) { _activeViewer.Close(); }
            if (_activeConduit != null) { _activeConduit.Enabled = false; _activeConduit = null; }

            ObjRef brepRef;
            Result rc = RhinoGet.GetOneObject(
                "Select closed solid for 3D design domain",
                false,
                ObjectType.Brep,
                out brepRef);

            if (rc != Result.Success || brepRef == null)
                return rc;

            Brep brep = brepRef.Brep();
            if (brep == null || !brep.IsSolid)
            {
                RhinoApp.WriteLine("Selected Brep must be a closed solid.");
                return Result.Failure;
            }

            BoundingBox bbox = brep.GetBoundingBox(true);
            if (!bbox.IsValid)
            {
                RhinoApp.WriteLine("Solid bounding box is invalid.");
                return Result.Failure;
            }

            double width  = bbox.Max.X - bbox.Min.X;
            double depth  = bbox.Max.Y - bbox.Min.Y;
            double height = bbox.Max.Z - bbox.Min.Z;

            if (width <= 0.0 || depth <= 0.0 || height <= 0.0)
            {
                RhinoApp.WriteLine("Solid bounding box is invalid.");
                return Result.Failure;
            }

            double longestDim      = Math.Max(width, Math.Max(depth, height));
            double defaultCellSize = Math.Max(longestDim / 30.0, doc.ModelAbsoluteTolerance * 10.0);

            var setupForm = new TopletSetupForm(width, depth, height, defaultCellSize);
            if (setupForm.ShowDialog() != System.Windows.Forms.DialogResult.OK)
                return Result.Cancel;

            double cellSize = setupForm.CellSize;
            LoadInputMode    loadMode    = setupForm.LoadIsFace    ? LoadInputMode.Face       : LoadInputMode.Point;
            SupportInputMode supportMode = setupForm.SupportIsFace ? SupportInputMode.Faces   : SupportInputMode.Points;

            // Round up to even so the brep centroid lands exactly on a cell boundary,
            // giving identical left/right (and front/back) half-grids.
            int nelx = Math.Max(2, (int)Math.Ceiling(width  / cellSize)); if (nelx % 2 != 0) nelx++;
            int nely = Math.Max(2, (int)Math.Ceiling(depth  / cellSize)); if (nely % 2 != 0) nely++;
            int nelz = Math.Max(2, (int)Math.Ceiling(height / cellSize)); if (nelz % 2 != 0) nelz++;

            double dx = cellSize;
            double dy = cellSize;
            double dz = cellSize;

            // Centre the grid on the brep centroid.
            double centX = (bbox.Min.X + bbox.Max.X) * 0.5;
            double centY = (bbox.Min.Y + bbox.Max.Y) * 0.5;
            double centZ = (bbox.Min.Z + bbox.Max.Z) * 0.5;
            double startX = centX - (nelx * 0.5) * cellSize;
            double startY = centY - (nely * 0.5) * cellSize;
            double startZ = centZ - (nelz * 0.5) * cellSize;

            BoundingBox gridBox = new BoundingBox(
                new Point3d(startX, startY, startZ),
                new Point3d(startX + nelx * cellSize, startY + nely * cellSize, startZ + nelz * cellSize));

            TopOptProblem3D problem = new TopOptProblem3D
            {
                NelX = nelx,
                NelY = nely,
                NelZ = nelz,
                VolumeFraction = 0.5,
                Penal = 3.0,
                FilterRadius = 1.5,
                MaxIterations = 15
            };

            int nnx = nelx + 1;
            int nny = nely + 1;
            int nnz = nelz + 1;
            int nodeCount = nnx * nny * nnz;
            int dofCount = nodeCount * 3;

            double tol = Math.Max(doc.ModelAbsoluteTolerance, cellSize * 0.35);

            bool[,,] rawMask = VoxelDomainBuilder3D.BuildMaskFromBrep(
                brep,
                gridBox,
                nelx,
                nely,
                nelz,
                dx,
                dy,
                dz,
                tol);

            int rawActiveCount = VoxelDomainBuilder3D.CountActiveElements(rawMask, nelx, nely, nelz);
            if (rawActiveCount == 0)
            {
                RhinoApp.WriteLine("No active voxels found inside the selected Brep.");
                return Result.Failure;
            }

            double[] forces = new double[dofCount];
            bool[] fixedDofs = new bool[dofCount];

            HashSet<int> supportNodes = new HashSet<int>();

            rc = BuildSupports(
                brepRef.ObjectId,
                supportMode,
                startX,
                startY,
                startZ,
                dx,
                dy,
                dz,
                nelx,
                nely,
                nelz,
                nny,
                nnz,
                fixedDofs,
                gridBox,
                tol,
                supportNodes);

            if (rc != Result.Success)
                return rc;

            bool[,,] designMask = MaskConnectivity3D.ExtractSupportConnectedMask(
                rawMask,
                supportNodes,
                bbox,
                nelx,
                nely,
                nelz,
                nny,
                nnz);

            int connectedActiveCount = VoxelDomainBuilder3D.CountActiveElements(designMask, nelx, nely, nelz);
            if (connectedActiveCount == 0)
            {
                RhinoApp.WriteLine("No support-connected active voxels remain.");
                return Result.Failure;
            }

            HashSet<int> activeNodes = BuildActiveNodeSet(designMask, nelx, nely, nelz, nny, nnz);

            rc = BuildLoad(
                brepRef.ObjectId,
                loadMode,
                startX,
                startY,
                startZ,
                dx,
                dy,
                dz,
                nelx,
                nely,
                nelz,
                nny,
                nnz,
                forces,
                gridBox,
                tol,
                activeNodes);

            if (rc != Result.Success)
                return rc;

            RhinoApp.WriteLine("Grid: {0} x {1} x {2}", nelx, nely, nelz);
            RhinoApp.WriteLine("Cell size used: {0}", cellSize);
            RhinoApp.WriteLine("DOF count: {0}", dofCount);
            RhinoApp.WriteLine("Raw active voxels: {0}", rawActiveCount);
            RhinoApp.WriteLine("Support-connected active voxels: {0}", connectedActiveCount);
            RhinoApp.WriteLine("Active node count: {0}", activeNodes.Count);
            RhinoApp.WriteLine("Load mode: {0}", loadMode);
            RhinoApp.WriteLine("Support mode: {0}", supportMode);
            RhinoApp.WriteLine("Support node count: {0}", supportNodes.Count);

            TopOptDomain3D domain = new TopOptDomain3D
            {
                DesignMask = designMask,
                Forces = forces,
                FixedDofs = fixedDofs
            };

            TopOptResult3D result;
            try
            {
                RhinoApp.WriteLine("Solver backend: Native C++ (Eigen)");
                result = NativeSolver3D.Solve(problem, domain);
            }
            catch (Exception ex)
            {
                RhinoApp.WriteLine("TopOpt3D solve failed.");
                RhinoApp.WriteLine(ex.Message);
                return Result.Failure;
            }

            if (result == null)
                return Result.Cancel;

            if (result.Density == null)
            {
                RhinoApp.WriteLine("Solver returned no result.");
                return Result.Failure;
            }

            _activeConduit = new TopletDisplayConduit(result, designMask, gridBox);
            _activeConduit.Enabled = true;

            _activeViewer = new ResultViewerForm(_activeConduit, doc);
            _activeViewer.Show();

            RhinoApp.WriteLine("TopOpt3D completed — use the Result Viewer to explore and bake.");
            RhinoApp.WriteLine("Iterations: {0}", result.Iterations);
            RhinoApp.WriteLine("Final compliance: {0:F4}", result.Compliance);

            doc.Views.Redraw();
            return Result.Success;
        }

        private static Result BuildLoad(
            Guid parentBrepId,
            LoadInputMode loadMode,
            double minX,
            double minY,
            double minZ,
            double dx,
            double dy,
            double dz,
            int nelx,
            int nely,
            int nelz,
            int nny,
            int nnz,
            double[] forces,
            BoundingBox bbox,
            double tol,
            HashSet<int> activeNodes)
        {
            if (loadMode == LoadInputMode.Point)
                return BuildPointLoad(minX, minY, minZ, dx, dy, dz, nelx, nely, nelz, nny, nnz, forces, activeNodes);

            return BuildFaceLoad(parentBrepId, bbox, nelx, nely, nelz, nny, nnz, tol, forces, activeNodes);
        }

        private static Result BuildPointLoad(
            double minX,
            double minY,
            double minZ,
            double dx,
            double dy,
            double dz,
            int nelx,
            int nely,
            int nelz,
            int nny,
            int nnz,
            double[] forces,
            HashSet<int> activeNodes)
        {
            GetPoint gp = new GetPoint();
            gp.SetCommandPrompt("Pick load point (global -Z point load)");
            GetResult gr = gp.Get();

            if (gr != GetResult.Point)
                return Result.Cancel;

            Point3d pt = gp.Point();

            int ix = ClampIndex((int)Math.Round((pt.X - minX) / dx), 0, nelx);
            int iy = ClampIndex((int)Math.Round((pt.Y - minY) / dy), 0, nely);
            int iz = ClampIndex((int)Math.Round((pt.Z - minZ) / dz), 0, nelz);

            int node = NodeIndex(ix, iy, iz, nny, nnz);

            if (!activeNodes.Contains(node))
            {
                RhinoApp.WriteLine("Selected load point mapped to a node that is not on the active support-connected structure.");
                return Result.Failure;
            }

            forces[3 * node + 2] = -1.0;

            RhinoApp.WriteLine("Point load node ijk: ({0}, {1}, {2})", ix, iy, iz);
            return Result.Success;
        }

        private static Result BuildFaceLoad(
            Guid parentBrepId,
            BoundingBox bbox,
            int nelx,
            int nely,
            int nelz,
            int nny,
            int nnz,
            double tol,
            double[] forces,
            HashSet<int> activeNodes)
        {
            ObjRef faceRef = GetSingleFaceFromBrep("Select ONE load face (Ctrl+Shift-click a face)", parentBrepId);
            if (faceRef == null)
                return Result.Cancel;

            BrepFace face = faceRef.Face();
            if (face == null)
            {
                RhinoApp.WriteLine("Could not read load face.");
                return Result.Failure;
            }

            List<int> candidateNodes = GetNodesNearFace(face, bbox, nelx, nely, nelz, nny, nnz, tol);
            List<int> nodes = FilterNodesToActive(candidateNodes, activeNodes);

            RhinoApp.WriteLine("Load face candidate node count: {0}", candidateNodes.Count);
            RhinoApp.WriteLine("Load face active node count: {0}", nodes.Count);

            if (nodes.Count == 0)
            {
                RhinoApp.WriteLine("No active FE nodes found near the selected load face.");
                RhinoApp.WriteLine("This usually means the face is not connected to the support-connected voxel structure at the current resolution.");
                return Result.Failure;
            }

            ApplyDistributedFaceLoad(forces, nodes, face, bbox, nelx, nely, nelz, nny, nnz, 1.0);

            RhinoApp.WriteLine("Load face node count: {0}", nodes.Count);
            return Result.Success;
        }

        private static Result BuildSupports(
            Guid parentBrepId,
            SupportInputMode supportMode,
            double minX,
            double minY,
            double minZ,
            double dx,
            double dy,
            double dz,
            int nelx,
            int nely,
            int nelz,
            int nny,
            int nnz,
            bool[] fixedDofs,
            BoundingBox bbox,
            double tol,
            HashSet<int> supportNodes)
        {
            if (supportMode == SupportInputMode.Points)
                return BuildPointSupports(minX, minY, minZ, dx, dy, dz, nelx, nely, nelz, nny, nnz, fixedDofs, supportNodes);

            return BuildFaceSupports(parentBrepId, bbox, nelx, nely, nelz, nny, nnz, tol, fixedDofs, supportNodes);
        }

        private static Result BuildPointSupports(
            double minX,
            double minY,
            double minZ,
            double dx,
            double dy,
            double dz,
            int nelx,
            int nely,
            int nelz,
            int nny,
            int nnz,
            bool[] fixedDofs,
            HashSet<int> supportNodes)
        {
            List<Point3d> supportPoints = new List<Point3d>();

            while (true)
            {
                GetPoint gp = new GetPoint();
                gp.SetCommandPrompt("Pick support point. Press Enter when done");
                gp.AcceptNothing(true);

                GetResult gr = gp.Get();

                if (gr == GetResult.Nothing)
                    break;

                if (gr != GetResult.Point)
                    return Result.Cancel;

                supportPoints.Add(gp.Point());
            }

            if (supportPoints.Count == 0)
            {
                RhinoApp.WriteLine("At least one support point is required.");
                return Result.Failure;
            }

            foreach (Point3d pt in supportPoints)
            {
                int ix = ClampIndex((int)Math.Round((pt.X - minX) / dx), 0, nelx);
                int iy = ClampIndex((int)Math.Round((pt.Y - minY) / dy), 0, nely);
                int iz = ClampIndex((int)Math.Round((pt.Z - minZ) / dz), 0, nelz);

                supportNodes.Add(NodeIndex(ix, iy, iz, nny, nnz));
            }

            ApplyFixedNodes(fixedDofs, supportNodes);
            return Result.Success;
        }

        private static Result BuildFaceSupports(
            Guid parentBrepId,
            BoundingBox bbox,
            int nelx,
            int nely,
            int nelz,
            int nny,
            int nnz,
            double tol,
            bool[] fixedDofs,
            HashSet<int> supportNodes)
        {
            ObjRef[] faceRefs = GetMultipleFacesFromBrep(
                "Select ONE OR MORE support faces (Ctrl+Shift-click faces, Enter when done)",
                parentBrepId);

            if (faceRefs == null || faceRefs.Length == 0)
            {
                RhinoApp.WriteLine("At least one support face is required.");
                return Result.Failure;
            }

            for (int i = 0; i < faceRefs.Length; i++)
            {
                BrepFace face = faceRefs[i].Face();
                if (face == null)
                    continue;

                List<int> nodes = GetNodesNearFace(face, bbox, nelx, nely, nelz, nny, nnz, tol);
                for (int j = 0; j < nodes.Count; j++)
                    supportNodes.Add(nodes[j]);
            }

            if (supportNodes.Count == 0)
            {
                RhinoApp.WriteLine("No FE nodes found near selected support faces.");
                return Result.Failure;
            }

            ApplyFixedNodes(fixedDofs, supportNodes);
            return Result.Success;
        }

        private static ObjRef GetSingleFaceFromBrep(string prompt, Guid parentBrepId)
        {
            GetObject go = new GetObject();
            go.SetCommandPrompt(prompt);
            go.GeometryFilter = ObjectType.Surface;
            go.GeometryAttributeFilter = GeometryAttributeFilter.SubSurface;
            go.SubObjectSelect = true;
            go.EnablePreSelect(false, true);
            go.DeselectAllBeforePostSelect = false;

            GetResult gr = go.Get();
            if (gr != GetResult.Object)
                return null;

            ObjRef faceRef = go.Object(0);
            if (faceRef == null)
                return null;

            if (faceRef.ObjectId != parentBrepId)
            {
                RhinoApp.WriteLine("Selected face must belong to the selected solid.");
                return null;
            }

            if (faceRef.Face() == null)
            {
                RhinoApp.WriteLine("You must select a Brep face.");
                return null;
            }

            return faceRef;
        }

        private static ObjRef[] GetMultipleFacesFromBrep(string prompt, Guid parentBrepId)
        {
            GetObject go = new GetObject();
            go.SetCommandPrompt(prompt);
            go.GeometryFilter = ObjectType.Surface;
            go.GeometryAttributeFilter = GeometryAttributeFilter.SubSurface;
            go.SubObjectSelect = true;
            go.AcceptNothing(true);
            go.EnablePreSelect(false, true);
            go.DeselectAllBeforePostSelect = false;

            GetResult gr = go.GetMultiple(1, 0);
            if (gr != GetResult.Object)
                return null;

            List<ObjRef> refs = new List<ObjRef>();

            for (int i = 0; i < go.ObjectCount; i++)
            {
                ObjRef faceRef = go.Object(i);
                if (faceRef == null)
                    continue;

                if (faceRef.ObjectId != parentBrepId)
                {
                    RhinoApp.WriteLine("All selected faces must belong to the selected solid.");
                    return null;
                }

                if (faceRef.Face() == null)
                {
                    RhinoApp.WriteLine("Selection must contain Brep faces only.");
                    return null;
                }

                refs.Add(faceRef);
            }

            return refs.ToArray();
        }

        private static List<int> GetNodesNearFace(
            BrepFace face,
            BoundingBox bbox,
            int nelx,
            int nely,
            int nelz,
            int nny,
            int nnz,
            double tol)
        {
            List<int> nodes = new List<int>();

            double dx = bbox.Diagonal.X / nelx;
            double dy = bbox.Diagonal.Y / nely;
            double dz = bbox.Diagonal.Z / nelz;

            double nodeTol = Math.Max(tol, 0.30 * Math.Min(dx, Math.Min(dy, dz)));

            for (int ix = 0; ix <= nelx; ix++)
            {
                for (int iy = 0; iy <= nely; iy++)
                {
                    for (int iz = 0; iz <= nelz; iz++)
                    {
                        Point3d p = new Point3d(
                            bbox.Min.X + ix * dx,
                            bbox.Min.Y + iy * dy,
                            bbox.Min.Z + iz * dz);

                        double u, v;
                        if (!face.ClosestPoint(p, out u, out v))
                            continue;

                        Point3d q = face.PointAt(u, v);
                        if (p.DistanceTo(q) <= nodeTol)
                            nodes.Add(NodeIndex(ix, iy, iz, nny, nnz));
                    }
                }
            }

            return nodes;
        }

        private static HashSet<int> BuildActiveNodeSet(
            bool[,,] mask,
            int nelx,
            int nely,
            int nelz,
            int nny,
            int nnz)
        {
            HashSet<int> activeNodes = new HashSet<int>();

            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    for (int ez = 0; ez < nelz; ez++)
                    {
                        if (!mask[ex, ey, ez])
                            continue;

                        activeNodes.Add(NodeIndex(ex, ey, ez, nny, nnz));
                        activeNodes.Add(NodeIndex(ex + 1, ey, ez, nny, nnz));
                        activeNodes.Add(NodeIndex(ex + 1, ey + 1, ez, nny, nnz));
                        activeNodes.Add(NodeIndex(ex, ey + 1, ez, nny, nnz));
                        activeNodes.Add(NodeIndex(ex, ey, ez + 1, nny, nnz));
                        activeNodes.Add(NodeIndex(ex + 1, ey, ez + 1, nny, nnz));
                        activeNodes.Add(NodeIndex(ex + 1, ey + 1, ez + 1, nny, nnz));
                        activeNodes.Add(NodeIndex(ex, ey + 1, ez + 1, nny, nnz));
                    }
                }
            }

            return activeNodes;
        }

        private static List<int> FilterNodesToActive(IEnumerable<int> nodes, HashSet<int> activeNodes)
        {
            List<int> filtered = new List<int>();

            foreach (int node in nodes)
            {
                if (activeNodes.Contains(node))
                    filtered.Add(node);
            }

            return filtered;
        }

        private static void ApplyDistributedFaceLoad(
            double[] forces,
            List<int> nodes,
            BrepFace face,
            BoundingBox bbox,
            int nelx,
            int nely,
            int nelz,
            int nny,
            int nnz,
            double totalLoadMagnitude)
        {
            if (forces == null || nodes == null || nodes.Count == 0 || face == null)
                return;

            double dx = bbox.Diagonal.X / nelx;
            double dy = bbox.Diagonal.Y / nely;
            double dz = bbox.Diagonal.Z / nelz;

            Vector3d avg = Vector3d.Zero;

            foreach (int node in nodes)
            {
                Point3d p = NodePoint(node, bbox, nny, nnz, dx, dy, dz);

                double u, v;
                if (!face.ClosestPoint(p, out u, out v))
                    continue;

                Vector3d n = face.NormalAt(u, v);
                if (face.OrientationIsReversed)
                    n.Reverse();

                if (n.Unitize())
                    avg += n;
            }

            if (!avg.Unitize())
            {
                Interval du = face.Domain(0);
                Interval dv = face.Domain(1);
                avg = face.NormalAt(0.5 * (du.T0 + du.T1), 0.5 * (dv.T0 + dv.T1));
                if (face.OrientationIsReversed)
                    avg.Reverse();
                avg.Unitize();
            }

            double fx = totalLoadMagnitude * avg.X / nodes.Count;
            double fy = totalLoadMagnitude * avg.Y / nodes.Count;
            double fz = totalLoadMagnitude * avg.Z / nodes.Count;

            foreach (int node in nodes)
            {
                forces[3 * node] += fx;
                forces[3 * node + 1] += fy;
                forces[3 * node + 2] += fz;
            }
        }

        private static void ApplyFixedNodes(bool[] fixedDofs, IEnumerable<int> nodes)
        {
            foreach (int node in nodes)
            {
                fixedDofs[3 * node] = true;
                fixedDofs[3 * node + 1] = true;
                fixedDofs[3 * node + 2] = true;
            }
        }

        private static Point3d NodePoint(
            int node,
            BoundingBox bbox,
            int nny,
            int nnz,
            double dx,
            double dy,
            double dz)
        {
            int ix = node / (nny * nnz);
            int rem = node % (nny * nnz);
            int iy = rem / nnz;
            int iz = rem % nnz;

            return new Point3d(
                bbox.Min.X + ix * dx,
                bbox.Min.Y + iy * dy,
                bbox.Min.Z + iz * dz);
        }

        private static int ClampIndex(int value, int min, int max)
        {
            if (value < min) return min;
            if (value > max) return max;
            return value;
        }

        private static int NodeIndex(int ix, int iy, int iz, int nny, int nnz)
        {
            return ix * (nny * nnz) + iy * nnz + iz;
        }

    }
}