using Toplet_v0_Alpha.TopOpt3D;
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

        private enum BoxFaceType
        {
            Unknown,
            XMin,
            XMax,
            YMin,
            YMax,
            ZMin,
            ZMax
        }

        public Toplet3DCommand()
        {
            Instance = this;
        }

        public static Toplet3DCommand Instance { get; private set; }

        public override string EnglishName => "TopletTopOpt3D";

        protected override Result RunCommand(RhinoDoc doc, RunMode mode)
        {
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

            double minX = bbox.Min.X;
            double minY = bbox.Min.Y;
            double minZ = bbox.Min.Z;
            double maxX = bbox.Max.X;
            double maxY = bbox.Max.Y;
            double maxZ = bbox.Max.Z;

            double width = maxX - minX;
            double depth = maxY - minY;
            double height = maxZ - minZ;

            if (width <= 0.0 || depth <= 0.0 || height <= 0.0)
            {
                RhinoApp.WriteLine("Solid bounding box is invalid.");
                return Result.Failure;
            }

            LoadInputMode loadMode = LoadInputMode.Point;
            SupportInputMode supportMode = SupportInputMode.Points;

            rc = GetBoundaryConditionModes(ref loadMode, ref supportMode);
            if (rc != Result.Success)
                return rc;

            double cellSize = 1.0;

            int maxElementsX = 20;
            int maxElementsY = 20;
            int maxElementsZ = 20;

            if (width / cellSize > maxElementsX ||
                depth / cellSize > maxElementsY ||
                height / cellSize > maxElementsZ)
            {
                double cellSizeX = width / maxElementsX;
                double cellSizeY = depth / maxElementsY;
                double cellSizeZ = height / maxElementsZ;
                cellSize = Math.Max(cellSizeX, Math.Max(cellSizeY, cellSizeZ));
            }

            int nelx = Math.Max(1, (int)Math.Ceiling(width / cellSize));
            int nely = Math.Max(1, (int)Math.Ceiling(depth / cellSize));
            int nelz = Math.Max(1, (int)Math.Ceiling(height / cellSize));

            double dx = cellSize;
            double dy = cellSize;
            double dz = cellSize;

            TopOptProblem3D problem = new TopOptProblem3D
            {
                NelX = nelx,
                NelY = nely,
                NelZ = nelz,
                VolumeFraction = 0.5,
                Penal = 3.0,
                FilterRadius = 1.5,
                MaxIterations = 30
            };

            int nnx = nelx + 1;
            int nny = nely + 1;
            int nnz = nelz + 1;
            int nodeCount = nnx * nny * nnz;
            int dofCount = nodeCount * 3;

            if (dofCount > 12000)
            {
                RhinoApp.WriteLine("Grid is too large for the current dense 3D solver.");
                RhinoApp.WriteLine("Try a larger cell size or a smaller solid.");
                RhinoApp.WriteLine("Current DOF count: {0}", dofCount);
                return Result.Failure;
            }

            bool[,,] designMask = new bool[nelx, nely, nelz];
            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    for (int ez = 0; ez < nelz; ez++)
                    {
                        designMask[ex, ey, ez] = true;
                    }
                }
            }

            double[] forces = new double[dofCount];
            bool[] fixedDofs = new bool[dofCount];

            double tol = Math.Max(doc.ModelAbsoluteTolerance, cellSize * 0.25);

            rc = BuildLoad(
                doc,
                brepRef.ObjectId,
                bbox,
                loadMode,
                nelx,
                nely,
                nelz,
                dx,
                dy,
                dz,
                minX,
                minY,
                minZ,
                nny,
                nnz,
                tol,
                forces);

            if (rc != Result.Success)
                return rc;

            rc = BuildSupports(
                doc,
                brepRef.ObjectId,
                bbox,
                supportMode,
                nelx,
                nely,
                nelz,
                dx,
                dy,
                dz,
                minX,
                minY,
                minZ,
                nny,
                nnz,
                tol,
                fixedDofs);

            if (rc != Result.Success)
                return rc;

            RhinoApp.WriteLine("Grid: {0} x {1} x {2}", nelx, nely, nelz);
            RhinoApp.WriteLine("Cell size used: {0}", cellSize);
            RhinoApp.WriteLine("DOF count: {0}", dofCount);
            RhinoApp.WriteLine("Load mode: {0}", loadMode);
            RhinoApp.WriteLine("Support mode: {0}", supportMode);

            TopOptDomain3D domain = new TopOptDomain3D
            {
                DesignMask = designMask,
                Forces = forces,
                FixedDofs = fixedDofs
            };

            TopOptResult3D result;
            try
            {
                TopOptSolver3D solver = new TopOptSolver3D();
                result = solver.Solve(problem, domain);
            }
            catch (Exception ex)
            {
                RhinoApp.WriteLine("TopOpt3D solve failed.");
                RhinoApp.WriteLine(ex.Message);
                return Result.Failure;
            }

            if (result == null || result.Density == null)
            {
                RhinoApp.WriteLine("Solver returned no result.");
                return Result.Failure;
            }

            Visualization3D.AddDensityVoxelsToDocument(
                doc,
                result,
                bbox,
                0.3,
                "Toplet3D_Result");

            RhinoApp.WriteLine("TopOpt3D completed.");
            RhinoApp.WriteLine("Iterations: {0}", result.Iterations);
            RhinoApp.WriteLine("Final compliance: {0:F4}", result.Compliance);

            doc.Views.Redraw();
            return Result.Success;
        }

        private static Result GetBoundaryConditionModes(
            ref LoadInputMode loadMode,
            ref SupportInputMode supportMode)
        {
            GetOption go = new GetOption();
            int loadPointOpt = go.AddOption("LoadPoint");
            int loadFaceOpt = go.AddOption("LoadFace");
            int supportPointsOpt = go.AddOption("SupportPoints");
            int supportFacesOpt = go.AddOption("SupportFaces");
            go.AcceptNothing(true);
            go.SetCommandPrompt("Choose boundary condition modes. Press Enter to accept current modes");

            while (true)
            {
                string current = string.Format(
                    "Current: Load={0}, Supports={1}",
                    loadMode,
                    supportMode);

                go.SetCommandPrompt("Choose boundary condition modes. " + current);

                GetResult gr = go.Get();

                if (gr == GetResult.Nothing)
                    return Result.Success;

                if (gr != GetResult.Option)
                    return Result.Cancel;

                int idx = go.OptionIndex();
                if (idx == loadPointOpt) loadMode = LoadInputMode.Point;
                if (idx == loadFaceOpt) loadMode = LoadInputMode.Face;
                if (idx == supportPointsOpt) supportMode = SupportInputMode.Points;
                if (idx == supportFacesOpt) supportMode = SupportInputMode.Faces;
            }
        }

        private static Result BuildLoad(
            RhinoDoc doc,
            Guid parentBrepId,
            BoundingBox bbox,
            LoadInputMode loadMode,
            int nelx,
            int nely,
            int nelz,
            double dx,
            double dy,
            double dz,
            double minX,
            double minY,
            double minZ,
            int nny,
            int nnz,
            double tol,
            double[] forces)
        {
            if (loadMode == LoadInputMode.Point)
                return BuildPointLoad(minX, minY, minZ, dx, dy, dz, nelx, nely, nelz, nny, nnz, forces);

            return BuildFaceLoad(parentBrepId, bbox, nelx, nely, nelz, nny, nnz, tol, forces);
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
            double[] forces)
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
            double[] forces)
        {
            ObjRef faceRef = GetSingleFaceFromBrep("Select ONE load face", parentBrepId);
            if (faceRef == null)
                return Result.Cancel;

            BrepFace face = faceRef.Face();
            if (face == null)
            {
                RhinoApp.WriteLine("Could not read load face.");
                return Result.Failure;
            }

            BoxFaceType faceType = GetBoxFaceType(face, bbox, tol);
            if (faceType == BoxFaceType.Unknown)
            {
                RhinoApp.WriteLine("Load face must lie on a bounding-box side of the selected solid.");
                return Result.Failure;
            }

            List<int> nodes = GetNodesOnFace(faceType, nelx, nely, nelz, nny, nnz);
            if (nodes.Count == 0)
            {
                RhinoApp.WriteLine("No FE nodes found on the selected load face.");
                return Result.Failure;
            }

            Vector3d normal = GetFaceNormal(face);
            if (!normal.Unitize())
            {
                RhinoApp.WriteLine("Could not compute a valid load face normal.");
                return Result.Failure;
            }

            ApplyDistributedFaceLoad(forces, nodes, normal, 1.0);

            RhinoApp.WriteLine("Load face type: {0}", faceType);
            RhinoApp.WriteLine("Load face node count: {0}", nodes.Count);
            RhinoApp.WriteLine("Load normal: ({0:F3}, {1:F3}, {2:F3})", normal.X, normal.Y, normal.Z);

            return Result.Success;
        }

        private static Result BuildSupports(
            RhinoDoc doc,
            Guid parentBrepId,
            BoundingBox bbox,
            SupportInputMode supportMode,
            int nelx,
            int nely,
            int nelz,
            double dx,
            double dy,
            double dz,
            double minX,
            double minY,
            double minZ,
            int nny,
            int nnz,
            double tol,
            bool[] fixedDofs)
        {
            if (supportMode == SupportInputMode.Points)
                return BuildPointSupports(minX, minY, minZ, dx, dy, dz, nelx, nely, nelz, nny, nnz, fixedDofs);

            return BuildFaceSupports(parentBrepId, bbox, nelx, nely, nelz, nny, nnz, tol, fixedDofs);
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
            bool[] fixedDofs)
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

                Point3d pt = gp.Point();
                supportPoints.Add(pt);

                RhinoApp.WriteLine(
                    "Support point added: ({0:F3}, {1:F3}, {2:F3})",
                    pt.X, pt.Y, pt.Z);
            }

            if (supportPoints.Count == 0)
            {
                RhinoApp.WriteLine("At least one support point is required.");
                return Result.Failure;
            }

            HashSet<int> nodes = new HashSet<int>();

            foreach (Point3d pt in supportPoints)
            {
                int ix = ClampIndex((int)Math.Round((pt.X - minX) / dx), 0, nelx);
                int iy = ClampIndex((int)Math.Round((pt.Y - minY) / dy), 0, nely);
                int iz = ClampIndex((int)Math.Round((pt.Z - minZ) / dz), 0, nelz);

                int node = NodeIndex(ix, iy, iz, nny, nnz);
                nodes.Add(node);
            }

            ApplyFixedFaceSupports(fixedDofs, nodes);
            RhinoApp.WriteLine("Support point node count: {0}", nodes.Count);

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
            bool[] fixedDofs)
        {
            ObjRef[] faceRefs = GetMultipleFacesFromBrep(
                "Select ONE OR MORE support faces. Press Enter when done",
                parentBrepId);

            if (faceRefs == null || faceRefs.Length == 0)
            {
                RhinoApp.WriteLine("At least one support face is required.");
                return Result.Failure;
            }

            HashSet<int> supportNodeSet = new HashSet<int>();

            for (int i = 0; i < faceRefs.Length; i++)
            {
                BrepFace face = faceRefs[i].Face();
                if (face == null)
                    continue;

                BoxFaceType faceType = GetBoxFaceType(face, bbox, tol);
                if (faceType == BoxFaceType.Unknown)
                {
                    RhinoApp.WriteLine("A support face must lie on a bounding-box side of the selected solid.");
                    return Result.Failure;
                }

                List<int> nodes = GetNodesOnFace(faceType, nelx, nely, nelz, nny, nnz);
                for (int n = 0; n < nodes.Count; n++)
                    supportNodeSet.Add(nodes[n]);
            }

            if (supportNodeSet.Count == 0)
            {
                RhinoApp.WriteLine("No FE nodes found on selected support faces.");
                return Result.Failure;
            }

            ApplyFixedFaceSupports(fixedDofs, supportNodeSet);
            RhinoApp.WriteLine("Support face node count: {0}", supportNodeSet.Count);

            return Result.Success;
        }

        private static ObjRef GetSingleFaceFromBrep(string prompt, Guid parentBrepId)
        {
            GetObject go = new GetObject();
            go.SetCommandPrompt(prompt);
            go.GeometryFilter = ObjectType.Surface;
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
                    RhinoApp.WriteLine("All support faces must belong to the selected solid.");
                    return null;
                }

                if (faceRef.Face() == null)
                {
                    RhinoApp.WriteLine("Support selection must contain Brep faces only.");
                    return null;
                }

                refs.Add(faceRef);
            }

            return refs.ToArray();
        }

        private static BoxFaceType GetBoxFaceType(BrepFace face, BoundingBox bbox, double tol)
        {
            BoundingBox fb = face.GetBoundingBox(true);

            bool atXMin = Math.Abs(fb.Min.X - bbox.Min.X) <= tol && Math.Abs(fb.Max.X - bbox.Min.X) <= tol;
            bool atXMax = Math.Abs(fb.Min.X - bbox.Max.X) <= tol && Math.Abs(fb.Max.X - bbox.Max.X) <= tol;

            bool atYMin = Math.Abs(fb.Min.Y - bbox.Min.Y) <= tol && Math.Abs(fb.Max.Y - bbox.Min.Y) <= tol;
            bool atYMax = Math.Abs(fb.Min.Y - bbox.Max.Y) <= tol && Math.Abs(fb.Max.Y - bbox.Max.Y) <= tol;

            bool atZMin = Math.Abs(fb.Min.Z - bbox.Min.Z) <= tol && Math.Abs(fb.Max.Z - bbox.Min.Z) <= tol;
            bool atZMax = Math.Abs(fb.Min.Z - bbox.Max.Z) <= tol && Math.Abs(fb.Max.Z - bbox.Max.Z) <= tol;

            if (atXMin) return BoxFaceType.XMin;
            if (atXMax) return BoxFaceType.XMax;
            if (atYMin) return BoxFaceType.YMin;
            if (atYMax) return BoxFaceType.YMax;
            if (atZMin) return BoxFaceType.ZMin;
            if (atZMax) return BoxFaceType.ZMax;

            return BoxFaceType.Unknown;
        }

        private static List<int> GetNodesOnFace(
            BoxFaceType faceType,
            int nelx,
            int nely,
            int nelz,
            int nny,
            int nnz)
        {
            List<int> nodes = new List<int>();

            switch (faceType)
            {
                case BoxFaceType.XMin:
                    for (int iy = 0; iy <= nely; iy++)
                        for (int iz = 0; iz <= nelz; iz++)
                            nodes.Add(NodeIndex(0, iy, iz, nny, nnz));
                    break;

                case BoxFaceType.XMax:
                    for (int iy = 0; iy <= nely; iy++)
                        for (int iz = 0; iz <= nelz; iz++)
                            nodes.Add(NodeIndex(nelx, iy, iz, nny, nnz));
                    break;

                case BoxFaceType.YMin:
                    for (int ix = 0; ix <= nelx; ix++)
                        for (int iz = 0; iz <= nelz; iz++)
                            nodes.Add(NodeIndex(ix, 0, iz, nny, nnz));
                    break;

                case BoxFaceType.YMax:
                    for (int ix = 0; ix <= nelx; ix++)
                        for (int iz = 0; iz <= nelz; iz++)
                            nodes.Add(NodeIndex(ix, nely, iz, nny, nnz));
                    break;

                case BoxFaceType.ZMin:
                    for (int ix = 0; ix <= nelx; ix++)
                        for (int iy = 0; iy <= nely; iy++)
                            nodes.Add(NodeIndex(ix, iy, 0, nny, nnz));
                    break;

                case BoxFaceType.ZMax:
                    for (int ix = 0; ix <= nelx; ix++)
                        for (int iy = 0; iy <= nely; iy++)
                            nodes.Add(NodeIndex(ix, iy, nelz, nny, nnz));
                    break;
            }

            return nodes;
        }

        private static Vector3d GetFaceNormal(BrepFace face)
        {
            Interval du = face.Domain(0);
            Interval dv = face.Domain(1);

            double u = 0.5 * (du.T0 + du.T1);
            double v = 0.5 * (dv.T0 + dv.T1);

            Vector3d normal = face.NormalAt(u, v);
            if (face.OrientationIsReversed)
                normal.Reverse();

            normal.Unitize();
            return normal;
        }

        private static void ApplyDistributedFaceLoad(
            double[] forces,
            List<int> nodes,
            Vector3d normal,
            double totalLoadMagnitude)
        {
            if (forces == null || nodes == null || nodes.Count == 0)
                return;

            double fx = totalLoadMagnitude * normal.X / nodes.Count;
            double fy = totalLoadMagnitude * normal.Y / nodes.Count;
            double fz = totalLoadMagnitude * normal.Z / nodes.Count;

            for (int i = 0; i < nodes.Count; i++)
            {
                int node = nodes[i];
                forces[3 * node] += fx;
                forces[3 * node + 1] += fy;
                forces[3 * node + 2] += fz;
            }
        }

        private static void ApplyFixedFaceSupports(bool[] fixedDofs, IEnumerable<int> nodes)
        {
            if (fixedDofs == null || nodes == null)
                return;

            foreach (int node in nodes)
            {
                fixedDofs[3 * node] = true;
                fixedDofs[3 * node + 1] = true;
                fixedDofs[3 * node + 2] = true;
            }
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