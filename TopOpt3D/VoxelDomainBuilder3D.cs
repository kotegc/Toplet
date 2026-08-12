using System;
using System.Threading.Tasks;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;

namespace Toplet_v0_Alpha.TopOpt3D
{
    public static class VoxelDomainBuilder3D
    {
        public static bool[,,] BuildMaskFromBrep(
            Brep brep,
            BoundingBox bbox,
            int nelx,
            int nely,
            int nelz,
            double dx,
            double dy,
            double dz,
            double tol)
        {
            if (brep == null)
                throw new ArgumentNullException(nameof(brep));
            if (!bbox.IsValid)
                throw new ArgumentException("Bounding box is invalid.", nameof(bbox));
            if (nelx < 1 || nely < 1 || nelz < 1)
                throw new ArgumentException("Grid dimensions must all be >= 1.");
            if (dx <= 0.0 || dy <= 0.0 || dz <= 0.0)
                throw new ArgumentException("Voxel size must be positive.");

            // Tessellate Brep once — max edge length half the smallest voxel dimension
            // so no ray can slip through a gap between triangles
            double minVoxel = Math.Min(Math.Min(dx, dy), dz);
            var mp = new MeshingParameters();
            mp.MaximumEdgeLength = minVoxel * 0.5;
            mp.MinimumEdgeLength = tol;

            Mesh[] parts = Mesh.CreateFromBrep(brep, mp);
            if (parts == null || parts.Length == 0)
                return FallbackIsPointInside(brep, bbox, nelx, nely, nelz, dx, dy, dz, tol);

            var combined = new Mesh();
            foreach (var part in parts) combined.Append(part);
            combined.Normals.ComputeNormals();
            combined.Compact();

            bool[,,] mask = new bool[nelx, nely, nelz];

            // One ray per (ex, ey) column — parallelised across columns
            Parallel.For(0, nelx * nely, idx =>
            {
                int ex = idx / nely;
                int ey = idx % nely;

                double cx = bbox.Min.X + (ex + 0.5) * dx;
                double cy = bbox.Min.Y + (ey + 0.5) * dy;
                double originZ = bbox.Min.Z - 1.0;

                var line = new Line(
                    new Point3d(cx, cy, originZ),
                    new Point3d(cx, cy, bbox.Max.Z + 1.0));

                Point3d[] hits = Intersection.MeshLine(combined, line, out _);

                if (hits == null || hits.Length == 0) return;

                double[] zHits = new double[hits.Length];
                for (int i = 0; i < hits.Length; i++)
                    zHits[i] = hits[i].Z;
                Array.Sort(zHits);

                // Parity test for each voxel in this column
                for (int ez = 0; ez < nelz; ez++)
                {
                    double cz = bbox.Min.Z + (ez + 0.5) * dz;
                    int crossings = 0;
                    foreach (double z in zHits)
                        if (z < cz) crossings++;
                    mask[ex, ey, ez] = (crossings & 1) == 1;
                }
            });

            return mask;
        }

        // Fallback if mesh tessellation fails
        private static bool[,,] FallbackIsPointInside(
            Brep brep, BoundingBox bbox,
            int nelx, int nely, int nelz,
            double dx, double dy, double dz, double tol)
        {
            bool[,,] mask = new bool[nelx, nely, nelz];
            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++)
            for (int ez = 0; ez < nelz; ez++)
            {
                double cx = bbox.Min.X + (ex + 0.5) * dx;
                double cy = bbox.Min.Y + (ey + 0.5) * dy;
                double cz = bbox.Min.Z + (ez + 0.5) * dz;
                mask[ex, ey, ez] = brep.IsPointInside(new Point3d(cx, cy, cz), tol, true);
            }
            return mask;
        }

        public static int CountActiveElements(bool[,,] mask, int nelx, int nely, int nelz)
        {
            int count = 0;
            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++)
            for (int ez = 0; ez < nelz; ez++)
                if (mask[ex, ey, ez]) count++;
            return count;
        }
    }
}
