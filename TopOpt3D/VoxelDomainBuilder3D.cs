using Rhino.Geometry;
using System;

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

            bool[,,] mask = new bool[nelx, nely, nelz];

            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    for (int ez = 0; ez < nelz; ez++)
                    {
                        double cx = bbox.Min.X + (ex + 0.5) * dx;
                        double cy = bbox.Min.Y + (ey + 0.5) * dy;
                        double cz = bbox.Min.Z + (ez + 0.5) * dz;

                        Point3d center = new Point3d(cx, cy, cz);
                        mask[ex, ey, ez] = brep.IsPointInside(center, tol, true);
                    }
                }
            }

            return mask;
        }

        public static int CountActiveElements(bool[,,] mask, int nelx, int nely, int nelz)
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
    }
}