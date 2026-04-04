using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace Toplet_v0_Alpha.TopOpt3D
{
    public static class MaskConnectivity3D
    {
        public static bool[,,] ExtractSupportConnectedMask(
            bool[,,] rawMask,
            IEnumerable<int> supportNodes,
            BoundingBox bbox,
            int nelx,
            int nely,
            int nelz,
            int nny,
            int nnz)
        {
            if (rawMask == null)
                throw new ArgumentNullException(nameof(rawMask));
            if (supportNodes == null)
                throw new ArgumentNullException(nameof(supportNodes));
            if (!bbox.IsValid)
                throw new ArgumentException("Bounding box is invalid.", nameof(bbox));

            bool[,,] connected = new bool[nelx, nely, nelz];
            Queue<(int ex, int ey, int ez)> queue = new Queue<(int ex, int ey, int ez)>();

            double dx = bbox.Diagonal.X / nelx;
            double dy = bbox.Diagonal.Y / nely;
            double dz = bbox.Diagonal.Z / nelz;

            foreach (int node in supportNodes)
            {
                int ix = node / (nny * nnz);
                int rem = node % (nny * nnz);
                int iy = rem / nnz;
                int iz = rem % nnz;

                AddSeedVoxel(ix - 1, iy - 1, iz - 1, rawMask, connected, nelx, nely, nelz, queue);
                AddSeedVoxel(ix - 1, iy - 1, iz, rawMask, connected, nelx, nely, nelz, queue);
                AddSeedVoxel(ix - 1, iy, iz - 1, rawMask, connected, nelx, nely, nelz, queue);
                AddSeedVoxel(ix - 1, iy, iz, rawMask, connected, nelx, nely, nelz, queue);
                AddSeedVoxel(ix, iy - 1, iz - 1, rawMask, connected, nelx, nely, nelz, queue);
                AddSeedVoxel(ix, iy - 1, iz, rawMask, connected, nelx, nely, nelz, queue);
                AddSeedVoxel(ix, iy, iz - 1, rawMask, connected, nelx, nely, nelz, queue);
                AddSeedVoxel(ix, iy, iz, rawMask, connected, nelx, nely, nelz, queue);
            }

            int[] nx = { -1, 1, 0, 0, 0, 0 };
            int[] ny = { 0, 0, -1, 1, 0, 0 };
            int[] nz = { 0, 0, 0, 0, -1, 1 };

            while (queue.Count > 0)
            {
                var cell = queue.Dequeue();

                for (int k = 0; k < 6; k++)
                {
                    int ex2 = cell.ex + nx[k];
                    int ey2 = cell.ey + ny[k];
                    int ez2 = cell.ez + nz[k];

                    if (ex2 < 0 || ex2 >= nelx ||
                        ey2 < 0 || ey2 >= nely ||
                        ez2 < 0 || ez2 >= nelz)
                        continue;

                    if (!rawMask[ex2, ey2, ez2])
                        continue;

                    if (connected[ex2, ey2, ez2])
                        continue;

                    connected[ex2, ey2, ez2] = true;
                    queue.Enqueue((ex2, ey2, ez2));
                }
            }

            return connected;
        }

        private static void AddSeedVoxel(
            int ex,
            int ey,
            int ez,
            bool[,,] rawMask,
            bool[,,] connected,
            int nelx,
            int nely,
            int nelz,
            Queue<(int ex, int ey, int ez)> queue)
        {
            if (ex < 0 || ex >= nelx ||
                ey < 0 || ey >= nely ||
                ez < 0 || ez >= nelz)
                return;

            if (!rawMask[ex, ey, ez])
                return;

            if (connected[ex, ey, ez])
                return;

            connected[ex, ey, ez] = true;
            queue.Enqueue((ex, ey, ez));
        }
    }
}