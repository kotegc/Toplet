using System;
using System.Collections.Generic;
using System.Drawing;
using Rhino.Display;
using Rhino.Geometry;

namespace Toplet_v0_Alpha.TopOpt3D
{
    public enum TopletDisplayMode { DensityHeatmap, SolidVoid, Ghosted }

    public sealed class TopletDisplayConduit : DisplayConduit
    {
        private readonly TopOptResult3D _result;
        private readonly bool[,,]       _mask;
        private readonly BoundingBox    _gridBox;

        private TopletDisplayMode _mode      = TopletDisplayMode.DensityHeatmap;
        private double            _threshold = 0.50;
        private Mesh              _mesh;
        private bool              _dirty     = true;

        // Face directions: +X -X +Y -Y +Z -Z
        private static readonly int[,] FaceDir = {
            { 1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}
        };

        // Corner offsets for each face, wound so the normal points outward
        private static readonly int[,,] FaceCorners = {
            {{1,0,0},{1,1,0},{1,1,1},{1,0,1}},  // +X
            {{0,1,0},{0,0,0},{0,0,1},{0,1,1}},  // -X
            {{1,1,0},{0,1,0},{0,1,1},{1,1,1}},  // +Y
            {{0,0,0},{1,0,0},{1,0,1},{0,0,1}},  // -Y
            {{0,0,1},{1,0,1},{1,1,1},{0,1,1}},  // +Z
            {{1,0,0},{0,0,0},{0,1,0},{1,1,0}},  // -Z
        };

        public TopletDisplayConduit(TopOptResult3D result, bool[,,] mask, BoundingBox gridBox)
        {
            _result  = result;
            _mask    = mask;
            _gridBox = gridBox;
        }

        public TopletDisplayMode Mode
        {
            get => _mode;
            set { if (_mode != value) { _mode = value; _dirty = true; } }
        }

        public double Threshold
        {
            get => _threshold;
            set { if (Math.Abs(_threshold - value) > 1e-9) { _threshold = value; _dirty = true; } }
        }

        protected override void DrawForeground(DrawEventArgs e)
        {
            if (_dirty) { RebuildMesh(); _dirty = false; }
            if (_mesh == null || !_mesh.IsValid) return;

            if (_mode == TopletDisplayMode.Ghosted)
            {
                var mat = new DisplayMaterial();
                mat.Diffuse      = Color.FromArgb(255, 80, 60, 220);
                mat.Transparency = 0.55;
                e.Display.DrawMeshShaded(_mesh, mat);
                e.Display.DrawMeshWires(_mesh, Color.FromArgb(50, 30, 80, 180));
            }
            else
            {
                e.Display.DrawMeshFalseColors(_mesh);
            }
        }

        private void RebuildMesh()
        {
            int nelx = _result.NelX, nely = _result.NelY, nelz = _result.NelZ;
            double dx = _gridBox.Diagonal.X / nelx;
            double dy = _gridBox.Diagonal.Y / nely;
            double dz = _gridBox.Diagonal.Z / nelz;
            double ox = _gridBox.Min.X, oy = _gridBox.Min.Y, oz = _gridBox.Min.Z;

            var mesh   = new Mesh();
            var colors = new List<Color>();

            for (int ex = 0; ex < nelx; ex++)
            for (int ey = 0; ey < nely; ey++)
            for (int ez = 0; ez < nelz; ez++)
            {
                if (!_mask[ex, ey, ez]) continue;

                double density = _result.Density[ex, ey, ez];

                bool isHeatmap = _mode == TopletDisplayMode.DensityHeatmap;
                if (!isHeatmap && density < _threshold) continue;

                Color faceColor = GetColor(density);

                for (int f = 0; f < 6; f++)
                {
                    int nx = ex + FaceDir[f, 0];
                    int ny = ey + FaceDir[f, 1];
                    int nz = ez + FaceDir[f, 2];

                    bool neighborVisible =
                        nx >= 0 && nx < nelx &&
                        ny >= 0 && ny < nely &&
                        nz >= 0 && nz < nelz &&
                        _mask[nx, ny, nz] &&
                        (isHeatmap || _result.Density[nx, ny, nz] >= _threshold);

                    if (neighborVisible) continue;

                    int baseIdx = mesh.Vertices.Count;
                    for (int c = 0; c < 4; c++)
                    {
                        mesh.Vertices.Add(new Point3f(
                            (float)(ox + (ex + FaceCorners[f, c, 0]) * dx),
                            (float)(oy + (ey + FaceCorners[f, c, 1]) * dy),
                            (float)(oz + (ez + FaceCorners[f, c, 2]) * dz)));
                        colors.Add(faceColor);
                    }
                    mesh.Faces.AddFace(baseIdx, baseIdx+1, baseIdx+2, baseIdx+3);
                }
            }

            mesh.VertexColors.SetColors(colors.ToArray());
            mesh.Normals.ComputeNormals();
            mesh.Compact();
            _mesh = mesh.IsValid ? mesh : null;
        }

        private Color GetColor(double density)
        {
            if (_mode == TopletDisplayMode.DensityHeatmap) return HeatmapColor(density);
            if (_mode == TopletDisplayMode.SolidVoid)      return Color.FromArgb(45, 45, 55);
            return Color.FromArgb(60, 130, 210);
        }

        private static Color HeatmapColor(double t)
        {
            t = Math.Max(0.0, Math.Min(1.0, t));
            double r, g, b;
            if      (t < 0.25) { r = 0;               g = t / 0.25;          b = 1; }
            else if (t < 0.50) { r = 0;               g = 1;                 b = 1 - (t - 0.25) / 0.25; }
            else if (t < 0.75) { r = (t - 0.50)/0.25; g = 1;                 b = 0; }
            else               { r = 1;               g = 1 - (t-0.75)/0.25; b = 0; }
            return Color.FromArgb((int)(r*255), (int)(g*255), (int)(b*255));
        }

        public Mesh BakeMesh()
        {
            if (_dirty) { RebuildMesh(); _dirty = false; }
            return _mesh;
        }
    }
}
