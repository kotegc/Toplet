using Rhino;
using Rhino.Display;
using Rhino.DocObjects;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Drawing;

namespace Toplet_v0_Alpha.TopOpt3D
{
    public static class Visualization3D
    {
        public static Guid AddDensityMeshToDocument(
            RhinoDoc doc,
            TopOptResult3D result,
            BoundingBox bbox,
            double threshold,
            string layerName = "Toplet3D_Result",
            bool[,,] activeMask = null)
        {
            if (doc == null)
                throw new ArgumentNullException(nameof(doc));
            if (result == null)
                throw new ArgumentNullException(nameof(result));
            if (result.Density == null)
                throw new ArgumentException("Result density cannot be null.", nameof(result));
            if (!bbox.IsValid)
                throw new ArgumentException("Bounding box is invalid.", nameof(bbox));

            int nelx = result.NelX;
            int nely = result.NelY;
            int nelz = result.NelZ;

            if (nelx < 1 || nely < 1 || nelz < 1)
                throw new ArgumentException("Result grid dimensions must all be >= 1.");

            double dx = bbox.Diagonal.X / nelx;
            double dy = bbox.Diagonal.Y / nely;
            double dz = bbox.Diagonal.Z / nelz;

            int layerIndex = EnsureLayer(doc, layerName);

            ObjectAttributes attr = new ObjectAttributes();
            attr.LayerIndex = layerIndex;
            attr.ColorSource = ObjectColorSource.ColorFromObject;
            attr.ObjectColor = Color.Blue;
            ApplyGhostedOverrideToAllViews(doc, attr);

            Mesh mesh = new Mesh();

            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    for (int ez = 0; ez < nelz; ez++)
                    {
                        if (activeMask != null && !activeMask[ex, ey, ez])
                            continue;

                        double rho = result.Density[ex, ey, ez];
                        if (rho < threshold)
                            continue;

                        double x0 = bbox.Min.X + ex * dx;
                        double y0 = bbox.Min.Y + ey * dy;
                        double z0 = bbox.Min.Z + ez * dz;

                        double x1 = x0 + dx;
                        double y1 = y0 + dy;
                        double z1 = z0 + dz;

                        AddBoxToMesh(mesh, x0, y0, z0, x1, y1, z1);
                    }
                }
            }

            if (mesh.Vertices.Count == 0)
                return Guid.Empty;

            mesh.Normals.ComputeNormals();
            mesh.Compact();

            Guid id = doc.Objects.AddMesh(mesh, attr);
            doc.Views.Redraw();
            return id;
        }

        private static void AddBoxToMesh(Mesh mesh, double x0, double y0, double z0, double x1, double y1, double z1)
        {
            int v0 = mesh.Vertices.Add(x0, y0, z0);
            int v1 = mesh.Vertices.Add(x1, y0, z0);
            int v2 = mesh.Vertices.Add(x1, y1, z0);
            int v3 = mesh.Vertices.Add(x0, y1, z0);
            int v4 = mesh.Vertices.Add(x0, y0, z1);
            int v5 = mesh.Vertices.Add(x1, y0, z1);
            int v6 = mesh.Vertices.Add(x1, y1, z1);
            int v7 = mesh.Vertices.Add(x0, y1, z1);

            mesh.Faces.AddFace(v0, v1, v2, v3); // bottom
            mesh.Faces.AddFace(v4, v5, v6, v7); // top
            mesh.Faces.AddFace(v0, v1, v5, v4); // front
            mesh.Faces.AddFace(v1, v2, v6, v5); // right
            mesh.Faces.AddFace(v2, v3, v7, v6); // back
            mesh.Faces.AddFace(v3, v0, v4, v7); // left
        }

        private static void ApplyGhostedOverrideToAllViews(RhinoDoc doc, ObjectAttributes attr)
        {
            if (doc == null || attr == null)
                return;

            DisplayModeDescription ghosted = DisplayModeDescription.FindByName("Ghosted");
            if (ghosted == null)
                return;

            foreach (RhinoView view in doc.Views)
                attr.SetDisplayModeOverride(ghosted, view.ActiveViewportID);
        }

        private static int EnsureLayer(RhinoDoc doc, string layerName)
        {
            int existing = doc.Layers.FindByFullPath(layerName, -1);
            if (existing >= 0)
                return existing;

            Layer layer = new Layer
            {
                Name = layerName
            };

            return doc.Layers.Add(layer);
        }
    }
}