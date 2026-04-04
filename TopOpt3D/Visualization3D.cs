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
        public static uint[] AddDensityVoxelsToDocument(
            RhinoDoc doc,
            TopOptResult3D result,
            BoundingBox bbox,
            double threshold,
            string layerName = "Toplet3D_Result")
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

            if (dx <= 0.0 || dy <= 0.0 || dz <= 0.0)
                throw new ArgumentException("Computed voxel size is invalid.");

            int layerIndex = EnsureLayer(doc, layerName);

            ObjectAttributes baseAttr = new ObjectAttributes();
            baseAttr.LayerIndex = layerIndex;
            baseAttr.ColorSource = ObjectColorSource.ColorFromObject;
            baseAttr.ObjectColor = Color.Blue;
            baseAttr.MaterialSource = ObjectMaterialSource.MaterialFromObject;

            List<uint> added = new List<uint>();

            for (int ex = 0; ex < nelx; ex++)
            {
                for (int ey = 0; ey < nely; ey++)
                {
                    for (int ez = 0; ez < nelz; ez++)
                    {
                        double rho = result.Density[ex, ey, ez];
                        if (rho < threshold)
                            continue;

                        double x0 = bbox.Min.X + ex * dx;
                        double y0 = bbox.Min.Y + ey * dy;
                        double z0 = bbox.Min.Z + ez * dz;

                        double x1 = x0 + dx;
                        double y1 = y0 + dy;
                        double z1 = z0 + dz;

                        Point3d[] corners = new Point3d[]
                        {
                            new Point3d(x0, y0, z0),
                            new Point3d(x1, y0, z0),
                            new Point3d(x1, y1, z0),
                            new Point3d(x0, y1, z0),
                            new Point3d(x0, y0, z1),
                            new Point3d(x1, y0, z1),
                            new Point3d(x1, y1, z1),
                            new Point3d(x0, y1, z1)
                        };

                        Brep brep = Brep.CreateFromBox(corners);
                        if (brep == null)
                            continue;

                        ObjectAttributes attr = baseAttr.Duplicate();
                        ApplyGhostedOverrideToAllViews(doc, attr);

                        Guid id = doc.Objects.AddBrep(brep, attr);
                        if (id == Guid.Empty)
                            continue;

                        RhinoObject obj = doc.Objects.FindId(id);
                        if (obj != null)
                            added.Add(obj.RuntimeSerialNumber);
                    }
                }
            }

            doc.Views.Redraw();
            return added.ToArray();
        }

        private static void ApplyGhostedOverrideToAllViews(RhinoDoc doc, ObjectAttributes attr)
        {
            if (doc == null || attr == null)
                return;

            DisplayModeDescription ghosted = DisplayModeDescription.FindByName("Ghosted");
            if (ghosted == null)
                return;

            foreach (RhinoView view in doc.Views)
            {
                attr.SetDisplayModeOverride(ghosted, view.ActiveViewportID);
            }
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