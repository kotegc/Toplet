using Rhino;
using Rhino.DocObjects;
using Rhino.Geometry;

namespace Toplet_v0_Alpha.TopOpt3D
{
    public static class Visualization3D
    {
        public static int EnsureLayer(RhinoDoc doc, string layerName)
        {
            int existing = doc.Layers.FindByFullPath(layerName, -1);
            if (existing >= 0)
                return existing;

            return doc.Layers.Add(new Layer { Name = layerName });
        }
    }
}
