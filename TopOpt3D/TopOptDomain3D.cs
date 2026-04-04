namespace Toplet_v0_Alpha.TopOpt3D
{
    public class TopOptDomain3D
    {
        // One bool per voxel element.
        // true = active design element
        // false = ignored / nonexistent element
        public bool[,,] DesignMask { get; set; }

        // One bool per global DOF.
        // true = fixed
        // false = free
        public bool[] FixedDofs { get; set; }

        // One force value per global DOF.
        public double[] Forces { get; set; }
    }
}