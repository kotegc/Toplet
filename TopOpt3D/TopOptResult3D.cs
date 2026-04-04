namespace Toplet_v0_Alpha.TopOpt3D
{
    public class TopOptResult3D
    {
        public int NelX { get; set; }
        public int NelY { get; set; }
        public int NelZ { get; set; }

        public double[,,] Density { get; set; }
        public double Compliance { get; set; }
        public int Iterations { get; set; }
    }
}