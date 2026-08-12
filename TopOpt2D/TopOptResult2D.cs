namespace Toplet_v0_Alpha.TopOpt2D
{
    public class TopOptResult2D
    {
        public int NelX { get; set; }
        public int NelY { get; set; }
        public double[,] Density { get; set; }
        public double Compliance { get; set; }
        public int Iterations { get; set; }
    }
}