namespace Anemone_V1_Alpha.TopOpt
{
    public class TopOptResult
    {
        public int NelX { get; set; }
        public int NelY { get; set; }
        public double[,] Density { get; set; }
        public double Compliance { get; set; }
        public int Iterations { get; set; }
    }
}