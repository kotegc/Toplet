namespace Anemone_V1_Alpha.TopOpt
{
    public class TopOptProblem
    {
        public int NelX { get; set; }
        public int NelY { get; set; }
        public double VolumeFraction { get; set; }
        public double Penal { get; set; }
        public double FilterRadius { get; set; }
        public int MaxIterations { get; set; }

        public double YoungsModulusSolid { get; set; } = 1.0;
        public double YoungsModulusMin { get; set; } = 1e-9;
        public double PoissonRatio { get; set; } = 0.3;
    }
}