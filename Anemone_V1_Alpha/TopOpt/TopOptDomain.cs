namespace Anemone_V1_Alpha.TopOpt
{
    public class TopOptDomain
    {
        // One bool per element/cell in the design grid.
        // true = this element is part of the design domain
        // false = this element does not exist / should be ignored
        public bool[,] DesignMask { get; set; }

        // One bool per global DOF.
        // true = fixed
        // false = free
        public bool[] FixedDofs { get; set; }

        // One force value per global DOF.
        public double[] Forces { get; set; }
    }
}