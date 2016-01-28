package bdsky;

/**
 * A piecewise constant segment of a skyline.
 */
public class BDSSkylineSegment extends SkylineSegment {

    public BDSSkylineSegment(double lambda, double mu, double psi, double r, double t1, double t2) {

        super(t1, t2, new double[]{lambda, mu, psi, r});
    }

    /**
     * @return the birth rate per unit time.
     */
    public double lambda() { return value[0]; };

    /**
     * @return the death rate per unit time.
     */
    public double mu() { return value[1]; };

    /**
     * @return the sampling rate per unit time.
     */
    public double psi() { return value[2]; };

    /**
     * @return the removal probability, i.e. the probability that sampling causes recovery/removal/death.
     */
    public double r() { return value[3]; };

    // TODO fold rho sampling events into the Skyline
    // public boolean hasRho();
}
