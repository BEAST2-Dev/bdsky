package bdsky;

import beast.core.Input;

/**
 * Created by alexeid on 7/12/15.
 */
public class CanonicalParameterization extends BDSParameterization {

    public Input<SimpleSkyline> birthRate =
            new Input<>("birthRate", "BirthRate = BirthRateVector * birthRateScalar, birthrate can change over time");
    public Input<SimpleSkyline> deathRate =
            new Input<>("deathRate", "The deathRate vector with birthRates between times");
    public Input<SimpleSkyline> samplingRate =
            new Input<>("samplingRate", "The sampling rate per individual");      // psi
    public Input<SimpleSkyline> removalProbability =
            new Input<>("removalProbability", "The probability of an individual to become noninfectious immediately after the sampling");

    @Override
    public void initAndValidate() {

        // Need to add defaults, not all inputs are required

        MultiSkyline multiSkyline = new MultiSkyline(
                birthRate.get(),
                deathRate.get(),
                samplingRate.get(),
                removalProbability.get()
        );
        setMultiSkyline(multiSkyline);
    }


    @Override
    public BDSSkylineSegment toCanonicalSegment(SkylineSegment segment) {

        double birth = segment.value[0]; // lambda = birth rate
        double death = segment.value[1]; // mu = death rate
        double psi = segment.value[2]; // psi = sampling rate
        double r = segment.value[3]; // removal probability

        return new BDSSkylineSegment(birth, death, psi, r, segment.t1, segment.t2);
    }

}
