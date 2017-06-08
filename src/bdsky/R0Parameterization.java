package bdsky;

import beast.core.Input;

/**
 * Created by alexeid on 7/12/15.
 */
public class R0Parameterization extends BDSParameterization {

    public Input<SimpleSkyline> R0 =
              new Input<>("R0",
                      "The skyline of the basic reproduction number");
    public Input<SimpleSkyline> becomeUninfectiousRate =
              new Input<>("becomeUninfectiousRate",
                      "Rate at which individuals become uninfectious (through recovery or sampling)");
    public Input<SimpleSkyline> samplingProportion =
              new Input<>("samplingProportion",
                      "The samplingProportion = samplingRate / becomeUninfectiousRate");
    public Input<SimpleSkyline> removalProbability =
            new Input<>("removalProbability",
                    "The probability of death/removal/recovery upon sampling. " +
                            "If 1.0 then no sampled ancestors are produced in that interval.");

    @Override
    public void initAndValidate() {

        // Need to add defaults, not all inputs are required

        MultiSkyline multiSkyline = new MultiSkyline(
                R0.get(),
                becomeUninfectiousRate.get(),
                samplingProportion.get(),
                removalProbability.get()
        );
        setMultiSkyline(multiSkyline);
    }


    @Override
    public BDSSkylineSegment toCanonicalSegment(SkylineSegment segment) {

        double R = segment.value[0]; // R
        double b = segment.value[1]; // become uninfectious
        double p = segment.value[2]; // sampling proportion
        double r = segment.value[3]; // removal probability

        double birth =  R * b;
        double psi = p * b;
        double death = b - psi*r;


        return new BDSSkylineSegment(birth, death, psi, r, segment.t1, segment.t2);
    }

}
