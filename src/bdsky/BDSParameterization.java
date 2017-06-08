package bdsky;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;

import java.util.ArrayList;
import java.util.List;

/**
 * A parameterization of the birth-death skyline model
 */
public abstract class BDSParameterization extends CalculationNode {

    MultiSkyline multiSkyline;

    public Input<RealParameter> origin =
            new Input<RealParameter>("origin", "The time from origin to last sample (must be larger than tree height)", Input.Validate.REQUIRED);

    public final void setMultiSkyline(MultiSkyline multiSkyline) {
        this.multiSkyline = multiSkyline;
    }

    /**
     * @return the canonical segments for this skyline model.
     */
    public final List<BDSSkylineSegment> canonicalSegments(){

        List<BDSSkylineSegment> canonical = new ArrayList<>();

        for (SkylineSegment seg : multiSkyline.getSegments()) {
            canonical.add(toCanonicalSegment(seg));
        }
        return canonical;
    }

    /**
     * @return the number of segments in this parameterization.
     */
    public final int size() {
        return multiSkyline.getSegments().size();
    }

    public final void populateCanonical(Double[] birth, Double[] death, Double[] psi, Double[] r, Double[] times) {
        int size = size();

        if (birth.length != size || death.length != size || psi.length != size || r.length != size) {
            throw new RuntimeException("array size unexpected!");
        }
        List<BDSSkylineSegment> canonicalSegments = canonicalSegments();
        for (int i = 0; i < size; i++) {
            BDSSkylineSegment seg  = canonicalSegments.get(i);
            birth[i] = seg.lambda();
            death[i] = seg.mu();
            psi[i] = seg.psi();
            r[i] = seg.r();
            times[i] = seg.start();
        }
    }

    /**
     * @return the time of the origin of the process before the present.
     */
    public double origin() {
        return origin.get().getValue();
    }

    /**
     * @return true if any segments have a removalProbability < 1
     */
    public final boolean isSampledAncestorModel() {
        for (BDSSkylineSegment seg : canonicalSegments()) {
            if (seg.r() < 1.0) return true;
        }
        return false;
    }

    public abstract BDSSkylineSegment toCanonicalSegment(SkylineSegment segment);

}
