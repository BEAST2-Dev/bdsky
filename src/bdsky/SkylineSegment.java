package bdsky;

import java.util.Arrays;

/**
 * A piecewise constant segment of a skyline.
 *
 * TODO next and prev are never set or used
 * (should be set by SimpleSkyline.getSegments() and MultiSkyline.getSegments())
 */
public class SkylineSegment {

    // the start time of this segment
    public double t1;

    // the end time of this segment
    public double t2;

    // the parameter values
    public double[] value;

    SkylineSegment next, prev = null;

    public SkylineSegment(double start, double end, double value) {
        this.t1 = start;
        this.t2 = end;
        this.value = new double[] {value};
    }

    public SkylineSegment(double start, double end, double[] value) {
        this.t1 = start;
        this.t2 = end;
        this.value = value;
    }

    void setNextSegment(SkylineSegment next) {
        if (next.t1 == t2) {
            this.next = next;
        } else {
            throw new RuntimeException("next.t1 must equal this.t2!");
        }
        if (next.prev != this) {
            next.setPreviousSegment(this);
        }
    }

    void setPreviousSegment(SkylineSegment prev) {
        if (prev.t2 == t1) {
            this.prev = prev;
        } else {
            throw new RuntimeException("prev.t2 must equal this.t1!");
        }
        if (prev.next != this) {
            prev.setNextSegment(this);
        }
    }

    /**
     * @return the start time of the segment.
     */
    public final double start() { return t1; }

    /**
     * @return the end time of the segment.
     */
    public final double end() { return t2; }

    public String toString() {
        return "segment(" + t1 + "," + t2 + ") = " + Arrays.toString(value);
    }


}
