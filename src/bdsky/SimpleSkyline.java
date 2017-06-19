package bdsky;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A skyline function for a parameter
 *
 * @author Alexei Drummond
 */
public class SimpleSkyline extends CalculationNode implements Skyline {

    // the interval times for the skyline function (e.g. "0 1 2 3")
    public Input<RealParameter> timesInput =
            new Input<>("times", "The times t_i specifying when the parameter changes occur. " +
                    "Times must be in ascending order.", (RealParameter) null);

    // the parameter values, must have the same length as times vector
    // (e.g. "-1 3 4 -1.5", means -1 between [0,1), 3 between [1,2), ..., -1.5 between [3,origin))
    public Input<RealParameter> parameterInput =
            new Input<>("parameter",
                    "The parameter values specifying the value for each piecewise constant segment of the skyline function. " +
                    "The first value is between t_0 and t_1, the last value is between t_n and infinity. " +
                    "Should be the same length as time vector", (RealParameter) null, Input.Validate.REQUIRED);

    public Input<RealParameter> originInput =
            new Input<>("origin",
                    "The 'origin' of the skyline. Either the beginning or the end, depending on the direction of time.", (RealParameter) null, Input.Validate.REQUIRED);

    public Input<Boolean> timesRelativeToOriginInput =
            new Input<>("timesRelativeToOrigin",
                    "If true, then the times are expressed at relative to the origin. default is false.", (Boolean) false, Input.Validate.OPTIONAL);


    @Override
    public void initAndValidate() {

        Double[] times = getTimes();

        double smallest = Double.NEGATIVE_INFINITY;
        for (Double time : times) {
            if (time < smallest) {
                throw new RuntimeException("Times must be in ascending order!");
            }
            smallest = time;
        }

        if (times[0] != 0) {
            throw new IllegalArgumentException("Skyline times must start with 0!");
        }
    }

    /**
     * Set the bounds for the skyline parameter
     */
    public void setBounds(Double lower, Double upper) {
        parameterInput.get().setBounds(upper, lower);
    }

    /**
     * Set lower bound for skyline parameter
     */
    public void setLower(Double lower) {
        parameterInput.get().setLower(lower);
    }

    /**
     * Set upper bound for skyline parameter
     */
    public void setUpper(Double upper) {
        parameterInput.get().setUpper(upper);
    }

    /**
     * Get lower bound for skyline parameter
     */

    public Double getLower() {
        return parameterInput.get().getLower();
    }

    /**
     * Get upper bound for skyline parameter
     */
    public Double getUpper() {
        return parameterInput.get().getUpper();
    }


    /**
     *
     * @return the times for this skyline function. Will pre-calculate absolute times using the origin parameter if the times input parameter provides relative times.
     */
    public Double[] getTimes() {
        Double[] times = timesInput.get().getValues();
        if (timesRelativeToOriginInput.get()) {
            for (int i = 0; i < times.length; i++) {
                times[i] *= originInput.get().getValue();
            }
        }
        return times;
    }

    /**
     *
     * @return the values for this skyline function
     */
    public List<double[]> getValues() {

        List<double[]> values = new ArrayList<>();
        for (double val : parameterInput.get().getValues()) {
            values.add(new double[] {val});
        }

        return values;
    }

    private Double[] rawValues() {
        return parameterInput.get().getValues();
    }

    public double[] getValue(double time) {
        Double[] times = getTimes();

        if (time < times[0]) {
            throw new RuntimeException("Time is smaller than smallest time in skyline function!");
        }

        int index = Arrays.binarySearch(times,time);

        Double[] values = rawValues();
        if (index < 0) {
            //returns (-(insertion point) - 1)
            int insertionPoint = -(index+1);
            return new double[] {values[insertionPoint-1]};
        } else {
            return new double[] {values[index]};
        }
    }

    /**
     * @param time1
     * @param time2
     * @return the segments of the skyline plot between the two times.
     */
    public List<SkylineSegment> getSegments(double time1, double time2) {

        Double[] times = getTimes();
        int insertionPoint;

        if (time1 < times[0] || time2 < times[0]) {
            throw new RuntimeException("Time is smaller than smallest time in skyline function!");
        }

        if (time1 > time2 || time1 == time2) {
            throw new RuntimeException("time1 must be smaller than time2!");
        }

        List<SkylineSegment> segments = new ArrayList<>();
        Double[] rawValues = rawValues();

        // Only one interval in skyline
        if (times.length == 1) {
            segments.add(new SkylineSegment(time1,time2,rawValues[0]));
            return segments;
        }

        int index1 = Arrays.binarySearch(times, time1);
        int index2 = Arrays.binarySearch(times, time2);


        // Get the first segment
        if (index1 == index2) {
            // Case 1: Both times fall within one interval and NOT on a boundary
            // (Since time1 != time2 this can only happen when the times are not in the skyline,
            //  which means both indices will be negative)
            insertionPoint = -(index1 + 1);
            segments.add(new SkylineSegment(time1, time2, rawValues[insertionPoint-1]));
        } else
        if (index1 < 0) {
            // Case 2: time1 is in the middle of a segment
            insertionPoint = -(index1 + 1);
            segments.add(new SkylineSegment(time1,times[insertionPoint], rawValues[insertionPoint-1]));
            index1 = insertionPoint;
        } else {
            // Case 3: time1 is at the start of a segment
            segments.add(new SkylineSegment(times[index1],times[index1+1], rawValues[index1]));
            index1 += 1;
        }

        // If there are more segments, add them
        if (index1 != index2) {
            // Does time2 end in the middle of an interval?
            insertionPoint = index2 < 0 ? -(index2 + 1) : index2+1;

            // Add all complete intervals between the first and time2
            for (int i = index1; i < insertionPoint-1; i++ ) {
                segments.add(new SkylineSegment(times[i],times[i+1], rawValues[i]));
            }

            // If time2 is in the middle of an interval add the last segment
            if (index2 < 0) {
                segments.add(new SkylineSegment(times[insertionPoint-1],time2, rawValues[insertionPoint-1]));
            }
        }

        return segments;
    }


    @Override
    public List<SkylineSegment> getSegments() {
        return getSegments(0, Double.POSITIVE_INFINITY);
    }

    @Override
    public int getDimension() {
        return 1;
    }
}
