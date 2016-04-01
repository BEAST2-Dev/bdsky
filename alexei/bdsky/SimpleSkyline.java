package bdsky;

import bdsky.Skyline;
import bdsky.SkylineSegment;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A skyline function for a parameter
 */
public class SimpleSkyline extends CalculationNode implements Skyline {

    // the interval times for the skyline function (e.g. "0 1 2 3")
    public Input<RealParameter> timesInput =
            new Input<RealParameter>("times", "The times t_i specifying when the parameter changes occur. " +
                    "Times must be in ascending order.", (RealParameter) null);

    // the parameter values, must have the same length as times vector
    // (e.g. "-1 3 4 -1.5", means -1 between [0,1), 3 between [1,2), ..., -1.5 between [3,infinity))
    public Input<RealParameter> parameterInput =
            new Input<RealParameter>("parameter",
                    "The parameter values specifying the value for each piecewise constant segment of the skyline function. " +
                    "The first value is between t_0 and t_1, the last value is between t_n and infinity. " +
                    "Should be the same length as time vector", (RealParameter) null);

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
    }

    /**
     *
     * @return the times for this skyline function
     */
    public Double[] getTimes() {
        return timesInput.get().getValues();
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

        if (time1 < times[0] || time2 < times[0]) {
            throw new RuntimeException("Time is smaller than smallest time in skyline function!");
        }

        if (time1 > time2 || time1 == time2) {
            throw new RuntimeException("time1 must be smaller than time2!");
        }

        List<SkylineSegment> segments = new ArrayList<>();

        int index1 = Arrays.binarySearch(times, time1);
        int index2 = Arrays.binarySearch(times, time2);

        Double[] rawValues = rawValues();

        // same insertion point
        if (index1 == index2) {
            int insertionPoint = -(index1 + 1);
            segments.add(new SkylineSegment(time1, time2, rawValues[insertionPoint-1]));
            return segments;
        }

        // not same insertion points
        if (index1 < 0) {
            int insertionPoint = -(index1 + 1);
            segments.add(new SkylineSegment(time1,times[insertionPoint], rawValues[insertionPoint-1]));
            index1 = insertionPoint;
            if (index1 == index2) return segments;
        } else {
            segments.add(new SkylineSegment(times[index1],times[index1+1], rawValues[index1]));
            index1 += 1;
            if (index1 == index2) return segments;
        }
        if (index2 < 0) {
            int insertionPoint = -(index2 + 1);
            for (int i = index1; i < insertionPoint-1; i++ ) {
                segments.add(new SkylineSegment(times[i],times[i+1], rawValues[i]));
            }
            segments.add(new SkylineSegment(times[insertionPoint-1],time2, rawValues[insertionPoint-1]));
            return segments;
        } else {
            for (int i = index1; i < index2; i++ ) {
                segments.add(new SkylineSegment(times[i],times[i+1], rawValues[i]));
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
