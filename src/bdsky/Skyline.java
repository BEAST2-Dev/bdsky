package bdsky;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A multivariate skyline interface
 *
 * @author Alexei Drummond
 */
public interface Skyline {

    /**
     * @return a list of segments in increasing time order.
     */
    List<SkylineSegment> getSegments();

    /**
     * @param time the time of interest
     * @return the value of this skyline at the given time
     */
    default double[] getValue(double time) {

        Double[] times = getTimes();

        if (time < times[0]) {
            throw new RuntimeException("Time is smaller than smallest time in skyline function!");
        }

        int index = Arrays.binarySearch(times,time);

        List<double[]> values = getValues();
        if (index < 0) {
            //returns (-(insertion point) - 1)
            int insertionPoint = -(index+1);
            return values.get(insertionPoint-1);
        } else {
            return values.get(index);
        }
    }

    /**
     * @return the absolute start times of the segments in index order.
     */
    default Double[] getTimes() {

        List<SkylineSegment> segments = getSegments();

        Double[] times = new Double[segments.size()];
        for (int i = 0; i < times.length; i++) {
            times[i] = segments.get(i).t1;
        }

        if (segments.get(segments.size()-1).t2 < Double.POSITIVE_INFINITY) {
            throw new RuntimeException("Last segment should extend to positive infinity!");
        }

        return times;
    }

    /**
     * @return the values of the segments in index order.
     */
    default List<double[]> getValues() {

        List<SkylineSegment> segments = getSegments();

        List<double[]> values = new ArrayList<>();
        for (int i = 0; i < segments.size(); i++) {
            values.add(segments.get(i).value);
        }
        return values;
    }

    /**
     * This is not the number of segments, but the dimension of each segment.
     * @return the dimension of the parameter in this skyline.
     */
    int getDimension();
}
