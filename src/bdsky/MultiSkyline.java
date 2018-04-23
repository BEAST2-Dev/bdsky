package bdsky;

import beast.core.CalculationNode;
import beast.core.Input;

import java.util.*;

/**
 * A multiskyline made up of simple skylines.
 * This skyline will have a number of segments equal to the union of the number of unique boundaries in the daughter skylines.
 */
public class MultiSkyline extends CalculationNode implements Skyline {

    public Input<List<SimpleSkyline>> skylineInput = new Input<>("skyline", "the simple skylines making up this multiple skyline", new ArrayList<>());

    public MultiSkyline(SimpleSkyline... skyline) {
        if (skyline.length > 0) {
            try {
                initByName("skyline", Arrays.asList(skyline));
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        } else {
            throw new IllegalArgumentException("The multiple skyline cannot be the empty list of simple skylines !");
        }
    }

    @Override
    public void initAndValidate() {}

    @Override
    public List<SkylineSegment> getSegments() {

        List<SimpleSkyline> skylines = skylineInput.get();


        List<Boundary2> boundaries = new ArrayList<>();
        for (int i = 0; i < skylines.size(); i++) {
            Skyline skyline = skylines.get(i);
            List<SkylineSegment> segments = skyline.getSegments();
            for (int j = 0; j < segments.size(); j++) {
                boundaries.add(new Boundary2(j, segments.get(j).t1, skyline, i));
            }
        }
        Collections.sort(boundaries, (o1, o2) -> Double.compare(o1.time, o2.time));

        //System.out.println(boundaries);

        //int[] index = new int[skylines.size()];
        List<SkylineSegment> segments = new ArrayList<>();


        //System.out.println("Boundaries.size = " + boundaries.size());

        int i = 0;
        double start = boundaries.get(0).time;
        while (i < boundaries.size()) {

            int j = i + 1;
            double end = Double.POSITIVE_INFINITY;
            if (j != boundaries.size()) {
                end = boundaries.get(j).time;

                while (j < boundaries.size() && end == start) {
                    j += 1;

                    if (j == boundaries.size()) {
                        end = Double.POSITIVE_INFINITY;
                    } else {
                        end = boundaries.get(j).time;
                    }
                }
                //System.out.println("next end = " + boundaries.get(j));
            }


            double[] value = new double[skylines.size()];
            for (int k = 0; k < skylines.size(); k++) {

                //int ind = index[k];
                //value[k] = skylines.get(k).getValues().get(ind)[0];

                value[k] = skylines.get(k).getValue(start)[0];
            }

            SkylineSegment segment = new SkylineSegment(start, end, value);
            segments.add(segment);
            //System.out.println("Added segment: " + segment);

            // This doesn't work when two skylines have segments ending at the same time
            // i.e. when two skylines have the same set of times
            /*
            if (j != boundaries.size()) {
                index[boundaries.get(j).skylineIndex] += 1;
                System.out.println("incremented index for skyline " + boundaries.get(j).skylineIndex);
            }
            */

            i = j;
            start = end;
        }

        return segments;
    }

    @Override
    public int getDimension() {

        int dim = 0;
        for (SimpleSkyline skyline : skylineInput.get()) {

            dim += skyline.getDimension();
        }
        return dim;
    }

    private class Boundary2 {

        // the index
        int index;

        // time of the boundary
        double time;

        // the skyline the boundary is in
        Skyline skyline;

        int skylineIndex;

        Boundary2(int index, double time, Skyline skyline, int skylineIndex) {
            this.index = index;
            this.time = time;
            this.skyline = skyline;
            this.skylineIndex = skylineIndex;
        }

        public String toString() {
            return "skyline[" + skylineIndex + "].time(" + index + ")=" + time;
        }
    }
}
