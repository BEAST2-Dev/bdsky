package bdsky;

import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.List;
import javax.swing.JFrame;

import beast.core.parameter.RealParameter;
import info.monitorenter.gui.chart.Chart2D;
import info.monitorenter.gui.chart.ITrace2D;
import info.monitorenter.gui.chart.traces.Trace2DSimple;

public class SkylinePlot {

    private SkylinePlot() {
        super();
    }

    public static void addTrace(Skyline skyline, Chart2D chart, Color color) {


        List<SkylineSegment> segments = skyline.getSegments();
        int size = segments.get(0).value.length;

        for (int i = 0; i < size; i++) {

            ITrace2D trace = new Trace2DSimple();
            // Add the trace to the chart. This has to be done before adding points (deadlock prevention):

            chart.addTrace(trace);
            trace.setColor(color);

            for (SkylineSegment segment : segments) {
                trace.addPoint(segment.t1, segment.value[0]);

                if (segment.t2 != Double.POSITIVE_INFINITY) {
                    trace.addPoint(segment.t2, segment.value[0]);
                } else {
                    double extra = 1.0;
                    if (segments.size() > 1) {
                        extra = (segment.t1 - trace.getMinX()) / (segments.size() - 1);
                    }
                    trace.addPoint(segment.t1 + extra, segment.value[0]);
                }
            }
        }
    }

    public static void main(String[] args) throws Exception {

        SimpleSkyline skyline = new SimpleSkyline();
        skyline.setInputValue("times", new RealParameter("0 1 2 3 4"));
        skyline.setInputValue("parameter", new RealParameter("0.5 0.0 3 2 5.5"));
        skyline.setID("skyline1");

        SimpleSkyline skyline2 = new SimpleSkyline();
        skyline2.setInputValue("times", new RealParameter("0 1.1 2.2 3.3 4.4"));
        skyline2.setInputValue("parameter", new RealParameter("1.5 3.2 0.2 5 2.5"));
        skyline2.setID("skyline2");

        // Create a chart:
        Chart2D chart = new Chart2D();
        // Create an ITrace:

        addTrace(skyline, chart, Color.red);
        addTrace(skyline2, chart, Color.blue);
        // Add all points, as it is static:




        // Make it visible:
        // Create a frame.
        JFrame frame = new JFrame("SkylinePlot");
        // add the chart to the frame:
        frame.getContentPane().add(chart);
        frame.setSize(800, 600);
        // Enable the termination button [cross on the upper right edge]:
        frame.addWindowListener(
                new WindowAdapter() {
                    public void windowClosing(WindowEvent e) {
                        System.exit(0);
                    }
                }
        );
        frame.setVisible(true);
    }
}