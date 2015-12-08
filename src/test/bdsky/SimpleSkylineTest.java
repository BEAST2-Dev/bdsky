package test.bdsky;

import bdsky.SkylineSegment;
import beast.core.parameter.RealParameter;
import bdsky.SimpleSkyline;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.*;

/**
 * Tests for simple skyline
 */
public class SimpleSkylineTest {

    @Test
    public void testGetValue() throws Exception {

        SimpleSkyline skyline = new SimpleSkyline();
        skyline.setInputValue("times", new RealParameter("0 1 2 3"));
        skyline.setInputValue("parameter", new RealParameter("-1 3 4 -1.5"));

        assertEquals(-1, skyline.getValue(0)[0], 0);
        assertEquals(-1, skyline.getValue(0.5)[0], 0);

        assertEquals(3, skyline.getValue(1)[0], 0);
        assertEquals(3, skyline.getValue(1.5)[0], 0);

        assertEquals(4, skyline.getValue(2)[0], 0);
        assertEquals(4, skyline.getValue(2.5)[0], 0);

        assertEquals(-1.5, skyline.getValue(3)[0], 0);
        assertEquals(-1.5, skyline.getValue(3.5)[0], 0);
    }

    @Test
    public void testGetSegments() throws Exception {

        SimpleSkyline skyline = new SimpleSkyline();
        skyline.setInputValue("times", new RealParameter("0 1 2 3"));
        skyline.setInputValue("parameter", new RealParameter("-1 3 4 -1.5"));

        List<SkylineSegment> segments = skyline.getSegments();

        assertEquals("Checking number of segments", 4, segments.size());

        assertEquals(0, segments.get(0).t1, 0.0);
        assertEquals(1, segments.get(1).t1, 0.0);
        assertEquals(2, segments.get(2).t1, 0.0);
        assertEquals(3, segments.get(3).t1, 0.0);

        assertEquals(1, segments.get(0).t2, 0.0);
        assertEquals(2, segments.get(1).t2, 0.0);
        assertEquals(3, segments.get(2).t2, 0.0);
        assertEquals(Double.POSITIVE_INFINITY, segments.get(3).t2, 0.0);

        assertEquals(-1, segments.get(0).value[0], 0.0);
        assertEquals(3, segments.get(1).value[0], 0.0);
        assertEquals(4, segments.get(2).value[0], 0.0);
        assertEquals(-1.5, segments.get(3).value[0], 0.0);
    }

    @Test
    public void testGetSegments1() throws Exception {

        SimpleSkyline skyline = new SimpleSkyline();
        skyline.setInputValue("times", new RealParameter("0 1 2 3"));
        skyline.setInputValue("parameter", new RealParameter("-1 3 4 -1.5"));

        List<SkylineSegment> segments = skyline.getSegments(0,1);

        assertEquals("Checking number of segments", 1, segments.size());

        assertEquals(0, segments.get(0).t1, 0.0);

        assertEquals(1, segments.get(0).t2, 0.0);

        assertEquals(-1, segments.get(0).value[0], 0.0);

    }

    @Test
    public void testGetSegments2() throws Exception {

        SimpleSkyline skyline = new SimpleSkyline();
        skyline.setInputValue("times", new RealParameter("0 1 2 3"));
        skyline.setInputValue("parameter", new RealParameter("-1 3 4 -1.5"));

        List<SkylineSegment> segments = skyline.getSegments(0.2,0.3);

        assertEquals("Checking number of segments", 1, segments.size());

        assertEquals(0.2, segments.get(0).t1, 0.0);

        assertEquals(0.3, segments.get(0).t2, 0.0);

        assertEquals(-1, segments.get(0).value[0], 0.0);

    }

    @Test
    public void testGetSegments3() throws Exception {

        SimpleSkyline skyline = new SimpleSkyline();
        skyline.setInputValue("times", new RealParameter("0 1 2 3"));
        skyline.setInputValue("parameter", new RealParameter("-1 3 4 -1.5"));

        List<SkylineSegment> segments = skyline.getSegments(0.2,1.3);

        assertEquals("Checking number of segments", 2, segments.size());

        assertEquals(0.2, segments.get(0).t1, 0.0);
        assertEquals(1, segments.get(1).t1, 0.0);

        assertEquals(1, segments.get(0).t2, 0.0);
        assertEquals(1.3, segments.get(1).t2, 0.0);

        assertEquals(-1, segments.get(0).value[0], 0.0);
        assertEquals(3, segments.get(1).value[0], 0.0);
    }

    @Test
    public void testGetSegments4() throws Exception {

        SimpleSkyline skyline = new SimpleSkyline();
        skyline.setInputValue("times", new RealParameter("0 1 2 3"));
        skyline.setInputValue("parameter", new RealParameter("-1 3 4 -1.5"));

        List<SkylineSegment> segments = skyline.getSegments(0.2,2.3);

        assertEquals("Checking number of segments", 3, segments.size());

        assertEquals(0.2, segments.get(0).t1, 0.0);
        assertEquals(1, segments.get(1).t1, 0.0);
        assertEquals(2, segments.get(2).t1, 0.0);

        assertEquals(1, segments.get(0).t2, 0.0);
        assertEquals(2, segments.get(1).t2, 0.0);
        assertEquals(2.3, segments.get(2).t2, 0.0);

        assertEquals(-1, segments.get(0).value[0], 0.0);
        assertEquals(3, segments.get(1).value[0], 0.0);
        assertEquals(4, segments.get(2).value[0], 0.0);
    }

    @Test
    public void testGetSegments5() throws Exception {

        SimpleSkyline skyline = new SimpleSkyline();
        skyline.setInputValue("times", new RealParameter("0 1 2 3"));
        skyline.setInputValue("parameter", new RealParameter("-1 3 4 -1.5"));

        List<SkylineSegment> segments = skyline.getSegments(1.1,5.3);

        assertEquals("Checking number of segments", 3, segments.size());

        assertEquals(1.1, segments.get(0).t1, 0.0);
        assertEquals(2, segments.get(1).t1, 0.0);
        assertEquals(3, segments.get(2).t1, 0.0);

        assertEquals(2, segments.get(0).t2, 0.0);
        assertEquals(3, segments.get(1).t2, 0.0);
        assertEquals(5.3, segments.get(2).t2, 0.0);

        assertEquals(3, segments.get(0).value[0], 0.0);
        assertEquals(4, segments.get(1).value[0], 0.0);
        assertEquals(-1.5, segments.get(2).value[0], 0.0);
    }

    @Test
    public void testGetSegments6() throws Exception {

        SimpleSkyline skyline = new SimpleSkyline();
        skyline.setInputValue("times", new RealParameter("0 1 2 3"));
        skyline.setInputValue("parameter", new RealParameter("-1 3 4 -1.5"));

        List<SkylineSegment> segments = skyline.getSegments(1,3);

        assertEquals("Checking number of segments", 2, segments.size());

        assertEquals(1, segments.get(0).t1, 0.0);
        assertEquals(2, segments.get(1).t1, 0.0);

        assertEquals(2, segments.get(0).t2, 0.0);
        assertEquals(3, segments.get(1).t2, 0.0);

        assertEquals(3, segments.get(0).value[0], 0.0);
        assertEquals(4, segments.get(1).value[0], 0.0);
    }


}