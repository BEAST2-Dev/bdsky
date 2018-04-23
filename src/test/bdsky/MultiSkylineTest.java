package test.bdsky;

import bdsky.MultiSkyline;
import bdsky.SimpleSkyline;
import beast.core.parameter.RealParameter;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;

/**
 * Tests for simple skyline
 */
public class MultiSkylineTest {

    @Test
    public void testGetValue() throws Exception {

        SimpleSkyline skyline = new SimpleSkyline();
        skyline.setInputValue("times", new RealParameter("0 1 2 3"));
        skyline.setInputValue("parameter", new RealParameter("-1 3 4 -1.5"));

        SimpleSkyline skyline2 = new SimpleSkyline();
        skyline2.setInputValue("times", new RealParameter("0 1.5 2.5 3.5"));
        skyline2.setInputValue("parameter", new RealParameter("2 1.5 -2.7 4.5"));

        List<SimpleSkyline> simpleSkylines = new ArrayList<>();
        simpleSkylines.add(skyline);
        simpleSkylines.add(skyline2);

        MultiSkyline multiSkyline = new MultiSkyline(simpleSkylines.toArray(new SimpleSkyline[0]));
//        multiSkyline.skylineInput.setValue(simpleSkylines,multiSkyline);
//        multiSkyline.initAndValidate();

        assertEquals(2, multiSkyline.getDimension());

        assertEquals(7, multiSkyline.getSegments().size());

        assertEquals(-1, multiSkyline.getValue(0.75)[0], 0);
        assertEquals(2, multiSkyline.getValue(0.75)[1], 0);

        assertEquals(3, multiSkyline.getValue(1.25)[0], 0);
        assertEquals(2, multiSkyline.getValue(1.25)[1], 0);

        assertEquals(3, multiSkyline.getValue(1.75)[0], 0);
        assertEquals(1.5, multiSkyline.getValue(1.75)[1], 0);

        assertEquals(4, multiSkyline.getValue(2.25)[0], 0);
        assertEquals(1.5, multiSkyline.getValue(2.25)[1], 0);

        assertEquals(4, multiSkyline.getValue(2.75)[0], 0);
        assertEquals(-2.7, multiSkyline.getValue(2.75)[1], 0);

        assertEquals(-1.5, multiSkyline.getValue(3.25)[0], 0);
        assertEquals(-2.7, multiSkyline.getValue(3.25)[1], 0);

        assertEquals(-1.5, multiSkyline.getValue(3.75)[0], 0);
        assertEquals(4.5, multiSkyline.getValue(3.75)[1], 0);
    }

}