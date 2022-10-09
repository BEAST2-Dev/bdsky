package test.beast.app.beauti;

import java.io.File;

import org.junit.Test;
import org.testfx.api.FxRobot;

import test.beastfx.app.beauti.BeautiBase;

public class BeautiBDSKYTest extends BeautiBase {

	@Test
	public void simpleTest(FxRobot robot) throws Exception {
		
		importAlignment("../beast2/examples/nexus", new File("anolis.nex"));
		robot.clickOn("#Priors");
		
		warning("Change Tree prior to Birth Death Skyline");
		selectFromCombobox(robot, "TreeDistribution", "Birth Death Skyline Contemporary");
		printBeautiState();

		warning("Change Tree prior to Birth Death Skyline");
		selectFromCombobox(robot, "TreeDistribution", "Birth Death Skyline Serial");
		printBeautiState();

		makeSureXMLParses();
	}
}
