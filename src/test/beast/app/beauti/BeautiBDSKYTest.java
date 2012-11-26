package test.beast.app.beauti;

import java.io.File;

import org.fest.swing.fixture.JTabbedPaneFixture;
import org.junit.Test;

public class BeautiBDSKYTest extends BeautiBase {

	@Test
	public void simpleTest() throws Exception {
		
		importAlignment("../beast2/examples/nexus", new File("anolis.nex"));
		JTabbedPaneFixture f = beautiFrame.tabbedPane();
		f.selectTab("Priors");
		
		warning("Change Tree prior to Birth Death Skyline");
		beautiFrame.comboBox("TreeDistribution").selectItem("Birth Death Skyline Contemporary");
		printBeautiState(f);

		warning("Change Tree prior to Birth Death Skyline");
		beautiFrame.comboBox("TreeDistribution").selectItem("Birth Death Skyline Serial");
		printBeautiState(f);

		makeSureXMLParses();
	}
}
