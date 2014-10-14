package beast.app.beauti;

import beast.app.beauti.BeautiDoc;
import beast.app.draw.InputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.evolution.speciation.BirthDeathSkylineModel;

import javax.swing.*;
import javax.swing.border.Border;
import java.awt.*;

/**
 * Created with IntelliJ IDEA.
 * User: Denise
 * Date: 15.09.14
 * Time: 17:29
 */
public class BirthDeathSkylineInputEditor extends InputEditor.Base {

    @Override
    public Class<?> type() {
        return BirthDeathSkylineModel.class;
    }

    public BirthDeathSkylineInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?>[] types() {
        return new Class<?>[0];  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {
        JLabel label = new JLabel("Greetings from Basel");
        add(label);
    }



}
