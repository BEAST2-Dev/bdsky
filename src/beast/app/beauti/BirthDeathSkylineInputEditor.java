//package beast.app.beauti;
//
//import beastfx.app.inputeditor.InputEditor;
//import beast.base.core.BEASTInterface;
//import beast.base.core.Input;
//import beast.evolution.speciation.BirthDeathSkylineModel;
//
//import javax.swing.*;
//
///**
// * User: Denise
// * Date: 15.09.14
// * Time: 17:29
// */
//public class BirthDeathSkylineInputEditor extends InputEditor.Base {
//
//    @Override
//    public Class<?> type() {
//        return BirthDeathSkylineModel.class;
//    }
//
//    public BirthDeathSkylineInputEditor(BeautiDoc doc) {
//        super(doc);
//    }
//
//    @Override
//    public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {
//        JLabel label = new JLabel("Greetings from Basel " + input.get().toString());
//        add(label);
//    }
//
//
//
//}
