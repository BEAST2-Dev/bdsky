package test.beast.evolution.tree.coalescent;

import junit.framework.TestCase;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.evolution.tree.Tree;
import beast.evolution.tree.Node;
import beast.evolution.speciation.BirthDeathSkylineModel;
import beast.core.parameter.RealParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.Description;
import beast.util.TreeParser;
import org.junit.Test;

import java.io.PrintStream;

/**
 * @author Denise Kuehnert
 * Date: Mar 29, 2011
 */

@Description("Test the BirthDeathSerialSkylineModel with a small tree example")

public class BirthDeathSkylineTest extends TestCase {


    @Test
    public void testRhoSasha() throws Exception {


        Tree tree = new TreeParser("(((((t1:0.4595008531,t25:0.4595008531):0.3373053072,t23:0.3567584538):0.007310819036,t16:0.3489190732):0.331009529,((t18:0.03315384045,t14:0.03315384045):0.5063451374,(t10:0.4211543131,t15:0.4211543131):0.1183446648):0.5956275305):0.1158090878,((t19:0.9429393194,((t6:0.363527235,t11:0.4417423167):0.01881829549,((((t3:0.3071904376,(((t24:0.01065209364,t13:0.01065209364):0.06076485145,t8:0.07141694509):0.123620245,(t22:0.1616119808,t2:0.1616119808):0.03342520927):0.1121532475):0.24520579,t9:0.5523962276):0.3852615426,(((t20:0.2935970782,(t17:0.06569090089,t4:0.06569090089):0.2279061773):0.08350780408,(t21:0.05109047139,t5:0.05109047139):0.3260144109):0.2298344132,t7:0.6069392955):0.3307184747):0.01206284377,t26:0.9497206139):0.05755333197):0.03290891884):0.07263755325,t12:1.112820418):0.1381151782);",false);

        PrintStream treeString = new PrintStream("out.tree");
        tree.log(1, treeString);


        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("2."));
        bdssm.setInputValue("conditionOnSurvival", false);

        bdssm.setInputValue("birthRate", new RealParameter("3."));
        bdssm.setInputValue("deathRate", new RealParameter("2.5"));
        bdssm.setInputValue("samplingRate", new RealParameter("2."));
        bdssm.setInputValue("rho", new RealParameter("0.0 0.05 0.01"));
        bdssm.setInputValue("rhoSamplingTimes","0. 1. 1.5");
        bdssm.setInputValue("reverseTimeArrays","false false false false");
        bdssm.initAndValidate();
        assertEquals(-124.96086690757612, bdssm.calculateTreeLogLikelihood(tree), 1e-4);     // this result is from BEAST, not double checked in R


        // now the same with reverse rhoSamplingTimes
        bdssm =  new BirthDeathSkylineModel();
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("2."));
        bdssm.setInputValue("conditionOnSurvival", false);

        bdssm.setInputValue("birthRate", new RealParameter("3."));
        bdssm.setInputValue("deathRate", new RealParameter("2.5"));
        bdssm.setInputValue("samplingRate", new RealParameter("2."));
        bdssm.setInputValue("rho", new RealParameter("0.0 0.05 0.01"));
        bdssm.setInputValue("rhoSamplingTimes","0. 0.5 1.");
        bdssm.setInputValue("reverseTimeArrays","false false false true");
        bdssm.initAndValidate();
        assertEquals(-124.96086690757612, bdssm.calculateTreeLogLikelihood(tree), 1e-4);     // this result is from BEAST, not double checked in R

    }


    @Test
    public void test3intsmitRho() throws Exception {


        Tree tree = new TreeParser("(((((t1:0.4595008531,t25:0.4595008531):0.3373053072,t23:0.3567584538):0.007310819036,t16:0.3489190732):0.331009529,((t18:0.03315384045,t14:0.03315384045):0.5063451374,(t10:0.4211543131,t15:0.4211543131):0.1183446648):0.5956275305):0.1158090878,((t19:0.9429393194,((t6:0.363527235,t11:0.4417423167):0.01881829549,((((t3:0.3071904376,(((t24:0.01065209364,t13:0.01065209364):0.06076485145,t8:0.07141694509):0.123620245,(t22:0.1616119808,t2:0.1616119808):0.03342520927):0.1121532475):0.24520579,t9:0.5523962276):0.3852615426,(((t20:0.2935970782,(t17:0.06569090089,t4:0.06569090089):0.2279061773):0.08350780408,(t21:0.05109047139,t5:0.05109047139):0.3260144109):0.2298344132,t7:0.6069392955):0.3307184747):0.01206284377,t26:0.9497206139):0.05755333197):0.03290891884):0.07263755325,t12:1.112820418):0.1381151782);",false);

        PrintStream treeString = new PrintStream("out.tree");
        tree.log(1, treeString);


        for (int i = 0; i<3; i++){


            switch (i){
                case 0:{
                    BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();
                    bdssm.setInputValue("tree", tree);
                    bdssm.setInputValue("origin", new RealParameter("2."));
                    bdssm.setInputValue("conditionOnSurvival", false);

                    bdssm.setInputValue("birthRate", new RealParameter("3. 2. 4."));
                    bdssm.setInputValue("deathRate", new RealParameter("2.5 1. .5"));
                    bdssm.setInputValue("samplingRate", new RealParameter("2. 0.5 1."));
                    bdssm.setInputValue("rho", new RealParameter("0.01"));
                    bdssm.setInputValue("birthRateChangeTimes", new RealParameter("0. 1. 1.5"));
                    bdssm.setInputValue("deathRateChangeTimes", new RealParameter("0. 1. 1.5"));
                    bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("0. 1. 1.5"));

                    bdssm.initAndValidate();
                    bdssm.printTempResults = true;

                    assertEquals(-87.59718586549747, bdssm.calculateTreeLogLikelihood(tree), 1e-4);

                }
                case 1:{
                    BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();
                    bdssm.setInputValue("tree", tree);
                    bdssm.setInputValue("origin", new RealParameter("2."));
                    bdssm.setInputValue("conditionOnSurvival", false);

                    bdssm.setInputValue("birthRate", new RealParameter("3. 2. 4."));
                    bdssm.setInputValue("deathRate", new RealParameter("2.5 1. .5"));
                    bdssm.setInputValue("samplingRate", new RealParameter("2. 0.5 1."));
                    bdssm.setInputValue("birthRateChangeTimes", new RealParameter("0. 1. 1.5"));
                    bdssm.setInputValue("deathRateChangeTimes", new RealParameter("0. 1. 1.5"));
                    bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("0. 1. 1.5"));

                    bdssm.setInputValue("rho", new RealParameter("0.05 0.01"));
                    bdssm.setInputValue("rhoSamplingTimes","0. 1.");
                    bdssm.initAndValidate();
                    assertEquals(-87.96488, bdssm.calculateTreeLogLikelihood(tree), 1e-4);
                }

                case 2:{
                    BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();
                    bdssm.setInputValue("tree", tree);
                    bdssm.setInputValue("origin", new RealParameter("2."));
                    bdssm.setInputValue("conditionOnSurvival", false);

                    bdssm.setInputValue("birthRate", new RealParameter("3. 2. 4. 4."));
                    bdssm.setInputValue("deathRate", new RealParameter("2.5 1. .5 0.5"));
                    bdssm.setInputValue("samplingRate", new RealParameter("2. 0.5 1. 2.0"));
                    bdssm.setInputValue("rho", new RealParameter("0.05 0.01"));
                    bdssm.setInputValue("rhoSamplingTimes","0. 1.");
                    bdssm.setInputValue("intervalTimes", new RealParameter("0. 0.5 1. 1.5"));
                    bdssm.initAndValidate();
                    assertEquals(-106.06555718977357, bdssm.calculateTreeLogLikelihood(tree), 1e-4);
                }
            }
        }
    }


    @Test
    public void test1int() throws Exception{
        Tree tree = new TreeParser("(((((t1:0.4595008531,t25:0.4595008531):0.3373053072,t23:0.3567584538):0.007310819036,t16:0.3489190732):0.331009529,((t18:0.03315384045,t14:0.03315384045):0.5063451374,(t10:0.4211543131,t15:0.4211543131):0.1183446648):0.5956275305):0.1158090878,((t19:0.9429393194,((t6:0.363527235,t11:0.4417423167):0.01881829549,((((t3:0.3071904376,(((t24:0.01065209364,t13:0.01065209364):0.06076485145,t8:0.07141694509):0.123620245,(t22:0.1616119808,t2:0.1616119808):0.03342520927):0.1121532475):0.24520579,t9:0.5523962276):0.3852615426,(((t20:0.2935970782,(t17:0.06569090089,t4:0.06569090089):0.2279061773):0.08350780408,(t21:0.05109047139,t5:0.05109047139):0.3260144109):0.2298344132,t7:0.6069392955):0.3307184747):0.01206284377,t26:0.9497206139):0.05755333197):0.03290891884):0.07263755325,t12:1.112820418):0.1381151782);", false);


        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", "2.");
        bdssm.setInputValue("conditionOnSurvival", true);
        bdssm.setInputValue("birthRate", new RealParameter("2."));
        bdssm.setInputValue("deathRate", new RealParameter("0.5"));
        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));

        bdssm.setInputValue("rho", new RealParameter("1."));
        bdssm.initAndValidate();
//        System.out.println("\na) Likelihood: " + bdssm.calculateTreeLogLikelihood(tree));
        assertEquals(-21.42666177086957, bdssm.calculateTreeLogLikelihood(tree), 1e-5);

    }

    @Test
    public void testRhoA() throws Exception{

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", true);
        //        bdssm.setInputValue("birthRate", new RealParameter("2."));
        //        bdssm.setInputValue("deathRate", new RealParameter("1."));
        //        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));
        bdssm.setInputValue("R0", new RealParameter(new Double[]{4./3.}));
        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{1./3.}));
        //bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("birthRateChangeTimes", new RealParameter("0."));
        bdssm.setInputValue("deathRateChangeTimes", new RealParameter("0."));
        bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("0."));


//        // a)
        bdssm.setInputValue("rho", new RealParameter("0.01"));
        bdssm.initAndValidate();
//        System.out.println("\na) Likelihood: " + bdssm.calculateTreeLogLikelihood(tree));
        assertEquals(-22.54937791357737, bdssm.calculateTreeLogLikelihood(tree), 1e-5);
//-22.54937791357737

    }

//    @Test
//    public void testRhoB() throws Exception{
//
//        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();
//
//        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
//        bdssm.setInputValue("tree", tree);
//        bdssm.setInputValue("orig_root", new RealParameter("1."));
//        bdssm.setInputValue("conditionOnSurvival", true);
//        //        bdssm.setInputValue("birthRate", new RealParameter("2."));
//        //        bdssm.setInputValue("deathRate", new RealParameter("1."));
//        //        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));
//        bdssm.setInputValue("R0", new RealParameter(new Double[]{4./3.}));
//        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
//        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{1./3.}));
//        //bdssm.setInputValue("intervalNumber", 1);
//        bdssm.setInputValue("birthRateChangeTimes", new RealParameter("0."));
//        bdssm.setInputValue("deathRateChangeTimes", new RealParameter("0."));
//        bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("0."));
//
//        // b)
//        bdssm.setInputValue("intervalNumber", 3);
//
//        bdssm.setInputValue("R0", new RealParameter(new Double[]{4./3., 4./3., 4./3.}));
//        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter(new Double[]{1.5,1.5,1.5}));
//        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{1./3., 1./3., 1./3.}));
//
//        bdssm.setInputValue("birthRateChangeTimes", new RealParameter("0. 2.5 3.5"));
//        bdssm.setInputValue("deathRateChangeTimes", new RealParameter("0. 2.5 3.5"));
//        bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("0. 2.5 3.5"));
//        bdssm.setInputValue("rhoSamplingTimes", new RealParameter("2.5 3.5 6.0"));
//
////        bdssm.setInputValue("rho", new RealParameter("0.01 0.02 0.05"));
//        bdssm.setInputValue("rho", new RealParameter("0.05 0.02 0.01"));
//        bdssm.initAndValidate();
////        System.out.println("\nb) Likelihood: " + bdssm.calculateTreeLogLikelihood(tree));
//        assertEquals(-28.19069830850966, bdssm.calculateTreeLogLikelihood(tree), 1e-5);
//        //-28.19069830850966
//
//    }

    @Test
    public void testRhoC1() throws Exception{

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", true);
        //        bdssm.setInputValue("birthRate", new RealParameter("2."));
        //        bdssm.setInputValue("deathRate", new RealParameter("1."));
        //        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));
        bdssm.setInputValue("R0", new RealParameter(new Double[]{4./3.}));
        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{1./3.}));
        //bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("birthRateChangeTimes", new RealParameter("0."));
        bdssm.setInputValue("deathRateChangeTimes", new RealParameter("0."));
        bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("0."));

        // c1)
        bdssm.setInputValue("R0", new RealParameter(new Double[]{4./3.}));
        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter(new Double[]{1.5}));
        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{1./3.}));
//        bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("intervalTimes", new RealParameter("0."));
        bdssm.setInputValue("birthRateChangeTimes", null);
        bdssm.setInputValue("deathRateChangeTimes", null);
        bdssm.setInputValue("samplingRateChangeTimes", null);
        bdssm.setInputValue("rho", new RealParameter("0.1"));
        bdssm.initAndValidate();
//        System.out.println("\nc) Likelihood: " + bdssm.calculateTreeLogLikelihood(tree));
        assertEquals(-20.74594782518312, bdssm.calculateTreeLogLikelihood(tree), 1e-5);

    }

//    @Test
//    public void testRhoC3() throws Exception{
//
//        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();
//
//        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
//        bdssm.setInputValue("tree", tree);
//        bdssm.setInputValue("orig_root", new RealParameter("1."));
//        bdssm.setInputValue("conditionOnSurvival", true);
//        //        bdssm.setInputValue("birthRate", new RealParameter("2."));
//        //        bdssm.setInputValue("deathRate", new RealParameter("1."));
//        //        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));
//        bdssm.setInputValue("R0", new RealParameter(new Double[]{4./3.}));
//        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
//        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{1./3.}));
//        //bdssm.setInputValue("intervalNumber", 1);
//        bdssm.setInputValue("birthRateChangeTimes", new RealParameter("0."));
//        bdssm.setInputValue("deathRateChangeTimes", new RealParameter("0."));
//        bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("0."));
//
//        //        // c3)
//        bdssm.setInputValue("intervalNumber", 1);
//        bdssm.setInputValue("rhoSamplingTimes", new RealParameter("2.7 3.7 6"));
//        bdssm.setInputValue("rho", new RealParameter("0 0 0.1"));
//        bdssm.initAndValidate();
////        System.out.println("\nc3) Likelihood: " + bdssm.calculateTreeLogLikelihood(tree));
//        assertEquals(-20.745947825183116, bdssm.calculateTreeLogLikelihood(tree), 1e-5);
//
//
//    }

    @Test
    public void testRhoD() throws Exception{

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", true);
        //        bdssm.setInputValue("birthRate", new RealParameter("2."));
        //        bdssm.setInputValue("deathRate", new RealParameter("1."));
        //        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));
        bdssm.setInputValue("R0", new RealParameter(new Double[]{4./3.}));
        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{1./3.}));
        //bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("birthRateChangeTimes", new RealParameter("0."));
        bdssm.setInputValue("deathRateChangeTimes", new RealParameter("0."));
        bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("0."));


        // d)
        bdssm.setInputValue("rho", null);
//        bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("intervalTimes", new RealParameter("0."));
        bdssm.initAndValidate();
//        System.out.println("\nd) Likelihood: " + bdssm.calculateTreeLogLikelihood(tree));
        assertEquals(-18.574104165089775, bdssm.calculateTreeLogLikelihood(tree), 1e-5);

    }
    @Test
    public void testRhoE() throws Exception{

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

//        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
//        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", true);
        //        bdssm.setInputValue("birthRate", new RealParameter("2."));
        //        bdssm.setInputValue("deathRate", new RealParameter("1."));
        //        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));
        bdssm.setInputValue("R0", new RealParameter(new Double[]{4./3.}));
        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
//        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{1./3.}));
        //bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("birthRateChangeTimes", new RealParameter("0."));
        bdssm.setInputValue("deathRateChangeTimes", new RealParameter("0."));
        bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("0."));


        //e) contemp tree:

        Tree tree = new TreeParser("((3:4,4:4):1,(1:2,2:2):3);", false);
        bdssm.setInputValue("tree", tree);

//        bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{0.}));
        //        bdssm.setInputValue("deathRate", new RealParameter("1.5"));
        //        bdssm.setInputValue("samplingRate", new RealParameter(new Double[]{0.}));
        bdssm.setInputValue("rho", new RealParameter("0.01"));
        bdssm.initAndValidate();

//        System.out.println("\ne) Contemp. TreeLikelihood: " + bdssm.calculateTreeLogLikelihood(tree));
        assertEquals(-8.130835517289412, bdssm.calculateTreeLogLikelihood(tree), 1e-5); //-8.130835517289412

    }

    @Test
    public void testLikelihoodCalculation1() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", false);

        // test without rate change
        bdssm.setInputValue("birthRate", new RealParameter("2."));
        bdssm.setInputValue("deathRate", new RealParameter("1."));
        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));

        bdssm.initAndValidate();
        bdssm.printTempResults = true;

        assertEquals(-19.0198, bdssm.calculateTreeLogLikelihood(tree), 1e-5);
    }

    @Test
    public void testLikelihoodCalculation2() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", false);

        // test with rate change outside tree range
        bdssm.setInputValue("birthRate", new RealParameter("2."));
        bdssm.setInputValue("deathRate", new RealParameter("1."));
        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));
        bdssm.setInputValue("forceRateChange", false);

        bdssm.initAndValidate();
        bdssm.printTempResults = true;

        assertEquals(-19.0198, bdssm.calculateTreeLogLikelihood(tree), 1e-5);
    }

    @Test
    public void testLikelihoodCalculation3() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
//        TreeIntervals intervals = new TreeIntervals();
//        intervals.init(tree);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", false);

        // test with rate change outside tree range
        //bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("birthRate", new RealParameter("2."));
        bdssm.setInputValue("deathRate", new RealParameter("1."));
        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));
        //bdssm.setInputValue("intervalTimes", new RealParameter("0."));
        bdssm.setInputValue("forceRateChange", false);

        bdssm.initAndValidate();
        bdssm.printTempResults = true;

        assertEquals(-19.0198, bdssm.calculateTreeLogLikelihood(tree), 1e-5);

    }


    @Test
    public void testLikelihoodCalculation4() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", false);


        // test with 1 rate change within interval
//        bdssm.setInputValue("intervalNumber", 2);
        bdssm.setInputValue("birthRate", new RealParameter("3. 2."));
        bdssm.setInputValue("deathRate", new RealParameter("2.5 1."));
        bdssm.setInputValue("samplingRate", new RealParameter("2. 0.5"));
        bdssm.setInputValue("intervalTimes", new RealParameter("0. 3."));

        bdssm.initAndValidate();
        bdssm.printTempResults = true;

        assertEquals(-33.7573, bdssm.calculateTreeLogLikelihood(tree), 1e-4);
    }

    @Test
    public void testLikelihoodCalculation5() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", false);


        PrintStream treeString = new PrintStream("out.tree");
        tree.log(1, treeString);

        // test with 2 rate changes
//        bdssm.setInputValue("intervalNumber", 3);
        bdssm.setInputValue("birthRate", new RealParameter("3. 2. 4."));
        bdssm.setInputValue("deathRate", new RealParameter("2.5 1. .5"));
        bdssm.setInputValue("samplingRate", new RealParameter("2. 0.5 1."));
        bdssm.setInputValue("intervalTimes", new RealParameter("0. 3. 4.5"));

        bdssm.initAndValidate();
        bdssm.printTempResults = true;

        assertEquals(-37.8056, bdssm.calculateTreeLogLikelihood(tree), 1e-4);

    }

    @Test
    public void testLikelihoodCalculation6() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));

        //same test with epi-parametrization

        bdssm.setInputValue("conditionOnSurvival", false);

//        bdssm.setInputValue("intervalNumber", 3);
        bdssm.setInputValue("intervalTimes", new RealParameter("0. 3. 4.5"));

        bdssm.setInputValue("R0", new RealParameter(new Double[]{2./3., 4./3., 8./3.}));
        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("4.5 1.5 1.5"));
        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{4./9., 1./3., 2./3.}));
        bdssm.initAndValidate();

        assertEquals(-37.8056, bdssm.calculateTreeLogLikelihood(tree), 1e-4);

    }

    @Test
    public void testMini() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("(1 : 1.5, 2 : 0.5);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("2.5"));
        bdssm.setInputValue("conditionOnSurvival", false);

//        bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("birthRate", new RealParameter("2."));
        bdssm.setInputValue("deathRate", new RealParameter("1."));
        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));

        bdssm.initAndValidate();

        assertEquals(-4.719294304452187, bdssm.calculateTreeLogLikelihood(tree), 1e-5);

    }

    @Test
    public void testLikelihoodCalculation4reverse() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", false);


        // test with 1 rate change within interval
//        bdssm.setInputValue("intervalNumber", 2);
        bdssm.setInputValue("birthRate", new RealParameter("3. 2."));
        bdssm.setInputValue("deathRate", new RealParameter("2.5 1."));
        bdssm.setInputValue("samplingRate", new RealParameter("2. 0.5"));
        bdssm.setInputValue("intervalTimes", new RealParameter("0. 3."));
        bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("3. 0."));

        bdssm.setInputValue("reverseTimeArrays", "false false true false");

        bdssm.initAndValidate();
        bdssm.printTempResults = true;

        assertEquals(-33.7573, bdssm.calculateTreeLogLikelihood(tree), 1e-4);
    }
    @Test

    public void testLikelihoodCalculation5reverse() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", false);


        PrintStream treeString = new PrintStream("out.tree");
        tree.log(1, treeString);

        // test with 2 rate changes
//        bdssm.setInputValue("intervalNumber", 3);
        bdssm.setInputValue("birthRate", new RealParameter("3. 2. 4."));
        bdssm.setInputValue("deathRate", new RealParameter("2.5 1. .5"));
        bdssm.setInputValue("samplingRate", new RealParameter("2. 0.5 1."));
        bdssm.setInputValue("intervalTimes", new RealParameter("0. 3. 4.5"));

        bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("1.5 3. 0."));
        bdssm.setInputValue("reverseTimeArrays", "false false true false");


        bdssm.initAndValidate();
        bdssm.printTempResults = true;

        assertEquals(-37.8056, bdssm.calculateTreeLogLikelihood(tree), 1e-4);

    }

    public void testLikelihoodCalculation5reverseAll() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", false);


        PrintStream treeString = new PrintStream("out.tree");
        tree.log(1, treeString);

        // test with 2 rate changes
//        bdssm.setInputValue("intervalNumber", 3);
        bdssm.setInputValue("birthRate", new RealParameter("3. 2. 4."));
        bdssm.setInputValue("deathRate", new RealParameter("2.5 1. .5"));
        bdssm.setInputValue("samplingRate", new RealParameter("2. 0.5 1."));
//         bdssm.setInputValue("intervalTimes", new RealParameter("0. 3. 4.5"));

        bdssm.setInputValue("birthRateChangeTimes", new RealParameter("1.5 3. 0."));
        bdssm.setInputValue("deathRateChangeTimes", new RealParameter("1.5 3. 0."));
        bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("1.5 3. 0."));
        bdssm.setInputValue("reverseTimeArrays", "true true true false");


        bdssm.initAndValidate();
        bdssm.printTempResults = true;

        assertEquals(-37.8056, bdssm.calculateTreeLogLikelihood(tree), 1e-4);

    }


    public void testTreeParser() throws Exception {

        TreeParser tree = new TreeParser();
        String newick = "(((1[&state='1']:1, (2[&state='0']:.5)[&state='1']:1.5)[&state='1']:2)[&state='0']:1, (3[&state='0']:1.5, (4[&state='1']:1.5)[&state='0']:1 )[&state='0']:2)[&state='0']:1;";
        tree.initByName("adjustTipHeights",false, "singlechild", true, "newick", newick);

        printNodeState(tree.getRoot());

    }

    void printNodeState(Node node){

        System.out.println("Node " + node.getNr() + " has colour " + node.metaDataString + "\t" + Integer.parseInt(node.metaDataString.split("=")[1].replaceAll("'","").replaceAll("\"","")));
        if (node.getLeft()!=null)
            printNodeState(node.getLeft());
        if (node.getRight()!=null)
            printNodeState(node.getRight());
    }

}


/* R code for test example

# calculate BDSSM treelikelihood for Tree('(((1:1,2:1):2,3:3):1,4:4);')


A <- function(b, g, psi){ return (sqrt((b - g - psi)^2 + 4*b*psi))}
B <- function(b, g, psi, p0){ return (-((1-p0)*b + g + psi)/A(b, g, psi))}
p_0 <- function(b, g, psi, A, B, t, ti){ return ((b + g + psi - A *((exp(A*(t - ti))*(1-B)-(1+B)))/(exp(A*(t-ti))*(1-B)+(1+B)))/(2*b))}
g <- function(A, B, t, ti){ return (4/(2*(1-B^2)+ exp(A*(t-ti))*(1-B)^2 + exp(A*(t-ti))*(1+B)^2  )) }

t0 = 0
t1 = 2
n1 = 3
x0=5
x1=4
x2=3
t1=2
x3=1
y=0

b0 = 1/200.
b1 = 1/300.
g0=1/4.
g1 = 1/7.
psi = 0.5

A0 = A(b0,g0,psi)
A1 = A(b1,g1,psi)

B0=B(b0, g0, psi, 1)

p0_0_t1=p_0(b0, g0, psi, A0, B0, t1, t1)
p0_0_x3=p_0(b0, g0, psi, A0, B0, 1, 0)
p0_0_y=p_0(b0, g0, psi, A0, B0, y, t0)

B1=B(b1, g1, psi, p0_0_t1)

g_x0 = g(A1,B1,x0,t1)
g_x1 = g(A1,B1,x1,t1)
g_y = g(A0,B0,y,t0)
g_x2 = g(A1,B1,x2,t1)
g_x3 = g(A0,B0,x3,t0)
g_0_t1 = g(A0,B0,t1,t0)

f_r1 = g_x0 * b1 * g_x1 * b1 * g_x2 * b0 * g_x3 * g_0_t1^3 * psi^4
f_r0 = f_r1 * p0_0_y^4

# test factors
T0 = g_x0 *psi^4
T1 = b1 * g_x1 * b1 * g_x2 * b0 * g_x3
T2 = p0_0_y^4
T3 = g_0_t1^3

*/