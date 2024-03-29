package beast.evolution.speciation;

import junit.framework.TestCase;
import beast.base.evolution.tree.TreeIntervals;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.core.Description;
import beast.base.evolution.tree.TreeParser;
import org.junit.Test;

import bdsky.evolution.speciation.BirthDeathSkylineModel;

import java.io.PrintStream;

/**
 * @author Denise Kuehnert
 * Date: Mar 29, 2011
 */

@Description("Test the BirthDeathSerialSkylineModel with a small tree example")

public class BirthDeathSkylineTest extends TestCase {

    @Test
    public void testmultiRho2tipsNoDeath() throws Exception {

        Tree tree = new TreeParser("(3[&type=0]: 4, 4[&type=0]: 4) ;", false);

        BirthDeathSkylineModel bdm = new BirthDeathSkylineModel();

        bdm.setInputValue("tree", tree);

        bdm.setInputValue("birthRate", new RealParameter("2.25"));
        bdm.setInputValue("deathRate", new RealParameter("0."));
        bdm.setInputValue("samplingRate", new RealParameter("0."));
//        bdm.setInputValue("rho", new RealParameter("0. 1."));
        bdm.setInputValue("rhoSamplingTimes", new RealParameter("0. 2.5"));
        bdm.setInputValue("reverseTimeArrays", "false false false true");
        bdm.setInputValue("origin", "5.");

        bdm.setInputValue("rho", new RealParameter("0.2 1."));

//        bdm.setInputValue("conditionOnSurvival", false);
//        bdm.initAndValidate();
//        assertEquals(-19.88536, bdm.calculateLogP(), 1e-4);   // this result is from R: LikBirthDeathMultiRho(2.25,0.,c(1,0.2),c(4.,5.),c(0.,0.),survival=0,root=0)

        bdm.setInputValue("rho", new RealParameter("0.2 0.6"));
        bdm.initAndValidate();
        assertEquals(-18.865767180278915, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R

        tree = new TreeParser("(3[&type=0]: 1.5, 4[&type=0]: 4) ;",false);
        bdm.setInputValue("tree", tree);

//        bdm.initAndValidate();
//        assertEquals(-15.646651247531981, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R

        bdm.setInputValue("conditionOnSurvival", false);
        bdm.setInputValue("rho", new RealParameter("0.2 0.6"));
        bdm.initAndValidate();
        assertEquals(-15.133091119955177, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R

        bdm.setInputValue("conditionOnSurvival", false);
        bdm.setInputValue("deathRate", "1.5");
        bdm.setInputValue("rho", new RealParameter("0.2 0.6"));
        bdm.initAndValidate();
        assertEquals(-8.637410990319223, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R

    }

    @Test
    public void testmultiRho2tips() throws Exception {

        Tree tree = new TreeParser("(3[&type=0]: 4, 4[&type=0]: 4) ;",false);

        BirthDeathSkylineModel bdm =  new BirthDeathSkylineModel();

        bdm.setInputValue("tree", tree);

        bdm.setInputValue("reproductiveNumber", new RealParameter("1.5"));
        bdm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
        bdm.setInputValue("samplingProportion", new RealParameter("0.") );
        bdm.setInputValue("rho", new RealParameter("0. 1.") );
        bdm.setInputValue("rhoSamplingTimes", new RealParameter("0. 2.5") );
        bdm.setInputValue("reverseTimeArrays", "false false false true");
        bdm.setInputValue("origin", "5.");

        bdm.setInputValue("conditionOnSurvival", true);
        bdm.initAndValidate();
        assertEquals(-9.151595, bdm.calculateLogP(), 1e-4);   // this result is from R: LikConstant(2.25,1.5,1,c(4.0,5.0),root=0,survival=1)

        bdm.setInputValue("conditionOnSurvival", false);
        bdm.initAndValidate();
        assertEquals(-10.2344, bdm.calculateLogP(), 1e-4);   // this result is from R: LikConstant(2.25,1.5,1,c(4.0,5.0),root=0,survival=0)

        bdm.setInputValue("rho", new RealParameter("0.2 1.") );

        bdm.setInputValue("conditionOnSurvival", true);
        bdm.initAndValidate();
        assertEquals(-10.569863754307026, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R

        bdm.setInputValue("conditionOnSurvival", false);
        bdm.initAndValidate();
        assertEquals(-11.616647517280333, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R


        tree = new TreeParser("(3[&type=0]: 1.5, 4[&type=0]: 4) ;",false);
        bdm.setInputValue("tree", tree);

        bdm.setInputValue("conditionOnSurvival", true);
        bdm.initAndValidate();
        assertEquals(-8.099631076932816, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R

        bdm.setInputValue("conditionOnSurvival", false);
        bdm.initAndValidate();
        assertEquals(-9.146414839906122, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R


    }

    @Test
    public void testmultiRho() throws Exception {

        Tree tree = new TreeParser("((3[&type=0]: 1.5, 4[&type=0]: 1.5)[&type=0]: 1 , (1[&type=0]: 2, 2[&type=0]: 2)[&type=0]: 3)[&type=0];",false);

        BirthDeathSkylineModel bdm =  new BirthDeathSkylineModel();

        bdm.setInputValue("tree", tree);

        bdm.setInputValue("reproductiveNumber", new RealParameter("1.5"));
        bdm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
        bdm.setInputValue("samplingProportion", new RealParameter("0.") );
        bdm.setInputValue("rho", new RealParameter("1. 1.") );
        bdm.setInputValue("rhoSamplingTimes", new RealParameter("0. 2.5") );
        bdm.setInputValue("reverseTimeArrays", "false false false true");
        bdm.setInputValue("conditionOnSurvival", false);
        bdm.setInputValue("origin", "10.");

        bdm.initAndValidate();

        assertEquals(Double.NEGATIVE_INFINITY, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R

        bdm.setInputValue("rho", new RealParameter("0.5 1.") );
        bdm.initAndValidate();
        assertEquals(-19.29773458054159, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R
//
        bdm.setInputValue("rho", new RealParameter("0.01 0.01") );
        bdm.initAndValidate();
        assertEquals(-14.676660714677519, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R

        bdm.setInputValue("conditionOnSurvival", true);

        //        bdm.setInputValue("rho", new RealParameter("0.5 1.") );
        //        bdm.initAndValidate();
        //        assertEquals(-14.37075396578577, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R
        //

//        bdm.setInputValue("rho", new RealParameter("0.") );
//        bdm.setInputValue("samplingProportion", new RealParameter("0.1") );
//        bdm.initAndValidate();
//        assertEquals(-13.738804989064041, bdm.calculateLogP(), 1e-4);   // this result is from BEAST, not double checked in R


    }
    @Test
    public void testSingleRho() throws Exception {
        BirthDeathSkylineModel bdm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((1: 4.5, 2: 4.5):1,3:5.5);",false);
        bdm.setInputValue("tree", tree);

        bdm.setInputValue("reproductiveNumber", new RealParameter("1.5"));
        bdm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
        bdm.setInputValue("samplingProportion", new RealParameter("0.") );
        bdm.setInputValue("conditionOnSurvival", false);
        bdm.setInputValue("origin", "10");
//        bdm.setInputValue("origin", "1e-100");
        // bdm.setInputValue("originIsRootEdge", true);
        bdm.setInputValue("rho", new RealParameter("0.01") );
        bdm.initAndValidate();
        assertEquals(-8.520565, bdm.calculateLogP(), 1e-4);   // this result is from R: LikConstant(2.25,1.5,0.01,c(4.5,5.5,10),root=0,survival=0)

        bdm.setInputValue("origin", "1e-100");
        bdm.setInputValue("originIsRootEdge", true);
        bdm.initAndValidate();
        assertEquals(-5.950979, bdm.calculateLogP(), 1e-4);   // this result is from R: LikConstant(2.25,1.5,0.01,c(4.5,5.5,5.5+1e-100),root=0,survival=0)

        bdm.setInputValue("conditionOnSurvival", true);
        bdm.initAndValidate();
        assertEquals(-4.431935, bdm.calculateLogP(), 1e-4);   // this result is from R: LikConstant(2.25,1.5,0.01,c(4.5,5.5,5.5+1e-100),root=0,survival=1)

        bdm.setInputValue("conditionOnSurvival", true);
        bdm.setInputValue("origin", "10");
        bdm.setInputValue("originIsRootEdge", false);
        bdm.initAndValidate();
        assertEquals(-7.404227, bdm.calculateLogP(), 1e-4);   // this result is from R: LikConstant(2.25,1.5,0.01,c(4.5,5.5,5.5+1e-100),root=0,survival=1)

    }

    @Test
    public void testYule() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree1 = new TreeParser("((A:1.0,B:1.0):1.0,C:2.0);",false);
        bdssm.setInputValue("tree", tree1);

        bdssm.setInputValue("birthRate", new RealParameter("10."));
        bdssm.setInputValue("deathRate", new RealParameter("0."));
        bdssm.setInputValue("samplingRate", new RealParameter("0."));
        bdssm.setInputValue("rho", new RealParameter("1."));
        bdssm.setInputValue("origin", new RealParameter("10."));
        bdssm.setInputValue("origin", new RealParameter("10."));

        bdssm.setInputValue("conditionOnSurvival", false);

        bdssm.initAndValidate();

        double logP1 = bdssm.calculateTreeLogLikelihood(tree1);



        Tree tree = new TreeParser("((A:1.0,B:1.0):2.0,C:3.0);",false);
        bdssm.setInputValue("tree", tree);

        bdssm.setInputValue("birthRate", new RealParameter("10."));
        bdssm.setInputValue("deathRate", new RealParameter("0."));
        bdssm.setInputValue("samplingRate", new RealParameter("0."));
        bdssm.setInputValue("rho", new RealParameter("1."));
        bdssm.setInputValue("origin", new RealParameter("10."));

        bdssm.setInputValue("conditionOnSurvival", false);

        bdssm.initAndValidate();

        double logP2 = bdssm.calculateTreeLogLikelihood(tree);

        assertEquals(logP1-logP2,10.0);

    }

    @Test
    public void testLikelihoodCalculationSimple() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("10."));
        bdssm.setInputValue("conditionOnSurvival", false);
        bdssm.setInputValue("removalProbability", "1");

        bdssm.setInputValue("reproductiveNumber", new RealParameter("1.5"));
        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
        bdssm.setInputValue("samplingProportion", new RealParameter("0.3") );

        bdssm.initAndValidate();
        bdssm.printTempResults = false;

        assertEquals(-26.105360134266082, bdssm.calculateTreeLogLikelihood(tree), 1e-4);
    }


    @Test
    public void testLikelihoodCalculationSimple2() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("10."));
        bdssm.setInputValue("conditionOnSurvival", false);
        bdssm.setInputValue("removalProbability", "1");


        bdssm.setInputValue("netDiversification", new RealParameter("0.75"));
        bdssm.setInputValue("turnOver", new RealParameter("0.6666667"));
        bdssm.setInputValue("samplingProportion", new RealParameter("0.3") );

        bdssm.initAndValidate();
        bdssm.printTempResults = false;

        assertEquals(-26.105360134266082, bdssm.calculateTreeLogLikelihood(tree), 1e-4);
    }



    @Test
    public void testLikelihoodCalculationSimpleForBDMM() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("10."));
        bdssm.setInputValue("conditionOnSurvival", true);
        bdssm.setInputValue("removalProbability", "0.9");


        bdssm.setInputValue("reproductiveNumber", new RealParameter("1.5"));
        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
        bdssm.setInputValue("samplingProportion", new RealParameter("0.3") );

        bdssm.initAndValidate();
        bdssm.printTempResults = false;

        assertEquals(-25.991511346557598, bdssm.calculateTreeLogLikelihood(tree), 1e-4);
    }

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
                    bdssm.setInputValue("rho", new RealParameter("0. 0. 0.01"));
                    bdssm.setInputValue("birthRateChangeTimes", new RealParameter("0. 1. 1.5"));
                    bdssm.setInputValue("deathRateChangeTimes", new RealParameter("0. 1. 1.5"));
                    bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("0. 1. 1.5"));
                    bdssm.setInputValue("rhoSamplingTimes", new RealParameter("0. 1. 1.5"));

                    bdssm.initAndValidate();
                    bdssm.printTempResults = false;

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
                case 3:{
                    BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();
                    bdssm.setInputValue("tree", tree);
                    bdssm.setInputValue("origin", new RealParameter("2"));
                    bdssm.setInputValue("conditionOnSurvival", false);

                    bdssm.setInputValue("birthRate", new RealParameter("3. 2. 4. 4."));
                    bdssm.setInputValue("deathRate", new RealParameter("2.5 1. .5 0.5"));
                    bdssm.setInputValue("samplingRate", new RealParameter("2. 0.5 1. 2.0"));
                    bdssm.setInputValue("rho", new RealParameter("0.05 0.01"));
                    bdssm.setInputValue("rhoSamplingTimes","0. 1.");
                    bdssm.setInputValue("intervalTimes", new RealParameter("0. 0.5 1. 1.1"));
                    bdssm.initAndValidate();
                    assertEquals(-109.6105592003162, bdssm.calculateTreeLogLikelihood(tree), 1e-4);
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
    public void test1intOrig() throws Exception{
        Tree tree = new TreeParser("(((((t1:0.4595008531,t25:0.4595008531):0.3373053072,t23:0.3567584538):0.007310819036,t16:0.3489190732):0.331009529,((t18:0.03315384045,t14:0.03315384045):0.5063451374,(t10:0.4211543131,t15:0.4211543131):0.1183446648):0.5956275305):0.1158090878,((t19:0.9429393194,((t6:0.363527235,t11:0.4417423167):0.01881829549,((((t3:0.3071904376,(((t24:0.01065209364,t13:0.01065209364):0.06076485145,t8:0.07141694509):0.123620245,(t22:0.1616119808,t2:0.1616119808):0.03342520927):0.1121532475):0.24520579,t9:0.5523962276):0.3852615426,(((t20:0.2935970782,(t17:0.06569090089,t4:0.06569090089):0.2279061773):0.08350780408,(t21:0.05109047139,t5:0.05109047139):0.3260144109):0.2298344132,t7:0.6069392955):0.3307184747):0.01206284377,t26:0.9497206139):0.05755333197):0.03290891884):0.07263755325,t12:1.112820418):0.1381151782);", false);


        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter(new Double[]{1e-15}));
        bdssm.setInputValue("originIsRootEdge", true);
//        bdssm.setInputValue("conditionOnSurvival", false);
        bdssm.setInputValue("birthRate", new RealParameter("2."));
        bdssm.setInputValue("deathRate", new RealParameter("0.5"));
        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));

        bdssm.setInputValue("rho", new RealParameter("1."));
        bdssm.initAndValidate();
//        System.out.println("\na) Likelihood: " + bdssm.calculateTreeLogLikelihood(tree));
        assertEquals(-19.74719012042364, bdssm.calculateTreeLogLikelihood(tree), 1e-5);   // conditionedOnSurvival = true
//        assertEquals(-19.946795863214046, bdssm.calculateTreeLogLikelihood(tree), 1e-5);   // conditionedOnSurvival = false

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
        bdssm.setInputValue("reproductiveNumber", new RealParameter(new Double[]{4./3.}));
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


    @Test
    public void testRhoA2() throws Exception{

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", true);
        //        bdssm.setInputValue("birthRate", new RealParameter("2."));
        //        bdssm.setInputValue("deathRate", new RealParameter("1."));
        //        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));
        bdssm.setInputValue("netDiversification", new RealParameter("0.5"));
        bdssm.setInputValue("turnOver", new RealParameter("0.75"));
        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{1./3.}));
        //bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("birthRateChangeTimes", new RealParameter("0."));
        bdssm.setInputValue("deathRateChangeTimes", new RealParameter("0."));
        bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("0."));

        bdssm.setInputValue("rho", new RealParameter("0.01"));
        bdssm.initAndValidate();
//        System.out.println("\na) Likelihood: " + bdssm.calculateTreeLogLikelihood(tree));
        assertEquals(-22.54937791357737, bdssm.calculateTreeLogLikelihood(tree), 1e-5);
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
//        bdssm.setInputValue("reproductiveNumber", new RealParameter(new Double[]{4./3.}));
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
//        bdssm.setInputValue("reproductiveNumber", new RealParameter(new Double[]{4./3., 4./3., 4./3.}));
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
        bdssm.setInputValue("reproductiveNumber", new RealParameter(new Double[]{4./3.}));
        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{1./3.}));
        //bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("birthRateChangeTimes", new RealParameter("0."));
        bdssm.setInputValue("deathRateChangeTimes", new RealParameter("0."));
        bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("0."));

        // c1)
        bdssm.setInputValue("reproductiveNumber", new RealParameter(new Double[]{4./3.}));
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


    @Test
    public void testRhoC2() throws Exception{

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", true);
        //        bdssm.setInputValue("birthRate", new RealParameter("2."));
        //        bdssm.setInputValue("deathRate", new RealParameter("1."));
        //        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));
        bdssm.setInputValue("netDiversification", new RealParameter("0.5"));
        bdssm.setInputValue("turnOver", new RealParameter("0.75"));
        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{1./3.}));
        //bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("birthRateChangeTimes", new RealParameter("0."));
        bdssm.setInputValue("deathRateChangeTimes", new RealParameter("0."));
        bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("0."));

        // c1)
        bdssm.setInputValue("netDiversification", new RealParameter("0.5"));
        bdssm.setInputValue("turnOver", new RealParameter("0.75"));
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
//        bdssm.setInputValue("reproductiveNumber", new RealParameter(new Double[]{4./3.}));
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
        bdssm.setInputValue("reproductiveNumber", new RealParameter(new Double[]{4./3.}));
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
    public void testRhoD2() throws Exception{

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", true);
        //        bdssm.setInputValue("birthRate", new RealParameter("2."));
        //        bdssm.setInputValue("deathRate", new RealParameter("1."));
        //        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));
        bdssm.setInputValue("netDiversification", new RealParameter("0.5"));
        bdssm.setInputValue("turnOver", new RealParameter("0.75"));
        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{1./3.}));
        //bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("birthRateChangeTimes", new RealParameter("0."));
        bdssm.setInputValue("deathRateChangeTimes", new RealParameter("0."));
        bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("0."));

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
        bdssm.setInputValue("reproductiveNumber", new RealParameter(new Double[]{4./3.}));
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
    public void testRhoE2() throws Exception{

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

//        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
//        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", true);
        //        bdssm.setInputValue("birthRate", new RealParameter("2."));
        //        bdssm.setInputValue("deathRate", new RealParameter("1."));
        //        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));
        bdssm.setInputValue("netDiversification", new RealParameter("0.5"));
        bdssm.setInputValue("turnOver", new RealParameter("0.75"));
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

        bdssm.setInputValue("reproductiveNumber", new RealParameter(new Double[]{2./3., 4./3., 8./3.}));
        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("4.5 1.5 1.5"));
        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{4./9., 1./3., 2./3.}));
        bdssm.initAndValidate();

        assertEquals(-37.8056, bdssm.calculateTreeLogLikelihood(tree), 1e-4);

    }


    @Test
    public void testLikelihoodCalculation7() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));

        //same test with epi-parametrization

        bdssm.setInputValue("conditionOnSurvival", false);

//        bdssm.setInputValue("intervalNumber", 3);
        bdssm.setInputValue("intervalTimes", new RealParameter("0. 3. 4.5"));

//        bdssm.setInputValue("reproductiveNumber", new RealParameter(new Double[]{2./3., 4./3., 8./3.}));
//        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("4.5 1.5 1.5"));
        bdssm.setInputValue("netDiversification", new RealParameter("-1.5 0.5 2.5"));
        bdssm.setInputValue("turnOver", new RealParameter("1.5 0.75 0.375"));
        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{4./9., 1./3., 2./3.}));
        bdssm.initAndValidate();

        assertEquals(-37.8056, bdssm.calculateTreeLogLikelihood(tree), 1e-4);
    }


    @Test
    public void testVeronika() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((((t73:0.01530208968,t75:1.156805242):0.7271576473,((t89:0.4480175162,t7:1.074752261):0.9766680278,((((t27:0.9767821943,t44:2.499842623):0.06922451595,t88:0.4686480382):3.874499393,((((t95:0.05986914016,t8:1.179246844):1.605580582,(((t23:0.5241800061,t24:1.907396834):0.7562073329,t72:1.108092025):0.1765734653,t94:0.5636128448):4.094115933):2.101413042,((((t4:0.4136749452,t68:2.118798524):1.546886938,t33:3.106435992):1.389941118,t39:0.9960299057):3.843767009,t41:2.556867517):0.7400688616):1.817360336,t10:0.256382727):3.241502636):0.1069251666,t61:0.9035797139):1.226269224):2.39211925):0.5107025246,((((((((t67:6.628374555,(t62:0.16943938,(t87:0.09916894415,(t57:2.422526259,t96:1.721313766):2.127590837):0.04973452396):0.04740879255):0.05006949562,t97:0.3717134135):1.767710383,t66:1.600025384):0.252854429,(t22:1.138167137,((t31:0.3295741843,(t99:2.877648776,(t42:0.7597191713,t28:2.726267134):0.389066579):3.382678787):3.40448748,((((t1:3.095420121,((t69:1.346146083,((t84:5.053722198,(t36:1.436074204,t16:0.04158946158):1.8259038):0.4341720219,(t79:0.9463279059,t6:5.062900921):1.675701577):1.373223424):0.2797454029,(((t15:0.9615474372,t80:0.9107358375):0.2921729326,t46:3.260674574):3.204694482,t100:0.3433167979):2.185334754):0.1344187902):0.966406539,t18:0.8853363796):1.26869644,(t17:4.540458062,t92:1.784855062):0.07509309778):0.04720669724,t81:1.494819847):0.5850014075):2.704587599):0.8940644426):0.8320247232,t60:2.83755279):0.7854434687,(t82:0.4053230926,(t51:0.3298092201,t30:2.098196195):0.04537985467):0.0499363901):1.085620168,(((t52:0.3203829029,(t71:0.04132246426,(((((t25:1.599591892,(t38:0.3831072194,t32:0.7282078865):1.182042631):3.050266335,t48:0.2952627683):2.046776454,(t14:0.0267379961,((t58:0.6469808535,(t12:1.795493868,t90:0.008599763368):0.5775253675):1.697312748,t35:3.275512845):0.4792751271):2.005515846):1.374856184,t98:1.224523058):0.2387958775,(t26:0.2286381189,((((t63:0.7896555658,(t49:0.03358752082,t21:5.113921775):0.03643435543):1.535988246,(t43:0.3065202928,((t64:4.921085652,t54:2.15858213):1.262544996,(t86:0.3518927006,t93:4.653646325):0.5668623644):0.9812764567):0.01211019273):0.6086295407,(((((t47:0.475461723,t50:0.123102006):0.008164226095,t83:0.2253500598):0.7998366182,t74:0.4107139442):0.8356302041,(t77:5.916373763,((t5:1.258774918,(t2:1.377399228,t9:1.589954493):0.8200025224):0.775524383,t37:0.7927213631):2.083755218):0.2804356569):0.7645776215,t13:1.572106645):0.6566036856):0.05619800187,((t40:0.4101428013,t59:0.009648414467):0.3950404599,t3:0.90401304):0.5449025378):0.4961783597):0.3148637741):1.07255348):0.2142166838):0.3037582735,t85:1.558286393):7.531358997,t70:0.9296960371):0.2269842114):1.532885028,((t20:0.3673285773,(t34:0.226258026,t29:2.869768924):1.893042559):1.166459247,((((t91:2.079325562,(t76:3.475974831,((t45:1.142879794,(t19:0.3007867251,t65:2.76581448):4.091195997):0.1588251617,t55:4.146246866):0.5259188711):0.9303453928):1.031361348,t56:0.7870142322):0.3289922425,t53:0.5056011922):1.883232687,t78:0.8842096268):0.5595520835):5.06033768):1.110615062):0.06633297341,t11:1.111858009):0.7862074588;",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter(new Double[]{0.7862+tree.getRoot().getHeight()}));
        bdssm.setInputValue("conditionOnSurvival", false);

        bdssm.setInputValue("reproductiveNumber", new RealParameter("7.767395"));
        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("0.1057396"));
//        bdssm.setInputValue("netDiversification", new RealParameter("0.7155816"));
//        bdssm.setInputValue("turnOver", new RealParameter("0.1287433"));
        bdssm.setInputValue("samplingProportion", new RealParameter("0.5"));

        bdssm.initAndValidate();

        assertEquals(-455.05623026907483, bdssm.calculateTreeLogLikelihood(tree), 1e-5);

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

    @Test
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


    @Test
    public void testTreeParser() throws Exception {

        TreeParser tree = new TreeParser();
        String newick = "(((1[&state='1']:1, (2[&state='0']:0.5)[&state='1']:1.5)[&state='1']:2)[&state='0']:1, (3[&state='0']:1.5, (4[&state='1']:1.5)[&state='0']:1 )[&state='0']:2)[&state='0']:1;";
        tree.initByName("adjustTipHeights",false, "singlechild", true, "newick", newick);

        printNodeState(tree.getRoot());

    }


    @Test
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