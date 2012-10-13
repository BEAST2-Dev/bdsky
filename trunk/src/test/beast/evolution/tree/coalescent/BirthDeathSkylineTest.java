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
    public void testRho() throws Exception{

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("orig_root", new RealParameter("1."));
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

        //
        //
//        // a)
        bdssm.setInputValue("rho", new RealParameter("0.01"));
        bdssm.setInputValue("rhoSamplingTimes", new RealParameter("6.0"));
        bdssm.initAndValidate();
//        System.out.println("\na) Likelihood: " + bdssm.calculateTreeLogLikelihood(tree));
        assertEquals(-22.54937791357737, bdssm.calculateTreeLogLikelihood(tree), 1e-5);
//-22.54937791357737

        // b)
        bdssm.setInputValue("intervalNumber", 3);

        bdssm.setInputValue("R0", new RealParameter(new Double[]{4./3., 4./3., 4./3.}));
        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter(new Double[]{1.5,1.5,1.5}));
        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{1./3., 1./3., 1./3.}));

        bdssm.setInputValue("birthRateChangeTimes", new RealParameter("0. 2.5 3.5"));
        bdssm.setInputValue("deathRateChangeTimes", new RealParameter("0. 2.5 3.5"));
        bdssm.setInputValue("samplingRateChangeTimes", new RealParameter("0. 2.5 3.5"));
        bdssm.setInputValue("rhoSamplingTimes", new RealParameter("2.5 3.5 6.0"));

//        bdssm.setInputValue("rho", new RealParameter("0.01 0.02 0.05"));
        bdssm.setInputValue("rho", new RealParameter("0.05 0.02 0.01"));
        bdssm.initAndValidate();
//        System.out.println("\nb) Likelihood: " + bdssm.calculateTreeLogLikelihood(tree));
        assertEquals(-28.19069830850966, bdssm.calculateTreeLogLikelihood(tree), 1e-5);
 //-28.19069830850966

        // c1)
        bdssm.setInputValue("R0", new RealParameter(new Double[]{4./3.}));
        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter(new Double[]{1.5}));
        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{1./3.}));
        bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("intervalTimes", new RealParameter("0."));
        bdssm.setInputValue("birthRateChangeTimes", null);
        bdssm.setInputValue("deathRateChangeTimes", null);
        bdssm.setInputValue("samplingRateChangeTimes", null);
        bdssm.setInputValue("rhoSamplingTimes", new RealParameter("6.0"));
        bdssm.setInputValue("rho", new RealParameter("0.1"));
        bdssm.initAndValidate();
//        System.out.println("\nc) Likelihood: " + bdssm.calculateTreeLogLikelihood(tree));
        assertEquals(-20.74594782518312, bdssm.calculateTreeLogLikelihood(tree), 1e-5);

        //        // c3)
        bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("rhoSamplingTimes", new RealParameter("2.7 3.7 6"));
        bdssm.setInputValue("rho", new RealParameter("0 0 0.1"));
        bdssm.initAndValidate();
//        System.out.println("\nc3) Likelihood: " + bdssm.calculateTreeLogLikelihood(tree));
        assertEquals(-20.745947825183116, bdssm.calculateTreeLogLikelihood(tree), 1e-5);
          
        // d)
        bdssm.setInputValue("rho", null);
        bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("intervalTimes", new RealParameter("0."));
        bdssm.initAndValidate();
//        System.out.println("\nd) Likelihood: " + bdssm.calculateTreeLogLikelihood(tree));
        assertEquals(-18.574104165089775, bdssm.calculateTreeLogLikelihood(tree), 1e-5);
        

        // --------------

        //e) contemp tree:

        tree = new Tree("((3:4,4:4):1,(1:2,2:2):3);");
        bdssm.setInputValue("tree", tree);

        bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{0.}));
        //        bdssm.setInputValue("deathRate", new RealParameter("1.5"));
        //        bdssm.setInputValue("samplingRate", new RealParameter(new Double[]{0.}));
        bdssm.setInputValue("intervalTimes", new RealParameter("0."));
        bdssm.setInputValue("rho", new RealParameter("0.01"));
        bdssm.initAndValidate();

//        System.out.println("\ne) Contemp. TreeLikelihood: " + bdssm.calculateTreeLogLikelihood(tree));
        assertEquals(-8.130835517289412, bdssm.calculateTreeLogLikelihood(tree), 1e-5); //-8.130835517289412
              

        //f) contemp tree:

        bdssm.setInputValue("intervalNumber", 3);
        bdssm.setInputValue("intervalTimes", new RealParameter("0. 2.5 3.5"));
        bdssm.setInputValue("rhoSamplingTimes", new RealParameter("2.5 3.5 6"));
        bdssm.setInputValue("R0", new RealParameter("3. .5 2."));
        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("1. 2. 1."));
        bdssm.setInputValue("samplingProportion", new RealParameter("0. 0. 0."));
        bdssm.setInputValue("rho", new RealParameter("0. 0. 0.01"));
        bdssm.initAndValidate();

//        System.out.println("\nf) Contemp. TreeLikelihood: " + bdssm.calculateTreeLogLikelihood(tree));
        assertEquals(-7.176847889608425, bdssm.calculateTreeLogLikelihood(tree), 1e-5);  //-7.176847889608425
               

    }

    @Test
    public void testLikelihoodCalculation() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        TreeIntervals intervals = new TreeIntervals();
        intervals.init(tree);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("orig_root", new RealParameter("1."));
        bdssm.setInputValue("conditionOnSurvival", false);


        PrintStream treeString = new PrintStream("out.tree");
        tree.log(1, treeString);

        // test without rate change
        //bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("birthRate", new RealParameter("2."));
        bdssm.setInputValue("deathRate", new RealParameter("1."));
        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));



//        bdssm.setInputValue("R0", new RealParameter(new Double[]{4./3.}));
//        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
//        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{1./3.}));

//        bdssm.setInputValue("intervalTimes", new RealParameter("0."));

        bdssm.initAndValidate();
        bdssm.printTempResults = true;

        assertEquals(-19.0198, bdssm.calculateTreeLogLikelihood(tree), 1e-5);


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


        // test with 1 rate change within interval
        bdssm.setInputValue("intervalNumber", 2);

        bdssm.setInputValue("birthRate", new RealParameter("3. 2."));
        //        bdssm.setInputValue("birthRate", new RealParameter("2.  3."));
        bdssm.setInputValue("deathRate", new RealParameter("2.5 1."));
        //        bdssm.setInputValue("deathRate", new RealParameter("1. 2.5"));
        bdssm.setInputValue("samplingRate", new RealParameter("2. 0.5"));
        //        bdssm.setInputValue("samplingRate", new RealParameter("0.5 2."));

        bdssm.setInputValue("intervalTimes", new RealParameter("0. 3."));

        bdssm.initAndValidate();
        bdssm.printTempResults = true;

        assertEquals(-33.7573, bdssm.calculateTreeLogLikelihood(tree), 1e-4);


        // test with 2 rate changes
        bdssm.setInputValue("intervalNumber", 3);
        bdssm.setInputValue("birthRate", new RealParameter("3. 2. 4."));
        //        bdssm.setInputValue("birthRate", new RealParameter("4. 2. 3."));
        bdssm.setInputValue("deathRate", new RealParameter("2.5 1. .5"));
        //        bdssm.setInputValue("deathRate", new RealParameter(".5 1. 2.5"));
        bdssm.setInputValue("samplingRate", new RealParameter("2. 0.5 1."));
        //        bdssm.setInputValue("samplingRate", new RealParameter("1. 0.5 2."));
        bdssm.setInputValue("intervalTimes", new RealParameter("0. 3. 4.5"));

        bdssm.initAndValidate();
        bdssm.printTempResults = true;

        assertEquals(-37.8056, bdssm.calculateTreeLogLikelihood(tree), 1e-4);


        //same test with epi-parametrization
        bdssm =  new BirthDeathSkylineModel();

        bdssm.setInputValue("conditionOnSurvival", false);

        intervals.init(tree);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("orig_root", new RealParameter("1."));
        bdssm.setInputValue("intervalNumber", 3);
        bdssm.setInputValue("intervalTimes", new RealParameter("0. 3. 4.5"));

        bdssm.setInputValue("R0", new RealParameter(new Double[]{2./3., 4./3., 8./3.}));
        bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("4.5 1.5 1.5"));
        bdssm.setInputValue("samplingProportion", new RealParameter(new Double[]{4./9., 1./3., 2./3.}));
        bdssm.initAndValidate();

        assertEquals(-37.8056, bdssm.calculateTreeLogLikelihood(tree), 1e-4);

        System.out.println("Test finished.");

    }

    @Test
    public void testMini() throws Exception {

        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        Tree tree = new TreeParser("(1 : 1.5, 2 : 0.5);",false);
        TreeIntervals intervals = new TreeIntervals();
        intervals.init(tree);
        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("orig_root", new RealParameter("1."));
        bdssm.setInputValue("conditionOnSurvival", false);

        bdssm.setInputValue("intervalNumber", 1);
        bdssm.setInputValue("birthRate", new RealParameter("2."));
        bdssm.setInputValue("deathRate", new RealParameter("1."));
        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));

        bdssm.initAndValidate();

        assertEquals(-4.719294304452187, bdssm.calculateTreeLogLikelihood(tree), 1e-5);

    }

    public void testTreeParser() throws Exception {

        TreeParser tree = new TreeParser();
        String newick = "(((1[&state='1']:1, (2[&state='0']:.5)[&state='1']:1.5)[&state='1']:2)[&state='0']:1, (3[&state='0']:1.5, (4[&state='1']:1.5)[&state='0']:1 )[&state='0']:2)[&state='0']:1;";
        tree.initByName("adjustTipHeights",false, "singlechild", true, "newick", newick);

        printNodeState(tree.getRoot());
        
    }

    void printNodeState(Node node){

        System.out.println("Node " + node.getNr() + " has colour " + node.m_sMetaData + "\t" + Integer.parseInt(node.m_sMetaData.split("=")[1].replaceAll("'","").replaceAll("\"","")));
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