package beast.evolution.speciation;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import junit.framework.TestCase;
import org.junit.Test;

import bdsky.evolution.speciation.BirthDeathSkylineDiversifiedSampling;
import bdsky.evolution.speciation.BirthDeathSkylineModel;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * @author Alexandra Gavryushkina
 */
public class SABirthDeathSkylineModelTest extends TestCase {

    @Test
    public void testLikelihoodCalculationNoRho() throws Exception {

        ArrayList<String> taxa = new ArrayList<String>(Arrays.asList("1", "2", "3"));
        Tree tree = new TreeParser(taxa, "((1:2.0,2:0.5):0.5,3:0.0):0.0", 1, false);

        BirthDeathSkylineModel model = new BirthDeathSkylineModel();
        model.setInputValue("tree", tree);
        model.setInputValue("origin", new RealParameter("3."));
        model.setInputValue("birthRate", new RealParameter("1.5 0.8"));
        model.setInputValue("deathRate", new RealParameter("0.4 0.3"));
        model.setInputValue("samplingRate", new RealParameter("1.2 0.4"));
        model.setInputValue("removalProbability", new RealParameter("0.8 0.7"));
        model.setInputValue("reverseTimeArrays", "true true true false true");
        model.setInputValue("birthRateChangeTimes", new RealParameter("1. 0."));
        model.setInputValue("deathRateChangeTimes", new RealParameter("1. 0."));
        model.setInputValue("samplingRateChangeTimes", new RealParameter("1. 0."));
        model.setInputValue("removalProbabilityChangeTimes", new RealParameter("1. 0."));
        model.setInputValue("conditionOnSurvival", true);
        model.initAndValidate();

        assertEquals(-7.22936943044961, model.calculateTreeLogLikelihood(tree), 1e-14);

    }

    @Test
    public void testLikelihoodCalculationWithRho() throws Exception {
        ArrayList<String> taxa1 = new ArrayList<String>(Arrays.asList("1", "2", "3"));
        Tree tree1 = new TreeParser(taxa1, "((1:1.0,2:0.0):1.0,3:0.0):0.0", 1, false);

        BirthDeathSkylineModel model = new BirthDeathSkylineModel();
        model.setInputValue("tree", tree1);
        model.setInputValue("origin", new RealParameter("3."));
        model.setInputValue("birthRate", new RealParameter("1.2 1."));
        model.setInputValue("deathRate", new RealParameter("0.4 0.3"));
        model.setInputValue("samplingRate", new RealParameter("0.1 0.4"));
        model.setInputValue("removalProbability", new RealParameter("0.8 0.7"));
        model.setInputValue("rho", new RealParameter("0.1 0.4"));
        model.setInputValue("reverseTimeArrays", "true true true true true");
        model.setInputValue("birthRateChangeTimes", new RealParameter("1. 0."));
        model.setInputValue("deathRateChangeTimes", new RealParameter("1. 0."));
        model.setInputValue("samplingRateChangeTimes", new RealParameter("1. 0."));
        model.setInputValue("rhoSamplingTimes", new RealParameter("1. 0."));
        model.setInputValue("removalProbabilityChangeTimes", new RealParameter("1. 0."));
        model.setInputValue("conditionOnSurvival", true);
        model.initAndValidate();

        assertEquals(-10.7246150935598, model.calculateTreeLogLikelihood(tree1), 1e-14);

        ArrayList<String> taxa2 = new ArrayList<String>(Arrays.asList("1", "2"));
        Tree tree2 = new TreeParser(taxa2, "(1:2.0,2:1.0):0.0", 1, false);
        model.setInputValue("tree", tree2);

        assertEquals(-6.03094565272802, model.calculateTreeLogLikelihood(tree2), 1e-14);

    }


    @Test
    public void testLikelihoodCalculationThreeIntervalsWithRho() throws Exception {

        ArrayList<String> taxa = new ArrayList<String>(Arrays.asList("1", "2", "3"));
        Tree tree = new TreeParser(taxa, "((1:1.5,3:0.0):1.0,2:1.5):0.0", 1, false);

        BirthDeathSkylineModel model = new BirthDeathSkylineModel();
        model.setInputValue("tree", tree);
        model.setInputValue("origin", new RealParameter("3."));
        model.setInputValue("birthRate", new RealParameter("1.5 0.5 0.6"));
        model.setInputValue("deathRate", new RealParameter("0.4 0.3 0.5"));
        model.setInputValue("samplingRate", new RealParameter("0.3 0.8 0.4"));
        model.setInputValue("rho", new RealParameter("0.0 0.5 0.9"));
        model.setInputValue("removalProbability", new RealParameter("0.8 0.5 0.9"));
        model.setInputValue("reverseTimeArrays", "true true true true true");
        model.setInputValue("rhoSamplingTimes", new RealParameter("1. 2. 0."));
        model.setInputValue("intervalTimes", new RealParameter("1. 2. 0."));
        model.setInputValue("conditionOnSurvival", true);
        model.initAndValidate();

        assertEquals(-7.64384832894999, model.calculateTreeLogLikelihood(tree), 1e-14);

        model.setInputValue("rho", new RealParameter("0.3 0.5 0.9"));

        assertEquals(-8.51026663900650, model.calculateTreeLogLikelihood(tree), 1e-14);

    }

    @Test
    public void testLikelihoodCalculationThreeIntervalsWithRhoConditionOnRhoSampling() throws Exception {

        ArrayList<String> taxa1 = new ArrayList<String>(Arrays.asList("1", "2", "3"));
        Tree tree1 = new TreeParser(taxa1, "((1:1.5,3:0.0):1.0,2:1.7):0.0", 1, false);

        BirthDeathSkylineModel model = new BirthDeathSkylineModel();
        model.setInputValue("tree", tree1);
        model.setInputValue("origin", new RealParameter("3."));
        model.setInputValue("birthRate", new RealParameter("1.5 0.5 0.6"));
        model.setInputValue("deathRate", new RealParameter("0.4 0.3 0.5"));
        model.setInputValue("samplingRate", new RealParameter("0.3 0.8 0.4"));
        model.setInputValue("rho", new RealParameter("0.0 0.9"));
        model.setInputValue("removalProbability", new RealParameter("0.0"));
        model.setInputValue("reverseTimeArrays", "true true true true true");
        model.setInputValue("birthRateChangeTimes", new RealParameter("1. 2. 0."));
        model.setInputValue("deathRateChangeTimes", new RealParameter("1. 2. 0."));
        model.setInputValue("samplingRateChangeTimes", new RealParameter("1. 2. 0."));
        model.setInputValue("rhoSamplingTimes", new RealParameter("1. 0."));
        model.setInputValue("conditionOnSurvival", false);
        model.setInputValue("conditionOnRhoSampling", true);
        model.initAndValidate();

        assertEquals(-7.49698456787236, model.calculateTreeLogLikelihood(tree1), 1e-14);

        ArrayList<String> taxa2 = new ArrayList<String>(Arrays.asList("1", "2", "3", "4"));
        Tree tree2 = new TreeParser(taxa2, "((1:0.8,3:0.0):1.5,(2:1.5,4:0.4):0.8):0.0", 1, false);
        model.setInputValue("tree", tree2);
        model.setInputValue("rho", new RealParameter("0.0 0.3"));
        model.initAndValidate();

        assertEquals(-8.80702177958906, model.calculateTreeLogLikelihood(tree2), 1e-14);
    }

    @Test
    public void testLikelihoodCalculationTwoIntervalsWithRhoConditionOnRhoSamplingAndRoot() throws Exception {

        ArrayList<String> taxa1 = new ArrayList<String>(Arrays.asList("1", "2", "3"));
        Tree tree1 = new TreeParser(taxa1, "((1:1.5,3:0.0):1.0,2:1.7):0.0", 1, false);

        BirthDeathSkylineModel model = new BirthDeathSkylineModel();
        model.setInputValue("tree", tree1);
        model.setInputValue("birthRate", new RealParameter("1.5 0.5"));
        model.setInputValue("deathRate", new RealParameter("0.4 0.3"));
        model.setInputValue("samplingRate", new RealParameter("0.3 0.8"));
        model.setInputValue("rho", new RealParameter("0.9"));
        model.setInputValue("removalProbability", new RealParameter("0.0"));
        model.setInputValue("reverseTimeArrays", "true true true true true");
        model.setInputValue("birthRateChangeTimes", new RealParameter("0. 1."));
        model.setInputValue("deathRateChangeTimes", new RealParameter("0. 1."));
        model.setInputValue("samplingRateChangeTimes", new RealParameter("0. 1."));
        model.setInputValue("rhoSamplingTimes", new RealParameter("0.0"));
        model.setInputValue("conditionOnRoot", true);
        model.setInputValue("conditionOnSurvival", false);
        model.setInputValue("conditionOnRhoSampling", true);
        model.initAndValidate();

        // the true value is calculated in R 'scripts/LikelihoodForConditionOnRootAndRhoTest.R'
        assertEquals(-8.55232499940742, model.calculateTreeLogLikelihood(tree1), 1e-14);
    }

    @Test
    public void testLikelihoodCalculationDiversifiedSampling() throws Exception {

        Tree tree = new TreeParser("(f4:12.78,((((f2:5.06,(t5:11.81,t1:11.81):1.25):3.89,t3:16.95):4.88,((f3:1.92,t2:10.92):9.08,f5:0.00):1.83):1.15,(t4:5.00,f1:0.00):17.98):4.80):2.22");

        BirthDeathSkylineDiversifiedSampling model = new BirthDeathSkylineDiversifiedSampling();
        model.setInputValue("tree", tree);
        model.setInputValue("netDiversification", new RealParameter("0.3 0.2"));
        model.setInputValue("turnOver", new RealParameter("0.5"));
        model.setInputValue("samplingProportion", new RealParameter("0.2 0.3 0.4"));
        model.setInputValue("rho", new RealParameter("0.0 0.5"));
        model.setInputValue("removalProbability", new RealParameter("0.0"));
        model.setInputValue("reverseTimeArrays", "true true true true true");
        model.setInputValue("birthRateChangeTimes", new RealParameter("0.0 14."));
        model.setInputValue("deathRateChangeTimes", new RealParameter("0.0"));
        model.setInputValue("samplingRateChangeTimes", new RealParameter("0.0 10. 18."));
        model.setInputValue("rhoSamplingTimes", new RealParameter("0.0 10."));
        model.setInputValue("conditionOnSurvival", true);  // conditionOnRhoSampling is false
        model.setInputValue("conditionOnRoot", true);
        model.initAndValidate();

        assertEquals(-69.62623542386, model.calculateTreeLogLikelihood(tree), 1e-8);

        model.setInputValue("origin", new RealParameter("30."));
        model.setInputValue("conditionOnRoot", false);  // condition on origin

        assertEquals(-71.73431504021, model.calculateTreeLogLikelihood(tree), 1e-8);
    }

    @Test
    public void testLikelihoodCalculationRandomSampling() throws Exception {

        Tree tree = new TreeParser("(f4:12.78,((((f2:5.06,(t5:11.81,t1:11.81):1.25):3.89,t3:16.95):4.88,((f3:1.92,t2:10.92):9.08,f5:0.00):1.83):1.15,(t4:5.00,f1:0.00):17.98):4.80):2.22");

        BirthDeathSkylineModel model = new BirthDeathSkylineModel();
        model.setInputValue("tree", tree);
        // model.setInputValue("origin", new RealParameter("30."));
        model.setInputValue("netDiversification", new RealParameter("0.3 0.2"));
        model.setInputValue("turnOver", new RealParameter("0.5"));
        model.setInputValue("samplingProportion", new RealParameter("0.2 0.3 0.4"));
        model.setInputValue("rho", new RealParameter("0.0 0.5"));
        model.setInputValue("removalProbability", new RealParameter("0.0"));
        model.setInputValue("reverseTimeArrays", "true true true true true");
        model.setInputValue("birthRateChangeTimes", new RealParameter("0.0 14."));
        model.setInputValue("deathRateChangeTimes", new RealParameter("0.0"));
        model.setInputValue("samplingRateChangeTimes", new RealParameter("0.0 10. 18."));
        model.setInputValue("rhoSamplingTimes", new RealParameter("0.0 10."));
        model.setInputValue("conditionOnSurvival", true);  // conditionOnRhoSampling is false
        model.setInputValue("conditionOnRoot", true);
        model.initAndValidate();

        assertEquals(-72.18435989275, model.calculateTreeLogLikelihood(tree), 1e-8);

        // mimic diversified sampling
        model.setInputValue("samplingProportion", new RealParameter("0.2 0.3 0.4 0.0"));
        model.setInputValue("rho", new RealParameter("0.0 1.0"));
        model.setInputValue("samplingRateChangeTimes", new RealParameter("0.0 4.75 10. 18."));
        model.initAndValidate();

        // by adding the correction term, should equal to
        assertEquals(-69.62623542386, model.calculateTreeLogLikelihood(tree)-1.37038559905, 1e-8);
    }
}

