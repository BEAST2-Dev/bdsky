package test.beast.evolution.tree.coalescent;

import beast.core.parameter.RealParameter;
import beast.evolution.speciation.BirthDeathSkylineModel;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import junit.framework.TestCase;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * @author Alexandra Gavryushkina
 */
public class SABirthDeathSkylineModel  extends TestCase {

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

}

