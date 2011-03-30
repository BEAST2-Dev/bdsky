package test.beast.evolution.tree.coalescent;

import junit.framework.TestCase;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.evolution.tree.Tree;
import beast.evolution.speciation.BirthDeathSerialSkylineModel;
import beast.core.parameter.RealParameter;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: dkuh004
 * Date: Mar 29, 2011
 * Time: 2:19:46 PM
 */
public class BirthDeathSerialSkylineTest extends TestCase {

    @Test
    public void testSkyline() throws Exception {

        RealParameter times = new RealParameter(new double[]{2.0});
        RealParameter birthRate = new RealParameter(new double[]{1/200., 1/300.});
        RealParameter deathRate = new RealParameter(new double[]{1/4., 1/7.});


        BirthDeathSerialSkylineModel bdssm = new BirthDeathSerialSkylineModel();

        bdssm.setInputValue("times", times);
        bdssm.setInputValue("birthRate", birthRate);
        bdssm.setInputValue("deathRate", deathRate);

        bdssm.setInputValue("serialSamplingRate", 0.5);
        bdssm.setInputValue("extantSamplingRate", 0.5);
        bdssm.setInputValue("relativeDeath", false);
        bdssm.setInputValue("sampledIndividualsRemainInfectious", true);
        bdssm.setInputValue("finalTimeInterval", 0.);

        bdssm.initAndValidate();
        
        Tree tree = new Tree("(((1:1,2:1):2,3:3):1,4:4);");
        TreeIntervals intervals = new TreeIntervals();
        intervals.init(tree);

//        for (int i = 0; i < size; i += 1) {
//            System.out.println("deathRate at time " + i + " is " + bdssm.deathRate(i));
//            System.out.println("p0 at time " + i + " is " + bdssm.p0(i, i));
//        }

        double logL = bdssm.calculateTreeLogLikelihood(tree);

        
        System.out.println("loglikelihood: " + logL + " " + Math.exp(logL));
        


    }


}
