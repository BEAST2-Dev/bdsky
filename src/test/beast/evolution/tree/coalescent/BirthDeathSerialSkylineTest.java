package test.beast.evolution.tree.coalescent;

import junit.framework.TestCase;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.evolution.tree.Tree;
import beast.evolution.speciation.BirthDeathSerialSkylineModel;
import beast.core.parameter.RealParameter;
import beast.core.Description;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Denise Kuehnert
 * Date: Mar 29, 2011
 */

@Description("Test the BirthDeathSerialSkylineModel with a small tree example")

public class BirthDeathSerialSkylineTest extends TestCase {

    @Test
    public void testLikelihoodCalculation() throws Exception {

        double PRECISION = 1e-12;
        RealParameter times = new RealParameter(new double[]{2.0});
        RealParameter birthRateVector = new RealParameter(new double[]{1/200., 1/300.});
        RealParameter deathRate = new RealParameter(new double[]{1/4., 1/7.});


        for (int i = 0; i <2; i++){

            Boolean r = (i!=0);

            BirthDeathSerialSkylineModel bdssm = new BirthDeathSerialSkylineModel();

            bdssm.setInputValue("times", times);
            bdssm.setInputValue("birthRateVector", birthRateVector);
            bdssm.setInputValue("birthRateScalar", 1.0);
            bdssm.setInputValue("deathRate", deathRate);

            bdssm.setInputValue("serialSamplingRate", 0.5);
            bdssm.setInputValue("extantSamplingRate", 0.5);
            bdssm.setInputValue("relativeDeath", false);

            bdssm.setInputValue("sampledIndividualsRemainInfectious", r);
            bdssm.setInputValue("finalTimeInterval", 0.);
            bdssm.setInputValue("origin", 5.0);

            bdssm.initAndValidate();

            Tree tree = new Tree("(((1:1,2:1):2,3:3):1,4:4);");
            TreeIntervals intervals = new TreeIntervals();
            intervals.init(tree);

    //        for (int i = 0; i < size; i += 1) {
    //            System.out.println("deathRate at time " + i + " is " + bdssm.deathRate(i));
    //            System.out.println("p0 at time " + i + " is " + bdssm.p0(i, i));
    //        }

            double logL = bdssm.calculateTreeLogLikelihood(tree);

            // parts of likelihood calculation
            assertEquals(0.5000000000000115, bdssm.p0(0,2,2), PRECISION); // p0_0_t1
            assertEquals(0.4111946067942296, bdssm.p0(0,1,0), PRECISION); // p0_0_x3
            assertEquals(0.5000000000000115, bdssm.p0(0,0,0), PRECISION); // p0_0_y
    //
            assertEquals(1.0, bdssm.g(0,0,0), PRECISION); // g_y
            assertEquals(0.1445844893653338, bdssm.g(1, 5, 2), PRECISION); // g_x0
            assertEquals(0.2754869293957173, bdssm.g(1, 4, 2), PRECISION); // g_x1
            assertEquals(0.5248860533202707, bdssm.g(1, 3, 2), PRECISION); // g_x2

            assertEquals(0.472130416276374, bdssm.g(0, 1, 0), PRECISION); // g_x3
            assertEquals(0.2227681214206327, bdssm.g(0, 2, 0), PRECISION); // g_0_t1

            //factors in likelihood formula
            assertEquals(0.00903653058533336, Math.exp(bdssm.TestFactor[0]), PRECISION); // T0
//            assertEquals(3.792761262897912e-09, Math.exp(bdssm.TestFactor[1]), PRECISION); // T1

            System.out.println("sampledIndividualsRemainInfectious = " + r);
            assertEquals((r?1.:0.06250000000000577), Math.exp(bdssm.TestFactor[2]), PRECISION); // T2
            assertEquals(0.01105500968848732, Math.exp(bdssm.TestFactor[3]), PRECISION); // T3



            // total likelihood
            assertEquals(1., (r ? 3.788928039364497e-13: 2.36808002460303e-14) / Math.exp(logL), PRECISION);
            System.out.println("loglikelihood: " + logL + " " + Math.exp(logL));
            System.out.println();
        }
    }


}


/* R code for test example

# calculate BDSSM treelikelihood for Tree("(((1:1,2:1):2,3:3):1,4:4);")


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