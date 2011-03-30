package beast.evolution.speciation;


import beast.evolution.tree.Tree;
import beast.evolution.alignment.Taxon;
import beast.core.parameter.*;
import beast.core.Input;
import beast.core.Description;

import java.util.*;

/**
 * @author Alexei Drummond, Denise KŸhnert
 */

@Description("Adaptation of Tanja Stadler's BirthDeathSerialSamplingModel, to allow for birth and death rates to change at times t_i")

public class BirthDeathSerialSkylineModel extends SpeciationLikelihood {

    public Input<RealParameter> times =
            new Input<RealParameter>("times", "times t_i specifying when rate changes can occur", Input.Validate.REQUIRED);
    public Input<RealParameter> birthRate =
            new Input<RealParameter>("birthRate", "birthRate vector with birthRates between times times", Input.Validate.REQUIRED);
    public Input<RealParameter> deathRate =
            new Input<RealParameter>("deathRate", "deathRate vector with birthRates between times times", Input.Validate.REQUIRED);
    public Input<Double> serialSamplingRate =
            new Input<Double>("serialSamplingRate", "serial sampling proportion", Input.Validate.REQUIRED);      // psi
    public Input<Double> extantSamplingRate =
            new Input<Double>("extantSamplingRate", "extant sampling proportion", Input.Validate.REQUIRED);

    public Input<Boolean> relativeDeath =
            new Input<Boolean>("relativeDeath", "is death relative?", Input.Validate.REQUIRED);

    public Input<Boolean> sampledIndividualsRemainInfectious =
            new Input<Boolean>("sampledIndividualsRemainInfectious", "stating whether sampled individuals remain infectious, or become non-infectious", Input.Validate.REQUIRED);

    public Input<Double> finalTimeInterval =
            new Input<Double>("finalTimeInterval", "final time interval", Input.Validate.REQUIRED);

    // the origin of the infection x0 which is greater than tree.getRoot();
    public Input<Double> origin =
            new Input<Double>("origin", "the origin of infection x0", Input.Validate.REQUIRED);

    // need TYPES?? 
//    public Input<String> m_pType =
//            new Input<String>("type", "tree type, should be one of " + Arrays.toString(TYPES)+" (default unscaled)",
//                    "unscaled", TYPES);

    protected double[] p0_iMinus1;
    protected double[] Ai;
    protected double[] Bi;
    protected int m;

    Boolean birthChanges;
    Boolean deathChanges;


    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        m = times.get().getDimension() + 1; // t0 separately: finalTimeInterval

        if (birthRate.get().getDimension() != 1 && birthRate.get().getDimension() != m)
            throw new RuntimeException("Length of birthRate parameter should be one or equal to the size of time parameter (size = " + (m+1) + ")");

        birthChanges = (birthRate.get().getDimension() == m) ;

        if (deathRate.get().getDimension() != 1 && deathRate.get().getDimension() != m)
            throw new RuntimeException("Length of mu parameter should be one or equal to the size of time parameter (size = " + (m+1) + ")");

        deathChanges = (deathRate.get().getDimension() == m) ;

    }

    /*    calculate and store Ai, Bi and p0_iMinus1        */
    public void preCalculation(double[] times){

        Ai = new double[m];
        Bi = new double[m];
        p0_iMinus1 = new double[m];

        for (int i = 0; i < m; i++){
            Ai[i] = Ai(birthRate.get().getArrayValue(birthChanges? i : 0), deathRate.get().getArrayValue(deathChanges? i : 0), serialSamplingRate.get());
        }

        Bi[0] = Bi(birthRate.get().getArrayValue(0), deathRate.get().getArrayValue(0), serialSamplingRate.get(), Ai[0], 1);
        for (int i = 1; i < m; i++){
            p0_iMinus1[i-1] = p0(birthRate.get().getArrayValue(birthChanges? (i-1) : 0), deathRate.get().getArrayValue(deathChanges? (i-1) : 0), serialSamplingRate.get(), Ai[i-1], Bi[i-1], times[i] , times[i]);
            Bi[i] = Bi(birthRate.get().getArrayValue(birthChanges? i : 0), deathRate.get().getArrayValue(deathChanges? i : 0), serialSamplingRate.get(), Ai[i], p0_iMinus1[i-1]);
        }
    }

    public double Ai(double b, double g, double psi){

        return Math.sqrt((b - g - psi)*(b - g - psi) + 4*b*psi);

    }

    public double Bi(double b, double g, double psi, double A, double p0){

        return (-((1-p0)*b + g + psi)/A);

    }

    public double p0(int index, double t, double ti){

        return p0(birthRate.get().getArrayValue(index), deathRate.get().getArrayValue(index), serialSamplingRate.get(), Ai[index], Bi[index], t, ti);

    }

    public double p0(double b, double g, double psi, double A, double B, double t, double ti){

        return ((b + g + psi - A *((Math.exp(A*(t - ti))*(1-B)-(1+B)))/(Math.exp(A*(t-ti))*(1-B)+(1+B)))/(2*b));

    }

    public double g(int index, double t, double ti){

        return (4 / (2*(1-Bi[index]*Bi[index]) + Math.exp(Ai[index]*(t - ti))*2*(1+Bi[index]*Bi[index])));

    }

    // get  array of interval times
    public double[] getEndTimes() {
        double[] endTimes = new double[m];
        endTimes[0] = finalTimeInterval.get();
        for (int i = 1; i < m; i++) {
            endTimes[i] = times.get().getArrayValue(i-1);
        }
        return endTimes;
    }

    // todo: make nice & efficient! 
    public int index(double t) {

        double[] endTime = getEndTimes();
        int epoch = Arrays.binarySearch(endTime, t);

        if (epoch < 0) {
            epoch = -epoch - 1;

        }

        return Math.max((epoch - 1), 0);
     }


    public double calculateTreeLogLikelihood(Tree tree) {

        double[] times = getEndTimes();
        double finalTime = times[0];

        // number of lineages at each time ti
        int[] n = new int[m];
        int nTips = tree.getLeafNodeCount();
        preCalculation(times);

        double x0 = tree.getRoot().getHeight() + finalTime;

        int index = m-1;      // x0 must be in last interval
        double print;

        // the first factor for origin
        print =nTips * Math.log(serialSamplingRate.get()) + Math.log(g(index, x0, times[index]));
        double logL = print;

//        double logL = nTips * Math.log(serialSamplingRate.get()) + Math.log(g(index, x0, times[index]));
        System.out.println("p0_0_t1: " + Math.exp(print));

         // first product term in f[T]
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {

            double x = tree.getNode(nTips+i).getHeight() + finalTime;

            index = index(x);
//            print = Math.log(birthRate.get().getArrayValue(index) * g(index, x, times[index]));
//            logL += print;
            logL += Math.log(birthRate.get().getArrayValue(index) * g(index, x, times[index]));

        }

         // middle product term in f[T]
        for (int i = 0; i < nTips; i++) {

            double y = tree.getNode(i).getHeight() + finalTime;

            index = index(y);


            // sampledIndividualsRemainInfectious is r_i, but we dont allow it to change here. can only be 1 or 0 for the whole time
            logL += Math.log( g(index, y, times[index]) ) + (sampledIndividualsRemainInfectious.get() ? Math.log(p0(index, y, times[index])) : 0);

        }

        // last product term in f[T], factorizing from 1 to m
        for (int j = 1; j < m; j++){
            double time = times[j];
            for (int i = 0; i < nTips; i++) {
                if (tree.getNode(i).getHeight() >= time)  n[j] += 1;
            }

            if (n[j] > 0) {
                logL += n[j] * Math.log(g(j-1, time, times[j]));     
            }
        }

        return logL;
    }

    public double calculateTreeLogLikelihood(Tree tree, Set<Taxon> exclude) {
        if (exclude.size() == 0) return calculateTreeLogLikelihood(tree);
        throw new RuntimeException("Not implemented!");
    }


}