package beast.evolution.speciation;


import beast.evolution.tree.Tree;
import beast.evolution.tree.Node;
import beast.evolution.alignment.Taxon;
import beast.core.parameter.*;
import beast.core.Input;
import beast.core.Description;

import java.util.*;

/**
 * @author Alexei Drummond, Denise Kuehnert
 */

@Description("Adaptation of Tanja Stadler's BirthDeathSerialSamplingModel, to allow for birth and death rates to change at times t_i")

public class BirthDeathSerialSkylineModel extends SpeciationLikelihood {

    public Input<RealParameter> times =
            new Input<RealParameter>("times", "The times t_i specifying when rate changes can occur", Input.Validate.REQUIRED);
    public Input<RealParameter> birthRateVector =
            new Input<RealParameter>("birthRateVector", "BirthRate = BirthRateVector * birthRateScalar, birthrate can change over time", Input.Validate.REQUIRED);
    public Input<Double> birthRateScalar =
                new Input<Double>("birthRateScalar", "BirthRate = BirthRateVector * birthRateScalar, birthrate can change over time", Input.Validate.REQUIRED);

    public Input<RealParameter> deathRate =
            new Input<RealParameter>("deathRate", "The deathRate vector with birthRates between times", Input.Validate.REQUIRED);
    public Input<Double> serialSamplingRate =
            new Input<Double>("serialSamplingRate", "The serial sampling proportion", Input.Validate.REQUIRED);      // psi
    public Input<Double> extantSamplingRate =
            new Input<Double>("extantSamplingRate", "The extant sampling proportion", Input.Validate.REQUIRED);

    public Input<Boolean> relativeDeath =
            new Input<Boolean>("relativeDeath", "Boolean, is death relative?", Input.Validate.REQUIRED);

    public Input<Boolean> sampledIndividualsRemainInfectious =
            new Input<Boolean>("sampledIndividualsRemainInfectious", "Boolean, stating whether sampled individuals remain infectious, or become non-infectious", Input.Validate.REQUIRED);

    public Input<Double> finalTimeInterval =
            new Input<Double>("finalTimeInterval", "The final time interval", Input.Validate.REQUIRED);

    // the origin of the infection x0 which is greater than tree.getRoot();
    public Input<Double> origin =
            new Input<Double>("origin", "The origin of infection x0", Input.Validate.REQUIRED);

    // need TYPES?? 
//    public Input<String> m_pType =
//            new Input<String>("type", "tree type, should be one of " + Arrays.toString(TYPES)+" (default unscaled)",
//                    "unscaled", TYPES);

    RealParameter birthRate;
    protected double[] p0_iMinus1;
    protected double[] Ai;
    protected double[] Bi;
    protected int m;

    Boolean birthChanges;
    Boolean deathChanges;

    //Testvector for BirthDeathSerialSkylineTest
    public double[] TestFactor = new double[4];

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        epiutil.ElementwiseMultiplication mult = new epiutil.ElementwiseMultiplication();
        mult.setInputValue("vector", birthRateVector.get());
        mult.setInputValue("scalar", birthRateScalar.get());
        birthRate = mult.get();


        m = times.get().getDimension() + 1; // t0 separately: finalTimeInterval

        if (birthRate.getDimension() != 1 && birthRate.getDimension() != m)
            throw new RuntimeException("Length of birthRate parameter should be one or equal to the size of time parameter (size = " + (m+1) + ")");

        birthChanges = (birthRate.getDimension() == m) ;

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
            Ai[i] = Ai(birthRate.getArrayValue(birthChanges? i : 0), deathRate.get().getArrayValue(deathChanges? i : 0), serialSamplingRate.get());
        }

        Bi[0] = Bi(birthRate.getArrayValue(0), deathRate.get().getArrayValue(0), serialSamplingRate.get(), Ai[0], 1);
        for (int i = 1; i < m; i++){
            p0_iMinus1[i-1] = p0(birthRate.getArrayValue(birthChanges? (i-1) : 0), deathRate.get().getArrayValue(deathChanges? (i-1) : 0), serialSamplingRate.get(), Ai[i-1], Bi[i-1], times[i] , times[i]);
            Bi[i] = Bi(birthRate.getArrayValue(birthChanges? i : 0), deathRate.get().getArrayValue(deathChanges? i : 0), serialSamplingRate.get(), Ai[i], p0_iMinus1[i-1]);
        }
    }

    public double Ai(double b, double g, double psi){

        return Math.sqrt((b - g - psi)*(b - g - psi) + 4*b*psi);

    }

    public double Bi(double b, double g, double psi, double A, double p0){

        return (-((1-p0)*b + g + psi)/A);

    }

    public double p0(int index, double t, double ti){

        double testp = p0(birthRate.getArrayValue(index), deathRate.get().getArrayValue(index), serialSamplingRate.get(), Ai[index], Bi[index], t, ti);
//        System.out.println("Calculating p0(" + index + ", " + t + ", " + ti + ") ... " + testp) ;

        return testp;

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

    public int nCurrentLineages(double time, Node root){
        int count = 0;
        if (root.getHeight() < time) return 0;
        count = getCurrentChildren(time, root, count);
        return count;
    }

    public int getCurrentChildren(double time, Node node, int count){

        if (node.m_left.getHeight() < time)
            count++;
        else count = getCurrentChildren(time, node.m_left, count);

        if (node.m_right.getHeight() < time)
            count++;
        else count = getCurrentChildren(time, node.m_left, count);

        return count;
    }


    public double calculateTreeLogLikelihood(Tree tree) {

        double[] times = getEndTimes();
        double finalTime = times[0];

        // number of lineages at each time ti
        int[] n = new int[m];
        int nTips = tree.getLeafNodeCount();
        preCalculation(times);

        double x0 = origin.get() + finalTime; //tree.getRoot().getHeight() + finalTime;

        int index = m-1;      // x0 must be in last interval

        // the first factor for origin
        TestFactor[0] = nTips * Math.log(serialSamplingRate.get()) + Math.log(g(index, x0, times[index]));
        double logL = nTips * Math.log(serialSamplingRate.get()) + Math.log(g(index, x0, times[index]));

        // first product term in f[T]
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {

            double x = tree.getNode(nTips+i).getHeight() + finalTime;

            index = index(x);

            TestFactor[1] += Math.log(birthRate.getArrayValue(index) * g(index, x, times[index]));

            logL += Math.log(birthRate.getArrayValue(index) * g(index, x, times[index]));

        }

        // middle product term in f[T]
        for (int i = 0; i < nTips; i++) {

            double y = tree.getNode(i).getHeight() + finalTime;

            index = index(y);

            // sampledIndividualsRemainInfectious is r_i, but we dont allow it to change here. can only be 1 or 0 for the whole time
            TestFactor[2] += (sampledIndividualsRemainInfectious.get() ? 0 : Math.log(p0(index, y, times[index]))) - Math.log(g(index, y, times[index]) ) ;

            logL += (sampledIndividualsRemainInfectious.get() ? 0 : Math.log(p0(index, y, times[index]))) - Math.log(g(index, y, times[index]) ) ;

        }

        // last product term in f[T], factorizing from 1 to m
        for (int j = 0; j < m-1; j++){
            double time = times[j+1];
            n[j] = nCurrentLineages(time, tree.getRoot());
            if (n[j] > 0) {
                TestFactor[3] +=  n[j] * Math.log(g(j, time, times[j]));
                logL +=  n[j] * Math.log(g(j, time, times[j]));
            }
        }
        logP = logL;
        return logL;
    }

    public double calculateTreeLogLikelihood(Tree tree, Set<Taxon> exclude) {
        if (exclude.size() == 0) return calculateTreeLogLikelihood(tree);
        throw new RuntimeException("Not implemented!");
    }


}