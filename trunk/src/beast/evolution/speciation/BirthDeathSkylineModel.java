package beast.evolution.speciation;


import beast.evolution.tree.Tree;
import beast.evolution.alignment.Taxon;
import beast.core.parameter.*;
import beast.core.Input;
import beast.core.Description;

import java.util.*;

/**
 * @author Denise Kuehnert
 *
 * maths: Tanja Stadler
 */

@Description("Adaptation of Tanja Stadler's BirthDeathSamplingModel, to allow for birth and death rates to change at times t_i")
public class BirthDeathSkylineModel extends SpeciesTreeDistribution {

    // assume equidistant intervals if intervaltimes are not specified
    public Input<RealParameter> intervalTimes =
            new Input<RealParameter>("intervalTimes", "The times t_i specifying when rate changes can occur", (RealParameter) null);
    public Input<RealParameter> orig_root =
            new Input<RealParameter>("orig_root", "The origin of infection x0", Input.Validate.REQUIRED);
    public Input<IntegerParameter> intervalNumber =
            new Input<IntegerParameter>("intervalNumber", "The number of intervals in which rates can change", Input.Validate.REQUIRED);

    public Input<RealParameter> birthRate =
            new Input<RealParameter>("birthRate", "BirthRate = BirthRateVector * birthRateScalar, birthrate can change over time", Input.Validate.REQUIRED);
    public Input<RealParameter> deathRate =
            new Input<RealParameter>("deathRate", "The deathRate vector with birthRates between times", Input.Validate.REQUIRED);
    public Input<RealParameter> samplingRate =
            new Input<RealParameter>("samplingRate", "The sampling rate per individual", Input.Validate.REQUIRED);      // psi

    public Input<Boolean> forceRateChange =
            new Input<Boolean>("forceRateChange", "If there is more than one interval and we estimate the time of rate change, do we enforce it to be within the tree interval? Default true", true);

    double t_root;
    protected double[] p0_iMinus1;
    protected double[] Ai;
    protected double[] Bi;
    protected int m;


    Double[] birth;
    Double[] death;
    Double[] psi;
    Boolean birthChanges;
    Boolean deathChanges;
    Boolean samplingChanges;
    Double[] times;
    Boolean timesFromXML;

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        m = intervalNumber.get().getValue();

        death = deathRate.get().getValues();
        psi = samplingRate.get().getValues();
        birth = new Double[m];

        if (birth.length != 1 && birth.length != m)
            throw new RuntimeException("Length of birthRate parameter ("+ birth.length + ") should be one or equal to intervalNumber (" + m + ")");

        birthChanges = (birthRate.get().getDimension() == m && m > 1)  ;

        if (death.length != 1 && death.length != m)
            throw new RuntimeException("Length of mu parameter should be one or equal to intervalNumber (" + m + ")");

        deathChanges = (death.length == m && m > 1) ;

        if (psi.length != 1 && psi.length != m)
            throw new RuntimeException("Length of birthRate parameter ("+ psi.length + ") should be one or equal to intervalNumber (" + m + ")");

        samplingChanges = (psi.length == m && m > 1) ;

        if (intervalTimes.get() != null){
            if (intervalTimes.get().getDimension() != m)
                throw new RuntimeException("Length of intervalTimes parameter should equal to intervalNumber (" + m + ")");

            if (intervalTimes.get().getValue() != 0)
                throw new RuntimeException("First entry in intervalTimes must be 0!");
            timesFromXML = true;
        }
        else {
            times  = new Double[m];
            for (int i = 0; i < m; i++) {
                timesFromXML = false;
            }
        }
    }

    /*    calculate and store Ai, Bi and p0_iMinus1        */
    public Double preCalculation(){

        death = deathRate.get().getValues();
        birth = birthRate.get().getValues();

        t_root = m_tree.get().getRoot().getHeight();

        if (timesFromXML){
            times = intervalTimes.get().getValues();

            // make sure intervalTimes are increasing
            for (int i = 1; i < times.length; i++){
                if (times[i] <= times[i-1]) return Double.NEGATIVE_INFINITY;
            }

            // if forceRateChange: force rate change to be within tree range
            if (forceRateChange.get()){
                if ( times[times.length-1] >= (t_root + orig_root.get().getValue()) )
                    return Double.NEGATIVE_INFINITY;
            }


            Double[] temp = new Double[m];
            temp[0] = 0.;
            for (int i = 1; i < m; i++) {
                temp[i] = Math.max(t_root + orig_root.get().getValue() - times[m-i], 0);
            }
            times = temp;

        } else {
            times  = new Double[m];
            for (int i = 0; i < m; i++) {
                times[i] = i * (t_root+orig_root.get().getValue()) / m;
            }
        }


        // changing parametrization to have everything forward..
        for(int i = 0; i < birth.length / 2; i++)
        {
            double tempBirth = birth[i];
            birth[i] = birth[birth.length - i - 1];
            birth[birth.length - i - 1] = tempBirth;

            if (i < death.length / 2){
                double tempDeath = death[i];
                death[i] = death[death.length - i - 1];
                death[death.length - i - 1] = tempDeath;
            }

            if (i < psi.length / 2){
                double tempPsi = psi[i];
                psi[i] = psi[psi.length - i - 1];
                psi[psi.length - i - 1] = tempPsi;
            }
        }


        Ai = new double[m];
        Bi = new double[m];
        p0_iMinus1 = new double[m];

        for (int i = 0; i < m; i++){
            Ai[i] = Ai(birth[birthChanges? i : 0], death[deathChanges? i : 0], psi[samplingChanges? i : 0]);
//            System.out.println("Ai[" + i + "] = " + Ai[i] + " " +  Math.log(Ai[i]));
        }

        Bi[0] = Bi(birth[0], death[0], psi[0], Ai[0], 1);  //  (p0_iMinus1[-1] = 1)
//        System.out.println("B0[0] = " + Bi[0] + " " + Math.log(Bi[0]));
        for (int i = 1; i < m; i++){
            p0_iMinus1[i-1] = p0(birth[birthChanges? (i-1) : 0], death[deathChanges? (i-1) : 0], psi[samplingChanges? (i-1) : 0], Ai[i-1], Bi[i-1], times[i] , times[i-1]);
            Bi[i] = Bi(birth[birthChanges? i : 0], death[deathChanges? i : 0], psi[samplingChanges? i : 0], Ai[i], p0_iMinus1[i-1]);
//            System.out.println("Bi[" + i + "] = " + Bi[i] + " " +  Math.log(Bi[i]));

        }

        return 0.;
    }

    public double Ai(double b, double g, double psi){

        return Math.sqrt((b - g - psi)*(b - g - psi) + 4*b*psi);
    }

    public double Bi(double b, double g, double psi, double A, double p0){

        return (-((1-2*p0)*b + g + psi)/A);
    }

    public double p0(int index, double t, double ti){

        return p0(birth[birthChanges? index : 0], death[deathChanges?index:0], psi[samplingChanges? index : 0], Ai[index], Bi[index], t, ti);
    }

    public double p0(double b, double g, double psi, double A, double B, double t, double ti){

        return ((b + g + psi - A *((Math.exp(A*(t - ti))*(1-B)-(1+B)))/(Math.exp(A*(t-ti))*(1-B)+(1+B)) ) / (2*b));
    }

    public double g(int index, double t, double ti){

        return (4 * Math.exp(-Ai[index]*(t - ti))) / Math.pow((Math.exp(-Ai[index]*(t - ti))*(1+Bi[index])+(1-Bi[index])),2);
    }

    public int index(double t) {

        int epoch = Arrays.binarySearch(times, t);

        if (epoch < 0) {
            epoch = -epoch - 1;
        }

        return Math.max((epoch - 1), 0);
    }


    /**
     * @param time the time
     * @param tree the tree
     * @return the number of lineages that exist at the given time in the given tree.
     */
    public int lineageCountAtTime(double time, Tree tree) {

        int count = 1;
        int tipCount = tree.getLeafNodeCount();
        for (int i = tipCount; i < tipCount + tree.getInternalNodeCount(); i++) {
            if (tree.getNode(i).getHeight() > time) count += 1;

        }
        for (int i = 0; i < tipCount; i++) {
            if (tree.getNode(i).getHeight() > time) count -= 1;
        }
        return count;
    }

    public double calculateTreeLogLikelihood(Tree tree) {

        m = intervalNumber.get().getValue();

        // number of lineages at each time ti
        int[] n = new int[m];
        int nTips = tree.getLeafNodeCount();
        if (preCalculation() < 0)
            return Double.NEGATIVE_INFINITY;

        double x0 = tree.getRoot().getHeight() + orig_root.get().getValue(); //tree.getRoot().getHeight() ;

        int index = m-1;      // x0 must be in last interval
        double temp;

        // the first factor for origin
        temp =  Math.log(g(index, x0, times[index]));
        logP = temp;
//        System.out.println("orig = " + temp);


        // first product term in f[T]
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {

            double x = tree.getNode(nTips+i).getHeight();
            index = index(x);

            temp = Math.log(birth[birthChanges? index:0] * g(index, x, times[index]));
            logP += temp;
//            System.out.println(temp);
        }

        // middle product term in f[T]
        for (int i = 0; i < nTips; i++) {

            double y = tree.getNode(i).getHeight();
            index = index(y);

            temp = Math.log(psi[samplingChanges? index : 0]) - Math.log(g(index, y, times[index]) ) ;
            logP += temp;
//            System.out.println(temp);
        }

        // last product term in f[T], factorizing from 1 to m
        for (int j = 0; j < m-1; j++){
            double time = times[j+1];
            n[j] = lineageCountAtTime(time, tree);

            if (n[j] > 0) {
                temp =  n[j] * Math.log(g(j, time, times[j]));
                logP += temp;
//                System.out.println(temp);
            }
        }
        return logP;
    }

    public double calculateTreeLogLikelihood(Tree tree, Set<Taxon> exclude) {
        if (exclude.size() == 0) return calculateTreeLogLikelihood(tree);
        throw new RuntimeException("Not implemented!");
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }


    public static int index(Double t, Double[] timeArray) {

        int epoch = Arrays.binarySearch(timeArray, t);

        if (epoch < 0) {
            epoch = -epoch - 1;
        }

        return Math.max((epoch - 1), 0);
    }


    public static void main(String[] args){

        Double[] times = new Double[]{0., 0., 2.};
        double t = 1.;

        System.out.println(index(t, times));
    }

}