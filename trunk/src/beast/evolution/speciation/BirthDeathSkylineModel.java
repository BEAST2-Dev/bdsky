package beast.evolution.speciation;


import beast.evolution.tree.Tree;
import beast.evolution.alignment.Taxon;
import beast.core.parameter.*;
import beast.core.Input;
import beast.core.Description;
import beast.core.util.ParameterConstrainer;

import java.util.*;

/**
 * @author Denise Kuehnert
 * @author Alexei Drummond
 *         <p/>
 *         maths: Tanja Stadler
 */

@Description("Adaptation of Tanja Stadler's BirthDeathSamplingModel, to allow for birth and death rates to change at times t_i")
public class BirthDeathSkylineModel extends SpeciesTreeDistribution {

    // the interval times for the birth rate
    public Input<RealParameter> birthRateChangeTimesInput =
            new Input<RealParameter>("birthRateChangeTimes", "The times t_i specifying when birth/R rate changes occur", (RealParameter) null);

    // the interval times for the death rate
    public Input<RealParameter> deathRateChangeTimesInput =
            new Input<RealParameter>("deathRateChangeTimes", "The times t_i specifying when death/becomeUninfectious rate changes occur", (RealParameter) null);

    // the interval times for sampling rate
    public Input<RealParameter> samplingRateChangeTimesInput =
            new Input<RealParameter>("samplingRateChangeTimes", "The times t_i specifying when sampling rate or sampling proportion changes occur", (RealParameter) null);

    public Input<Boolean> birthRateChangeTimesRelativeInput =
            new Input<Boolean>("birthRateTimesRelative", "True if birth rate change times specified relative to tree height? Default true", true);

    public Input<Boolean> deathRateChangeTimesRelativeInput =
            new Input<Boolean>("deathRateTimesRelative", "True if death rate change times specified relative to tree height? Default true", true);

    public Input<Boolean> samplingRateChangeTimesRelativeInput =
            new Input<Boolean>("samplingRateTimesRelative", "True if sampling rate times specified relative to tree height? Default true", true);

    // the times for rho sampling
    public Input<RealParameter> rhoSamplingTimes =
            new Input<RealParameter>("rhoSamplingTimes", "The times t_i specifying when rho-sampling occurs", (RealParameter) null);


    public Input<RealParameter> orig_root =
            new Input<RealParameter>("orig_root", "The origin of infection x0", Input.Validate.REQUIRED);

    public Input<RealParameter> birthRate =
            new Input<RealParameter>("birthRate", "BirthRate = BirthRateVector * birthRateScalar, birthrate can change over time");
    public Input<RealParameter> deathRate =
            new Input<RealParameter>("deathRate", "The deathRate vector with birthRates between times");
    public Input<RealParameter> samplingRate =
            new Input<RealParameter>("samplingRate", "The sampling rate per individual");      // psi

    public Input<RealParameter> m_rho =
            new Input<RealParameter>("rho", "The proportion of lineages sampled at rho-sampling times (default 0.)");
    public Input<Boolean> contemp =
            new Input<Boolean>("contemp", "Only contemporaneous sampling (i.e. all tips are from same sampling time, default false)", false);


    public Input<RealParameter> R0 =
            new Input<RealParameter>("R0", "The basic reproduction number", Input.Validate.XOR, birthRate);
    public Input<RealParameter> becomeUninfectiousRate =
            new Input<RealParameter>("becomeUninfectiousRate", "Rate at which individuals become uninfectious (throuch recovery or sampling)", Input.Validate.XOR, deathRate);
    public Input<RealParameter> samplingProportion =
            new Input<RealParameter>("samplingProportion", "The samplingProportion = samplingRate / becomeUninfectiousRate", Input.Validate.XOR, samplingRate);


    public Input<Boolean> forceRateChange =
            new Input<Boolean>("forceRateChange", "If there is more than one interval and we estimate the time of rate change, do we enforce it to be within the tree interval? Default true", true);
    public Input<Boolean> conditionOnSurvival =
            new Input<Boolean>("conditionOnSurvival", "condition on at least one survival? Default true.", true);

    public Input<ParameterConstrainer> psiConstrainer = new Input<ParameterConstrainer>("psiConstrainer", "psiConstrainer, constrain samplingRate to period of sampling");

    public Input<IntegerParameter> S0_input =
            new Input<IntegerParameter>("S0", "The numbers of susceptible individuals");

    double t_root;
    protected double[] p0;
    protected double[] Ai;
    protected double[] Bi;
    protected int[] N;   // number of leaves sampled at each time t_i

    // these four arrays are totalIntervals in length
    Double[] birth;
    Double[] death;
    Double[] psi;
    Double[] rho;

    // true if the node of the given index occurs at the time of a rho-sampling event
    boolean[] isRhoTip;

    /**
     * The number of change points in the birth rate
     */
    int birthChanges;

    /**
     * The number of change points in the death rate
     */
    int deathChanges;

    /**
     * The number of change points in the sampling rate
     */
    int samplingChanges;

    /**
     * The number of times rho-sampling occurs
     */
    int rhoSamplingCount;

    /**
     * Total interval count
     */
    int totalIntervals;

    private List<Double> birthRateChangeTimes = new ArrayList<Double>();
    private List<Double> deathRateChangeTimes = new ArrayList<Double>();
    private List<Double> samplingRateChangeTimes = new ArrayList<Double>();

    /**
     * map of birth rate parameters keyed by times in timesSet
     */
//    private Map<Double, Double> birthRates = new HashMap<Double, Double>();
//    private Map<Double, Double> deathRates = new HashMap<Double, Double>();
//    private Map<Double, Double> samplingRates = new HashMap<Double, Double>();

    Boolean contempData;
    //List<Interval> intervals = new ArrayList<Interval>();
    SortedSet<Double> timesSet = new TreeSet<Double>();

    Double[] times;

    Boolean transform;
    Boolean m_forceRateChange;

    Boolean birthRateTimesRelative = true;
    Boolean deathRateTimesRelative = true;
    Boolean samplingRateTimesRelative = true;


    Boolean printTempResults;

/************************************************************************************************/
    /*              "constructor"                                                                   */

    /**
     * ********************************************************************************************
     */

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        m_forceRateChange = forceRateChange.get();
        birthRateTimesRelative = birthRateChangeTimesRelativeInput.get();
        deathRateTimesRelative = deathRateChangeTimesRelativeInput.get();
        samplingRateTimesRelative = samplingRateChangeTimesRelativeInput.get();

        contempData = contemp.get();
        rhoSamplingCount = 0;

        collectTimes();

//        if (birthRate.get() != null && deathRate.get() != null && samplingRate.get() != null) {
//
//            transform = false;
//            death = deathRate.get().getValues();
//            psi = samplingRate.get().getValues();
//            birth = birthRate.get().getValues();
//
//        } else if (R0.get() != null && becomeUninfectiousRate.get() != null && samplingProportion.get() != null) {
//
//            transform = true;
//            transformParameters(S0_input.get() == null ? 1 : S0_input.get().getValue());
//
//        } else {
//            throw new RuntimeException("Either specify birthRate, deathRate and samplingRate OR specify R0, becomeUninfectiousRate and samplingProportion!");
//        }
//
//        if (transform) {
//            if (R0.get().getDimension() != birthRateChangeTimesInput.get().getDimension())
//                throw new RuntimeException("Length of R0 parameter should be equal to length of interval times (" +
//                        birthRateChangeTimesInput.get().getDimension() + ")");
//
//
//            if (becomeUninfectiousRate.get().getDimension() != deathRateChangeTimesInput.get().getDimension())
//                throw new RuntimeException("Length of becomeUninfectiousRate parameter should be equal to intervalNumber (" +
//                        deathRateChangeTimesInput.get().getDimension() + ")");
//
//            if (samplingProportion.get().getDimension() != samplingRateChangeTimesInput.get().getDimension())
//                throw new RuntimeException("Length of samplingProportion parameter should be equal to intervalNumber ("
//                        + samplingRateChangeTimesInput.get().getDimension() + ")");
//
//            birthChanges = R0.get().getDimension() - 1;
//            samplingChanges = samplingProportion.get().getDimension() - 1;
//            deathChanges = becomeUninfectiousRate.get().getDimension() - 1;
//        } else {
//
//            birthChanges = birthRate.get().getDimension() - 1;
//            deathChanges = deathRate.get().getDimension() - 1;
//            samplingChanges = samplingRate.get().getDimension() - 1;
//        }

        if (m_rho.get() != null) {

            if (m_rho.get().getDimension() == 1) {
                if (!contempData && ((samplingProportion.get() != null && samplingProportion.get().getDimension() == 1 && samplingProportion.get().getValue() == 0.) ||
                        (samplingRate.get() != null && samplingRate.get().getDimension() == 1 && samplingRate.get().getValue() == 0.))) {
                    contempData = true;
                    System.out.println("Parameters were chosen for contemporaneously sampled data. Setting contemp=true.");
                }
            }

            if (contempData) {
                if (m_rho.get().getDimension() != 1)
                    throw new RuntimeException("when contemp=true, rho must have dimension 1");

                else {
                    rho = new Double[totalIntervals];
                    Arrays.fill(rho, 0.);
                    rho[totalIntervals - 1] = m_rho.get().getValue();
                    rhoSamplingCount = 1;
                }
            } else {

                rho = new Double[totalIntervals];

                RealParameter rhoSampling = rhoSamplingTimes.get();

                for (int i = 0; i < rhoSampling.getDimension(); i++) {
                    rho[index(rhoSampling.getValue(i))] = m_rho.get().getValue(i);
                }
                rhoSamplingCount = rho.length;
            }
        } else {
            rho = new Double[totalIntervals];
            Arrays.fill(rho, 0.);
        }
        isRhoTip = new boolean[m_tree.get().getLeafNodeCount()];


        printTempResults = false;
    }

    /**
     * @return a list of intervals
     */
    public void getChangeTimes(List<Double> changeTimes, RealParameter intervalTimes, int numChanges, boolean relative) {
        changeTimes.clear();

        double maxTime = m_tree.get().getRoot().getHeight() + orig_root.get().getValue();

        if (intervalTimes == null) {
            //equidistant

            double intervalWidth = maxTime / (numChanges + 1);

            double start = 0.0;
            double end;
            for (int i = 1; i <= numChanges; i++) {
                end = (intervalWidth) * i;
                changeTimes.add(start);
                start = end;
            }
            changeTimes.add(start);
            //end = maxTime;

        } else {

            if (intervalTimes.getValue(0) != 0.0) {
                throw new RuntimeException("First time in interval times parameter should always be zero.");
            }

            if (numChanges > 0 && intervalTimes.getDimension() != numChanges + 1) {
                throw new RuntimeException("The time interval parameter should be numChanges + 1 long (" + (numChanges + 1));
            }


            double start = 0.0;
            double end;
            for (int i = 1; i < intervalTimes.getDimension(); i++) {
                end = intervalTimes.getValue(i);
                if (relative) end *= maxTime;
                changeTimes.add(start);
                start = end;
            }
            changeTimes.add(start);
            //end = maxTime;
        }
    }

    /*
    * Counts the number of tips at each of the contemporaneous sampling times ("rho" sampling time)
    * @return negative infinity if tips are found at a time when rho is zero, zero otherwise.
    */
    private double computeN(Tree tree) {

        isRhoTip = new boolean[tree.getLeafNodeCount()];

        N = new int[totalIntervals];

        int tipCount = tree.getLeafNodeCount();

        double[] dates = new double[tipCount];

        for (int i = 0; i < tipCount; i++) {
            dates[i] = tree.getNode(i).getDate();
        }

        for (int k = 0; k < totalIntervals; k++) {


            for (int i = 0; i < tipCount; i++) {

                if (Math.abs((times[totalIntervals - 1] - times[k]) - dates[i]) < 1e-10) {
                    if (rho[k] == 0)
                        return Double.NEGATIVE_INFINITY;
                    N[k] += 1;
                    isRhoTip[i] = true;
                }
            }
        }
        return 0.;
    }

    /**
     * Collect all the times of parameter value changes and rho-sampling events
     */
    private void collectTimes() {

        getChangeTimes(birthRateChangeTimes, birthRateChangeTimesInput.get(), birthChanges, true);
        getChangeTimes(deathRateChangeTimes, deathRateChangeTimesInput.get(), deathChanges, true);
        getChangeTimes(samplingRateChangeTimes, samplingRateChangeTimesInput.get(), samplingChanges, true);

        //eventsSet.clear();
        for (Double time : birthRateChangeTimes) {
            //eventsSet.add(new BDSEvent(BDSEvent.Type.birth,time));
            timesSet.add(time);
        }
        for (Double time : deathRateChangeTimes) {
            //eventsSet.add(new BDSEvent(BDSEvent.Type.death,time));
            timesSet.add(time);

        }
        for (Double time : samplingRateChangeTimes) {
            //eventsSet.add(new BDSEvent(BDSEvent.Type.sampling,time));
            timesSet.add(time);
        }

        RealParameter rhoSampling = rhoSamplingTimes.get();
        for (int i = 0; i < rhoSampling.getDimension(); i++) {
            //eventsSet.add(new BDSEvent(BDSEvent.Type.rhoSampling, rhoSampling.getValue(i)));
            timesSet.add(rhoSampling.getValue(i));
        }

        times = timesSet.toArray(times);
        totalIntervals = times.length;
    }

    private Double updateRatesAndTimes(Tree tree) {

        collectTimes();

        t_root = tree.getRoot().getHeight();

        if (m_forceRateChange && timesSet.last() > t_root) {
            return Double.NEGATIVE_INFINITY;
        }

        if (transform)
            transformParameters(S0_input.get() == null ? 1 : S0_input.get().getValue());
        else {

            Double[] birthRates = birthRate.get().getValues();
            Double[] deathRates = deathRate.get().getValues();
            Double[] samplingRates = samplingRate.get().getValues();

            for (int i = 0; i < totalIntervals; i++) {
                birth[i] = birthRates[index(times[i], birthRateChangeTimes)];
                death[i] = deathRates[index(times[i], deathRateChangeTimes)];
                psi[i] = samplingRates[index(times[i], samplingRateChangeTimes)];

            }
        }

        Double[] rhos = m_rho.get().getValues();
        RealParameter rhoSampling = rhoSamplingTimes.get();

        for (int i = 0; i < totalIntervals; i++) {
            for (int j = 0; i < rhos.length; i++) {
                if (times[i].equals(rhoSampling.getValue(j))) rho[i] = rhos[j];
            }
        }

        return 0.;
    }

/************************************************************************************************/
    /*                   calculations for likelihood                                                */

    /**
     * ********************************************************************************************
     */

    /*    calculate and store Ai, Bi and p0        */
    public Double preCalculation(Tree tree) {

        // updateRatesAndTimes must be called before calls to index() below
        if (updateRatesAndTimes(tree) < 0) return Double.NEGATIVE_INFINITY;

        if (m_rho.get() != null) {
            if (contempData) {
                Arrays.fill(rho, 0.);
                rho[totalIntervals - 1] = m_rho.get().getValue();
            } else {
                rho = new Double[totalIntervals];

                RealParameter rhoSampling = rhoSamplingTimes.get();

                for (int i = 0; i < rhoSampling.getDimension(); i++) {
                    rho[index(rhoSampling.getValue(i))] = m_rho.get().getValue(i);
                }
                rhoSamplingCount = rho.length;
            }
        } else rho = new Double[]{0.};

        if (m_rho.get() != null)
            if (computeN(tree) < 0)
                return Double.NEGATIVE_INFINITY;  // todo: check if it's enough to do this once at the beginning (only if interval times and sample dates don't change)

        int intervalCount = times.length;

        Ai = new double[intervalCount];
        Bi = new double[intervalCount];
        p0 = new double[intervalCount];

        for (int i = 0; i < intervalCount; i++) {

            Ai[i] = Ai(birth[i], death[i], psi[i]);

            if (printTempResults) System.out.println("Ai[" + i + "] = " + Ai[i] + " " + Math.log(Ai[i]));
        }

        Bi[totalIntervals - 1] = Bi(
                birth[totalIntervals - 1],
                death[totalIntervals - 1],
                psi[totalIntervals - 1],
                rho[totalIntervals - 1],
                Ai[totalIntervals - 1], 1.);  //  (p0[m-1] = 1)

        if (printTempResults)
            System.out.println("Bi[m-1] = " + Bi[totalIntervals - 1] + " " + Math.log(Bi[totalIntervals - 1]));
        for (int i = totalIntervals - 2; i >= 0; i--) {

            p0[i + 1] = p0(birth[i], death[i], psi[i], Ai[i + 1], Bi[i + 1], times[i + 1], times[i]);
            if (Math.abs(p0[i + 1] - 1) < 1e-10) {
                return Double.NEGATIVE_INFINITY;
            }
            if (printTempResults) System.out.println("p0[" + (i + 1) + "] = " + p0[i + 1]);

            Bi[i] = Bi(birth[i], death[i], psi[i], rho[i], Ai[i], p0[i + 1]);

            if (printTempResults) System.out.println("Bi[" + i + "] = " + Bi[i] + " " + Math.log(Bi[i]));
        }

        if (printTempResults) {
            System.out.println("g(0, x0, 0):" + g(0, times[0], 0));
            System.out.println("g(index(1),times[index(1)],1.) :" + g(index(1), times[index(1)], 1.));
            System.out.println("g(index(2),times[index(2)],2.) :" + g(index(2), times[index(2)], 2));
            System.out.println("g(index(4),times[index(4)],4.):" + g(index(4), times[index(4)], 4));

        }

        return 0.;
    }

    public double Ai(double b, double g, double psi) {

        return Math.sqrt((b - g - psi) * (b - g - psi) + 4 * b * psi);
    }

    public double Bi(double b, double g, double psi, double r, double A, double p0) {

        return ((1 - 2 * p0 * (1 - r)) * b + g + psi) / A;
    }

    public double p0(int index, double t, double ti) {

        return p0(birth[index], death[index], psi[index], Ai[index], Bi[index], t, ti);
    }

    public double p0(double b, double g, double psi, double A, double B, double ti, double t) {

        if (printTempResults)
            System.out.println("in p0: b = " + b + "; g = " + g + "; psi = " + psi + "; A = " + A + " ; B = " + B + "; ti = " + ti + "; t = " + t);
//        return ((b + g + psi - A *((Math.exp(A*(ti - t))*(1+B)-(1-B)))/(Math.exp(A*(ti - t))*(1+B)+(1-B)) ) / (2*b));
        // formula from manuscript slightly rearranged for numerical stability
        return ((b + g + psi - A * ((1 + B) - (1 - B) * (Math.exp(A * (t - ti)))) / ((1 + B) + Math.exp(A * (t - ti)) * (1 - B))) / (2 * b));

    }

    public double g(int index, double ti, double t) {

//        return (Math.exp(Ai[index]*(ti - t))) / (0.25*Math.pow((Math.exp(Ai[index]*(ti - t))*(1+Bi[index])+(1-Bi[index])),2));
        // formula from manuscript slightly rearranged for numerical stability
        return (4 * Math.exp(Ai[index] * (t - ti))) / (Math.exp(Ai[index] * (t - ti)) * (1 - Bi[index]) + (1 + Bi[index])) / (Math.exp(Ai[index] * (t - ti)) * (1 - Bi[index]) + (1 + Bi[index]));
    }

    /**
     * @param t the time in question
     * @return the index of the given time in the event list, or if the time is not in the event list, the index of the event
     *         with the next smallest time
     */
//    public int index(double t, List<BDSEvent> events) {
//
//        int epoch = Collections.binarySearch(events, new BDSEvent(events.get(0).type,t));
//
//        if (epoch < 0) {
//            epoch = -epoch - 1;
//        }
//
//        return epoch;
//    }

    /**
     * @param t the time in question
     * @return the index of the given time in the list of times, or if the time is not in the list, the index of the
     *         next smallest time
     */
    public int index(double t, List<Double> times) {

        int epoch = Collections.binarySearch(times, t);

        if (epoch < 0) {
            epoch = -epoch - 1;
        }

        return epoch;
    }


    /**
     * @param t the time in question
     * @return the index of the given time in the times array, or if the time is not in the array the index of the time
     *         next smallest
     */
    public int index(double t) {

        int epoch = Arrays.binarySearch(times, t);

        if (epoch < 0) {
            epoch = -epoch - 1;
        }

        // TODO this minimum should not be activated unless the time is greater than the root of the tree?
        return Math.min(epoch, totalIntervals - 1); //Math.max((epoch - 1), 0);
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
            if (tree.getNode(i).getHeight() >= time) count -= 1;
        }
        return count;
    }

    private void transformParameters(int S0) {

        Double[] R = R0.get().getValues();
        Double[] b = becomeUninfectiousRate.get().getValues();
        Double[] p = samplingProportion.get().getValues();

        for (int i = 0; i < totalIntervals; i++) {
            birth[i] = R[index(times[i], birthRateChangeTimes)] * b[index(times[i], deathRateChangeTimes)] / S0;
            psi[i] = p[index(times[i], samplingRateChangeTimes)] * b[index(times[i], deathRateChangeTimes)];
            death[i] = b[index(times[i], deathRateChangeTimes)] - psi[i];

        }
    }

    @Override
    public double calculateTreeLogLikelihood(Tree tree) {

        int nTips = tree.getLeafNodeCount();

        if (preCalculation(tree) < 0)
            return Double.NEGATIVE_INFINITY;

        // number of lineages at each time ti
        int[] n = new int[totalIntervals];

        double x0 = 0;
        int index = 0;

        double temp;

        // the first factor for origin
        if (!conditionOnSurvival.get())
            temp = Math.log(g(index, times[index], x0));  // NOT conditioned on at least one sampled individual
        else {
            temp = p0(index, times[index], x0);
            if (temp == 1)
                return Double.NEGATIVE_INFINITY;
            temp = Math.log(g(index, times[index], x0) / (1 - temp));   // DEFAULT: conditioned on at least one sampled individual
        }

        logP = temp;
        if (Double.isInfinite(logP))
            return logP;

        if (printTempResults) System.out.println("first factor for origin = " + temp);

        // first product term in f[T]
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {

            double x = times[totalIntervals - 1] - tree.getNode(nTips + i).getHeight();
            index = index(x);

            temp = Math.log(birth[index] * g(index, times[index], x));
            logP += temp;
            if (printTempResults) System.out.println("1st pwd" +
                    " = " + temp + "; interval = " + i);
            if (Double.isInfinite(logP))
                return logP;

        }

        // middle product term in f[T]
        for (int i = 0; i < nTips; i++) {

            if (!isRhoTip[i] || m_rho.get() == null) {
                double y = times[totalIntervals - 1] - tree.getNode(i).getHeight();
                index = index(y);

                temp = Math.log(psi[index]) - Math.log(g(index, times[index], y));
                logP += temp;
                if (printTempResults) System.out.println("2nd PI = " + temp);
                if (Double.isInfinite(logP))
                    return logP;

            }
        }

        // last product term in f[T], factorizing from 1 to m
        double time;
        for (int j = 0; j < totalIntervals; j++) {
            time = j < 1 ? 0 : times[j - 1];
            n[j] = (j == (0) ? 0 : lineageCountAtTime(times[totalIntervals - 1] - time, tree));

            if (n[j] > 0) {
                temp = n[j] * Math.log(g(j, times[j], time));
                logP += temp;
                if (printTempResults)
                    System.out.println("3rd factor (nj loop) = " + temp + "; interval = " + j + "; n[j] = " + n[j]);//+ "; Math.log(g(j, times[j], time)) = " + Math.log(g(j, times[j], time)));
                if (Double.isInfinite(logP))
                    return logP;

            }
            if (rho[j] > 0 && N[j] > 0) {
                temp = N[j] * Math.log(rho[j]);    // term for contemporaneous sampling
                logP += temp;
                if (printTempResults)
                    System.out.println("3rd factor (Nj loop) = " + temp + "; interval = " + j + "; N[j] = " + N[j]);
                if (Double.isInfinite(logP))
                    return logP;

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
}