package beast.evolution.speciation;


import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.util.*;

/**
 * @author Denise Kuehnert
 * @author Alexei Drummond
 *         <p/>
 *         maths: Tanja Stadler
 */

@Description("Adaptation of Tanja Stadler's BirthDeathSamplingModel, " +
        "to allow for birth and death rates to change at times t_i")
@Citation("Stadler, T., Kuehnert, D., Bonhoeffer, S., and Drummond, A. J. (2013). “Birth-death skyline" +
        "plot reveals temporal changes of epidemic spread in HIV and hepatitis C virus (HCV).” Proc" +
        "Natl Acad Sci U S A, 110(1): 228–33.")
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

    public Input<Integer> intervalNumber =
            new Input<Integer>("intervalNumber", "The number of intervals in which rates can change", 1);

    public Input<RealParameter> intervalTimes =
            new Input<RealParameter>("intervalTimes", "The time t_i for all parameters if they are the same", (RealParameter) null);

    public Input<Boolean> birthRateChangeTimesRelativeInput =
            new Input<Boolean>("birthRateTimesRelative", "True if birth rate change times specified relative to tree height? Default false", false);

    public Input<Boolean> deathRateChangeTimesRelativeInput =
            new Input<Boolean>("deathRateTimesRelative", "True if death rate change times specified relative to tree height? Default false", false);

    public Input<Boolean> samplingRateChangeTimesRelativeInput =
            new Input<Boolean>("samplingRateTimesRelative", "True if sampling rate times specified relative to tree height? Default false", false);

    BooleanParameter bp = new BooleanParameter();

    {
        try {
            bp.initByName("value", "false false false false");
        } catch (Exception e) { //ignore
        }
    }

    public Input<BooleanParameter> reverseTimeArrays =
            new Input<BooleanParameter>("reverseTimeArrays", "True if the time arrays are given in backwards time (from the present back to root). Order: 1) birth 2) death 3) sampling 4) rho. Default false." +
                    "Careful, rate array must still be given in FORWARD time (root to tips). If rhosamplingTimes given, they should be backwards and this should be true.",
                    bp);

    // the times for rho sampling
    public Input<RealParameter> rhoSamplingTimes =
            new Input<RealParameter>("rhoSamplingTimes", "The times t_i specifying when rho-sampling occurs", (RealParameter) null);


    public Input<RealParameter> origin =
            new Input<RealParameter>("origin", "The time from origin to last sample (must be larger than tree height)", Input.Validate.REQUIRED);

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


    double t_root;
    protected double[] p0;
    protected double[] Ai;
    protected double[] Bi;
    protected int[] N;   // number of leaves sampled at each time t_i

    // these four arrays are totalIntervals in length
    protected Double[] birth;
    Double[] death;
    Double[] psi;
    Double[] rho;

    // true if the node of the given index occurs at the time of a rho-sampling event
    boolean[] isRhoTip;

    /**
     * The number of change points in the birth rate
     */
    protected int birthChanges;

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
    Boolean constantRho;

    /**
     * Total interval count
     */
    protected int totalIntervals;

    protected List<Double> birthRateChangeTimes = new ArrayList<Double>();
    protected List<Double> deathRateChangeTimes = new ArrayList<Double>();
    protected List<Double> samplingRateChangeTimes = new ArrayList<Double>();

    Boolean contempData;
    //List<Interval> intervals = new ArrayList<Interval>();
    SortedSet<Double> timesSet = new TreeSet<Double>();

    protected Double[] times = new Double[]{0.};

    protected Boolean transform;
    Boolean m_forceRateChange;

    Boolean birthRateTimesRelative = false;
    Boolean deathRateTimesRelative = false;
    Boolean samplingRateTimesRelative = false;

    public Boolean printTempResults;

/************************************************************************************************/
    /*              "constructor"                                                                   */

    /**
     * ********************************************************************************************
     */

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        birth = null;
        death = null;
        psi = null;
        rho = null;
        birthRateChangeTimes.clear();
        deathRateChangeTimes.clear();
        samplingRateChangeTimes.clear();
        totalIntervals = 0;

        m_forceRateChange = forceRateChange.get();
        birthRateTimesRelative = birthRateChangeTimesRelativeInput.get();
        deathRateTimesRelative = deathRateChangeTimesRelativeInput.get();
        samplingRateTimesRelative = samplingRateChangeTimesRelativeInput.get();

        contempData = contemp.get();
        rhoSamplingCount = 0;
        printTempResults = false;


        if (birthRate.get() != null && deathRate.get() != null && samplingRate.get() != null) {

            transform = false;
            death = deathRate.get().getValues();
            psi = samplingRate.get().getValues();
            birth = birthRate.get().getValues();

        } else if (R0.get() != null && becomeUninfectiousRate.get() != null && samplingProportion.get() != null) {

            transform = true;

        } else {
            throw new RuntimeException("Either specify birthRate, deathRate and samplingRate OR specify R0, becomeUninfectiousRate and samplingProportion!");
        }

        if (transform) {

            if (birthChanges < 1) birthChanges = R0.get().getDimension() - 1;
            samplingChanges = samplingProportion.get().getDimension() - 1;
            deathChanges = becomeUninfectiousRate.get().getDimension() - 1;

        } else {

            if (birthChanges < 1) birthChanges = birthRate.get().getDimension() - 1;
            deathChanges = deathRate.get().getDimension() - 1;
            samplingChanges = samplingRate.get().getDimension() - 1;
        }

        collectTimes();

        if (m_rho.get() != null) {

            constantRho = !(m_rho.get().getDimension() > 1);

            if (m_rho.get().getDimension() == 1) {
                if (!contempData && ((samplingProportion.get() != null && samplingProportion.get().getDimension() == 1 && samplingProportion.get().getValue() == 0.) ||
                        (samplingRate.get() != null && samplingRate.get().getDimension() == 1 && samplingRate.get().getValue() == 0.))) {
                    contempData = true;
                    if (printTempResults)
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
                if (rhoSampling != null) {
                    for (int i = 0; i < rhoSampling.getDimension(); i++) {
                        rho[index(reverseTimeArrays.get().getValue(3) ? (times[totalIntervals - 1] - rhoSampling.getValue(rhoSampling.getDimension() - i - 1)) : rhoSampling.getValue(i))]
                                = m_rho.get().getValue(constantRho ? 0 : i);
                    }
                    rhoSamplingCount = rho.length;
                }
            }
        } else {
            rho = new Double[totalIntervals];
            Arrays.fill(rho, 0.);
        }
        isRhoTip = new boolean[treeInput.get().getLeafNodeCount()];

        printTempResults = false;
    }

    /**
     * @return a list of intervals
     */
    public void getChangeTimes(List<Double> changeTimes, RealParameter intervalTimes, int numChanges, boolean relative, boolean reverse) {
        changeTimes.clear();

        if (printTempResults) System.out.println("relative = " + relative);

        double maxTime = origin.get().getValue(); // m_tree.get().getRoot().getHeight() + orig_root.get().getValue();

        if (intervalTimes == null) {
            //equidistant

            double intervalWidth = maxTime / (numChanges + 1);

            double start = 0.0;
            double end;
            for (int i = 1; i <= numChanges; i++) {
                end = (intervalWidth) * i;
                changeTimes.add(end);
                start = end;
            }
            end = maxTime;
            changeTimes.add(end);

        } else {

            if (!reverse && intervalTimes.getValue(0) != 0.0) {
                throw new RuntimeException("First time in interval times parameter should always be zero.");
            }

            if (reverse && intervalTimes.getValue(intervalTimes.getDimension() - 1) != 0.0) {
                throw new RuntimeException("Last time in backward interval times parameter should always be zero");
            }

            if ((!isBDSIR()) && numChanges > 0 && intervalTimes.getDimension() != numChanges + 1) {
                throw new RuntimeException("The time interval parameter should be numChanges + 1 long (" + (numChanges + 1) + ").");
            }

//            if (numChanges>0){
            int dim = intervalTimes.getDimension();
            double start = 0.0;
            double end;
            for (int i = 1; i < dim; i++) {
                end = reverse ? (maxTime - intervalTimes.getValue(dim - i - 1)) : intervalTimes.getValue(i);
                if (relative) end *= maxTime;
                changeTimes.add(end);
                start = end;
            }
            end = maxTime;
            changeTimes.add(end);
        }
//        }
    }

    /*
    * Counts the number of tips at each of the contemporaneous sampling times ("rho" sampling time)
    * @return negative infinity if tips are found at a time when rho is zero, zero otherwise.
    */
    private double computeN(TreeInterface tree) {

        isRhoTip = new boolean[tree.getLeafNodeCount()];

        N = new int[totalIntervals];

        int tipCount = tree.getLeafNodeCount();

        double[] dates = new double[tipCount];

        for (int i = 0; i < tipCount; i++) {
            dates[i] = tree.getNode(i).getHeight();
        }

        for (int k = 0; k < totalIntervals; k++) {


            for (int i = 0; i < tipCount; i++) {

                if (Math.abs((times[totalIntervals - 1] - times[k]) - dates[i]) < 1e-10) {
                    if (rho[k] == 0 && psi[k] == 0) {
                        return Double.NEGATIVE_INFINITY;
                    }
                    if (rho[k] > 0) {
                        N[k] += 1;
                        isRhoTip[i] = true;
                    }
                }
            }
        }
        return 0.;
    }

    /**
     * Collect all the times of parameter value changes and rho-sampling events
     */
    private void collectTimes() {

        timesSet.clear();

        if (isBDSIR() && intervalNumber.get() != null) {
            birthChanges = getSIRdimension() - 1;
//            deathChanges = birthChanges;
//            samplingChanges = birthChanges;

        }

        getChangeTimes(birthRateChangeTimes,
                birthRateChangeTimesInput.get() != null && !isSeasonalBDSIR() ? birthRateChangeTimesInput.get() : intervalTimes.get(),
                birthChanges, birthRateTimesRelative, reverseTimeArrays.get().getValue(0));

        getChangeTimes(deathRateChangeTimes,
                deathRateChangeTimesInput.get() != null ? deathRateChangeTimesInput.get() : intervalTimes.get(),
                deathChanges, deathRateTimesRelative, reverseTimeArrays.get().getValue(1));

        getChangeTimes(samplingRateChangeTimes,
                samplingRateChangeTimesInput.get() != null ? samplingRateChangeTimesInput.get() : intervalTimes.get(),
                samplingChanges, samplingRateTimesRelative, reverseTimeArrays.get().getValue(2));

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
        if (rhoSampling != null) {

            double maxTime = origin.get().getValue();
            int dim = rhoSampling.getDimension();

            for (int i = 0; i < dim; i++) {
                //eventsSet.add(new BDSEvent(BDSEvent.Type.rhoSampling, rhoSampling.getValue(i)));
                timesSet.add(reverseTimeArrays.get().getValue(3) ? (maxTime - rhoSampling.getValue(dim - i - 1)) : rhoSampling.getValue(i));
            }
        }

        if (printTempResults) System.out.println("times = " + timesSet);

        times = timesSet.toArray(new Double[timesSet.size()]);
        totalIntervals = times.length;

        if (printTempResults) System.out.println("total intervals = " + totalIntervals);

    }

    protected Double updateRatesAndTimes(TreeInterface tree) {

        collectTimes();

        t_root = tree.getRoot().getHeight();

        if (m_forceRateChange && timesSet.last() > origin.get().getValue()) {
            return Double.NEGATIVE_INFINITY;
        }

        if (transform)
            transformParameters();
        else {

            Double[] birthRates = birthRate.get().getValues();
            Double[] deathRates = deathRate.get().getValues();
            Double[] samplingRates = samplingRate.get().getValues();

            birth = new Double[totalIntervals];
            death = new Double[totalIntervals];
            psi = new Double[totalIntervals];

            birth[0] = birthRates[0];

            for (int i = 0; i < totalIntervals; i++) {
                if (!isBDSIR()) birth[i] = birthRates[index(times[i], birthRateChangeTimes)];
                death[i] = deathRates[index(times[i], deathRateChangeTimes)];
                psi[i] = samplingRates[index(times[i], samplingRateChangeTimes)];

                if (printTempResults) {
                    if (!isBDSIR()) System.out.println("birth[" + i + "]=" + birth[i]);
                    System.out.println("death[" + i + "]=" + death[i]);
                    System.out.println("psi[" + i + "]=" + psi[i]);
                }
            }
        }

        if (m_rho.get() != null && rhoSamplingTimes.get() != null) {

            Double[] rhos = m_rho.get().getValues();
            RealParameter rhoSampling = rhoSamplingTimes.get();

            for (int i = 0; i < totalIntervals; i++) {

                for (int j = 0; i < rhos.length; i++) {
                    if (times[i].equals(reverseTimeArrays.get().getValue(3) ? (origin.get().getValue() - rhoSampling.getValue(rhoSampling.getDimension() - j - 1)) : rhoSampling.getValue(j)))
                        rho[i] = rhos[j];
                }
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
    public Double preCalculation(TreeInterface tree) {

        if (tree.getRoot().getHeight() >= origin.get().getValue()) {
            return Double.NEGATIVE_INFINITY;
        }

        // updateRatesAndTimes must be called before calls to index() below
        if (updateRatesAndTimes(tree) < 0) {
            return Double.NEGATIVE_INFINITY;
        }

        if (printTempResults) System.out.println("After update rates and times");

        if (m_rho.get() != null) {
            if (contempData) {
                rho = new Double[totalIntervals];
                Arrays.fill(rho, 0.);
                rho[totalIntervals - 1] = m_rho.get().getValue();
            } else {
                rho = new Double[totalIntervals];
                Arrays.fill(rho, 0.);

                RealParameter rhoSampling = rhoSamplingTimes.get();

                for (int i = 0; i < rhoSampling.getDimension(); i++) {
                    rho[index(reverseTimeArrays.get().getValue(3) ? (times[totalIntervals - 1] - rhoSampling.getValue(rhoSampling.getDimension() - i - 1)) : rhoSampling.getValue(i))]
                            = m_rho.get().getValue(constantRho ? 0 : i);

                }
                rhoSamplingCount = rho.length;
            }
        } else {
            rho = new Double[totalIntervals];
            Arrays.fill(rho, 0.0);
        }

        if (m_rho.get() != null)
            if (computeN(tree) < 0)
                return Double.NEGATIVE_INFINITY;

        int intervalCount = times.length;

        //System.out.println("intervalCount=" + intervalCount);

        Ai = new double[intervalCount];
        Bi = new double[intervalCount];
        p0 = new double[intervalCount];

        for (int i = 0; i < intervalCount; i++) {

            Ai[i] = Ai(birth[i], death[i], psi[i]);

            if (printTempResults) System.out.println("Ai[" + i + "] = " + Ai[i] + " " + Math.log(Ai[i]));
        }

        if (printTempResults) {
            System.out.println("birth[m-1]=" + birth[totalIntervals - 1]);
            System.out.println("death[m-1]=" + death[totalIntervals - 1]);
            System.out.println("psi[m-1]=" + psi[totalIntervals - 1]);
            System.out.println("rho[m-1]=" + rho[totalIntervals - 1]);
            System.out.println("Ai[m-1]=" + Ai[totalIntervals - 1]);
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

            p0[i + 1] = p0(birth[i + 1], death[i + 1], psi[i + 1], Ai[i + 1], Bi[i + 1], times[i + 1], times[i]);
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

        if (t >= times[totalIntervals - 1])
            return totalIntervals - 1;

        int epoch = Arrays.binarySearch(times, t);

        if (epoch < 0) {
            epoch = -epoch - 1;
        }

        return epoch; //Math.min(epoch, totalIntervals - 1); //Math.max((epoch - 1), 0);
    }


    /**
     * @param time the time
     * @param tree the tree
     * @return the number of lineages that exist at the given time in the given tree.
     */
    public int lineageCountAtTime(double time, TreeInterface tree) {

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

    protected void transformParameters() {

        Double[] R = R0.get().getValues();
        Double[] b = becomeUninfectiousRate.get().getValues();
        Double[] p = samplingProportion.get().getValues();

        birth = new Double[totalIntervals];
        death = new Double[totalIntervals];
        psi = new Double[totalIntervals];

        if (isBDSIR()) birth[0] = R[0] * b[0]; // the rest will be done in BDSIR class

        for (int i = 0; i < totalIntervals; i++) {
            if (!isBDSIR())
                birth[i] = R[birthChanges > 0 ? index(times[i], birthRateChangeTimes) : 0] * b[deathChanges > 0 ? index(times[i], deathRateChangeTimes) : 0];
            psi[i] = p[samplingChanges > 0 ? index(times[i], samplingRateChangeTimes) : 0] * b[deathChanges > 0 ? index(times[i], deathRateChangeTimes) : 0];
            death[i] = b[deathChanges > 0 ? index(times[i], deathRateChangeTimes) : 0] - psi[i];

        }
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {

        int nTips = tree.getLeafNodeCount();

        if (preCalculation(tree) < 0) {
            return Double.NEGATIVE_INFINITY;
        }

        // number of lineages at each time ti
        int[] n = new int[totalIntervals];

        double x0 = 0;
        int index = 0;

        if (times[index] < 0.)
            index = index(0.);

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
                if (psi[index] == 0 || Double.isInfinite(logP))
                    return logP;

            }
        }

        // last product term in f[T], factorizing from 1 to m
        double time;
        for (int j = 0; j < totalIntervals; j++) {
            time = j < 1 ? 0 : times[j - 1];
            n[j] = ((j == 0) ? 0 : lineageCountAtTime(times[totalIntervals - 1] - time, tree));

            if (n[j] > 0) {
                temp = n[j] * (Math.log(g(j, times[j], time)) + Math.log(1 - rho[j]));
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

    @Override
    public boolean canHandleTipDates() {
        return (m_rho.get() == null);
    }


    public Boolean isBDSIR() {
        return false;
    }

    public Boolean isSeasonalBDSIR() {
        return false;
    }

    public int getSIRdimension() {
        throw new NotImplementedException();
    }
}