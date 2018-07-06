package bdsky;


import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.speciation.SpeciesTreeDistribution;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

import java.util.*;

/**
 * @author Alexei Drummond
 * @author Denise Kuehnert
 * @author Alexandra Gavryushkina
 *         <p/>
 *         maths: Tanja Stadler, sampled ancestor extension Alexandra Gavryushkina
 */

@Description("BirthDeathSkylineModel with generalized parameterizations")
@Citation(value = "Stadler, T., Kuehnert, D., Bonhoeffer, S., and Drummond, A. J. (2013):\n Birth-death skyline " +
        "plot reveals temporal changes of\n epidemic spread in HIV and hepatitis C virus (HCV). PNAS 110(1): 228–33.\n" +
        "If sampled ancestors are used then please also cite: Gavryushkina A, Welch D, Stadler T, Drummond AJ (2014) \n" +
        "Bayesian inference of sampled ancestor trees for epidemiology and fossil calibration. \n" +
        "PLoS Comput Biol 10(12): e1003919. doi:10.1371/journal.pcbi.1003919",
        year = 2013, firstAuthorSurname = "Stadler", DOI="10.1073/pnas.1207965110")
public class ParameterizedBirthDeathSkylineModel extends SpeciesTreeDistribution {

    public Input<BDSParameterization> parameterizationInput = new Input<>("parameterization", "The parameterization to use.", Input.Validate.REQUIRED);

    // the times for rho sampling
    public Input<RealParameter> rhoSamplingTimes =
            new Input<RealParameter>("rhoSamplingTimes", "The times t_i specifying when rho-sampling occurs", (RealParameter) null);

    // the rho parameter, one for each rho sampling time
    public Input<RealParameter> rhoInput =
            new Input<RealParameter>("rho", "The proportion of lineages sampled at rho-sampling times (default 0.)");

    public Input<Boolean> originIsRootEdge =
            new Input<>("originIsRootEdge", "The origin is only the length of the root edge", false);

    public Input<Boolean> contemp =
            new Input<Boolean>("contemp", "Only contemporaneous sampling (i.e. all tips are from same sampling time, default false)", false);

    public Input<Boolean> conditionOnSurvival =
            new Input<Boolean>("conditionOnSurvival", "if is true then condition on sampling at least one individual (psi-sampling).", true);
    public Input<Boolean> conditionOnRhoSampling =
            new Input<Boolean> ("conditionOnRhoSampling","if is true then condition on sampling at least one individual at present.", false);

    double t_root;
    protected double[] p0, p0hat;
    protected double[] Ai, Aihat;
    protected double[] Bi, Bihat;
    protected int[] N;   // number of leaves sampled at each time t_i

    // these four arrays are totalIntervals in length
    protected Double[] birth;
    Double[] death;
    Double[] psi;
    Double[] rho;
    Double[] r;

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
    int rhoChanges;

    /**
     * The number of change point in the removal probability
     */
    int rChanges;

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
    protected List<Double> rhoSamplingChangeTimes = new ArrayList<Double>();
    protected List<Double> rChangeTimes = new ArrayList<Double>();

    Boolean contempData;
    //List<Interval> intervals = new ArrayList<Interval>();
    SortedSet<Double> timesSet = new TreeSet<Double>();

    protected Double[] times = new Double[]{0.};

    protected Boolean transform;
    Boolean m_forceRateChange;

    Boolean birthRateTimesRelative = false;
    Boolean deathRateTimesRelative = false;
    Boolean samplingRateTimesRelative = false;
    Boolean rTimesRelative = false;
    //Boolean[] reverseTimeArrays;

    public boolean SAModel;

    enum ConditionOn {NONE, SURVIVAL, RHO_SAMPLING};
    protected ConditionOn conditionOn= ConditionOn.SURVIVAL;

    public Boolean printTempResults;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if (!originIsRootEdge.get() && treeInput.get().getRoot().getHeight() >= origin())
            throw new RuntimeException("Origin parameter ("+ origin() +" ) must be larger than tree height("+treeInput.get().getRoot().getHeight()+" ). Please change initial origin value!");

        // check if this is a sampled ancestor model
        if (parameterizationInput.get().isSampledAncestorModel()) SAModel = true;

        birth = null;
        death = null;
        psi = null;
        rho = null;
        r = null;
        birthRateChangeTimes.clear();
        deathRateChangeTimes.clear();
        samplingRateChangeTimes.clear();
        if (SAModel) rChangeTimes.clear();
        totalIntervals = 0;

        contempData = contemp.get();
        rhoSamplingCount = 0;
        printTempResults = false;

        //if (SAModel) rChanges = removalProbability.get().getDimension() -1;

        if (rhoInput.get()!=null) {
            rho = rhoInput.get().getValues();
            rhoChanges = rhoInput.get().getDimension() - 1;
        }

        collectTimes();

        if (rhoInput.get() != null) {

            constantRho = !(rhoInput.get().getDimension() > 1);

            if (rhoInput.get().getDimension() == 1 && rhoSamplingTimes.get()==null || rhoSamplingTimes.get().getDimension() < 2) {
                // TODO figure this out!
                //if (!contempData && ((samplingProportion.get() != null && samplingProportion.get().getDimension() == 1 && samplingProportion.get().getValue() == 0.) ||
                //        (samplingRate.get() != null && samplingRate.get().getDimension() == 1 && samplingRate.get().getValue() == 0.))) {
                //    contempData = true;
                //    if (printTempResults)
                //        System.out.println("Parameters were chosen for contemporaneously sampled data. Setting contemp=true.");
                //}
            }

            if (contempData) {
                if (rhoInput.get().getDimension() != 1)
                    throw new RuntimeException("when contemp=true, rho must have dimension 1");

                else {
                    rho = new Double[totalIntervals];
                    Arrays.fill(rho, 0.);
                    rho[totalIntervals - 1] = rhoInput.get().getValue();
                    rhoSamplingCount = 1;
                }
            }

        } else {
            rho = new Double[totalIntervals];
            Arrays.fill(rho, 0.);
        }
        isRhoTip = new boolean[treeInput.get().getLeafNodeCount()];

        if (conditionOnSurvival.get()) {
            conditionOn = ConditionOn.SURVIVAL;
            if (conditionOnRhoSampling.get()) {
                throw new RuntimeException("conditionOnSurvival and conditionOnRhoSampling can not be both true at the same time." +
                        "Set one of them to true and another one to false.");
            }
        } else if (conditionOnRhoSampling.get()) {
            if (!rhoSamplingConditionHolds()) {
                throw new RuntimeException("Conditioning on rho-sampling is only available for sampled ancestor analyses where r " +
                        "is set to zero and all except the last rho are zero");
            }
            conditionOn = ConditionOn.RHO_SAMPLING;
        } else {
            conditionOn = ConditionOn.NONE;
        }

        printTempResults = false;
    }

    private double origin() {
        return parameterizationInput.get().origin();
    }


    /**
     * checks if r is zero, all elements of rho except the last one are
     * zero and the last one is not zero
     * @return
     */
    private boolean rhoSamplingConditionHolds() {

        if (SAModel) {
            for (BDSSkylineSegment segment : parameterizationInput.get().canonicalSegments()) {
                if (segment.r() != 0.0) {
                    return false;
                }
            }
        } else return false;

        for (int i=0; i<rho.length-1; i++) {
            if (rho[i] != 0.0) {
                return false;
            }
        }

        return (rho[rho.length-1] != 0.0);
    }

    /**
     * @return a list of intervals
     */
    public void getChangeTimes(List<Double> changeTimes, RealParameter intervalTimes, int numChanges, boolean relative, boolean reverse) {
        changeTimes.clear();

        if (printTempResults) System.out.println("relative = " + relative);

        double maxTime = originIsRootEdge.get()? treeInput.get().getRoot().getHeight() + origin() : origin();

        if (intervalTimes == null) { //equidistant

            double intervalWidth = maxTime / (numChanges + 1);

            double end;
            for (int i = 1; i <= numChanges; i++) {
                end = (intervalWidth) * i;
                changeTimes.add(end);
            }
            end = maxTime;
            changeTimes.add(end);

        } else {

            int dim = intervalTimes.getDimension();

            ArrayList<Double> sortedIntervalTimes = new ArrayList<>();
            for (int i=0; i< dim; i++) {
                sortedIntervalTimes.add(intervalTimes.getValue(i));
            }
            Collections.sort(sortedIntervalTimes);

            if (!reverse && sortedIntervalTimes.get(0) != 0.0) {
                throw new RuntimeException("First time in interval times parameter should always be zero.");
            }

//            if(intervalTimes.getValue(dim-1)==maxTime) changeTimes.add(0.); //rhoSampling

            double end;
            for (int i = (reverse?0:1); i < dim; i++) {
                end = reverse ? (maxTime - sortedIntervalTimes.get(dim - i - 1)) : sortedIntervalTimes.get(i);
                if (relative) end *= maxTime;
                if (end != maxTime) changeTimes.add(end);
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
     * Collect all the times of multiskyline parameterization and the rho-sampling events
     */
    private void collectTimes() {

        timesSet.clear();

        for (BDSSkylineSegment seg : parameterizationInput.get().canonicalSegments()) {
            timesSet.add(seg.start());
        }

        getChangeTimes(rhoSamplingChangeTimes, rhoSamplingTimes.get(), rhoChanges, false, false);

        if (printTempResults) System.out.println("times = " + timesSet);

        times = timesSet.toArray(new Double[timesSet.size()]);
        totalIntervals = times.length;

        if (printTempResults) System.out.println("total intervals = " + totalIntervals);

    }

    protected Double updateRatesAndTimes(TreeInterface tree) {

        collectTimes();

        t_root = tree.getRoot().getHeight();

        int size = ((BDSParameterization)parameterizationInput.get()).size();

        if (birth == null || birth.length != size) birth = new Double[size];
        if (death == null || death.length != size) death = new Double[size];
        if (psi == null || psi.length != size) psi = new Double[size];
        if (r == null || r.length != size) r = new Double[size];

        parameterizationInput.get().populateCanonical(birth, death, psi, r, times);
//            for (int i = 0; i < totalIntervals; i++) {
//                death[i] = deathRates[index(times[i], deathRateChangeTimes)];
//                psi[i] = samplingRates[index(times[i], samplingRateChangeTimes)];
//                if (SAModel) r[i] = removalProbabilities[index(times[i], rChangeTimes)];
//
//                if (printTempResults) {
//                    System.out.println("death[" + i + "]=" + death[i]);
//                    System.out.println("psi[" + i + "]=" + psi[i]);
//                    if (SAModel) System.out.println("r[" + i + "]=" + r[i]);
//                }
//            }

        if (rhoInput.get() != null && (rhoInput.get().getDimension()==1 ||  rhoSamplingTimes.get() != null)) {

            Double[] rhos = rhoInput.get().getValues();
            rho = new Double[totalIntervals];

//            rho[totalIntervals-1]=rhos[rhos.length-1];
            for (int i = 0; i < totalIntervals; i++) {

                rho[i]= //rhoSamplingChangeTimes.contains(times[i]) ? rhos[rhoSamplingChangeTimes.indexOf(times[i])] : 0.;
                        rhoChanges>0?
                                rhoSamplingChangeTimes.contains(times[i]) ? rhos[rhoSamplingChangeTimes.indexOf(times[i])] : 0.
                                : rhos[0];
            }
        }

        return 0.;
    }


    /*    calculate and store Ai, Bi and p0        */
    public Double preCalculation(TreeInterface tree) {

        if (!originIsRootEdge.get() && tree.getRoot().getHeight() >= parameterizationInput.get().origin()) {
            return Double.NEGATIVE_INFINITY;
        }

        // updateRatesAndTimes must be called before calls to index() below
        if (updateRatesAndTimes(tree) < 0) {
            return Double.NEGATIVE_INFINITY;
        }

        if (printTempResults) System.out.println("After update rates and times");

        if (rhoInput.get() != null) {
            if (contempData) {
                rho = new Double[totalIntervals];
                Arrays.fill(rho, 0.);
                rho[totalIntervals-1] = rhoInput.get().getValue();
            }

        } else {
            rho = new Double[totalIntervals];
            Arrays.fill(rho, 0.0);
        }

        if (rhoInput.get() != null)
            if (computeN(tree) < 0)
                return Double.NEGATIVE_INFINITY;

        int intervalCount = times.length;

        Ai = new double[intervalCount];
        Bi = new double[intervalCount];
        p0 = new double[intervalCount];

        if (conditionOn == ConditionOn.RHO_SAMPLING) {
            Aihat = new double[intervalCount];
            Bihat = new double[intervalCount];
            p0hat = new double[intervalCount];
        }

        for (int i = 0; i < intervalCount; i++) {

            Ai[i] = Ai(birth[i], death[i], psi[i]);

            if (conditionOn == ConditionOn.RHO_SAMPLING) {
                Aihat[i] = Ai(birth[i], death[i], 0.0);
            }

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

        if (conditionOn == ConditionOn.RHO_SAMPLING) {
            Bihat[totalIntervals - 1] = Bi(
                    birth[totalIntervals - 1],
                    death[totalIntervals - 1],
                    0.0,
                    rho[totalIntervals - 1],
                    Aihat[totalIntervals - 1], 1.);  //  (p0[m-1] = 1)
        }

        if (printTempResults)
            System.out.println("Bi[m-1] = " + Bi[totalIntervals - 1] + " " + Math.log(Bi[totalIntervals - 1]));
        for (int i = totalIntervals - 2; i >= 0; i--) {

            p0[i + 1] = p0(birth[i + 1], death[i + 1], psi[i + 1], Ai[i + 1], Bi[i + 1], times[i + 1], times[i]);
            if (Math.abs(p0[i + 1] - 1) < 1e-10) {
                return Double.NEGATIVE_INFINITY;
            }
            if (conditionOn == ConditionOn.RHO_SAMPLING) {
                p0hat[i + 1] = p0(birth[i + 1], death[i + 1], 0.0, Aihat[i + 1], Bihat[i + 1], times[i + 1], times[i]);
                if (Math.abs(p0hat[i + 1] - 1) < 1e-10) {
                    return Double.NEGATIVE_INFINITY;
                }
            }
            if (printTempResults) System.out.println("p0[" + (i + 1) + "] = " + p0[i + 1]);

            Bi[i] = Bi(birth[i], death[i], psi[i], rho[i], Ai[i], p0[i + 1]);
            if (conditionOn == ConditionOn.RHO_SAMPLING) {
                Bihat[i] = Bi(birth[i], death[i], 0.0, rho[i], Aihat[i], p0hat[i + 1]);
            }

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

    public double p0hat(int index, double t, double ti) {

        return p0(birth[index], death[index], 0.0, Aihat[index], Bihat[index], t, ti);
    }


    public double g(int index, double ti, double t) {

//        return (Math.exp(Ai[index]*(ti - t))) / (0.25*Math.pow((Math.exp(Ai[index]*(ti - t))*(1+Bi[index])+(1-Bi[index])),2));
        // formula from manuscript slightly rearranged for numerical stability
        return (4 * Math.exp(Ai[index] * (t - ti))) / (Math.exp(Ai[index] * (t - ti)) * (1 - Bi[index]) + (1 + Bi[index])) / (Math.exp(Ai[index] * (t - ti)) * (1 - Bi[index]) + (1 + Bi[index]));
    }

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

        return epoch;
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
            if (tree.getNode(i).getHeight() >= time) count -= 1;
        }
        return count;
    }

    /**
     * @param time the time
     * @param tree the tree
     * @param k count the number of sampled ancestors at the given time
     * @return the number of lineages that exist at the given time in the given tree.
     */
    public int lineageCountAtTime(double time, TreeInterface tree, int[] k) {

        int count = 1;
        k[0]=0;
        int tipCount = tree.getLeafNodeCount();
        for (int i = tipCount; i < tipCount + tree.getInternalNodeCount(); i++) {
            if (tree.getNode(i).getHeight() >= time) count += 1;

        }
        for (int i = 0; i < tipCount; i++) {
            if (tree.getNode(i).getHeight() > time) count -= 1;
            if (Math.abs(tree.getNode(i).getHeight() - time) < 1e-10) {
                count -= 1;
                if (tree.getNode(i).isDirectAncestor()) {
                    count -= 1;
                    k[0]++;
                }

            }
        }
        return count;
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

        double temp=0;

        switch (conditionOn) {
            case NONE:
                temp = Math.log(g(index, times[index], x0));
                break;
            case SURVIVAL:
                temp = p0(index, times[index], x0);
                if (temp == 1)
                    return Double.NEGATIVE_INFINITY;
                temp = Math.log(g(index, times[index], x0) / (1 - temp));
                break;
            case RHO_SAMPLING:
                temp = p0hat(index, times[index], x0);
                if (temp == 1)
                    return Double.NEGATIVE_INFINITY;
                temp = Math.log(g(index, times[index], x0) / (1 - temp));
                break;
            default:
                break;
        }

        logP = temp;
        if (Double.isInfinite(logP))
            return logP;

        if (printTempResults) System.out.println("first factor for origin = " + temp);

        // first product term in f[T]
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {

            double x = times[totalIntervals - 1] - tree.getNode(nTips + i).getHeight();
            index = index(x);
            if (!(tree.getNode(nTips + i)).isFake()) {
                temp = Math.log(birth[index] * g(index, times[index], x));
                logP += temp;
                if (printTempResults) System.out.println("1st pwd" +
                        " = " + temp + "; interval = " + i);
                if (Double.isInfinite(logP))
                    return logP;
            }
        }

        // middle product term in f[T]
        for (int i = 0; i < nTips; i++) {

            if (!isRhoTip[i] || rhoInput.get() == null) {
                double y = times[totalIntervals - 1] - tree.getNode(i).getHeight();
                index = index(y);

                if (!(tree.getNode(i)).isDirectAncestor()) {
                    if (!SAModel) {
                        temp = Math.log(psi[index]) - Math.log(g(index, times[index], y));
                    } else {
                        temp = Math.log(psi[index] * (r[index] + (1 - r[index]) * p0(index, times[index], y))) - Math.log(g(index, times[index], y));
                    }
                    logP += temp;
                    if (printTempResults) System.out.println("2nd PI = " + temp);
                    if (psi[index] == 0 || Double.isInfinite(logP))
                        return logP;
                } else {
                    if (r[index] != 1) {
                        logP += Math.log((1 - r[index])*psi[index]);
                        if (Double.isInfinite(logP)) {
                            return logP;
                        }
                    } else {
                        //throw new Exception("There is a sampled ancestor in the tree while r parameter is 1");
                        System.out.println("There is a sampled ancestor in the tree while r parameter is 1");
                        System.exit(0);
                    }
                }
            }
        }

        // last product term in f[T], factorizing from 1 to m //
        double time;
        for (int j = 0; j < totalIntervals; j++) {
            time = j < 1 ? 0 : times[j - 1];
            int[] k = {0};
            if (!SAModel) {
                n[j] = ((j == 0) ? 0 : lineageCountAtTime(times[totalIntervals - 1] - time, tree));
            } else {
                n[j] = ((j == 0) ? 0 : lineageCountAtTime(times[totalIntervals - 1] - time, tree, k));
            }
            if (n[j] > 0) {
                temp = n[j] * (Math.log(g(j, times[j], time)) + Math.log(1 - rho[j-1]));
                logP += temp;
                if (printTempResults)
                    System.out.println("3rd factor (nj loop) = " + temp + "; interval = " + j + "; n[j] = " + n[j]);//+ "; Math.log(g(j, times[j], time)) = " + Math.log(g(j, times[j], time)));
                if (Double.isInfinite(logP))
                    return logP;

            }

            if (SAModel && j>0 && N != null) { // term for sampled leaves and two-degree nodes at time t_i
                logP += k[0] * (Math.log(g(j, times[j], time)) + Math.log(1-r[j])) + //here g(j,..) corresponds to q_{i+1}, r[j] to r_{i+1},
                        (N[j-1]-k[0])*(Math.log(r[j]+ (1-r[j])*p0(j, times[j], time))); //N[j-1] to N_i, k[0] to K_i,and thus N[j-1]-k[0] to M_i
                if (Double.isInfinite(logP)) {
                    return logP;
                }
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

        if (SAModel) {
            int internalNodeCount = tree.getLeafNodeCount() - ((Tree)tree).getDirectAncestorNodeCount()- 1;
            logP +=  Math.log(2)*internalNodeCount;
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
        return (rhoInput.get() == null);
    }
}
