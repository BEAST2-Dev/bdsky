package beast.evolution.speciation;


import beast.evolution.tree.Tree;
import beast.evolution.alignment.Taxon;
import beast.core.parameter.*;
import beast.core.Input;
import beast.core.Description;
import beast.core.util.ParameterConstrainer;
//import beast.util.ParameterConstrainer;

import java.util.*;

/**
 * @author Denise Kuehnert
 *
 * maths: Tanja Stadler
 */

@Description("Adaptation of Tanja Stadler's BirthDeathSamplingModel, to allow for birth and death rates to change at times t_i")
public class BirthDeathSkylineModel extends SpeciesTreeDistribution{

    // assume equidistant intervals if intervaltimes are not specified
    public Input<RealParameter> intervalTimes =
            new Input<RealParameter>("intervalTimes", "The times t_i specifying when rate changes can occur", (RealParameter) null);
    public Input<RealParameter> orig_root =
            new Input<RealParameter>("orig_root", "The origin of infection x0", Input.Validate.REQUIRED);
    public Input<Integer> intervalNumber =
            new Input<Integer>("intervalNumber", "The number of intervals in which rates can change", Input.Validate.REQUIRED);

    public Input<RealParameter> birthRate =
            new Input<RealParameter>("birthRate", "BirthRate = BirthRateVector * birthRateScalar, birthrate can change over time");
    public Input<RealParameter> deathRate =
            new Input<RealParameter>("deathRate", "The deathRate vector with birthRates between times");
    public Input<RealParameter> samplingRate =
            new Input<RealParameter>("samplingRate", "The sampling rate per individual");      // psi

    public Input<RealParameter> m_rho =
            new Input<RealParameter>("rho", "The proportion of samples taken at times m_i (default 0.)");
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
    protected int   [] N;   // number of leaves sampled at each time t_i
    protected int m;


    Double[] birth;
    Double[] death;
    Double[] psi;
    Double[] rho;
    boolean isRhoTip[];
    Boolean birthChanges;
    Boolean deathChanges;
    Boolean samplingChanges;
    Boolean rhoChanges;
    Boolean contempData;
    Double[] times;
    Boolean timesFromXML;
    Boolean transform;
    Boolean m_forceRateChange;

    Boolean printTempResults;


/************************************************************************************************/
    /*              "constructor"                                                                   */
/************************************************************************************************/

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        m = intervalNumber.get();
        m_forceRateChange = forceRateChange.get();

        contempData =  contemp.get();
        rhoChanges = false;
        if (m_rho.get()!=null){

            if (m_rho.get().getDimension() == 1 ){
                if ( !contempData && ((samplingProportion.get()!=null && samplingProportion.get().getDimension()==1 && samplingProportion.get().getValue()==0.) ||
                      (samplingRate.get()!=null && samplingRate.get().getDimension()==1 && samplingRate.get().getValue()==0.) )) {
                    contempData = true;
                    System.out.println("Parameters were chosen for contemporaneously sampled data. Setting contemp=true.");
                }
            }

            if (contempData){
                if (m_rho.get().getDimension()!=1) throw new RuntimeException("when contemp=true, rho must have dimension 1");

                else if (m>1){
                    rho = new Double[m];
                    Arrays.fill(rho, 0.);
                    rho[m-1] = m_rho.get().getValue();
                    rhoChanges = true;
                }
                else rho = m_rho.get().getValues();
            }
            else
                rho = m_rho.get().getValues();
                if (rho.length != 1 && rho.length != m)
                    throw new RuntimeException("Length of rho parameter ("+ rho.length + ") should be one or equal to intervalNumber (" + m + ")");
                rhoChanges = (rho.length == m && m > 1)  ;
        }
        else rho = new Double[]{0.};
        isRhoTip = new boolean[m_tree.get().getLeafNodeCount()];

        
        if (birthRate.get() != null && deathRate.get() != null && samplingRate.get() != null){

            transform = false;
            death = deathRate.get().getValues();
            psi = samplingRate.get().getValues();
            birth = birthRate.get().getValues();

        }
        else if (R0.get() != null && becomeUninfectiousRate.get() != null && samplingProportion.get() != null){

            transform = true;
            transformParameters(S0_input.get()==null?1:S0_input.get().getValue());

        }

        else{
            throw new RuntimeException("Either specify birthRate, deathRate and samplingRate OR specify R0, becomeUninfectiousRate and samplingProportion!");
        }


        if (transform){
            if (R0.get().getDimension() != 1 && R0.get().getDimension() != m)
                throw new RuntimeException("Length of R0 parameter should be one or equal to intervalNumber (" + m + ")");


            if (becomeUninfectiousRate.get().getDimension() != 1 && becomeUninfectiousRate.get().getDimension() != m)
                throw new RuntimeException("Length of becomeUninfectiousRate parameter should be one or equal to intervalNumber (" + m + ")");

            if (samplingProportion.get().getDimension() != 1 && samplingProportion.get().getDimension() != m)
                throw new RuntimeException("Length of samplingProportion parameter should be one or equal to intervalNumber (" + m + ")");

            birthChanges = (R0.get().getDimension() == m && m > 1)  ;
            samplingChanges = (samplingProportion.get().getDimension() == m && m > 1) ;
            deathChanges = (becomeUninfectiousRate.get().getDimension() == m && m > 1) ;

        }
        else {

            if (birth.length != 1 && birth.length != m)
                throw new RuntimeException("Length of birthRate parameter ("+ birth.length + ") should be one or equal to intervalNumber (" + m + ")");


            if (death.length != 1 && death.length != m)
                throw new RuntimeException("Length of mu parameter should be one or equal to intervalNumber (" + m + ")");


            if (psi.length != 1 && psi.length != m)
                throw new RuntimeException("Length of birthRate parameter ("+ psi.length + ") should be one or equal to intervalNumber (" + m + ")");

            birthChanges = (birthRate.get().getDimension() == m && m > 1)  ;
            deathChanges = (death.length == m && m > 1) ;
            samplingChanges = (psi.length == m && m > 1) ;
        }


        if (intervalTimes.get() != null){
            if (intervalTimes.get().getDimension() != m)
                throw new RuntimeException("Length of intervalTimes parameter should equal to intervalNumber (" + m + ")");

            if (intervalTimes.get().getValue() != 0)
                throw new RuntimeException("First entry in intervalTimes must be 0!");
            timesFromXML = true;
            if (m_forceRateChange && intervalTimes.get().getArrayValue(m-1) >  m_tree.get().getRoot().getHeight())
                throw new RuntimeException("IntervalTimes have to be between [0, tree + origroot] ([0," + (orig_root.get().getValue() + m_tree.get().getRoot().getHeight()) + "])!");

        }
        else {
            times  = new Double[m];
            for (int i = 0; i < m; i++) {
                timesFromXML = false;
            }
        }

        printTempResults = false;
    }

    private double computeN(Tree tree){

        isRhoTip = new boolean[tree.getLeafNodeCount()];

        N = new int[m];

        int tipCount = tree.getLeafNodeCount();

        double[] dates= new double[tipCount];

        for (int i = 0; i<tipCount; i++){
            dates[i] = tree.getNode(i).getDate();
        }

        for (int k = 0; k<m; k++){


            for (int i = 0; i<tipCount; i++){

                if (Math.abs((times[m-1] - times[k]) - dates[i]) < 1e-10){
                    if (rho[rhoChanges?k:0] ==0)
                        return Double.NEGATIVE_INFINITY;
                    N[k] += 1;
                    isRhoTip[i] = true;
                }
            }
        }
        return 0.;
    }

    private Double getTimes(Tree tree){

        times  = new Double[m];
        double orig = orig_root.get().getValue();
        t_root = tree.getRoot().getHeight();

        if (timesFromXML){
            Double[] tempTimes = intervalTimes.get().getValues();

            // if forceRateChange: force rate change to be within tree range
            if (m_forceRateChange && tempTimes[m-1] > t_root )
                return Double.NEGATIVE_INFINITY;


            // make sure intervalTimes are increasing
            for (int i = 1; i < tempTimes.length; i++){
                if (tempTimes[i] <= tempTimes[i-1])
                    return Double.NEGATIVE_INFINITY;
                times[i-1]=tempTimes[i];
            }

        } else {
            for (int i = 1; i < m; i++) {
                times[i-1] = i * (t_root+orig) / m;
            }
        }
        times[m-1] = t_root + orig;
        

        return 0.;
    }

/************************************************************************************************/
    /*                   calculations for likelihood                                                */
/************************************************************************************************/

    /*    calculate and store Ai, Bi and p0        */
    public Double preCalculation(Tree tree){


        updateRates(tree);
        
        if (m_rho.get()!=null) {
            if (contempData){
                Arrays.fill(rho, 0.);
                rho[m-1] = m_rho.get().getValue();
            }
            else
                rho = m_rho.get().getValues();
        }
        else rho = new Double[]{0.};


        if (getTimes(tree) < 0)
            return Double.NEGATIVE_INFINITY;

        if (psiConstrainer.get()!=null)
            psi = psiConstrainer.get().get(psi, m) ;

        if (m_rho.get()!=null)
            if (computeN(tree) < 0)
                return Double.NEGATIVE_INFINITY;  // todo: check if it's enough to do this once at the beginning (only if intervaltimes and sample dates don't change)

        Ai = new double[m];
        Bi = new double[m];
        p0 = new double[m];

        for (int i = 0; i < m; i++){
            Ai[i] = Ai(birth[birthChanges? i : 0], death[deathChanges? i : 0], psi[samplingChanges? i : 0]);
            if (printTempResults) System.out.println("Ai[" + i + "] = " + Ai[i] + " " +  Math.log(Ai[i]));
        }

        Bi[m-1] = Bi(birth[birthChanges?m-1:0], death[deathChanges?m-1:0], psi[samplingChanges?m-1:0],rho[rhoChanges?m-1:0], Ai[m-1], 1.);  //  (p0[m-1] = 1)
        if (printTempResults) System.out.println("Bi[m-1] = " + Bi[m-1] + " " + Math.log(Bi[m-1]));
        for (int i = m-2; i>=0; i--){
            p0[i+1] = p0(birth[birthChanges? (i+1) : 0], death[deathChanges? (i+1) : 0], psi[samplingChanges? (i+1) : 0], Ai[i+1], Bi[i+1], times[i+1], times[i] );
            if (Math.abs(p0[i+1]-1)<1e-10) {
                return Double.NEGATIVE_INFINITY;
            }
            if (printTempResults) System.out.println("p0[" + (i+1) + "] = " + p0[i+1] );
            Bi[i] = Bi(birth[birthChanges? i : 0], death[deathChanges? i : 0], psi[samplingChanges? i : 0],rho[rhoChanges? i : 0], Ai[i], p0[i+1]);
            if (printTempResults) System.out.println("Bi[" + i + "] = " + Bi[i] + " " +  Math.log(Bi[i]));
        }

        if (printTempResults) {
            System.out.println("g(index, x0, times[index]):" +g(0, times[0],0));
            System.out.println("g 1 :" +g(index(1), times[index(1)],1.));
            System.out.println("g 2 :" +g(index(2), times[index(2)],2));
            System.out.println("g 4 :" +g(index(4), times[index(4)],4));
        }

        return 0.;
    }

    public double Ai(double b, double g, double psi){

        return Math.sqrt((b - g - psi)*(b - g - psi) + 4*b*psi);
    }

    public double Bi(double b, double g, double psi, double r, double A, double p0){

        return ((1-2*p0*(1-r))*b + g + psi)/A;
    }

    public double p0(int index, double t, double ti){

        return p0(birth[birthChanges? index : 0], death[deathChanges?index:0], psi[samplingChanges? index : 0], Ai[index], Bi[index], t, ti);
    }

    public double p0(double b, double g, double psi, double A, double B, double ti, double t){

        if (printTempResults) System.out.println("in p0: b = " + b + "; g = " + g  + "; psi = " +psi+ "; A = " + A + " ; B = " + B + "; ti = " + ti + "; t = " + t);
//        return ((b + g + psi - A *((Math.exp(A*(ti - t))*(1+B)-(1-B)))/(Math.exp(A*(ti - t))*(1+B)+(1-B)) ) / (2*b));
        // formula from manuscript slightly rearranged for numerical stability
        return ((b + g + psi - A *((1+B)-(1-B)*(Math.exp(A*(t - ti))))/((1+B)+Math.exp(A*(t - ti))*(1-B)) ) / (2*b));

    }

    public double g(int index, double ti, double t){

//        return (Math.exp(Ai[index]*(ti - t))) / (0.25*Math.pow((Math.exp(Ai[index]*(ti - t))*(1+Bi[index])+(1-Bi[index])),2));
        // formula from manuscript slightly rearranged for numerical stability
        return (4*Math.exp(Ai[index]*(t - ti))) /(Math.exp(Ai[index]*(t - ti))*(1-Bi[index])+(1+Bi[index])) / (Math.exp(Ai[index]*(t - ti))*(1-Bi[index])+(1+Bi[index]));
    }

    public int index(double t) {

        int epoch = Arrays.binarySearch(times, t);

        if (epoch < 0) {
            epoch = -epoch - 1;
        }

        return Math.min(epoch, m-1); //Math.max((epoch - 1), 0);
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

    public void transformParameters(int S0){


        birth = new Double[m];
        death = new Double[m];
        psi = new Double[m];
        
        Double[] p = samplingProportion.get().getValues();
        Double[] b = becomeUninfectiousRate.get().getValues();
        Double[] R = R0.get().getValues();

        for (int i = 0; i < m; i++){
            birth[i] = R[R.length > 1 ? i : 0] * b[b.length > 1 ? i : 0] / S0;
            psi[i] = p[p.length > 1 ? i : 0]  * b[b.length > 1 ? i : 0] ;
            death[i] = b[b.length > 1 ? i : 0] - psi[i];
        }
    }


    public double updateRates(Tree tree){

        if (transform) transformParameters(S0_input.get()==null?1:S0_input.get().getValue());
        else {
            death = deathRate.get().getValues();
            birth = birthRate.get().getValues();
            psi = samplingRate.get().getValues();
        }

        return 0.;
    }

    @Override
    public double calculateTreeLogLikelihood(Tree tree){

        m = intervalNumber.get();

        // number of lineages at each time ti
        int[] n = new int[m];
        int nTips = tree.getLeafNodeCount();
        if (preCalculation(tree) < 0)
            return Double.NEGATIVE_INFINITY;

        double x0 = 0;
        int index = 0;

        double temp;

        // the first factor for origin
        if (!conditionOnSurvival.get())
            temp =  Math.log(g(index, times[index], x0)) ;  // NOT conditioned on at least one sampled individual
        else{
            temp = p0(index, times[index], x0);
            if (temp == 1)
                return Double.NEGATIVE_INFINITY;
            temp =  Math.log(g(index, times[index], x0)/(1 - temp));   // DEFAULT: conditioned on at least one sampled individual
        }

        logP = temp;
        if (Double.isInfinite(logP))
            return logP;

        if (printTempResults) System.out.println("first factor for origin = " + temp);

        // first product term in f[T]
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {

            double x = times[m-1] - tree.getNode(nTips+i).getHeight();
            index = index(x);

            temp = Math.log(birth[birthChanges? index:0] * g(index, times[index], x));
            logP += temp;
            if (printTempResults) System.out.println("1st pwd" +
                    " = " +temp+ "; interval = " + i);
            if (Double.isInfinite(logP))
                return logP;

        }

        // middle product term in f[T]
        for (int i = 0; i < nTips; i++) {

            if (!isRhoTip[i] || m_rho.get()==null ) {
                double y = times[m-1] - tree.getNode(i).getHeight();
                index = index(y);

                temp = Math.log(psi[samplingChanges? index : 0]) - Math.log(g(index, times[index], y) ) ;
                logP += temp;
                if (printTempResults) System.out.println("2nd PI = " + temp);
                if (Double.isInfinite(logP))
                    return logP;

            }
        }

        // last product term in f[T], factorizing from 1 to m
        double time;
        for (int j = 0; j < m; j++){
            time = j<1?0:times[j-1];
            n[j] = (j==(0)? 0 : lineageCountAtTime(times[m-1]-time, tree));

            if (n[j] > 0) {
                temp =  n[j] * Math.log(g(j, times[j], time));
                logP += temp;
                if (printTempResults) System.out.println("3rd factor (nj loop) = " +temp + "; interval = " + j+ "; n[j] = " + n[j] );//+ "; Math.log(g(j, times[j], time)) = " + Math.log(g(j, times[j], time)));
                if (Double.isInfinite(logP))
                    return logP;

            }
            if (rho[rhoChanges? j : 0] > 0  && N[j] > 0 ){
                temp = N[j] * Math.log(rho[rhoChanges? j : 0]);    // term for contemporaneous sampling
                logP += temp;
                if (printTempResults) System.out.println("3rd factor (Nj loop) = " +temp + "; interval = " + j+ "; N[j] = " + N[j]);
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