package beast.evolution.speciation;

import beast.core.Citation;
import beast.core.Description;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

/**
 * @author Chi Zhang
 */

@Description("Extension of the birth-death skyline model to account for diversified sampling of extant taxa.")
@Citation(value = "Zhang, C., Stadler, T., Klopfstein, S., Heath, T. A., & Ronquist, F. (2015). \n" +
        "Total-evidence dating under the fossilized birth–death process. \n" +
        "Systematic biology, 65(2), 228-249.",
        year = 2015, firstAuthorSurname = "Zhang", DOI="10.1093/sysbio/syv080")
@Citation(value = "Stadler, T., & Smrckova, J. (2016). \n" +
        "Estimating shifts in diversification rates based on higher-level phylogenies. \n" +
        "Biology letters, 12(10), 20160273.",
        year = 2016, firstAuthorSurname = "Stadler", DOI="10.1098/rsbl.2016.0273")
public class BirthDeathSkylineDiversifiedSampling extends BirthDeathSkylineModel {

    private double samplingProp;  // extant sampling proportion
    private int numExtant;        // number of extant taxa
    private double x_cut;

    @Override
    protected Double updateRatesAndTimes(TreeInterface tree) {
        Double value = super.updateRatesAndTimes(tree);

        /* the cut-off time (x_cut) should be smaller than the youngest internal node and fossil, as we assume sampling
           exactly one representative extant species per clade and no fossil sampling between x_cut and the present */
        double minTime = tree.getRoot().getHeight();
        numExtant = 0;
        for(Node node: tree.getNodesAsArray()) {
            if (node.isLeaf()) {
                if (node.getHeight() > 0.0 && node.getHeight() < minTime)
                    minTime = node.getHeight();
                else if (node.getHeight() == 0.0)
                    numExtant++;
            } else {
                if(!node.isFake() && node.getHeight() < minTime)
                    minTime = node.getHeight();
            }
        }
        /* further more, we assume no rate shifting between x_cut and the present (for convenience) */
        double maxTime;
        if (origin.get() == null)
            maxTime = treeInput.get().getRoot().getHeight();
        else
            maxTime = originIsRootEdge.get()? treeInput.get().getRoot().getHeight() + origin.get().getValue() :origin.get().getValue();

        if (totalIntervals > 1 && minTime > maxTime - times[totalIntervals - 2])
            minTime = maxTime - times[totalIntervals - 2];

        // set x_cut to be slightly smaller than minTime
        x_cut = minTime * 0.95;

        // insert x_cut to rate change times
        timesSet.add(maxTime - x_cut);
        times = timesSet.toArray(new Double[timesSet.size()]);

        // increase dimension of birth rate
        Double[] tempb = new Double[totalIntervals + 1];
        System.arraycopy(birth, 0, tempb, 0, totalIntervals);
        tempb[totalIntervals] = birth[totalIntervals - 1];
        birth = tempb;

        // increase dimension of death rate
        Double[] tempd = new Double[totalIntervals + 1];
        System.arraycopy(death, 0, tempd, 0, totalIntervals);
        tempd[totalIntervals] = death[totalIntervals - 1];
        death = tempd;

        // increase dimension of sampling rate
        Double[] tempp = new Double[totalIntervals + 1];
        System.arraycopy(psi, 0, tempp, 0, totalIntervals);
        tempp[totalIntervals] = 0.0;  // add 0 to the last entry
        psi = tempp;

        // set extant sampling proportion here as rho will be overwritten
        samplingProp = rho[totalIntervals - 1];
        // increase dimension of rho rate, add 0 to the last entry
        Double[] tempr = new Double[totalIntervals + 1];
        System.arraycopy(rho, 0, tempr, 0, totalIntervals - 1);
        tempr[totalIntervals - 1] = 0.0;
        tempr[totalIntervals] = 1.0;  // assuming complete sampling
        rho = tempr;

        // increase dimension of r rate
        if (SAModel) {
            Double[] temp = new Double[totalIntervals + 1];
            System.arraycopy(r, 0, temp, 0, totalIntervals);
            temp[totalIntervals] = r[totalIntervals - 1];
            r = temp;
        }

        // increase interval length
        totalIntervals = times.length;

        if (printTempResults) {
            System.out.println("times = " + timesSet);
            System.out.println("total intervals = " + totalIntervals);
            for (int i = 0; i < totalIntervals; i++) {
                if (!isBDSIR()) System.out.println("birth[" + i + "]=" + birth[i]);
                System.out.println("death[" + i + "]=" + death[i]);
                System.out.println("psi[" + i + "]=" + psi[i]);
                if (SAModel) System.out.println("r[" + i + "]=" + r[i]);
                System.out.println("rho[" + i + "]=" + rho[i]);
            }
            System.out.println("samplingProp=" + samplingProp);
        }

        return value;
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {
        logP = super.calculateTreeLogLikelihood(tree);

        // if (logP == Double.POSITIVE_INFINITY) throw new RuntimeException("logP is positive infinity");
        if (Double.isInfinite(logP))
            return Double.NEGATIVE_INFINITY;

        final double lambda = birth[totalIntervals - 1];
        final double mu = death[totalIntervals - 1];
        final double term;
        if ((lambda - mu) * x_cut > 1e-6)
            term = Math.log(lambda * (1 - Math.exp((mu - lambda) * x_cut)))
                 - Math.log(lambda - mu * Math.exp((mu - lambda) * x_cut));
        else if ((mu - lambda) * x_cut > 1e-6)
            term = Math.log(lambda * (1 - Math.exp((lambda - mu) * x_cut)))
                 - Math.log(mu - lambda * Math.exp((lambda - mu) * x_cut));
        else  // for numerical stability
            term = Math.log(lambda / (mu + 1/x_cut));

        // number of extant taxa not sampled
        final int ita_ext = (int)Math.round(numExtant/samplingProp) - numExtant;

        logP += ita_ext * term;

        return logP;
    }
}
