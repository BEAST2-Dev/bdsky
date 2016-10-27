package beast.evolution.speciation;

import beast.core.Citation;
import beast.core.Description;
import beast.evolution.tree.TreeInterface;

/**
 * @author Chi Zhang
 */

@Description("Extension of the birth-death skyline model to account for diversified sampling of extant taxa.")
@Citation("Zhang C., Stadler T., Klopfstein S., Heath T.A., Ronquist F. 2016. Total-evidence dating under the fossilized " +
        "birth-death process. Syst. Biol. 65:228–249.  See also: Stadler T., Smrckova J. 2016. Estimating shifts in " +
        "diversification rates based on higher-level phylogenies. Biol. Lett. 12:20160273.")
public class BirthDeathSkylineDiversifiedSampling extends BirthDeathSkylineModel {


    @Override
    protected Double updateRatesAndTimes(TreeInterface tree) {

        // Insert cut time to sampling rate change times
        // samplingRateChangeTimes = ...

        Double d = super.updateRatesAndTimes(tree);

        // Increase dimension of sampling rate, add value 0 to the last entry
        // Double[] samplingRates = samplingRate.get().getValues();

        return d;

    }





}
