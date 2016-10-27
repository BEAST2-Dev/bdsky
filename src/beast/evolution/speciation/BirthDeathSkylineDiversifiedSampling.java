package beast.evolution.speciation;

import beast.evolution.tree.TreeInterface;

/**
 * Created by dlouis on 27/10/16.
 */
public class BirthDeathSkylineModelPreferentialSampling extends BirthDeathSkylineModel {


    protected Double updateRatesAndTimes(TreeInterface tree) {

        // Insert cut time to sampling rate change times
        // samplingRateChangeTimes = ...

        Double d = super.updateRatesAndTimes(tree);

        // Increase dimension of sampling rate, add value 0 to the last entry
        // Double[] samplingRates = samplingRate.get().getValues();

        return d;

    }





    }
