package util;

import beast.core.parameter.RealParameter;
import beast.core.*;

/**
 * @author: Denise Kuehnert
 * Date: Apr 12, 2011
 * Time: 2:41:12 PM
 */

@Description("Multiply each element of a vector with the given scalar")
public class ElementwiseMultiplication extends BEASTObject{

    public Input<RealParameter> vector =
                new Input<RealParameter>("vector", "vector to be multiplied by the scalar", Input.Validate.REQUIRED);

    public Input<Double> scalar =
                    new Input<Double>("scalar", "scalar to multiply with vector", Input.Validate.REQUIRED);


    RealParameter multipliedVector;


    /* Methods */
    public void initAndValidate() throws Exception {
        Double[] mult = new Double[vector.get().getDimension()];
        for (int i =0; i<mult.length; i++){
            mult[i] =   vector.get().getArrayValue(i) * scalar.get();
        }

        multipliedVector = new RealParameter(mult);
    }

    public RealParameter get() throws Exception {

        return multipliedVector;

    }

}
