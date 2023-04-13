package bdsky;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

@Description("Operator necessary for inferring negative diversification rates.")
public class JointDiversificationTurnoverOperator extends Operator {

    public Input<RealParameter> divRateInput = new Input<>("netDiversification",
            "Net diversification rate parameter.", Input.Validate.REQUIRED);

    public Input<RealParameter> turnOverInput = new Input<>("turnOver",
            "Turnover parameter.", Input.Validate.REQUIRED);

    RealParameter divRate, turnOver;

    @Override
    public void initAndValidate() {
        divRate = divRateInput.get();
        turnOver = turnOverInput.get();

        if (divRate.getDimension() != turnOver.getDimension())
            throw new IllegalArgumentException("JointDiversificationTurnoverOperator " +
                    "requires that the dimension of diversification rate and " +
                    "turnover are identical.");
    }

    @Override
    public double proposal() {

        int idx = Randomizer.nextInt(divRate.getDimension());

        double newDiv = -divRate.getValue(idx);
        double newTurnOver = 1.0/turnOver.getValue(idx);

        if (newDiv < divRate.getLower() || newDiv > divRate.getUpper()
                || newTurnOver < turnOver.getLower() || newTurnOver > turnOver.getUpper())
            return Double.NEGATIVE_INFINITY;

        divRate.setValue(idx, newDiv);
        turnOver.setValue(idx, newTurnOver);

        return 2*Math.log(newTurnOver);
    }
}
