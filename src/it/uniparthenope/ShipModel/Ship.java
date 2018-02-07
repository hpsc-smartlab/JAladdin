package it.uniparthenope.ShipModel;

import it.uniparthenope.Const;

public interface Ship {
    public void shipResistance(double wHeight, Const constants);
    public void vesselResponse(Const constants, String outdir);
}
