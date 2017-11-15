package it.uniparthenope;

public class TemporalGrid {
    private double Dt;
    private int minNt;

    public TemporalGrid(){
        this.Dt = 1.00;//h
        this.minNt = 24;
    }

    public double getDt() {
        return Dt;
    }

    public void setDt(double dt) {
        Dt = dt;
    }

    public int getMinNt() {
        return minNt;
    }

    public void setMinNt(int minNt) {
        this.minNt = minNt;
    }
}
