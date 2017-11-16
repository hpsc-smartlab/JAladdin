package it.uniparthenope;

public class TemporalGrid {
    private double Dt;
    private long minNt;

    public TemporalGrid(){
        this.Dt = 1.00;//h
        this.minNt = 24;
    }

    public double getDt() {
        return Dt;
    }

    public long getMinNt() {
        return minNt;
    }

}
