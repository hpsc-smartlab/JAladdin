package it.uniparthenope.EnvFieldsIn;

import java.io.IOException;

public interface waveInputData {
    public boolean open(String filename) throws IOException;
    public int[] getTimestep();
    public double[] getLongitude();
    public double[] getLatitude();
    public double[][][] getDirection();
    public double[][][] getHeight();
    public double[][][] getPeriod();
    public void dispatch();
}
