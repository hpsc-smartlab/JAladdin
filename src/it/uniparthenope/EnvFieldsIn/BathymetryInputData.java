package it.uniparthenope.EnvFieldsIn;

import java.io.IOException;
import java.util.ArrayList;

public interface BathymetryInputData {
    public boolean open(String filename) throws IOException;
    public ArrayList<Double> getLongitude();
    public ArrayList<Double> getLatitude();
    public double[][] getBathymetry();
    public void dispatch();
}
