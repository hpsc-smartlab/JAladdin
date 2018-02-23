package it.uniparthenope.EnvFieldsIn;

import it.uniparthenope.Polygon;
import org.json.simple.parser.ParseException;

import java.io.IOException;
import java.util.ArrayList;

public interface VirtualFencesInputData {
    public boolean open(String filename) throws IOException, ParseException;
    public ArrayList<Polygon> getPolygons();
}
