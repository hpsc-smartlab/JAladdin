package it.uniparthenope.EnvFieldsIn;

import it.uniparthenope.Polygon;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.ParseException;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class GeoJSONAMPInputData implements AMPInputData {
    private String outDir;
    private Object polyFile;

    public GeoJSONAMPInputData(String outDir){
        this.outDir = outDir;
    }

    @Override
    public boolean open(String filename) throws IOException, ParseException {
        org.json.simple.parser.JSONParser parser = new org.json.simple.parser.JSONParser();
        this.polyFile = parser.parse(new FileReader(filename));
        return false;
    }

    @Override
    public ArrayList<Polygon> getPolygons() {
        ArrayList<Polygon> virtualFances = new ArrayList<>();
        JSONArray features = ((JSONArray)((JSONObject) polyFile).get("features"));
        for(Object feature : features){
            Polygon vFence = new Polygon();
            JSONObject jsonFeature = (JSONObject) feature;//POLYGON
            JSONArray coordinates = ((JSONArray)((JSONArray)((JSONObject)jsonFeature.get("geometry")).get("coordinates")).get(0));
            for(int i=0;i<coordinates.size();++i){
                JSONArray elem = (JSONArray) coordinates.get(i);
                vFence.addCorner((double) elem.get(0), (double) elem.get(1));
            }
            virtualFances.add(vFence);
        }
        return virtualFances;
    }

}
