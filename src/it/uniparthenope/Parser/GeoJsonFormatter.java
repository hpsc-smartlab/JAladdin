package it.uniparthenope.Parser;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import java.io.FileWriter;
import java.io.IOException;

public class GeoJsonFormatter {


    public static void writeGeoJson(String filename, double[][] coords) throws IOException{
        JSONObject geojson = new JSONObject();
        geojson.put("type", "FeatureCollection");
        JSONArray features = new JSONArray();
        JSONObject feature = new JSONObject();
        feature.put("type", "Feature");
        feature.put("properties", new JSONObject());
        //geometry
        JSONObject geometry = new JSONObject();
        geometry.put("type", "LineString");
        //coordinates
        JSONArray coordinates = new JSONArray();
        for(int i=0;i<coords.length;++i){
            JSONArray point = new JSONArray();
            point.add(0,coords[i][0]);
            point.add(1,coords[i][1]);
            coordinates.add(point);
        }
        geometry.put("coordinates", coordinates);
        feature.put("geometry", geometry);
        features.add(feature);
        geojson.put("features", features);
        FileWriter file = new FileWriter("Output/"+filename+".geojson");
        file.write(geojson.toJSONString());
        file.flush();
        file.close();
    }

}
