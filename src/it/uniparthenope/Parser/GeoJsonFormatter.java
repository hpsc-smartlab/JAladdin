package it.uniparthenope.Parser;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import java.io.FileWriter;
import java.io.IOException;

public class GeoJsonFormatter {


    public static void writeGeoJson(String filename, double[][] coords) throws IOException{
        JSONObject geojson = new JSONObject();
        //Adding line
        geojson.put("type", "FeatureCollection");
        JSONArray features = new JSONArray();
        JSONObject feature = new JSONObject();
        feature.put("type", "Feature");
        //feature.put("properties", new JSONObject());
        /*TESTING*/
        JSONObject properties = new JSONObject();
        properties.put("stroke","#fb0000");
        properties.put("stroke-width", 3);
        properties.put("stroke-opacity", 1);
        feature.put("properties", properties);
        //END TESTING
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

        //ADDING STARTER MARKER:
        feature = new JSONObject();
        feature.put("type", "Feature");
        properties = new JSONObject();
        properties.put("marker-color","#fffe00");
        properties.put("marker-size","medium");
        properties.put("marker-symbol","circle");
        feature.put("properties", properties);
        geometry = new JSONObject();
        geometry.put("type", "Point");
        coordinates = new JSONArray();
        coordinates.add(0,coords[0][0]);
        coordinates.add(1,coords[0][1]);
        geometry.put("coordinates",coordinates);
        feature.put("geometry", geometry);
        features.add(feature);

        //ADDING ENDING MARKER:
        feature = new JSONObject();
        feature.put("type", "Feature");
        properties = new JSONObject();
        properties.put("marker-color","#2dc702");
        properties.put("marker-size","medium");
        properties.put("marker-symbol","star");
        feature.put("properties", properties);
        geometry = new JSONObject();
        geometry.put("type", "Point");
        coordinates = new JSONArray();
        coordinates.add(0,coords[coords.length-1][0]);
        coordinates.add(1,coords[coords.length-1][1]);
        geometry.put("coordinates",coordinates);
        feature.put("geometry", geometry);
        features.add(feature);
        //END
        geojson.put("features", features);

        FileWriter file = new FileWriter("Output/"+filename+".geojson");
        file.write(geojson.toJSONString());
        file.flush();
        file.close();
    }

}
