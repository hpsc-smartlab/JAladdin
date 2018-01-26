package it.uniparthenope.Parser;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import java.io.FileWriter;
import java.io.IOException;

public class GeoJsonFormatter {

    public static  void writeGeoJson(double[][]... coordsList) throws IOException{
        JSONObject geojson = new JSONObject();
        geojson.put("type", "FeatureCollection");
        JSONArray features = new JSONArray();
        String[] colors = new String[]{"#ff0000","#00ff00","#ff00ff"};
        int i = 0;
        for(double[][] coords : coordsList){
            if(i==0){// Adding markers
                features.add(getMarker(coords[0][0], coords[0][1], "startPoint"));
                features.add(getMarker(coords[coords.length-1][0], coords[coords.length-1][1], "endPoint"));
            }
            JSONObject feature = new JSONObject();
            feature.put("type","Feature");
            JSONObject properties = new JSONObject();
            properties.put("stroke",colors[i]);
            properties.put("stroke-width", 2.5);
            properties.put("stroke-opacity", 1);
            feature.put("properties", properties);
            JSONObject geometry = new JSONObject();
            geometry.put("type", "LineString");
            //coordinates
            JSONArray coordinates = new JSONArray();
            for(int ix=0;ix<coords.length;++ix){
                JSONArray point = new JSONArray();
                point.add(0,coords[ix][0]);
                point.add(1,coords[ix][1]);
                coordinates.add(point);
            }
            geometry.put("coordinates", coordinates);
            feature.put("geometry", geometry);
            features.add(feature);
            ++i;
        }
        geojson.put("features", features);
        FileWriter file = new FileWriter("Output/shipRoute.geojson");
        file.write(geojson.toJSONString());
        file.flush();
        file.close();
    }


    private static JSONObject getMarker(double lon, double lat, String type){
        JSONObject marker = new JSONObject();
        marker.put("type", "Feature");
        JSONObject properties = new JSONObject();
        properties.put("marker-color","#ffff00");
        properties.put("marker-size","medium");
        if(type == "startPoint")
            properties.put("marker-symbol","ferry");
        else
            properties.put("marker-symbol","college");
        marker.put("properties", properties);
        JSONObject geometry = new JSONObject();
        geometry.put("type", "Point");
        JSONArray coordinates = new JSONArray();
        coordinates.add(0, lon);
        coordinates.add(1, lat);
        geometry.put("coordinates",coordinates);
        marker.put("geometry", geometry);
        return  marker;
    }





}
