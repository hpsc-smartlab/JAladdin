package it.uniparthenope.Parser;

import it.uniparthenope.Boxing.RouteInfo;
import it.uniparthenope.Polygon;
import it.uniparthenope.Utility;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.ParseException;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class GeoJsonFormatter {

    public static ArrayList<Polygon> getPolygons(String filename) throws IOException, ParseException{
        ArrayList<Polygon> virtualFances = new ArrayList<>();
        org.json.simple.parser.JSONParser parser = new org.json.simple.parser.JSONParser();
        Object polyFile = parser.parse(new FileReader(filename));
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

    public static void writeGeoJson(RouteInfo gdtRoute, RouteInfo optimalRoute, double[][] gdtRouteCoords, double[][] optimalRouteCoords, String outdir) throws IOException{
        JSONObject geojson = new JSONObject();
        geojson.put("type", "FeatureCollection");
        JSONArray features = new JSONArray();
        String RED = "#ff0000";
        String GREEN = "#00ff00";
        //ADDING ROUTES:
        features.add(getLineString(gdtRoute, gdtRouteCoords,RED));// GDT
        features.add(getLineString(optimalRoute, optimalRouteCoords, GREEN)); //OPTIMAL
        //END
        geojson.put("features", features);
        FileWriter file = new FileWriter(outdir+"/shipRoute.geojson");
        file.write(geojson.toJSONString());
        file.flush();
        file.close();

        geojson = new JSONObject();
        geojson.put("type", "FeatureCollection");
        features = new JSONArray();
        for(int i=0;i<optimalRouteCoords.length;++i){
            features.add(getMarker(optimalRoute, optimalRouteCoords[i][0], optimalRouteCoords[i][1], i));
        }
        geojson.put("features", features);
        file = new FileWriter(outdir+"/waypointInfo.geojson");
        file.write(geojson.toJSONString());
        file.flush();
        file.close();
    }

    private static JSONObject getMarker(RouteInfo optRoute, double lon, double lat, int waypointIdx){
        JSONObject marker = new JSONObject();
        marker.put("type", "Feature");
        JSONObject properties = new JSONObject();
        if((waypointIdx==0) || (waypointIdx==optRoute.getPath().size()-1)){//start point or end point
            properties.put("marker-color","#ffff00");//YELLOW
            if(waypointIdx==0){
                properties.put("marker-symbol","ferry");
                properties.put("nav-time", "0.0 sec");
                properties.put("avg-speed", "0.0 kts");
                properties.put("distance", "0.0 NM");
            } else {
                properties.put("marker-symbol","college");
                properties.put("nav-time", Utility.secs2hms((optRoute.getPartialTimes()[waypointIdx]*3600.0)));
                properties.put("avg-speed", (optRoute.getDr_cum()[waypointIdx]/optRoute.getPartialTimes()[waypointIdx])+" kts");
                properties.put("distance", (optRoute.getDr_cum()[waypointIdx])+" NM");
            }
        } else {
            properties.put("marker-symbol","circle");
            properties.put("marker-color","#7e7e7e");
            properties.put("nav-time", Utility.secs2hms((optRoute.getPartialTimes()[waypointIdx]*3600.0)));
            properties.put("avg-speed", (optRoute.getDr_cum()[waypointIdx]/optRoute.getPartialTimes()[waypointIdx])+" kts");
            properties.put("distance", (optRoute.getDr_cum()[waypointIdx])+" NM");
        }
        properties.put("marker-size","medium");
        marker.put("properties", properties);
        JSONObject geometry = new JSONObject();
        geometry.put("type", "Point");
        JSONArray coordinates = new JSONArray();
        coordinates.add(0, lon);
        coordinates.add(1, lat);
        geometry.put("coordinates",coordinates);
        marker.put("geometry", geometry);
        return marker;
    }


    private static JSONObject getLineString(RouteInfo route, double[][] coords, String lineColor){
        JSONObject feature = new JSONObject();
        feature.put("type","Feature");
        JSONObject properties = new JSONObject();
        properties.put("stroke",lineColor);
        properties.put("stroke-width", 2.5);
        properties.put("stroke-opacity", 1);
        JSONObject property = new JSONObject();
        //ADDING SUMMARY ROUTE INFO TO GEOJSON:
        double time = 3600.0*route.getPartialTimes()[route.getPartialTimes().length-1];
        double avgV = route.getDr_cum()[route.getDr_cum().length-1] / route.getPartialTimes()[route.getPartialTimes().length-1];
        if(route.getType()==0){//GDT Route
            property.put("route-type","geodetic");
            property.put("route-distance", Utility.forDecimalPts(route.getCost()) + " NM");
            if((time == 0) || Double.isNaN(avgV)){
                property.put("error: ", "the route crosses a virtual fence!");
            } else{
                property.put("navigation-time", Utility.secs2hms(time));
                property.put("average-speed", Utility.forDecimalPts(avgV) + " kts");
            }
        } else {
            double navDist = route.getDr_cum()[route.getDr_cum().length-1];
            property.put("route-distance", Utility.forDecimalPts(navDist) + " NM");
            if(route.getType() == 1){//STATIC Route
                property.put("route-type","static-optimal");
            } else { //DYNAMIC Route
                property.put("route-type","dynamic-optimal");
            }
            property.put("navigation-time", Utility.secs2hms(time));
            property.put("average-speed", Utility.forDecimalPts(avgV) + " kts");
        }
        properties.put("summary", property);
        //END

        feature.put("properties", properties);
        JSONObject geometry = new JSONObject();
        geometry.put("type", "LineString");
        //coordinates:

        JSONArray coordinates = new JSONArray();
        for(int ix=0;ix<coords.length;++ix){
            coordinates.add(getPoint(coords[ix][0], coords[ix][1]));

        }


        geometry.put("coordinates", coordinates);
        feature.put("geometry", geometry);
        return feature;
    }

    private static JSONArray getPoint(double lon, double lat){
        JSONArray point = new JSONArray();
        point.add(0, lon);
        point.add(1, lat);
        return point;
    }








}
