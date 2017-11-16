package it.uniparthenope;

import java.io.FileReader;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;

public class MyJSONParser {
    private JSONObject object;

    public MyJSONParser(String path){
        JSONParser parser = new JSONParser();
        try{
            Object obj = parser.parse(new FileReader("inputFiles/"+path));
            object = (JSONObject) obj;
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public String getValueAsString(String key){
        String value = (String) object.get(key);
        if(value == null){
            value = "NOT FOUND";
        }
        return value;
    }

    public long getValueAsLong(String key){
        long value = (long) object.get(key);
        return value;
    }

    public double getValueAsDouble(String key){
        double value = (double) object.get(key);
        return value;
    }
}
