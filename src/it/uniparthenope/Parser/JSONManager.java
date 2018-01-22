package it.uniparthenope.Parser;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Reader;

public class JSONManager {
    private String filename;
    private Reader reader;
    private Object jsonObj;
    private JSONObject toWrite;

    public int retrieveInteger(String varName){
        long value;
        JSONObject jsonObject = (JSONObject) jsonObj;
        value = (Long) jsonObject.get(varName);
        return (int) value;
    }

    public long[] retrieveLongArray(String varName){
        long[] values;
        JSONArray vals = (JSONArray) ((JSONObject) jsonObj).get(varName);
        values = new long[vals.size()];
        for(int i=0;i<values.length;++i)
            values[i] = (Long) vals.get(i);
        return values;
    }

    public double[][] retrieveDouble2D(String varName){
        JSONArray rows = (JSONArray) ((JSONObject) jsonObj).get(varName);
        int nRows = rows.size();
        int nCols = ((JSONArray) rows.get(0)).size();
        double[][] matrix = new double[nRows][nCols];
        for(int i=0;i<rows.size();++i){
            JSONArray col = (JSONArray) rows.get(i);
            for(int j=0;j<col.size();++j){
                matrix[i][j] = (double) col.get(j);
            }
        }
        return matrix;
    }

    public void initReading(String filename) throws IOException, org.json.simple.parser.ParseException{
        this.filename = filename;
        reader = new FileReader(filename);
        org.json.simple.parser.JSONParser parser = new org.json.simple.parser.JSONParser();
        jsonObj = parser.parse(reader);
    }


    public void initWriting(String filename){
        toWrite = new JSONObject();
        this.filename = filename;
    }

    public void putInteger(String varName, int value){
        if(toWrite != null)
            toWrite.put(varName, value);
    }

    public void putLongArray(String varname, long[] values){
        if(toWrite != null){
            JSONArray longValues = new JSONArray();
            for(long item : values)
                longValues.add(item);
            toWrite.put(varname, longValues);
        }
    }

    public void putDouble2D(String varname, double[][] matrix){
        if(toWrite != null){
            JSONArray rows = new JSONArray();
            for(double[] row: matrix){
                JSONArray cols = new JSONArray();
                for(double col : row){
                    cols.add(col);
                }
                rows.add(cols);
            }
            toWrite.put(varname, rows);
        }
    }

    public void dispose() throws IOException{
        if(reader != null)
            reader.close();
        if(toWrite != null){
            FileWriter file = new FileWriter(filename);
            file.write(toWrite.toJSONString());
            file.flush();
            file.close();
        }
    }
}
