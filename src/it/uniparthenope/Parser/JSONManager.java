package it.uniparthenope.Parser;

import com.sun.org.apache.xpath.internal.operations.Bool;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Reader;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Locale;

public class JSONManager {
    private String filename;
    private Reader reader;
    private Object jsonObj;
    private JSONObject toWrite;
    NumberFormat nf = NumberFormat.getInstance(Locale.US);

    public char retrieveChar(String varName){
        JSONObject jsonObject = (JSONObject) jsonObj;
        return getAschar(jsonObject.get(varName));
        //return (Character) jsonObject.get(varName);
    }

    public int retrieveInteger(String varName){
        long value;
        JSONObject jsonObject = (JSONObject) jsonObj;
        return getAsint(jsonObject.get(varName));
//        value = (Long) jsonObject.get(varName);
//        return (int) value;
    }

    public long retrieveLong(String varName){
        JSONObject jsonObject = (JSONObject) jsonObj;
        return getAslong(jsonObject.get(varName));
        //return (Long) jsonObject.get(varName);
    }

    public double retrieveDouble(String varName)throws ParseException{
        JSONObject jsonObject = (JSONObject) jsonObj;
        return getAsdouble(jsonObject.get(varName));
        //return (Double) jsonObject.get(varName);
    }

    public String retrieveString(String varName){
        JSONObject jsonObject = (JSONObject) jsonObj;
        return (String) jsonObject.get(varName);
    }

    public boolean[] retrieveBooleanArray(String varName){
        boolean[] values;
        JSONArray vals = (JSONArray) ((JSONObject) jsonObj).get(varName);
        values = new boolean[vals.size()];
        for(int i=0;i<values.length;++i)
            values[i] = getAsboolean(vals.get(i));
            //values[i] = (Boolean) vals.get(i);
        return values;
    }

    public int[] retrieveIntArray(String varName){
        int[] values;
        JSONArray vals = (JSONArray) ((JSONObject) jsonObj).get(varName);
        values = new int[vals.size()];
        for(int i=0;i<values.length;++i)
            values[i] = getAsint(vals.get(i));
            //values[i] = ((Long) vals.get(i)).intValue();
        return values;
    }

    public long[] retrieveLongArray(String varName){
        long[] values;
        JSONArray vals = (JSONArray) ((JSONObject) jsonObj).get(varName);
        values = new long[vals.size()];
        for(int i=0;i<values.length;++i)
            values[i] = getAslong(vals.get(i));
            //values[i] = (Long) vals.get(i);
        return values;
    }

    public double[] retrieveDoubleArray(String varName)throws ParseException{
        double[] values;
        JSONArray vals = (JSONArray) ((JSONObject) jsonObj).get(varName);
        values = new double[vals.size()];
        for(int i=0;i<values.length;++i)
            values[i] = getAsdouble(vals.get(i));
            //values[i] = (Double) vals.get(i);
        return values;
    }

    public ArrayList<Double> retrieveDoubleArrayList(String varName)throws ParseException{
        ArrayList<Double> values = new ArrayList<>();
        JSONArray vals = (JSONArray) ((JSONObject) jsonObj).get(varName);
        for(Object val : vals){
            values.add(getAsDouble(val));
            //values.add((Double) val);
        }
        return values;
    }

    public ArrayList<Integer> retrieveIntegerArrayList(String varName){
        ArrayList<Integer> values = new ArrayList<>();
        JSONArray vals = (JSONArray) ((JSONObject) jsonObj).get(varName);
        for(Object val : vals){
            values.add(getAsInteger(val));
            //values.add(((Long) val).intValue());
        }
        return values;
    }

    public ArrayList<Long> retrieveLongArrayList(String varName){
        ArrayList<Long> values = new ArrayList<>();
        JSONArray vals = (JSONArray) ((JSONObject) jsonObj).get(varName);
        for(Object val : vals){
            values.add(getAsLong(val));
            //values.add((Long) val);
        }
        return values;
    }

    public int[][] retrieveInt2D(String varName){
        JSONArray rows = (JSONArray) ((JSONObject) jsonObj).get(varName);
        int nRows = rows.size();
        int nCols = ((JSONArray) rows.get(0)).size();
        int[][] matrix = new int[nRows][nCols];
        for(int i=0;i<rows.size();++i){
            JSONArray col = (JSONArray) rows.get(i);
            for(int j=0;j<col.size();++j){
                matrix[i][j] = getAsint(col.get(j));
                //matrix[i][j] = ((Long) col.get(j)).intValue();
            }
        }
        return matrix;
    }

    public double[][] retrieveDouble2D(String varName)throws ParseException{
        JSONArray rows = (JSONArray) ((JSONObject) jsonObj).get(varName);
        int nRows = rows.size();
        int nCols = ((JSONArray) rows.get(0)).size();
        double[][] matrix = new double[nRows][nCols];
        for(int i=0;i<rows.size();++i){
            JSONArray col = (JSONArray) rows.get(i);
            for(int j=0;j<col.size();++j){
                matrix[i][j] = getAsdouble(col.get(j));
                //matrix[i][j] = (double) col.get(j);
            }
        }
        return matrix;
    }

    public int[][][] retrieveInt3D(String varName){
        JSONArray matrix = (JSONArray) ((JSONObject) jsonObj).get(varName);
        int nZ = matrix.size();
        int nRows = ((JSONArray) matrix.get(0)).size();
        int nCols = ((JSONArray) ((JSONArray) matrix.get(0)).get(0)).size();
        int[][][] mat3d = new int[nZ][nRows][nCols];
        for(int k=0;k<nZ;++k){
            JSONArray k_matrix = (JSONArray) matrix.get(k);
            for(int i=0;i<nRows;++i){
                JSONArray row = (JSONArray) k_matrix.get(i);
                for(int j=0;j<nCols;++j){
                    mat3d[k][i][j] = getAsint(row.get(j));
                    //mat3d[k][i][j] = ((Long) row.get(j)).intValue();
                }
            }
        }
        return mat3d;
    }

    public double[][][] retrieveDouble3D(String varName)throws ParseException{
        JSONArray matrix = (JSONArray) ((JSONObject) jsonObj).get(varName);
        int nZ = matrix.size();
        int nRows = ((JSONArray) matrix.get(0)).size();
        int nCols = ((JSONArray) ((JSONArray) matrix.get(0)).get(0)).size();
        double[][][] mat3d = new double[nZ][nRows][nCols];
        for(int k=0;k<nZ;++k){
            JSONArray k_matrix = (JSONArray) matrix.get(k);
            for(int i=0;i<nRows;++i){
                JSONArray row = (JSONArray) k_matrix.get(i);
                for(int j=0;j<nCols;++j){
                    mat3d[k][i][j] = getAsdouble(row.get(j));
                    //mat3d[k][i][j] = (double) row.get(j);
                }
            }
        }
        return mat3d;
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

    public void putCharacter(String varName, char value){
        if(toWrite != null)
            toWrite.put(varName, getAsString(value));
    }

    public void putInteger(String varName, int value){
        if(toWrite != null)
            toWrite.put(varName, getAsString(value));
    }

    public void putLong(String varName, long value){
        if(toWrite != null)
            toWrite.put(varName, getAsString(value));
    }

    public void putDouble(String varName, double value)throws ParseException{
        if(toWrite != null)
            toWrite.put(varName, getAsString(value));
    }

    public void putString(String varName, String value){
        if(toWrite != null)
            toWrite.put(varName, value);
    }

    public void putLongArray(String varname, long[] values){
        if(toWrite != null){
            JSONArray longValues = new JSONArray();
            for(long item : values)
                longValues.add(getAsString(item));
            toWrite.put(varname, longValues);
        }
    }

    public void putBooleanArray(String varname, boolean[] values){
        if(toWrite != null){
            JSONArray longValues = new JSONArray();
            for(boolean item : values)
                longValues.add(getAsString(item));
            toWrite.put(varname, longValues);
        }
    }

    public void putIntArray(String varname, int[] values){
        if(toWrite != null){
            JSONArray longValues = new JSONArray();
            for(int item : values)
                longValues.add(getAsString(item));
            toWrite.put(varname, longValues);
        }
    }

    public void putDoubleArray(String varname, double[] values)throws ParseException{
        if(toWrite != null){
            JSONArray longValues = new JSONArray();
            for(double item : values)
                longValues.add(getAsString(item));
            toWrite.put(varname, longValues);
        }
    }

    public void putIntegerArrayList(String varName, ArrayList<Integer> list){
        if(toWrite != null){
            JSONArray intValues = new JSONArray();
            for(Integer item : list)
                intValues.add(getAsString(item));
            //intValues.addAll(list);
            toWrite.put(varName, intValues);
        }
    }

    public void putLongArrayList(String varName, ArrayList<Long> list){
        if(toWrite != null){
            JSONArray intValues = new JSONArray();
            for(Long item : list)
                intValues.add(getAsString(item));
            //intValues.addAll(list);
            toWrite.put(varName, intValues);
        }
    }

    public void putDoubleArrayList(String varName, ArrayList<Double> list){
        if(toWrite != null){
            JSONArray intValues = new JSONArray();
            for(Double item : list)
                intValues.add(getAsString(item));
            //intValues.addAll(list);
            toWrite.put(varName, intValues);
        }
    }

    public void putInt2D(String varname, int[][] matrix){
        if(toWrite != null){
            JSONArray rows = new JSONArray();
            for(int[] row: matrix){
                JSONArray cols = new JSONArray();
                for(int col : row){
                    cols.add(getAsString(col));
                }
                rows.add(cols);
            }
            toWrite.put(varname, rows);
        }
    }

    public void putDouble2D(String varname, double[][] matrix)throws ParseException{
        if(toWrite != null){
            JSONArray rows = new JSONArray();
            for(double[] row: matrix){
                JSONArray cols = new JSONArray();
                for(double col : row){
                    cols.add(getAsString(col));
                }
                rows.add(cols);
            }
            toWrite.put(varname, rows);
        }
    }

    public void putInt3D(String varname, int[][][] matrix3d){
        if(toWrite != null){
            JSONArray matrices = new JSONArray();
            for(int[][] matrix : matrix3d){
                JSONArray rows  = new JSONArray();
                for(int[] row : matrix){
                    JSONArray cols = new JSONArray();
                    for(int col: row){
                        cols.add(getAsString(col));
                    }
                    rows.add(cols);
                }
                matrices.add(rows);
            }
            toWrite.put(varname, matrices);
        }
    }

    public void putDouble3D(String varname, double[][][] matrix3d) throws ParseException{
        if(toWrite != null){
            JSONArray matrices = new JSONArray();
            for(double[][] matrix : matrix3d){
                JSONArray rows  = new JSONArray();
                for(double[] row : matrix){
                    JSONArray cols = new JSONArray();
                    for(double col: row){
                        cols.add(getAsString(col));
                    }
                    rows.add(cols);
                }
                matrices.add(rows);
            }
            toWrite.put(varname, matrices);
        }
    }

    private String getAsString(double val){
        return val+"";
    }

    private double getAsdouble(Object val)throws ParseException{
        //return Double.parseDouble(((String) val));
        //return nf.parse((String) val).doubleValue();
        try{
            //return nf.parse((String) val).doubleValue();
            Number n = nf.parse((String) val);
            return Double.valueOf(n.toString());
        } catch (ParseException ex){
            return Double.parseDouble((String) val);
        }
    }

    private String getAsString(Double val){
        return val+"";
    }

    private Double getAsDouble(Object val) throws ParseException{
        //return nf.parse((String) val).doubleValue();
        //return Double.parseDouble(((String) val));
        try{
            Number n = nf.parse((String) val);
            return Double.valueOf(n.toString());
        } catch (ParseException ex){
            return Double.parseDouble((String) val);
        }
    }

    private String getAsString(long val){
        return val+"";
    }

    private long getAslong(Object val){
        return Long.parseLong(((String) val));
    }

    private String getAsString(int val){
        return val+"";
    }

    private int getAsint(Object val){
        return Integer.parseInt((String) val);
    }

    private String getAsString(Integer val){
        return val+"";
    }

    private Integer getAsInteger(Object val){
        return Integer.parseInt((String) val);
    }

    private String getAsString(Long val){
        return val+"";
    }

    private Long getAsLong(Object val){
        return Long.parseLong((String) val);
    }

    private String getAsString(boolean val){
        return val+"";
    }

    private boolean getAsboolean(Object val){
        return Boolean.parseBoolean((String) val);
    }

    private char getAschar(Object val){
        return ((String) val).charAt(0);
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
