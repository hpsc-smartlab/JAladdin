package it.uniparthenope.Parser;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;

public class JSONManager {
    private String filename;
    private Reader reader;
    private Object jsonObj;
    private JSONObject toWrite;

    public char retrieveChar(String varName){
        JSONObject jsonObject = (JSONObject) jsonObj;
        return (Character) jsonObject.get(varName);
    }

    public int retrieveInteger(String varName){
        long value;
        JSONObject jsonObject = (JSONObject) jsonObj;
        value = (Long) jsonObject.get(varName);
        return (int) value;
    }

    public long retrieveLong(String varName){
        JSONObject jsonObject = (JSONObject) jsonObj;
        return (Long) jsonObject.get(varName);
    }

    public double retrieveDouble(String varName){
        JSONObject jsonObject = (JSONObject) jsonObj;
        return (Double) jsonObject.get(varName);
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
            values[i] = (Boolean) vals.get(i);
        return values;
    }

    public int[] retrieveIntArray(String varName){
        int[] values;
        JSONArray vals = (JSONArray) ((JSONObject) jsonObj).get(varName);
        values = new int[vals.size()];
        for(int i=0;i<values.length;++i)
            values[i] = ((Long) vals.get(i)).intValue();
        return values;
    }

    public long[] retrieveLongArray(String varName){
        long[] values;
        JSONArray vals = (JSONArray) ((JSONObject) jsonObj).get(varName);
        values = new long[vals.size()];
        for(int i=0;i<values.length;++i)
            values[i] = (Long) vals.get(i);
        return values;
    }

    public double[] retrieveDoubleArray(String varName){
        double[] values;
        JSONArray vals = (JSONArray) ((JSONObject) jsonObj).get(varName);
        values = new double[vals.size()];
        for(int i=0;i<values.length;++i)
            values[i] = (Double) vals.get(i);
        return values;
    }

    public ArrayList<Double> retrieveDoubleArrayList(String varName){
        ArrayList<Double> values = new ArrayList<>();
        JSONArray vals = (JSONArray) ((JSONObject) jsonObj).get(varName);
        for(Object val : vals){
            values.add((Double) val);
        }
        return values;
    }

    public ArrayList<Integer> retrieveIntegerArrayList(String varName){
        ArrayList<Integer> values = new ArrayList<>();
        JSONArray vals = (JSONArray) ((JSONObject) jsonObj).get(varName);
        for(Object val : vals){
            values.add(((Long) val).intValue());
        }
        return values;
    }

    public ArrayList<Long> retrieveLongArrayList(String varName){
        ArrayList<Long> values = new ArrayList<>();
        JSONArray vals = (JSONArray) ((JSONObject) jsonObj).get(varName);
        for(Object val : vals){
            values.add((Long) val);
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
                matrix[i][j] = ((Long) col.get(j)).intValue();
            }
        }
        return matrix;
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
                    mat3d[k][i][j] = ((Long) row.get(j)).intValue();
                }
            }
        }
        return mat3d;
    }

    public double[][][] retrieveDouble3D(String varName){
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
                    mat3d[k][i][j] = (double) row.get(j);
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
            toWrite.put(varName, value);
    }

    public void putInteger(String varName, int value){
        if(toWrite != null)
            toWrite.put(varName, value);
    }

    public void putLong(String varName, long value){
        if(toWrite != null)
            toWrite.put(varName, value);
    }

    public void putDouble(String varName, double value){
        if(toWrite != null)
            toWrite.put(varName, value);
    }

    public void putString(String varName, String value){
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

    public void putBooleanArray(String varname, boolean[] values){
        if(toWrite != null){
            JSONArray longValues = new JSONArray();
            for(boolean item : values)
                longValues.add(item);
            toWrite.put(varname, longValues);
        }
    }

    public void putIntArray(String varname, int[] values){
        if(toWrite != null){
            JSONArray longValues = new JSONArray();
            for(int item : values)
                longValues.add(item);
            toWrite.put(varname, longValues);
        }
    }

    public void putDoubleArray(String varname, double[] values){
        if(toWrite != null){
            JSONArray longValues = new JSONArray();
            for(double item : values)
                longValues.add(item);
            toWrite.put(varname, longValues);
        }
    }

    public void putIntegerArrayList(String varName, ArrayList<Integer> list){
        if(toWrite != null){
            JSONArray intValues = new JSONArray();
            intValues.addAll(list);
            toWrite.put(varName, intValues);
        }
    }

    public void putLongArrayList(String varName, ArrayList<Long> list){
        if(toWrite != null){
            JSONArray intValues = new JSONArray();
            intValues.addAll(list);
            toWrite.put(varName, intValues);
        }
    }

    public void putDoubleArrayList(String varName, ArrayList<Double> list){
        if(toWrite != null){
            JSONArray intValues = new JSONArray();
            intValues.addAll(list);
            toWrite.put(varName, intValues);
        }
    }

    public void putInt2D(String varname, int[][] matrix){
        if(toWrite != null){
            JSONArray rows = new JSONArray();
            for(int[] row: matrix){
                JSONArray cols = new JSONArray();
                for(int col : row){
                    cols.add(col);
                }
                rows.add(cols);
            }
            toWrite.put(varname, rows);
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

    public void putInt3D(String varname, int[][][] matrix3d){
        if(toWrite != null){
            JSONArray matrices = new JSONArray();
            for(int[][] matrix : matrix3d){
                JSONArray rows  = new JSONArray();
                for(int[] row : matrix){
                    JSONArray cols = new JSONArray();
                    for(int col: row){
                        cols.add(col);
                    }
                    rows.add(cols);
                }
                matrices.add(rows);
            }
            toWrite.put(varname, matrices);
        }
    }

    public void putDouble3D(String varname, double[][][] matrix3d){
        if(toWrite != null){
            JSONArray matrices = new JSONArray();
            for(double[][] matrix : matrix3d){
                JSONArray rows  = new JSONArray();
                for(double[] row : matrix){
                    JSONArray cols = new JSONArray();
                    for(double col: row){
                        cols.add(col);
                    }
                    rows.add(cols);
                }
                matrices.add(rows);
            }
            toWrite.put(varname, matrices);
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
