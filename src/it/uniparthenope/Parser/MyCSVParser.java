package it.uniparthenope.Parser;

import java.io.*;
import java.util.ArrayList;

public class MyCSVParser {
    private InputStream inputStream;
    private ArrayList<String> parsedCSV;
    private static final char DEFAULT_SEPARATOR = ',';
    private FileWriter writer;

    public MyCSVParser(InputStream inputStream){
        this.inputStream = inputStream;
        this.parsedCSV = new ArrayList<String>();
        this.Parse();
    }

    public MyCSVParser(String csvFile){
        try {
            this.writer = new FileWriter(csvFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeCSV(double[][] toWrite) throws IOException {
        StringBuilder sb = new StringBuilder();
        int nRows = toWrite.length;
        int nCols = toWrite[0].length;
        for(int i=0; i<nRows; i++){
            boolean first = true;
            for(int j=0;j<nCols;j++){
                if(!first){
                    sb.append(DEFAULT_SEPARATOR);
                }
                sb.append(toWrite[i][j]);
                first=false;
            }
            sb.append("\n");
            this.writer.append(sb.toString());
        }
        this.writer.flush();
        this.writer.close();
    }

    public ArrayList<String> getParsedCSV() {
        return parsedCSV;
    }

    public ArrayList<Double> getDoubleValues(){
        ArrayList<Double> values = new ArrayList<Double>();
        for(String element: this.parsedCSV){
            values.add(Double.parseDouble(element));
        }
        return  values;
    }

    private void Parse(){
        BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
        try {
            String csvLine;
            while ((csvLine = reader.readLine()) != null) {
                String[] row = csvLine.split(",");
                for(int i=0;i<row.length;i++){
                    this.parsedCSV.add(row[i]);
                }
            }
        }
        catch (IOException ex) {
            throw new RuntimeException("Error in reading CSV file: "+ex);
        }
        finally {
            try {
                inputStream.close();
            }
            catch (IOException e) {
                throw new RuntimeException("Error while closing input stream: "+e);
            }
        }
    }

    public long getnElements() {
        return this.parsedCSV.size();
    }

}
