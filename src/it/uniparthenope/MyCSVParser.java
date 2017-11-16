package it.uniparthenope;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class MyCSVParser {
    private InputStream inputStream;
    private ArrayList<String> parsedCSV;

    public MyCSVParser(InputStream inputStream){
        this.inputStream = inputStream;
        this.parsedCSV = new ArrayList<String>();
        this.Parse();
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
