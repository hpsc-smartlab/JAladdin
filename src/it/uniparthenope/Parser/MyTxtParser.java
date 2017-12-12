package it.uniparthenope.Parser;

import it.uniparthenope.Debug.MyFileWriter;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;

public class MyTxtParser {
    private String fileName;
    private ArrayList<String> lines;
    private BufferedReader file;

    public MyTxtParser(String fileName){
        this.fileName = fileName;
        try {
            this.file = new BufferedReader(new FileReader(fileName));
            this.lines = new ArrayList<>();
            String line;
            while((line = this.file.readLine()) != null){
                lines.add(line);
            }
            this.file.close();
        } catch (Exception ex){
            ex.printStackTrace();
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("MyTxtParser: "+ex.getMessage());
            debug.CloseFile();
            System.exit(0);
        }
    }

    public String getLine(int i){
        return this.lines.get(i-1);
    }

    public ArrayList<String> getLines() {
        return lines;
    }

    public ArrayList<String> tail(int n){
        int index =(this.lines.size() -1);
        ArrayList<String> out = new ArrayList<>();
        for(int i =0;i<n;i++){
            out.add(this.lines.get(index));
            index --;
        }
        Collections.reverse(out);
        return out;
    }
}
