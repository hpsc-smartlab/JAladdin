package it.uniparthenope.Debug;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class MyFileWriter {
    private String fileName;
    private BufferedWriter bw;

    public MyFileWriter(String fileName){
        this.fileName = fileName;
        try {
            this.bw = new BufferedWriter(new FileWriter("DebugFiles/"+this.fileName+".txt"));
        } catch (Exception e){
            e.printStackTrace();
        }
    }

    public void WriteLine(String line){
        try {
            this.bw.write(line);
            this.bw.newLine();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void CloseFile(){
        try {
            if(bw!=null){
                bw.close();
            }
        } catch (Exception e){
            e.printStackTrace();
        }
    }


}