package it.uniparthenope.Debug;

import com.sun.org.apache.xpath.internal.operations.Bool;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

public class MyFileWriter {
    private String fileName;
    private BufferedWriter bw;
    private static final DateFormat sdf = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");

    public MyFileWriter(String fileName){
        this.fileName = fileName;
        try {
            this.bw = new BufferedWriter(new FileWriter("DebugFiles/"+this.fileName+".txt"));
        } catch (Exception e){
            e.printStackTrace();
        }
    }
    public MyFileWriter(String fileName,String mode, boolean append){
        //this.fileName = fileName;
        try {
            if(mode =="debug"){
                this.bw = new BufferedWriter(new FileWriter("DebugFiles/errors.log",append));
                this.fileName="errors";
            }
            else{
                this.bw = new BufferedWriter(new FileWriter("Output/exec.log",append));
                this.fileName="exec";
            }
            if(!append){
                WriteLine("Current timestamp: "+this.sdf.format(new Date()));
            }
        } catch (Exception e){
            e.printStackTrace();
        }
    }

    private String dateTime(){
        return this.sdf.format(new Date());
    }

    public void WriteLog(String message){
        WriteLine("************************");
        WriteLine(this.dateTime()+":\t"+message);
        WriteLine("************************");
    }

    public void WriteValue(String varName, Double value){
        WriteLine(varName+": "+value);
    }

    public void WriteValue(String varName, Long value){
        WriteLine(varName+": "+value);
    }

    public void WriteValue(String varName, Integer value){
        WriteLine(varName+": "+value);
    }

    public void WriteDoubleArray(String varName, ArrayList<Double> array){
        WriteLine(varName+":");
        for(Double element: array){
            WriteLine(element.toString());
        }
    }

    public void WriteLongArray(String varName, ArrayList<Long> array){
        WriteLine(varName+":");
        for(Long element: array){
            WriteLine(element.toString());
        }
    }

    public void WriteIntegerArray(String varName, ArrayList<Integer> array){
        WriteLine(varName+":");
        for(Integer element: array){
            WriteLine(element.toString());
        }
    }

    public void WriteBooleanArray(String varName, ArrayList<Boolean> array){
        WriteLine(varName+":");
        for(Boolean element: array){
            WriteLine(element.toString());
        }
    }

    public void WriteDoubleArray(String varName, Double[] array){
        WriteLine(varName+":");
        for(int i=0;i<array.length;i++){
            WriteLine(array[i].toString());
        }
    }

    public void WriteLongArray(String varName, Long[] array){
        WriteLine(varName+":");
        for(int i=0;i<array.length;i++){
            WriteLine(array[i].toString());
        }
    }

    public void WriteLongArray(String varName, long[] array){
        WriteLine(varName+":");
        for(int i=0;i<array.length;i++){
            WriteLine(""+array[i]);
        }
    }

    public void WriteIntegerArray(String varName, Integer[] array){
        WriteLine(varName+":");
        for(int i=0;i<array.length;i++){
            WriteLine(array[i].toString());
        }
    }

    public void WriteDoubleMatrix(String varName, Double[][] matrix){
        WriteLine(varName+":");
        for(int i=0;i<matrix.length;i++){
            String row="";
            for(int j=0;j<matrix[0].length;j++){
                row+=matrix[i][j]+" ";
            }
            WriteLine(row);
            WriteLine("");
        }
    }

    public void WriteIntMatrix(String varName, Integer[][] matrix){
        WriteLine(varName+":");
        for(int i=0;i<matrix.length;i++){
            String row="";
            for(int j=0;j<matrix[0].length;j++){
                row+=matrix[i][j]+" ";
            }
            WriteLine(row);
        }
    }

    public void WriteLongMatrix(String varName, Long[][] matrix){
        WriteLine(varName+":");
        for(int i=0;i<matrix.length;i++){
            String row="";
            for(int j=0;j<matrix[0].length;j++){
                row+=matrix[i][j]+" ";
            }
            WriteLine(row);
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