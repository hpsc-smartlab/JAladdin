package it.uniparthenope.Parser;

import it.uniparthenope.Utility;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

public class MyBinaryParser {
    private String fileName;
    private File file;
    private DataInputStream din;

    public MyBinaryParser(String fileName){
        this.fileName = fileName;
        try{
            this.file = new File(fileName);
            FileInputStream fin = new FileInputStream(this.file);
            BufferedInputStream bin = new BufferedInputStream(fin);
            this.din = new DataInputStream(bin);

        } catch(Exception e){
            e.printStackTrace();
        }
    }

    private int[] readAsInt32(){
        int count = (int) file.length()/4;
        int[] values = new int[count];
        try{
            for(int i=0;i<count;i++) {
                values[i] = this.din.readInt();
            }
        }catch (Exception e){
            e.printStackTrace();
        }
        return values;
    }


    public long[] readAsUInt32(){
        int[] values = readAsInt32();
        long[] output = new long[values.length];
        for(int i=0;i<values.length;i++){
            output[i]= Utility.getUnsignedInt(Utility.BigToLittleEndian(values[i]));
        }
        return output;
    }

    public int[] readAsUInt16(){
        int[] values = readAsInt32();
        int[] output = new int[values.length];
        for(int i=0;i<values.length;i++){
            output[i]=(int) Utility.getUnsignedInt(Utility.BigToLittleEndian(values[i]));
        }
        return output;
    }
}
