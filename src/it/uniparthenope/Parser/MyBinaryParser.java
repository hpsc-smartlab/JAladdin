package it.uniparthenope.Parser;

import java.io.*;

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

    public long[] readAsUint32(){
        int count = (int) file.length()/4;
        long[] values = new long[count];
        try{
            for(int i=0;i<count;i++) {
                values[i] = (long) this.din.readInt();
            }
        }catch (Exception e){
            e.printStackTrace();
        }
        return values;
    }
}
