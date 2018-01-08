package it.uniparthenope.Debug;

import java.io.*;

public class ObjectSerializer {

    public static void Serialize(String filename, Object obj) {
        String path = "DebugFiles/SerializedObjects/"+filename+".ser";
        ObjectOutputStream oos = null;
        FileOutputStream fout = null;
        try{
            fout = new FileOutputStream(path,false);
            oos = new ObjectOutputStream(fout);
            oos.writeObject(obj);
            fout.close();
            oos.close();
        } catch (Exception e){
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("ObjectSerializer - Serialization: "+e.getMessage());
            debug.CloseFile();
            e.printStackTrace();
        }
    }

    public static Object Deserialize(String filename){
        String path = "DebugFiles/SerializedObjects/"+filename+".ser";
        Object obj = null;
        try {
            FileInputStream fis = new FileInputStream(path);
            ObjectInputStream ois = new ObjectInputStream(fis);
            obj = ois.readObject();
            fis.close();
            ois.close();
        } catch (Exception e){
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("ObjectSerializer - Deserialization: "+e.getMessage());
            debug.CloseFile();
            e.printStackTrace();
        } finally {
            return obj;
        }
    }
}
