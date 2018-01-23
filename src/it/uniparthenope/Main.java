package it.uniparthenope;


//NOTE:
//visualization (in settings.m file) class not provided in this version.
public class Main {

    public static void main(String[] args) {
        //System.out.println(Runtime.getRuntime().maxMemory());
        //System.out.println(System.getProperty("sun.arch.data.model") );
        //mode=0: normal execution without serialize objects
        //mode=1: normal execution serializing objects
        //mode=2: debugging mode: load data deserializing objects
        JVisirModel v = new JVisirModel(2,2, 1);//Same input provided for testing .m files
        v.Start();
    }
}
