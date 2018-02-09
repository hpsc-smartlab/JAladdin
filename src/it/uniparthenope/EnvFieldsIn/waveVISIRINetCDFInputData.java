package it.uniparthenope.EnvFieldsIn;

import it.uniparthenope.Debug.MyFileWriter;
import org.jsoup.helper.StringUtil;
import ucar.ma2.ArrayFloat;
import ucar.ma2.ArrayInt;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

import java.io.IOException;

public class waveVISIRINetCDFInputData implements waveInputData {
    private String filename;
    private NetcdfFile dataFile = null;
    private String outDir;

    public waveVISIRINetCDFInputData(String outDir){
        this.outDir = outDir;
    }

    @Override
    public boolean open(String filename) throws IOException {
        this.filename = filename;
        try{
            this.dataFile = NetcdfFile.open(this.filename, null);
            return true;
        } catch (Exception e){
            e.printStackTrace();
            return false;
        }
    }

    @Override
    public double[] getLongitude() {
        //Retrive the variable named "longitude"
        Variable longitude = dataFile.findVariable("lon");
        if(longitude == null){
            System.out.println("Can't find the variable named lon");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
            debug.WriteLog("waveVISIRINetCDFInputData: Can't find the variable named lon");
            debug.CloseFile();
            System.exit(0);
        }
        int[] longitudeShape = longitude.getShape();
        int[] longitudeOrigin = new int[1];
        try{
            double[] parsedLongitude = new double[longitudeShape[0]];
            ArrayFloat.D1 longitudeArray = (ArrayFloat.D1) longitude.read(longitudeOrigin, longitudeShape);
            for(int i=0;i<longitudeShape[0];++i)
                parsedLongitude[i] = longitudeArray.get(i);
            return parsedLongitude;
        }catch (Exception ex){
            errLog(ex);
            return new double[0];
        }
    }

    @Override
    public double[] getLatitude() {
        //Retrive the variable named "latitude"
        Variable latitude = dataFile.findVariable("lat");
        if(latitude == null){
            System.out.println("Can't find the variable named lat");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
            debug.WriteLog("waveVISIRINetCDFInputData: Can't find the variable named lat");
            debug.CloseFile();
            System.exit(0);
        }
        int[] latitudeShape = latitude.getShape();
        int[] latitudeOrigin = new int[1];
        try {
            double[] parsedLatitude = new double[latitudeShape[0]];
            ArrayFloat.D1 latitudeArray = (ArrayFloat.D1) latitude.read(latitudeOrigin, latitudeShape);
            for(int i=0;i<latitudeShape[0];++i){
                parsedLatitude[i] = latitudeArray.get(i);
            }
            return parsedLatitude;
        } catch (Exception ex){
            errLog(ex);
            return new double[0];
        }
    }

    @Override
    public double[][][] getDirection() {
        Variable VDIR = dataFile.findVariable("VDIR");
        if(VDIR==null){
            System.out.println("Can't find the variable named VDIR or dirmn (wave direction)");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
            debug.WriteLog("waveVISIRINetCDFInputData: Can't find the variable named VDIR");
            debug.CloseFile();
            System.exit(0);
        }
        int[] VDIRShape = VDIR.getShape();
        int[] VDIROrigin = new int[3];
        try {
            ArrayFloat.D3 vdir3DMatrix = (ArrayFloat.D3) VDIR.read(VDIROrigin, VDIRShape);
            return ArrayFloatD3ToDouble3D(vdir3DMatrix, VDIRShape);
        } catch (Exception ex){
            errLog(ex);
            return new double[0][][];
        }
    }

    @Override
    public double[][][] getHeight() {
        //Retrive the variable named "VTDH"
        Variable VTDH = dataFile.findVariable("VTDH");
        if(VTDH == null){
            System.out.println("Can't find the variable named VTDH or hs (wave height)");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
            debug.WriteLog("waveVISIRINetCDFInputData: Can't find the variable named VTDH");
            debug.CloseFile();
            System.exit(0);
        }
        int[] VTDHShape = VTDH.getShape();
        int[] VTDHOrigin = new int[3];
        try {
            ArrayFloat.D3 vtdh3DMatrix = (ArrayFloat.D3) VTDH.read(VTDHOrigin, VTDHShape);
            return ArrayFloatD3ToDouble3D(vtdh3DMatrix, VTDHShape);
        } catch (Exception ex){
            errLog(ex);
            return new double[0][][];
        }
    }

    @Override
    public double[][][] getPeriod() {
        //Retrive the variable named "VTPK"
        Variable VTPK = dataFile.findVariable("VTPK");
        if(VTPK==null) {
            System.out.println("Can't find the variable named VTPK or tmn (wave period)");
            MyFileWriter debug = new MyFileWriter("", "debug", false, this.outDir);
            debug.WriteLog("waveVISIRINetCDFInputData: Can't find the variable named VTPK");
            debug.CloseFile();
            System.exit(0);
        }
        int[] VTPKShape = VTPK.getShape();
        int[] VTPKOrigin = new int[3];
        try {
            ArrayFloat.D3 vtpk3DMatrix = (ArrayFloat.D3) VTPK.read(VTPKOrigin, VTPKShape);
            return ArrayFloatD3ToDouble3D(vtpk3DMatrix, VTPKShape);
        } catch (Exception ex){
            errLog(ex);
            return new double[0][][];
        }
    }

    @Override
    public int[] getTimestep() {
        //Retrive the variable named "time"
        Variable time = dataFile.findVariable("time");
        if(time==null){
            System.out.println("Can't find the variable named time");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
            debug.WriteLog("waveVISIRINetCDFInputData: Can't find the variable named time");
            debug.CloseFile();
            System.exit(0);
        }
        int[] timeShape = time.getShape();
        int[] timeOrigin = new int[1];
        try {
            int[] parsedTime = new int[timeShape[0]];
            ArrayInt.D1 timeArray = (ArrayInt.D1) time.read(timeOrigin,timeShape);
            for(int i=0;i<timeShape[0];++i)
                parsedTime[i] = timeArray.get(i);
            return parsedTime;
        } catch (Exception ex){
            errLog(ex);
            return new int[0];
        }
    }

    @Override
    public void dispatch() {
        if(this.dataFile != null){
            try {
                this.dataFile.close();
            } catch (Exception ex){
                ex.printStackTrace();
            }
        }
    }

    private void errLog(Exception ex){
        MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
        debug.WriteLog("waveVISIRINetCDFInputData: "+ex.getMessage());
        debug.CloseFile();
        System.exit(0);
    }

    private double[][][] ArrayFloatD3ToDouble3D(ArrayFloat.D3 Matrix, int[] shape){
        double[][][] out = new double[shape[2]][shape[0]][shape[1]];
        for(int k=0;k<shape[2];++k){
            for(int i=0;i<shape[0];++i){
                for(int j=0;j<shape[1];++j){
                    if(Matrix.get(i,j,k) >= 1.0E20)
                        out[k][i][j] = Double.NaN;
                    else
                        out[k][i][j] = Matrix.get(i,j,k);
                }
            }
        }
        return out;
    }

}
