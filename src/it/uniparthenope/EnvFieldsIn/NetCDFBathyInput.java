package it.uniparthenope.EnvFieldsIn;

import it.uniparthenope.Debug.MyFileWriter;
import ucar.ma2.ArrayDouble;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

import java.io.IOException;
import java.util.ArrayList;

public class NetCDFBathyInput implements BathymetryInputData{
    private String filename;
    private NetcdfFile dataFile = null;
    private String outDir;

    public NetCDFBathyInput(String outDir){
        this.outDir=outDir;
    }

    @Override
    public ArrayList<Double> getLongitude() {
        //Retrive the variable named "longitude"
        Variable longitude = dataFile.findVariable("longitude");
        if(longitude == null){
            System.out.println("Can't find the variable named longitude");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
            debug.WriteLog("NetCDFBathyInput: Can't find the variable named longitude");
            debug.CloseFile();
            System.exit(0);
        }
        int[] longitudeShape = longitude.getShape();
        int[] longitudeOrigin = new int[1];
        try{
            ArrayDouble.D1 longitudeArray = (ArrayDouble.D1) longitude.read(longitudeOrigin, longitudeShape);
            ArrayList<Double> lon = new ArrayList<>();
            for(int i=0; i<longitudeShape[0];++i){
                lon.add(longitudeArray.get(i));
            }
            return lon;
        }catch (Exception ex){
            MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
            debug.WriteLog("NetCDFBathyInput-Longitude: "+ex.getMessage());
            debug.CloseFile();
            System.exit(0);
            return new ArrayList<>();
        }
    }

    @Override
    public ArrayList<Double> getLatitude() {
        //Retrive the variable named "latitude"
        Variable latitude = dataFile.findVariable("latitude");
        if(latitude==null){
            System.out.println("Can't find the variable named latitude");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
            debug.WriteLog("NetCDFBathyInput: Can't find the variable named longitude");
            debug.CloseFile();
            System.exit(0);
        }
        int[] latitudeShape = latitude.getShape();
        int[] latitudeOrigin = new int[1];
        try{
            ArrayDouble.D1 latitudeArray = (ArrayDouble.D1) latitude.read(latitudeOrigin,latitudeShape);
            ArrayList<Double> lat = new ArrayList<>();
            for(int i=0; i<latitudeShape[0]; ++i){
                lat.add(latitudeArray.get(i));
            }
            return lat;
        }catch (Exception ex){
            MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
            debug.WriteLog("NetCDFBathyInput-Latitude: "+ex.getMessage());
            debug.CloseFile();
            System.exit(0);
            return new ArrayList<>();
        }
    }

    @Override
    public double[][] getBathymetry() {
        Variable z = dataFile.findVariable("z");
        if(z == null){
            System.out.println("Can't find the variable named z");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
            debug.WriteLog("NetCDFBathyInput-Bathy: Can't find the variable named z");
            debug.CloseFile();
            System.exit(0);
        }
        int[] zShape = z.getShape();
        int[] zOrigin = new int[2];
        try {
            ArrayDouble.D2 zArray = (ArrayDouble.D2) z.read(zOrigin, zShape);
            double[][] depth = new double[zShape[0]][zShape[1]];
            for(int i=0;i<zShape[0];++i){
                for(int j=0;j<zShape[1];++j){
                    depth[i][j] = zArray.get(i, j);
                    if(depth[i][j]<=0)
                        depth[i][j] = Double.NaN;
                }
            }
            return depth;
        } catch (Exception ex){
            MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
            debug.WriteLog("NetCDFBathyInput-Bathy: "+ex.getMessage());
            debug.CloseFile();
            System.exit(0);
            return new double[0][];
        }
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
    public void dispatch() {
        if(this.dataFile != null){
            try {
                this.dataFile.close();
            } catch (Exception ex){
                ex.printStackTrace();
            }
        }
    }
}
