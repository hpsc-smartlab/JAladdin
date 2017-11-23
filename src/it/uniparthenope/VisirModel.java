package it.uniparthenope;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Scanner;

public class VisirModel {
    //input data
    private long bar_flag;
    private long timedep_flag;
    private Const constants;
    private Optim optim;
    private Ship ship;
    private Visualization visualization;
    private SpatialGrid sGrid;
    private TemporalGrid tGrid;
    private ExtremePoints extreme_pts;
    private DepartureParameters dep_datetime;
    private SafetyParameters safety;
    private EnvironmentalFields forcing;
    private Double[][] vel_LUT;
    private ArrayList<Double> H_array_m;
    private Double minWind;
    private Double maxWind;
    private Double polar_twa;
    private Double polar_tws;
    private ArrayList<Double> lon_ext;
    private ArrayList<Double> lat_ext;
    private ArrayList<Double> lon_int;
    private ArrayList<Double> lat_int;

    //approssimating computation time
    private long startTime;
    private long estimatedTime;


    //Constructor
    public VisirModel(long bar_flag, long timedep_flag){//Initialize with standard values defined in settings.m
        this.Tic();//Start the "timer"

        //bar_flag = 1: fresh compution of edges not crossing coastline (mode 1 of GMD-D paper)
        //bar_flag = 2: edges not crossing coastline read out from DB   (mode 2 of GMD-D paper)
        if((bar_flag!=1) && (bar_flag!=2)) { //If the input incorrect, set bar_flag = 2 as default value.
            this.bar_flag = 2;
        } else {
            this.bar_flag = bar_flag;
        }
        // timedep_flag=0: static algorithm with fields at time step #1
        // timedep_flag=2: time-dependent method
        if((timedep_flag!=0) && (timedep_flag!=2)) { //Same thing for timedep_flag.
            this.timedep_flag = 2;
        } else {
            this.timedep_flag = timedep_flag;
        }
        this.constants = new Const();
        this.optim = new Optim();
        this.ship = new Ship();
        this.sGrid = new SpatialGrid();
        this.tGrid = new TemporalGrid();
        this.visualization = new Visualization();
        this.lon_ext = new ArrayList<>();
        this.lon_int = new ArrayList<>();
        this.lat_ext = new ArrayList<>();
        this.lat_int = new ArrayList<>();
    }

    public void LoadData(){//Loading data parsing them from json file defined in inputFiles folder
        this.extreme_pts = new ExtremePoints();
        this.dep_datetime = new DepartureParameters();
        this.ship.LoadVesselParameters();
        this.safety = new SafetyParameters();
        this.optim.OptimizationParameters();
        this.visualization.VisualizationParameters();
    }

    public void CalculateParameters(){//Set some parameters based on other parameters
        this.ship.setStepsInPowerReduction(this.optim.getIntentional_speed_red());
        this.sGrid.setInvStepFields(this.optim.getWaveModel());
        this.forcing = new EnvironmentalFields(this.ship.getVessType(), this.ship.getSailType(),this.optim.getWindModel());
        this.safety.setCriteria();
        this.visualization.setData();
        if(this.forcing.getAnalytic() == 1){
            this.extreme_pts.setStart_lat(41.0);
            this.extreme_pts.setStart_lon(10.0);
            this.extreme_pts.setEnd_lat(40.5);
            this.extreme_pts.setEnd_lon(11.0);
            this.extreme_pts.setCycType("id");
            this.dep_datetime.setYear(14);
            this.dep_datetime.setMonth(1);
            this.dep_datetime.setDay(1);
            this.dep_datetime.setHour(0);
            this.dep_datetime.setMin(0);
            this.visualization.setEnv_forcing(1);
            this.visualization.setWaypoint_info(0);
            this.visualization.setSafegram(0);
            this.visualization.setScientific_mode(1);
            this.visualization.setH_s(1);
            this.visualization.setCustomMax(Double.NaN);
            this.ship.setVessType(1);
            this.timedep_flag = 0;
        }
    }

    private class Output{
        public ArrayList<Double> H_array_m;
        public Double[][] vel_LUT;
    }

    private void ship_Model(){//called from vessel_Response.m, ship_model.m implementation
        //Look-up table for involuntary ship speed reduction in waves.
        this.ship.setFn_max(this.constants.getMs2kts(), this.constants.getG0());
        //LUT independent variables:
        long Nh = 25; //40
        long Nl = 40; //30
        double Hmax = 8.0; //[m]
        long ih2 = (long) Math.floor(Nh/2);
        long il2 = (long) Math.floor(Nl/2);

        //Significant wave height
        //Adding to H_array_m (Nh-1) elements between 10^-1 and 10^log10(Hmax)
        ArrayList<Double> H_array_m = Utility.logspace(-1.0, Hmax, Nh-1);
        H_array_m.add(0, 0.0);//adding 0 at first element of H_array_m
        //preallocations:
        Double[][] vel_LUT = Utility.zeros((int) Nh+1, this.ship.getP_level_hp().size()+1);
        Double[][] Rc_LUT = Utility.zeros((int) Nh+1, this.ship.getP_level_hp().size()+1);
        Double[][] Raw_LUT = Utility.zeros((int) Nh+1, this.ship.getP_level_hp().size()+1);
        ArrayList<Double> P_level_thro = new ArrayList<Double>();
        Double max = Collections.max(this.ship.getP_level_hp());
        for(Double element : this.ship.getP_level_hp()){
            P_level_thro.add((100*element) / max);
        }
        //pars for v_Bowditch:
        Double m2ft = 3.2808399;
        Double a3_ref = 0.0248; //[kts/ ft^2]
        Double a2_ref = 0.0165; //[kts/ ft^2]
        Double a1_ref = 0.0083; //[kts/ ft^2]
        Double v0_ref = 18.0; //kts
        //LUT computation
        for(int ih = 1; ih< Nh; ih++){
            this.ship.ship_resistance(H_array_m.get(ih), this.constants);
            for(int j = 1; j<this.ship.getP_level_hp().size(); j++){
                vel_LUT[ih][j] = this.ship.getV_out_kts().get(j);
                Rc_LUT[ih][j] = this.ship.getR_c().get(j);
                Raw_LUT[ih][j] = this.ship.getR_aw().get(j);
                //v_Bowditch(ih)= max( 0, vel_LUT(1,1)/v0_ref* ( v0_ref - a2_ref *  (m2ft*H_array_m(ih)).^2) );
            }
        }
        this.ship.ship_resistance(0.0,this.constants);
        ArrayList<Double> v0 = this.ship.getV_out_kts();
        ArrayList<Double> Rc0 = this.ship.getR_c();
        ArrayList<Double> Raw0 = this.ship.getR_aw();
//        for ip=1: numel(P_level_hp)
//          disp([ min(vel_LUT(:,ip)), max(vel_LUT(:,ip)) ])
//        end
        //graphical output not implemented.
        this.H_array_m = H_array_m;
        this.vel_LUT = vel_LUT;
    }

    public void vessel_Response(){//vessel_Response.m implementation
        this.ship.shipPowerLevels();
        this.ship_Model();
        this.maxWind = Double.NaN;
        this.minWind = Double.NaN;
        this.polar_twa = Double.NaN;
        this.polar_tws = Double.NaN;
    }

    public void Grid_definition(){//Grid_definition.m implementation
        String freeedges_DB_filename = "freeedges_DB.dat";
        String freenodes_DB_filename = "freeNodes_DB.dat";
        // Bounding boxes:
        // there are 2 grids:
        // a reference grid (read from a DB) and an inset grid (defined by the user via a namelist)
        // following coordinate components may refer to each of them, depending on bar_flag (s. if-loop below):
        if(this.forcing.getAnalytic() == 1){
            this.extreme_pts.setBbox_deltaLat_u(0.1);
            this.extreme_pts.setBbox_deltaLon_l(0.1);
            this.extreme_pts.setBbox_deltaLat_d(0.1);
            this.extreme_pts.setBbox_deltaLon_r(0.1);
        }
        Double lat_max = this.extreme_pts.getBbox_deltaLat_u() + Math.max(this.extreme_pts.getStart_lat(), this.extreme_pts.getEnd_lat());
        Double lon_min = this.extreme_pts.getBbox_deltaLon_l() + Math.min(this.extreme_pts.getStart_lon(), this.extreme_pts.getEnd_lon());
        Double lat_min = this.extreme_pts.getBbox_deltaLat_d() + Math.min(this.extreme_pts.getStart_lat(), this.extreme_pts.getEnd_lat());
        Double lon_max = this.extreme_pts.getBbox_deltaLon_r() + Math.max(this.extreme_pts.getStart_lon(), this.extreme_pts.getEnd_lon());

        if(this.bar_flag == 1) {//fresh DB computation
            //bounding box of reference grid from namelist:
            this.sGrid.setDB_bbox__lat_max(lat_max);
            this.sGrid.setDB_bbox__lon_min(lon_min);
            this.sGrid.setDB_bbox__lat_min(lat_min);
            this.sGrid.setDB_bbox__lon_max(lon_max);
        } else if(this.bar_flag == 2){ //re-use DB information
            //bounding box of reference grid from free-edges DB header file:
            this.sGrid.setDB_bbox__lat_max(this.extreme_pts.getUL_lat());
            this.sGrid.setDB_bbox__lon_min(this.extreme_pts.getUL_lon());
            this.sGrid.setDB_bbox__lat_min(this.extreme_pts.getDR_lat());
            this.sGrid.setDB_bbox__lon_max(this.extreme_pts.getDR_lon());
            //bounding box of inset grid from namelist :
            this.sGrid.setBbox__lat_max(lat_max);
            this.sGrid.setBbox__lon_min(lon_min);
            this.sGrid.setBbox__lat_min(lat_min);
            this.sGrid.setBbox__lon_max(lon_max);
        }

        //Coastline:
        this.readout_coast();
        ArrayList<Double> y_coast = new ArrayList<>();
        ArrayList<Double> x_coast = new ArrayList<>();
        for(Double element : this.lat_int){//y_islands
            y_coast.add(element);
        }
        y_coast.add(Double.NaN);
        for(Double element : this.lat_ext){//y_continent
            y_coast.add(element);
        }
        for(Double element : this.lon_int){//x_islands
            x_coast.add(element);
        }
        x_coast.add(Double.NaN);
        for(Double element : this.lon_ext){//x_continent
            x_coast.add(element);
        }

        //Bathymetry:
        //read-out database:
        int bathy_code = 1;
        //#GM: here add code of compare_bathys if needed
        String bathy_title = "bathy: MedOneMin";
        if(bathy_code == 1){
            bathy_title = "bathy: MedOneMin";
        } else if(bathy_code == 2){
            bathy_title = "bathy: GEBCO-08";//warning: so far, I just downloaded a small subset of the full DB!
        } else if(bathy_code == 3) {
            bathy_title = "bathy: EMODnet";
        }//4: Adriatic 7.5 (not yet active)
        ArrayList<Object> bathymetry = this.readout_bathy(bathy_code);
        ArrayList<Double> lat_bathy = (ArrayList<Double>) bathymetry.get(0);
        ArrayList<Double> lon_bathy = (ArrayList<Double>) bathymetry.get(1);
        Double[][] z_bathy = (Double[][]) bathymetry.get(2);
        this.sGrid.setInv_step(1.0/Math.abs(lon_bathy.get(1)-lon_bathy.get(0)));
    }

    private void readout_coast(){
        //reads out coastline (_ext) and islands (_int) database - (source:  GSHHS via MEDSLIK-II)
        //This file contains a list of geographical coordinates of successive
        //points on the digitised coastline in a format similar to that used for
        //blanking files by the SURFER software. The format of each line in this
        //file is not important as long as adjacent entries are separated by a
        //space or comma. The first line contains the total number of closed
        //contours, counting all islands as well as the external coastline. The
        //second line contains the number of points on the external coastline
        //followed by the digit "0". There follows a list of longitudes and
        //latitudes of successive points on the external coastline, in decimal
        //degrees. These should be specified with adequate precision and it is
        //suggested that at least format (2f11.5) be used but for finely resolved
        //coastline maps the format (2f12.6) would be preferred. The external
        //coast must be described anti-clockwise and the last point must be
        //identical with the first. After this each island is listed: on the first
        //line the number of points on the island's coast and the digit "1", then
        //a list of longitudes and latitudes of points, this time described
        //clockwise. Again, for each island, the first point must be repeated as
        //the last one.
        try{
            int Nexternal = 0;
            ArrayList<Double> lat_tmp = new ArrayList<>();
            ArrayList<Double> lon_tmp = new ArrayList<>();
            //Opening file
            FileReader file = new FileReader("inputFiles/coast/medf.map");
            Scanner coastFile = new Scanner(file);
            //first line has a single column (i don't need)
            if(coastFile.hasNextInt()){
                coastFile.nextInt();
            }
            //in the second line, i need only 1st element, the 2nd element is discharged
            if(coastFile.hasNextInt()){
                Nexternal = coastFile.nextInt();
                int x = coastFile.nextInt();
            }
            int varType = 0;
            int totalLines = 0;//Total file lines (starting from 3rd line)
            while(coastFile.hasNext()){//while EOF
                if(varType%2==0){//1st element
                    lon_tmp.add(coastFile.nextDouble());//adding in temp longitude array
                }
                else{//2nd element
                    lat_tmp.add(coastFile.nextDouble());//adding in temp latitude array
                    totalLines++;
                }
                varType++;
            }
            for(int i=0;i<totalLines;i++){
                if(i<Nexternal){//adding in external coastline logitude and latitude array
                    //external coastline:  anti-clockwise and the last point is identical with the first
                    lon_ext.add(lon_tmp.get(i));
                    lat_ext.add(lat_tmp.get(i));
                }else {
                    //internal coastlines:  for each island, clockwise and the last point is identical with the first
                    lon_int.add(lon_tmp.get(i));
                    lat_int.add(lat_tmp.get(i));
                }
            }
            //closing file
            file.close();
        } catch(Exception e){
            e.printStackTrace();
        }
    }

    private ArrayList<Object> readout_bathy(int bathy_code){
        String filename ="";
        String varname = "";
        String lonname = "";
        String latname = "";
        switch(bathy_code){
            case 1:
                //MedOneMin:
                filename = "inputFiles/bathy/MedOneMin/med_one_min_single_vs2.nc";
                varname = "z";
                lonname = "longitude";
                latname = "latitude";
                break;
            case 2:
                //GEBCO_08:
                filename = "inputFiles/bathy/GEBCO/gebco_08_10_34_15_38.nc";
                varname = "z";
                lonname = "x_range";
                latname = "y_range";
                break;
            case 3:
                //EMODnet:
                filename = "inputFiles/bathy/EMODnet/Adriatic_Ionian_central__MedSea_mean.nc";
                varname = "depth_average";
                lonname = "lon";
                latname = "lat";
                break;
            default:
                System.out.println("unknow bathy DB code!");
                break;
        }
        //Parsing file:
        MyNetCDFParser test = new MyNetCDFParser(filename);
        ArrayList<Object> out = test.parseMedOneMin();
        if(out == null){
            System.out.println("Parsing fail!");
            return null;
        }
        ArrayList<Double> latTmp = (ArrayList<Double>) out.get(0);
        ArrayList<Double> lonTmp = (ArrayList<Double>) out.get(1);
        Double[][] zTmp = (Double[][]) out.get(2);
        //if sea depth >=0, setting as NaN
        for(int i =0 ;i<zTmp.length;i++){
            for(int j=0;j<zTmp[0].length;j++){
                if(zTmp[i][j]>=0){
                    zTmp[i][j]= Double.NaN;
                }
            }
        }
        //reduction within Med. Sea basin:
        Double lon_min= -6.5;
        Double lon_max= 37.0;
        Double lat_min= 30.0;
        Double lat_max= 46.0;
        //Check if coords are in the med. sea area
        ArrayList<Integer> latOkIndexes = new ArrayList<>();
        ArrayList<Integer> lonOkIndexes = new ArrayList<>();
        for(int i=0;i<latTmp.size();i++){
            if(latTmp.get(i)>=lat_min && latTmp.get(i)<=lat_max){
                latOkIndexes.add(i);
            }
        }
        for(int i=0;i<lonTmp.size();i++){
            if(lonTmp.get(i)>=lon_min && lonTmp.get(i)<=lon_max){
                lonOkIndexes.add(i);
            }
        }

        ArrayList<Double> lat = new ArrayList<>();
        ArrayList<Double> lon = new ArrayList<>();
        Double[][] z = new Double[latOkIndexes.size()][lonOkIndexes.size()];
        for(int i=0;i<latOkIndexes.size();i++){
            lat.add(latTmp.get(i));
        }
        for(int i=0;i<lonOkIndexes.size();i++){
            lon.add(lonTmp.get(i));
        }
        for(int i=0;i<latOkIndexes.size();i++){
            for(int j=0;j<lonOkIndexes.size();j++){
                z[i][j] = zTmp[latOkIndexes.get(i)][lonOkIndexes.get(j)];
            }
        }
        ArrayList<Object> output = new ArrayList<>();
        output.add((Object) lat);
        output.add((Object) lon);
        output.add((Object) z);
        return output;
    }

    private void Tic(){//Equivalent to MATLAB tic function
        this.startTime = System.nanoTime();
    }

    private void Toc(){
        this.estimatedTime = System.nanoTime() - this.startTime;
    }

    /*****************Getter methods*********************/
    public Const getConstants() {
        return this.constants;
    }

    public Optim getOptim() {
        return this.optim;
    }

    public Ship getShip() {
        return this.ship;
    }

    public SpatialGrid getsGrid() {
        return this.sGrid;
    }

    public TemporalGrid gettGrid() {
        return this.tGrid;
    }

    public ExtremePoints getExtreme_pts() {
        return extreme_pts;
    }

    public SafetyParameters getSafety() {
        return safety;
    }

    public DepartureParameters getDep_datetime() {
        return dep_datetime;
    }

    public Visualization getVisualization() {
        return visualization;
    }

    public long getBar_flag() {
        return bar_flag;
    }

    public long getTimedep_flag() {
        return timedep_flag;
    }

    public long getEstimatedTime() {
        return estimatedTime;
    }

    public Double[][] getVel_LUT() {
        return vel_LUT;
    }

    public ArrayList<Double> getH_array_m() {
        return H_array_m;
    }

    /******************************************************/
}
