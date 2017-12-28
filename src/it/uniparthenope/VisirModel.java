package it.uniparthenope;

import java.io.FileReader;
import java.sql.Time;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;
//Boxing classes
import com.sun.org.apache.xpath.internal.operations.Bool;
import it.uniparthenope.Boxing.*;
import it.uniparthenope.Debug.MyFileWriter;
import it.uniparthenope.Parser.MyBinaryParser;
import it.uniparthenope.Parser.MyCSVParser;
import it.uniparthenope.Parser.MyNetCDFParser;
import it.uniparthenope.Parser.MyTxtParser;

import javax.rmi.CORBA.Util;

public class VisirModel {
    //input data
    private long bar_flag;
    private long timedep_flag;
    private Fstats fstats;
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
    private double[][] vel_LUT;
    private ArrayList<Double> H_array_m;
    private double minWind;
    private double maxWind;
    private double polar_twa;
    private double polar_tws;
    private ArrayList<Double> lon_ext;
    private ArrayList<Double> lat_ext;
    private ArrayList<Double> lon_int;
    private ArrayList<Double> lat_int;
    private double[][] xy;
    private double[][] xg;
    private double[][] yg;
    private double[] xg_array;
    private double[] yg_array;
    private double[][] xy_DB;
    private ArrayList<Double> lat_bathy_Inset;
    private ArrayList<Double> lon_bathy_Inset;
    private double[][] bathy_Inset;
    private double[][] lsm_mask;
    private double[][] J_mask;
    ArrayList<Double> x_islands;
    ArrayList<Double> y_islands;
    ArrayList<Double> x_continent;
    ArrayList<Double> y_continent;
    private double estGdtDist;
    private MyFileWriter logFile;

    //approssimating computation time
    private long startTime;
    private long estimatedTime;


    //Constructor
    public VisirModel(long bar_flag, long timedep_flag){//Initialize with standard values defined in settings.m
        this.logFile = new MyFileWriter("","",false);
        this.logFile.WriteLine("");
        this.logFile.WriteLog("System initialization...");
        this.logFile.CloseFile();
        this.startTime = Utility.Tic();//Start the "timer"

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
        this.fstats = new Fstats();
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
        this.logFile = new MyFileWriter("","",true);
        this.logFile.WriteLog("Loading data...");
        this.logFile.CloseFile();
        this.extreme_pts = new ExtremePoints();
        this.dep_datetime = new DepartureParameters();
        this.ship.LoadVesselParameters();
        this.safety = new SafetyParameters();
        this.optim.OptimizationParameters();
        this.visualization.VisualizationParameters();
    }

    public void CalculateParameters(){//Set some parameters based on other parameters
        this.logFile = new MyFileWriter("","",true);
        this.logFile.WriteLog("Calculating parameters...");
        this.logFile.CloseFile();
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

    private void ship_Model(){//called from vessel_Response.m, ship_model.m implementation
        //Look-up table for involuntary ship speed reduction in waves.
        this.logFile = new MyFileWriter("","",true);
        this.logFile.WriteLog("Calculating ship model...");
        this.ship.setFn_max(this.constants.getMs2kts(), this.constants.getG0());
        this.logFile.CloseFile();
        //LUT independent variables:
        long Nh = 25; //40
        long Nl = 40; //30
        double Hmax = 8.0; //[m]
        long ih2 = (long) Math.floor(Nh/2);
        long il2 = (long) Math.floor(Nl/2);

        //Significant wave height
        //Adding to H_array_m (Nh-1) elements between 10^-1 and 10^log10(Hmax)
        ArrayList<Double> H_array_m = Utility.logspace(-1.0, Math.log10(Hmax), Nh-1);
        H_array_m.add(0, 0.0);//adding 0 at first element of H_array_m

        //preallocations:
        double[][] vel_LUT = new double[(int)Nh][this.ship.getP_level_hp().size()];
        double[][] Rc_LUT =  new double[(int)Nh][this.ship.getP_level_hp().size()];
        double[][] Raw_LUT =  new double[(int)Nh][this.ship.getP_level_hp().size()];
        ArrayList<Double> P_level_thro = new ArrayList<Double>();
        double max = Collections.max(this.ship.getP_level_hp());
        for(double element : this.ship.getP_level_hp()){
            P_level_thro.add((100*element) / max);
        }

        //pars for v_Bowditch:
        double m2ft = 3.2808399;
        double a3_ref = 0.0248; //[kts/ ft^2]
        double a2_ref = 0.0165; //[kts/ ft^2]
        double a1_ref = 0.0083; //[kts/ ft^2]
        double v0_ref = 18.0; //kts
        //LUT computation
        for(int ih = 0; ih< Nh; ih++){
            this.ship.ship_resistance(H_array_m.get(ih), this.constants);
            for(int j = 0; j<this.ship.getP_level_hp().size(); j++){
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
        this.logFile = new MyFileWriter("","",true);
        this.logFile.WriteLog("Calculating ship power levels...");
        this.ship.shipPowerLevels();
        this.logFile.CloseFile();
        this.ship_Model();
        this.maxWind = Double.NaN;
        this.minWind = Double.NaN;
        this.polar_twa = Double.NaN;
        this.polar_tws = Double.NaN;
    }

    public void Grid_definition(){//Grid_definition.m implementation
        this.logFile = new MyFileWriter("","",true);
        this.logFile.WriteLog("Grid definition: ");
        this.logFile.CloseFile();
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
        double lat_max = this.extreme_pts.getBbox_deltaLat_u() + Math.max(this.extreme_pts.getStart_lat(), this.extreme_pts.getEnd_lat());
        double lon_min = -this.extreme_pts.getBbox_deltaLon_l() + Math.min(this.extreme_pts.getStart_lon(), this.extreme_pts.getEnd_lon());

        double lat_min = -this.extreme_pts.getBbox_deltaLat_d() + Math.min(this.extreme_pts.getStart_lat(), this.extreme_pts.getEnd_lat());
        double lon_max = this.extreme_pts.getBbox_deltaLon_r() + Math.max(this.extreme_pts.getStart_lon(), this.extreme_pts.getEnd_lon());
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
        this.logFile = new MyFileWriter("","",true);
        this.logFile.WriteLog("\tLoading coastline data...");
        this.logFile.CloseFile();
        //Coastline:
        this.readout_coast();
        ArrayList<Double> y_coast = new ArrayList<>();
        ArrayList<Double> x_coast = new ArrayList<>();
        for(double element : this.lat_int){//y_islands
            y_coast.add(element);
        }
        y_coast.add(Double.NaN);
        for(double element : this.lat_ext){//y_continent
            y_coast.add(element);
        }
        for(double element : this.lon_int){//x_islands
            x_coast.add(element);
        }
        x_coast.add(Double.NaN);
        for(double element : this.lon_ext){//x_continent
            x_coast.add(element);
        }
        this.logFile = new MyFileWriter("","",true);
        this.logFile.WriteLog("\tLoading bathymetry data...");
        this.logFile.CloseFile();
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
        readout_bathyResults bathymetry = this.readout_bathy(bathy_code);
        ArrayList<Double> lat_bathy = bathymetry.getLat();
        ArrayList<Double> lon_bathy = bathymetry.getLon();
        double[][] z_bathy = bathymetry.getZ();
        this.sGrid.setInv_step( Math.round(1.0/Math.abs(lon_bathy.get(1)-lon_bathy.get(0)))+0.0 );
        this.logFile = new MyFileWriter("","",true);
        this.logFile.WriteLog("\tCalculating target grid and reduced bathy field...");
        this.logFile.CloseFile();
        //Target grid and reduced bathy field:
        //grid_extreme_coords
        grid_extreme_coordsResults insets = this.grid_extreme_coords(lat_bathy, lon_bathy,z_bathy);
        this.lat_bathy_Inset = insets.getLat_red();
        this.lon_bathy_Inset = insets.getLon_red();
        this.bathy_Inset = insets.getField_out();

        ArrayList<Boolean> x_bool=new ArrayList<>();
        ArrayList<Boolean> y_bool=new ArrayList<>();
        ArrayList<Boolean> inset_bool=new ArrayList<>();
        ArrayList<Double> x_coast_Inset = new ArrayList<>();
        ArrayList<Double> y_coast_Inset = new ArrayList<>();
        if(this.bar_flag == 2){
            // coastline excerpt within Inset:
            this.logFile = new MyFileWriter("","",true);
            this.logFile.WriteLog("\tCalculating coastline excerpt within inset...");
            this.logFile.CloseFile();
            double min_lon_bathy = Collections.min(this.lon_bathy_Inset);
            double max_lon_bathy = Collections.max(this.lon_bathy_Inset);
            double min_lat_bathy = Collections.min(this.lat_bathy_Inset);
            double max_lat_bathy = Collections.max(this.lat_bathy_Inset);
            for(double element : x_coast){
                if((element >= min_lon_bathy) && (element <= max_lon_bathy)){
                    x_bool.add(true);
                } else { x_bool.add(false); }
            }
            for(double element : y_coast){
                if((element >= min_lat_bathy) && (element <= max_lat_bathy)){
                    y_bool.add(true);
                } else { y_bool.add(false); }
            }
            for(int i =0 ;i<x_bool.size();i++){
                if(x_bool.get(i) && y_bool.get(i)){
                    inset_bool.add(true);
                } else { inset_bool.add(false); }
            }
            for(int i =0 ;i<inset_bool.size();i++){
                if(inset_bool.get(i)){
                    x_coast_Inset.add(x_coast.get(i));
                    y_coast_Inset.add(y_coast.get(i));
                }
            }
        }

        //ref. grid coordinates:
        ArrayList<Double> lat_bathy_DB = Utility.linspace(this.sGrid.getDB_yi(), this.sGrid.getDB_yf(), this.sGrid.getDB_Ny());
        ArrayList<Double> lon_bathy_DB = Utility.linspace(this.sGrid.getDB_xi(), this.sGrid.getDB_xf(), this.sGrid.getDB_Nx());


        mdata_gridResults data_gridOut = mdata_grid(lat_bathy_DB,lon_bathy_DB);
        this.xy_DB = data_gridOut.getXy();
        double[][] xg_DB = data_gridOut.getXg();
        double[][] yg_DB = data_gridOut.getYg();

        //inset grid plaid coordinates:
        this.logFile = new MyFileWriter("","",true);
        this.logFile.WriteLog("\tinset grid plaid coordinates...");
        this.logFile.CloseFile();
        mdata_gridResults data_gridOut2 = mdata_grid(this.lat_bathy_Inset, this.lon_bathy_Inset);
        this.xy = data_gridOut2.getXy();
        this.xg = data_gridOut2.getXg();
        this.yg = data_gridOut2.getYg();
        if(this.visualization.getGraphData()==1){
            //csv_write xy
            this.logFile = new MyFileWriter("","",true);
            this.logFile.WriteLog("\tWriting graph coords in GRAPH.node_LonLat.csv...");
            this.logFile.CloseFile();
            try{
                MyCSVParser csv = new MyCSVParser("Output/GRAPH.node_LonLat.csv");
                csv.writeCSV(this.xy);
            } catch (Exception e){
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("grid_definition, csv parser: "+e.getMessage());
                debug.CloseFile();
                e.printStackTrace();
            }
        }

        //--------------------------------------------------------------------
        //Joint coast-vessel safety mask:
        this.logFile = new MyFileWriter("","",true);
        this.logFile.WriteLog("\tcomputation of a joint coast-vessel safety mask...");
        this.logFile.CloseFile();
        double[][] UKC_mask = Utility.NaNmatrix(bathy_Inset.length, bathy_Inset[0].length);
        for(int i=0;i<UKC_mask.length;i++){
            for(int j=0;j<UKC_mask[0].length;j++){
                if(bathy_Inset[i][j]>this.ship.getDraught()){
                    UKC_mask[i][j] = 1.0;
                }
            }
        }
        double[][] xg_Jmasked;
        double[][] yg_Jmasked;
        double[][] J_bathy_Inset;
        this.xg_array = new double[0];
        this.yg_array = new double[0];
        double[][] xy_g;
        double[][] dist_mask;
        double[][] min_coast_dist;
        if(this.bar_flag == 2){
            //readout free nodes DB:
            this.logFile = new MyFileWriter("","",true);
            this.logFile.WriteLog("\treadout freenodes DB...");
            this.logFile.CloseFile();
            MyBinaryParser datFile = new MyBinaryParser("inputFiles/graph/freeNodes_DB.dat");
            long[] free_nodes_DB = datFile.readAsUInt32();
            //remapping free nodes:
            this.logFile = new MyFileWriter("","",true);
            this.logFile.WriteLog("\tremapping free nodes...");
            this.logFile.CloseFile();
            long lambda = Math.round((this.sGrid.getXi()-this.sGrid.getDB_xi())*this.sGrid.getInv_step());
            long mu = Math.round((this.sGrid.getYi()-this.sGrid.getDB_yi())*this.sGrid.getInv_step());
            idx_ref2inset_gridResults out = this.idx_ref2inset_grid(free_nodes_DB,this.sGrid.getDB_Nx(), this.sGrid.getDB_Ny(),this.sGrid.getInset_Nx(),this.sGrid.getInset_Ny(),lambda,mu);
            long[] row = out.getRow();
            long[] col = out.getCol();
            long[] free_nodes_extended = out.getIdx();
            ArrayList<Long> tmp = new ArrayList<>();
            for(int i=0;i<free_nodes_extended.length;i++){
                if(free_nodes_extended[i]>=0){
                    tmp.add(free_nodes_extended[i]);
                }
            }
            free_nodes_extended = new long[tmp.size()];
            for(int i=0;i<tmp.size();i++){
                free_nodes_extended[i]=tmp.get(i);
            }
            long[] free_nodes = free_nodes_extended;

            //checks:
            long free_nodes_Number = free_nodes.length;
            if(free_nodes_Number==0){
                System.out.println("no free nodes found in graph");
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("Grid_definition: no free nodes found in graph");
                debug.CloseFile();
                System.exit(0);
            } else if(free_nodes_Number > this.sGrid.getNodesLargeN()){
                System.out.println("too large graph");
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("Grid_definition: too large graph");
                debug.CloseFile();
                System.exit(0);
            }
            this.sGrid.setFreenodes(free_nodes_Number);
            //lsm, created using coastline DB (NaNs on landmass):
            this.logFile = new MyFileWriter("","",true);
            this.logFile.WriteLog("\tCreating lsm_mask using coastline DB...");
            this.logFile.CloseFile();
            lsm_mask = Utility.NaNmatrix(xg.length,xg[0].length);
            int nRows = xg.length;
            int ii=0;
            int jj=0;
            for(int k=0;k<free_nodes.length;k++){
                int index =(int) free_nodes[k]-1;
                jj=index/nRows;
                ii=(index%nRows);
                lsm_mask[ii][jj]=1.0;
            }
            int[] dim = new int[2];
            dim[0]=xg.length;
            dim[1]=xg[0].length;

            //Safe distance from coastline:
            this.logFile = new MyFileWriter("","",true);
            this.logFile.WriteLog("\tCalculating safe distance from coastline...");
            this.logFile.CloseFile();
            xg_array = Utility.reshape(xg,dim[0]*dim[1]);
            yg_array = Utility.reshape(yg, yg.length*yg[0].length);
            xy_g = new double[xg_array.length][2];
            for(int i=0;i<xg_array.length;i++){
                xy_g[i][0] = xg_array[i];
                xy_g[i][1] = yg_array[i];
            }
            int nC = x_coast_Inset.size();
            int cols = dim[0]*dim[1];
            if(nC>0){
                //Double[][] coast_dist = Utility.zeros(nC, cols);
                double[][] coast_dist = new double[nC][cols];
                double[][] P_b = new double[1][2];
                for(int i=0;i<nC;i++){
                    P_b[0][0]=x_coast_Inset.get(i);
                    P_b[0][1]=y_coast_Inset.get(i);
                    double[] hor_dist = hor_distance("s",xy_g, P_b);
                    for(int j=0; j<cols;j++){
                        coast_dist[i][j]= hor_dist[j];
                    }
                }
                double[] min_coast_distTmp = Utility.min(coast_dist,1);
                min_coast_dist = Utility.reshape(min_coast_distTmp, dim);
                dist_mask = Utility.NaNmatrix(dim[0],dim[1]);
                for(int i=0;i<dist_mask.length;i++){
                    for(int j=0;j<dist_mask[0].length;j++){
                        if(min_coast_dist[i][j]>=this.extreme_pts.getMinCoastDist()){
                            dist_mask[i][j]=1.0;
                        }
                    }
                }

            } else{
                dist_mask = Utility.ones(dim[0],dim[1]);
            }
            this.logFile = new MyFileWriter("","",true);
            this.logFile.WriteLog("\tCalculating Joint safe mask...");
            this.logFile.CloseFile();
            // %###############################################
            // %
            // %  Joint safe mask:
            // %
            // J_mask = lsm_mask' .* UKC_mask .* dist_mask' ;
            // %
            // %###############################################
            J_mask = Utility.MatrixComponentXcomponent(Utility.MatrixComponentXcomponent(Utility.transposeMatrix(lsm_mask),UKC_mask),Utility.transposeMatrix(dist_mask));

            xg_Jmasked = Utility.MatrixComponentXcomponent(xg,Utility.transposeMatrix(J_mask));
            yg_Jmasked = Utility.MatrixComponentXcomponent(yg,Utility.transposeMatrix(J_mask));
            J_bathy_Inset = Utility.MatrixComponentXcomponent(bathy_Inset,J_mask);
            // % Note: It is not necessary to pass [xg_Jmasked,yg_Jmasked] to the MAIN.m in place of [xg,yg].
            // % Indeed, shortest path search already accounts (s. Edges_definition.m)
            // % for coastline (free_edges) and bathymetry/lsm/dist-from-coastline (nogo_edges)
            // % [xg_Jmasked,yg_Jmasked] are used here just for the sake of searching departure and arrival nodes (s. below).
            //-----------------------------------------------------------------------------------------------------
            //Start and end nodes:
            xg_array = Utility.reshape(xg_Jmasked,xg_Jmasked.length*xg_Jmasked[0].length);
            yg_array = Utility.reshape(yg_Jmasked, yg_Jmasked.length*yg_Jmasked[0].length);
            xy_g = new double[xg_Jmasked.length*xg_Jmasked[0].length][2];
            for(int i =0 ;i<xg_Jmasked.length*xg_Jmasked[0].length;i++){
                xy_g[i][0]=xg_array[i];
                xy_g[i][1]=yg_array[i];
            }
            this.logFile = new MyFileWriter("","",true);
            this.logFile.WriteLog("\tCalculating distances from start/end nodes...");
            this.logFile.CloseFile();
            //distance from start/end nodes:
            double[][] tmpP_b=new double[1][2];
            tmpP_b[0][0]=this.extreme_pts.getStart_lon();
            tmpP_b[0][1]=this.extreme_pts.getStart_lat();
            double[] start_dist_matrix = hor_distance("s",xy_g,tmpP_b);
            tmpP_b[0][0]=this.extreme_pts.getEnd_lon();
            tmpP_b[0][1]=this.extreme_pts.getEnd_lat();
            double[] end_dist_matrix = hor_distance("s",xy_g,tmpP_b);

            this.sGrid.setMin_start_dist(Utility.min(start_dist_matrix));
            this.sGrid.setMin_end_dist(Utility.min(end_dist_matrix));

            if(Double.isNaN(this.sGrid.getMin_start_dist()) || Double.isNaN(this.sGrid.getMin_end_dist())){
                System.out.println("departure or arrival point not compliant with safety specifications");
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("Grid_definition: departure or arrival point not compliant with safety specifications");
                debug.CloseFile();
                System.exit(0);
            }
            ArrayList<Integer> tmpIndexes = new ArrayList<>();
            for(int i=0;i<start_dist_matrix.length;i++){
                if(start_dist_matrix[i]==this.sGrid.getMin_start_dist()){
                    tmpIndexes.add(i);
                }
            }
            this.sGrid.setNode_start(Collections.min(tmpIndexes));
            tmpIndexes = new ArrayList<>();
            for(int i=0;i<end_dist_matrix.length;i++){
                if(end_dist_matrix[i]==this.sGrid.getMin_end_dist()){
                    tmpIndexes.add(i);
                }
            }
            this.sGrid.setNode_end(Collections.min(tmpIndexes));

            this.sGrid.setNode_start_lat(yg_array[(int) this.sGrid.getNode_start()]);
            this.sGrid.setNode_start_lon(xg_array[(int) this.sGrid.getNode_start()]);
            this.sGrid.setNode_end_lat(yg_array[(int) this.sGrid.getNode_end()]);
            this.sGrid.setNode_end_lon(xg_array[(int) this.sGrid.getNode_end()]);

            double[] p_a = new double[2];
            double[] p_b = new double[2];
            p_a[0] = this.sGrid.getNode_start_lon();
            p_a[1] = this.sGrid.getNode_start_lat();
            p_b[0] = this.sGrid.getNode_end_lon();
            p_b[1] = this.sGrid.getNode_end_lat();
            estGdtDist = hor_distance("s", p_a, p_b);

            p_a[0] = this.sGrid.getNode_end_lon();
            p_a[1] = this.sGrid.getNode_end_lat();
            p_b[0] = this.sGrid.getNode_end_lon();
            p_b[1] = this.sGrid.getNode_start_lat();
            double delta_y = hor_distance("s", p_a, p_b);
            int sign = Utility.sign(this.sGrid.getNode_end_lat()-this.sGrid.getNode_start_lat());
            delta_y=delta_y*sign;

            p_a[0] = this.sGrid.getNode_end_lon();
            p_a[1] = this.sGrid.getNode_start_lat();
            p_b[0] = this.sGrid.getNode_start_lon();
            p_b[1] = this.sGrid.getNode_start_lat();
            double delta_x = hor_distance("s", p_a, p_b);
            sign = Utility.sign(this.sGrid.getNode_end_lon()-this.sGrid.getNode_start_lon());
            delta_x=delta_x*sign;

            //orientation of the gdt route
            this.logFile = new MyFileWriter("","",true);
            this.logFile.WriteLog("\tCalculating orientation of the geodetic route...");
            this.logFile.CloseFile();
            double atan2 = Double.NaN;
            atan2=Math.atan2(delta_x,delta_y);
            this.sGrid.setTheta_gdt(atan2);
            //th= Sgrid.theta_gdt/const.deg2rad

        } else if(this.bar_flag==1){
            lsm_mask = Utility.NaNmatrix(xg.length,xg[0].length);
            double[][] tmp = Utility.transposeMatrix(lsm_mask);
            J_mask = Utility.NaNmatrix(tmp.length,tmp[0].length);
            //-----------------------------------------------------------------------------------------------------
            //Start and end nodes:
            xg_array = new double[0];
            yg_array = new double[0];
            estGdtDist = Double.NaN;
        }
        this.x_islands = this.lon_int;
        this.y_islands = this.lat_int;
        this.x_continent = this.lon_ext;
        this.y_continent = this.lat_ext;
        this.logFile = new MyFileWriter("","",true);
        this.logFile.WriteLog("Done.");
        this.logFile.CloseFile();
    }

    private double changeDirRule(double inField){
        // % change of wind/current directional convention:
        // %
        // % from atan2 output concention to WAM-like convention, i.e.:
        // % from:  [-pi, pi] counterclockwise with 0 at  3:00 o'clock
        // %   to:  [0, 2*pi]        clockwise with 0 at 12:00 o'clock
        // %
        // % (inField and outField must be in radians)
        // %
        //---------------------------------------------
        double inField1 = inField;
        double inField1_ = 0;
        if(inField1 < 0){
            inField1_ = inField1;
            inField1_ = inField1_ + 2*Math.PI;
            inField1 = inField1_;
        }
        //---------------------------------------------
        double inField2 = inField1 -Math.PI/2;
        double inField2_ = 0;
        if(inField2 < 0){
            inField2_ = inField2;
            inField2_ = inField2_ + 2*Math.PI;
            inField2 = inField2_;
        }
        //---------------------------------------------
        return (-inField2+(2*Math.PI));
    }

    private double[] changeDirRule(double[] inField){
        // % change of wind/current directional convention:
        // %
        // % from atan2 output concention to WAM-like convention, i.e.:
        // % from:  [-pi, pi] counterclockwise with 0 at  3:00 o'clock
        // %   to:  [0, 2*pi]        clockwise with 0 at 12:00 o'clock
        // %
        // % (inField and outField must be in radians)
        // %
        //---------------------------------------------
        double[] inField1 = inField;
        for(int i=0;i<inField1.length;i++){
            if(inField1[i]<0){
                inField1[i]+=(2*Math.PI);
            }
        }

        //---------------------------------------------
        double[] inField2 = inField1;
        for(int i=0;i<inField2.length;i++){
            inField2[i]-=(Math.PI/2);
            if(inField2[i]<0){
                inField2[i]+=(2*Math.PI);
            }
        }
        //---------------------------------------------
        double[] outFiled = new double[inField2.length];
        for(int i=0;i<inField2.length;i++){
            outFiled[i] = -inField2[i]*(2*Math.PI);
        }
        return outFiled;
    }

    public double[][][] changeDirRule(double[][][] inField){
        // % change of wind/current directional convention:
        // %
        // % from atan2 output concention to WAM-like convention, i.e.:
        // % from:  [-pi, pi] counterclockwise with 0 at  3:00 o'clock
        // %   to:  [0, 2*pi]        clockwise with 0 at 12:00 o'clock
        // %
        // % (inField and outField must be in radians)
        // %
        //---------------------------------------------
        double[][][] inField1 = inField;
        for(int i=0;i<inField.length;i++){
            for(int j=0;j<inField[0].length;j++){
                for(int k=0;k<inField[0][0].length;k++){
                    if(inField1[i][j][k]<0){
                        inField1[i][j][k]+=(2*Math.PI);
                    }
                }
            }
        }
        //---------------------------------------------
        double[][][] inField2 = inField1;
        for(int i=0;i<inField2.length;i++){
            for(int j=0;j<inField2[0].length;j++){
                for(int k=0;k<inField2[0][0].length;k++){
                    inField2[i][j][k] -=(Math.PI/2);
                    if(inField2[i][j][k]<0){
                        inField2[i][j][k] += (2*Math.PI);
                    }
                }
            }
        }
        //---------------------------------------------
        double[][][] outField = new double[inField2.length][inField2[0].length][inField2[0][0].length];
        for(int i=0;i<inField.length;i++){
            for(int j=0;j<inField[0].length;j++){
                for(int k=0;k<inField[0][0].length;k++){
                    outField[i][j][k] = -inField2[i][j][k] * (2*Math.PI);
                }
            }
        }
        return outField;
    }

    private double[] hor_distance(String method, double[][] P_a, double[][] P_b, double... varargin){
        // % horizontal distance between pair of points (expressed in NM)
        // % either on the plane or the sphere
        // % works also with arrays
        // % optionally rescales distances by grid_step_in_NM
        // %
        // % P_a and P_b must be linear array of the same size
        // % (build via meshgrid and then reshape)
        // %
        // % P_a: [lon, lat]
        // % P_b: [lon, lat]
        // %
        // % method='p' for planar geometry
        // % method='s' for spherical geometry
        this.constants.setDeg2rad(Math.PI/180);
        double grid_step_in_NM = varargin.length > 0 ? varargin[0] : 1.0;

        double[] xa = new double[P_a.length];
        double[] ya = new double[P_a.length];
        for(int i=0;i<P_a.length;i++){
            xa[i]=P_a[i][0];
            ya[i]=P_a[i][1];
        }

        double[] xb = new double[P_b.length];
        double[] yb = new double[P_b.length];
        for(int i=0;i<P_b.length;i++){
            xb[i]=P_b[i][0];
            yb[i]=P_b[i][1];
        }


        double[] dd=new double[P_a.length];
        if(method=="plane"||method=="p"){
            for(int i=0;i<P_a.length;i++){
                dd[i]=grid_step_in_NM*Math.sqrt(Math.pow((xa[i]-xb[0]),2) + Math.pow((ya[i]-yb[0]),2));
            }
        } else if(method=="sphere" || method=="s"){
            // % from : http://mathworld.wolfram.com/GreatCircle.html
            // % x and y must be in degree.
            // % output in Nautical Miles (NM)
            double E_radius = 3444.0; //NM
            for(int i=0;i<P_a.length;i++){
                dd[i]=E_radius * Math.acos(Math.cos(this.constants.getDeg2rad()*ya[i])*Math.cos(this.constants.getDeg2rad()*yb[0])*
                        Math.cos(this.constants.getDeg2rad()*(xa[i]-xb[0]))+Math.sin(this.constants.getDeg2rad()*ya[i])*
                        Math.sin(this.constants.getDeg2rad()*yb[0]));
            }
        }

        return dd;
    }

    private double hor_distance(String method, double[] P_a, double[] P_b, double... varargin){
        // % horizontal distance between pair of points (expressed in NM)
        // % either on the plane or the sphere
        // % works also with arrays
        // % optionally rescales distances by grid_step_in_NM
        // %
        // % P_a and P_b must be linear array of the same size
        // % (build via meshgrid and then reshape)
        // %
        // % P_a: [lon, lat]
        // % P_b: [lon, lat]
        // %
        // % method='p' for planar geometry
        // % method='s' for spherical geometry
        this.constants.setDeg2rad(Math.PI/180);
        double grid_step_in_NM = varargin.length > 0 ? varargin[0] : 1.0;

        double xa = P_a[0];
        double ya = P_a[1];

        double xb = P_b[0];
        double yb = P_b[1];


        double dd=0;
        if(method=="plane"||method=="p"){
            dd=grid_step_in_NM*Math.sqrt(Math.pow((xa-xb),2) + Math.pow((ya-yb),2));
        } else if(method=="sphere" || method=="s"){
            // % from : http://mathworld.wolfram.com/GreatCircle.html
            // % x and y must be in degree.
            // % output in Nautical Miles (NM)
            double E_radius = 3444.0; //NM
            dd=E_radius * Math.acos(Math.cos(this.constants.getDeg2rad()*ya)*Math.cos(this.constants.getDeg2rad()*yb)*
                    Math.cos(this.constants.getDeg2rad()*(xa-xb))+Math.sin(this.constants.getDeg2rad()*ya)*
                    Math.sin(this.constants.getDeg2rad()*yb));
        }

        return dd;
    }

    private idx_ref2inset_gridResults idx_ref2inset_grid(long[] idx_big, long nx_big, long ny_big, long nx, long ny, long lambda, long mu){
        // % remap grid index
        // % from a nx_big columns grid (reference) to a nx columns grid (inset).
        // % (lambda, mu) are the offset coordinates of idx=1 gridpoint of the inset grid.
        // % (row, col) are Cartesian coordinates in the inset grid.
        // %
        // % --> Rectangular and equally spaced grids are assumed! <--
        // % In picture below, indexes between brackets refer to the inset grid.
        // %
        // % WARNING: idx outside inset grid are set to -1!
        // %   ny_big- -------------------------------------------------
        // %           |                                               |
        // %           |     reference grid                            |
        // %           |                                               |
        // %           |                                               |
        // %           |                                               |
        // %           |                                               |
        // %    (ny) - |..   .....    ......    .. |----------|        |
        // %           |                           |          |        |
        // %           |                           |  inset   |        |
        // %           |                           |          |        |
        // %(1) 1+mu - |..   .....    ......    .. |----------|        |
        // %           |                           .          .        |
        // %           |                           .          .        |
        // %       1 - ------------------------------------------------|
        // %           |                           |          .        .
        // %
        // %           1                      (1) 1+lambda   (nx)     nx_big
        // %
        //checks:
        if(nx > nx_big || ny> ny_big || lambda+nx > nx_big || mu+ny > ny_big){
            System.out.println("inset grid not within reference grid!");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("idx_ref2inset_grid: inset grid not within reference grid!");
            debug.CloseFile();
            System.exit(0);
        }
        if((Utility.any(idx_big,"<",1)) || Utility.any(idx_big,">",(nx_big*ny_big))){
            System.out.println("grid index not within reference grid!");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("idx_ref2inset_grid: grid index not within reference grid!");
            debug.CloseFile();
            System.exit(0);
        }
        //col idx (ref):
        long[] col_big = new long[idx_big.length];
        for(int i =0; i<idx_big.length;i++){
            col_big[i]=Math.floorMod(idx_big[i],nx_big);
            if(col_big[i]==0){
                col_big[i]=nx_big;
            }
        }
        //row idx (ref):
        long[] row_big = new long[idx_big.length];
        for(int i =0;i<idx_big.length;i++){
            row_big[i] =Math.round(1+(idx_big[i]-col_big[i])/nx_big);
        }
        //col idx (inset):
        long[] col = new long[col_big.length];
        for(int i=0;i<col.length;i++){
            col[i]=col_big[i]-lambda;
        }
        //row idx (inset):
        long[] row = new long[row_big.length];
        for(int i=0;i<row.length;i++){
            row[i]=row_big[i]-mu;
        }
        //idx (inset grid):
        // and pad with zero indexes out of inset grid:
        long[] idx = new long[col.length];
        for(int i=0;i<idx.length;i++){
            idx[i]=col[i]+nx*(row[i]-1);
            if( col[i] < 1 || col[i] > nx || row[i] < 1 || row[i] > ny){
                idx[i] = -1;
            }
        }

        return  new idx_ref2inset_gridResults(row, col, idx);
    }

    private mdata_gridResults mdata_grid(ArrayList<Double> lat, ArrayList<Double> lon){
        //converts 2 lat-lon 1-dimensional arrays (Nlat,1) and (Nlon,1) into :
        //a) a list of node coordinates (Nlat*Nlon, 2)
        //b) a meshgrid matrix       (Nlat,   Nlon)
        int NN = lat.size() * lon.size();
        meshgridResults out = Utility.meshgrid(lat,lon);
        double[][] yg = out.getX();
        double[][] xg = out.getY();
        int[] aDim = new int[2];
        aDim[0]=NN;
        aDim[1]=1;
        double[][] xa = Utility.reshape(xg, aDim);
        double[][] ya = Utility.reshape(yg, aDim);
        double[][] xy = new double[NN][2];
        for(int i=0;i<NN;i++){
            xy[i][0]=xa[i][0];
            xy[i][1]=ya[i][0];
        }
        return new mdata_gridResults(xy,xg,yg);
    }

    private grid_extreme_coordsResults grid_extreme_coords(ArrayList<Double> lat, ArrayList<Double> lon, double[][] field_in){
        //compute coordinates of grid extreme nodes,
        //both for the reference grid (grid read from DB)
        //and the inset grid

        //reference grid:
        //Finding indexes
        ArrayList<Integer> DB_lat_row = new ArrayList<>();
        ArrayList<Integer> DB_lon_row = new ArrayList<>();
        for(int i=0;i<lat.size();i++){
            if( (lat.get(i) >= this.sGrid.getDB_bbox__lat_min()) && (lat.get(i) <= this.sGrid.getDB_bbox__lat_max()) ){
                DB_lat_row.add(i);
            }
        }

        for(int i=0;i<lon.size();i++){
            if( (lon.get(i) >= this.sGrid.getDB_bbox__lon_min()) && (lon.get(i) <= this.sGrid.getDB_bbox__lon_max()) ){
                DB_lon_row.add(i);
            }
        }

        this.sGrid.setDB_xi(lon.get(DB_lon_row.get(0)));
        this.sGrid.setDB_yi(lat.get(DB_lat_row.get(0)));
        this.sGrid.setDB_xf(lon.get(DB_lon_row.get(DB_lon_row.size()-1)));
        this.sGrid.setDB_yf(lat.get(DB_lat_row.get(DB_lat_row.size()-1)));

        this.sGrid.setDB_Nx(1 + DB_lon_row.get(DB_lon_row.size()-1) - DB_lon_row.get(0));
        this.sGrid.setDB_Ny(1 + DB_lat_row.get(DB_lat_row.size()-1) - DB_lat_row.get(0));

        //double meshRes = 1/this.sGrid.getInvStepFields();

        //*****not sure of this*****
        ArrayList<Double> lat_red = new ArrayList<>();
        ArrayList<Double> lon_red = new ArrayList<>();
        double[][] field_out = null;
        //*****not sure of this*****
        //if barflag ==1, we are in creating graph DB mode, so we return empty lat_red, lon_red and field_out
        if(bar_flag == 2){//reading DB mode
            //Insert grid
            ArrayList<Integer> lat_row = new ArrayList<>();
            ArrayList<Integer> lon_row = new ArrayList<>();
            for(int i=0;i<lat.size();i++){
                if( (this.sGrid.getBbox__lat_min() <= lat.get(i)) && (lat.get(i) <= this.sGrid.getBbox__lat_max()) ){
                    lat_row.add(i);
                }
            }
            for(int i=0;i<lon.size();i++){
                if((this.sGrid.getBbox__lon_min() <= lon.get(i)) && (lon.get(i) <= this.sGrid.getBbox__lon_max())){
                    lon_row.add(i);
                }
            }
            this.sGrid.setXi(lon.get(lon_row.get(0)));
            this.sGrid.setYi(lat.get(lat_row.get(0)));
            this.sGrid.setXf(lon.get(lon_row.get(lon_row.size()-1)));
            this.sGrid.setYf(lat.get(lat_row.get(lat_row.size()-1)));
            this.sGrid.setInset_area( Math.abs( (this.sGrid.getXf()-this.sGrid.getXi()) * (this.sGrid.getYf()-this.sGrid.getYi()) ) );

            this.sGrid.setInset_Nx(1 + lon_row.get(lon_row.size()-1) - lon_row.get(0));
            this.sGrid.setInset_Ny(1+ lat_row.get(lat_row.size()-1) - lat_row.get(0));

            //reduced field coordinates and values:
            for(int element : lat_row){
                lat_red.add(lat.get(element));
            }
            for(int element : lon_row){
                lon_red.add(lon.get(element));
            }
            field_out = new double[lat_row.size()][lon_row.size()];
            for(int i=0;i<lat_row.size();i++){
                for(int j=0;j<lon_row.size();j++){
                    field_out[i][j] = field_in[lat_row.get(i)][lon_row.get(j)];
                }
            }
        }

        return new grid_extreme_coordsResults(lat_red,lon_red,field_out);
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
        int totalLines = 0;//Total file lines (starting from 3rd line)
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
            //int totalLines = 0;//Total file lines (starting from 3rd line)
            while(coastFile.hasNext()){//while EOF
                if(varType%2==0){//1st element
                    lon_tmp.add(Double.parseDouble(coastFile.next()));//adding in temp longitude array
                }
                else{//2nd element
                    lat_tmp.add(Double.parseDouble(coastFile.next()));//adding in temp latitude array
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
            System.out.println("Lines readed: "+totalLines);
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("readout_coast: "+e.getMessage());
            debug.CloseFile();
            e.printStackTrace();
        }
    }

    private readout_bathyResults readout_bathy(int bathy_code){
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
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("readout_bathy: unknow bathy DB code!");
                debug.CloseFile();
                break;
        }
        //Parsing file:
        MyNetCDFParser test = new MyNetCDFParser(filename);
        parseMedOneMinResults out = test.parseMedOneMin();
        if(out == null){
            System.out.println("Parsing fail!");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("readout_bathy: Parsing fail!");
            debug.CloseFile();
            return null;
        }
        ArrayList<Double> latTmp = out.getLat();
        ArrayList<Double> lonTmp = out.getLon();
        double[][] zTmp = out.getDepth();
        //if sea depth <=0, setting as NaN
        for(int i =0 ;i<zTmp.length;i++){
            for(int j=0;j<zTmp[0].length;j++){
                if(zTmp[i][j]<=0){
                    zTmp[i][j]= Double.NaN;
                }
            }
        }
        //reduction within Med. Sea basin:
        double lon_min= -6.5;
        double lon_max= 37.0;
        double lat_min= 30.0;
        double lat_max= 46.0;
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
        double[][] z = new double[latOkIndexes.size()][lonOkIndexes.size()];
        for(int i=0;i<latOkIndexes.size();i++) {
            lat.add(latTmp.get(latOkIndexes.get(i)));
        }
        for(int i=0;i<lonOkIndexes.size();i++){
            lon.add(lonTmp.get(lonOkIndexes.get(i)));
        }
        for(int i=0;i<latOkIndexes.size();i++){
            for(int j=0;j<lonOkIndexes.size();j++){
                z[i][j] = zTmp[latOkIndexes.get(i)][lonOkIndexes.get(j)];
            }
        }
        return new readout_bathyResults(lat, lon, z);
    }

    private nodes_free_form_barrierResults nodes_free_from_barrier(int[] nodes){
        // % finds nodes not on the landmass (both continent and islands)
        // %
        // % removal of islands outside of bbox:
        // % computation of nodes not on the landmasses:
        int Nn = nodes.length;

        //removal of islands outside of bbox:
        int n_prima = this.y_islands.size();

        boolean[] bool_in = new boolean[this.y_islands.size()];
        for(int i=0;i<this.y_islands.size();i++){
            if((this.sGrid.getDB_xi()<= this.x_islands.get(i)) && (this.x_islands.get(i)<=this.sGrid.getDB_xf())
                    && (this.sGrid.getDB_yi()<=this.y_islands.get(i) && (this.y_islands.get(i)<=this.sGrid.getDB_yf())) ){
                bool_in[i]=true;
            }else{
                bool_in[i]=false;
            }
        }
        boolean[] bool_nan = new boolean[this.y_islands.size()];
        for(int i=0;i<this.y_islands.size();i++){
            if(Double.isNaN(this.x_islands.get(i)) && Double.isNaN(this.y_islands.get(i))){
                bool_nan[i] = true;
            } else {
                bool_nan[i] = false;
            }
        }
        ArrayList<Integer> isl_in_bb = new ArrayList<>();
        for(int i=0;i<this.y_islands.size();i++){
            if(bool_in[i] || bool_nan[i]){
                isl_in_bb.add(i);
            }
        }
        ArrayList<Double> x_islands_tmp = new ArrayList<>();
        ArrayList<Double> y_islands_tmp = new ArrayList<>();
        for(Integer element : isl_in_bb){
            x_islands_tmp.add(this.x_islands.get(element));
            y_islands_tmp.add(this.y_islands.get(element));
        }
        this.x_islands = x_islands_tmp;
        this.y_islands = y_islands_tmp;

        int n_dopo = this.y_islands.size();//includes NaNs separators between islands (i.e., 5608 elements)
        int red_fact = Math.round((1 - n_dopo/n_prima)*100);

        ArrayList<Double> y_coast = new ArrayList<>();
        ArrayList<Double> x_coast = new ArrayList<>();
        for(int i=0; i<this.y_islands.size();i++){
            y_coast.add(this.y_islands.get(i));
            x_coast.add(this.x_islands.get(i));
        }
        y_coast.add(Double.NaN);
        x_coast.add(Double.NaN);

        for(int i=0;i<this.y_continent.size();i++){
            y_coast.add(this.y_continent.get(i));
            x_coast.add(this.x_continent.get(i));
        }

        //--------------------------------------------------------------------------------------------
        //computation of nodes not on the landmasses:
        ArrayList<Integer> free_nodes = new ArrayList<>();
        MyFileWriter fid = new MyFileWriter(true,"freeNodes_DB.dat.permil");
        for(int ie=0;ie<Nn;ie++){
            int frac_done = new Double(Math.floor(1000*ie/Nn)).intValue();
            fid.WriteLine(""+frac_done);

            double xP = this.xy[ie][0];
            double yP = this.xy[ie][1];

            inpolygonResults tmp = Utility.inpolygon(xP, yP, this.x_islands, this.y_islands);
            boolean IN_i = tmp.getIn();
            boolean ON_i = tmp.getOn();
            tmp = Utility.inpolygon(xP, yP, this.x_continent, this.y_continent);
            boolean IN_c = tmp.getIn();
            boolean ON_c = tmp.getOn();

            boolean canAdd = (IN_i == false) && (ON_i == false) && (IN_c == true) && (ON_c == false);

            if(canAdd){
                free_nodes.add(nodes[ie]);
            }
        }
        if(free_nodes.size() == 0){
            System.out.println("no free nodes found in graph");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("nodes_free_from_barrier: no free nodes found in graph");
            debug.CloseFile();
            System.exit(0);
        }
        fid.CloseFile();
        return new nodes_free_form_barrierResults(free_nodes.size(),free_nodes,red_fact);
    }

    public Fields_regriddingResults Fields_regridding(){
        this.logFile = new MyFileWriter("","",true);
        this.logFile.WriteLog("processing environmental fields...");
        this.logFile.CloseFile();
        //Spatial and temporal Grid:
        //Time parameters preprocessing:
        // % Fields reading:
        // % Spatial subsetting:
        // % Time processing:
        // % seaOverLand:

        //Time parameters preprocessing:
        double wave_t1= 12.5; // wave file start time is 1230 UTC (WW3 4e)
        double wind_t1 = Double.NaN;
        if (this.optim.getWindModel()==11 || this.optim.getWindModel()==12){//EXMWF
            wind_t1 = 12.0;//ECMWF wind file (analysis) start time is 1200 UTC
        } else if (this.optim.getWindModel() == 2){
            wind_t1 = 15.0; //COSMO-ME wind file (forecast) start time is 1500 UTC  *** use also analysis in the future!
        } else {
            System.out.println("unknown wind model code!");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("Fields_regridding: unknown wind model code!");
            debug.CloseFile();
            System.exit(0);
        }
        String dateListFile = "inputFiles/fields/an_dates_DB.txt";
        MyTxtParser file = new MyTxtParser(dateListFile);
        this.tGrid.setLatest_date(file.tail(1).get(0));

        String l_date_str = this.tGrid.getLatest_date()+"1200";
        long l_num = Utility.datenum(l_date_str,"yyyyMMddHHmm"); // taken at 1200 UTC (analysis time)

        this.tGrid.setDepDateTime(Utility.datenum(2000+this.dep_datetime.getYear(),this.dep_datetime.getMonth(),
                this.dep_datetime.getDay(),this.dep_datetime.getHour(),this.dep_datetime.getMin()));
        //Adding 6 hours delay to original departure date time:
        double delayed = this.tGrid.getDepDateTime()-0.25;
        //hrs elapsed between latest analysis and departure time
        int cento = 100;
        long deltaHr_anls=0;
        if(this.forcing.getAnalytic()==1){
            deltaHr_anls = 1;
        } else{
            deltaHr_anls = Math.round(cento*this.constants.getTwentyfour()*(delayed - l_num)/cento);
        }

        this.tGrid.setWave_dep_TS(Math.round(1+ deltaHr_anls - (wave_t1 - this.constants.getTwelve())));

        this.tGrid.setWind_dep_TS(Math.round(1+ deltaHr_anls - (wave_t1 - this.constants.getTwelve())));

        this.logFile = new MyFileWriter("","",true);
        this.logFile.WriteLog("\tlatest analysis date and time :"+Utility.datestr(l_num));
        this.logFile.WriteLog("\tdeparture date and time :"+Utility.datestr(this.tGrid.getDepDateTime()));
        this.logFile.CloseFile();

        //-----------------------------------------------------------------------------------------------------
        //Fields reading:
        if(this.forcing.getAnalytic()==1){
            this.tGrid.setNt(this.tGrid.getMinNt());
            prepare_sqrtY_fieldResults out = prepare_sqrtY_field();
            double[][][] VTDH_out = out.getVTDH_b();
            double[][][] VDIR_out = out.getVDIR_b();
            double[][][] VTPK_out = Utility.NaN3Dmatrix(VTDH_out.length,VTDH_out[0].length,VTDH_out[0][0].length);
            double[][][] windMAGN_out = Utility.NaN3Dmatrix(VTDH_out.length,VTDH_out[0].length,VTDH_out[0][0].length);
            double[][][] windDIR_out = Utility.NaN3Dmatrix(VTDH_out.length,VTDH_out[0].length,VTDH_out[0][0].length);
            int[] time_steps = null;
            return new Fields_regriddingResults(time_steps, VTDH_out, VTPK_out, VDIR_out, windMAGN_out, windDIR_out);
        } else { //Reading of environmental Fields:
            readout_envFieldsResults envFieldsResults = readout_envFields();
            //this check must be run after input file time size has been determined
            check_start_timestep(deltaHr_anls);

            //-----------------------------------------------------------------------------------------------------
            //Spatial subsetting:

            //Wave data reduction to inset grid:
            mFields_reductionResults waveDataReduction =mFields_reduction(envFieldsResults.getLat_wave(),envFieldsResults.getLon_wave(), envFieldsResults.getVTDH(), envFieldsResults.getVTPK(), envFieldsResults.getVDIR());
            double[] lat_wave_Inset = waveDataReduction.getLat_red();
            double[] lon_wave_Inset = waveDataReduction.getLon_red();
            double[][][] VTDH_Inset = waveDataReduction.getOut1();
            double[][][] VTPK_Inset = waveDataReduction.getOut2();
            double[][][] VDIR_Inset = waveDataReduction.getOut3();
            meshgridResults meshGWave = Utility.meshgrid(lon_wave_Inset, lat_wave_Inset);
            double[][] lon_wave_m = meshGWave.getX();
            double[][] lat_wave_m = meshGWave.getY();

            //Wave direction into Cartesian components:
            double[][][] X_Inset = new double[VDIR_Inset.length][VDIR_Inset[0].length][VDIR_Inset[0][0].length];
            double[][][] Y_Inset = new double[VDIR_Inset.length][VDIR_Inset[0].length][VDIR_Inset[0][0].length];
            for(int i=0;i<VDIR_Inset.length; i++){
                for(int j=0;j<VDIR_Inset[0].length; j++){
                    for(int k=0;k<VDIR_Inset[0][0].length; k++){
                        X_Inset[i][j][k] = Math.sin(VDIR_Inset[i][j][k]);
                        Y_Inset[i][j][k] = Math.cos(VDIR_Inset[i][j][k]);
                    }
                }
            }
            //clear VDIR_Inset
            //Wind  data reduction to inset grid:
            mFields_reductionResults ecmwfDataReduction = mFields_reduction(envFieldsResults.getEcmwf_lat_wind(), envFieldsResults.getEcmwf_lon_wind(), envFieldsResults.getEcmwf_U10m(), envFieldsResults.getEcmwf_V10m(), null);
            double[] ecmwf_lat_wind = ecmwfDataReduction.getLat_red();
            double[] ecmwf_lon_wind = ecmwfDataReduction.getLon_red();
            double[][][] ecmwf_U10M_Inset = ecmwfDataReduction.getOut1();
            double[][][] ecmwf_V10M_Inset = ecmwfDataReduction.getOut2();
            mFields_reductionResults cosmoDataReduction = mFields_reduction(envFieldsResults.getCosmo_lat_wind(), envFieldsResults.getCosmo_lon_wind(), envFieldsResults.getCosmo_U10m(), envFieldsResults.getCosmo_V10m(), null);
            double[] cosmof_lat_wind = cosmoDataReduction.getLat_red();
            double[] cosmo_lon_wind = cosmoDataReduction.getLon_red();
            double[][][] cosmo_U10M_Inset = cosmoDataReduction.getOut1();
            double[][][] cosmo_V10M_Inset = cosmoDataReduction.getOut2();

            //Wind unit conversion
            for(int i=0;i<ecmwf_U10M_Inset.length;i++){
                for(int j=0;j<ecmwf_U10M_Inset[0].length;j++){
                    for(int k=0;k<ecmwf_U10M_Inset[0][0].length;k++){
                        ecmwf_U10M_Inset[i][j][k] = this.constants.getMs2kts()*ecmwf_U10M_Inset[i][j][k];
                    }
                }
            }
            for(int i=0;i<ecmwf_V10M_Inset.length;i++){
                for(int j=0;j<ecmwf_V10M_Inset[0].length;j++){
                    for(int k=0;k<ecmwf_V10M_Inset[0][0].length;k++){
                        ecmwf_V10M_Inset[i][j][k] = this.constants.getMs2kts()*ecmwf_V10M_Inset[i][j][k];
                    }
                }
            }
            for(int i=0;i<cosmo_U10M_Inset.length;i++){
                for(int j=0;j<cosmo_U10M_Inset[0].length;j++){
                    for(int k=0;k<cosmo_U10M_Inset[0][0].length;k++){
                        cosmo_U10M_Inset[i][j][k] = this.constants.getMs2kts()*cosmo_U10M_Inset[i][j][k];
                    }
                }
            }
            for(int i=0;i<cosmo_V10M_Inset.length;i++){
                for(int j=0;j<cosmo_V10M_Inset[0].length;j++){
                    for(int k=0;k<cosmo_V10M_Inset[0][0].length;k++){
                        cosmo_V10M_Inset[i][j][k] = this.constants.getMs2kts()*cosmo_V10M_Inset[i][j][k];
                    }
                }
            }

            //-----------------------------------------------------------------------------------------------------
            //Time processing:
            fieldStatsResults field_res = fieldStats(this.bathy_Inset, VTDH_Inset, VTPK_Inset, ecmwf_U10M_Inset, ecmwf_V10M_Inset, cosmo_U10M_Inset, cosmo_V10M_Inset);
            double ecmwf_dir_avg = field_res.getEcmwf_dir_avg();
            double ecmwf_dir_std = field_res.getEcmwf_dir_std();
            double cosmo_dir_avg = field_res.getCosmo_dir_avg();
            double cosmo_dir_std = field_res.getCosmo_dir_std();
            //Just for GMD paper's route #2 (1589_c) uncomment following line:
            //this.fstats.setWlenght_max(85);

            this.estim_Nt(estGdtDist, Math.max(this.tGrid.getWave_dep_TS(), this.tGrid.getWind_dep_TS()), this.H_array_m, this.vel_LUT);

            int delta_hr_dep_int = Math.round(deltaHr_anls);
            int[] interpTimes = new int[(int)this.tGrid.getNt()];
            for(int i=0;i<interpTimes.length;i++){
                interpTimes[i] = delta_hr_dep_int+(i-1);
            }
            int[] time_steps = new int[(int) this.tGrid.getMaxNt()];
            //time_steps always counted starting from step #1:
            for(int i=0;i<time_steps.length; i++){
                time_steps[i] = i+1;
            }

            double[] wind_origTimes = new double[0];
            double[][][] U10M_Inset = new double[0][0][0];
            double[][][] V10M_Inset = new double[0][0][0];//FIN QUI OK
            if(this.optim.getWindModel() == 11 || this.optim.getWindModel() == 12){//ecmwf
                wind_origTimes = envFieldsResults.getEcmwf_wind_origTimes();
                U10M_Inset = ecmwf_U10M_Inset;
                V10M_Inset = ecmwf_V10M_Inset;
            } else if(this.optim.getWindModel() == 2){//cosmo-me
                wind_origTimes = envFieldsResults.getCosmo_wind_origTimes();
                U10M_Inset = cosmo_U10M_Inset;
                V10M_Inset = cosmo_V10M_Inset;
            }

            double[][][] U10M_at_TS = Utility.interp1(wind_origTimes,U10M_Inset, interpTimes);
            double[][][] V10M_at_TS = Utility.interp1(wind_origTimes,V10M_Inset, interpTimes);
            //(wave_t1-const.twelve)

            double[][][] VTDH_times = new double[VTDH_Inset.length][interpTimes.length][VTDH_Inset[0][0].length];
            for(int i=0;i<VTDH_Inset.length;i++){
                for(int j=0;j<interpTimes.length;j++){
                    for(int k=0; k<VTDH_Inset[0][0].length; k++){
                        VTDH_times[i][j][k]=VTDH_Inset[i][interpTimes[j]][k];
                    }
                }
            }

            double[][][] VTPK_times = new double[VTPK_Inset.length][interpTimes.length][VTPK_Inset[0][0].length];
            for(int i=0;i<VTPK_Inset.length;i++){
                for(int j=0;j<interpTimes.length;j++){
                    for(int k=0; k<VTPK_Inset[0][0].length; k++){
                        VTPK_times[i][j][k]=VTPK_Inset[i][interpTimes[j]][k];
                    }
                }
            }


            double[][][] X_times = new double[X_Inset.length][interpTimes.length][X_Inset[0][0].length];
            for(int i=0;i<X_Inset.length;i++){
                for(int j=0;j<interpTimes.length;j++){
                    for(int k=0; k<X_Inset[0][0].length; k++){
                        X_times[i][j][k]=X_Inset[i][interpTimes[j]][k];
                    }
                }
            }


            double[][][] Y_times = new double[Y_Inset.length][interpTimes.length][Y_Inset[0][0].length];
            for(int i=0;i<Y_Inset.length;i++){
                for(int j=0;j<interpTimes.length;j++){
                    for(int k=0; k<Y_Inset[0][0].length; k++){
                        Y_times[i][j][k]=Y_Inset[i][interpTimes[j]][k];
                    }
                }
            }


            //-----------------------------------------------------------------------------------------------------
            // seaOverLand:
            this.logFile = new MyFileWriter("","",true);//DA DEBUGGARE DA QUI IN POI!
            this.logFile.WriteLog("\tseaOverLand extrapolation...");
            this.logFile.CloseFile();

            //Wave fields processing:
            ArrayList<double[][][]> seaOverLand_3stepsOut = seaOverLand_3steps(this.lon_bathy_Inset, this.lat_bathy_Inset, this.lsm_mask, lon_wave_m, lat_wave_m,VTDH_times, VTPK_times, X_times, Y_times);
            double[][][] VTDH_out = seaOverLand_3stepsOut.get(0);
            double[][][] VTPK_out = seaOverLand_3stepsOut.get(1);
            X_times = seaOverLand_3stepsOut.get(2);
            Y_times = seaOverLand_3stepsOut.get(3);

            double[][][] VDIR_times = new double[X_times.length][X_times[0].length][X_times[0][0].length];
            for(int i=0;i<VDIR_times.length;i++){
                for(int j=0;j<VDIR_times[0].length;j++){
                    for(int k=0;k<VDIR_times[0][0].length;k++){
                        VDIR_times[i][j][k] = Math.atan2(Y_times[i][j][k], X_times[i][j][k]);
                    }
                }
            }

            double[][][] VDIR_out = changeDirRule(VDIR_times);
            for(int i=0;i<VDIR_out.length;i++){
                for(int j=0;j<VDIR_out[0].length;j++){
                    for(int k=0;k<VDIR_out[0][0].length;k++){
                        VDIR_out[i][j][k]=VDIR_out[i][j][k]/this.constants.getDeg2rad();
                    }
                }
            }
            double[] lat_wind_Inset = new double[0];
            double[] lon_wind_Inset = new double[0];
            if (this.optim.getWindModel() == 11 || this.optim.getWindModel() == 12) {//ecmwf
                lat_wind_Inset = ecmwfDataReduction.getLat_red();
                lon_wind_Inset = ecmwfDataReduction.getLon_red();
            } else {
                if(this.optim.getWindModel() == 2){ //cosmo-me
                    lat_wind_Inset = cosmoDataReduction.getLat_red();
                    lon_wind_Inset = cosmoDataReduction.getLon_red();
                }
            }

            meshgridResults wind_m = Utility.meshgrid(lon_wind_Inset, lat_wind_Inset);
            double[][] lon_wind_m = wind_m.getX();
            double[][] lat_wind_m = wind_m.getY();

            ArrayList<double[][][]> seaOut = seaOverLand_3steps(lon_bathy_Inset,lat_bathy_Inset,lsm_mask,lon_wind_m,lat_wind_m, U10M_at_TS,V10M_at_TS);
            U10M_at_TS = seaOut.get(0);
            V10M_at_TS = seaOut.get(1);

            double[][][] windDir_times = new double[V10M_at_TS.length][V10M_at_TS[0].length][V10M_at_TS[0][0].length];
            double[][][] windDIR_out = new double[V10M_at_TS.length][V10M_at_TS[0].length][V10M_at_TS[0][0].length];
            double[][][] windMAGN_out =new double[V10M_at_TS.length][V10M_at_TS[0].length][V10M_at_TS[0][0].length];
            for(int i=0;i<windDir_times.length;i++){
                for(int j=0;j<windDir_times[0].length;j++){
                    for(int k=0;k<windDir_times[0][0].length;k++){
                        windDir_times[i][j][k] = Math.atan2(V10M_at_TS[i][j][k],U10M_at_TS[i][j][k]);
                    }
                }
            }
            windDIR_out = changeDirRule(windDir_times);
            for(int i=0;i<windDIR_out.length;i++){
                for(int j=0;j<windDIR_out[0].length;j++){
                    for(int k=0;k<windDIR_out[0][0].length;k++){
                        windDIR_out[i][j][k] = windDIR_out[i][j][k]/this.constants.getDeg2rad();
                    }
                }
            }

            for(int i=0;i<windMAGN_out.length;i++){
                for(int j=0;j<windMAGN_out[0].length;j++){
                    for(int k=0;k<windMAGN_out[0][0].length;k++){
                        windMAGN_out[i][j][k] = Math.sqrt(Math.pow(U10M_at_TS[i][j][k],2)+Math.pow(V10M_at_TS[i][j][k],2));
                    }
                }
            }
            return new Fields_regriddingResults(time_steps, VTDH_out, VTPK_out, VDIR_out, windMAGN_out, windDIR_out);
        }
    }

    private ArrayList<double[][][]> seaOverLand_3steps(ArrayList<Double> lon_bathy, ArrayList<Double> lat_bathy, double[][] lsm_mask, double[][] lon_f, double[][] lat_f, double[][][]... varargs){
        //% "sea over land" 3-step process:
        //% (1) extrapolation:
        //% (2) regridding to bathy-grid:
        //% (3) masking landmass on target grid:
        //%
        //%--------------------------------------------------------------------------
        if(varargs.length < 1 || varargs.length > 4){
            System.out.println("seaOverLand_3steps: varargs must be between 1 and 4");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("\tseaOverLand_3steps: varargs must be between 1 and 4");
            debug.CloseFile();
            System.exit(0);
        }

        if(Math.min(lon_f.length, lon_f[0].length) < this.sGrid.getMinNoGridPoints()){
            System.out.println("seaOverLand_3steps: too small or too narrow bounding box");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("\tseaOverLand_3steps: too small or too narrow bounding box");
            debug.CloseFile();
            System.exit(0);
        }

        int n_loops = 50;
//        double[] lon_fT = new double[lon_f[0].length];
//        for(int i=0;i<lon_f[0].length;i++){
//            lon_fT[i] = lon_f[0][i];
//        }
//        double[] lat_fT = new double[lat_f.length];
//        for(int i=0;i<lat_f.length;i++){
//            lat_fT[i]=lat_f[i][0];
//        }
        double[][][] myfield_bathy = new double[lon_bathy.size()][(int)this.tGrid.getNt()][lat_bathy.size()];
        ArrayList<double[][][]> out = new ArrayList<>();
        for(double[][][] myfield : varargs){
            //(1) extrapolation - % #GM: check also mdata_EWeights.m:
            for(int it=0;it<(int)this.tGrid.getNt(); it++){
                double[][] myfield_mat = new double[0][0];
                double[][][] tmp = new double[myfield.length][1][myfield[0][0].length];
                for(int z=0;z<myfield.length;z++){
                    for(int j =0; j<myfield[0][0].length; j++){
                        tmp[z][0][j] = myfield[z][it][j];
                    }
                }
                if(this.forcing.getWind()!=1){
                    myfield_mat = SeaOverLand(Utility.squeeze(tmp),n_loops);//DEBUG DA QUI
                } else {
                    myfield_mat = Utility.squeeze(tmp);
                }

                // (2) regridding to bathy-grid:
                //myfield_mat= squeeze(myfield_mat); ????

                //check that bbox is larger than minimum allowed area:
                if(it==0){
                    int size_min = Math.min(myfield_mat.length, myfield_mat[0].length);
                    int nan_rank = Utility.rank(Utility.convertDouble(Utility.isnan(myfield_mat)));
                    if((size_min-nan_rank) < 2){
                        System.out.println("seaOverLand_3steps: too few sea grid points");
                        MyFileWriter debug = new MyFileWriter("","debug",false);
                        debug.WriteLog("\tseaOverLand_3steps: too few sea grid points");
                        debug.CloseFile();
                    }
                }
                double[] lon_bathyArray = new double[lon_bathy.size()];
                for(int i=0;i<lon_bathy.size();i++){
                    lon_bathyArray[i] = lon_bathy.get(i);
                }
                double[] lat_bathyArray = new double[lat_bathy.size()];
                for(int i=0;i<lat_bathy.size();i++){
                    lat_bathyArray[i] = lat_bathy.get(i);
                }
                double[][] tmpMtx=Utility.interp2(lon_f, lat_f, myfield_mat, lat_bathyArray, lon_bathyArray);
//                double[][] tmpMtx = Utility.interp2(lon_fT, lat_fT, myfield_mat,lat_bathyArray, lon_bathyArray);
                for(int i=0;i<tmpMtx.length;i++){
                    for(int j=0;j<tmpMtx[0].length;j++){
                        myfield_bathy[i][it][j]=tmpMtx[i][j];
                    }
                }
                //FORSE??
//                for(int i=0;i<tmpMtx.length;i++){
//                    for(int j=0;j<tmpMtx[0].length;j++){
//                        myfield[i][it][j]=tmpMtx[i][j];
//                    }
//                }
            }
            double[][][] myfield_Inset = new double[myfield_bathy.length][myfield_bathy[0].length][myfield_bathy[0][0].length];
            if(this.visualization.getScientific_mode()==0){
                /*
                * % (3) masking landmass on target grid:
                %
                % ok also for wind 10mt.
                % warning: for ocean currents it will become necessary to set field=0 on the
                % landmass for reproducing the no-slip boundary condition!*/
                if(this.forcing.getWind()!=1){
                    for(int it=0;it<(int) this.tGrid.getNt(); it++){
                        double[][][] tmp = new double[myfield_bathy.length][1][myfield_bathy[0][0].length];
                        for(int row=0;row<myfield_bathy.length;row++){
                            for(int col=0;col<myfield_bathy[0][0].length; col++){
                                tmp[row][0][col] = myfield_bathy[row][0][col];
                            }
                        }
                        double[][] vv1 = Utility.MatrixComponentXcomponent(Utility.squeeze(tmp),Utility.transposeMatrix(this.lsm_mask));
                        for(int i=0;i<myfield_bathy.length;i++){
                            for(int j=0;j<myfield_bathy[0][0].length;j++){
                                myfield_bathy[i][it][j] = vv1[i][j];
                            }
                        }
                    }
                } else {
                    myfield_Inset = myfield_bathy;
                }
            } else{
                myfield_Inset = myfield_bathy;
            }
            //varargout{ia} = myfield_Inset;
            out.add(myfield_Inset);
        }
        return out;
    }

    private static double[][] SeaOverLand(double[][] matrice_in, int loop){
        //source: Nicoletta Fabroni (SINCEM, Ravenna), developed for Relocatable HOPS package
        //% versione ricevuta da I.Federico via S. Falchetti in data 24/6/2013

        double[][] dummy = new double[matrice_in.length+2][matrice_in[0].length+2];
        for(int i=0;i<dummy.length;i++){
            for(int j=0;j<dummy[0].length;j++){
                if(i==0 || j==0 || i==(dummy.length-1) || j==(dummy[0].length-1)){
                    dummy[i][j] = 9999;
                }else{
                    dummy[i][j] = matrice_in[i-1][j-1];
                }
            }
        }

        ArrayList<Integer> idx = Utility.find(dummy, Double.NaN);

        for(int i=0;i<dummy.length;i++){
            for(int j=0;j<dummy[0].length;j++){
                if(i==0 || j==0 || i==(dummy.length-1) || j==(dummy[0].length-1)){
                    dummy[i][j] = Double.NaN;
                }
            }
        }
        //matrice dei vicini che voglio analizzare
        int M = dummy.length;
        double[] neighbor_offsets = new double[8];
        neighbor_offsets[0] = M;
        neighbor_offsets[1] = M+1;
        neighbor_offsets[2] = 1;
        neighbor_offsets[3] = -M+1;
        neighbor_offsets[4] = -M;
        neighbor_offsets[5] = -M-1;
        neighbor_offsets[6] = -1;
        neighbor_offsets[7] = M-1;
        int count=0;
        while(count<loop && idx.size()>0){

            //Creo matrice dove ogni elemento è la somma di idx con neighbor_offset
            double[][] neighbors = new double[neighbor_offsets.length][idx.size()];
            for(int i=0;i<neighbors.length;i++){
                for(int j=0;j<neighbors[0].length;j++){
                    neighbors[i][j] = idx.get(j)+neighbor_offsets[i];
                }
            }

            double[][] mat = new double[neighbors.length][neighbors[0].length];
            for(int i=0;i<mat.length;i++){
                for(int j=0;j<mat[0].length;j++){
                    double currentElem = neighbors[i][j];
                    int _colIndex = ((int) Math.floor(currentElem/dummy.length));
                    int _rowIndex = ((int) currentElem%dummy.length);
                    mat[i][j] = dummy[_rowIndex][_colIndex];
                }
            }

            boolean[][] nans = Utility.isnan(mat);
            int[] snn = new int[nans[0].length];
            for(int i=0;i<snn.length;i++)
                snn[i]=-1;

            for(int i=0;i<nans.length;i++){
                for(int j=0;j<nans[0].length;j++){
                    if(nans[i][j]){
                        mat[i][j] = 0;
                    } else{
                        snn[j]++;
                    }
                }
            }

            double[] media = Utility.sum(mat);
            for(int i=0;i<snn.length;i++){
                media[i]= media[i]/(snn[i]+1);
            }

            for(int i=0;i<idx.size();i++){
                double currentElem = media[i];
                int currentIndex = idx.get(i);
                int _colIndex = ((int) Math.floor(currentIndex/dummy.length));
                int _rowIndex = (currentIndex%dummy.length);
                dummy[_rowIndex][_colIndex] = currentElem;
            }

            int i=0;
            ArrayList<Integer> tmp=new ArrayList<>();
            while(i<idx.size() && i<snn.length){
                if(snn[i]==-1)
                    tmp.add(idx.get(i));
                i++;
            }
            idx = tmp;

            count++;
        }

        double[][] matrice_out = new double[dummy.length-2][dummy[0].length-2];
        for(int i=0;i<matrice_out.length;i++){
            for(int j=0;j<matrice_out[0].length;j++){
                matrice_out[i][j] = dummy[i+1][j+1];
            }
        }
        return matrice_out;
    }



    private void estim_Nt(double estGdtDist, long start_timestep, ArrayList<Double> H_array_m, double[][] ship_v_LUT){
        // % estimates Tgrid.Nt (number of time steps)
        // % - key quantity for RAM allocation and perfomance
        // % ( "temporal bbox" )
        // %
        // %
        // % fract parameter has a great impact on CPU times and is defined as follows:
        // % 0[more conservative]:  v_kts_signif= v_kts_min;
        // % 1[very optimistic]  :  v_kts_signif= v_kts_max

        double fract = 0.2; //RAM saving setting
        double v_kts_min = Double.NaN;
        double v_kts_max = Double.NaN;
        if(this.ship.getVessType() != this.ship.getSailType()){
            double[] ship_v_LUT_1stCol = new double[ship_v_LUT.length];
            for(int i=0;i<ship_v_LUT.length;i++){
                ship_v_LUT_1stCol[i] = ship_v_LUT[i][0];
            }
            v_kts_min = Utility.interp1(H_array_m, ship_v_LUT_1stCol, this.fstats.getWheight_max());
            v_kts_max = Utility.interp1(H_array_m, ship_v_LUT_1stCol, this.fstats.getWheight_min());
        } else{//sailboat
            v_kts_min=Utility.min(ship_v_LUT);
            v_kts_max=Utility.max(ship_v_LUT);
        }

        double v_kts_signif = v_kts_min*(1-fract) + v_kts_max*fract;

        double Nt_1 = this.tGrid.getMaxNt() - start_timestep + 1;
        double Nt_2b = Math.ceil(estGdtDist/v_kts_signif);

        double Nt_2 = Math.max(this.tGrid.getMinNt(), Nt_2b);
        this.tGrid.setNt(Math.min(Nt_1, Nt_2));

        System.out.println("# time steps of forecast file employed: "+this.tGrid.getNt());
        this.logFile = new MyFileWriter("","",true);
        this.logFile.WriteLog("\t# time steps of forecast file employed:"+this.tGrid.getNt());
        this.logFile.CloseFile();
    }

    private fieldStatsResults fieldStats(double[][] bathy, double[][][] VTDH, double[][][] VTPK, double[][][] ecmwf_U10M, double[][][] ecmwf_V10M, double[][][] cosmo_U10M, double[][][] cosmo_V10M){
        //fields statistics:
        // Please note that computed max and min values
        // refer to the whole temporal window selected via timeprocess_* routines

        double piccolo = 0.1;
        //----------------------------------------------------------
        //Bathy
        this.fstats.setBathy_min(Utility.minNotNaN(bathy));
        this.fstats.setBathy_max(Utility.maxNotNaN(bathy));

        //----------------------------------------------------------
        //Wave
        this.fstats.setWheight_min((1.0-piccolo)* Utility.min3d(VTDH));
        this.fstats.setWheight_max((1.0+piccolo)* Utility.max3d(VTDH));
        this.fstats.setWheight_avg(nanmean2(VTDH));

        this.fstats.setWperiod_min((1-piccolo)* Utility.min3d(VTPK));
        this.fstats.setWperiod_max((1+piccolo)*Utility.max3d(VTPK));

        double[][][] Lambda = wave_dispersion(VTPK);
        this.fstats.setWlength_min((1-piccolo)*Utility.min3d(Lambda));
        this.fstats.setWlenght_max(Utility.max3d(Lambda));
        //----------------------------------------------------------
        //Wind
        double[][][] ecmwf_wind10M_magn = new double[ecmwf_U10M.length][ecmwf_U10M[0].length][ecmwf_U10M[0][0].length];
        for(int i=0;i<ecmwf_wind10M_magn.length; i++){
            for(int j=0;j<ecmwf_wind10M_magn[0].length; j++){
                for(int k=0;k<ecmwf_wind10M_magn[0][0].length; k++){
                    ecmwf_wind10M_magn[i][j][k] = Math.sqrt(Math.pow(ecmwf_U10M[i][j][k],2) + Math.pow(ecmwf_V10M[i][j][k],2));
                }
            }
        }
        double ecmwf_wind10M_max = (1+piccolo) * Utility.max3d(ecmwf_wind10M_magn);
        double ecmwf_wind10M_min = (1-piccolo) * Utility.min3d(ecmwf_wind10M_magn);

        double[][][] cosmo_wind10M_magn = new double[cosmo_U10M.length][cosmo_U10M[0].length][cosmo_U10M[0][0].length];
        for(int i=0;i<cosmo_wind10M_magn.length; i++){
            for(int j=0;j<cosmo_wind10M_magn[0].length; j++){
                for(int k=0;k<cosmo_wind10M_magn[0][0].length; k++){
                    cosmo_wind10M_magn[i][j][k] = Math.sqrt(Math.pow(cosmo_U10M[i][j][k],2) + Math.pow(cosmo_V10M[i][j][k],2));
                }
            }
        }
        double cosmo_wind10M_max = (1+piccolo) * Utility.max3d(cosmo_wind10M_magn);
        double cosmo_wind10M_min = (1-piccolo) * Utility.min3d(cosmo_wind10M_magn);

        //wind direction:
        double[][][] TmpEcmwf_cos_avg = new double[ecmwf_U10M.length][ecmwf_U10M[0].length][ecmwf_U10M[0][0].length];
        double[][][] TmpEcmwf_sin_avg = new double[ecmwf_V10M.length][ecmwf_V10M[0].length][ecmwf_V10M[0][0].length];
        for(int i=0;i<TmpEcmwf_cos_avg.length; i++){
            for(int j=0;j<TmpEcmwf_cos_avg[0].length; j++){
                for(int k=0;k<TmpEcmwf_cos_avg[0][0].length; k++){
                    TmpEcmwf_cos_avg[i][j][k] = ecmwf_U10M[i][j][k]/ecmwf_wind10M_magn[i][j][k];
                }
            }
        }
        for(int i=0;i<TmpEcmwf_sin_avg.length; i++){
            for(int j=0;j<TmpEcmwf_sin_avg[0].length; j++){
                for(int k=0;k<TmpEcmwf_sin_avg[0][0].length; k++){
                    TmpEcmwf_sin_avg[i][j][k] = ecmwf_V10M[i][j][k]/ecmwf_wind10M_magn[i][j][k];
                }
            }
        }
        double ecmwf_cos_avg = nanmean2(TmpEcmwf_cos_avg);
        double ecmwf_sin_avg = nanmean2(TmpEcmwf_sin_avg);

        double ecmwf_dir_avg = Math.atan2(ecmwf_sin_avg, ecmwf_cos_avg);
        ecmwf_dir_avg = changeDirRule(ecmwf_dir_avg)/this.constants.getDeg2rad();


        double[][][] TmpCosmo_cos_avg = new double[cosmo_U10M.length][cosmo_U10M[0].length][cosmo_U10M[0][0].length];
        double[][][] TmpCosmo_sin_avg = new double[cosmo_V10M.length][cosmo_V10M[0].length][cosmo_V10M[0][0].length];
        for(int i=0;i<TmpCosmo_cos_avg.length; i++){
            for(int j=0;j<TmpCosmo_cos_avg[0].length; j++){
                for(int k=0;k<TmpCosmo_cos_avg[0][0].length; k++){
                    TmpCosmo_cos_avg[i][j][k] = cosmo_U10M[i][j][k]/cosmo_wind10M_magn[i][j][k];
                }
            }
        }
        for(int i=0;i<TmpCosmo_sin_avg.length; i++){
            for(int j=0;j<TmpCosmo_sin_avg[0].length; j++){
                for(int k=0;k<TmpCosmo_sin_avg[0][0].length; k++){
                    TmpCosmo_sin_avg[i][j][k] = cosmo_V10M[i][j][k]/cosmo_wind10M_magn[i][j][k];
                }
            }
        }

        double cosmo_cos_avg = nanmean2(TmpCosmo_cos_avg);
        double cosmo_sin_avg = nanmean2(TmpCosmo_sin_avg);

        double cosmo_dir_avg = Math.atan2(cosmo_sin_avg, cosmo_cos_avg);
        cosmo_dir_avg = changeDirRule(cosmo_dir_avg)/this.constants.getDeg2rad();


        // % Yamartino's method for wind std:
        // % http://journals.ametsoc.org/doi/pdf/10.1175/1520-0450%281984%29023%3C1362%3AACOSPE%3E2.0.CO%3B2
        double bb = -1 + 2/Math.sqrt(3.0);

        //Numero complesso... ma non lo uso ?
        double ecmwf_epsilon = Math.sqrt(1- (Math.pow(ecmwf_cos_avg, 2) + Math.pow(ecmwf_sin_avg, 2) ));
        double ecmwf_dir_std = Math.asin(ecmwf_epsilon) * (1+bb+Math.pow(ecmwf_epsilon, 2))/this.constants.getDeg2rad();


        //Numero complesso... ma non lo uso ?
        double cosmo_epsilon = Math.sqrt(Math.pow(cosmo_cos_avg, 2) + Math.pow(cosmo_sin_avg, 2));
        double cosmo_dir_std = Math.asin(cosmo_epsilon) * (1+bb+Math.pow(cosmo_epsilon,2))/this.constants.getDeg2rad();


        return new fieldStatsResults(ecmwf_dir_avg, ecmwf_dir_std, cosmo_dir_avg, cosmo_dir_std);
    }

    private double[][][] wave_dispersion(double[][][] wave_period){
        //    dispersion relation for monocromatic ocean waves
        //      -) If just wave period provided --> deep water approximation
        double twopi = 2*Math.PI;

        //deep water approximation
        double[][][] lambda = new double[wave_period.length][wave_period[0].length][wave_period[0][0].length];
        for(int i=0;i<lambda.length; i++){
            for(int j=0;j<lambda[0].length;j++){
                for(int k=0;k<lambda[0][0].length;k++){
                    lambda[i][j][k] = this.constants.getG0()/twopi* Math.pow(wave_period[i][j][k],2); //[m]
                }
            }
        }
        return lambda;
    }

    private double[][][] wave_dispersion(double[][][] wave_period, double[][][] depth){
        double[][][] lambda = wave_dispersion(wave_period);
        double[][][] fmk = Fenton_McKee_factor(wave_period,depth);
        // %   -) If also depth provided       --> generic depth formula after:
        // %   Fenton, JD and McKee, WD (1990)
        // %
        // %   wave_period [s]
        // %   depth       [m]
        // %   lambda      [m]
        // %
        if(this.forcing.getDeepWaterApprox()!=1){
            for(int i=0;i<lambda.length;i++){
                for(int j=0;j<lambda[0].length;j++){
                    for(int k=0;k<lambda[0][0].length;k++){
                        lambda[i][j][k] = lambda[i][j][k]*fmk[i][j][k];
                    }
                }
            }
        }
        return lambda;
    }

    private double[][][] Fenton_McKee_factor(double[][][] wave_period, double[][][] depth){
        // % Factor multiplying deep-water wavelength in
        // % Fenton, JD and McKee, WD (1990)
        // %
        // % (it leads to a reduced wavelength in shallow water)
        double twopi = Math.PI *2;
        double[][][] lambda_0 = new double[wave_period.length][wave_period[0].length][wave_period[0][0].length];
        for(int i=0;i<wave_period.length;i++){
            for(int j=0;j<wave_period[0].length;j++){
                for(int k=0;k<wave_period[0][0].length;k++){
                    lambda_0[i][j][k] = this.constants.getG0()/twopi * Math.pow(wave_period[i][j][k],2);
                }
            }
        }
        double[][][] fmk = new double[wave_period.length][wave_period[0].length][wave_period[0][0].length];
        for(int i=0;i<wave_period.length;i++){
            for(int j=0;j<wave_period[0].length;j++){
                for(int k=0;k<wave_period[0][0].length;k++){
                    Math.pow(Math.tanh(Math.pow(twopi*depth[i][j][k] / lambda_0[i][j][k],(3/4))),(2/3));
                }
            }
        }
        return fmk;
    }

//    private double[][] nanmean2(double[][][] A_mat){
//        //Abstract: compute mean of matrix elements of A_mat, even in presence of NaNs.
//        double[][][] matrix = A_mat;
//        for(int i=0;i<matrix.length;i++){
//            for(int j=0;j<matrix[0].length;j++){
//                for(int k=0;k<matrix[0][0].length;k++){
//                    if(Double.isNaN(matrix[i][j][k])){
//                        matrix[i][j][k] = 0.0;
//                    }
//                }
//            }
//        }
//        return Utility.mean3d(matrix);
//    }

    private double nanmean2(double[][][] A_mat){
        //Abstract: compute mean of matrix elements of A_mat, even in presence of NaNs.
        ArrayList<Double> elements = new ArrayList<>();
        for(int k=0;k<A_mat.length;k++){
            for(int i=0;i<A_mat[0].length;i++){
                for(int j=0;j<A_mat[0][0].length;j++){
                    if(!Double.isNaN(A_mat[k][i][j])){
                        elements.add(A_mat[k][i][j]);
                    }
                }
            }
        }
        int nElements = elements.size();
        double cumSum = 0.0;
        for(int i=0;i<elements.size();i++){
            cumSum+=elements.get(i);
        }
        return (cumSum/nElements);
    }

    private mFields_reductionResults mFields_reduction(double[] lat, double[] lon, double[][][] VTDH,
                                                       double[][][] VTPK, double [][][] VDIR){
        mFields_reductionResults retThis = new mFields_reductionResults();
        double meshRes = 1.0/this.sGrid.getInvStepFields();

        //reduction to given bounding box
        findResults latTmp = Utility.find(lat, "<=", this.sGrid.getBbox__lat_min()-meshRes, "<=", this.sGrid.getBbox__lat_max()+meshRes);
        int[] row_lat = latTmp.getIndexes();
        findResults lonTmp = Utility.find(lon, "<=", this.sGrid.getBbox__lon_min()-meshRes, "<=", this.sGrid.getBbox__lon_max()+meshRes);
        int[] row_lon = lonTmp.getIndexes();

        if(row_lon.length<this.sGrid.getMinNoGridPoints() || row_lat.length<this.sGrid.getMinNoGridPoints()){
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("check_start_timestep: too small or too narrow bounding box");
            debug.CloseFile();
            System.exit(0);
        }
        double[] lat_red = new double[row_lat.length];
        double[] lon_red = new double[row_lon.length];
        for(int i=0;i<row_lat.length;i++){
            lat_red[i]=lat[row_lat[i]];
        }
        for(int i=0;i<row_lon.length;i++){
            lon_red[i]=lon[row_lon[i]];
        }
        retThis.setLat_red(lat_red);
        retThis.setLon_red(lon_red);
        double[][][] out1 = null;
        double[][][] out2 = null;
        double[][][] out3 = null;
        if(VTDH != null){

            out1=new double[row_lon.length][VTDH[0].length][row_lat.length];
            for(int k=0;k<row_lon.length;k++){
                for(int i=0;i<out1[0].length;i++){
                    for(int j=0;j<row_lat.length;j++){
                        out1[k][i][j] = VTDH[row_lon[k]][i][row_lat[j]];
                    }
                }
            }
        }
        if(VTPK != null){
            out2 = new double[row_lon.length][VTPK[0].length][row_lat.length];
            for(int k=0;k<row_lon.length;k++){
                for(int i=0;i<VTPK[0].length;i++){
                    for(int j=0;j<row_lat.length;j++){
                        out2[k][i][j] = VTPK[row_lon[k]][i][row_lat[j]];
                    }
                }
            }
        }
        if(VDIR != null){
            out3 = new double[row_lon.length][VDIR[0].length][row_lat.length];
            for(int k=0;k<row_lon.length;k++){
                for(int i=0;i<VDIR[0].length;i++){
                    for(int j=0;j<row_lat.length;j++){
                        out3[k][i][j] = VDIR[row_lon[k]][i][row_lat[j]];
                    }
                }
            }
        }
        retThis.setOut1(out1);
        retThis.setOut2(out2);
        retThis.setOut3(out3);

        //percent data reduction:
        int n_elements=out1.length * out1[0].length * out1[0][0].length;
        int red_percent = ((1- n_elements)/n_elements)*100;
        this.logFile = new MyFileWriter("","",true);
        this.logFile.WriteLog("\t\tData reduction: "+red_percent+"%");
        this.logFile.CloseFile();
        return retThis;
    }

    private void check_start_timestep(long deltaHr_anls){
        if(deltaHr_anls <= 0){
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("mFields_reduction: departure time too far back in the past!");
            debug.CloseFile();
            System.exit(0);
        }

        if(this.forcing.getWave()==1){
            if(this.tGrid.getWave_dep_TS() < 1){
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("check_start_timestep: departure time too far back in the past!");
                debug.CloseFile();
                System.exit(0);
            }

            if(this.tGrid.getWave_dep_TS() > this.tGrid.getMaxNt()){
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("check_start_timestep: departure time too far in future!");
                debug.CloseFile();
                System.exit(0);
            }
        }

        if(this.forcing.getWind() == 1){
            if(this.tGrid.getWind_dep_TS() < 1){
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("check_start_timestep: departure time too far back in the past!");
                debug.CloseFile();
                System.exit(0);
            }

            if(this.tGrid.getWind_dep_TS() > this.tGrid.getMaxNt()){
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("check_start_timestep: departure time too far back in the past!");
                debug.CloseFile();
                System.exit(0);
            }
        }
    }

    private readout_envFieldsResults readout_envFields(){
        //reads wave and wind model forecast and (when relevant) analysis fields
        double[][][] VTDH;
        double[][][] VTPK;
        double[][][] VDIR;
        double[] lat_wave;
        double[] lon_wave;
        int[] wave_origtimes;
        double wave_maxNt = 0.0;
        if(this.forcing.getWave() == 1){
            //------------------------------------------------------------------------
            //wave model
            String model_str = "";
            if(this.optim.getWaveModel() < this.optim.getRelocAnlsCode()){
                model_str = "WW3";
            } else {
                if(this.optim.getWaveModel() >= this.optim.getRelocAnlsCode()){
                    model_str = "SWAN (reloc)";
                } else {
                    System.out.println("readout_envFields: unknown waveModel!");
                    MyFileWriter debug = new MyFileWriter("","debug",false);
                    debug.WriteLog("readout_envFields: unknown waveModel!");
                    debug.CloseFile();
                    System.exit(0);
                }
            }
            this.logFile = new MyFileWriter("","",true);
            this.logFile.WriteLog("\tattempt reading "+model_str+" wave forecast data...");
            this.logFile.CloseFile();

            readout_mWaveResults readout_mWave = readout_mWave();
            if(!readout_mWave.isWave_status()){
                System.out.println("readout_envFields: "+ readout_mWave.getWave_filename() + "wave forecast not found!");
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("readout_envFields: "+ readout_mWave.getWave_filename() + "wave forecast not found!");
                debug.CloseFile();
                System.exit(0);
            }
            lat_wave = readout_mWave.getLat();
            lon_wave = readout_mWave.getLon();
            VTDH = readout_mWave.getVTHD();
            VTPK = readout_mWave.getVTPK();
            VDIR = readout_mWave.getVDIR();
            wave_maxNt = VTDH[0].length;
            wave_origtimes = new int[(int)wave_maxNt];
            for(int i=0;i<(int)wave_maxNt;i++)
                wave_origtimes[i]=(i+1);
            this.logFile = new MyFileWriter("","",true);
            this.logFile.WriteLog("\tdeparture at time step # "+this.tGrid.getWave_dep_TS()+" of hourly interpolated file");
            this.logFile.CloseFile();
        } else{
            //fake wave fields in case of sailboat:
            fake_waveFieldsResults waveFields = fake_waveFields(192,192, 10, 30);
            wave_origtimes = waveFields.getWave_origtimes();
            lon_wave = waveFields.getLon_wave();
            lat_wave = waveFields.getLat_wave();
            VTDH = waveFields.getVTDH();
            VTPK = waveFields.getVTPK();
            VDIR = waveFields.getVDIR();
        }
        double[] ecmwf_wind_origTimes;
        double[] cosmo_wind_origTimes;
        double[] ecmwf_lat_wind;
        double[] cosmo_lat_wind;
        double[] ecmwf_lon_wind;
        double[] cosmo_lon_wind;
        double[][][] ecmwf_U10m;
        double[][][] cosmo_U10m;
        double[][][] ecmwf_V10m;
        double[][][] cosmo_V10m;
        double wind_maxNt =0.0;
        if(this.forcing.getWind() == 1){
            String star_str = " (**using this one**)";
            String read_str1 = " attempt reading ";
            String read_str2 = " wind forecast data ";
            String ecmwf_str = "";
            String cosmo_str = "";
            if(this.optim.getWindModel() == 11){ //ecmwf
                ecmwf_str = read_str1 + "ECMWF_1/4" + read_str2+ star_str;
                cosmo_str = read_str1 + "COSMO-ME" + read_str2;
            } else {
                if(this.optim.getWindModel() == 12){ //ecmwf
                    ecmwf_str = read_str1 + "ECMWF_1/8" + read_str2+ star_str;
                    cosmo_str = read_str1 + "COSMO-ME" + read_str2;
                } else {
                    if(this.optim.getWindModel() == 2){ //cosmo-me
                        ecmwf_str = read_str1 + "ECMWF " + read_str2+ star_str;
                        cosmo_str = read_str1 + "COSMO-ME_1/16" + read_str2;
                    } else {
                        System.out.println("readout_envFields: unknown windModel!");
                        MyFileWriter debug = new MyFileWriter("","debug",false);
                        debug.WriteLog("readout_envFields: unknown windModel!");
                        debug.CloseFile();
                        System.exit(0);
                    }
                }
            }
            //------------------------------------------------------------------------
            this.logFile = new MyFileWriter("","",true);
            this.logFile.WriteLog("\t"+ecmwf_str);
            this.logFile.CloseFile();

            ecmwf_lat_wind = new double[0];
            ecmwf_lon_wind = new double[0];
            cosmo_lat_wind = new double[0];
            cosmo_lon_wind = new double[0];
            ecmwf_U10m = new double[0][0][0];
            cosmo_U10m = new double[0][0][0];
            ecmwf_V10m = new double[0][0][0];
            cosmo_V10m = new double[0][0][0];
            ecmwf_wind_origTimes = new double[0];
            cosmo_wind_origTimes = new double[0];
        } else { //Fake wind fields in case of motorboat:
            fake_windFieldsResults res = fake_windFields(192, 45, 10, 30);
            ecmwf_wind_origTimes = res.getWind_origTimes();
            cosmo_wind_origTimes = res.getWind_origTimes();
            ecmwf_lat_wind = res.getLat_wind();
            cosmo_lat_wind = res.getLat_wind();
            ecmwf_lon_wind = res.getLon_wind();
            cosmo_lon_wind = res.getLon_wind();
            ecmwf_U10m = res.getU10m();
            cosmo_U10m = res.getU10m();
            ecmwf_V10m = res.getV10m();
            cosmo_V10m = res.getV10m();
        }
        //------------------------------------------------------------------------
        //Tgrid.maxNt:
        if(this.forcing.getWave()!=1){
            wave_maxNt = Double.NaN;
        }
        if(this.forcing.getWind()!=1){
            wind_maxNt = Double.NaN;
        }
        this.tGrid.setMaxNt(Utility.min(wave_maxNt, wind_maxNt));

        return new readout_envFieldsResults(lat_wave, lon_wave, ecmwf_lat_wind, ecmwf_lon_wind, cosmo_lat_wind, cosmo_lon_wind, VTDH, VTPK, VDIR,
                ecmwf_U10m, ecmwf_V10m, cosmo_U10m, cosmo_V10m, wave_origtimes, ecmwf_wind_origTimes, cosmo_wind_origTimes);
    }

    private fake_windFieldsResults fake_windFields(int ntwind, int nsteps, int nlat, int nlon){
        ArrayList<Double> lat_windTmp = Utility.linspace(this.sGrid.getBbox__lat_min(), this.sGrid.getBbox__lat_max(), nlat);
        double[] lat_wind = new double[lat_windTmp.size()];
        for(int i=0;i<lat_windTmp.size(); i++)
            lat_wind[i] = lat_windTmp.get(i);
        ArrayList<Double> lon_windTmp = Utility.linspace(this.sGrid.getBbox__lon_min(), this.sGrid.getBbox__lon_max(), nlon);
        double[] lon_wind = new double[lon_windTmp.size()];
        for(int i=0;i<lon_windTmp.size(); i++)
            lon_wind[i] = lon_windTmp.get(i);
        //FORSE
//        double[][][] U10m = Utility.ones3Dmatrix(nsteps, nlat, nlon);
//        double[][][] V10m = Utility.ones3Dmatrix(nsteps, nlat, nlon);
        double[][][] U10m = Utility.ones3Dmatrix(nlon, nsteps, nlat);
        double[][][] V10m = Utility.ones3Dmatrix(nlon, nsteps, nlat);
        ArrayList<Double> wind_origTimesTmp = Utility.linspace(0, ntwind, nsteps);
        double[] wind_origTimes = new double[wind_origTimesTmp.size()];
        for(int i=0;i<wind_origTimesTmp.size();i++)
            wind_origTimes[i]=wind_origTimesTmp.get(i);
        return new fake_windFieldsResults(wind_origTimes, lat_wind, lon_wind, U10m, V10m);
    }

    private fake_waveFieldsResults fake_waveFields(int mtwave, int nsteps, int nlat, int nlon){
        ArrayList<Double> lat_waveTmp = Utility.linspace(this.sGrid.getBbox__lat_min(), this.sGrid.getBbox__lat_max(), nlat);
        double[] lat_wave = new double[lat_waveTmp.size()];
        for(int i=0;i<lat_waveTmp.size();i++){
            lat_wave[i]=lat_waveTmp.get(i);
        }
        ArrayList<Double> lon_waveTmp = Utility.linspace(this.sGrid.getBbox__lon_min(), this.sGrid.getBbox__lon_max(), nlon);
        double[] lon_wave = new double[lon_waveTmp.size()];
        for(int i=0;i<lon_waveTmp.size();i++){
            lon_wave[i]=lon_waveTmp.get(i);
        }
        //FORSE
//        double[][][] VTDH = Utility.ones3Dmatrix(nsteps, nlat, nlon);
//        double[][][] VTPK = Utility.ones3Dmatrix(nsteps, nlat, nlon);
//        double[][][] VDIR = Utility.ones3Dmatrix(nsteps, nlat, nlon);
        double[][][] VTDH = Utility.ones3Dmatrix(nlon,nsteps,nlat);
        double[][][] VTPK = Utility.ones3Dmatrix(nlon,nsteps,nlat);
        double[][][] VDIR = Utility.ones3Dmatrix(nlon,nsteps,nlat);
        int[] wave_origtimes = new int[mtwave];
        for(int i=0;i<mtwave; i++)
            wave_origtimes[i]=(i+1);
        return new fake_waveFieldsResults(wave_origtimes, lon_wave, lat_wave, VTDH, VTPK, VDIR);
    }

    private readout_mWaveResults readout_mWave(){
        // % reads out Waves forecast data
        // % physical fields from either WAM or WW3 model:
        String Stagein_path = "";
        if(this.optim.getWaveModel() == 10){
            Stagein_path = "inputFiles/wave/WW3/analysis";
        } else {
            if(this.optim.getWaveModel() == 1){
                Stagein_path = "inputFiles/wave/WW3/forecast";
            } else {
                if(this.optim.getWaveModel() >= 20){//relocatable
                    Stagein_path = "inputFiles/wave/SWAN";
                }
            }
        }

        String wave_filename = Stagein_path+"/start__"+this.tGrid.getLatest_date()+".nc";
        MyNetCDFParser parser = new MyNetCDFParser(wave_filename);
        boolean wave_status = parser.isFileExists();
        double[] lat;
        double[] lon;
        double[][][] VTDH;
        double[][][] VDIR;
        double[][][] VTPK;
        if(wave_status){
            waveForecastResults waveForecast = parser.parseWaveForecastData();
            lat= waveForecast.getLatitude();
            lon = waveForecast.getLongitude();

            VTDH = waveForecast.getVTDH();
            VDIR = waveForecast.getVDIR();
            VTPK = waveForecast.getVTPK();

            //############################ WAM convention used in VISIR!! ###################
            // if optim.waveModel<optim.relocAnlsCode
            //     % converting WW3 wave directions to WAM directions convention:
            //     if numel(regexp(hostflag, 'okeanos'))==0
            //         VDIR= VDIR + 180* sign(180- VDIR);
            //         VDIR(VDIR==180) = 0;   % this accounts for the fact that, in Matlab, sign(0)=0
            //     end
            // elseif optim.waveModel>=optim.relocAnlsCode
            //     VY=  cos(const.deg2rad*VDIR);
            //     VX=  sin(const.deg2rad*VDIR);
            //     VDIR= atan2(VY, VX)/const.deg2rad;
            // end

            if(this.optim.getWaveModel() >= this.optim.getRelocAnlsCode()){
                double[][][] VY = new double[VDIR.length][VDIR[0].length][VDIR[0][0].length];
                double[][][] VX = new double[VDIR.length][VDIR[0].length][VDIR[0][0].length];
                for(int i=0;i<VDIR.length; i++){
                    for(int j=0;j<VDIR[0].length;j++){
                        for(int k=0;k<VDIR[0][0].length; k++){
                            VY[i][j][k]=Math.cos(this.constants.getDeg2rad()*VDIR[i][j][k]);
                            VX[i][j][k]=Math.sin(this.constants.getDeg2rad()*VDIR[i][j][k]);
                            VDIR[i][j][k] = Math.atan2(VY[i][j][k], VX[i][j][k])/this.constants.getDeg2rad();
                        }
                    }
                }
            }

            //######################################################################

            // %     GM's note 23/5/2016
            // %     WAM (WW3 pre13Jun2016) convention is the oceanographic (meteorological) convention
            // %
            // %     starting since 13/6/2016, INGV will provide VDIR with oceanographic convention
            // %     Thus, the conversion:
            // %          VDIR= VDIR + 180* sign(180- VDIR);
            // %         VDIR(VDIR==180) = 0;   % this accounts for the fact that, in Matlab, sign(0)=0
            // %     will have to be dismissed.

        } else{
            lat = new double[0];
            lon = new double[0];
            VTDH = new double[0][0][0];
            VTPK = new double[0][0][0];
            VDIR = new double[0][0][0];
        }
        return new readout_mWaveResults(lat, lon, wave_status, wave_filename, VTDH, VTPK, VDIR);
    }

    private prepare_sqrtY_fieldResults prepare_sqrtY_field(){
        // % preparation of a pseudo- waveheight field
        // % depending on sqrt of Y coordinate
        // % to be used for testing of analytical solution (=cycloid)

        //x, y
        double[] tmpLat = new double[1];
        tmpLat[0]=this.extreme_pts.getStart_lat();
        double[] tmpLon = new double[1];
        tmpLon[0]=this.extreme_pts.getStart_lon();

        deg2utmResults tmp = deg2utm(tmpLat, tmpLon);
        double x_start = tmp.getX()[0];
        double y_start = tmp.getY()[0];
        String[] utmzone_start = tmp.getUtmzone()[0];

        tmpLat[0]=this.extreme_pts.getEnd_lat();
        tmpLon[0]=this.extreme_pts.getEnd_lon();

        tmp = deg2utm(tmpLat,tmpLon);
        double x_end = tmp.getX()[0];
        double y_end = tmp.getY()[0];
        String[] utm_zone_end = tmp.getUtmzone()[0];

        int zoneN_start = Integer.parseInt(utmzone_start[0]);

        meshgridResults res = Utility.meshgrid(this.lat_bathy_Inset, this.lon_bathy_Inset);
        double[][] lat_gr = res.getX();
        double[][] lon_gr = res.getY();
        int Np = lat_gr.length * lat_gr[0].length;

        int[] utmzone_number = new int[Np];
        Arrays.fill(utmzone_number, zoneN_start);
        degzone2utmResults dz2uTmp = degzone2utm(Utility.reshape(lat_gr,Np), Utility.reshape(lon_gr, Np), utmzone_number);
        double[] xx=dz2uTmp.getX();
        double[] yy=dz2uTmp.getY();

        int Ny = this.lat_bathy_Inset.size();//43
        int Nx = this.lon_bathy_Inset.size();//72

        int[] dim = new int[2];
        dim[0]=Nx;
        dim[1]=Ny;
        double[][] Nyy = Utility.reshape(yy,dim);

        //Dy
        double[][] DeltaY = new double[Nx][Ny];
        if(this.extreme_pts.getCycType() == "id"){//dd-type cycloid:
            for(int ix = 0;ix<Nx; ix++){
                double[] y_diff = new double[Ny];
                for(int j =0;j<Ny;j++){
                    y_diff[j] = Math.abs(Nyy[ix][j] - y_start); //m
                }
                minResults min = Utility.minWithIndex(y_diff);
                double y_val = min.getElement();
                int y_idx = min.getIndex();
                for(int j=0;j<Ny;j++){
                    DeltaY[ix][j] = Nyy[ix][y_idx] - Nyy[ix][j]; //m
                }
            }
        } else { //id-type cycloid:
            for(int ix = 0;ix<Nx; ix++){
                double[] y_diff = new double[Ny];
                for(int j =0;j<Ny;j++){
                    y_diff[j] = Math.abs(Nyy[ix][j] - y_end); //m
                }
                minResults min = Utility.minWithIndex(y_diff);
                double y_val = min.getElement();
                int y_idx = min.getIndex();
                for(int j=0;j<Ny;j++){
                    DeltaY[ix][j] = Nyy[ix][j] - Nyy[ix][y_idx]; //m
                }
            }
        }
        for(int u =0 ;u<Nx; u++){
            for(int v =0;v<Ny;v++){
                if(DeltaY[u][v]<0){
                    DeltaY[u][v] = 0;
                }
            }
        }

        //speed
        this.extreme_pts.setPseudoG(0.001);
        double[][][] VTDH_b = new double[Ny][(int) this.tGrid.getNt()][Nx];
        double[][] DeltaYTransposed = Utility.transposeMatrix(DeltaY);
        for(int i=0;i<Ny;i++){
            for(int it=0;it<(int) this.tGrid.getNt(); it++){
                for(int j=0;j<Nx;j++){
                    VTDH_b[i][it][j] = this.constants.getMs2kts() * Math.sqrt(2*this.extreme_pts.getPseudoG()*DeltaYTransposed[i][j]);//kts
                }
            }
        }
        double[][][] VDIR_b = Utility.NaN3Dmatrix(Ny, (int) this.tGrid.getNt(), Nx);
        //FORSE??
//        double[][][] VTDH_b = new double[(int) this.tGrid.getNt()][Ny][Nx];
//        double[][] DeltaYTransposed = Utility.transposeMatrix(DeltaY);
//        for(int i=0;i<(int) this.tGrid.getNt();i++){
//            for(int j=0;j<Ny;j++){
//                for(int k=0;k<Nx;k++){
//                    VTDH_b[i][j][k] = this.constants.getMs2kts() * Math.sqrt(2*this.extreme_pts.getPseudoG() * DeltaYTransposed[j][k]); //kts
//                }
//            }
//        }
//        double[][][] VDIR_b = Utility.NaN3Dmatrix((int) this.tGrid.getNt(), Ny, Nx);

        //graphical pars
        double smallFract = 0.5;
        this.fstats.setWheight_min((1-smallFract)*Utility.min3d(VTDH_b));
        this.fstats.setWheight_max((1+smallFract)*Utility.max3d(VTDH_b));
        return new prepare_sqrtY_fieldResults(VTDH_b, VDIR_b);
    }

    private degzone2utmResults degzone2utm(double[] Lat, double[] Lon, int[] utmzone_number){
        // % utmzone_number can be forced to be different from the proper one.
        // % based on deg2utm; modified by G.Mannarini on dec.15, 2008

        //Argument checking
        int n1 = Lat.length;
        int n2 = Lon.length;
        int n3 = utmzone_number.length;

        if(n1!=n2){
            System.out.println("degzone2utm: Lat and Lon vectors should have the same length");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("degzone2utm: Lat and Lon vectors should have the same length");
            debug.CloseFile();
            System.exit(0);
        } else{
            if(n1!=n3){
                System.out.println("degzone2utm: Lat and Lon vectors should have the same length of utmzone_number vector");
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("degzone2utm: Lat and Lon vectors should have the same length of utmzone_number vector");
                debug.CloseFile();
                System.exit(0);
            }
        }

        //Memory pre-allocation
        double[] x = new double[n1];
        double[] y = new double[n1];

        //Main loop
        for(int i=0;i<n1;i++){
            double la = Lat[i];
            double lo = Lon[i];

            double sa = 6378137.000000 ;
            double sb = 6356752.314245;

            double e2 = ( Math.pow((Math.pow(sa,2) - Math.pow(sb,2)),0.5) )/sb;
            double e2cuadrada = Math.pow(e2,2);
            double c = Math.pow(sa,2)/sb;

            double lat = la * (Math.PI/180);
            double lon = lo * (Math.PI/180);

            int Huso = utmzone_number[i];
            int S =( (Huso * 6) - 183);
            double deltaS = lon - ( S* (Math.PI/180));

            double a = Math.cos(lat) * Math.sin(deltaS);
            double epsilon = 0.5 * Math.log((1+a)/(1-a));
            double nu = Math.atan(Math.tan(lat)/Math.cos(deltaS))-lat;
            double v = (c/(Math.pow((Math.pow((e2cuadrada*Math.cos(lat)),2)),0.5)))*0.9996;
            double ta = (e2cuadrada/2) * Math.pow(epsilon,2) + Math.pow(Math.cos(lat),2);
            double a1 = Math.sin(2*lat);
            double a2 = a1 * Math.pow(Math.cos(lat),2);
            double j2 = lat + (a1/2);
            double j4 = ((3*j2) + a2)/4;
            double j6 = (5*j4) + (a2*Math.pow(Math.cos(lat),2)) /3;
            double alfa = (3/4)*e2cuadrada;
            double beta = (5/3)*Math.pow(alfa,2);
            double gama = (35/27)*Math.pow(alfa,3);
            double Bm = 0.9996 * c * (lat - alfa * j2 + beta * j4 - gama * j6);
            double xx = epsilon * v * (1+(ta/3))+500000;
            double yy = nu*v*(1+ta)+Bm;

            if(yy<0)
                yy+=9999999;

            x[i] = xx;
            y[i] = yy;
        }

        return new degzone2utmResults(x,y);
    }

    private deg2utmResults deg2utm(double[] Lat, double[] Lon){
        // % Description: Function to convert lat/lon vectors into UTM coordinates (WGS84).
    // % Some code has been extracted from UTM.m function by Gabriel Ruiz Martinez.
    // %
    // % Inputs:
    // %    Lat: Latitude vector.   Degrees.  +ddd.ddddd  WGS84
    // %    Lon: Longitude vector.  Degrees.  +ddd.ddddd  WGS84
    // %
    // % Outputs:
    // %    x, y , utmzone.   See example
    // %
    // % Example 1:
    // %    Lat=[40.3154333; 46.283900; 37.577833; 28.645650; 38.855550; 25.061783];
    // %    Lon=[-3.4857166; 7.8012333; -119.95525; -17.759533; -94.7990166; 121.640266];
    // %    [x,y,utmzone] = deg2utm(Lat,Lon);
    // %    fprintf('%7.0f ',x)
    // %       458731  407653  239027  230253  343898  362850
    // %    fprintf('%7.0f ',y)
    // %      4462881 5126290 4163083 3171843 4302285 2772478
    // %    utmzone =
    // %       30 T
    // %       32 T
    // %       11 S
    // %       28 R
    // %       15 S
    // %       51 R
    // %
    // % Example 2: If you have Lat/Lon coordinates in Degrees, Minutes and Seconds
    // %    LatDMS=[40 18 55.56; 46 17 2.04];
    // %    LonDMS=[-3 29  8.58;  7 48 4.44];
    // %    Lat=dms2deg(mat2dms(LatDMS)); %convert into degrees
    // %    Lon=dms2deg(mat2dms(LonDMS)); %convert into degrees
    // %    [x,y,utmzone] = deg2utm(Lat,Lon)
    // %
    // % Author:
    // %   Rafael Palacios
    // %   Universidad Pontificia Comillas
    // %   Madrid, Spain
    // % Version: Apr/06, Jun/06, Aug/06, Aug/06
    // % Aug/06: fixed a problem (found by Rodolphe Dewarrat) related to southern
    // %    hemisphere coordinates.
    // % Aug/06: corrected m-Lint warnings
    // %-------------------------------------------------------------------------

        //Argument checking
        int n1 = Lat.length;
        int n2 = Lon.length;

        if(n1!=n2){
            System.out.println("Lat and Lon vectors should have the same length");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("deg2utm: Lat and Lon vectors should have the same length");
            debug.CloseFile();
            System.exit(0);
        }

        //Memory pre-allocation
        double[] x = new double[n1];
        double[] y = new double[n1];
        String[][] utmzone = new String[n1][2];

        for(int i=0;i<n1;i++){
            double la = Lat[i];
            double lo = Lon[i];

            double sa = 6378137.000000;
            double sb = 6356752.314245;

            double e2 = ( Math.pow((Math.pow(sa,2) - Math.pow(sb,2)),0.5) )/sb;
            double e2cuadrada = Math.pow(e2,2);
            double c = Math.pow(sa,2)/sb;

            double lat = la * (Math.PI/180);
            double lon = lo * (Math.PI/180);

            int Huso = Utility.fix( (lo/6) + 31);
            int S =( (Huso * 6) - 183);
            double deltaS = lon - ( S* (Math.PI/180));

            String Letra = "X";

            if(la<-72)
                Letra = "C";
            else if(la<-64)
                Letra = "D";
            else if(la<-56)
                Letra = "E";
            else if(la<-48)
                Letra = "F";
            else if(la<-40)
                Letra = "G";
            else if(la<-32)
                Letra = "H";
            else if(la<-24)
                Letra = "J";
            else if(la<-16)
                Letra = "K";
            else if(la<-8)
                Letra = "L";
            else if(la<0)
                Letra = "M";
            else if(la<8)
                Letra = "N";
            else if(la<16)
                Letra = "P";
            else if(la<24)
                Letra = "Q";
            else if(la<32)
                Letra = "R";
            else if(la<40)
                Letra = "S";
            else if(la<48)
                Letra = "T";
            else if(la<56)
                Letra = "U";
            else if(la<64)
                Letra = "V";
            else if(la<72)
                Letra = "W";
            else
                Letra ="X";

            double a = Math.cos(lat) * Math.sin(deltaS);
            double epsilon = 0.5 * Math.log((1+a)/(1-a));
            double nu = Math.atan(Math.tan(lat)/Math.cos(deltaS))-lat;
            double v = (c/(Math.pow((Math.pow((e2cuadrada*Math.cos(lat)),2)),0.5)))*0.9996;
            double ta = (e2cuadrada/2) * Math.pow(epsilon,2) + Math.pow(Math.cos(lat),2);
            double a1 = Math.sin(2*lat);
            double a2 = a1 * Math.pow(Math.cos(lat),2);
            double j2 = lat + (a1/2);
            double j4 = ((3*j2) + a2)/4;
            double j6 = (5*j4) + (a2*Math.pow(Math.cos(lat),2)) /3;
            double alfa = (3/4)*e2cuadrada;
            double beta = (5/3)*Math.pow(alfa,2);
            double gama = (35/27)*Math.pow(alfa,3);
            double Bm = 0.9996 * c * (lat - alfa * j2 + beta * j4 - gama * j6);
            double xx = epsilon * v * (1+(ta/3))+500000;
            double yy = nu*v*(1+ta)+Bm;

            if(yy<0)
                yy+=9999999;

            x[i] = xx;
            y[i] = yy;
            utmzone[i][0] = ""+Huso;
            utmzone[i][1] = Letra;
        }

        return new deg2utmResults(x,y,utmzone);
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

    public double[][] getVel_LUT() {
        return vel_LUT;
    }

    public ArrayList<Double> getH_array_m() {
        return H_array_m;
    }

    /******************************************************/
}
