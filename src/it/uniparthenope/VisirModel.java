package it.uniparthenope;

import java.io.FileReader;
import java.util.*;
//Boxing classes
import it.uniparthenope.Boxing.mdata_gridResults;
import it.uniparthenope.Boxing.meshgridResults;
import it.uniparthenope.Boxing.grid_extreme_coordsResults;
import it.uniparthenope.Boxing.idx_ref2inset_gridResults;
import it.uniparthenope.Boxing.readout_bathyResults;
import it.uniparthenope.Boxing.parseMedOneMinResults;
import it.uniparthenope.Parser.MyBinaryParser;
import it.uniparthenope.Parser.MyCSVParser;
import it.uniparthenope.Parser.MyNetCDFParser;

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
    private Double[][] xy;
    private Double[][] xg;
    private Double[][] yg;
    private Double[] xg_array;
    private Double[] yg_array;
    private Double[][] xy_DB;
    private ArrayList<Double> lat_bathy_Inset;
    private ArrayList<Double> lon_bathy_Inset;
    private Double[][] bathy_Inset;
    private Double[][] lsm_mask;
    private Double[][] J_mask;
    ArrayList<Double> x_islands;
    ArrayList<Double> y_islands;
    ArrayList<Double> x_continent;
    ArrayList<Double> y_continent;
    private Double[] estGdtDist;


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
        readout_bathyResults bathymetry = this.readout_bathy(bathy_code);
        ArrayList<Double> lat_bathy = bathymetry.getLat();
        ArrayList<Double> lon_bathy = bathymetry.getLon();
        Double[][] z_bathy = bathymetry.getZ();
        this.sGrid.setInv_step(1.0/Math.abs(lon_bathy.get(1)-lon_bathy.get(0)));

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
            Double min_lon_bathy = Collections.min(this.lon_bathy_Inset);
            Double max_lon_bathy = Collections.max(this.lon_bathy_Inset);
            Double min_lat_bathy = Collections.min(this.lat_bathy_Inset);
            Double max_lat_bathy = Collections.max(this.lat_bathy_Inset);
            for(Double element : x_coast){
                if((element >= min_lon_bathy) && (element <= max_lon_bathy)){
                    x_bool.add(true);
                } else { x_bool.add(false); }
            }
            for(Double element : y_coast){
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
        Double[][] xg_DB = data_gridOut.getXg();
        Double[][] yg_DB = data_gridOut.getYg();

        //inset grid plaid coordinates:

        mdata_gridResults data_gridOut2 = mdata_grid(this.lat_bathy_Inset, this.lon_bathy_Inset);
        this.xy = data_gridOut2.getXy();
        this.xg = data_gridOut2.getXg();
        this.yg = data_gridOut2.getYg();

        if(this.visualization.getGraphData()==1){
            //csv_write xy
            try{
                MyCSVParser csv = new MyCSVParser("output/GRAPH.node_LonLat.csv");
                csv.writeCSV(this.xy);
            } catch (Exception e){
                e.printStackTrace();
            }
        }

        //--------------------------------------------------------------------
        //Joint coast-vessel safety mask:
        System.out.println("computation of a joint coast-vessel safety mask...");
        Double[][] UKC_mask = Utility.ones(bathy_Inset.length, bathy_Inset[0].length);
        for(int i=0;i<UKC_mask.length;i++){
            for(int j=0;j<UKC_mask[0].length;j++){
                if(bathy_Inset[i][j] < this.ship.getDraught()){
                    UKC_mask[i][j] = Double.NaN;
                }
            }
        }
        Double[][] xg_Jmasked;
        Double[][] yg_Jmasked;
        Double[][] J_bathy_Inset;
        this.xg_array = new Double[0];
        this.yg_array = new Double[0];
        Double[][] xy_g;
        Double[][] dist_mask;
        Double[][] min_coast_dist;
        if(this.bar_flag == 2){
            //readout free nodes DB:
            MyBinaryParser datFile = new MyBinaryParser("inputFiles/graph/freeNodes_DB.dat");
            long[] free_nodes_DB = datFile.readAsUInt32();
            //remapping free nodes:
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
                System.exit(0);
            } else if(free_nodes_Number > this.sGrid.getNodesLargeN()){
                System.out.println("too large graph");
                System.exit(0);
            }
            this.sGrid.setFreenodes(free_nodes_Number);
            //lsm, created using coastline DB (NaNs on landmass):
            lsm_mask = Utility.NaNmatrix(xg.length,xg[0].length);
            for(int i=0;i<free_nodes.length;i++){
                lsm_mask[0][(int) free_nodes[i]] = 1.0;
            }
            int[] dim = new int[2];
            dim[0]=xg.length;
            dim[1]=xg[0].length;
            lsm_mask = Utility.reshape(lsm_mask,dim);

            //Safe distance from coastline:
            xg_array = Utility.reshape(xg,dim[0]*dim[1]);
            yg_array = Utility.reshape(yg, yg.length*yg[0].length);
            xy_g = new Double[xg_array.length][2];

            int nC = x_coast_Inset.size();
            int cols = dim[0]*dim[1];
            if(nC>0){
                Double[][] coast_dist = Utility.zeros(nC, cols);
                Double[][] P_b = new Double[1][2];
                for(int i=0;i<nC;i++){
                    P_b[0][0]=x_coast_Inset.get(i);
                    P_b[0][1]=y_coast_Inset.get(i);
                    Double[] hor_dist = hor_distance("s",xy_g, P_b);
                    for(int j=0; j<cols;j++){
                        coast_dist[i][j]= hor_dist[j];
                    }
                }
                Double[] min_coast_distTmp = Utility.min(coast_dist,1);
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

            // %###############################################
            // %
            // %  Joint safe mask:
            // %
            // J_mask = lsm_mask' .* UKC_mask .* dist_mask' ;
            // %
            // %###############################################
            J_mask = Utility.MatrixComponentXcomponent(Utility.MatrixComponentXcomponent(Utility.transposeMatrix(lsm_mask),Utility.transposeMatrix(UKC_mask)),Utility.transposeMatrix(dist_mask));

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
            xy_g = new Double[xg_Jmasked.length*xg_Jmasked[0].length][2];
            for(int i =0 ;i<xg_Jmasked.length*xg_Jmasked[0].length;i++){
                xy_g[i][0]=xg_array[i];
                xy_g[i][1]=yg_array[i];
            }
            //distance from start/end nodes:
            Double[][] tmpP_b=new Double[1][2];
            tmpP_b[0][0]=this.extreme_pts.getStart_lon();
            tmpP_b[0][1]=this.extreme_pts.getStart_lat();
            Double[] start_dist_matrix = hor_distance("s",xy_g,tmpP_b);
            tmpP_b[0][0]=this.extreme_pts.getEnd_lon();
            tmpP_b[0][1]=this.extreme_pts.getEnd_lat();
            Double[] end_dist_matrix = hor_distance("s",xy_g,tmpP_b);

            this.sGrid.setMin_start_dist(Utility.min(start_dist_matrix));
            this.sGrid.setMin_end_dist(Utility.min(end_dist_matrix));

            if(Double.isNaN(this.sGrid.getMin_start_dist()) || Double.isNaN(this.sGrid.getMin_end_dist())){
                System.out.println("departure or arrival point not compliant with safety specifications");
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
            this.sGrid.setNode_end_lon(yg_array[(int) this.sGrid.getNode_end()]);

            Double[][] tmpP_a = new Double[1][2];
            tmpP_a[0][0] = this.sGrid.getNode_start_lon();
            tmpP_a[0][1] = this.sGrid.getNode_start_lat();
            tmpP_b[0][0] = this.sGrid.getNode_end_lon();
            tmpP_b[0][1] = this.sGrid.getNode_end_lat();
            estGdtDist = hor_distance("s", tmpP_a, tmpP_b);

            tmpP_a[0][0] = this.sGrid.getNode_end_lon();
            tmpP_a[0][1] = this.sGrid.getNode_end_lat();
            tmpP_b[0][0] = this.sGrid.getNode_end_lon();
            tmpP_b[0][1] = this.sGrid.getNode_start_lat();
            Double[] delta_y = hor_distance("s", tmpP_a, tmpP_b);
            int sign = Utility.sign(this.sGrid.getNode_end_lat()-this.sGrid.getNode_start_lat());
            for(int i=0;i<delta_y.length;i++){
                delta_y[i]=delta_y[i]*sign;
            }

            tmpP_a[0][0] = this.sGrid.getNode_end_lon();
            tmpP_a[0][1] = this.sGrid.getNode_start_lat();
            tmpP_b[0][0] = this.sGrid.getNode_start_lon();
            tmpP_b[0][1] = this.sGrid.getNode_start_lat();
            Double[] delta_x = hor_distance("s", tmpP_a, tmpP_b);
            sign = Utility.sign(this.sGrid.getNode_end_lon()-this.sGrid.getNode_start_lon());
            for(int i=0;i<delta_y.length;i++){
                delta_y[i]=delta_y[i]*sign;
            }

            //orientation of the gdt route
            Double[] atan2 = new Double[delta_x.length];
            for(int i=0;i<delta_x.length;i++){
                atan2[i]=Math.atan2(delta_y[i],delta_x[i]);
            }
            this.sGrid.setTheta_gdt(atan2);
            //th= Sgrid.theta_gdt/const.deg2rad

        } else if(this.bar_flag==1){
            lsm_mask = Utility.NaNmatrix(xg.length,xg[0].length);
            Double[][] tmp = Utility.transposeMatrix(lsm_mask);
            J_mask = Utility.NaNmatrix(tmp.length,tmp[0].length);
            //-----------------------------------------------------------------------------------------------------
            //Start and end nodes:
            xg_array = new Double[0];
            yg_array = new Double[0];
            estGdtDist = new Double[0];
        }
        this.x_islands = this.lon_int;
        this.y_islands = this.lat_int;
        this.x_continent = this.lon_ext;
        this.y_continent = this.lat_ext;
    }

    private Double[] changeDirRule(Double[] inField){
        // % change of wind/current directional convention:
        // %
        // % from atan2 output concention to WAM-like convention, i.e.:
        // % from:  [-pi, pi] counterclockwise with 0 at  3:00 o'clock
        // %   to:  [0, 2*pi]        clockwise with 0 at 12:00 o'clock
        // %
        // % (inField and outField must be in radians)
        // %
        //---------------------------------------------
        Double[] inField1 = inField;
        for(int i=0;i<inField1.length;i++){
            if(inField1[i]<0){
                inField1[i]+=(2*Math.PI);
            }
        }

        //---------------------------------------------
        Double[] inField2 = inField1;
        for(int i=0;i<inField2.length;i++){
            inField2[i]-=(Math.PI/2);
            if(inField2[i]<0){
                inField2[i]+=(2*Math.PI);
            }
        }
        //---------------------------------------------
        Double[] outFiled = new Double[inField2.length];
        for(int i=0;i<inField2.length;i++){
            outFiled[i] = -inField2[i]*(2*Math.PI);
        }
        return outFiled;
    }

    private Double[] hor_distance(String method, Double[][] P_a, Double[][] P_b, Double... varargin){
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
        Double grid_step_in_NM = varargin.length > 0 ? varargin[0] : 1.0;
        Double[] xa = new Double[P_a.length];
        Double[] ya = new Double[P_a.length];
        for(int i=0;i<P_a.length;i++){
            xa[i]=P_a[i][0];
            ya[i]=P_a[i][1];
        }

        Double[] xb = new Double[P_b.length];
        Double[] yb = new Double[P_b.length];
        for(int i=0;i<P_b.length;i++){
            xb[i]=P_b[i][0];
            yb[i]=P_b[i][1];
        }

        Double[] dd=new Double[P_a.length];
        if(method=="plane"||method=="p"){
            for(int i=0;i<P_a.length;i++){
                dd[i]=grid_step_in_NM*Math.sqrt(Math.pow((xa[i]-xb[i]),2) + Math.pow((ya[i]-yb[i]),2));
            }
        } else if(method=="sphere" || method=="s"){
            // % from : http://mathworld.wolfram.com/GreatCircle.html
            // % x and y must be in degree.
            // % output in Nautical Miles (NM)
            Double E_radius = 3444.0; //NM
            for(int i=0;i<P_a.length;i++){
                dd[i]=E_radius*Math.acos(Math.cos(this.constants.getDeg2rad()*ya[i])*Math.cos(this.constants.getDeg2rad()*yb[i])*
                Math.cos(this.constants.getDeg2rad()*(xa[i]-xb[i]))*Math.sin(this.constants.getDeg2rad()*ya[i])*
                Math.sin(this.constants.getDeg2rad()*yb[i]));
            }
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
            System.out.println("inset grid not within reference grid !");
            System.exit(0);
        }
        if((Utility.any(idx_big,"<",1)) || Utility.any(idx_big,">",(nx_big*ny_big))){
            System.out.println("grid index not within reference grid !");
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
        Double[][] yg = out.getX();
        Double[][] xg = out.getY();
        int[] aDim = new int[2];
        aDim[0]=NN;
        aDim[1]=1;
        Double[][] xa = Utility.reshape(xg, aDim);
        Double[][] ya = Utility.reshape(yg, aDim);
        Double[][] xy = new Double[NN][2];
        for(int i=0;i<NN;i++){
            xy[i][0]=xa[i][0];
            xy[i][1]=ya[i][0];
        }
        return new mdata_gridResults(xy,xg,yg);
    }

    private grid_extreme_coordsResults grid_extreme_coords(ArrayList<Double> lat, ArrayList<Double> lon, Double[][] field_in){
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
        Double[][] field_out = null;
        //*****not sure of this*****
        //if barflag ==1, we are in creating graph DB mode, so we return empty lat_red, lon_red and field_out
        if(bar_flag == 2){//reading DB mode
            //Insert grid
            ArrayList<Integer> lat_row = new ArrayList<>();
            ArrayList<Integer> lon_row = new ArrayList<>();
            for(int i=0;i<lat.size();i++){
                if( (lat.get(i) >= this.sGrid.getBbox__lat_min()) && (lat.get(i) <= this.sGrid.getDB_bbox__lat_max())){
                    lat_row.add(i);
                }
            }
            for(int i=0;i<lon.size();i++){
                if( (lon.get(i) >= this.sGrid.getBbox__lon_min()) && (lon.get(i) <= this.sGrid.getBbox__lon_max()) ){
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
            field_out = new Double[lat_row.size()][lon_row.size()];
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
                break;
        }
        //Parsing file:
        MyNetCDFParser test = new MyNetCDFParser(filename);
        parseMedOneMinResults out = test.parseMedOneMin();
        if(out == null){
            System.out.println("Parsing fail!");
            return null;
        }
        ArrayList<Double> latTmp = out.getLat();
        ArrayList<Double> lonTmp = out.getLon();
        Double[][] zTmp = (Double[][]) out.getDepth();
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
        return new readout_bathyResults(lat, lon, z);
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
