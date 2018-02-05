package it.uniparthenope;

import java.io.FileReader;
import java.util.*;
import it.uniparthenope.Boxing.*;
import it.uniparthenope.Debug.MyFileWriter;
import it.uniparthenope.Parser.*;

public class JVisirModel {
    //input data
    private long bar_flag;
    private long timedep_flag;
    private int mode;
    private InputPaths paths;

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
    private MyFileWriter logFile;

    //approssimating computation time
    private double seconds;


    //Constructor
    public JVisirModel(long bar_flag, long timedep_flag, InputPaths paths){
        //bar_flag = 1: fresh compution of edges not crossing coastline (mode 1 of GMD-D paper)
        //bar_flag = 2: edges not crossing coastline read out from DB   (mode 2 of GMD-D paper)
        long tic = Utility.Tic();
        if((bar_flag!=1) && (bar_flag!=2)) { //If the input incorrect, set bar_flag = 2 as default value.
            this.bar_flag = 2;
        } else {
            this.bar_flag = bar_flag;
        }
        if((timedep_flag!=0) && (timedep_flag!=2)) { //Same thing for timedep_flag.
            this.timedep_flag = 2;
        } else {
            this.timedep_flag = timedep_flag;
        }
        this.mode=0;
        this.fstats = new Fstats();
        this.constants = new Const();
        this.optim = new Optim();
        this.ship = new Ship();
        this.sGrid = new SpatialGrid();
        this.tGrid = new TemporalGrid();
        this.visualization = new Visualization();

        if(paths==null)
            this.paths = new InputPaths();//default
        else
            this.paths = new InputPaths(paths.getDep_parameters(), paths.getExtr_parameters(),
                    paths.getOptim_parameters(), paths.getSafety_parameters(), paths.getShip_parameters(),
                    paths.getVisualization_parameters(), paths.getOutDir(), paths.getFreeNodesDB(),
                    paths.getCoastlineDB(), paths.getBathymetryDB(), paths.getAnalysisDB(), paths.getForecastFile(), paths.getFreeEdgesDB());
        this.logFile = new MyFileWriter("","",false, this.paths.getOutDir());
        this.logFile.WriteLine("");
        this.logFile.WriteLog("System initialization...");
        this.logFile.CloseFile();
        this.seconds = Utility.nanosecToSec(Utility.Toc(tic));
    }

    public JVisirModel(long bar_flag, long timedep_flag, int mode){//Initialize with standard values defined in settings.m
        this.logFile = new MyFileWriter("","",false, this.paths.getOutDir());
        this.logFile.WriteLine("");
        this.logFile.WriteLog("System initialization...");
        this.logFile.CloseFile();
        long tic = Utility.Tic();

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

        //mode=0: normal execution without serialize objects
        //mode=1: normal execution serializing objects
        //mode=2: debugging mode: load data deserializing objects
        this.mode = mode;
        if(this.mode==0 || this.mode==1){
            if(this.mode==1){
                this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
                this.logFile.WriteLog("\t[mode flag=1]: normal execution with serialization...");
                this.logFile.CloseFile();
            }
            this.fstats = new Fstats();
            this.constants = new Const();
            this.optim = new Optim();
            this.ship = new Ship();
            this.sGrid = new SpatialGrid();
            this.tGrid = new TemporalGrid();
            this.visualization = new Visualization();
        } else if(this.mode==2){
            this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
            this.logFile.WriteLog("\t[DEBUG MODE] loading serialized objects...");
            this.logFile.CloseFile();
            LoadState();
        } else {
            this.mode=0;
            this.fstats = new Fstats();
            this.constants = new Const();
            this.optim = new Optim();
            this.ship = new Ship();
            this.sGrid = new SpatialGrid();
            this.tGrid = new TemporalGrid();
            this.visualization = new Visualization();
        }
        this.seconds = Utility.nanosecToSec(Utility.Toc(tic));
    }

    public void Start(){
        long tic = Utility.Tic();
        LoadData();
        CalculateParameters();
        vessel_ResponseResults vesselResponse;
        Grid_definitionResults gridDefinitionResults;
        Fields_regriddingResults fieldsRegriddingResults;
        Edges_definitionResults edgesDefinitionResults;
        this.seconds += Utility.nanosecToSec(Utility.Toc(tic));
        double stopTime;
        tic = Utility.Tic();
        vesselResponse = vessel_Response();
        this.seconds += Utility.nanosecToSec(Utility.Toc(tic));
        tic = Utility.Tic();
        gridDefinitionResults = Grid_definition();
        stopTime = Utility.nanosecToSec(Utility.Toc(tic));
        seconds += stopTime;
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("Done. Grid definition tooks "+stopTime+" sec.");
        this.logFile.CloseFile();
        tic = Utility.Tic();
        fieldsRegriddingResults = Fields_regridding(gridDefinitionResults.getLat_bathy_Inset(), gridDefinitionResults.getLon_bathy_Inset(), gridDefinitionResults.getBathy_Inset(),
                gridDefinitionResults.getLsm_mask(), vesselResponse.getShip_v_LUT(), vesselResponse.getH_array_m(), gridDefinitionResults.getEstGdtDist());
        stopTime = Utility.nanosecToSec(Utility.Toc(tic));
        seconds += stopTime;
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("Done. Fields regridding tooks "+stopTime+" sec.");
        this.logFile.CloseFile();
        tic = Utility.Tic();
        edgesDefinitionResults = Edges_definition(gridDefinitionResults.getXy(), gridDefinitionResults.getXg_array(), gridDefinitionResults.getYg_array(), fieldsRegriddingResults.getVTDH_Inset(),
                fieldsRegriddingResults.getVTPK_Inset(), fieldsRegriddingResults.getVDIR_Inset(), fieldsRegriddingResults.getWindMAGN_Inset(), fieldsRegriddingResults.getWindDIR_Inset(),
                gridDefinitionResults.getBathy_Inset(), gridDefinitionResults.getJ_mask(), vesselResponse.getShip_v_LUT(), vesselResponse.getH_array_m());
        stopTime = Utility.nanosecToSec(Utility.Toc(tic));
        seconds += stopTime;
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("Done. Edges definition tooks "+stopTime+" sec.");
        this.logFile.CloseFile();

        RouteInfo gdtRoute = Gdt_route(gridDefinitionResults.getXy(), edgesDefinitionResults.getFree_edges(),edgesDefinitionResults.getEdge_lenght(),edgesDefinitionResults.getI_bool(),edgesDefinitionResults.getI_ord(),edgesDefinitionResults.getI_point(),
                this.sGrid.getNode_start(), this.sGrid.getNode_end(), edgesDefinitionResults.getNogo_edges(), edgesDefinitionResults.getWaveHeight_edges(), vesselResponse.getShip_v_LUT(), vesselResponse.getH_array_m());

        RouteInfo requestedRoute;
        if(this.timedep_flag==0){
            requestedRoute = Static_algorithm(gridDefinitionResults.getXy(), edgesDefinitionResults.getFree_edges(),edgesDefinitionResults.getSh_delay(),
                    (int) (this.tGrid.getWave_dep_TS()-1),edgesDefinitionResults.getI_bool(),edgesDefinitionResults.getI_ord(),
                    edgesDefinitionResults.getI_point(),this.sGrid.getNode_start(), this.sGrid.getNode_end());
        } else{
            requestedRoute = Dynamic_algorithm(gridDefinitionResults.getXy(), edgesDefinitionResults.getFree_edges(),edgesDefinitionResults.getSh_delay(),
                    (int) (this.tGrid.getWave_dep_TS()-1),edgesDefinitionResults.getI_bool(),edgesDefinitionResults.getI_ord(),
                    edgesDefinitionResults.getI_point(),this.sGrid.getNode_start(), this.sGrid.getNode_end());
        }

        seconds = seconds + gdtRoute.getComputationTime() + requestedRoute.getComputationTime();
        if(this.mode==1){//Serialize data
            this.SaveState(vesselResponse, gridDefinitionResults, fieldsRegriddingResults, edgesDefinitionResults);
        }
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("Route informations:");
        this.logFile.CloseFile();
        gettingRouteInfo(gdtRoute);
        gettingRouteInfo(requestedRoute);
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("Generating GeoJSON...");
        this.logFile.CloseFile();
        try{
            GeoJsonFormatter.writeGeoJson(gdtRoute, requestedRoute, getRouteCoords(gridDefinitionResults.getXy(), gdtRoute.getPath()),
                    getRouteCoords(gridDefinitionResults.getXy(), requestedRoute.getPath()), this.paths.getOutDir());
        } catch (Exception e){
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("GeoJsonFormatter: "+e.getMessage());
            debug.CloseFile();
        }
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("Done. Total execution time: "+Utility.secondsToMins(seconds)+" Min.");
        this.logFile.CloseFile();
    }

    private RouteInfo Gdt_route(double[][] xy, int[][] free_edges, double[] edge_costs, boolean[] I_bool,
                                        int[] I_ord, int[] I_point, long SID, long FID, ArrayList<Integer> nogo_edges,
                                double[][] waveHeight_edges, double[][] ship_v_LUT, ArrayList<Double> H_array_m){
        long t0 = Utility.Tic();
        //Joint mask- constraint:
        for(Integer element : nogo_edges)
            edge_costs[element] = Double.POSITIVE_INFINITY;

        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("Calculating geodetic route...");
        this.logFile.CloseFile();
        long dt = Utility.Toc(t0);
        long tic = Utility.Tic();
        Dijkstra2DResults geoRoute = Algorithms.Dijkstra(free_edges, edge_costs, I_bool, I_ord, I_point, SID, FID);
        long toc = Utility.Toc(tic);
        t0 = Utility.Tic();
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("Done. Geodetic route found in "+ Utility.nanosecToSec(toc)+" seconds.");
        this.logFile.CloseFile();
        double[] partialTimes = get_NodeLabel(free_edges, edge_costs, waveHeight_edges, geoRoute.getPath(), ship_v_LUT, H_array_m);
        get_Info_at_NodeResults routeInfo  = get_Info_at_Node(xy, partialTimes, geoRoute.getPath());
        return new RouteInfo(geoRoute, Utility.nanosecToSec(toc+dt+Utility.Toc(t0)), partialTimes, routeInfo, 0);
    }

    private RouteInfo Static_algorithm(double[][] xy, int[][] free_edges, double[][] sh_delay, int time_step, boolean[] I_bool,
                                       int[] I_ord, int[] I_point, long SID, long FID){
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("Calculating static route at "+time_step+" time step...");
        this.logFile.CloseFile();
        long tic = Utility.Tic();
        //Ignoring current time step, why?
        Dijkstra2DResults staticRoute = Algorithms.Dijkstra(free_edges, sh_delay, 0, I_bool, I_ord, I_point, SID, FID);
        long toc = Utility.Toc(tic);
        tic = Utility.Tic();
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("Done. Static route found in  "+ Utility.nanosecToSec(toc)+" seconds.");
        this.logFile.CloseFile();
        double[] partial_cost = path_labels(free_edges, sh_delay, 0, staticRoute.getPath());
        get_Info_at_NodeResults routeInfo = get_Info_at_Node(xy, partial_cost, staticRoute.getPath());
        return new RouteInfo(staticRoute, Utility.nanosecToSec(toc+Utility.Toc(tic)), partial_cost, routeInfo, 1);

    }

    private RouteInfo Dynamic_algorithm(double[][] xy, int[][] free_edges, double[][] sh_delay, int time_step, boolean[] I_bool,
                                        int[] I_ord, int[] I_point, long SID, long FID){
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("Calculating dynamic route...");
        this.logFile.CloseFile();
        long tic = Utility.Tic();
        DijkstraTimeResults dynamicRoute = Algorithms.DijkstraTime(free_edges, sh_delay, I_bool, I_ord, I_point, SID, FID);
        long toc = Utility.Toc(tic);
        tic = Utility.Tic();
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("Done. Dynamic route found in  "+ Utility.nanosecToSec(toc)+" seconds.");
        this.logFile.CloseFile();
        get_Info_at_NodeResults routeInfo = get_Info_at_Node(xy, dynamicRoute.getPartial_times(), dynamicRoute.getPath());
        return new RouteInfo(dynamicRoute, Utility.nanosecToSec(toc+Utility.Toc(tic)), routeInfo, 2);
    }

    private double[] path_labels(int[][] edge_list, double[][] edge_weights, int colIDX, LinkedList<Integer> path_nodes){
        double[] Labels = new double[path_nodes.size()];
        for(int i=0; i<Labels.length-1; ++i){
            int edge_no = MemorySaver.getIdx(edge_list, path_nodes.get(i), path_nodes.get(i+1));
            Labels[i+1] = Labels[i] + edge_weights[edge_no][colIDX];
        }
        return Labels;
    }

    private void gettingRouteInfo(RouteInfo gdt){
        //type flag:
        //0 = geodetic
        //1 = static
        //2 = dynamic
        String tempString = "";
        switch(gdt.getType()){
            case 0:{
                tempString = "\tGeodetic nav. distance: "+gdt.getCost()+" NM";
                break;
            }
            case 1: {
                double navDist = gdt.getDr_cum()[gdt.getDr_cum().length-1];
                tempString = "\tOptimal nav. distance: " + navDist + " NM";
                break;
            }
            case 2: {
                double navDist = gdt.getDr_cum()[gdt.getDr_cum().length - 1];
                tempString = "\tOptimal nav. distance: " + navDist + " NM";
                break;
            }
        }
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog(tempString);
        this.logFile.CloseFile();
        double gdt_time = 3600.0*gdt.getPartialTimes()[gdt.getPartialTimes().length-1];
        double gdt_avgV = gdt.getDr_cum()[gdt.getDr_cum().length-1] / gdt.getPartialTimes()[gdt.getPartialTimes().length-1];
        switch (gdt.getType()){
            case 0:{
                tempString = "\tGeodetic nav.time: "+Utility.secs2hms(gdt_time);
                break;
            }
            case 1:{
                tempString = "\tOptimal nav.time: "+Utility.secs2hms(gdt_time);
                break;
            }
            case 2:{
                tempString = "\tOptimal nav.time: "+Utility.secs2hms(gdt_time);
                break;
            }
        }
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog(tempString);
        this.logFile.CloseFile();
        switch (gdt.getType()){
            case 0:{
                tempString = "\tGeodetic average speed: "+gdt_avgV+" kts";
                break;
            }
            case 1:{
                tempString = "\tOptimal average speed: "+gdt_avgV+" kts";
                break;
            }
            case 2:{
                tempString = "\tOptimal average speed: "+gdt_avgV+" kts";
                break;
            }
        }
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog(tempString);
        this.logFile.CloseFile();
    }



    private double[] get_NodeLabel(int[][] free_edges, double[] edge_costs, double[][] waveHeight_edges, LinkedList<Integer> path_vector,
                               double[][] ship_v_LUT, ArrayList<Double> H_array_m){
        // finds Node labels (=cumulated time since departure node) along a
        // specified (topological) path in presence of a time dependent environmental field
        int Nwp = path_vector.size();
        double[] label = new double[Nwp];
        int itime = 0;
        double sh_vel =0;

        double const_vel = Utility.maxNotNaN(ship_v_LUT);

        for(int ip=0; ip<(Nwp-1); ++ip){ //waypoint index
            int edge_no = MemorySaver.getIdx(free_edges, path_vector.get(ip), path_vector.get(ip+1));
            if(this.forcing.getAnalytic() == 1)
                sh_vel = waveHeight_edges[edge_no][0];
            else{
                if(this.ship.getVessType() != this.ship.getSailType()){ //motorboat
                    sh_vel = Utility.interp1(H_array_m, ship_v_LUT ,0, waveHeight_edges[edge_no][itime],"extrap", this.paths.getOutDir());
                } else{ //sailboat
                    sh_vel = const_vel;
                }
            }

            label[ip+1] = label[ip] + edge_costs[edge_no] / sh_vel;

            if(this.timedep_flag == 2) {//Dynamic algorithm
                //itime = Math.min(1+ (int) Math.floor(label[ip+1]/this.tGrid.getDt()) , (int) this.tGrid.getNt());
                itime = Math.min((int) Math.floor(label[ip+1]/this.tGrid.getDt()) , (int) this.tGrid.getNt());
            }
        }
        return label;
    }

    private get_Info_at_NodeResults get_Info_at_Node(double[][] xy, double[] partial_cum_times, LinkedList<Integer> path){
        //computes velocity at waypoints of given path
        int nw = partial_cum_times.length; //n of waypoints (start + end node included)
        double[] Dr_cum = new double[nw]; //optimal edge length
        double[] v_opt = new double[nw]; //optimal ship velocity (average along edge)
        double[] theta_opt = new double[nw]; //optimal ship course
        double[] theta_VCM = new double[nw]; //course in the direction of target

        double[] P_VMC = new double[]{xy[path.getLast()][0], xy[path.getLast()][1]};

        double[] edge_delay = Utility.diff(partial_cum_times);

        for(int iw=1; iw<nw ; ++iw){
            double[] P_a = new double[]{xy[path.get(iw-1)][0], xy[path.get(iw-1)][1]};
            double[] P_b = new double[]{xy[path.get(iw)][0], xy[path.get(iw)][1]};

            double edge_length = Haversine_distance(P_a, P_b);

            Dr_cum[iw] = Dr_cum[iw-1] + edge_length;

            v_opt[iw-1] = edge_length/edge_delay[iw-1];

            double theta_sign = Math.atan2(P_b[0] - P_a[0], P_b[1] - P_a[1]) / this.constants.getDeg2rad();
            double alpha_VMC_sign = Math.atan2(P_VMC[0] - P_a[0], P_VMC[1] - P_a[1]) / this.constants.getDeg2rad();

            if(theta_sign<0)
                theta_opt[iw-1] = 360.0 + theta_sign;
            else
                theta_opt[iw-1] = theta_sign;

            if(alpha_VMC_sign < 0)
                theta_VCM[iw -1] = 360.0 + alpha_VMC_sign;
            else
                theta_VCM[iw -1] = alpha_VMC_sign;
        }

        return new get_Info_at_NodeResults(Dr_cum, v_opt, MemorySaver.addNaNAtTheEnd(edge_delay), theta_opt, theta_VCM);
    }


    private double[][] getRouteCoords(double[][] nodes, LinkedList<Integer> nodeIDX){
        double[][] coords = new double[nodeIDX.size()][2];
        for(int i=0; i<nodeIDX.size(); ++i){
            coords[i][0] = nodes[nodeIDX.get(i)][0];
            coords[i][1] = nodes[nodeIDX.get(i)][1];
        }
        return coords;
    }

    private void SaveState(vessel_ResponseResults vesselResponse, Grid_definitionResults gridDefinitionResults, Fields_regriddingResults fieldsRegriddingResults, Edges_definitionResults edgesDefinitionResults){
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("Saving state...");
        this.logFile.CloseFile();
        try{
            this.fstats.saveState();
            this.optim.saveState();
            this.ship.saveState();
            this.visualization.saveState();
            this.sGrid.saveState();
            this.tGrid.saveState();
            this.extreme_pts.saveState();
            this.dep_datetime.saveState();
            this.safety.saveState();
            this.forcing.saveState();
            vesselResponse.saveState();
            gridDefinitionResults.saveState();
            fieldsRegriddingResults.saveState();
            edgesDefinitionResults.saveState();
        } catch (Exception e){
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("SaveState: "+e.getMessage());
            debug.CloseFile();
            e.printStackTrace();
        }
    }

    private void LoadState(){
        try{
            this.fstats = new Fstats(true);
            this.optim = new Optim(true);
            this.ship = new Ship(true);
            this.visualization = new Visualization(true);
            this.constants = new Const();
            this.sGrid = new SpatialGrid(true);
            this.tGrid = new TemporalGrid(true);
            this.extreme_pts = new ExtremePoints(true);
            this.dep_datetime = new DepartureParameters(true);

        } catch (Exception e){
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("LoadState: "+e.getMessage());
            debug.CloseFile();
            e.printStackTrace();
        }
    }



    private void LoadData(){//Loading data parsing them from json file defined in inputFiles folder
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("Processing namelists...");
        this.logFile.CloseFile();
//        this.extreme_pts = new ExtremePoints();
//        this.dep_datetime = new DepartureParameters();
//        this.ship.LoadVesselParameters();
//        this.safety = new SafetyParameters();
//        this.optim.OptimizationParameters();
//        this.visualization.VisualizationParameters();
        this.extreme_pts = new ExtremePoints(this.paths.getExtr_parameters(), this.paths.getOutDir());
        this.dep_datetime = new DepartureParameters(this.paths.getDep_parameters());
        this.ship.LoadVesselParameters(this.paths.getShip_parameters(), this.paths.getOutDir());
        this.safety = new SafetyParameters(this.paths.getSafety_parameters());
        this.optim.OptimizationParameters(this.paths.getOptim_parameters());
        this.visualization.VisualizationParameters(this.paths.getVisualization_parameters());
    }

    private void CalculateParameters(){//namelist_postproc.m
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
            this.dep_datetime.setHour(6);
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


    private ship_ModelResults ship_Model(){//called from vessel_Response.m, ship_model.m implementation
        //Look-up table for involuntary ship speed reduction in waves.
        this.ship.setFn_max(this.constants.getMs2kts(), this.constants.getG0());
        //LUT independent variables:
        long Nh = 25; //40
        long Nl = 40; //30
        double Hmax = 8.0; //[m]
        //long ih2 = (long) Math.floor(Nh/2);
        //long il2 = (long) Math.floor(Nl/2);

        //Significant wave height
        //Adding to H_array_m (Nh-1) elements between 10^-1 and 10^log10(Hmax)
        ArrayList<Double> H_array_m = Utility.logspace(-1.0, Math.log10(Hmax), Nh-1);
        H_array_m.add(0, 0.0);//adding 0 at first element of H_array_m

        //preallocations:
        double[][] vel_LUT = new double[(int)Nh][this.ship.getP_level_hp().size()];
        //double[][] Rc_LUT =  new double[(int)Nh][this.ship.getP_level_hp().size()];
        //double[][] Raw_LUT =  new double[(int)Nh][this.ship.getP_level_hp().size()];
        //ArrayList<Double> P_level_thro = new ArrayList<Double>();
        double max = Collections.max(this.ship.getP_level_hp());
//        for(double element : this.ship.getP_level_hp()){
//            P_level_thro.add((100*element) / max);
//        }

        //pars for v_Bowditch:
//        double m2ft = 3.2808399;
//        double a3_ref = 0.0248; //[kts/ ft^2]
//        double a2_ref = 0.0165; //[kts/ ft^2]
//        double a1_ref = 0.0083; //[kts/ ft^2]
//        double v0_ref = 18.0; //kts
        //LUT computation
        for(int ih = 0; ih<Nh; ++ih){
            this.ship.ship_resistance(H_array_m.get(ih), this.constants);
            for(int j = 0; j<this.ship.getP_level_hp().size(); ++j){
                vel_LUT[ih][j] = this.ship.getV_out_kts().get(j);
                //Rc_LUT[ih][j] = this.ship.getR_c().get(j);
                //Raw_LUT[ih][j] = this.ship.getR_aw().get(j);
            }
        }
        this.ship.ship_resistance(0.0,this.constants);
        //ArrayList<Double> v0 = this.ship.getV_out_kts();
        //ArrayList<Double> Rc0 = this.ship.getR_c();
        //ArrayList<Double> Raw0 = this.ship.getR_aw();

        return new ship_ModelResults(vel_LUT, H_array_m);
    }



    private vessel_ResponseResults vessel_Response(){//vessel_Response.m implementation
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        if(this.forcing.getAnalytic()==1){
            this.logFile.WriteLog("Analitical benchmark...");
        } else{
            this.logFile.WriteLog("Ship model response (motorboat)...");

        }
        this.logFile.CloseFile();
        this.ship.shipPowerLevels();
        ship_ModelResults tmp = this.ship_Model();
        this.ship.setMaxWind(Double.NaN);
        this.ship.setMinWind(Double.NaN);
        return new vessel_ResponseResults(tmp.getShip_v_LUT(), tmp.getH_array_m(), new ArrayList<Double>(), new ArrayList<Double>());
    }

    private boolean pointInPolygon(double x, double y, Polygon polygon){
        boolean isInside = false;
        int nVert = polygon.getCornersNumber();
        int i, j;
        for(i=0, j=nVert-1; i<nVert; j = i++){
            if(((polygon.getCorner(i).getY() >= y) != (polygon.getCorner(j).getY() >= y)) &&
                    (x <= (polygon.getCorner(j).getX()-polygon.getCorner(i).getX()) * (y-polygon.getCorner(i).getY()) /
                            (polygon.getCorner(j).getY()-polygon.getCorner(i).getY()) + polygon.getCorner(i).getX() ) ){
                isInside = !isInside;
            }
        }

        return isInside;
    }

    private Grid_definitionResults Grid_definition(){//Grid_definition.m implementation
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("Bounding boxes and bathymetry postprocessing (grid definition): ");
        this.logFile.CloseFile();
        Grid_definitionResults GridOut = new Grid_definitionResults();
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
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("\tLoading coastline data...");
        this.logFile.CloseFile();
        //Coastline:
        readout_coastResults readoutCoastResults = this.readout_coast();
        ArrayList<Double> y_coast = new ArrayList<>();
        ArrayList<Double> x_coast = new ArrayList<>();
        y_coast.addAll(readoutCoastResults.getLat_int());
        y_coast.add(Double.NaN);
        y_coast.addAll(readoutCoastResults.getLat_ext());
        x_coast.addAll(readoutCoastResults.getLon_int());
        x_coast.add(Double.NaN);
        x_coast.addAll(readoutCoastResults.getLon_ext());
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("\tCreating target grid and loading bathymetry data...");
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
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("\tCalculating target grid and reduced bathy field...");
        this.logFile.CloseFile();
        //Target grid and reduced bathy field:
        //grid_extreme_coords
        grid_extreme_coordsResults insets = this.grid_extreme_coords(lat_bathy, lon_bathy,z_bathy);
        ArrayList<Double> lat_bathy_Inset = insets.getLat_red();
        ArrayList<Double> lon_bathy_Inset = insets.getLon_red();
        double[][] bathy_Inset = insets.getField_out();

        ArrayList<Boolean> x_bool=new ArrayList<>();
        ArrayList<Boolean> y_bool=new ArrayList<>();
        ArrayList<Boolean> inset_bool=new ArrayList<>();
        ArrayList<Double> x_coast_Inset = new ArrayList<>();
        ArrayList<Double> y_coast_Inset = new ArrayList<>();
        if(this.bar_flag == 2){
            // coastline excerpt within Inset:
            this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
            this.logFile.WriteLog("\tCalculating coastline excerpt within inset...");
            this.logFile.CloseFile();
            double min_lon_bathy = Collections.min(lon_bathy_Inset);
            double max_lon_bathy = Collections.max(lon_bathy_Inset);
            double min_lat_bathy = Collections.min(lat_bathy_Inset);
            double max_lat_bathy = Collections.max(lat_bathy_Inset);
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
        double[][] xy_DB = data_gridOut.getXy();
       // double[][] xg_DB = data_gridOut.getXg();
        //double[][] yg_DB = data_gridOut.getYg();

        //inset grid plaid coordinates:
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("\tinset grid plaid coordinates...");
        this.logFile.CloseFile();
        mdata_gridResults data_gridOut2 = mdata_grid(lat_bathy_Inset, lon_bathy_Inset);
        GridOut.setXy(data_gridOut2.getXy());
        GridOut.setXg(data_gridOut2.getXg());
        GridOut.setYg(data_gridOut2.getYg());
        if(this.visualization.getGraphData()==1){
            //csv_write xy
            this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
            this.logFile.WriteLog("\tWriting graph coords in GRAPH.node_LonLat.csv...");
            this.logFile.CloseFile();
            try{
                MyCSVParser csv = new MyCSVParser(this.paths.getOutDir()+"GRAPH.node_LonLat.csv");
                csv.writeCSV(GridOut.getXy());
            } catch (Exception e){
                MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
                debug.WriteLog("grid_definition, csv parser: "+e.getMessage());
                debug.CloseFile();
                e.printStackTrace();
            }
        }

        //--------------------------------------------------------------------
        //Joint coast-vessel safety mask:
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
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
        //double[][] J_bathy_Inset;//VALORE MAI USATO
        double[] xg_array;
        double[] yg_array;
        double[][] xy_g;
        double[][] dist_mask;
        double[][] min_coast_dist;
        if(this.bar_flag == 2){
            //readout free nodes DB:
            this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
            this.logFile.WriteLog("\treadout freenodes DB...");
            this.logFile.CloseFile();
            MyBinaryParser datFile = new MyBinaryParser(this.paths.getFreeNodesDB());
            long[] free_nodes_DB = datFile.readAsUInt32();
            //remapping free nodes:
            this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
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
                MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
                debug.WriteLog("Grid_definition: no free nodes found in graph");
                debug.CloseFile();
                System.exit(0);
            } else if(free_nodes_Number > this.sGrid.getNodesLargeN()){
                System.out.println("too large graph");
                MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
                debug.WriteLog("Grid_definition: too large graph");
                debug.CloseFile();
                System.exit(0);
            }
            this.sGrid.setFreenodes(free_nodes_Number);
            //lsm, created using coastline DB (NaNs on landmass):
            this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
            this.logFile.WriteLog("\tCreating lsm_mask using coastline DB...");
            this.logFile.CloseFile();
            double[][] lsm_mask = Utility.NaNmatrix(GridOut.getXg().length,GridOut.getXg()[0].length);
            int nRows = GridOut.getXg().length;
            int ii;
            int jj;
            for(int k=0;k<free_nodes.length;++k){
                int index =(int) free_nodes[k]-1;
                jj=index/nRows;
                ii=(index%nRows);
                lsm_mask[ii][jj]=1.0;
            }
            int[] dim = new int[2];
            dim[0]=GridOut.getXg().length;
            dim[1]=GridOut.getXg()[0].length;

            //Safe distance from coastline:
            this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
            this.logFile.WriteLog("\tCalculating safe distance from coastline...");
            this.logFile.CloseFile();
            xg_array = Utility.reshape(GridOut.getXg(),dim[0]*dim[1], this.paths.getOutDir());
            yg_array = Utility.reshape(GridOut.getYg(), GridOut.getYg().length*GridOut.getYg()[0].length, this.paths.getOutDir());
            xy_g = new double[xg_array.length][2];
            for(int i=0;i<xg_array.length;i++){
                xy_g[i][0] = xg_array[i];
                xy_g[i][1] = yg_array[i];
            }
            int nC = x_coast_Inset.size();
            int cols = dim[0]*dim[1];
            if(nC>0){
                double[][] coast_dist = new double[nC][cols];
                double[][] P_b = new double[1][2];
                for(int i=0;i<nC;++i){
                    P_b[0][0]=x_coast_Inset.get(i);
                    P_b[0][1]=y_coast_Inset.get(i);
                    coast_dist[i] = Haversine_distanceOneP_b(xy_g, P_b);
//                    double[] hor_dist = Haversine_distanceOneP_b(xy_g, P_b);
//                    for(int j=0; j<cols;j++){
//                        coast_dist[i][j]= hor_dist[j];
//                    }
                }
                double[] min_coast_distTmp = Utility.min(coast_dist,1);
                min_coast_dist = Utility.reshape(min_coast_distTmp, dim, this.paths.getOutDir());
                dist_mask = Utility.NaNmatrix(dim[0],dim[1]);
                for(int i=0;i<dist_mask.length;++i){
                    for(int j=0;j<dist_mask[0].length;++j){
                        if(min_coast_dist[i][j]>=this.extreme_pts.getMinCoastDist()){
                            dist_mask[i][j]=1.0;
                        }
                    }
                }

            } else{
                dist_mask = Utility.ones(dim[0],dim[1]);
            }
            this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
            this.logFile.WriteLog("\tCalculating Joint safe mask...");
            this.logFile.CloseFile();
            // %###############################################
            // %
            // %  Joint safe mask:
            // %
            // J_mask = lsm_mask' .* UKC_mask .* dist_mask' ;
            // %
            // %###############################################
            double[][] J_mask = Utility.MatrixComponentXcomponent(Utility.MatrixComponentXcomponent(Utility.transposeMatrix(lsm_mask),UKC_mask, this.paths.getOutDir()),Utility.transposeMatrix(dist_mask), this.paths.getOutDir());

            xg_Jmasked = Utility.MatrixComponentXcomponent(GridOut.getXg(),Utility.transposeMatrix(J_mask), this.paths.getOutDir());
            yg_Jmasked = Utility.MatrixComponentXcomponent(GridOut.getYg(),Utility.transposeMatrix(J_mask), this.paths.getOutDir());
            //J_bathy_Inset = Utility.MatrixComponentXcomponent(bathy_Inset,J_mask);//VALORE MAI USATO

            // % Note: It is not necessary to pass [xg_Jmasked,yg_Jmasked] to the MAIN.m in place of [xg,yg].
            // % Indeed, shortest path search already accounts (s. Edges_definition.m)
            // % for coastline (free_edges) and bathymetry/lsm/dist-from-coastline (nogo_edges)
            // % [xg_Jmasked,yg_Jmasked] are used here just for the sake of searching departure and arrival nodes (s. below).
            //-----------------------------------------------------------------------------------------------------
            //Start and end nodes:
            xg_array = Utility.reshape(xg_Jmasked,xg_Jmasked.length*xg_Jmasked[0].length, this.paths.getOutDir());
            yg_array = Utility.reshape(yg_Jmasked, yg_Jmasked.length*yg_Jmasked[0].length, this.paths.getOutDir());
            xy_g = new double[xg_Jmasked.length*xg_Jmasked[0].length][2];
            for(int i =0 ;i<xg_Jmasked.length*xg_Jmasked[0].length;++i){
                xy_g[i][0]=xg_array[i];
                xy_g[i][1]=yg_array[i];
            }
            this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
            this.logFile.WriteLog("\tCalculating distances from start/end nodes...");
            this.logFile.CloseFile();
            //distance from start/end nodes:
            double[][] tmpP_b=new double[1][2];
            tmpP_b[0][0]=this.extreme_pts.getStart_lon();
            tmpP_b[0][1]=this.extreme_pts.getStart_lat();
            //double[] start_dist_matrix = hor_distance("s",xy_g,tmpP_b);
            double[] start_dist_matrix = Haversine_distanceOneP_b(xy_g, tmpP_b);
            tmpP_b[0][0]=this.extreme_pts.getEnd_lon();
            tmpP_b[0][1]=this.extreme_pts.getEnd_lat();
            //double[] end_dist_matrix = hor_distance("s",xy_g,tmpP_b);
            double[] end_dist_matrix = Haversine_distanceOneP_b(xy_g,tmpP_b);

            this.sGrid.setMin_start_dist(Utility.min(start_dist_matrix));
            this.sGrid.setMin_end_dist(Utility.min(end_dist_matrix));

            if(Double.isNaN(this.sGrid.getMin_start_dist()) || Double.isNaN(this.sGrid.getMin_end_dist())){
                System.out.println("departure or arrival point not compliant with safety specifications");
                MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
                debug.WriteLog("Grid_definition: departure or arrival point not compliant with safety specifications");
                debug.CloseFile();
                System.exit(0);
            }
            ArrayList<Integer> tmpIndexes = new ArrayList<>();
            for(int i=0;i<start_dist_matrix.length;++i){
                if(start_dist_matrix[i]==this.sGrid.getMin_start_dist()){
                    tmpIndexes.add(i);
                }
            }
            this.sGrid.setNode_start(Collections.min(tmpIndexes));
            tmpIndexes = new ArrayList<>();
            for(int i=0;i<end_dist_matrix.length;++i){
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

            double estGdtDist = Haversine_distance(p_a, p_b);

            p_a[0] = this.sGrid.getNode_end_lon();
            p_a[1] = this.sGrid.getNode_end_lat();
            p_b[0] = this.sGrid.getNode_end_lon();
            p_b[1] = this.sGrid.getNode_start_lat();
            double delta_y = Haversine_distance(p_a, p_b);
            int sign = Utility.sign(this.sGrid.getNode_end_lat()-this.sGrid.getNode_start_lat());
            delta_y=delta_y*sign;

            p_a[0] = this.sGrid.getNode_end_lon();
            p_a[1] = this.sGrid.getNode_start_lat();
            p_b[0] = this.sGrid.getNode_start_lon();
            p_b[1] = this.sGrid.getNode_start_lat();
            double delta_x = Haversine_distance(p_a, p_b);
            sign = Utility.sign(this.sGrid.getNode_end_lon()-this.sGrid.getNode_start_lon());
            delta_x=delta_x*sign;

            //orientation of the gdt route
            this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
            this.logFile.WriteLog("\tCalculating orientation of the geodetic route...");
            this.logFile.CloseFile();
            double atan2 = Math.atan2(delta_x,delta_y);
            this.sGrid.setTheta_gdt(atan2);
            //th= Sgrid.theta_gdt/const.deg2rad


            /**NEW*/
            ArrayList<Polygon> vFences;
            try{
                vFences = GeoJsonFormatter.getPolygons("inputFiles/Fences/virtualFences.geojson");
                this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
                this.logFile.WriteLog("\tGetting virtual fences info from DB...");
                this.logFile.CloseFile();
            } catch (Exception e){
                this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
                this.logFile.WriteLog("\tFences Info not found, skipping...");
                this.logFile.CloseFile();
                vFences = new ArrayList<>();
            }
            /**NEW*/
            if(vFences.size()>0){
                this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
                this.logFile.WriteLog("\tRemoving virtual fences from grid...");
                this.logFile.CloseFile();
                for(int i=0;i<lsm_mask.length;++i){
                    for(int j=0;j<lsm_mask[0].length;++j){
                        for(Polygon polygon : vFences){
                            if(pointInPolygon(lon_bathy_Inset.get(i), lat_bathy_Inset.get(j), polygon)){
                                lsm_mask[i][j] = Double.NaN;
                            }
                        }
                    }
                }
            }

            GridOut.setXg_array(xg_array);
            GridOut.setYg_array(yg_array);
            GridOut.setXy_DB(xy_DB);
            GridOut.setLat_bathy_Inset(lat_bathy_Inset);
            GridOut.setLon_bathy_Inset(lon_bathy_Inset);
            GridOut.setBathy_Inset(bathy_Inset);
            GridOut.setLsm_mask(lsm_mask);
            GridOut.setJ_mask(J_mask);
            GridOut.setEstGdtDist(estGdtDist);

        } else if(this.bar_flag==1){
            double[][] lsm_mask = Utility.NaNmatrix(GridOut.getXg().length,GridOut.getXg()[0].length);
            double[][] tmp = Utility.transposeMatrix(lsm_mask);
            double[][] J_mask = Utility.NaNmatrix(tmp.length,tmp[0].length);
            //-----------------------------------------------------------------------------------------------------
            //Start and end nodes:
            xg_array = new double[0];
            yg_array = new double[0];
            double estGdtDist = Double.NaN;
            GridOut.setXg_array(xg_array);
            GridOut.setYg_array(yg_array);
            GridOut.setXy_DB(xy_DB);
            GridOut.setLat_bathy_Inset(lat_bathy_Inset);
            GridOut.setLon_bathy_Inset(lon_bathy_Inset);
            GridOut.setBathy_Inset(bathy_Inset);
            GridOut.setLsm_mask(lsm_mask);
            GridOut.setJ_mask(J_mask);
            GridOut.setEstGdtDist(estGdtDist);
        }
        GridOut.setX_islands(readoutCoastResults.getLon_int());
        GridOut.setY_islands(readoutCoastResults.getLat_int());
        GridOut.setX_continent(readoutCoastResults.getLon_ext());
        GridOut.setY_continent(readoutCoastResults.getLat_ext());
        return GridOut;
    }

    private double changeDirRule(double inField){
        double out;
        double tmp = inField;
        if(tmp<0)
            tmp+=(2*Math.PI);
        tmp -= (Math.PI/2);
        if(tmp < 0)
            tmp+=(2*Math.PI);
        out = -tmp + (2*Math.PI);
        return out;
    }

//    private double changeDirRule(double inField){
//        // % change of wind/current directional convention:
//        // %
//        // % from atan2 output concention to WAM-like convention, i.e.:
//        // % from:  [-pi, pi] counterclockwise with 0 at  3:00 o'clock
//        // %   to:  [0, 2*pi]        clockwise with 0 at 12:00 o'clock
//        // %
//        // % (inField and outField must be in radians)
//        // %
//        //---------------------------------------------
//        double inField1 = inField;
//        double inField1_ = 0;
//        if(inField1 < 0){
//            inField1_ = inField1;
//            inField1_ = inField1_ + 2*Math.PI;
//            inField1 = inField1_;
//        }
//        //---------------------------------------------
//        double inField2 = inField1 -Math.PI/2;
//        double inField2_ = 0;
//        if(inField2 < 0){
//            inField2_ = inField2;
//            inField2_ = inField2_ + 2*Math.PI;
//            inField2 = inField2_;
//        }
//        //---------------------------------------------
//        return (-inField2+(2*Math.PI));
//    }

//    private double[] changeDirRule(double[] inField){
//        // % change of wind/current directional convention:
//        // %
//        // % from atan2 output concention to WAM-like convention, i.e.:
//        // % from:  [-pi, pi] counterclockwise with 0 at  3:00 o'clock
//        // %   to:  [0, 2*pi]        clockwise with 0 at 12:00 o'clock
//        // %
//        // % (inField and outField must be in radians)
//        // %
//        //---------------------------------------------
//        double[] inField1 = Utility.deepCopy(inField);
//        for(int i=0;i<inField1.length;i++){
//            if(inField1[i]<0){
//                inField1[i]+=(2*Math.PI);
//            }
//        }
//
//        //---------------------------------------------
//        double[] inField2 = inField1;
//        for(int i=0;i<inField2.length;i++){
//            inField2[i]-=(Math.PI/2);
//            if(inField2[i]<0){
//                inField2[i]+=(2*Math.PI);
//            }
//        }
//        //---------------------------------------------
//        double[] outFiled = new double[inField2.length];
//        for(int i=0;i<inField2.length;i++){
//            outFiled[i] = -inField2[i]+(2*Math.PI);
//        }
//        return outFiled;
//    }

    private double[][] changeDirRule(double[][] inField){
        double[][] out = new double[inField.length][inField[0].length];
            for(int i=0;i<inField.length;++i){
                for(int j=0;j<inField[0].length;++j){
                    double tmp = inField[i][j];
                    if(tmp<0)
                        tmp+=(2*Math.PI);
                    tmp -= (Math.PI/2);
                    if(tmp < 0)
                        tmp+=(2*Math.PI);

                    out[i][j] = -tmp + (2*Math.PI);
                }
            }
        return out;
    }

//    private double[][] changeDirRule(double[][] inField){
//        // % change of wind/current directional convention:
//        // %
//        // % from atan2 output concention to WAM-like convention, i.e.:
//        // % from:  [-pi, pi] counterclockwise with 0 at  3:00 o'clock
//        // %   to:  [0, 2*pi]        clockwise with 0 at 12:00 o'clock
//        // %
//        // % (inField and outField must be in radians)
//        // %
//        //---------------------------------------------
//        double[][] inField1 = Utility.deepCopy(inField);
//        for(int i=0;i<inField.length;i++){
//            for(int j=0;j<inField[0].length;j++){
//                if(inField1[i][j]<0){
//                    inField1[i][j]+=(2*Math.PI);
//                }
//            }
//        }
//        //---------------------------------------------
//        double[][] inField2 = Utility.deepCopy(inField1);
//        for(int i=0;i<inField2.length;i++){
//            for(int j=0;j<inField2[0].length;j++){
//                inField2[i][j] -=(Math.PI/2);
//                if(inField2[i][j]<0){
//                    inField2[i][j] += (2*Math.PI);
//                }
//            }
//        }
//        //---------------------------------------------
//        double[][] outField = new double[inField2.length][inField2[0].length];
//        for(int i=0;i<inField.length;i++){
//            for(int j=0;j<inField[0].length;j++){
//                outField[i][j] = -inField2[i][j] + (2*Math.PI);
//            }
//        }
//        return outField;
//    }

    private double[][][] changeDirRule(double[][][] inField){
        double[][][] out = new double[inField.length][inField[0].length][inField[0][0].length];
        for(int k=0;k<inField.length;++k){
            for(int i=0;i<inField[0].length;++i){
                for(int j=0;j<inField[0][0].length;++j){
                    double tmp = inField[k][i][j];
                    if(tmp<0)
                        tmp+=(2*Math.PI);
                    tmp -= (Math.PI/2);
                    if(tmp < 0)
                        tmp+=(2*Math.PI);

                    out[k][i][j] = -tmp + (2*Math.PI);
                }
            }
        }
        return out;
    }

//    private double[][][] changeDirRule(double[][][] inField){
//        // % change of wind/current directional convention:
//        // %
//        // % from atan2 output concention to WAM-like convention, i.e.:
//        // % from:  [-pi, pi] counterclockwise with 0 at  3:00 o'clock
//        // %   to:  [0, 2*pi]        clockwise with 0 at 12:00 o'clock
//        // %
//        // % (inField and outField must be in radians)
//        // %
//        //---------------------------------------------
//        double[][][] inField1 = Utility.deepCopy(inField);
//        for(int i=0;i<inField.length;i++){
//            for(int j=0;j<inField[0].length;j++){
//                for(int k=0;k<inField[0][0].length;k++){
//                    if(inField1[i][j][k]<0){
//                        inField1[i][j][k]+=(2*Math.PI);
//                    }
//                }
//            }
//        }
//        //---------------------------------------------
//        double[][][] inField2 = Utility.deepCopy(inField1);
//        for(int i=0;i<inField2.length;i++){
//            for(int j=0;j<inField2[0].length;j++){
//                for(int k=0;k<inField2[0][0].length;k++){
//                    inField2[i][j][k] -=(Math.PI/2);
//                    if(inField2[i][j][k]<0){
//                        inField2[i][j][k] += (2*Math.PI);
//                    }
//                }
//            }
//        }
//        //---------------------------------------------
//        double[][][] outField = new double[inField2.length][inField2[0].length][inField2[0][0].length];
//        for(int i=0;i<inField.length;i++){
//            for(int j=0;j<inField[0].length;j++){
//                for(int k=0;k<inField[0][0].length;k++){
//                    outField[i][j][k] = -inField2[i][j][k] + (2*Math.PI);
//                }
//            }
//        }
//        return outField;
//    }

    public double[] Haversine_distance(double[][] xy, int[][] free_edges){

        this.constants.setDeg2rad(Math.PI/180);
        double E_radius = 3444.0; //Earth radius in NM
        double[] dd = new double[free_edges.length];
        for(int i=0;i<free_edges.length;i++){
            double a = Math.sqrt(Math.pow( Math.sin(this.constants.getDeg2rad()*(xy[ free_edges[i][1] ][1] - xy[ free_edges[i][0] ][1])/2) , 2) + Math.cos(this.constants.getDeg2rad()*xy[ free_edges[i][0] ][1])
                    * Math.cos(this.constants.getDeg2rad()*xy[ free_edges[i][1] ][1]) * Math.pow( Math.sin( this.constants.getDeg2rad()*(xy[ free_edges[i][1] ][0]-xy[ free_edges[i][0] ][0])/2 ) , 2));
            dd[i] = E_radius * 2 * Math.asin(a);
        }
        return dd;
    }

//    private double[] hor_distanceForEdges(String method, double[][]xy, int[][]free_edges, double... varargin){
//        /****
//         //              P_a(:,0)                    P_b(:,0)
//         xm[i] = (xy[ free_edges[i][0] ][0] + xy[ free_edges[i][1] ][0])/2;
//         //              P_a(:,1)                    P_b(:,1)
//         ym[i] = (xy[ free_edges[i][0] ][1] + xy[ free_edges[i][1] ][1])/2;
//         ***/
//        this.constants.setDeg2rad(Math.PI/180);
//        double grid_step_in_NM = varargin.length > 0 ? varargin[0] : 1.0;
//        //P_a = P_1
//        //P_b = P_2
//        double[] dd = new double[free_edges.length];
//        if(method=="plane"||method=="p"){
//            for(int i=0;i<free_edges.length;i++){
//                dd[i]=grid_step_in_NM*Math.sqrt(Math.pow(xy[ free_edges[i][0] ][0] - xy[ free_edges[0][1] ][0], 2) +
//                        Math.pow(xy[ free_edges[i][0] ][1] - xy[ free_edges[0][1] ][1], 2));
//            }
//        } else{
//            if(method=="sphere" || method=="s"){
//                double E_radius = 3444.0; //NM
//                for(int i=0;i<free_edges.length;i++){
//                    /*
//                    * dd1[i]=E_radius * Math.acos(Math.cos(this.constants.getDeg2rad()*P_a[i][1])*Math.cos(this.constants.getDeg2rad()*P_b[0][1])*
//                        Math.cos(this.constants.getDeg2rad()*(P_a[i][0]-P_b[0][0]))+Math.sin(this.constants.getDeg2rad()*P_a[i][1])*
//                        Math.sin(this.constants.getDeg2rad()*P_b[0][1]));
//                    * */
//                    dd[i]=E_radius *
//                            Math.acos( Math.cos(this.constants.getDeg2rad()*xy[ free_edges[i][0] ][1]) *
//                            Math.cos(this.constants.getDeg2rad()*xy[ free_edges[0][1] ][1]) *
//                            Math.cos(this.constants.getDeg2rad()*(xy[ free_edges[i][0] ][0]-xy[ free_edges[0][1] ][0])) +
//                            Math.sin(this.constants.getDeg2rad()*xy[ free_edges[i][0] ][1]) *
//                            Math.sin(this.constants.getDeg2rad()*xy[ free_edges[0][1] ][1]) );
//                }
//            }
//        }
//        return dd;
//    }

    public double[] Haversine_distanceOneP_b(double[][] P_a, double[][] P_b){//P_B è un solo punto, per questo non scorro

        this.constants.setDeg2rad(Math.PI/180);
        double E_radius = 3444.0; //Earth radius in NM
        double[] dd=new double[P_a.length];
        for(int i=0;i<P_a.length;++i){
            double a = Math.sqrt(Math.pow( Math.sin(this.constants.getDeg2rad()*(P_b[0][1] - P_a[i][1])/2) , 2) + Math.cos(this.constants.getDeg2rad()*P_a[i][1])
                    * Math.cos(this.constants.getDeg2rad()*P_b[0][1]) * Math.pow( Math.sin( this.constants.getDeg2rad()*(P_b[0][0]-P_a[i][0])/2 ) , 2));
            dd[i] = E_radius * 2 * Math.asin(a);
        }
        return dd;
    }

//    public double[] Haversine_distance(double[][] P_a, double[][] P_b){
//
//        this.constants.setDeg2rad(Math.PI/180);
//        double E_radius = 3444.0; //Earth radius in NM
//        double[] dd=new double[P_a.length];
//        for(int i=0;i<P_a.length;++i){
//            double a = Math.sqrt(Math.pow( Math.sin(this.constants.getDeg2rad()*(P_b[i][1] - P_a[i][1])/2) , 2) + Math.cos(this.constants.getDeg2rad()*P_a[i][1])
//                    * Math.cos(this.constants.getDeg2rad()*P_b[i][1]) * Math.pow( Math.sin( this.constants.getDeg2rad()*(P_b[i][0]-P_a[i][0])/2 ) , 2));
//            dd[i] = E_radius * 2 * Math.asin(a);
//        }
//        return dd;
//    }

//    private double[] hor_distance(String method, double[][] P_a, double[][] P_b, double... varargin){
//        // % horizontal distance between pair of points (expressed in NM)
//        // % either on the plane or the sphere
//        // % works also with arrays
//        // % optionally rescales distances by grid_step_in_NM
//        // %
//        // % P_a and P_b must be linear array of the same size
//        // % (build via meshgrid and then reshape)
//        // %
//        // % P_a: [lon, lat]
//        // % P_b: [lon, lat]
//        // %
//        // % method='p' for planar geometry
//        // % method='s' for spherical geometry
//        this.constants.setDeg2rad(Math.PI/180);
//        double grid_step_in_NM = varargin.length > 0 ? varargin[0] : 1.0;
//
////        double[] xa = new double[P_a.length];
////        double[] ya = new double[P_a.length];
////        for(int i=0;i<P_a.length;i++){
////            xa[i]=P_a[i][0];
////            ya[i]=P_a[i][1];
////        }
////
////        double[] xb = new double[P_b.length];
////        double[] yb = new double[P_b.length];
////        for(int i=0;i<P_b.length;i++){
////            xb[i]=P_b[i][0];
////            yb[i]=P_b[i][1];
////        }
//
//        double[] dd1=new double[P_a.length];
//        if(method=="plane"||method=="p"){
//            for(int i=0;i<P_a.length;i++){
//                dd1[i]=grid_step_in_NM*Math.sqrt(Math.pow((P_a[i][0]-P_b[0][0]),2) + Math.pow((P_a[i][1]-P_b[0][1]),2));
//            }
//        } else if(method=="sphere" || method=="s"){
//            // % from : http://mathworld.wolfram.com/GreatCircle.html
//            // % x and y must be in degree.
//            // % output in Nautical Miles (NM)
//            double E_radius = 3444.0; //NM
//            for(int i=0;i<P_a.length;i++){
//                dd1[i]=E_radius * Math.acos( Math.cos( this.constants.getDeg2rad()*P_a[i][1])*Math.cos(this.constants.getDeg2rad()*P_b[0][1])*
//                        Math.cos(this.constants.getDeg2rad()*(P_a[i][0]-P_b[0][0]))+Math.sin(this.constants.getDeg2rad()*P_a[i][1])*
//                        Math.sin(this.constants.getDeg2rad()*P_b[0][1]));
//            }
//        }
//
////        double[] dd=new double[P_a.length];
////        if(method=="plane"||method=="p"){
////            for(int i=0;i<P_a.length;i++){
////                dd[i]=grid_step_in_NM*Math.sqrt(Math.pow((xa[i]-xb[0]),2) + Math.pow((ya[i]-yb[0]),2));
////            }
////        } else if(method=="sphere" || method=="s"){
////            // % from : http://mathworld.wolfram.com/GreatCircle.html
////            // % x and y must be in degree.
////            // % output in Nautical Miles (NM)
////            double E_radius = 3444.0; //NM
////            for(int i=0;i<P_a.length;i++){
////                dd[i]=E_radius * Math.acos(Math.cos(this.constants.getDeg2rad()*ya[i])*Math.cos(this.constants.getDeg2rad()*yb[0])*
////                        Math.cos(this.constants.getDeg2rad()*(xa[i]-xb[0]))+Math.sin(this.constants.getDeg2rad()*ya[i])*
////                        Math.sin(this.constants.getDeg2rad()*yb[0]));
////            }
////        }
//
//        return dd1;
//    }

    public double Haversine_distance(double[] P_a, double[] P_b){
        this.constants.setDeg2rad(Math.PI/180);

        double E_radius = 3444.0; //Earth radius in NM
        double a = Math.sqrt(Math.pow( Math.sin(this.constants.getDeg2rad()*(P_b[1] - P_a[1])/2) , 2) + Math.cos(this.constants.getDeg2rad()*P_a[1])
                * Math.cos(this.constants.getDeg2rad()*P_b[1]) * Math.pow( Math.sin( this.constants.getDeg2rad()*(P_b[0]-P_a[0])/2 ) , 2));
        return E_radius * 2 * Math.asin(a);
    }

//    private double hor_distance(String method, double[] P_a, double[] P_b, double... varargin){
//        // % horizontal distance between pair of points (expressed in NM)
//        // % either on the plane or the sphere
//        // % works also with arrays
//        // % optionally rescales distances by grid_step_in_NM
//        // %
//        // % P_a and P_b must be linear array of the same size
//        // % (build via meshgrid and then reshape)
//        // %
//        // % P_a: [lon, lat]
//        // % P_b: [lon, lat]
//        // %
//        // % method='p' for planar geometry
//        // % method='s' for spherical geometry
//        this.constants.setDeg2rad(Math.PI/180);
//        double grid_step_in_NM = varargin.length > 0 ? varargin[0] : 1.0;
//
//        double xa = P_a[0];
//        double ya = P_a[1];
//
//        double xb = P_b[0];
//        double yb = P_b[1];
//
//
//        double dd=0;
//        if(method=="plane"||method=="p"){
//            dd=grid_step_in_NM*Math.sqrt(Math.pow((xa-xb),2) + Math.pow((ya-yb),2));
//        } else if(method=="sphere" || method=="s"){
//            // % from : http://mathworld.wolfram.com/GreatCircle.html
//            // % x and y must be in degree.
//            // % output in Nautical Miles (NM)
//            double E_radius = 3444.0; //NM
//            dd=E_radius * Math.acos(Math.cos(this.constants.getDeg2rad()*ya)*Math.cos(this.constants.getDeg2rad()*yb)*
//                    Math.cos(this.constants.getDeg2rad()*(xa-xb))+Math.sin(this.constants.getDeg2rad()*ya)*
//                    Math.sin(this.constants.getDeg2rad()*yb));
//        }
//
//        return dd;
//    }

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
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("idx_ref2inset_grid: inset grid not within reference grid!");
            debug.CloseFile();
            System.exit(0);
        }
        if((Utility.any(idx_big,"<",1)) || Utility.any(idx_big,">",(nx_big*ny_big))){
            System.out.println("grid index not within reference grid!");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("idx_ref2inset_grid: grid index not within reference grid!");
            debug.CloseFile();
            System.exit(0);
        }

        long[] col = new long[idx_big.length];
        long[] row = new long[idx_big.length];
        long[] idx = new long[idx_big.length];
        for(int i=0;i<idx_big.length;++i){
            //col idx (ref):
            long col_big = Math.floorMod(idx_big[i],nx_big);
            if(col_big==0)
                col_big=0;
            //row idx (ref):
            long row_big = Math.round(1+(idx_big[i]-col_big)/nx_big);
            //col idx (inset):
            col[i] = col_big-lambda;
            //row idx (inset):
            row[i] = row_big-mu;
            //idx (inset grid):
            // and pad with zero indexes out of inset grid:
            idx[i]=col[i]+nx*(row[i]-1);
            if( col[i] < 1 || col[i] > nx || row[i] < 1 || row[i] > ny){
                idx[i] = -1;
            }
        }
        return  new idx_ref2inset_gridResults(row, col, idx);
    }

    private int[][] idx_ref2inset_gridCompact(int[][] idx_big, long nx_big, long ny_big, long nx, long ny, long lambda, long mu){
        idx_ref2inset_gridV2Results res = idx_ref2inset_gridV2(idx_big, nx_big, ny_big, nx, ny, lambda, mu);
        int[][] idx= new int[res.getNewDim()][2];
        int counter = 0;
        //removeing elements from free_edges_extended
        //NOTE: sub 1 because MATLAB indexes starts from 1s, java from 0s
        for(int i=0;i<res.getIdx().length;++i){
            if(res.getIdx()[i][0]!=-1 && res.getIdx()[i][1] != -1){
                idx[counter][0] = res.getIdx()[i][0] -1;
                idx[counter][1] = res.getIdx()[i][1] -1;
                counter++;
            }
        }
        return idx;
    }


    private idx_ref2inset_gridV2Results idx_ref2inset_gridV2(int[][] idx_big, long nx_big, long ny_big, long nx, long ny, long lambda, long mu){
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
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("idx_ref2inset_grid: inset grid not within reference grid!");
            debug.CloseFile();
            System.exit(0);
        }
        if((Utility.any(idx_big,"<",1)) || Utility.any(idx_big,">",(int)(nx_big*ny_big))){
            System.out.println("grid index not within reference grid!");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("idx_ref2inset_grid: grid index not within reference grid!");
            debug.CloseFile();
            System.exit(0);
        }
        int[][] idx = new int[idx_big.length][idx_big[0].length];
        int minusOne=0;
        for(int i=0;i<idx_big.length;++i){
            for(int j=0;j<idx_big[0].length;j++){
                int _col_big = (int)Math.floorMod(idx_big[i][j],nx_big);
                if(_col_big == 0)
                    _col_big = (int)nx_big;
                int _row_big = Math.round(1+(idx_big[i][j]-_col_big)/nx_big);
                int _col = _col_big-(int)lambda;
                int _row = _row_big-(int)mu;
                idx[i][j] = _col+(int)nx*(_row-1);
                if( _col < 1 || _col > nx || _row < 1 || _row > ny){
                    idx[i][j] = -1;
                }
            }
            boolean found=false;
            for(int j=0;j<idx_big[0].length;++j){
                if(idx[i][j]==-1)
                    found=true;
            }
            if(found)
                minusOne++;
        }
        int newDim = idx.length-minusOne;
        return new idx_ref2inset_gridV2Results(idx, newDim);
    }

    private mdata_gridResults mdata_grid(ArrayList<Double> lat, ArrayList<Double> lon){
        //converts 2 lat-lon 1-dimensional arrays (Nlat,1) and (Nlon,1) into :
        //a) a list of node coordinates (Nlat*Nlon, 2)
        //b) a meshgrid matrix       (Nlat,   Nlon)
        int NN = lat.size() * lon.size();
        meshgridResults out = Utility.meshgrid(lat,lon);
//        double[][] yg = out.getX();
//        double[][] xg = out.getY();
        int[] aDim = new int[2];
        aDim[0]=NN;
        aDim[1]=1;
        double[][] xa = Utility.reshape(out.getY(), aDim, this.paths.getOutDir());
        double[][] ya = Utility.reshape(out.getX(), aDim, this.paths.getOutDir());
//        double[][] xa = Utility.reshape(xg, aDim);
//        double[][] ya = Utility.reshape(yg, aDim);
        double[][] xy = new double[NN][2];
        for(int i=0;i<NN;++i){
            xy[i][0]=xa[i][0];
            xy[i][1]=ya[i][0];
        }
        return new mdata_gridResults(xy,out.getY(),out.getX());
    }

    private grid_extreme_coordsResults grid_extreme_coords(ArrayList<Double> lat, ArrayList<Double> lon, double[][] field_in){
        //compute coordinates of grid extreme nodes,
        //both for the reference grid (grid read from DB)
        //and the inset grid

        //reference grid:
        //Finding indexes
        ArrayList<Integer> DB_lat_row = new ArrayList<>();
        ArrayList<Integer> DB_lon_row = new ArrayList<>();
        for(int i=0;i<lat.size();++i){
            if( (lat.get(i) >= this.sGrid.getDB_bbox__lat_min()) && (lat.get(i) <= this.sGrid.getDB_bbox__lat_max()) ){
                DB_lat_row.add(i);
            }
        }

        for(int i=0;i<lon.size();++i){
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

        ArrayList<Double> lat_red = new ArrayList<>();
        ArrayList<Double> lon_red = new ArrayList<>();
        double[][] field_out = null;

        //if barflag ==1, we are in creating graph DB mode, so we return empty lat_red, lon_red and field_out
        if(bar_flag == 2){//reading DB mode
            //Insert grid
            ArrayList<Integer> lat_row = new ArrayList<>();
            ArrayList<Integer> lon_row = new ArrayList<>();
            for(int i=0;i<lat.size();++i){
                if( (this.sGrid.getBbox__lat_min() <= lat.get(i)) && (lat.get(i) <= this.sGrid.getBbox__lat_max()) ){
                    lat_row.add(i);
                }
            }
            for(int i=0;i<lon.size();++i){
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
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("\t\tGrid size: "+this.sGrid.getInset_Nx()+"x"+this.sGrid.getInset_Ny()+" = "+(this.sGrid.getInset_Ny()*this.sGrid.getInset_Nx()));
        this.logFile.WriteLog("\t\t[min, max] grid depth (m): ["+Utility.min2d(field_out)+", "+Utility.max2d(field_out)+"]");
        this.logFile.CloseFile();
        return new grid_extreme_coordsResults(lat_red,lon_red,field_out);
    }


    private readout_coastResults readout_coast(){
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
            int Nexternal=0;
            ArrayList<Double> lat_tmp = new ArrayList<>();
            ArrayList<Double> lon_tmp = new ArrayList<>();
            //Opening file
            FileReader file = new FileReader(this.paths.getCoastlineDB());
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
            ArrayList<Double> lon_ext = new ArrayList<>();
            ArrayList<Double> lat_ext = new ArrayList<>();
            ArrayList<Double> lon_int = new ArrayList<>();
            ArrayList<Double> lat_int = new ArrayList<>();
            for(int i=0;i<totalLines;++i){
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
            return new readout_coastResults(lon_ext, lat_ext, lon_int, lat_int);
        } catch(Exception e){
            System.out.println("Lines readed: "+totalLines);
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("readout_coast: "+e.getMessage());
            debug.CloseFile();
            e.printStackTrace();
            return null;
        }
    }

    private readout_bathyResults readout_bathy(int bathy_code){
        String filename = this.paths.getBathymetryDB();
//        String varname = "";
//        String lonname = "";
//        String latname = "";
//        switch(bathy_code){
//            case 1:
//                //MedOneMin:
//                filename = "inputFiles/bathy/MedOneMin/med_one_min_single_vs2.nc";
//                varname = "z";
//                lonname = "longitude";
//                latname = "latitude";
//                break;
//            case 2:
//                //GEBCO_08:
//                filename = "inputFiles/bathy/GEBCO/gebco_08_10_34_15_38.nc";
//                varname = "z";
//                lonname = "x_range";
//                latname = "y_range";
//                break;
//            case 3:
//                //EMODnet:
//                filename = "inputFiles/bathy/EMODnet/Adriatic_Ionian_central__MedSea_mean.nc";
//                varname = "depth_average";
//                lonname = "lon";
//                latname = "lat";
//                break;
//            default:
//                System.out.println("unknow bathy DB code!");
//                MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
//                debug.WriteLog("readout_bathy: unknow bathy DB code!");
//                debug.CloseFile();
//                break;
//        }
        //Parsing file:
//        MyNetCDFParser test = new MyNetCDFParser(filename);
        MyNetCDFParser test = new MyNetCDFParser(filename);
        parseMedOneMinResults out = test.parseMedOneMin(this.paths.getOutDir());
        if(out == null){
            System.out.println("Parsing fail!");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("readout_bathy: Parsing fail!");
            debug.CloseFile();
            return null;
        }
//        ArrayList<Double> latTmp = out.getLat();
//        ArrayList<Double> lonTmp = out.getLon();
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
        for(int i=0;i<out.getLat().size();i++){
            if(out.getLat().get(i)>=lat_min && out.getLat().get(i)<=lat_max){
                latOkIndexes.add(i);
            }
        }
        for(int i=0;i<out.getLon().size();i++){
            if(out.getLon().get(i)>=lon_min && out.getLon().get(i)<=lon_max){
                lonOkIndexes.add(i);
            }
        }

        ArrayList<Double> lat = new ArrayList<>();
        ArrayList<Double> lon = new ArrayList<>();
        double[][] z = new double[latOkIndexes.size()][lonOkIndexes.size()];
        for(int i=0;i<latOkIndexes.size();++i) {
            lat.add(out.getLat().get(latOkIndexes.get(i)));
        }
        for(int i=0;i<lonOkIndexes.size();++i){
            lon.add(out.getLon().get(lonOkIndexes.get(i)));
        }
        for(int i=0;i<latOkIndexes.size();++i){
            for(int j=0;j<lonOkIndexes.size();++j){
                z[i][j] = zTmp[latOkIndexes.get(i)][lonOkIndexes.get(j)];
            }
        }
        return new readout_bathyResults(lat, lon, z);
    }

    //NEVER USED!
//    private nodes_free_form_barrierResults nodes_free_from_barrier(int[] nodes, double[][] xy,ArrayList<Double> x_islands, ArrayList<Double> y_islands, ArrayList<Double> x_continent, ArrayList<Double> y_continent){
//        // % finds nodes not on the landmass (both continent and islands)
//        // %
//        // % removal of islands outside of bbox:
//        // % computation of nodes not on the landmasses:
//        int Nn = nodes.length;
//
//        //removal of islands outside of bbox:
//        int n_prima = y_islands.size();
//
//        boolean[] bool_in = new boolean[y_islands.size()];
//        for(int i=0;i<y_islands.size();i++){
//            if((this.sGrid.getDB_xi()<= x_islands.get(i)) && (x_islands.get(i)<=this.sGrid.getDB_xf())
//                    && (this.sGrid.getDB_yi()<=y_islands.get(i) && (y_islands.get(i)<=this.sGrid.getDB_yf())) ){
//                bool_in[i]=true;
//            }else{
//                bool_in[i]=false;
//            }
//        }
//        boolean[] bool_nan = new boolean[y_islands.size()];
//        for(int i=0;i<y_islands.size();i++){
//            if(Double.isNaN(x_islands.get(i)) && Double.isNaN(y_islands.get(i))){
//                bool_nan[i] = true;
//            } else {
//                bool_nan[i] = false;
//            }
//        }
//        ArrayList<Integer> isl_in_bb = new ArrayList<>();
//        for(int i=0;i<y_islands.size();i++){
//            if(bool_in[i] || bool_nan[i]){
//                isl_in_bb.add(i);
//            }
//        }
//        ArrayList<Double> x_islands_tmp = new ArrayList<>();
//        ArrayList<Double> y_islands_tmp = new ArrayList<>();
//        for(Integer element : isl_in_bb){
//            x_islands_tmp.add(x_islands.get(element));
//            y_islands_tmp.add(y_islands.get(element));
//        }
//        x_islands = x_islands_tmp;
//        y_islands = y_islands_tmp;
//
//        int n_dopo = y_islands.size();//includes NaNs separators between islands (i.e., 5608 elements)
//        int red_fact = Math.round((1 - n_dopo/n_prima)*100);
//
//        ArrayList<Double> y_coast = new ArrayList<>();
//        ArrayList<Double> x_coast = new ArrayList<>();
//        for(int i=0; i<y_islands.size();i++){
//            y_coast.add(y_islands.get(i));
//            x_coast.add(x_islands.get(i));
//        }
//        y_coast.add(Double.NaN);
//        x_coast.add(Double.NaN);
//
//        for(int i=0;i<y_continent.size();i++){
//            y_coast.add(y_continent.get(i));
//            x_coast.add(x_continent.get(i));
//        }
//
//        //--------------------------------------------------------------------------------------------
//        //computation of nodes not on the landmasses:
//        ArrayList<Integer> free_nodes = new ArrayList<>();
//        MyFileWriter fid = new MyFileWriter(true,"freeNodes_DB.dat.permil");
//        for(int ie=0;ie<Nn;ie++){
//            int frac_done = new Double(Math.floor(1000*ie/Nn)).intValue();
//            fid.WriteLine(""+frac_done);
//
//            double xP = xy[ie][0];
//            double yP = xy[ie][1];
//
//            inpolygonResults tmp = Utility.inpolygon(xP, yP, x_islands, y_islands);
//            boolean IN_i = tmp.getIn();
//            boolean ON_i = tmp.getOn();
//            tmp = Utility.inpolygon(xP, yP, x_continent, y_continent);
//            boolean IN_c = tmp.getIn();
//            boolean ON_c = tmp.getOn();
//
//            boolean canAdd = (IN_i == false) && (ON_i == false) && (IN_c == true) && (ON_c == false);
//
//            if(canAdd){
//                free_nodes.add(nodes[ie]);
//            }
//        }
//        if(free_nodes.size() == 0){
//            System.out.println("no free nodes found in graph");
//            MyFileWriter debug = new MyFileWriter("","debug",false);
//            debug.WriteLog("nodes_free_from_barrier: no free nodes found in graph");
//            debug.CloseFile();
//            System.exit(0);
//        }
//        fid.CloseFile();
//        return new nodes_free_form_barrierResults(free_nodes.size(),free_nodes,red_fact);
//    }


    private Fields_regriddingResults Fields_regridding(ArrayList<Double> lat_bathy_Inset, ArrayList<Double> lon_bathy_Inset, double[][] bathy_Inset, double[][] lsm_mask,
                                                        double[][] ship_v_LUT, ArrayList<Double> H_array_m, double estGdtDist){
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
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
        double wind_t1;
        if (this.optim.getWindModel()==11 || this.optim.getWindModel()==12){//ECMWF
            wind_t1 = 12.0;//ECMWF wind file (analysis) start time is 1200 UTC
        } else if (this.optim.getWindModel() == 2){
            wind_t1 = 15.0; //COSMO-ME wind file (forecast) start time is 1500 UTC  *** use also analysis in the future!
        } else {
            System.out.println("unknown wind model code!");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("Fields_regridding: unknown wind model code!");
            debug.CloseFile();
            System.exit(0);
        }
        //String dateListFile = "inputFiles/fields/an_dates_DB.txt";
//        String dateListFile = this.paths.getAnalysisDB();
//        MyTxtParser file = new MyTxtParser(dateListFile, this.paths.getOutDir());
//        this.tGrid.setLatest_date(file.tail(1).get(0));
        boolean isGRIB = false;
        String noExtension = "";
        if(this.paths.getForecastFile().substring(this.paths.getForecastFile().length()-3).equals("grb")){
            //grib file
            isGRIB = true;
            noExtension = this.paths.getForecastFile().substring(0, paths.getForecastFile().length()-4);
        } else {
            //netcdf
            noExtension = this.paths.getForecastFile().substring(0, paths.getForecastFile().length()-3);
        }
        this.tGrid.setLatest_date(noExtension.substring(noExtension.length()-8));

        String l_date_str = this.tGrid.getLatest_date()+"1200";
        //String l_date_str = this.tGrid.getLatest_date()+"0000";//TEST
        long l_num = Utility.datenum(l_date_str,"yyyyMMddHHmm", this.paths.getOutDir()); // taken at 1200 UTC (analysis time)

        this.tGrid.setDepDateTime(Utility.datenum(2000+this.dep_datetime.getYear(),this.dep_datetime.getMonth(),
                this.dep_datetime.getDay(),this.dep_datetime.getHour(),this.dep_datetime.getMin(), this.paths.getOutDir()));

        //hrs elapsed between latest analysis and departure time
        int cento = 100;
        long deltaHr_anls;
        if(this.forcing.getAnalytic()==1){
            deltaHr_anls = 1;
        } else{
            deltaHr_anls = Utility.hoursBetween("yyyyMMddHHmm",l_date_str,dep_datetime.getDepDateTime(), this.paths.getOutDir());
        }


        this.tGrid.setWave_dep_TS(Math.round(1+ deltaHr_anls - (wave_t1 - this.constants.getTwelve())));

        this.tGrid.setWind_dep_TS(Math.round(1+ deltaHr_anls - (wave_t1 - this.constants.getTwelve())));

        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("\tlatest analysis date and time: "+Utility.datestr(l_num,12,00, this.paths.getOutDir()));
//        this.logFile.WriteLog("\tlatest analysis date and time: "+Utility.datestr(l_num,00,00, this.paths.getOutDir()));//TEST
        this.logFile.WriteLog("\tdeparture date and time: "+Utility.datestr(this.tGrid.getDepDateTime(),this.dep_datetime.getHour(),this.dep_datetime.getMin(), this.paths.getOutDir()));
        this.logFile.CloseFile();

        //-----------------------------------------------------------------------------------------------------
        //Fields reading:
        if(this.forcing.getAnalytic()==1){
            this.tGrid.setNt(this.tGrid.getMinNt());
            prepare_sqrtY_fieldResults out = prepare_sqrtY_field(lat_bathy_Inset, lon_bathy_Inset);
            return new Fields_regriddingResults(null, out.getVTDH_b(),
                    Utility.NaN3Dmatrix(out.getVTDH_b().length,out.getVTDH_b()[0].length,out.getVTDH_b()[0][0].length),
                    out.getVDIR_b(),Utility.NaN3Dmatrix(out.getVTDH_b().length,out.getVTDH_b()[0].length,out.getVTDH_b()[0][0].length),
                    Utility.NaN3Dmatrix(out.getVTDH_b().length,out.getVTDH_b()[0].length,out.getVTDH_b()[0][0].length));
        } else { //Reading of environmental Fields:
            readout_envFieldsResults envFieldsResults = readout_envFields();
            //this check must be run after input file time size has been determined
            check_start_timestep(deltaHr_anls);

            //-----------------------------------------------------------------------------------------------------
            //Spatial subsetting:

            //Wave data reduction to inset grid:
            mFields_reductionResults waveDataReduction =mFields_reduction(envFieldsResults.getLat_wave(),envFieldsResults.getLon_wave(), envFieldsResults.getVTDH(), envFieldsResults.getVTPK(), envFieldsResults.getVDIR());

            meshgridResults meshGWave = Utility.meshgrid(waveDataReduction.getLon_red(), waveDataReduction.getLat_red());


            //Wave direction into Cartesian components:
            double[][][] X_Inset = new double[waveDataReduction.getOut(2).length][waveDataReduction.getOut(2)[0].length][waveDataReduction.getOut(2)[0][0].length];
            double[][][] Y_Inset = new double[waveDataReduction.getOut(2).length][waveDataReduction.getOut(2)[0].length][waveDataReduction.getOut(2)[0][0].length];
            for(int i=0;i<waveDataReduction.getOut(2).length; ++i){
                for(int j=0;j<waveDataReduction.getOut(2)[0].length; ++j){
                    for(int k=0;k<waveDataReduction.getOut(2)[0][0].length; ++k){
                        X_Inset[i][j][k] = Math.sin(waveDataReduction.getOut(2)[i][j][k]*this.constants.getDeg2rad());
                        Y_Inset[i][j][k] = Math.cos(waveDataReduction.getOut(2)[i][j][k]*this.constants.getDeg2rad());
                    }
                }
            }
            //clear VDIR_Inset
            waveDataReduction.setOut(2, null);
            //Wind  data reduction to inset grid:
            mFields_reductionResults ecmwfDataReduction = mFields_reduction(envFieldsResults.getEcmwf_lat_wind(), envFieldsResults.getEcmwf_lon_wind(), envFieldsResults.getEcmwf_U10m(), envFieldsResults.getEcmwf_V10m());
            mFields_reductionResults cosmoDataReduction = mFields_reduction(envFieldsResults.getCosmo_lat_wind(), envFieldsResults.getCosmo_lon_wind(), envFieldsResults.getCosmo_U10m(), envFieldsResults.getCosmo_V10m());

            //Wind unit conversion
            ecmwfDataReduction.multiplyElementFor(0, this.constants.getMs2kts());
            ecmwfDataReduction.multiplyElementFor(1, this.constants.getMs2kts());

            cosmoDataReduction.multiplyElementFor(0, this.constants.getMs2kts());
            cosmoDataReduction.multiplyElementFor(1, this.constants.getMs2kts());

            //-----------------------------------------------------------------------------------------------------
            //Time processing:
            fieldStatsResults field_res = fieldStats(bathy_Inset, waveDataReduction.getOut(0), waveDataReduction.getOut(1),
                    ecmwfDataReduction.getOut(0), ecmwfDataReduction.getOut(1),
                    cosmoDataReduction.getOut(0), cosmoDataReduction.getOut(1));
            //double ecmwf_dir_avg = field_res.getEcmwf_dir_avg();
            //double ecmwf_dir_std = field_res.getEcmwf_dir_std();
            //double cosmo_dir_avg = field_res.getCosmo_dir_avg();
            //double cosmo_dir_std = field_res.getCosmo_dir_std();
            //Just for GMD paper's route #2 (1589_c) uncomment following line:
            //this.fstats.setWlenght_max(85);

            this.estim_Nt(estGdtDist, Math.max(this.tGrid.getWave_dep_TS(), this.tGrid.getWind_dep_TS()), H_array_m, ship_v_LUT);

            int delta_hr_dep_int = Math.round(deltaHr_anls);
            int[] interpTimes = new int[(int)this.tGrid.getNt()];
            for(int i=0;i<interpTimes.length;++i)
                interpTimes[i] = delta_hr_dep_int+(i-1);

            int[] time_steps = new int[(int) this.tGrid.getMaxNt()];
            //time_steps always counted starting from step #1:
            for(int i=0;i<time_steps.length; ++i)
                time_steps[i] = i+1;
            this.tGrid.setDayHrs(myMod(wave_t1, interpTimes,this.constants.getTwentyfour()));
            double[][][] U10M_at_TS = new double[0][0][0];
            double[][][] V10M_at_TS = new double[0][0][0];
            if(this.optim.getWindModel() == 11 || this.optim.getWindModel() == 12){//ecmwf

                U10M_at_TS = Utility.interp1(envFieldsResults.getEcmwf_wind_origTimes(),
                        ecmwfDataReduction.getOut(0), interpTimes, this.paths.getOutDir());

                V10M_at_TS = Utility.interp1(envFieldsResults.getEcmwf_wind_origTimes(),
                        ecmwfDataReduction.getOut(1), interpTimes, this.paths.getOutDir());
            } else if(this.optim.getWindModel() == 2){//cosmo-me

                U10M_at_TS = Utility.interp1(envFieldsResults.getCosmo_wind_origTimes(),
                        cosmoDataReduction.getOut(0), interpTimes, this.paths.getOutDir());

                V10M_at_TS = Utility.interp1(envFieldsResults.getCosmo_wind_origTimes(),
                        cosmoDataReduction.getOut(1), interpTimes, this.paths.getOutDir());
            }


            double[][][] VTDH_times = new double[waveDataReduction.getOut(0).length][interpTimes.length][waveDataReduction.getOut(0)[0][0].length];
            double[][][] VTPK_times = new double[waveDataReduction.getOut(1).length][interpTimes.length][waveDataReduction.getOut(1)[0][0].length];
            double[][][] X_times = new double[X_Inset.length][interpTimes.length][X_Inset[0][0].length];
            double[][][] Y_times = new double[Y_Inset.length][interpTimes.length][Y_Inset[0][0].length];

            for(int i=0;i<VTDH_times.length;++i){
                for(int j=0;j<VTDH_times[0].length;++j){
                    for(int k=0;k<VTDH_times[0][0].length;++k){
                        VTDH_times[i][j][k]=waveDataReduction.getOut(0)[i][interpTimes[j]][k];
                        VTPK_times[i][j][k]=waveDataReduction.getOut(1)[i][interpTimes[j]][k];
                        X_times[i][j][k]=X_Inset[i][interpTimes[j]][k];
                        Y_times[i][j][k]=Y_Inset[i][interpTimes[j]][k];
                    }
                }
            }


            //-----------------------------------------------------------------------------------------------------
            // seaOverLand:
            this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
            this.logFile.WriteLog("\tseaOverLand extrapolation...");
            this.logFile.CloseFile();
            //Wave fields processing:
            ArrayList<double[][][]> seaOverLand_3stepsOut = seaOverLand_3steps(lon_bathy_Inset, lat_bathy_Inset,
                    lsm_mask, meshGWave.getX(), meshGWave.getY(),VTDH_times, VTPK_times, X_times, Y_times);
            X_times = seaOverLand_3stepsOut.get(2);
            Y_times = seaOverLand_3stepsOut.get(3);
            double[][][] VDIR_out = new double[X_times.length][X_times[0].length][X_times[0][0].length];
            for(int k=0;k<VDIR_out.length;++k){
                for(int i=0;i<VDIR_out[0].length;++i){
                    for(int j=0;j<VDIR_out[0][0].length;++j){
                        double tmp = Math.atan2(Y_times[k][i][j], X_times[k][i][j]);
                        tmp = changeDirRule(tmp);
                        VDIR_out[k][i][j] = tmp/this.constants.getDeg2rad();
                    }
                }
            }

            X_times = null;
            Y_times = null;

            meshgridResults wind_m;
            if (this.optim.getWindModel() == 11 || this.optim.getWindModel() == 12) {//ecmwf
                wind_m = Utility.meshgrid(ecmwfDataReduction.getLon_red(), ecmwfDataReduction.getLat_red());
            } else {
                wind_m = Utility.meshgrid(cosmoDataReduction.getLon_red(), cosmoDataReduction.getLat_red());
            }


            ArrayList<double[][][]> seaOut = seaOverLand_3steps(lon_bathy_Inset,lat_bathy_Inset,lsm_mask,
                    wind_m.getX(),wind_m.getY(), U10M_at_TS,V10M_at_TS);
            U10M_at_TS = seaOut.get(0);
            V10M_at_TS = seaOut.get(1);

            double[][][] windDIR_out = new double[V10M_at_TS.length][V10M_at_TS[0].length][V10M_at_TS[0][0].length];
            double[][][] windMAGN_out =new double[V10M_at_TS.length][V10M_at_TS[0].length][V10M_at_TS[0][0].length];

            for(int k=0;k<windDIR_out.length;++k){
                for(int i=0;i<windDIR_out[0].length;++i){
                    for(int j=0;j<windDIR_out[0][0].length;++j){
                        double tmp = Math.atan2(V10M_at_TS[k][i][j], U10M_at_TS[k][i][j]);
                        tmp = changeDirRule(tmp);
                        windDIR_out[k][i][j] = tmp/this.constants.getDeg2rad();
                        windMAGN_out[k][i][j] = Math.sqrt(Math.pow(U10M_at_TS[k][i][j],2)+Math.pow(V10M_at_TS[k][i][j],2));
                    }
                }
            }

            return new Fields_regriddingResults(time_steps, seaOverLand_3stepsOut.get(0),
                    seaOverLand_3stepsOut.get(1), VDIR_out, windMAGN_out, windDIR_out);
        }
    }

    private double[] myMod(double value, int[] array, double under){
        double[] out = new double[array.length];
        for(int i=0;i<out.length;++i)
            out[i] = (value+array[i])%under;
        return out;
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
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("\tseaOverLand_3steps: varargs must be between 1 and 4");
            debug.CloseFile();
            System.exit(0);
        }

        if(Math.min(lon_f.length, lon_f[0].length) < this.sGrid.getMinNoGridPoints()){
            System.out.println("seaOverLand_3steps: too small or too narrow bounding box");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("\tseaOverLand_3steps: too small or too narrow bounding box");
            debug.CloseFile();
            System.exit(0);
        }

        int n_loops = 50;

        double[][][] myfield_bathy = new double[lon_bathy.size()][(int)this.tGrid.getNt()][lat_bathy.size()];
        ArrayList<double[][][]> out = new ArrayList<>();
        for(double[][][] myfield : varargs){
            //(1) extrapolation - % #GM: check also mdata_EWeights.m:
            for(int it=0;it<(int)this.tGrid.getNt(); ++it){
                double[][] myfield_mat;
                double[][][] tmp = new double[myfield.length][1][myfield[0][0].length];
                for(int z=0;z<myfield.length;++z){
                    for(int j =0; j<myfield[0][0].length; ++j){
                        tmp[z][0][j] = myfield[z][it][j];
                    }
                }
                if(this.forcing.getWind()!=1){
                    myfield_mat = SeaOverLand(Utility.squeeze(tmp, this.paths.getOutDir()),n_loops);
                } else {
                    myfield_mat = Utility.squeeze(tmp, this.paths.getOutDir());
                }

                // (2) regridding to bathy-grid:

                //check that bbox is larger than minimum allowed area:
                if(it==0){
                    int size_min = Math.min(myfield_mat.length, myfield_mat[0].length);
                    int nan_rank = Utility.rank(Utility.convertDouble(Utility.isnan(myfield_mat)));
                    if((size_min-nan_rank) < 2){
                        System.out.println("seaOverLand_3steps: too few sea grid points");
                        MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
                        debug.WriteLog("\tseaOverLand_3steps: too few sea grid points");
                        debug.CloseFile();
                    }
                }
                double[] lon_bathyArray = Utility.arrayfy(lon_bathy);
                double[] lat_bathyArray = Utility.arrayfy(lat_bathy);
                double[][] transposedTmpMtx = Utility.transposeMatrix(Utility.interp2(lon_f, lat_f, myfield_mat, lon_bathyArray, lat_bathyArray));
                for(int i=0;i<transposedTmpMtx.length;++i){
                    for(int j=0;j<transposedTmpMtx[0].length;++j){
                        myfield_bathy[i][it][j]=transposedTmpMtx[i][j];
                    }
                }
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
                    for(int it=0;it<(int) this.tGrid.getNt(); ++it){
                        double[][] tmp = Utility.squeeze(myfield_bathy,1,it, this.paths.getOutDir());
                        double[][] vv1 = Utility.MatrixComponentXcomponent(tmp,Utility.transposeMatrix(lsm_mask), this.paths.getOutDir());
                        for(int i=0;i<vv1.length;i++){
                            for(int j=0;j<vv1[0].length;j++){
                                myfield_Inset[j][it][i] = vv1[i][j];
                            }
                        }
                    }
                } else {
                    myfield_Inset = myfield_bathy;
                }
            } else{
                myfield_Inset = myfield_bathy;
            }
            out.add(myfield_Inset);
        }
        return out;
    }

    private static double[][] SeaOverLand(double[][] matrice_in, int loop){
        //source: Nicoletta Fabroni (SINCEM, Ravenna), developed for Relocatable HOPS package
        //% versione ricevuta da I.Federico via S. Falchetti in data 24/6/2013

        double[][] dummy = new double[matrice_in.length+2][matrice_in[0].length+2];
        for(int i=0;i<dummy.length;++i){
            for(int j=0;j<dummy[0].length;++j){
                if(i==0 || j==0 || i==(dummy.length-1) || j==(dummy[0].length-1)){
                    dummy[i][j] = 9999;
                }else{
                    dummy[i][j] = matrice_in[i-1][j-1];
                }
            }
        }

        ArrayList<Integer> idx = Utility.find(dummy, Double.NaN);

        for(int i=0;i<dummy.length;++i){
            for(int j=0;j<dummy[0].length;++j){
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
            for(int i=0;i<neighbors.length;++i){
                for(int j=0;j<neighbors[0].length;++j){
                    neighbors[i][j] = idx.get(j)+neighbor_offsets[i];
                }
            }

            double[][] mat = new double[neighbors.length][neighbors[0].length];
            for(int i=0;i<mat.length;++i){
                for(int j=0;j<mat[0].length;++j){
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

            for(int i=0;i<nans.length;++i){
                for(int j=0;j<nans[0].length;++j){
                    if(nans[i][j]){
                        mat[i][j] = 0;
                    } else{
                        snn[j]++;
                    }
                }
            }

            double[] media = Utility.sum(mat);
            for(int i=0;i<snn.length;++i){
                media[i]= media[i]/(snn[i]+1);
            }

            for(int i=0;i<idx.size();++i){
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
        for(int i=0;i<matrice_out.length;++i){
            for(int j=0;j<matrice_out[0].length;++j){
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
            for(int i=0;i<ship_v_LUT.length;++i){
                ship_v_LUT_1stCol[i] = ship_v_LUT[i][0];
            }
            v_kts_min = Utility.interp1(H_array_m, ship_v_LUT_1stCol, this.fstats.getWheight_max(), this.paths.getOutDir());
            v_kts_max = Utility.interp1(H_array_m, ship_v_LUT_1stCol, this.fstats.getWheight_min(), this.paths.getOutDir());
        } else{//sailboat
            v_kts_min=Utility.min(ship_v_LUT);
            v_kts_max=Utility.max(ship_v_LUT);
        }

        double v_kts_signif = v_kts_min*(1-fract) + v_kts_max*fract;

        double Nt_1 = this.tGrid.getMaxNt() - start_timestep + 1;
        double Nt_2b = Math.ceil(estGdtDist/v_kts_signif);

        double Nt_2 = Math.max(this.tGrid.getMinNt(), Nt_2b);
        this.tGrid.setNt(Math.min(Nt_1, Nt_2));

        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("\t# time steps of forecast file employed: "+this.tGrid.getNt());
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

        //WIND
        double[][][] ecmwf_wind10M_magn = new double[ecmwf_U10M.length][ecmwf_U10M[0].length][ecmwf_U10M[0][0].length];
        double[][][] cosmo_wind10M_magn = new double[cosmo_U10M.length][cosmo_U10M[0].length][cosmo_U10M[0][0].length];
        double[][][] TmpEcmwf_cos_avg = new double[ecmwf_U10M.length][ecmwf_U10M[0].length][ecmwf_U10M[0][0].length];
        double[][][] TmpEcmwf_sin_avg = new double[ecmwf_V10M.length][ecmwf_V10M[0].length][ecmwf_V10M[0][0].length];
        //WIND direction
        for(int i=0;i<ecmwf_U10M.length;++i){
            for(int j=0;j<ecmwf_U10M[0].length;++j){
                for(int k=0;k<ecmwf_U10M[0][0].length;++k){
                    ecmwf_wind10M_magn[i][j][k] = Math.sqrt(Math.pow(ecmwf_U10M[i][j][k],2) + Math.pow(ecmwf_V10M[i][j][k],2));
                    cosmo_wind10M_magn[i][j][k] = Math.sqrt(Math.pow(cosmo_U10M[i][j][k],2) + Math.pow(cosmo_V10M[i][j][k],2));
                    TmpEcmwf_cos_avg[i][j][k] = ecmwf_U10M[i][j][k]/ecmwf_wind10M_magn[i][j][k];
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
        for(int i=0;i<TmpCosmo_cos_avg.length; ++i){
            for(int j=0;j<TmpCosmo_cos_avg[0].length; ++j){
                for(int k=0;k<TmpCosmo_cos_avg[0][0].length; ++k){
                    TmpCosmo_cos_avg[i][j][k] = cosmo_U10M[i][j][k]/cosmo_wind10M_magn[i][j][k];
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
        for(int i=0;i<lambda.length; ++i){
            for(int j=0;j<lambda[0].length;++j){
                for(int k=0;k<lambda[0][0].length;++k){
                    lambda[i][j][k] = this.constants.getG0()/twopi* Math.pow(wave_period[i][j][k],2); //[m]
                }
            }
        }
        return lambda;
    }


    //NEVER USED
//    private double[][][] wave_dispersion(double[][][] wave_period, double[][][] depth){
//        double[][][] lambda = wave_dispersion(wave_period);
//        double[][][] fmk = Fenton_McKee_factor(wave_period,depth);
//        // %   -) If also depth provided       --> generic depth formula after:
//        // %   Fenton, JD and McKee, WD (1990)
//        // %
//        // %   wave_period [s]
//        // %   depth       [m]
//        // %   lambda      [m]
//        // %
//        if(this.forcing.getDeepWaterApprox()!=1){
//            for(int i=0;i<lambda.length;i++){
//                for(int j=0;j<lambda[0].length;j++){
//                    for(int k=0;k<lambda[0][0].length;k++){
//                        lambda[i][j][k] = lambda[i][j][k]*fmk[i][j][k];
//                    }
//                }
//            }
//        }
//        return lambda;
//    }


    //NEVER USED
//    private double[][][] Fenton_McKee_factor(double[][][] wave_period, double[][][] depth){
//        // % Factor multiplying deep-water wavelength in
//        // % Fenton, JD and McKee, WD (1990)
//        // %
//        // % (it leads to a reduced wavelength in shallow water)
//        double twopi = Math.PI *2;
//        double[][][] lambda_0 = new double[wave_period.length][wave_period[0].length][wave_period[0][0].length];
//        for(int i=0;i<wave_period.length;i++){
//            for(int j=0;j<wave_period[0].length;j++){
//                for(int k=0;k<wave_period[0][0].length;k++){
//                    lambda_0[i][j][k] = this.constants.getG0()/twopi * Math.pow(wave_period[i][j][k],2);
//                }
//            }
//        }
//        double[][][] fmk = new double[wave_period.length][wave_period[0].length][wave_period[0][0].length];
//        for(int i=0;i<wave_period.length;i++){
//            for(int j=0;j<wave_period[0].length;j++){
//                for(int k=0;k<wave_period[0][0].length;k++){
//                    Math.pow( Math.tanh( Math.pow( twopi*depth[i][j][k] / lambda_0[i][j][k],(3/4) ) ), (2/3) );
//                }
//            }
//        }
//        return fmk;
//    }

    private double[] Fenton_McKee_factor(double[][] wave_period, int colIndex, double[] depth){
//         % Factor multiplying deep-water wavelength in
//         % Fenton, JD and McKee, WD (1990)
//         %
//         % (it leads to a reduced wavelength in shallow water)
        double twopi = Math.PI *2;
        double[] fmk = new double[wave_period.length];
        for(int i=0;i<fmk.length;++i){
            double lambda_0 = this.constants.getG0()/twopi * Math.pow(wave_period[i][colIndex],2);
            fmk[i]=Math.pow( Math.tanh( Math.pow( twopi*depth[i] / lambda_0, (3.0/4.0) ) ) , (2.0/3.0));
        }
//        for(int i=0;i<fmk.length;++i){
//            double lambda_0 = this.constants.getG0()/twopi * Math.pow(wave_period[i][colIndex],2);
//            double num = depth[i]*twopi;
//            fmk[i]=Math.pow(Math.pow(Math.tanh(num/lambda_0) , (3.0/4.0)), (2.0/3.0));
//        }
        return fmk;
    }


    private double nanmean2(double[][][] A_mat){
        //Abstract: compute mean of matrix elements of A_mat, even in presence of NaNs.
        ArrayList<Double> elements = new ArrayList<>();
        for(int k=0;k<A_mat.length;++k){
            for(int i=0;i<A_mat[0].length;++i){
                for(int j=0;j<A_mat[0][0].length;++j){
                    if(!Double.isNaN(A_mat[k][i][j])){
                        elements.add(A_mat[k][i][j]);
                    }
                }
            }
        }
        int nElements = elements.size();
        double cumSum = 0.0;
        for(int i=0;i<elements.size();++i){
            cumSum+=elements.get(i);
        }
        return (cumSum/nElements);
    }

    private mFields_reductionResults mFields_reduction(double[] lat, double[] lon, double[][][]... varargin){
        mFields_reductionResults retThis = new mFields_reductionResults();
        double meshRes = 1.0/this.sGrid.getInvStepFields();

        //reduction to given bounding box
        findResults latTmp = Utility.find(lat, "<=", this.sGrid.getBbox__lat_min()-meshRes, "<=", this.sGrid.getBbox__lat_max()+meshRes);
        int[] row_lat = latTmp.getIndexes();
        findResults lonTmp = Utility.find(lon, "<=", this.sGrid.getBbox__lon_min()-meshRes, "<=", this.sGrid.getBbox__lon_max()+meshRes);
        int[] row_lon = lonTmp.getIndexes();

        if(row_lon.length<this.sGrid.getMinNoGridPoints() || row_lat.length<this.sGrid.getMinNoGridPoints()){
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("mFields_reduction: check_start_timestep: too small or too narrow bounding box");
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
        int n_in =0;
        double[][][] last=new double[0][0][0];
        for(double[][][] mat3d : varargin){
            n_in++;
            last=mat3d;
            double[][][] out=new double[row_lon.length][mat3d[0].length][row_lat.length];
            for(int k=0;k<row_lon.length;++k){
                for(int i=0;i<out[0].length;++i){
                    for(int j=0;j<row_lat.length;++j){
                        out[k][i][j] = mat3d[row_lon[k]][i][row_lat[j]];
                    }
                }
            }
            retThis.addOut(out);
        }
        long n_elementsIn = Utility.numel(last);
        long n_elementsOut = Utility.numel(retThis.getOut(n_in-1));
        double red_percent = (n_elementsOut+0.0)/(n_elementsIn+0.0);
        if(red_percent<1.0){
            this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
            this.logFile.WriteLog("\t\tData reduction: "+(red_percent*100)+"%");
            this.logFile.CloseFile();
        }
        return retThis;
    }


    private void check_start_timestep(long deltaHr_anls){
        if(deltaHr_anls <= 0){
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("mFields_reduction: departure time too far back in the past!");
            debug.CloseFile();
            System.exit(0);
        }

        if(this.forcing.getWave()==1){
            if(this.tGrid.getWave_dep_TS() < 1){
                MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
                debug.WriteLog("check_start_timestep: departure time too far back in the past!");
                debug.CloseFile();
                System.exit(0);
            }

            if(this.tGrid.getWave_dep_TS() > this.tGrid.getMaxNt()){
                MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
                debug.WriteLog("check_start_timestep: departure time too far in future!");
                debug.CloseFile();
                System.exit(0);
            }
        }

        if(this.forcing.getWind() == 1){
            if(this.tGrid.getWind_dep_TS() < 1){
                MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
                debug.WriteLog("check_start_timestep: departure time too far back in the past!");
                debug.CloseFile();
                System.exit(0);
            }

            if(this.tGrid.getWind_dep_TS() > this.tGrid.getMaxNt()){
                MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
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
                    MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
                    debug.WriteLog("readout_envFields: unknown waveModel!");
                    debug.CloseFile();
                    System.exit(0);
                }
            }
            this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
            this.logFile.WriteLog("\tattempt reading "+model_str+" wave forecast data...");
            this.logFile.CloseFile();

            readout_mWaveResults readout_mWave = readout_mWave();
            if(!readout_mWave.isWave_status()){
                System.out.println("readout_envFields: "+ readout_mWave.getWave_filename() + "wave forecast not found!");
                MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
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
            this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
            this.logFile.WriteLog("\tdeparture at time step #"+(this.tGrid.getWave_dep_TS()-1)+" of hourly interpolated file");
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
                        MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
                        debug.WriteLog("readout_envFields: unknown windModel!");
                        debug.CloseFile();
                        System.exit(0);
                    }
                }
            }
            //------------------------------------------------------------------------
            this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
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
        for(int i=0;i<lat_windTmp.size(); ++i)
            lat_wind[i] = lat_windTmp.get(i);
        ArrayList<Double> lon_windTmp = Utility.linspace(this.sGrid.getBbox__lon_min(), this.sGrid.getBbox__lon_max(), nlon);
        double[] lon_wind = new double[lon_windTmp.size()];
        for(int i=0;i<lon_windTmp.size(); ++i)
            lon_wind[i] = lon_windTmp.get(i);
        double[][][] U10m = Utility.ones3Dmatrix(nlon, nsteps, nlat);
        double[][][] V10m = Utility.ones3Dmatrix(nlon, nsteps, nlat);
        ArrayList<Double> wind_origTimesTmp = Utility.linspace(0, ntwind, nsteps);
        double[] wind_origTimes = new double[wind_origTimesTmp.size()];
        for(int i=0;i<wind_origTimesTmp.size();++i)
            wind_origTimes[i]=wind_origTimesTmp.get(i);
        return new fake_windFieldsResults(wind_origTimes, lat_wind, lon_wind, U10m, V10m);
    }

    private fake_waveFieldsResults fake_waveFields(int mtwave, int nsteps, int nlat, int nlon){
        ArrayList<Double> lat_waveTmp = Utility.linspace(this.sGrid.getBbox__lat_min(), this.sGrid.getBbox__lat_max(), nlat);
        double[] lat_wave = Utility.arrayfy(lat_waveTmp);
        ArrayList<Double> lon_waveTmp = Utility.linspace(this.sGrid.getBbox__lon_min(), this.sGrid.getBbox__lon_max(), nlon);
        double[] lon_wave = Utility.arrayfy(lon_waveTmp);
        double[][][] VTDH = Utility.ones3Dmatrix(nlon,nsteps,nlat);
        double[][][] VTPK = Utility.ones3Dmatrix(nlon,nsteps,nlat);
        double[][][] VDIR = Utility.ones3Dmatrix(nlon,nsteps,nlat);
        int[] wave_origtimes = new int[mtwave];
        for(int i=0;i<mtwave; ++i)
            wave_origtimes[i]=(i+1);
        return new fake_waveFieldsResults(wave_origtimes, lon_wave, lat_wave, VTDH, VTPK, VDIR);
    }

    private readout_mWaveResults readout_mWave(){
        // % reads out Waves forecast data
        // % physical fields from either WAM or WW3 model:
//        String Stagein_path = "";
//        if(this.optim.getWaveModel() == 10){
//            Stagein_path = "inputFiles/wave/WW3/analysis";
//        } else {
//            if(this.optim.getWaveModel() == 1){
//                Stagein_path = "inputFiles/wave/WW3/forecast";
//            } else {
//                if(this.optim.getWaveModel() >= 20){//relocatable
//                    Stagein_path = "inputFiles/wave/SWAN";
//                }
//            }
//        }
//
//        String wave_filename = Stagein_path+"/start__"+this.tGrid.getLatest_date()+".nc";

        if(this.paths.getForecastFile().substring(this.paths.getForecastFile().length()-3).equals("grb")){//if is grib
            try{
                waveForecastResults waveForecast =  MyGribParser.parseGrib(this.paths.getForecastFile(), this.paths.getOutDir());
                double[][][] VDIR = waveForecast.getVDIR();

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


                if(this.optim.getWaveModel() < this.optim.getRelocAnlsCode()){
                    for(int i=0;i<VDIR.length; ++i){
                        for(int j=0;j<VDIR[0].length;++j){
                            for(int k=0;k<VDIR[0][0].length; ++k){
                                VDIR[i][j][k] = VDIR[i][j][k] + (180*Math.signum((180.0-VDIR[i][j][k])));
                                if(VDIR[i][j][k]==180.0)
                                    VDIR[i][j][k] = 0.0;
                            }
                        }
                    }
                } else {
                    for(int i=0;i<VDIR.length; ++i){
                        for(int j=0;j<VDIR[0].length; ++j){
                            for(int k=0;k<VDIR[0][0].length; ++k){
                                double VY = Math.cos(this.constants.getDeg2rad()*VDIR[i][j][k]);
                                double VX = Math.sin(this.constants.getDeg2rad()*VDIR[i][j][k]);
                                VDIR[i][j][k] = Math.atan2(VY, VX)/this.constants.getDeg2rad();
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

                return new readout_mWaveResults(waveForecast.getLatitude(), waveForecast.getLongitude(), true, this.paths.getForecastFile(),
                        waveForecast.getVTDH(), waveForecast.getVTPK(), VDIR);
            } catch (Exception ex){
                return null;
            }
        } else {//NETCDF
            String wave_filename = this.paths.getForecastFile();
            MyNetCDFParser parser = new MyNetCDFParser(wave_filename);
            boolean wave_status = parser.isFileExists();

            if(wave_status){
                waveForecastResults waveForecast = parser.parseWaveForecastData(this.paths.getOutDir());

                double[][][] VDIR = waveForecast.getVDIR();

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


                if(this.optim.getWaveModel() < this.optim.getRelocAnlsCode()){
                    for(int i=0;i<VDIR.length; ++i){
                        for(int j=0;j<VDIR[0].length;++j){
                            for(int k=0;k<VDIR[0][0].length; ++k){
                                VDIR[i][j][k] = VDIR[i][j][k] + (180*Math.signum((180.0-VDIR[i][j][k])));
                                if(VDIR[i][j][k]==180.0)
                                    VDIR[i][j][k] = 0.0;
                            }
                        }
                    }
                } else {
                    for(int i=0;i<VDIR.length; ++i){
                        for(int j=0;j<VDIR[0].length; ++j){
                            for(int k=0;k<VDIR[0][0].length; ++k){
                                double VY = Math.cos(this.constants.getDeg2rad()*VDIR[i][j][k]);
                                double VX = Math.sin(this.constants.getDeg2rad()*VDIR[i][j][k]);
                                VDIR[i][j][k] = Math.atan2(VY, VX)/this.constants.getDeg2rad();
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

                return new readout_mWaveResults(waveForecast.getLatitude(), waveForecast.getLongitude(), wave_status, wave_filename,
                        waveForecast.getVTDH(), waveForecast.getVTPK(), VDIR);
            } else{
//            lat = new double[0];
//            lon = new double[0];
//            VTDH = new double[0][0][0];
//            VTPK = new double[0][0][0];
//            VDIR = new double[0][0][0];
                return new readout_mWaveResults(new double[0], new double[0], wave_status,
                        wave_filename, new double[0][0][0], new double[0][0][0], new double[0][0][0]);
            }
        }

    }


    private prepare_sqrtY_fieldResults prepare_sqrtY_field(ArrayList<Double> lat_bathy_Inset, ArrayList<Double> lon_bathy_Inset){
        // % preparation of a pseudo- waveheight field
        // % depending on sqrt of Y coordinate
        // % to be used for testing of analytical solution (=cycloid)

        //x, y
        double[] tmpLat = new double[1];
        tmpLat[0]=this.extreme_pts.getStart_lat();
        double[] tmpLon = new double[1];
        tmpLon[0]=this.extreme_pts.getStart_lon();

        deg2utmResults tmp = deg2utm(tmpLat, tmpLon);
        //double x_start = tmp.getX()[0];
        double y_start = tmp.getY()[0];
        String[] utmzone_start = tmp.getUtmzone()[0];

        tmpLat[0]=this.extreme_pts.getEnd_lat();
        tmpLon[0]=this.extreme_pts.getEnd_lon();

        tmp = deg2utm(tmpLat,tmpLon);
        double x_end = tmp.getX()[0];
        double y_end = tmp.getY()[0];
        //String[] utm_zone_end = tmp.getUtmzone()[0];

        int zoneN_start = Integer.parseInt(utmzone_start[0]);

        meshgridResults res = Utility.meshgrid(lat_bathy_Inset, lon_bathy_Inset);
        double[][] lat_gr = res.getX();
        double[][] lon_gr = res.getY();
        int Np = lat_gr.length * lat_gr[0].length;

        int[] utmzone_number = new int[Np];
        Arrays.fill(utmzone_number, zoneN_start);
        degzone2utmResults dz2uTmp = degzone2utm(Utility.reshape(lat_gr,Np, this.paths.getOutDir()), Utility.reshape(lon_gr, Np, this.paths.getOutDir()), utmzone_number);
        double[] xx=dz2uTmp.getX();
        double[] yy=dz2uTmp.getY();

        int Ny = lat_bathy_Inset.size();//43
        int Nx = lon_bathy_Inset.size();//72

        int[] dim = new int[2];
        dim[0]=Nx;
        dim[1]=Ny;
        double[][] Nyy = Utility.reshape(yy,dim, this.paths.getOutDir());

        //Dy
        double[][] DeltaY = new double[Nx][Ny];
        if(this.extreme_pts.getCycType() == "id"){//dd-type cycloid:
            for(int ix = 0;ix<Nx; ++ix){
                double[] y_diff = new double[Ny];
                for(int j =0;j<Ny; ++j){
                    y_diff[j] = Math.abs(Nyy[ix][j] - y_start); //m
                }
                minResults min = Utility.minWithIndex(y_diff);
                //double y_val = min.getElement();
                int y_idx = min.getIndex();
                for(int j=0;j<Ny; ++j){
                    DeltaY[ix][j] = Nyy[ix][y_idx] - Nyy[ix][j]; //m
                }
            }
        } else { //id-type cycloid:
            for(int ix = 0;ix<Nx; ++ix){
                double[] y_diff = new double[Ny];
                for(int j =0;j<Ny; ++j){
                    y_diff[j] = Math.abs(Nyy[ix][j] - y_end); //m
                }
                minResults min = Utility.minWithIndex(y_diff);
                //double y_val = min.getElement();
                int y_idx = min.getIndex();
                for(int j=0;j<Ny;j++){
                    DeltaY[ix][j] = Nyy[ix][j] - Nyy[ix][y_idx]; //m
                }
            }
        }
        for(int u =0 ;u<Nx; ++u){
            for(int v =0;v<Ny; ++v){
                if(DeltaY[u][v]<0){
                    DeltaY[u][v] = 0;
                }
            }
        }

        //speed
        this.extreme_pts.setPseudoG(0.001);
        double[][][] VTDH_b = new double[Ny][(int) this.tGrid.getNt()][Nx];
        double[][] DeltaYTransposed = Utility.transposeMatrix(DeltaY);
        for(int i=0;i<Ny; ++i){
            for(int it=0;it<(int) this.tGrid.getNt(); ++it){
                for(int j=0;j<Nx;j++){
                    VTDH_b[i][it][j] = this.constants.getMs2kts() * Math.sqrt(2*this.extreme_pts.getPseudoG()*DeltaYTransposed[i][j]);//kts
                }
            }
        }
        double[][][] VDIR_b = Utility.NaN3Dmatrix(Ny, (int) this.tGrid.getNt(), Nx);

        //graphical pars
        double smallFract = 0.5;
        this.fstats.setWheight_min((1-smallFract)*Utility.min3d(VTDH_b));
        this.fstats.setWheight_max((1+smallFract)*Utility.max3d(VTDH_b));
        return new prepare_sqrtY_fieldResults(VTDH_b, VDIR_b);
    }
    private Edges_definitionResults Edges_definition(double[][] xy, double[] xg_array, double[] yg_array, double[][][] VTDH_Inset, double[][][] VTPK_Inset,
                                  double[][][] VDIR_Inset, double[][][] windMAGN_Inset, double[][][] windDIR_Inset, double[][] bathy_Inset,
                                  double[][] J_mask, double[][] ship_v_LUT, ArrayList<Double>... varargin){

        /* readout DB:
         *  remapping free edges:
         *  pointers:
         *  edge lengths and angles:
         *  edge weights:
         *  edge delays:
         *  */
        boolean motorboat = false;
        if(varargin.length == 1) {
            motorboat = true; //so, varargin = H_array_m
            //else sailboat, so varargin(0) = twa_array, varargin(1) = tws_array
        }
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("Defining edges...\n");
        this.logFile.WriteLog("\tReading out free edges from DB...");
        this.logFile.CloseFile();
        MyBinaryParser edgesDBParser = new MyBinaryParser(this.paths.getFreeEdgesDB());
        int[] free_edges_DB_array = edgesDBParser.readAsUInt16();
        int[] free_edges_DB_size = new int[2];
        free_edges_DB_size[0]=(int) Math.floor(free_edges_DB_array.length/2);
        free_edges_DB_size[1]=2;
        int[][] free_edges_DB = Utility.reshape(free_edges_DB_array,free_edges_DB_size, this.paths.getOutDir());

        //remapping free edges:
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("\tremapping free edges to inset Sgrid...");
        this.logFile.CloseFile();
        long ll = Math.round((this.sGrid.getXi()-this.sGrid.getDB_xi())*this.sGrid.getInv_step());
        long mm =  Math.round((this.sGrid.getYi()-this.sGrid.getDB_yi())*this.sGrid.getInv_step());
        int[][] free_edges = idx_ref2inset_gridCompact(free_edges_DB, this.sGrid.getDB_Nx(), this.sGrid.getDB_Ny(), this.sGrid.getInset_Nx(), this.sGrid.getInset_Ny(), ll, mm);
        int free_edges_Number = free_edges.length;
        if(free_edges_Number==0){
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("Edges_definition: no free edges found in graph!");
            debug.CloseFile();
            System.exit(0);
        }
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("\tnumber of free nodes: "+this.sGrid.getFreenodes());
        this.logFile.WriteLog("\tnumber of free edges: "+free_edges_Number);
        this.logFile.CloseFile();
        doPointerResults ptrResults = doPointer(free_edges);
        edge__lenghts_anglesResults edgeResults = edge__lenghts_angles(xy,free_edges);
        //edge weights:
        fields_node2edgeResults fields_node2edgeRes = fields_node2edge(free_edges,VTDH_Inset,VTPK_Inset,
                MemorySaver.getYInset(VDIR_Inset, this.constants.getDeg2rad()), MemorySaver.getXInset(VDIR_Inset, this.constants.getDeg2rad()),
                bathy_Inset, J_mask, windMAGN_Inset, MemorySaver.getYInset(windDIR_Inset, this.constants.getDeg2rad()), MemorySaver.getXInset(windDIR_Inset, this.constants.getDeg2rad()));

        double[][] waveDir_edges = changeDirRule( Utility.atan2(fields_node2edgeRes.getYwave_edges(), fields_node2edgeRes.getXwave_edges(), this.paths.getOutDir()) );
        for(int i=0;i<waveDir_edges.length;++i)
            for(int j=0;j<waveDir_edges[0].length;++j)
                waveDir_edges[i][j] = waveDir_edges[i][j]/this.constants.getDeg2rad();

        double[][] windDir_edges = changeDirRule( Utility.atan2(fields_node2edgeRes.getYwind_edges(), fields_node2edgeRes.getXwind_edges(), this.paths.getOutDir()) );
        for(int i=0;i<windDir_edges.length;++i)
            for(int j=0;j<windDir_edges[0].length;++j)
                windDir_edges[i][j] = windDir_edges[i][j]/this.constants.getDeg2rad();

        ArrayList<Integer> nogo_edges = Utility.findNaNs(fields_node2edgeRes.getJ_edges());

        //edge delays:
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("\tcomputing edge delays...");
        this.logFile.CloseFile();
        //edge delays
        if(this.ship.getVessType() != this.ship.getSailType()){//motorboats
            edge_delaysResults edge_delays = edge_delays(free_edges_Number, edgeResults.getTheta_grid(), edgeResults.getEdge_lenght(),
                    fields_node2edgeRes.getWaveHeight_edges(), fields_node2edgeRes.getWavePeriod_edges(), waveDir_edges,
                    fields_node2edgeRes.getWindMAGN_edges(), windDir_edges, nogo_edges, fields_node2edgeRes.getBathy_edges(),
                    ship_v_LUT, varargin[0]);
            return new Edges_definitionResults(free_edges, nogo_edges, edgeResults.getTheta_grid(), edgeResults.getEdge_lenght(),
                    edge_delays.getSh_delay(), edge_delays.getVarargout().get(0), edge_delays.getSafe_indexes(), edge_delays.getTdep_danger_idx(),
                    fields_node2edgeRes.getWaveHeight_edges(), fields_node2edgeRes.getWavePeriod_edges(), fields_node2edgeRes.getWaveLenght_edges(),
                    waveDir_edges, fields_node2edgeRes.getWindMAGN_edges(), windDir_edges, fields_node2edgeRes.getBathy_edges(),
                    ptrResults.getHead_bool(), ptrResults.getHead_ord(), ptrResults.getPointer());
        } else { //sailboats
            /**TODO: Sailboats*/
            return null;
        }
    }

    private edge_delaysResults edge_delays(int free_edges_Number,double[] theta_grid, double[] edge_length, double[][] waveHeight_edges, double[][] wavePeriod_edges, double[][] waveDir_edges,
                             double[][] windMAGN_edges, double[][] windDir_edges, ArrayList<Integer> nogo_edges,
                             double[] bathy_edges, double[][] ship_v_LUT, ArrayList<Double>... varargin){
        //%
        //% computation of graph edge delays a_{jk} :
        //%
        //%   a_{jk}= |r_k - r_j| / v(P, H_{jk})
        //%
        //% and definition of safe edge indexes.
        //%
        //%
        //          if nargin==nMin
        //              H_array_m= varargin{1};
        //          else
        //              twa_array= varargin{1};
        //              tws_array= varargin{2};
        //          end
        //ArrayList<DangerIndexes> tdep_danger_idx = new ArrayList<>();
        edge_delaysResults out = new edge_delaysResults();
        double[] H_array_m = new double[0];
        //double[] twa_array = new double[0];
        //double[] tws_array = new double[0];
        if(varargin.length==1){
            //motorboat
            H_array_m = Utility.arrayfy(varargin[0]);
        } else {
            //sailboat
            //twa_array = Utility.arrayfy(varargin[0]);
            //tws_array = Utility.arrayfy(varargin[1]);
        }
        double[][] sh_delay = new double[free_edges_Number][(int) this.tGrid.getNt()];
        int[][] gear_idx = new int[free_edges_Number][(int) this.tGrid.getNt()];

        int[][][] safe_indexes;
        if((int) this.visualization.getSafegram() == 1){
            safe_indexes = new int[(int) this.tGrid.getNt()][(int) this.ship.getNvel()][free_edges_Number];
            Arrays.fill(safe_indexes,Double.NaN);
            //safe_indexes = Utility.NaN((int) this.ship.getNvel(), free_edges_Number, (int) this.tGrid.getNt());
        } else {
            safe_indexes = new int[1][1][1];
            safe_indexes[0][0][0] = Integer.MIN_VALUE;//Instead of nan
        }
        //etc

        //-------------------------------------------------------------------------
        //PARAMETERS:
        //surfriding and pure loss of stability pars after Belenky et.al. 2010:
        double steep_min1 = 1.0/40.0;
        double steep_max1 = 1.0/5.0;
        double steep_min2 = 1.0/25.0;
        double steep_min3 = 1.0/20.0;
        this.ship.setMin_lambda_L(0.8);
        this.ship.setMax_lambda_L(2.0);
        double min_vTw = 1.3;
        double max_vTw = 2.0;
        //rolling resonance linewidth parameters after:
        //Benedict,Baldauf,Kirchhoff (Hochschule Wismar - University of Technology,
        //Business and Design; Department of Maritime Studies)
        //"Development of ARROW software for decision support to avoid roll resonance and wave impact for ship operation in heavy seas"
        double ONE_inf = 0.8;
        double ONE_sup = 1.1;
        double TWO_sup = 2.1;
        double TWO_inf = 1.8;
        //double TWO_sup = 2.1;
        //cp. MAN propulsion manual:
        //LPP = 0.97 x LWL
        double fncrit1=0.2324;
        double fncrit2=0.07364;
        //-------------------------------------------------------------------------
        double[][][] sh_delay_gear = new double[(int) this.tGrid.getNt()][(int) this.ship.getNvel()][edge_length.length];
        boolean motorboat = false;
        if(this.ship.getVessType() < this.ship.getSailType()){ //motorboats
            motorboat=true;
//            if forcing.analytic
//                    sh_vel = waveHeight_edges(:,1);
//            end
            do_Fr_crit_LUTResults do_fr_crit_lutRes = do_Fr_crit_LUT(fncrit1, fncrit2, steep_min1, steep_max1);

            for(int it=0;it<(int) this.tGrid.getNt(); ++it){

                //wave direction
                DangerIndexes danger_idx = new DangerIndexes();
                double[] alpha = wave2ship_reldir(theta_grid, waveDir_edges, it);
                boolean[] follseas_bol = new boolean[alpha.length];
                double[] IMO_secans = new double[alpha.length];
                for(int i=0;i<alpha.length;++i){
                    double tmp = Math.abs(Math.abs(alpha[i])-180);
                    if(tmp <= 45)
                        follseas_bol[i] = true;
                    else
                        follseas_bol[i] = false;
                    tmp = 1.0/Math.abs(Math.cos(this.constants.getDeg2rad()*(180.0-alpha[i])));
                    IMO_secans[i] = Math.min(Math.sqrt(2.0), tmp);
                }

                //wave period
                double[] waveLenght_edges = wave_dispersion(wavePeriod_edges, it, bathy_edges);
                boolean[] lambda_bol = new boolean[waveLenght_edges.length];
                boolean[] lambda_bol2 = new boolean[waveLenght_edges.length];
                for(int i=0;i<lambda_bol.length; ++i){
                    double lambda_L = waveLenght_edges[i]/this.ship.getLength();
                    if(this.ship.getMin_lambda_L() <= lambda_L){
                        lambda_bol2[i] = true;
                        if(lambda_L <= this.ship.getMax_lambda_L())
                            lambda_bol[i] = true;
                        else
                            lambda_bol[i] = false;
                    } else {
                        lambda_bol2[i] = false;
                    }
                }

                //wave height
                boolean[] HL_bol2 = new boolean[waveHeight_edges.length];
                boolean[] HL_bol3 = new boolean[waveHeight_edges.length];
                double[] steepness = new double[waveLenght_edges.length];
                for(int i=0;i<waveHeight_edges.length; ++i){
                    double tmp = waveHeight_edges[i][it]/this.ship.getLength();
                    if(tmp >= steep_min2)
                        HL_bol2[i]=true;
                    else
                        HL_bol2[i]=false;
                    if(tmp >= steep_min3)
                        HL_bol3[i]=true;
                    else
                        HL_bol3[i]=false;

                    steepness[i]=waveHeight_edges[i][it]/waveLenght_edges[i];
                }

                double[] Fr_crit = Utility.interp1(do_fr_crit_lutRes.getSteep_var(), do_fr_crit_lutRes.getFr_LUT(), steepness, Double.POSITIVE_INFINITY, this.paths.getOutDir());
                //-------------------------------------------------------------------------------------------------
                //ship speed and ship edge delays:
                for( int jv = 0; jv<(int) this.ship.getNvel(); ++jv){
                    double[] sh_vel;
                    if(this.forcing.getAnalytic() != 1){//FALSO
                        sh_vel = Utility.interp1(H_array_m, ship_v_LUT, jv, waveHeight_edges, it, this.paths.getOutDir());
                    } else {
                        sh_vel = new double[waveHeight_edges.length];
                        for(int i=0;i<waveHeight_edges.length;++i)
                            sh_vel[i]=waveHeight_edges[i][0];
                    }
                    for(int j=0;j<edge_length.length;++j)
                        sh_delay_gear[it][jv][j] = edge_length[j]/sh_vel[j];
                    //-------------------------------------------------------------------------------------------------
                    // SAFETY constraint #1) parametric rolling
                    double[] absT_E = enc_wavePeriod(wavePeriod_edges, it, bathy_edges, sh_vel, alpha);
                    for(int i=0;i<absT_E.length;++i)
                        absT_E[i] = Math.abs(absT_E[i]);

                    boolean[] parRoll_bol = new boolean[absT_E.length];
                    boolean[] syncRoll_bol = new boolean[absT_E.length];
                    boolean[] res2_bol = new boolean[absT_E.length];
                    boolean[] res1_bol = new boolean[absT_E.length];

                    for(int i=0;i<parRoll_bol.length;++i){

                        if(absT_E[i] * ONE_inf <= this.ship.getRoll_period() && this.ship.getRoll_period() <= absT_E[i] * ONE_sup)
                            parRoll_bol[i] = true;
                        else
                            parRoll_bol[i]= false;

                        if(absT_E[i] * TWO_inf <= this.ship.getRoll_period() && this.ship.getRoll_period() <= absT_E[i]* TWO_sup)
                            syncRoll_bol[i] = true;
                        else
                            syncRoll_bol[i] = false;

                        if(lambda_bol[i] && HL_bol3[i] && parRoll_bol[i])
                            res2_bol[i] = true;
                        else
                            res2_bol[i] = false;

                        if(lambda_bol[i] && HL_bol3[i] && syncRoll_bol[i])
                            res1_bol[i] = true;
                        else
                            res1_bol[i] = false;

                    }

                    danger_idx.setResonance2(Utility.find(res2_bol));
                    danger_idx.setResonance1(Utility.find(res1_bol));
//                    -------------------------------------------------------------------------------------------------
//                    SAFETY constraint #2) pure loss of stability
                    boolean[] vTw_bol = new boolean[sh_vel.length];
                    boolean[] pureLossStab_bol = new boolean[sh_vel.length];
                    for(int i=0;i<vTw_bol.length;++i){
                        if(sh_vel[i] >= min_vTw * IMO_secans[i] * wavePeriod_edges[i][it] && sh_vel[i] <= max_vTw * IMO_secans[i] * wavePeriod_edges[i][it])
                            vTw_bol[i]=true;
                        else
                            vTw_bol[i] = false;
                        if(lambda_bol2[i] && HL_bol2[i] && follseas_bol[i] && vTw_bol[i])
                            pureLossStab_bol[i] = true;
                        else
                            pureLossStab_bol[i] = false;
                    }

                    danger_idx.setPureLossStab(Utility.find(pureLossStab_bol));

                    //-------------------------------------------------------------------------------------------------
                    //             SAFETY constraint #3) surfriding/broaching-to

                    double[] Fr = new double[sh_vel.length];
                    boolean[] Fr_bol = new boolean[sh_vel.length];
                    boolean[] surfRid_bol = new boolean[sh_vel.length];
                    for(int i=0;i<Fr.length;++i){
                        Fr[i] = sh_vel[i]/(this.constants.getMs2kts() * Math.sqrt(this.constants.getG0()*this.ship.getLength()));
                        if(Fr[i] >= Fr_crit[i] * IMO_secans[i])
                            Fr_bol[i] = true;
                        else
                            Fr_bol[i] = false;
                        if(lambda_bol[i] && follseas_bol[i] && Fr_bol[i])
                            surfRid_bol[i] = true;
                        else
                            surfRid_bol[i] = false;
                    }

                    danger_idx.setSurfRiding(Utility.find(surfRid_bol));


                    //-------------------------------------------------------------------------------------------------

                    if(jv==0)
                        out.addDangerIdx(danger_idx);
                        //tdep_danger_idx.add(danger_idx);

                    if (this.safety.getCriteria().get(0)==0){//parametric rolling means checking for both Te=Tr and 2*Te=Tr conditions:
                        danger_idx.setResonance1(new ArrayList<>());
                        danger_idx.setResonance2(new ArrayList<>());
                    }
                    if(this.safety.getCriteria().get(1)==0)
                        danger_idx.setPureLossStab(new ArrayList<>());
                    if(this.safety.getCriteria().get(2)==0)
                        danger_idx.setSurfRiding(new ArrayList<>());



                    //setting  edge_weight=Inf for edges with dangerous conditions:
                    ArrayList<Integer> danger_indexes = new ArrayList<>();
                    danger_indexes.addAll(danger_idx.getResonance1());
                    danger_indexes.addAll(danger_idx.getResonance2());
                    danger_indexes.addAll(danger_idx.getPureLossStab());
                    danger_indexes.addAll(danger_idx.getSurfRiding());
                    danger_indexes.addAll(nogo_edges);

                    //setting  edge_weight=Inf for edges with dangerous conditions:
                    int[] all_indexes = new int[free_edges_Number];
                    for(int i=0; i<free_edges_Number; ++i)
                        all_indexes[i] = i;
                    for(int index : danger_indexes) {
                        sh_delay_gear[it][jv][index] = Double.POSITIVE_INFINITY;
                        all_indexes[index] = Integer.MIN_VALUE; //used instead of NaN
                    }

                    if(this.visualization.getSafegram() == 1){
                        for(int j=0;j<all_indexes.length;++j){
                            safe_indexes[it][jv][j] = all_indexes[j]; //to be used for adjacency plot - danger edges are removed from graph!
                        }
                    }

                }//Nvel
                //following finds optimal delay by intentional speed reduction:

                //search min along dimension of Nvel
                for(int j=0;j<sh_delay_gear[0][0].length;++j){
                    double minVal = Double.MAX_VALUE;
                    int min_idx = -1;
                    for(int i=0;i<sh_delay_gear[0].length;++i){
                        if(!Double.isNaN(sh_delay_gear[it][i][j])){
                            if(sh_delay_gear[it][i][j]<minVal){
                                minVal = sh_delay_gear[it][i][j];
                                min_idx = j;
                            }
                        }
                    }
                    if(min_idx==-1){//All NaNs along Nvel dimension
                        min_idx = 0;
                        minVal = Double.NaN;
                    }
                    sh_delay[j][it] = minVal;
                    gear_idx[j][it] = min_idx;
                }

            }//tgrid.nt
        } else {//sailboat
            /** TODO: SAILBOAT*/

        }
        //edge_delaysResults out = new edge_delaysResults(sh_delay,safe_indexes,tdep_danger_idx);
        out.setSh_delay(sh_delay);
        out.setSafe_indexes(safe_indexes);
        if(motorboat){
            out.AddVarargout(gear_idx);
        }
        return out;
    }

    private double[] enc_wavePeriod(double[][] wave_period, int wpIndex, double[] depth, double[] v_ship, double[] alpha){
//        %   encounter wave period relation for monocromatic ocean waves
//                %     source: IMO 1228, sect 1.6 (holds just in deep waters

        if((v_ship.length != wave_period.length) || (alpha.length != wave_period.length)){
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("enc_wavePeriod: all arguments must have same dimensions!");
            debug.CloseFile();
            System.exit(0);
        }

        double cc = this.constants.getMs2kts()*this.constants.getG0()/(2*Math.PI);
        double[] v0 = new double[wave_period.length];
        double[] T_E = new double[wave_period.length];
        double[] fmk = null;
        if(this.forcing.getDeepWaterApprox() == 1)
            fmk=Fenton_McKee_factor(wave_period,wpIndex,depth);
        for(int i=0;i<wave_period.length;++i){
            v0[i] = cc * wave_period[i][wpIndex];
            if(this.forcing.getDeepWaterApprox() != 1)
                v0[i] *= fmk[i];
        }
//        % correct, general purpose formula:
//      % WARNING: depth may be NaN in some parts of the bbox!
        for(int i=0;i<T_E.length;++i)
            T_E[i] = wave_period[i][wpIndex]/( 1.0 + v_ship[i]/v0[i]* Math.cos(this.constants.getDeg2rad() * alpha[i] )  );
        return T_E;
    }

    private double[] wave2ship_reldir(double[] ship_dir, double[][] envField_dir, int colIndex){
//         % angles are expressed in degrees (!!):
//          % 0  = head sea
//           % 90 = beam sea
//            % 180= following sea
//          % in agreement with alpha in IMO circ. no 1228
//            % and sailboat PolarPlots convention

        if(ship_dir.length != envField_dir.length){
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("wave2ship_reldir: ship_dir and envField_dir have different sizes.");
            debug.CloseFile();
            System.exit(0);
        }

        double[] alpha = new double[ship_dir.length];
        for(int i=0;i<alpha.length;++i){
            double delta = envField_dir[i][colIndex] - ship_dir[i];
            if(delta==0.0)
                alpha[i] = 180.0;
            else
                alpha[i] = delta - (180.0* (double) Utility.sign(delta));
        }
        return alpha;
    }

    private do_Fr_crit_LUTResults do_Fr_crit_LUT(double fncrit1, double fncrit2, double steep_min, double steep_max){
        //         % LUT of critical Fr for surfriding/broaching-to according to Belenky
        //         % et.al. 2010 (cp. Fig.4.22)
        //         % xvar is wave steepness: H/lambda
        int nx=20;
        ArrayList<Double> xVarTmp = Utility.linspace(steep_min, steep_max, nx);
        double[] xvar = new double[xVarTmp.size()];
        double[] inv_xvar = new double[xVarTmp.size()];
        for(int i=0;i<xvar.length;i++){
            xvar[i] = xVarTmp.get(i);
            inv_xvar[i] = 1/xvar[i];
        }

        double[] Fr_LUT = new double[inv_xvar.length];
        for(int i=0;i<inv_xvar.length;i++)
            Fr_LUT[i] = fncrit1 * Utility.nthroot(inv_xvar[i], 3, this.paths.getOutDir()) - fncrit2 * Math.sqrt(inv_xvar[i]);

        return new do_Fr_crit_LUTResults(xvar, Fr_LUT);
    }

    private doPointerResults doPointer(int[][] head){
        //Considering only first column
//      % (values of "head" array must be nonnull integer numbers)
//      % head_bool: logical index of j-th node  (1/0: present/absent in "head" array)
//      % head_ord: ordinal position of head nodes in head array (may contain repetions if not
//                       % all head_bool = 1)
//      % pointer: first position of j-th node in head array

        int Nhead = head.length;
        int[] pointer = new int[Nhead+1];
        for(int i=0;i<pointer.length;i++){
            pointer[i]=-1;
        }

        //initialization
        int ip=0;
        pointer[ip]=0;
        int ic=0; // generic counter
        int head_val = head[ic][0];

        //iteration:
        while(ic<(Nhead-1)){//length of array to be pointicized
            ic++;
            if(head[ic][0] > head_val){
                ip++;
                pointer[ip] = ic;
                head_val = head[ic][0];
            }
        }
        pointer[Nhead]=Nhead;

        //post-proc:
        pointer = Utility.removeElements(pointer,-1);
        //int[] pointer_short = Utility.deepCopy(pointer, -1);
        int[] head_vals = getElements(head,Utility.deepCopy(pointer, -1));
        int max = Utility.max(head,0);
        int[] head_argm = new int[max];
        int counter=0;
        for(int i=0;i<max;i++){
            head_argm[i]=counter;
            counter++;
        }

        boolean[] head_bool = Utility.ismember(head_argm, head_vals);
        int[] head_ord = Utility.cumsum(head_bool);
        return new doPointerResults(head_bool, head_ord, pointer);
    }

    private int[] getElements(int[][] array, int[] indexes){
        int[] out = new int[indexes.length];
        for(int i=0;i<out.length;i++)
            out[i]=array[indexes[i]][0];
        return out;
    }

    private edge__lenghts_anglesResults edge__lenghts_angles(double[][] xy, int[][] free_edges){
        double[] xm = new double[free_edges.length];
        double[] ym = new double[free_edges.length];
        for(int i=0;i<free_edges.length;i++){
            //              P_1(:,1)                    P_2(:,1)
            xm[i] = (xy[ free_edges[i][0] ][0] + xy[ free_edges[i][1] ][0])/2;
            //              P_1(:,2)                    P_2(:,2)
            ym[i] = (xy[ free_edges[i][0] ][1] + xy[ free_edges[i][1] ][1])/2;
        }
        //dipending on number og neighbors, we can call dir_on_grid(this.sGrid.numberOfNeighbors); for now, we have numberOfNeighbors=24
        double[] theta_grid = dir_on_grid24n(free_edges);
        //double[] edge_lenght = hor_distanceForEdges("s",xy,free_edges);
        double[] edge_lenght = Haversine_distance(xy, free_edges);
        return new edge__lenghts_anglesResults(theta_grid,edge_lenght,xm,ym);
    }

    private fields_node2edgeResults fields_node2edge(int[][] free_edges, double[][][] VTDH_Inset, double[][][] VTPK_Inset, double[][][] Ywave_Inset, double[][][] Xwave_Inset,
                                  double[][] bathy_Inset, double[][] J_mask, double[][][]... varargin){
        //varargin[0] = windMAGN_Inset
        //varargin[1] = Ywind_Inset
        //varargin[2] = Xwind_Inset
        this.logFile = new MyFileWriter("","",true, this.paths.getOutDir());
        this.logFile.WriteLog("\tcomputing edge weights from model data");
        this.logFile.CloseFile();
        fields_node2edgeResults res = new fields_node2edgeResults();
        //Bathy and Joint-mask:
        //double[] bathy_edges = mStaticF_EWeights(free_edges, bathy_Inset,"min").getMin();
        res.setBathy_edges(mStaticF_EWeights(free_edges, bathy_Inset, "min"));
        //double J_edges = mStaticF_EWeights(free_edges, J_mask, "mean").getMean();
        res.setJ_edges(mStaticF_EWeights(free_edges, J_mask, "mean"));

        ArrayList<double[][]> out;
        if(varargin != null){//wave and wind fields
            out = mDynamicF_EWeights(free_edges, VTDH_Inset, VTPK_Inset, Ywave_Inset, Xwave_Inset, varargin[0], varargin[1], varargin[2]);
            res.setValues(out);
            //[waveHeight_edges, wavePeriod_edges, Ywave_edges,Xwave_edges, windMAGN_edges, Ywind_edges,Xwind_edges]
//            varargout{1}= windMAGN_edges;
//            varargout{2}= Ywind_edges;
//            varargout{3}= Xwind_edges;

            //lambda
            //double[][] waveLenght_edges = Utility.NaN2Dmatrix(out.get(1).length, out.get(1)[0].length);
            double[][] waveLenght_edges = Utility.NaN2Dmatrix(res.getWavePeriod_edges().length, res.getWavePeriod_edges()[0].length);
            for(int it=0; it<(int)this.tGrid.getNt();++it){
                //waveLength_edges(:,it) = wave_dispersion(wavePeriod_edges(:,it), bathy_edges);
                //double[] waveDisp = wave_dispersion(out.get(1),it,res.getBathy_edges());
                double[] waveDisp = wave_dispersion(res.getWavePeriod_edges(),it,res.getBathy_edges());
                for(int i=0;i<waveLenght_edges.length;++i){
                    waveLenght_edges[i][it]=waveDisp[i];
                }
            }
            res.setWaveLenght_edges(waveLenght_edges);
        } else { //just wave fields (for case forcing.analytic)
            out = mDynamicF_EWeights(free_edges, VTDH_Inset, Ywave_Inset, Xwave_Inset);
            res.setValues(out);
            res.setWavePeriod_edges(Utility.NaN2Dmatrix(res.getWaveHeight_edges().length,res.getWaveHeight_edges()[0].length));
            res.setWaveLenght_edges(Utility.NaN2Dmatrix(res.getWaveHeight_edges().length,res.getWaveHeight_edges()[0].length));
//            % resting fictious fields:
//            wavePeriod_edges = NaN*ones(size(waveHeight_edges));
//            waveLength_edges = NaN*ones(size(waveHeight_edges));
        }
        return res;
    }

    private double[] wave_dispersion(double[][] wave_period, int colIndex, double[] depth){

//     dispersion relation for monocromatic ocean waves
//                %
//%
//%   -) If just wave period provided --> deep water approximation
//%
//%   -) If also depth provided       --> generic depth formula after:
//%   Fenton, JD and McKee, WD (1990)
//                %
//                %   wave_period [s]
//                %   depth       [m]
//                %   lambda      [m]

        double twopi = 2*Math.PI;
        //deep water approximation
        double[] lambda = new double[wave_period.length];
        for(int i=0;i<lambda.length;i++){
            lambda[i]=this.constants.getG0()/twopi*Math.pow(wave_period[i][colIndex],2);
        }

        if(this.forcing.getDeepWaterApprox()!=1){
            if(depth!=null){
                //generic depth formula
                double[] fmk = Fenton_McKee_factor(wave_period, colIndex, depth);
                for(int i=0;i<lambda.length;i++){
                    lambda[i]=lambda[i]*fmk[i];
                }
            }
        }
        return lambda;
    }

    private ArrayList<double[][]> mDynamicF_EWeights(int[][] free_edges, double[][][]... varargin){

        ArrayList<double[][]> vargout = new ArrayList<>();
        for(int i=0;i<varargin.length;++i){
            int[] dim = {(int)this.tGrid.getNt(),(int) (Utility.numel(varargin[i])/(int)this.tGrid.getNt())};
            //condense the spatial dimensions into a single dimension:
            double[][] model_field = Utility.reshape(Utility.permute3DColsWithZdim(varargin[i]), dim);
            varargin[i]=null;
            //function max(a,b) outputs the not-NaN value among a and b, if existing:
            double[][] Weights = new double[free_edges.length][(int) this.tGrid.getNt()];
            for(int it=0;it<(int) this.tGrid.getNt(); ++it){
//                phi_I   = model_field(it,free_I_nodes);
//                phi_J   = model_field(it,free_J_nodes);
//                phi_IJ2 = ( phi_I + phi_J ) /2.;
//
//                Weights(:,it)= phi_IJ2; % to be used with sea-over-land

                for(int row=0;row<free_edges.length; ++row){
                    Weights[row][it] = (model_field[it][free_edges[row][0]] + model_field[it][free_edges[row][1]]) / 2.0;
                }
            }
            vargout.add(Weights);
        }
        return vargout;
    }

    private double[] mStaticF_EWeights(int[][] free_edges, double[][] input_field, String operation){
        //condense the spatial dimensions into a single dimension:
        double[] model_field=Utility.reshape(Utility.transposeMatrix(input_field), Utility.numel(input_field), this.paths.getOutDir());
        //conservative approach: minimum edge depth
        double[] phi_I = new double[free_edges.length];
        double[] phi_J = new double[free_edges.length];
        for(int i=0;i<free_edges.length;++i){
            //QUI AVEVO MESSO -1
            phi_I[i]=model_field[free_edges[i][0]];
            phi_J[i]=model_field[free_edges[i][1]];
        }
        double[] output_field = null;
        switch (operation){
            case "min":
                output_field = Utility.min(phi_I, phi_J, this.paths.getOutDir());
                break;
            case "mean":
                //out_field = Utility.mean(phi_I, phi_J);
                output_field = Utility.mean(phi_I, phi_J, this.paths.getOutDir());
                break;
                default:
                    MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
                    debug.WriteLog("mStaticF_EWeights: unknown operation!");
                    debug.CloseFile();
                    System.exit(0);
                    break;
        }

        return output_field;
    }

    private double[] dir_on_grid24n(int[][] free_edges){
        double[] theta = new double[free_edges.length];
        double theta_27= Math.atan(1.0/2.0)/this.constants.getDeg2rad();
        int Nx =(int) this.sGrid.getInset_Nx();
        for(int ie=0;ie<free_edges.length;++ie){
            int J_I = free_edges[ie][1] - free_edges[ie][0];
            double angolo=0;
            if((J_I == Nx) || J_I==(2*Nx))
                angolo=0.0;
            else{
                if((J_I == (Nx+1)) || J_I == (2*Nx+2))
                    angolo=45.0;
                else{
                    if((J_I == 1) || (J_I == 2))
                        angolo=90.0;
                    else{
                        if((J_I == (-Nx+1)) || (J_I == ((-2*Nx)+2)))
                            angolo=135.0;
                        else{
                            if((J_I == (-Nx)) || J_I == (-2*Nx))
                                angolo=180.0;
                            else{
                                if((J_I == (-Nx-1)) || J_I == (-2*Nx)-2)
                                    angolo=225.0;
                                else{
                                    if((J_I == (-1)) || J_I == -2)
                                        angolo=270.0;
                                    else{
                                        if((J_I == (Nx-1)) || J_I == (2*Nx)-2)
                                            angolo=315.0;
                                        else{
                                            if(J_I == (2*Nx+1))
                                                angolo=theta_27;
                                            else{
                                                if(J_I == Nx+2)
                                                    angolo=90-theta_27;
                                                else{
                                                    if(J_I == (-Nx + 2))
                                                        angolo=90+theta_27;
                                                    else{
                                                        if(J_I==(-2*Nx + 1))
                                                            angolo=180.0-theta_27;
                                                        else{
                                                            if(J_I == (-2*Nx) -1)
                                                                angolo=180.0+theta_27;
                                                            else{
                                                                if(J_I == -Nx-2)
                                                                    angolo=270.0-theta_27;
                                                                else{
                                                                    if(J_I == -2+Nx)
                                                                        angolo = 270.0+theta_27;
                                                                    else{
                                                                        if(J_I == (2*Nx)-1)
                                                                            angolo=360.0-theta_27;
                                                                        else{
                                                                            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
                                                                            debug.WriteLog("dirOnGrid24N: I and J ("+free_edges[ie][0]+" and "+free_edges[ie][1]+") are not next or second next neighbours of a squared grid with "+Nx+" horizontal nodes!");
                                                                            debug.CloseFile();
                                                                            System.exit(0);
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            theta[ie]=angolo;
        }
        return theta;
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
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
            debug.WriteLog("degzone2utm: Lat and Lon vectors should have the same length");
            debug.CloseFile();
            System.exit(0);
        } else{
            if(n1!=n3){
                System.out.println("degzone2utm: Lat and Lon vectors should have the same length of utmzone_number vector");
                MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
                debug.WriteLog("degzone2utm: Lat and Lon vectors should have the same length of utmzone_number vector");
                debug.CloseFile();
                System.exit(0);
            }
        }

        //Memory pre-allocation
        double[] x = new double[n1];
        double[] y = new double[n1];

        //Main loop
        for(int i=0;i<n1;++i){
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
            MyFileWriter debug = new MyFileWriter("","debug",false, this.paths.getOutDir());
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

}
