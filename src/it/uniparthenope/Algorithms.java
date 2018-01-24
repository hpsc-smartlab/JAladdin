package it.uniparthenope;

import com.google.common.primitives.Ints;
import com.mchange.v1.util.ArrayUtils;
import it.uniparthenope.Boxing.Dijkstra2DResults;

import java.util.*;

public class Algorithms {

    public static Dijkstra2DResults Dijkstra(int Nnodes, int[][] free_edges, double[] edge_costs, boolean[] I_bool, int[] I_ord, int[] I_point, long SID, long FID){
        /*
        * nodes = lon and lat for each node in graph
        * free_edges = free edges
        * edge_cost = weights
        * I_bool, I_ord, I_point = for "forward star" data structure
        * SID = Starting node
        * FID = Ending node
        * */

        //Initializations:
        double costs = 0;
        ArrayList<Integer> paths = new ArrayList<>();
        //int Nnodes = nodes.length;
        HashMap<Integer, Double> TBL = new HashMap<>(); // Used as a sparse array. First is index, second is value
        double[] min_cost = new double[Nnodes];
        Arrays.fill(min_cost, Double.POSITIVE_INFINITY);
        boolean[] settled = new boolean[Nnodes];
        Arrays.fill(settled, false);
        ArrayList<Integer>[] path = new ArrayList[Nnodes];
        int I = (int) SID;
        min_cost[I] = 0.0;
        TBL.put(I,0.0);
        settled[I] = false;
        if(path[I]==null)
            path[I] = new ArrayList<>();
        path[I].add(I);
        boolean isEmpty = false;

        while(!settled[(int) FID] || !isEmpty){
            HashMap<Integer, Double> TAB = (HashMap<Integer, Double>) TBL.clone();
            TBL.put(I, 0.0);

            int[] nids = {0};
            if(I_bool[I]){
                int head_ordinal = I_ord[I];
                int block_start = I_point[head_ordinal];
                int block_end = I_point[head_ordinal +1] -1;
                nids = new int[(block_end-block_start)+1];
                int tmp = block_start;
                for(int i=0;i<nids.length;++i) {
                    nids[i] = tmp;
                    tmp++;
                }
            }

            //Calculate the Costs to the Neighbor Points and Record Paths
            for(int kk=0; kk<nids.length;++kk){
                int row = nids[kk];
                int J = free_edges[row][1];
                if(!settled[J]){
                   double edge_cost = edge_costs[row];
                   boolean empty = false;
                   if(TAB.get(J) == null){//value not found
                       empty = true;
                   }
                   if( empty || (TAB.get(J) > (TAB.get(I) + edge_cost)) ){
                       TAB.put(J, TAB.get(I) + edge_cost);
                       if(path[J]==null)
                           path[J] = new ArrayList<>();
                       path[J].addAll(path[I]);
                       path[J].add(J);
                   } else {
                       TBL.put(J, TAB.get(J));
                   }
                }
            }

            //Set<Integer> keys = TBL.keySet();
            //int[] K = keys.stream().mapToInt(Integer::intValue).toArray();
            //Find the Minimum Value in the Table% typical step of label setting methods (Dijsktra)
            Map.Entry<Integer, Double> min = Collections.min(TBL.entrySet(), Comparator.comparingDouble(Map.Entry::getValue));
            if(min != null){
                I = min.getKey();
                min_cost[I] = TBL.get(I);
                settled[I] = true;
            } else {
                isEmpty = true;
            }
        }

        //Store Costs and Paths
        return new Dijkstra2DResults(min_cost[(int) FID], path[(int) FID]);
    }


}
