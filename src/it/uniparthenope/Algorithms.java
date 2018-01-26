package it.uniparthenope;

import com.google.common.primitives.Ints;
import com.mchange.v1.util.ArrayUtils;
import it.uniparthenope.Boxing.Dijkstra2DResults;

import java.util.*;

public class Algorithms {

    public static Dijkstra2DResults Dijkstra(int Nnodes, int[][] free_edges, double[] edge_costs, boolean[] I_bool,
                                             int[] I_ord, int[] I_point, long SID, long FID){
        /*
         * nodes = lon and lat for each node in graph
         * free_edges = free edges
         * edge_cost = weights
         * I_bool, I_ord, I_point = for "forward star" data structure
         * SID = Starting node
         * FID = Ending node
         * */

        //IDXs
        Set<Integer> unSettledNodes = new HashSet<>();
        Map<Integer, Integer> predecessors = new HashMap<>();
        Map<Integer, Double> distance = new HashMap<>();

        distance.put((int) SID, 0.0);
        unSettledNodes.add((int) SID);

        while(unSettledNodes.size() > 0){
            int node = getMinimum(unSettledNodes, distance);
            unSettledNodes.remove(node);
            findMinimalDistances(node, (int) FID,I_bool, I_ord, I_point, distance,free_edges, edge_costs, predecessors, unSettledNodes);
        }
        return new Dijkstra2DResults(distance.get((int) FID), getPath((int) FID, predecessors, free_edges, edge_costs));
        //return getPath((int) FID, predecessors, free_edges, edge_costs);
    }

    public static LinkedList<Integer> getPath(int target, Map<Integer, Integer> predecessors, int[][] free_edges, double[] edge_costs){
        LinkedList<Integer> path = new LinkedList<>();
        double cost = 0.0;
        int step = target;
        if(predecessors.get(step) == null)
            return null;
        path.add(step);
        while(predecessors.get(step) != null){
            int next = predecessors.get(step);
//            cost += edge_costs[free_edges[next][step]];
            //step = predecessors.get(step);
            step = next;
            path.add(step);
        }
        Collections.reverse(path);
//        System.out.println("COSTO: "+cost);
        return path;
    }

    private static void findMinimalDistances(int source, int destination, boolean[] I_bool,int[] I_ord, int[] I_point, Map<Integer, Double> distance,int[][] free_edges, double[] edge_cost, Map<Integer, Integer> predecessors, Set<Integer> unSettledNodes){
        int[] adjacentNodes = getNeighbors(I_bool, I_ord, I_point, source, free_edges);
        for(int target: adjacentNodes){
            if(target != -1) {
                if (getShortestDistance(target, distance) > (getShortestDistance(source, distance) + getDistance(target, edge_cost))) {
                    distance.put(target, getShortestDistance(source, distance) + getDistance(target, edge_cost));
                    predecessors.put(target, source);
                    unSettledNodes.add(target);
                }
                if(target == destination)
                    unSettledNodes.removeAll(unSettledNodes);
            }
        }
    }

    private static double getDistance(int destination, double[] edge_cost){
        return edge_cost[destination];
    }

    private static int getMinimum(Set<Integer> src, Map<Integer, Double> distance){
        Integer minimum = null;
        for(int vertex : src){
            if(minimum == null)
                minimum = vertex;
            else{
                if( getShortestDistance(vertex, distance) < getShortestDistance(minimum, distance))
                    minimum = vertex;
            }
        }
        return minimum;
    }

    private static double getShortestDistance(int vertex, Map<Integer, Double> distance){
        Double dist = distance.get(vertex);
        if(dist == null)
            return Double.POSITIVE_INFINITY;
        else
            return dist;
    }

    private static int[] getNeighbors(boolean[] I_bool, int[] I_ord, int[] I_point, int I, int[][] free_edges){
        if(I>=I_bool.length)
            return new int[]{-1};
        if(I_bool[I]){
            int head_ordinal = I_ord[I]-1;//Necessary because java indexes are from 0!
            int block_start = I_point[head_ordinal];
            int block_end = I_point[head_ordinal +1] -1;
            int[] destinations = new int[(block_end-block_start)+1];
            int tmp = block_start;
            for(int i=0;i<destinations.length; ++i) {
                destinations[i] = free_edges[tmp][1];
                ++tmp;
            }
            return destinations;
        } else
            return new int[]{-1};
    }

}
