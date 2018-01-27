package it.uniparthenope;

import it.uniparthenope.Boxing.Dijkstra2DResults;

import java.util.*;

public class Algorithms {

    //For GDT route
    public static Dijkstra2DResults Dijkstra(int[][] free_edges, double[] edge_costs, boolean[] I_bool,
                                             int[] I_ord, int[] I_point, long SID, long FID){
        /*
         * nodes = lon and lat for each node in graph
         * free_edges = free edges
         * edge_cost = weights
         * I_bool, I_ord, I_point = for "forward star" data structure
         * SID = Starting node
         * FID = Ending node
         * */
        Set<Integer> unSettledNodes = new HashSet<>();
        Set<Integer> settledNodes = new HashSet<>();
        Map<Integer, Integer> predecessors = new HashMap<>();
        Map<Integer, Double> distance = new HashMap<>();

        distance.put((int) SID, 0.0);
        unSettledNodes.add((int) SID);

        while(unSettledNodes.size() > 0){
            int node = getMinimum(unSettledNodes, distance);
            settledNodes.add(node);
            unSettledNodes.remove(node);
            findMinimalDistances(node,I_bool, I_ord, I_point, distance,free_edges, edge_costs, predecessors, unSettledNodes, settledNodes);
        }
        return new Dijkstra2DResults(distance.get((int) FID), getPath((int) FID, predecessors));
    }

    public static LinkedList<Integer> getPath(int target, Map<Integer, Integer> predecessors){
        LinkedList<Integer> path = new LinkedList<>();
        int step = target;
        if(predecessors.get(step) == null)
            return null;
        path.add(step);
        while(predecessors.get(step) != null){
            int next = predecessors.get(step);
            step = next;
            path.add(step);
        }
        Collections.reverse(path);
        return path;
    }

    private static void findMinimalDistances(int source, boolean[] I_bool,int[] I_ord, int[] I_point, Map<Integer, Double> distance,int[][] free_edges, double[] edge_cost, Map<Integer, Integer> predecessors, Set<Integer> unSettledNodes,Set<Integer> settledNodes){
        int[] adjacentNodes = getNeighbors(settledNodes, I_bool, I_ord, I_point, source, free_edges);
        for(int target: adjacentNodes){
            if(target != -1) {
                if (getShortestDistance(target, distance) > (getShortestDistance(source, distance) + getDistance(target, edge_cost))) {
                    distance.put(target, getShortestDistance(source, distance) + getDistance(target, edge_cost));
                    predecessors.put(target, source);
                    unSettledNodes.add(target);
                }
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

    private static int[] getNeighbors(Set<Integer> settledNodes, boolean[] I_bool, int[] I_ord, int[] I_point, int I, int[][] free_edges){
        if(I>=I_bool.length)
            return new int[]{-1};
        if(I_bool[I]){
            int head_ordinal = I_ord[I]-1;//Necessary because java indexes are from 0!
            int block_start = I_point[head_ordinal];
            int block_end = I_point[head_ordinal +1] -1;
//            int[] destinations = new int[(block_end-block_start)+1];
//            int tmp = block_start;
//            for(int i=0;i<destinations.length; ++i) {
//                destinations[i] = free_edges[tmp][1];
//                ++tmp;
//            }
//            return destinations;
            ArrayList<Integer> tmp_dest = new ArrayList<>();
            int tmp = block_start;
            for(int i = 0; i<(block_end-block_start)+1; ++i ){
                if(!settledNodes.contains(free_edges[tmp][1])){
                    tmp_dest.add(free_edges[tmp][1]);
                }
                ++tmp;
            }
            int[] neighbors = new int[tmp_dest.size()];
            for(int i=0;i<tmp_dest.size();++i){
                neighbors[i] = tmp_dest.get(i);
            }
            return neighbors;
        } else
            return new int[]{-1};
    }

    //FOR STATIC ALGORITHM
    public static Dijkstra2DResults Dijkstra(int[][] free_edges, double[][] sh_delay, int time_step, boolean[] I_bool,
                                             int[] I_ord, int[] I_point, long SID, long FID){
        Set<Integer> unSettledNodes = new HashSet<>();
        Map<Integer, Integer> predecessors = new HashMap<>();
        Map<Integer, Double> distance = new HashMap<>();
        Set<Integer> settledNodes = new HashSet<>();

        distance.put((int) SID, 0.0);
        unSettledNodes.add((int) SID);

        while(unSettledNodes.size() > 0){
            int node = getMinimum(unSettledNodes, distance);
            settledNodes.add(node);
            unSettledNodes.remove(node);
            findMinimalDistances(node, (int) FID,I_bool, I_ord, I_point, distance,free_edges, sh_delay, time_step, predecessors, unSettledNodes, settledNodes);
        }
        return new Dijkstra2DResults(distance.get((int) FID), getPath((int) FID, predecessors));
    }

    private static void findMinimalDistances(int source, int destination, boolean[] I_bool,int[] I_ord, int[] I_point, Map<Integer, Double> distance,
                                             int[][] free_edges, double[][] edge_cost, int time_step, Map<Integer, Integer> predecessors, Set<Integer> unSettledNodes, Set<Integer> settledNodes){
        int[] adjacentNodes = getNeighbors(settledNodes, I_bool, I_ord, I_point, source, free_edges);
        for(int target : adjacentNodes){
            if(target != -1){
                if(getShortestDistance(target, distance) > (getShortestDistance(source, distance) + getDistance(target, edge_cost, time_step))){
                    distance.put(target, (getShortestDistance(source, distance) + getDistance(target, edge_cost, time_step)));
                    predecessors.put(target, source);
                    unSettledNodes.add(target);
                }
            }
        }
    }

    private static double getDistance(int destination, double[][] edge_cost, int time_step){
        return edge_cost[destination][time_step];
    }
}
