package it.uniparthenope;

import it.uniparthenope.Boxing.Dijkstra2DResults;
import it.uniparthenope.Boxing.getNeighborsResults;

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
        getNeighborsResults adjacentNodes = getNeighbors(settledNodes, I_bool, I_ord, I_point, source, free_edges);
        for(int i=0;i<adjacentNodes.getRows().length;++i){
            int target = adjacentNodes.getNeighbors()[i];
            int row = adjacentNodes.getRows()[i];
            if(target != -1){
                if (getShortestDistance(target, distance) > (getShortestDistance(source, distance) + getDistance(row, edge_cost))){
                    distance.put(target, (getShortestDistance(source, distance) + getDistance(row, edge_cost)));
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

    private static getNeighborsResults getNeighbors(Set<Integer> settledNodes, boolean[] I_bool, int[] I_ord, int[] I_point, int I, int[][] free_edges){
        if(I>=I_bool.length)
            return new getNeighborsResults(new int[]{-1}, new int[]{-1});
        if(I_bool[I]){
            int head_ordinal = I_ord[I]-1;//Necessary because java indexes are from 0!
            int block_start = I_point[head_ordinal];
            int block_end = I_point[head_ordinal +1] -1;
            int[] neighbors = new int[(block_end-block_start)+1];
            int[] nids = new int[(block_end-block_start)+1];
            int tmp = block_start;
            for(int i=0;i<(block_end-block_start)+1; ++i){
                nids[i] = tmp;
                if(!settledNodes.contains(free_edges[tmp][1]))
                    neighbors[i] = free_edges[tmp][1];
                else
                    neighbors[i] = -1;
                ++tmp;
            }
            return new getNeighborsResults(nids, neighbors);
        } else
            return new getNeighborsResults(new int[]{-1}, new int[]{-1});
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
            findMinimalDistances(node,I_bool, I_ord, I_point, distance,free_edges, sh_delay, time_step, predecessors, unSettledNodes, settledNodes);
        }
        return new Dijkstra2DResults(distance.get((int) FID), getPath((int) FID, predecessors));
    }

    private static void findMinimalDistances(int source, boolean[] I_bool,int[] I_ord, int[] I_point, Map<Integer, Double> distance,
                                             int[][] free_edges, double[][] edge_cost, int time_step, Map<Integer, Integer> predecessors, Set<Integer> unSettledNodes, Set<Integer> settledNodes){

        getNeighborsResults adjacentNodes = getNeighbors(settledNodes, I_bool, I_ord, I_point, source, free_edges);
        for(int i=0;i<adjacentNodes.getRows().length;++i){
            int target = adjacentNodes.getNeighbors()[i];
            int row = adjacentNodes.getRows()[i];
            if(target != -1){
                if (getShortestDistance(target, distance) > (getShortestDistance(source, distance) + getDistance(row, edge_cost, time_step))){
                    distance.put(target, (getShortestDistance(source, distance) + getDistance(row, edge_cost, time_step)));
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
