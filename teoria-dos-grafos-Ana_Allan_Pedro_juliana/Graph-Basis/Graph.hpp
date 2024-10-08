#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "Node.hpp"
#include "defines.hpp"
#include <string>
#include <unordered_set>
#include <unordered_map>

class Graph
{
public:
    /*Assinatura dos métodos básicos para o funcionamento da classe*/

    Graph(std::ifstream& instance);
    Graph();
    ~Graph();

    void remove_node(size_t node_id);
    void remove_edge(size_t node_id_1, size_t node_id_2);
    void add_node(size_t node_id, float weight = 0);
    void add_edge(size_t node_id_1, size_t node_id_2, float weight = 0);
    void print_graph(std::ofstream& output_file);
    void print_graph();
    void dist_min_Djkstra(size_t node_id_1, size_t node_id_2);
    Node* get_node(size_t id);
    int conected(size_t node_id_1, size_t node_id_2);
    void floyd_warshall();
    std::vector<size_t> get_shortest_path(size_t node_id_1, size_t node_id_2);
    int get_radius();
    int get_diameter();
    std::vector<size_t> get_center();
    std::vector<size_t> get_periphery();
    std::vector<size_t> direct_closure(size_t node_id);
    std::vector<size_t> indirect_closure(size_t node_id);
    std::string prim(std::vector<size_t> subgraph);
    std::string kruskal(std::vector<size_t>subgraph);
    void depth_first_search(size_t start_node);
    void find_articulation_points();
    bool has_edges(size_t u);
    int count_connected_components();
    void dfs_count_components(size_t u, std::vector<bool>& visited);
    void restore_node(size_t node_position);
    
private:
    size_t _number_of_nodes;
    size_t _number_of_edges;
    bool   _directed;
    bool   _weighted_edges;
    bool   _weighted_nodes;
    Node  *_first;
    Node  *_last;
    std::vector<std::vector<float>> dist;
    std::vector<std::vector<int>> next;
    //std::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>> _removed_edges;
};

#endif  //GRAPH_HPP
