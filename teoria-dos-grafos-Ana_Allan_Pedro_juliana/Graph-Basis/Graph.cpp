#include "Graph.hpp"
#include <limits>
#include <unordered_set>
#include <stack>
#include <set>

const float INF = std::numeric_limits<float>::infinity();

Graph::Graph() : _number_of_nodes(0), _number_of_edges(0), _directed(false), _weighted_edges(false), _weighted_nodes(false), _first(nullptr), _last(nullptr) {
}


Graph::~Graph()
{
}

void Graph::remove_node(size_t node_position)
{
    Node *node = this->_first;

    while (node != nullptr)
    {
        if (node->_id == node_position)
        {
            if (node == this->_first)
            {
                this->_first = node->_next_node;
                this->_first->_previous_node = nullptr;
            }
            else if (node == this->_last)
            {
                this->_last = node->_previous_node;
                this->_last->_next_node = nullptr;
            }
            else
            {
                Node *previous_node = node->_previous_node;
                Node *next_node = node->_next_node;
                previous_node->_next_node = next_node;
                next_node->_previous_node = previous_node;
            }

            //deletar as arestas ligadas nesse nó:
            Node *nodeCompared = this->_first;
            while (nodeCompared != nullptr)
            {
                if (conected(nodeCompared->_id, node_position) == 1)
                {
                    //_removed_edges[node_position].push_back({nodeCompared->_id, node_position});
                    remove_edge(nodeCompared->_id, node_position);
                }

                nodeCompared = nodeCompared->_next_node;
            }

            delete node;
            this->_number_of_nodes--;
            return;
        }

        node = node->_next_node;
    }
}

void Graph::remove_edge(size_t node_position_1, size_t node_position_2)
{
    Node *node1 = this->_first;
    Node *node2 = this->_first;

    while (node1 != nullptr)
    {
        if (node1->_id == node_position_1)
        {
            break;
        }
        node1 = node1->_next_node;
    }

    while (node2 != nullptr)
    {
        if (node2->_id == node_position_2)
        {
            break;
        }
        node2 = node2->_next_node;
    }

    Edge *edge = node1->_first_edge;
    while (&edge->_target_id != nullptr)
    {
        if (edge->_target_id == node2->_id)
        {
            break;
        }
        edge = edge->_next_edge;
    }

    if (edge == node1->_first_edge)
    {
        node1->_first_edge = edge->_next_edge;
    }
    else
    {
        Edge *aux = node1->_first_edge;
        while (aux->_next_edge != edge)
        {
            aux = aux->_next_edge;
        }
        aux->_next_edge = edge->_next_edge;
    }

    delete edge;
    node1->_number_of_edges = node1->_number_of_edges - 1;
    node2->_number_of_edges = node2->_number_of_edges - 1;
    
}

void Graph::add_node(size_t node_id, float weight)
{
    // para verificar se já não existe um nó com esse id
    for(Node* currentNode = _first; currentNode != nullptr; currentNode = currentNode->_next_node){
        if(currentNode->_id == node_id){
            return;
        } // um ponteiro criado dentro de uma estrutura de repetição é automaticamente deletado após o fim da mesma
    }
    // Se chegou aqui, é porque o nó não existe ainda, então vamos criá-lo!
    Node* newNode = new Node;
    newNode->_number_of_edges = 0;
    newNode->_id = node_id;
    newNode->_weight = weight;
    newNode->_first_edge = nullptr;
    newNode->_next_node = nullptr;
    newNode->_previous_node = _last;

    // agora, vamos ver se esse é o primeiro nó criado ou se ele deve ser posicionado após o último nó criado
    if(_first == nullptr){
        _first = newNode;
        _last = newNode;
    } 
    else{
        _last->_next_node = newNode;
        newNode->_previous_node = _last;
        _last = newNode;
    }
    _number_of_nodes++;
}



void Graph::add_edge(size_t source_node_id, size_t target_node_id, float weight) {
    
    Node* no_saida = nullptr;
    //Aqui ele vai procurar se o no de onde sai a aresta existe
    for (Node* current = _first; current != nullptr; current = current->_next_node) {
        if (current->_id == source_node_id) {
            no_saida = current; //Se ele encontrar
            break;
        }
    }

    // Se o nó de origem não existe, ele manda uma mensagem
    if (!no_saida) {
         std::cerr <<"nao existe um no com esse id "<< source_node_id << "para saida"<< std::endl;
            return;
    }

    Node* no_entrada = nullptr;
    //Aqui ele vai procurar se o no de entrada existe
    for (Node* current = _first; current != nullptr; current = current->_next_node) {
        if (current->_id == target_node_id) {
            no_entrada = current;//Se ele encontrar
            break;
        }
    }
    // Se o nó de origem não existe, ele manda uma mensagem
    if (!no_entrada) {
        std::cerr <<"nao existe um no com esse id"<< target_node_id << "para entrada" << std::endl;
            return;
    }

    // Cria uma nova aresta e declara tudo
    Edge* new_edge = new Edge{no_saida->_first_edge, weight, target_node_id};
    no_saida->_first_edge = new_edge;
    no_saida->_number_of_edges++;

}

void Graph::print_graph()
{
    Node *node = this->_first;
    while (node != nullptr)
    {
        std::cout << node->_id << " -- ";
        Edge *edge = node->_first_edge;
        while (edge != nullptr)
        {
            std::cout << edge->_target_id << " [" << edge->_weight << "] ";
            edge = edge->_next_edge;
        }
        std::cout << std::endl;
        node = node->_next_node;
    }
}


void Graph::print_graph(std::ofstream& output_file)
{
    Node *node = this->_first;
    while (node != nullptr)
    {
        output_file << node->_id << " -- ";
        Edge *edge = node->_first_edge;
        while (edge != nullptr)
        {
            output_file << edge->_target_id << edge->_weight;
            edge = edge->_next_edge;
        }
        output_file << std::endl;
        node = node->_next_node;
    }
    output_file.close();
}

int Graph::conected(size_t node_id_1, size_t node_id_2)
{
    Node *node1 = this->_first;
    Node *node2 = this->_first;

    while (node1 != nullptr)
    {
        if (node1->_id == node_id_1)
        {
            break;
        }
        node1 = node1->_next_node;
    }

    while (node2 != nullptr)
    {
        if (node2->_id == node_id_2)
        {
            break;
        }
        node2 = node2->_next_node;
    }

    Edge *edge = node1->_first_edge;
    while (edge != nullptr)
    {
        if (edge->_target_id == node2->_id)
        {
            return 1;
        }
        edge = edge->_next_edge;
    }
    
    return 0;
}
Node* Graph::get_node(size_t id) {
    Node* current = _first;
    while (current) {
        if (current->_id == id) return current;
        current = current->_next_node;
    }
    return nullptr;
}

int Graph::get_radius(){
        
    floyd_warshall();
    
    int radius = std::numeric_limits<int>::max();
    
    for (size_t i = 0; i < _number_of_nodes; ++i) {
       
        int eccentricity = 0;
        for (size_t j = 0; j < _number_of_nodes; ++j) {
            if (i != j) {
                
                if (dist[i][j] == std::numeric_limits<float>::infinity()) {
                    eccentricity = std::numeric_limits<int>::max();
                    break;
                }
                eccentricity = std::max(eccentricity, static_cast<int>(dist[i][j]));
            }
        }
        
    
        radius = std::min(radius, eccentricity);
    }
    
 
    return (radius == std::numeric_limits<int>::max()) ? -1 : radius;
}

int Graph::get_diameter(){
    floyd_warshall();
    
    int diameter = 0;
    
    for (size_t i = 0; i < _number_of_nodes; ++i) {
        for (size_t j = 0; j < _number_of_nodes; ++j) {
            if (i != j) {
                if (dist[i][j] == std::numeric_limits<float>::infinity()) {
                    return -1;
                }
                diameter = std::max(diameter, static_cast<int>(dist[i][j]));
            }
        }
    }
    
    return diameter;
}

std::vector<size_t> Graph::get_center() {
    floyd_warshall();
    
    std::vector<size_t> centers;
    int min_eccentricity = std::numeric_limits<int>::max();
    
    for (size_t i = 0; i < _number_of_nodes; ++i) {
        int eccentricity = 0;
        for (size_t j = 0; j < _number_of_nodes; ++j) {
            if (i != j) {
                if (dist[i][j] == std::numeric_limits<float>::infinity()) {
                    eccentricity = std::numeric_limits<int>::max();
                    break;
                }
                eccentricity = std::max(eccentricity, static_cast<int>(dist[i][j]));
            }
        }
        
        if (eccentricity < min_eccentricity) {
            min_eccentricity = eccentricity;
            centers.clear();
            centers.push_back(i);
        } else if (eccentricity == min_eccentricity) {
            centers.push_back(i);
        }
    }
    
    return centers;
}

std::vector<size_t> Graph::get_periphery() {
    floyd_warshall();
    
    std::vector<size_t> periphery;
    int max_eccentricity = 0;
    
    for (size_t i = 0; i < _number_of_nodes; ++i) {
        int eccentricity = 0;
        for (size_t j = 0; j < _number_of_nodes; ++j) {
            if (i != j) {
                if (dist[i][j] == std::numeric_limits<float>::infinity()) {
                    eccentricity = std::numeric_limits<int>::max();
                    break;
                }
                eccentricity = std::max(eccentricity, static_cast<int>(dist[i][j]));
            }
        }
        
        if (eccentricity > max_eccentricity) {
            max_eccentricity = eccentricity;
            periphery.clear();
            periphery.push_back(i);
        } else if (eccentricity == max_eccentricity) {
            periphery.push_back(i);
        }
    }
    
    return periphery;
}

std::vector<size_t> Graph::direct_closure(size_t node_id){

    Node* start_node = get_node(node_id); 
    if (!start_node) return {};  

    std::vector<size_t> closure;  

    // Percorre todas as arestas do nó inicial
    Edge* edge = start_node->_first_edge;
    while (edge != nullptr) {
      
        closure.push_back(edge->_target_id);
        edge = edge->_next_edge;
    }

    return closure;
}

std::vector<size_t> Graph::indirect_closure(size_t node_id){
    
    Node* firstNode = get_node(node_id);

    std::vector<size_t> closure;
    std::stack<Node*> stack;     

    stack.push(firstNode);

    while (!stack.empty()) {
        Node* current_node = stack.top();
        stack.pop();

        if (std::find(closure.begin(), closure.end(), current_node->_id) != closure.end()) {
            continue;
        }

        closure.push_back(current_node->_id);

        Edge* edge = current_node->_first_edge;
        while (edge != nullptr) {

            Node* target_node = get_node(edge->_target_id); 
            if (std::find(closure.begin(), closure.end(), target_node->_id) == closure.end()) {
                stack.push(target_node); 
            }
            edge = edge->_next_edge;
        }
    }

    return closure;  

}

void Graph::dist_min_Djkstra(size_t node_id_1, size_t node_id_2) {
    std::unordered_map<size_t, bool> visited;
    std::unordered_map<size_t, size_t> predecessors; // Mapeamento de predecessores
    auto compare = [](std::pair<float, size_t> a, std::pair<float, size_t> b) { return a.first > b.first; };
    std::priority_queue<std::pair<float, size_t>, std::vector<std::pair<float, size_t>>, decltype(compare)> pq(compare);

    Node* startNode = get_node(node_id_1);
    if (!startNode) {
        std::cerr << "esse no de partida nao existe" << std::endl;
        return;
    }

    pq.push({0.0, node_id_1});

    while (!pq.empty()) {
        size_t current_id = pq.top().second;
        pq.pop();
        if (visited[current_id]) continue;
        visited[current_id] = true;

        Node* current = get_node(current_id);
        if (current_id == node_id_2) break; // Para a busca ao alcançar o nó de destino

        for (Edge* edge = current->_first_edge; edge; edge = edge->_next_edge) {
            size_t neighbor_id = edge->_target_id;
            if (!visited[neighbor_id]) {
                pq.push({edge->_weight, neighbor_id});
                predecessors[neighbor_id] = current_id; // Armazena o predecessor
            }
        }
    }

    if (visited[node_id_2]) {
        // Reconstruindo o caminho
        std::vector<size_t> path;
        for (size_t at = node_id_2; at != node_id_1; at = predecessors[at]) {
            path.push_back(at);
        }
        path.push_back(node_id_1);

        std::reverse(path.begin(), path.end()); // Reverte o caminho para a ordem correta

        std::cout << "Caminho: ";
        for (size_t node : path) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    } else {
        std::cout << "nao ha caminho entre esses nos" << std::endl;
    }
}

void Graph::floyd_warshall() {
    std::unordered_map<size_t, size_t> node_to_index;
    std::unordered_map<size_t, size_t> index_to_node;
    size_t index = 0;

    //  garantindo que os ids dos nós sejam contíguos e comecem em 0
    for (Node* current = _first; current; current = current->_next_node) {
        node_to_index[current->_id] = index;
        index_to_node[index] = current->_id;
        index++;
    }

    // Inicializar as matrizes dist e next
    dist.assign(_number_of_nodes, std::vector<float>(_number_of_nodes, INF));// armazena a menor distância conhecida do nó i ao nó j
    next.assign(_number_of_nodes, std::vector<int>(_number_of_nodes, -1));// rastreia os próximos nós no caminho mínimo entre dois nós

    // Preenche as matrizes dist e next com base nos índices mapeados
    Node *current = _first;
    while (current) {
        size_t u = node_to_index[current->_id];
        dist[u][u] = 0;
        for (Edge *edge = current->_first_edge; edge; edge = edge->_next_edge) {
            size_t v = node_to_index[edge->_target_id];
            dist[u][v] = edge->_weight;
            next[u][v] = v;
        }
        current = current->_next_node;
    }

    // algoritmo Floyd-Warshall
    for (size_t k = 0; k < _number_of_nodes; ++k) {
        for (size_t i = 0; i < _number_of_nodes; ++i) {
            for (size_t j = 0; j < _number_of_nodes; ++j) {
                if (dist[i][k] < INF && dist[k][j] < INF) {
                    if (dist[i][j] > dist[i][k] + dist[k][j]) {
                        dist[i][j] = dist[i][k] + dist[k][j];
                        next[i][j] = next[i][k];
                    }
                }
            }
        }
    }


    for (size_t i = 0; i < _number_of_nodes; ++i) {
        for (size_t j = 0; j < _number_of_nodes; ++j) {     
        }
       
    }
}

//matriz de predecessores gerada pela floyd vai achar o caminho min entre dois nós
std::vector<size_t> Graph::get_shortest_path(size_t node_id_1, size_t node_id_2) {
    std::unordered_map<size_t, size_t> node_to_index;
    std::unordered_map<size_t, size_t> index_to_node;
    size_t index = 0;

    for (Node* current = _first; current; current = current->_next_node) {
        node_to_index[current->_id] = index;
        index_to_node[index] = current->_id;
        index++;
    }

    size_t u = node_to_index[node_id_1];
    size_t v = node_to_index[node_id_2];
    std::vector<size_t> path;

    if (next[u][v] == -1) return path;

    while (u != v) {
        path.push_back(index_to_node[u]);
        u = next[u][v];
    }
    path.push_back(index_to_node[v]);

    return path;
}

std::string Graph::prim(std::vector<size_t> subgraph) {
    
    Graph* graphLocal =  new Graph();
    int orderAux = 0;

    graphLocal->_weighted_edges = true;

    for(int i=0; i<subgraph.size(); i++){
        Node *nodeAux = get_node(subgraph[i]);
        Edge *edgeAux = nodeAux->_first_edge;
        orderAux+=1;
        while(edgeAux != nullptr){

            graphLocal->add_node(nodeAux->_id); // cria o no
            graphLocal->add_node(edgeAux->_target_id);
            graphLocal->add_edge(nodeAux->_id, edgeAux->_target_id, edgeAux->_weight);
            edgeAux = edgeAux->_next_edge;
            orderAux += 1;
        }
    }

    graphLocal->_number_of_nodes = orderAux;

    std::vector<bool> inAGM(graphLocal->_number_of_nodes, false);
    std::vector<int> key(graphLocal->_number_of_nodes, INT32_MAX);
    std::vector<int> dad(graphLocal->_number_of_nodes, -1);

    int source = graphLocal->_first->_id;
    key[source-1] = 0;

    std::priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::greater<std::pair<int,int>>> priorityqueue;
    priorityqueue.push(std::make_pair(0,source));

    while (!priorityqueue.empty()){
        int u = priorityqueue.top().second;
        priorityqueue.pop();

        if(inAGM[u-1]){
            continue;
        }

        inAGM[u-1] = true;

        Edge *edge = graphLocal->get_node(u)->_first_edge;

        while(edge != nullptr){
            int v = edge->_target_id;
            int weight = edge->_weight;

            if(!inAGM[v-1] && weight < key[v-1]){
                key[v-1] = weight;
                dad[v-1] = u;

                priorityqueue.push(std::make_pair(weight,v));
            }

            edge = edge->_next_edge;

        }

    }

    std::string result = "Arestas da Arvore Geradora Minima obtida pelo metodo de PRIM: \n";

    int totalWeight = 0;

    for(int i=0; i< graphLocal->_number_of_nodes; i++){
            
        if(dad[i] != -1){
            int u = dad[i];
            int v = i+1;
            int weight = key[i];

            result += std::to_string(u) + " -- " + std::to_string(v) + " [" + std::to_string(weight) +"]" "\n";

            totalWeight+= weight;
        }
    }

    result+= "Peso total da Arvore Geradora Minima: " + std::to_string(totalWeight) + "\n";

    return result;
}

std::string Graph::kruskal(std::vector<size_t>subgraph){

    Graph* localGraph =  new Graph();
    int orderAux = 0;

    localGraph->_weighted_edges = true;

    for(int i=0; i<subgraph.size(); i++){
        Node *nodeAux = this->get_node(subgraph[i]);
        Edge *edgeAux = nodeAux->_first_edge;
        orderAux+=1;
        while(edgeAux != nullptr){

            localGraph->add_node(nodeAux->_id); // cria o no
            localGraph->add_node(edgeAux->_target_id);
            localGraph->add_edge(nodeAux->_id, edgeAux->_target_id, edgeAux->_weight);
            edgeAux = edgeAux->_next_edge;
            orderAux += 1;
        }
    }

    localGraph->_number_of_nodes = orderAux;

    std::vector<std::pair<int, std::pair<int,int>>> orderedEdges;

    for(Node *currentNode = localGraph->_first; currentNode != nullptr; currentNode = currentNode->_next_node){

        Edge *edge =  currentNode->_first_edge;
        while(edge != nullptr){
            int sourceNode = currentNode->_id;
            int targetNode = edge->_target_id;
            int weight = edge->_weight;
            orderedEdges.push_back({weight, {sourceNode, targetNode}});
            edge = edge->_next_edge;
        }
    }

    std::sort(orderedEdges.begin(), orderedEdges.end());

    std::vector<int> group(localGraph->_number_of_nodes);

    for(int i = 0; i< localGraph->_number_of_nodes; i++){
        group[i] = i+1;
    }

    std::vector<std::pair<int,std::pair<int,int>>> msc;

    for(int i=0; i< orderedEdges.size(); i++){
        int origin = orderedEdges[i].second.first;
        int destiny = orderedEdges[i].second.second;

        std::cout << origin << std::endl;
        std::cout << destiny << std::endl;

        if(group[origin - 1] != group[destiny -1 ]){
            msc.push_back(orderedEdges[i]);

            int originSet = group[origin-1];
            int destinySet = group[destiny-1];
            for(int j=0; j<localGraph->_number_of_nodes; j++){
                if(group[j] == destinySet){
                    group[j] = originSet;
                }
            }

        }

    }

    std::string result = "Arestas da Arvore Geradora Minima obtida pelo metodo de Kruskal: \n";

    int totalWeight = 0;

    for(auto &edge : msc){
        int origin = edge.second.first;
        int destiny = edge.second.second;
        int weight = edge.first;
        std::cout << weight << std::endl;
        totalWeight += weight;
        result += std::to_string(origin) + " -- " + std::to_string(destiny) + " [" + std::to_string(weight)+ "]" "\n";
    }

    result +="Peso total da Arvore Geradora Minima obtida pelo metodo de Kruskal: " + std::to_string(totalWeight) + "\n";

    return result;
}

void Graph::depth_first_search(size_t start_node) {
    std::unordered_set<size_t> visited;
    std::unordered_map<size_t, size_t> parent;
    std::stack<size_t> stack;
    stack.push(start_node);
    parent[start_node] = -1; //nó inicial não tem pai

    while (!stack.empty()) {
        size_t current_node = stack.top();
        stack.pop();

        if (visited.find(current_node) != visited.end()) {
            continue;
        }

        visited.insert(current_node);
        std::cout << current_node << " ";

        Node* node = get_node(current_node);
        if (node) {
            for (Edge* edge = node->_first_edge; edge; edge = edge->_next_edge) {
                size_t target_id = edge->_target_id;

                if (visited.find(target_id) == visited.end()) {
                    stack.push(target_id);
                    parent[target_id] = current_node;
                } else if (target_id != parent[current_node]) {
                    std::cout << "\nAresta de retorno: " << current_node << " -> " << target_id << std::endl;
                }
            }
        }
    }
    std::cout << std::endl;
}

bool Graph::has_edges(size_t u) {
    Node* node = get_node(u);
    if (!node) {
        return false; // O vértice não existe
    }
    return node->_first_edge != nullptr; // Retorna true se houver pelo menos uma aresta conectada
}

void Graph::find_articulation_points() {
    std::unordered_set<size_t> articulation_points;
    int initial_components = count_connected_components();

    for (size_t u = 0; u < _number_of_nodes; ++u) {
        if (has_edges(u)) {
            remove_node(u);
            int new_components = count_connected_components();
            if (new_components < initial_components) {
                articulation_points.insert(u);
            }
            restore_node(u); // Restaura o vértice removido
        }
    }

    std::cout << "Pontos de articulação: ";
    for (size_t point : articulation_points) {
        std::cout << point << " ";
    }
    std::cout << std::endl;
}

int Graph::count_connected_components() {
    std::vector<bool> visited(_number_of_nodes, false);
    int component_count = 0;

    for (size_t i = 0; i < _number_of_nodes; ++i) {
        if (!visited[i] && has_edges(i)) {
            dfs_count_components(i, visited);
            ++component_count;
        }
    }

    return component_count;
}

void Graph::dfs_count_components(size_t u, std::vector<bool>& visited) {
    visited[u] = true;

    Node* node = get_node(u);
    if (!node) return;

    for (Edge* edge = node->_first_edge; edge != nullptr; edge = edge->_next_edge) {
        size_t v = edge->_target_id;
        if (!visited[v]) {
            dfs_count_components(v, visited);
        }
    }
}

void Graph::restore_node(size_t node_position) {
    // if (_removed_edges.find(node_position) == _removed_edges.end()) {
    // //     std::cerr << "Erro: Não há informações para restaurar o nó " << node_position << std::endl;
    // //     return;
    // // }

    // // Restaura o nó
    // Node* new_node = new Node();
    // new_node->_id = node_position;
    
    // // Insere o nó de volta na lista de nós do grafo
    // if (this->_first == nullptr) {
    //     this->_first = new_node;
    //     this->_last = new_node;
    // } else {
    //     this->_last->_next_node = new_node;
    //     new_node->_previous_node = this->_last;
    //     this->_last = new_node;
    // }

    // // Restaura as arestas removidas
    // for (const auto& edge : _removed_edges[node_position]) {
    //     // Arestas armazenadas como pares {source, target}
    //     add_edge(edge.first, edge.second);
    // }

    // // Limpa as arestas armazenadas após restaurar
    // _removed_edges.erase(node_position);

    // this->_number_of_nodes++;
}
