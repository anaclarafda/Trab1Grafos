#include "Graph.hpp"



int main(int argc, char* argv[])
{
    if (argc != 5) { 
        std::cerr << "Uso: " << argv[0] << "<arquivo_entrada> <arquivo_saida> <Op_Direc> <Op_PesoAresta>\n";
        return 1;
    }

    // Lê os argumentos da linha de comando
    const char* entrada = argv[1];
    const char* saida = argv[2];
    int direcionado = std::stoi(argv[3]);
    int ponderado = std::stoi(argv[4]);

    int valor = -1;
    std::ifstream arquivoEntrada(entrada);
    std::ifstream arquivoSaida(saida);
    
    
    if(!arquivoEntrada){
         std::cerr << "Erro ao abrir o arquivo de instancia." << std::endl;
         return 1;
    }

     Graph *graph = new Graph();

    size_t number_of_nodes;
    arquivoEntrada >> number_of_nodes; // Lê o número de nós do arquivo de instância

    if(ponderado == 0 ){
        for (size_t i = 0; i < number_of_nodes; ++i) {
            size_t node_id1;
            size_t node_id2;  
            arquivoEntrada >> node_id1;// Lê o no de saida
            graph->add_node(node_id1); // cria o no
            arquivoEntrada >> node_id2; // Lê o no alvo
            graph->add_node(node_id2); // cria o no
            float weight;
            arquivoEntrada >> weight; // Lê o edge weight
            graph->add_edge(node_id1,node_id2,0);
     } 
    }
    else{
        for (size_t i = 0; i < number_of_nodes; ++i) {
            size_t node_id1;
            size_t node_id2;  
            arquivoEntrada >> node_id1;// Lê o no de saida
            graph->add_node(node_id1); // cria o no
            arquivoEntrada >> node_id2; // Lê o no alvo
            graph->add_node(node_id2); // cria o no
            float weight;
            arquivoEntrada >> weight; // Lê o edge weight
            graph->add_edge(node_id1,node_id2,weight);
     } 
    }
    std::cerr << "O seu grafo foi criado! Veja-o abaixo:" << std::endl;
    graph->print_graph();

    valor = 1;



    while (valor != 0 )
    {
        std::cerr << "FUNCIONALIDADES: O que voce quer ver?" << std::endl;
        std::cerr << "Selecione um dos valores abaixo:" << std::endl;
        std::cerr << "(1) Caminho min entre dois vertices (alg Djkstra)       (2) Caminho min entre dois vertices (alg Floyd)" << std::endl;
        std::cerr << "(3) Arvore em ordem de caminhamento em profundidade     (4) O raio, o diametro, o centro e a periferia do grafo" << std::endl;
        std::cerr << "(5) O conjunto de vertices de articulacao               " << std::endl;

        if(ponderado == 1) //Só aparecem para grafos ponderados
        {
            std::cerr << "(6) Arvore Geradora Minima (alg Prim)                   (7) Arvore Geradora Minima (alg Kruskal)" << std::endl;
        }
        if(direcionado == 1) // Essas duas funcio. so aparecem p/ grafos direcionados
        {
            std::cerr << "(8) Fecho transitivo direto de um vertice "<< std::endl << "(9) Fecho transitivo indireto de um vertice " << std::endl;
        }

        std::cerr << "Clique no 0 para finalizar" << std::endl;
        std::cin >> valor;
        
        if(valor == 1) //Caminho min entre dois vertices (alg Djkstra)
        {
            int vertice1, vertice2; // a funcao usa dois parametros (vertice inicial e vertice final)

            std::cerr << "Qual o vertice inicial?" << std::endl;
            std::cin >> vertice1;
            std::cerr << "Qual o vertice final?" << std::endl;
            std::cin >> vertice2; 
            graph->dist_min_Djkstra(vertice1, vertice2);//  chama a funcao 
        }
        else if(valor == 2) //Caminho min entre dois vertices (alg Floyd)
        {
            int vertice1, vertice2; // a funcao usa dois parametros (vertice inicial e vertice final)
            std::cerr << "Qual o vertice inicial?" << std::endl;
            std::cin >> vertice1;
            std::cerr << "Qual o vertice final?" << std::endl;
            std::cin >> vertice2;

            graph->floyd_warshall();//esta dando algum erro dentro da funcao


            std::vector<size_t> path = graph->get_shortest_path(vertice1, vertice2);
            std::cout << "Caminho minimo de " << vertice1 << " para " << vertice2 << ":" << std::endl;
            for (size_t node : path) {
                std::cout << node << " ";
            }
            std::cout << std::endl;

            delete graph;
           
        }
        else if(valor == 3) //Arvore em ordem de caminhamento em profundidade
        {
            int vertice_inicial;
            std::cerr << "Qual o vértice inicial para o caminhamento em profundidade?" << std::endl;
            std::cin >> vertice_inicial;

            graph->depth_first_search(vertice_inicial);
        }
        else if(valor == 4) //O raio, o diametro, o centro e a periferia do grafo
        {
            int radius = graph->get_radius();
            int diameter = graph->get_diameter();
            std::vector<size_t> centerNodes = graph->get_center();
            std::vector<size_t> peripheryNodes = graph->get_periphery();

            std::cerr << "O raio do grafo é  " << radius << std::endl << "O diametro do grafo é " << diameter << std::endl;


            std::cerr << "O centro do grafo é compostos pelo nos: ";  
            for (const size_t& value : centerNodes) {
                std::cerr << value << " ";
            }

            std::cerr << std::endl;
             
            std::cerr << "A periferia do grafo é composta pelos nos: ";
            for (const size_t& value : peripheryNodes) {
                std::cout << value << " ";
            }

            std::cerr << std::endl;
        }
        else if(valor == 5) //Conjunto de vertices de articulação
        {
            graph->find_articulation_points();
        }
        else if(valor == 6) //Arvore Geradora Minima (alg Prim) 
        {
            bool readSubGraph = true;
            std::vector<size_t> vectorAux = {};
            size_t no;

            while(readSubGraph){
                
                std::cout << "Seu sub-grafo é: ";
                for(int i=0; i< vectorAux.size(); i++){
                    std::cout << vectorAux[i] << " ";
                }

                std::cout << std::endl;

                std::cout << "Insira um novo no no sub-grafo ou digite -1 para envia-lo a funcao:" << std::endl;

                std::cin>> no;

                if(no == -1){
                    readSubGraph = false;
                }else{
                    bool flag = false;
                    for(int i;i< vectorAux.size(); i++){
                        if(no == vectorAux[i]){
                            flag = true;
                        }     
                    }
                    if(flag == false){
                        if(graph->get_node(no) == nullptr){
                            std::cout << "Nó nao existente no grafo, favor inserir um nó valido."<< std::endl;
                        }else{
                            vectorAux.push_back(no);
                        }
                    }else{
                        std::cout<< "Vértice repetido, favor inserir um nó válido." << std::endl;
                    }
                }
            }

            if(vectorAux.empty()){
                std::cout<< "Sub-grafo vazio." << std::endl;
            }else{
                std::string result = graph->prim(vectorAux);

                std::cout << result << std::endl;
            }
            
        }
        else if(valor == 7) //Arvore Geradora Minima (alg Kruskal)
        {   
            bool readSubGraph = true;
            std::vector<size_t> vectorAux = {};
            size_t no;

            while(readSubGraph){
                
                std::cout << "Seu sub-grafo e: ";
                for(int i=0; i< vectorAux.size(); i++){
                    std::cout << vectorAux[i] << " ";
                }

                std::cout << std::endl;

                std::cout << "Insira um novo no no sub-grafo ou digite -1 para envia-lo a funcao:" << std::endl;

                std::cin>> no;

                if(no == -1){
                    readSubGraph = false;
                }else{
                    bool flag = false;
                    for(int i;i< vectorAux.size(); i++){
                        if(no == vectorAux[i]){
                            flag = true;
                        }     
                    }
                    if(flag == false){
                        if(graph->get_node(no) == nullptr){
                            std::cout << "No nao existente no grafo, favor inserir um no valido."<< std::endl;
                        }else{
                            vectorAux.push_back(no);
                        }
                    }else{
                        std::cout<< "Vertice repetido, favor inserir um no valido." << std::endl;
                    }
                }
            }

            if(vectorAux.empty()){
                std::cout<< "Sub-grafo vazio." << std::endl;
            }else{
                std::string result = graph->kruskal(vectorAux); // -> IMPLEMENTAR AINDA!

                std::cout << result << std::endl;
            }

        }
        else if(valor == 8) //Fecho transitivo direto de um vertice
        {
            int vertice;
            std::cerr << "Qual o vertice desejado?" << std::endl;
            std::cin >> vertice;

            std::vector<size_t> closure = graph->direct_closure(vertice);

            std::cerr << "O fecho é composta pelos nos: ";
            for (const size_t& value : closure) {
                std::cout << value << " ";
            }

            std::cerr << std::endl;

        }
        else if(valor == 9) //Fecho transitivo indireto de um vertice
        {
            int vertice;
            std::cerr << "Qual o vertice desejado?" << std::endl;
            std::cin >> vertice;

            std::vector<size_t> closure = graph->indirect_closure(vertice);

            std::cerr << "O fecho é composta pelos nos: ";
            for (const size_t& value : closure) {
                std::cout << value << " ";
            }

            std::cerr << std::endl;
        }
    }
    
    
    return 0;
}