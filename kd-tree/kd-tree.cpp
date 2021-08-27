#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <cmath>
#include <map>

using namespace std;

struct Node{
    vector<double> point;
    int dim;
    int k;
    Node* left;
    Node* right; 
    map<string, double> info;
    
    Node(int new_k, int new_dim, vector<double> new_point, map<string, double> new_info){
        k = new_k;
		dim = new_dim;
		point = new_point;
		info = new_info;
        left = nullptr;
        right = nullptr;
    }
    
};


struct fake_bpq{
	vector<Node*> nodes;
	vector<double> dists;	
	double biggest_dist;
	int max_nodes;
	
	fake_bpq(){
		biggest_dist = -1;
	}
	
	reset(int max_points){
		max_nodes = max_points;
		nodes.clear();
		dists.clear();
		biggest_dist = -1;
	}
	
	bool is_not_full(){
		return nodes.size() < max_nodes;
	}
	
	/* 
	 * Function will check if the distance is smaller than the distances,
	 * if it found one, it will replace the distance and the node on the vectors.
	 */
	void add_node(Node* pNode, double new_dist){
		if(is_not_full()){
			nodes.push_back(pNode);
			dists.push_back(new_dist);
			return;
		}
		
		vector<double>::iterator it = dists.begin();
		vector<Node*>::iterator it_n = nodes.begin();
		
		while(new_dist > *it){
			if(*it > biggest_dist) biggest_dist = *it;
			it++; it_n++;
			if(it == dists.end()){
				return;
			}
		}
		*it = new_dist;
		*it_n = pNode;	
		
		while(it != dists.end()){
			if(*it > biggest_dist) biggest_dist = *it;
			it++;
		}
	}
};


class KDTree{
	public: 
    Node* pRoot;
    int k;
    vector<double> coord_weights;
    fake_bpq candidates;
    
    KDTree(int new_k){
    	k = new_k;
        pRoot = nullptr;
    }
    
    void print(){
    	print(pRoot);
	}
	
	void print(Node * pNode){
		if(pNode == nullptr) return;
		cout << "Node point:";
		vector<double>::iterator it = pNode->point.begin();
		for(; it != pNode->point.end(); it++){
			cout << *it << " ";
		}
		cout << endl;
		print(pNode->left);
		print(pNode->right);
	}

    void insert(vector<double> new_point, map<string, double> info){
        Node** pNode = &pRoot;
        int dim = 0;
        while((*pNode) != nullptr){
            vector<double>::iterator it = (*pNode)->point.begin();
            vector<double>::iterator it_new = new_point.begin();
            for(int i = 0; i < dim; i++){
            	it++; it_new++;
			}
			
            if((*it_new) < (*it)){
                pNode = &((*pNode)->left);
            }else{
                pNode = &((*pNode)->right);
            }
            dim = (dim + 1)%3;
        }
        
        (*pNode) = new Node(k, dim, new_point, info);
    }
	
    void k_neighbors(int max_points, vector<double> point, string filename){
		Node *pNode = pRoot;
		int dim = 0;
		
		//Setting the bpq
		candidates.reset(max_points);
	
		//Start recursive search	
		recursive_search(point, pNode, dim);
		
		//Saving the result in a file
		ofstream output (filename);
		vector<Node*>::iterator it = candidates.nodes.begin();
		for(; it != candidates.nodes.end(); it++){
			vector<double>::iterator it2 = (*it)->point.begin();
			for(; it2 != (*it)->point.end(); it2++){
				output << (*it2) << " ";
			}
			
			for (auto& t : (*it)->info){
				output << t.second << " ";
			}
			output << endl;
		}
		output.close();
	}
	
	void recursive_search(vector<double> point, Node* pNode, int dim){
		if(pNode == nullptr) return;
		
		//Compute distance between searched point and current point
		double dist = 0;
		vector<double>::iterator it = pNode->point.begin();
        vector<double>::iterator it_point = point.begin();
        vector<double>::iterator it_weight = coord_weights.begin();
        
        for(int i = 0; i < k; i++){
        	dist = dist + (*it - *it_point) * (*it - *it_point) * (*it_weight);
        	it++; it_point++; it_weight++;
		}
		
		//Update bounded priority queue
		candidates.add_node(pNode, dist);
		
		//Check what side to search recursively
		it = pNode->point.begin();
        it_point = point.begin();
        for(int i = 0; i < dim; i++){
        	it++; it_point++;
		}
		bool search_left = false;
		if(*it < *it_point){
			search_left = true;
			recursive_search(point, pNode->left, dim+1);
		}else{
			recursive_search(point, pNode->right, dim+1);
		}
		
		
		//Check if should verify also the other side
		if(candidates.is_not_full() | 
		  abs(*it - *it_point) < candidates.biggest_dist){
			if(search_left){
		  		recursive_search(point, pNode->right, dim+1);
			}else{
				recursive_search(point, pNode->left, dim+1);
			}
		}
	}
    
};


void create_random_farms(string filename, int n){
	ofstream output (filename);
	for(int i = 0; i < n; i++){
		double strawberry_prod = 500 + rand()%400;
		double property_size = 40 + rand()%60;
		double sell_value = (100 + rand()%400);
		sell_value = sell_value/100;
		double x = rand()%10000;
		double y = rand()%10000;
		output << strawberry_prod << " " << property_size << " " << sell_value << " " << x << " " << y << endl;
	}
	output.close();
}

void create_first_tree(KDTree * tree, string filename){
	ifstream input(filename);
	string point_line;
	while(getline(input, point_line)){
		//Read from line and transform in doubles
        stringstream line_stream(point_line);
        string strawberry_prod_s, property_size_s, sell_value_s, x_s, y_s;
        double strawberry_prod, property_size, sell_value, x, y;
        line_stream >> strawberry_prod_s >> property_size_s >> sell_value_s >> x_s >> y_s;
        strawberry_prod = stod(strawberry_prod_s);
        property_size = stod(property_size_s);
        sell_value = stod(sell_value_s);
        x = stod(x_s);
        y = stod(y_s);
        
        //Add new node to tree
        vector<double> farm = {strawberry_prod, property_size, sell_value};
        map<string, double> info = {{"x", x}, {"y", y}};
        tree->insert(farm, info);
	}
	tree->coord_weights.push_back(0.01);
	tree->coord_weights.push_back(0.1);
	tree->coord_weights.push_back(1);
	input.close();	
}

void create_second_tree(KDTree * tree, string filename){
	ifstream input(filename);
	string point_line;
	while(getline(input, point_line)){
		//Read from line and transform in doubles
        stringstream line_stream(point_line);
        string strawberry_prod_s, property_size_s, sell_value_s, x_s, y_s;
        double strawberry_prod, property_size, sell_value, x, y;
        line_stream >> strawberry_prod_s >> property_size_s >> sell_value_s >> x_s >> y_s;
        strawberry_prod = stod(strawberry_prod_s);
        property_size = stod(property_size_s);
        sell_value = stod(sell_value_s);
        x = stod(x_s);
        y = stod(y_s);
        
        //Add new node to tree
        vector<double> farm = {x, y};
        map<string, double> info = {{"strawberry_prod", strawberry_prod}, 
									{"property_size", property_size},
									{"sell_value", sell_value}};
        tree->insert(farm, info);
	}
	tree->coord_weights.push_back(1);
	tree->coord_weights.push_back(1);
	input.close();	
	
}

int main(){
	create_random_farms("random_farms.txt", 20);
    KDTree tree_3d(3);
    //vector<double> farm1 = {600, 60, 2000000};
    //tree.insert(farm1);
    create_first_tree(&tree_3d,"random_farms.txt");
    tree_3d.print();
    vector<double> farm1 = {600, 60, 2};
    tree_3d.k_neighbors(6, farm1, "result_closest.txt");
    
    KDTree tree_2d(2);
    create_second_tree(&tree_2d, "result_closest.txt");
	tree_2d.print();
	vector<double> farm1_coord = {650, 200};
	tree_2d.k_neighbors(2, farm1_coord, "result_closest_2.txt");
	
    return 0;
}

