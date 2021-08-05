#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <algorithm>
#include <random>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point; 
typedef Kernel::Triangle_2 Triangle;

struct Halfedge; struct TriangleNode;

struct Vertex{
  Point p;
  Halfedge* he;
  int index, r, g, b;

  Vertex(){}

  Vertex(Point new_p){
    p = new_p;
  }

  /**
   * @brief Create a copy of the vertex and return.
   * 
   */
  Vertex copy(){
    Vertex * v = new Vertex();
    v->p = p;
    v->he = he; 
    v->index = index;
    return *v;
  }

  /**
   * @brief Print the information of the vertex.
   * 
   */
  void print(){
    std::cout << " p(" << index << ")";
    return;
  }

};

struct Halfedge{
  Halfedge* twin;
  Halfedge* next;
  Vertex* origin;
  TriangleNode* tri;

  Halfedge(){}
};

struct TriangleNode{
  std::vector<TriangleNode*> childs;
  Triangle t;
  Vertex nodeVertex[3];
  Halfedge* he;
  int depth;
  
  TriangleNode(){ }
  
  /** 
   * @brief Create a new TriangleNode with the link to the halfedge, and set the nodeVertex and t attributes.
   * 
   */
  TriangleNode(Halfedge* new_he){
    he = new_he;
    nodeVertex[0] = new_he->origin->copy();
    nodeVertex[1] = new_he->next->origin->copy();
    nodeVertex[2] = new_he->next->next->origin->copy();
    Triangle new_triangle (nodeVertex[0].p, nodeVertex[1].p, nodeVertex[2].p);
    t = new_triangle;
  }

  /**
   * @brief Print the vertex index of this TriangleNode.
   * 
   */
  void print(){
    std::cout << "T: "; 
    for(int i = 0; i < 3; i++){
      nodeVertex[i].print();
    }
    std::cout << " d: " << depth << std::endl;
  }


  /**
   * @brief Verify if a point is inside the triangle of the TriangleNode.
   * 
   * @param p 
   * @return true 
   * @return false 
   */
  bool contains_point(Point p){
    return t.bounded_side(p) == CGAL::ON_BOUNDED_SIDE | t.bounded_side(p) == CGAL::ON_BOUNDARY;
  }
};

class DCEL{
  public:
  std::vector<Vertex*> vertices;
  std::vector<Halfedge*> halfedges;

  /**
   * @brief Print the vertices of each edge (the origin of the halfedge and the origin of the twin).
   * 
   */
  void print_edges(){
    std::cout << std::endl << "DCEL halfedges: " << std::endl;
    std::vector<Halfedge*>::iterator it = halfedges.begin();
    for(; it != halfedges.end(); it++){
      (*it)->origin->print();
      std::cout << "------";
      (*it)->twin->origin->print();
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  /**
   * @brief Print all halfedges that don't contain helper vertex.
   * 
   */
  void print_final_edges(){
    std::cout << std::endl << "DCEL halfedges: " << std::endl;
    std::vector<Halfedge*>::iterator it = halfedges.begin();
    for(; it != halfedges.end(); it++){
      if(((*it)->origin->index >= 0) & ((*it)->next->origin->index >= 0)){
        (*it)->origin->print();
        std::cout << "------";
        (*it)->twin->origin->print();
        std::cout << std::endl;
      }
    }
    std::cout << std::endl;
  }

  /**
   * @brief Save the triangulation to a file.
   * 
   */
  void save_triangulation(){
    std::ofstream output ("output.json");

    output << "[" << std::endl;
    std::vector<Halfedge*>::iterator it = halfedges.begin();
    for(; it != halfedges.end(); it++){
      Vertex *v1, *v2, *v3;
      v1 = (*it)->origin;
      v2 = (*it)->next->origin;
      v3 = (*it)->next->next->origin;
      if((v1->index >= 0) && (v2->index >= 0) & (v3->index>= 0)){
        output << "{\"points\": [" ;
        output << "[" << v1->p.x() << "," << v1->p.y() << "],";
        output << "[" << v2->p.x() << "," << v2->p.y() << "],";
        output << "[" << v3->p.x() << "," << v3->p.y() << "]";
        output << "],\"colors\": [";
        output << "[" << v1->r << "," << v1->g << "," << v1->b << "],";
        output << "[" << v2->r << "," << v2->g << "," << v2->b << "],";
        output << "[" << v3->r << "," << v3->g << "," << v3->b << "]";
        output << "]}," << std::endl;
      }
    }
    output << "]" << std::endl;
    output.close();
    return;
  }
  

  /**
   * @brief Print the 3 vertex of the triangle that contain the halfedge.
   * 
   * @param he halfedge.
   */
  void print_triangle(Halfedge *he){
    std::cout << "Triangle: ";
    he->origin->print();
    he->next->origin->print();
    he->next->next->origin->print();
    std::cout << std::endl;
  }

  /**
   * @brief Add vertex to DCEL.
   * 
   * @param v 
   */
  void add_vertex(Vertex* v){
    v->index = vertices.back()->index + 1;
    vertices.push_back(v);
    return;
  }

  void set_vertex_color(Vertex* v, int r, int g, int b){
    v->r = r;
    v->g = g;
    v->b = b;
    return;
  }

  /**
   * @brief Create two halfedges (twins) between vertex v1 and v2,
   * will only have origin and twin attributes.
   * 
   * @param v1 First vertex of halfedge.
   * @param v2 Second vertex of halfedge.
   * @return Halfedge* one of the two created, the other can be acessed by twin attribute.
   */
  Halfedge* halfedges_between(Vertex* v1, Vertex* v2){
    Halfedge* he1 = new Halfedge();
    Halfedge* he2 = new Halfedge();
    he1->origin = v1;   he2->origin = v2;
    he1->twin = he2;    he2->twin = he1;
    return he1;
  }

  /**
   * @brief Create a triangle structure between the three halfedges recieved
   * in the clockwise direction.
   * 
   * @param h1 First halfedge
   * @param h2 Second halfedge
   * @param h3 Third halfedge
   */
  void create_triangle(Halfedge* h1, Halfedge* h2, Halfedge* h3){
    h1->next = h2;
    h2->next = h3;
    h3->next = h1;
    return;
  }

  /**
   * @brief Create the initial triangle of the structure with helper points 
   * that contain all points that will be used and link all created halfedges to a new TriangleNode.
   * 
   * @param M value for biggest coordinate
   * @param tRoot Pointer to PointLocation TriangleNode
   * @return Halfedge* pointer to be add to PointLocation.
   */
  Halfedge* start_structure(double M, TriangleNode* tRoot){
    Vertex * v_minus_3 = new Vertex();
    v_minus_3->index = -3;
    set_vertex_color(v_minus_3, 0, 0, 0);
    v_minus_3->p = Point (-3*M, -3*M);
    Vertex* v_minus_2 = new Vertex();
    v_minus_2->index = -2;
    set_vertex_color(v_minus_2, 0, 0, 0);
    v_minus_2->p = Point (0, 3*M);
    Vertex* v_minus_1 = new Vertex();
    v_minus_1->index = -1;
    set_vertex_color(v_minus_1, 0, 0, 0);
    v_minus_1->p = Point (3*M, 0);

    vertices.push_back(v_minus_3);
    vertices.push_back(v_minus_2);
    vertices.push_back(v_minus_1);

    Halfedge* h1 = halfedges_between(v_minus_3, v_minus_2);
    Halfedge* h2 = halfedges_between(v_minus_2, v_minus_1);
    Halfedge* h3 = halfedges_between(v_minus_1, v_minus_3);

    create_triangle(h1, h2, h3);
    create_triangle(h3->twin, h2->twin, h1->twin);

    h1->tri = tRoot;  h1->twin->tri = tRoot;
    h2->tri = tRoot;  h2->twin->tri = tRoot;
    h3->tri = tRoot;  h3->twin->tri = tRoot;

    halfedges.push_back(h1);   
    halfedges.push_back(h2); 
    halfedges.push_back(h3); 

    return h1;
  }

  /**
   * @brief Add a vertex inside the triangle that contains the halfedge he,
   * and for each vertex of this triangle, create a new halfedge to it.
   * 
   * @param vr Vertex that will be linked to every other vertex inside this triangle.
   * @param he1 One halfedge that is inside of the triangle.
   * @return std::vector<Halfedge*> pointer to halfedges of the new triangles.
   */
  std::vector<Halfedge*>  add_center_vertex(Vertex* vr, Halfedge* he1){
    add_vertex(vr);
    Halfedge* he2 = he1->next;
    Halfedge* he3 = he2->next;

    Vertex* v1 = he1->origin;
    Vertex* v2 = he2->origin;
    Vertex* v3 = he3->origin;

    Halfedge* new_he1 = halfedges_between(vr, v1);
    Halfedge* new_he2 = halfedges_between(vr, v2);
    Halfedge* new_he3 = halfedges_between(vr, v3);

    create_triangle(new_he1, he1, new_he2->twin);
    create_triangle(new_he2, he2, new_he3->twin);
    create_triangle(new_he3, he3, new_he1->twin);

    halfedges.push_back(new_he1);
    halfedges.push_back(new_he2); 
    halfedges.push_back(new_he3); 

    std::vector<Halfedge*> newTri_he {new_he1, new_he2, new_he3};
    return newTri_he;
  }

  /**
   * @brief Checks if the edge i-j, that belon go the triangles ijk and irj,
   * with r being the newest point added, is illegal. In the general case, to check if is illegal,
   * a circle is created with the vertex v_i, v_r and v_j (origin of edges
   * i-r, r-j and j-i respectively), and is check if the vertex v_k is inside this circle.
   * 
   * @param vi Vertex* origin of edge i-r;
   * @param vj Vertex* origin of edge j-i;
   * @param vr Vertex* origin of edge r-j;
   * @param vk Vertex* origin of edge k-i;
   * @return true 
   * @return false 
   */
  bool is_illegal(Vertex* vi, Vertex* vj, Vertex* vr, Vertex*vk){
    if((vi->index < 0) & (vj->index < 0)){
      return false;
    }else{
      Point origin = CGAL::circumcenter(vi->p, vj->p, vr->p);
      double radium = (vi->p.x() - origin.x())*(vi->p.x() - origin.x()) + 
                      (vi->p.y() - origin.y())*(vi->p.y() - origin.y());
      double dist = (vk->p.x() - origin.x())*(vk->p.x() - origin.x()) + 
                    (vk->p.y() - origin.y())*(vk->p.y() - origin.y());
      return dist < radium;
    }
  }

  /**
   * @brief Flip the edge i-j (both halfedges),
   * fix the pointers of the others halfedges and add new childs to PointLocation.
   * Return the changed halfedges that are needed to check if is illegal.
   *        k                k
   *      /   \            / | \
   *     /     \          /  |  \
   *    /       \        /   |   \
   *  j --------- i    j     |    i
   *    \       /        \   |   /
   *     \     /          \  |  /
   *      \   /            \ | /
   *        r                r
   * 
   * @param he_ij 
   * @return std::vector<Halfedge*> halfedges that are needed to be updated.
   */
  std::vector<Halfedge*> flip(Halfedge* he_ji){
    Halfedge *he_ir, *he_rj, *he_ij, *he_jk, *he_ki;
    he_ir = he_ji->next;
    he_rj = he_ir->next;
    he_ij = he_ji->twin;
    he_jk = he_ij->next;
    he_ki = he_jk->next;

    TriangleNode *tNode1, *tNode2, *tNode1_new, *tNode2_new;
    tNode1 = he_ji->tri;
    tNode2 = he_ij->tri;

    he_ji->origin = he_rj->origin;
    he_ij->origin = he_ki->origin;
    create_triangle(he_ji, he_ki, he_ir);
    create_triangle(he_ij, he_rj, he_jk);

    tNode1_new = new TriangleNode(he_ji);
    tNode2_new = new TriangleNode(he_ij);

    int max_depth = tNode1->depth > tNode2->depth ? tNode1->depth : tNode2->depth;

    tNode1_new->depth = max_depth + 1;
    tNode2_new->depth = max_depth + 1;
    tNode1->childs.push_back(tNode1_new);
    tNode1->childs.push_back(tNode2_new);
    tNode2->childs.push_back(tNode1_new);
    tNode2->childs.push_back(tNode2_new); 

    he_ji->tri = tNode1_new;    he_ij->tri = tNode2_new;
    he_ki->tri = tNode1_new;    he_rj->tri = tNode2_new;
    he_ir->tri = tNode1_new;    he_jk->tri = tNode2_new;

    std::vector<Halfedge*> updated_he = {he_ki, he_jk};

    return updated_he;
  }

  /**
   * @brief Check if the halfedge is legal, and if it isn't, call the function to flip the halfedge.
   * 
   * @param vr Vertex of the new created point.
   * @param he Halfedge to be checked if is legal.
   */
  void legalize_edge(Vertex* vr, Halfedge* he){
    Vertex* vi = he->origin;
    Vertex* vj = he->twin->origin;
    Vertex* vk = he->twin->next->next->origin;
    if(is_illegal(vi, vj, vr, vk)){
        std::vector<Halfedge*> updated_he = flip(he);
        std::vector<Halfedge*>::iterator it = updated_he.begin();
        for(; it != updated_he.end(); it++){
          legalize_edge(vr, *it);
        }
    }
  }
};

class PointLocation{
  public:
    TriangleNode* root = new TriangleNode();

    /**
     * @brief Recursive print function of the nodes of this structure.
     * 
     */
    void print(){
      std::cout << std::endl;
      std::cout << "PointLocation:"  << std::endl;
      root->print();
      if(root->childs.empty()){
        std::cout << "/\\" << std::endl;
      }else{
        std::vector<TriangleNode*>::iterator it = root->childs.begin();
        for(; it != root->childs.end(); it++){
          print(*it);
        }
      }
      std::cout << std::endl;
      return;
    }

    /**
     * @brief Recursive print function of the nodes of this structure.
     * 
     * @param tNode node to print and recurse for childs.
     */
    void print(TriangleNode* tNode){
      tNode->print();
      if(tNode->childs.empty()){
        std::cout << "/\\" << std::endl;
      }else{
        std::vector<TriangleNode*>::iterator it = tNode->childs.begin();
        for(; it != tNode->childs.end(); it++){
          print(*it);
        }
      }
    }

    /**
     * @brief Update the root TriangleNode with the halfedge, the array of vertex and the depth.
     * 
     * @param he 
     */
    void update_root(Halfedge* he){
      root->depth = 0;
      root->he = he;
      root->nodeVertex[0] = he->origin->copy();
      root->nodeVertex[1] = he->next->origin->copy();
      root->nodeVertex[2] = he->next->next->origin->copy();
      root->t = Triangle(root->nodeVertex[0].p,
                           root->nodeVertex[1].p,
                           root->nodeVertex[2].p);
      return;
    }

    /**
     * Search for triangle node that contain point and returns it,
     * if there isn't a triangle that contains, return a null pointer.
     *
     * @param tNode TriangleNode to search childs.
     * @param p Point to seach.
     * @return Pointer to TriangleNode or null pointer.
     */
    TriangleNode* search(TriangleNode* tNode, Point p){
      if(tNode->childs.empty()){
        if(tNode->contains_point(p)){
          return tNode;
        }else{
          return nullptr;
        }
      }
      
      std::vector<TriangleNode*>::iterator it = tNode->childs.begin();
      for (; it != tNode->childs.end(); ++it){
        if((*it)->contains_point(p)){
            return search(*it, p);
        }
      }
      return nullptr;
    }

    /**
     * Search for triangle node that contain point and returns it,
     * if there isn't a triangle that contains, return a null pointer.
     *
     * @param p Point to seach.
     * @return Pointer to TriangleNode or null pointer.
     */
    TriangleNode* search(Point p){
      return search(root, p);
    }    

    /**
     *  Append childs to a TriangleNode.
     * 
     * @param tNode TriangleNode that will recieve childs.
     * @param childtNodes Vector of TriangleNodes* that are the new childs.
     */
    void insert(TriangleNode* tNode, std::vector<TriangleNode*> childtNodes){
      std::vector<TriangleNode*>::iterator it = childtNodes.begin();
      for (; it != childtNodes.end(); ++it){
        tNode->childs.push_back(&**it);
      }
      return;
    }

    /**
     * Append childs to the root node.
     * 
     * @param childtNodes Vector of TriangleNodes* that are the new childs.
     */
    void insert(std::vector<TriangleNode*> childtNodes){
      std::vector<TriangleNode*>::iterator it = childtNodes.begin();
      for (; it != childtNodes.end(); ++it){
        root->childs.push_back(&**it);
      }
      return;
    }

    /**
     * @brief Set the childs of a TriangleNode from the vector of halfedges.
     * 
     * @param leafTri leaf TriangleNode to get more childs.
     * @param newTri_he vector of halfedges of the new childs.
     */
    void add_childs(TriangleNode* leafTri, std::vector<Halfedge*> newTri_he){
      std::vector<TriangleNode*> childs;
      std::vector<Halfedge*>::iterator it = newTri_he.begin();
      TriangleNode* curTri;
      for(; it != newTri_he.end(); it++){
        curTri = new TriangleNode(*it);
        (*it)->tri = curTri;
        (*it)->next->tri = curTri;
        (*it)->next->next->tri = curTri;
        curTri->depth = leafTri->depth + 1;
        childs.push_back(curTri);
      }
      leafTri->childs = childs;
      return; 
    }
};

class Delaunay{
  public:
    std::vector<Point> points;
    std::vector<int> r_v, g_v, b_v;
    DCEL T;
    PointLocation D;

    /** 
     * @brief Load points from a file with two coordinates by line separated by space.
     * 
     */
    void load_file(std::string filename){
      std::ifstream points_file(filename);
      std::string point_line;
      while(getline(points_file, point_line)){
        std::stringstream line_stream(point_line);
        std::string x_s, y_s, r_s, g_s, b_s;
        double x, y, r, g, b;
        line_stream >> x_s >> y_s >> r_s >> g_s >> b_s;
        x = std::stod(x_s);   
        y = std::stod(y_s);
        r = std::stod(r_s);
        g = std::stod(g_s);
        b = std::stod(b_s);
        points.push_back(Point (x, y));
        r_v.push_back(r);
        g_v.push_back(g);
        b_v.push_back(b);
      }
      points_file.close();
      return;
    }

    /**
     * @brief Find the biggest value of the coordinates of points.
     * 
     * @return M double with value.
     */
    double biggest_coordinate(){
      std::vector<Point>::iterator it = points.begin();
      double M = std::abs((*it).x());
      for(; it != points.end(); ++it){
        M = M <= std::abs((*it).x()) ? std::abs((*it).x()) : M;
        M = M <= std::abs((*it).y()) ? std::abs((*it).y()) : M;
      }
      return M;
    }

    /**
     * @brief Start DCEL and PointLocation structures with the initial point.
     * 
     * @param M biggest coordinate.
     */
    void start_structures(double M){
      TriangleNode* tRoot = D.root;
      Halfedge* he = T.start_structure(M, tRoot);
      D.update_root(he);
      return;
    }

    void run(){
      double M = biggest_coordinate();
      start_structures(M);
      
      auto rng = std::default_random_engine {};
      std::shuffle(std::begin(points), std::end(points), rng);
      std::vector<Point>::iterator it = points.begin();
      std::vector<Point>::iterator it_r = r_v.begin();
      std::vector<Point>::iterator it_g = g_v.begin();
      std::vector<Point>::iterator it_b = b_v.begin();

      for(; it != points.end(); ++it){
        TriangleNode* leafTri = D.search(*it);
      
        Vertex* vr = new Vertex(*it);
        vr->r = (*it_r);
        vr->g = (*it_g);
        vr->b = (*it_b);
        std::vector<Halfedge*> newTri_he;

        newTri_he = T.add_center_vertex(vr, leafTri->he);
        D.add_childs(leafTri, newTri_he);

        std::vector<Halfedge*>::iterator it2 = newTri_he.begin();
        
        for(; it2 != newTri_he.end(); it2++){
          T.legalize_edge(vr, (*it2)->next);
        }
        it_r++; it_g++; it_b++;        
      }

      T.print_final_edges();
      T.save_triangulation();
      return;
    }
};

int main() {
  Delaunay de;
  de.load_file("points.txt");
  de.run();
  return 0;
}
