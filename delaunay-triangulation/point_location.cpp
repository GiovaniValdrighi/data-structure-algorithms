#include <iostream>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/HalfedgeDS_vertex_base.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <algorithm>
#include <random>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2; 
typedef Kernel::Triangle_2 Triangle_2;

struct Halfedge; struct TriangleNode;

struct Vertex{
  Point_2 p;
  Halfedge* he;
  int index;

  Vertex(){}

  Vertex(Point_2 new_p){
    p = new_p;
  }

  /**
   * @brief Print the information of the vertex.
   * 
   */
  void print(){
    if(index == -2){
      std::cout << " Left p_{-2} vertex ";
    }else if(index == -1){
      std::cout << " Right p_{-1} vertex ";
    }else{
      std::cout << " (" << p << ")" << " vertex";
    }
    return;
  }

};

struct Halfedge{
  Halfedge* twin;
  Halfedge* prev;
  Halfedge* next;
  Vertex* origin;
  TriangleNode* tri;

  Halfedge(){}
};

struct TriangleNode{
  std::vector<TriangleNode*> childs;
  Triangle_2 t;
  Vertex nodeVertex [3];
  Halfedge* he;
  Point_2 p0;
  int depth;
  
  TriangleNode(){ }

  TriangleNode(Halfedge* new_he){
    he = new_he;
    Point_2 p1 = new_he->origin->p;
    Point_2 p2 = new_he->next->origin->p;
    Point_2 p3 = new_he->next->next->origin->p;

    Triangle_2 new_triangle (p1, p2, p3);
    set_triangle(new_triangle);
    
  }

  void set_triangle(Triangle_2 new_triangle){
    t = new_triangle;
  }

  /**
   * @brief Print the vertex index of this TriangleNode.
   * 
   */
  void print(){
    if(he != nullptr){
      std::cout << "T: " << he->origin->index << " ";
      std::cout << he->next->origin->index << " ";
      std::cout << he->next->next->origin->index; 
      std::cout << " d: " << depth << std::endl;
    }else{
      std::cout << "T: ---" << std::endl;
    }
    return;
  }

  /**
   * @brief Return true if p is on the right side of the segment p1-p2.
   * 
   * @param v1 
   * @param v2 
   * @param p 
   * @return true 
   * @return false 
   */
  bool isright(Vertex* v1, Vertex* v2, Point_2 p){
    if((v1->index >= 0) & (v2->index >= 0)){
      return CGAL::right_turn(v1->p, v2->p, p);
    }else if((v1->index == -1) & (v2->index == -2)){
      return true;
    }else if(v1->index == -1){
      return (v2->p.x() < p.x()) & (v2->p.y() <= p.y());
    }else if(v2->index == -1){
      return !((v1->p.x() < p.x()) & (v1->p.y() <= p.y()));
    }else if(v1->index == -2){
      return !((v2->p.x() < p.x()) & (v2->p.y() <= p.y()));
    }else if(v2->index == -2){
      return (v1->p.x() < p.x()) & (v1->p.y() <= p.y());
    }
    return false;
  }

  /**
   * @brief Verify if a point is inside the triangle of the TriangleNode.
   * 
   * @param p 
   * @return true 
   * @return false 
   */
  bool contains_point(Point_2 p){
    Vertex *v1, *v2, *v3;
    v1 = he->origin;
    v2 = he->next->origin;
    v3 = he->next->next->origin;
    if((v1->index >= 0) & 
       (v2->index >= 0) & 
       (v3->index >= 0)){
      return t.bounded_side(p) == CGAL::ON_BOUNDED_SIDE;
    }else{
      bool isright_1 = isright(v1, v2, p);
      bool isright_2 = isright(v2, v3, p);
      bool isright_3 = isright(v3, v1, p);
      return isright_1 & isright_2 & isright_3;
    }
  }
};


class DCEL{
  public:
  std::vector<Vertex*> vertices;
  std::vector<Halfedge*> halfedges;

  /**
   * @brief Print the vertexs of each edge (the origin of the halfedge and the origin of the twin).
   * 
   */
  void print_edges(){
    std::cout << std::endl;
    std::cout << "DCEL halfedges: " << std::endl;
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
   * @brief Create two halfedges (twins) between vertex p1 and p2,
   * they will not contain prev e next attributes, only twin and origin.
   * 
   * @param p1 First vertex of halfedge.
   * @param p2 Second vertex of halfedge.
   * @return Halfedge* one of the two created, the other can be acessed by twin attribute.
   */
  Halfedge* halfedges_between(Vertex* p1, Vertex* p2){
    Halfedge* he1 = new Halfedge();
    Halfedge* he2 = new Halfedge();
    he1->origin = p1;   he2->origin = p2;
    he1->twin = he2;    he2->twin = he1;
    return he1;
  }

  /**
   * @brief Set the prev and next attribute of two halfedges.
   * 
   * @param hep Previous halfedge.
   * @param hen Next halfedge.
   */
  void set_prev_next(Halfedge *hep, Halfedge* hen){
    hep->next = hen;    //hen->next = hep;
    return;
  }

  /**
   * @brief Create a triangle structure between the three halfedges recieve
   * in the clockwise direction
   * 
   * @param h1 First halfedge
   * @param h2 Second halfedge
   * @param h3 Third halfedge
   */
  void create_triangle(Halfedge* h1, Halfedge* h2, Halfedge* h3){
    set_prev_next(h1, h2);
    set_prev_next(h2, h3);
    set_prev_next(h3, h1);
    return;
  }

  /**
   * @brief Create the initial triangle of the structure with
   * the point p_left and p_right that are points with coordinates 
   * (-inf, max(points.y)), (inf, min(points.y)) and the rightmost highest
   * point between points and link all halfedges to a TriangleNode.
   * 
   * @param p Rightmost highest Point_2.
   * @param tRoot Pointer to PointLocation TriangleNode
   * @return Halfedge* pointer to be add to PointLocation.
   */
  Halfedge* start_structure(Point_2 p, TriangleNode* tRoot){
    Vertex* v_0 = new Vertex(p);
    v_0->index = 0;
    Vertex* v_minus_2 = new Vertex();
    v_minus_2->index = -2;
    Vertex* v_minus_1 = new Vertex();
    v_minus_1->index = -1;

    vertices.push_back(v_minus_2);
    vertices.push_back(v_minus_1);
    vertices.push_back(v_0);

    Halfedge* h1 = halfedges_between(v_minus_2, v_0);
    Halfedge* h2 = halfedges_between(v_0, v_minus_1);
    Halfedge* h3 = halfedges_between(v_minus_1, v_minus_2);

    std::cout << "h1 origin: " << h1->origin->index << std::endl;
    std::cout << "h2 origin: " << h2->origin->index << std::endl;
    std::cout << "h3 origin: " << h3->origin->index << std::endl;

    Halfedge* h1_twin = h1->twin;
    Halfedge* h2_twin = h2->twin;
    Halfedge* h3_twin = h3->twin;

    create_triangle(h1, h2, h3);
    create_triangle(h3_twin, h2_twin, h1_twin);

    std::cout << "check if triangle is right:" << std::endl;
    std::cout << "first next: " << h1->next->origin->index << h2->origin->index << std::endl;
    std::cout << "second next: " << h1->next->next->origin->index << h3->origin->index << std::endl;


    h1->tri = tRoot;
    h2->tri = tRoot;
    h3->tri = tRoot;

    halfedges.push_back(h1);
    halfedges.push_back(h2);
    halfedges.push_back(h3);

    std::cout << "after some operations" << std::endl;
    std::cout << "h1 origin: " << h1->origin->index << std::endl;
    std::cout << "h2 origin: " << h2->origin->index << std::endl;
    std::cout << "h3 origin: " << h3->origin->index << std::endl;

    return h1;
  }

  /**
   * @brief dd a point inside the triangle that contains the halfedge he,
   * and for each vertex of this triangle, create a new halfedge to the center point.
   * 
   * @param pr Vertex that will be linked to every other vertex inside this triangle.
   * @param he1 One halfedge that is inside of the triangle.
   * @return std::vector<Halfedge*> pointer to halfedges of the new triangles.
   */
  std::vector<Halfedge*>  add_center_point(Vertex* pr, Halfedge* he1){
    pr->index = vertices.back()->index + 1;
    vertices.push_back(pr);
    Halfedge* he2 = he1->next;
    Halfedge* he3 = he2->next;

    Vertex* p1 = he1->origin;
    Vertex* p2 = he2->origin;
    Vertex* p3 = he3->origin;

    std::cout << "p1 index: " << p1->index << std::endl;
    std::cout << "p2 index: " << p2->index << std::endl;
    std::cout << "p3 index: " << p3->index << std::endl;

    Halfedge* new_he1 = halfedges_between(pr, p1);
    Halfedge* new_he2 = halfedges_between(pr, p2);
    Halfedge* new_he3 = halfedges_between(pr, p3);

    halfedges.push_back(new_he1);
    halfedges.push_back(new_he2);
    halfedges.push_back(new_he3);

    std::cout << "Created edges between." << std::endl;

    create_triangle(new_he1, he1, new_he2->twin);
    create_triangle(new_he2, he2, new_he3->twin);
    create_triangle(new_he3, he3, new_he1->twin);
    
    std::cout << "Created triangles of edges." << std::endl;

    std::vector<Halfedge*> newTri_he;

    newTri_he.push_back(new_he1);
    newTri_he.push_back(new_he2);
    newTri_he.push_back(new_he3);

    return newTri_he;

  }

  /**
   * @brief Checks if the edge i-j, that belon go the triangles ijk and irj,
   * with r being the newest point added, is illegal. In the general case, to check if is illegal,
   * a circle is created with the vertex p_i, p_r and p_j (origin of edges
   * i-r, r-j and j-i respectively), and is check if the vertex p_j is inside this circle.
   * 
   * @param pi Vertex* origin of edge i-r;
   * @param pj Vertex* origin of edge j-i;
   * @param pr Vertex* origin of edge r-j;
   * @param pk Vertex* origin of edge k-i;
   * @return true 
   * @return false 
   */
  bool is_illegal(Vertex* pi, Vertex* pj, Vertex* pr, Vertex*pk){
    if(
      ((pi->index == -2) & (pj->index == 0)) |
      ((pi->index == 0) & (pj->index == -1)) |
      ((pi->index == -1) & (pj->index == -2))
    ){
      return false;
    }else if(
      (pi->index >= 0) & 
      (pj->index >= 0) &
      (pr->index >= 0) &
      (pk->index >= 0)
    ){
      Point_2 origin = CGAL::circumcenter(pi->p, pj->p, pr->p);
      double radium = (pi->p.x() - origin.x())*(pi->p.x() - origin.x()) + (pi->p.y() - origin.y())*(pi->p.y() - origin.y());
      double dist = (pk->p.x() - origin.x())*(pk->p.x() - origin.x()) + (pk->p.y() - origin.y())*(pk->p.y() - origin.y());
      return dist < radium;
    }else{
      int min_kr = pk->index;
      if(pr->index < pk->index) min_kr = pr->index;

      int min_ij = pj->index;
      if(pi->index < pj->index) min_ij = pi->index;

      return min_kr < min_ij;
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

    //create two new trainglenodes
    //link the old to the new triangle nodes
    //link the halfedges to the new trianglenodes
    create_triangle(he_ji, he_ki, he_ir);
    create_triangle(he_ij, he_rj, he_jk);

    tNode1_new = new TriangleNode(he_ji);
    tNode2_new = new TriangleNode(he_ij);
    std::vector<TriangleNode*> childs {tNode1_new, tNode2_new};
    tNode1->childs = childs;
    tNode2->childs = childs;    

    he_ji->tri = tNode1_new;
    he_ki->tri = tNode1_new;
    he_ir->tri = tNode1_new;

    he_ij->tri = tNode2_new;
    he_rj->tri = tNode2_new;
    he_jk->tri = tNode2_new;

    std::vector<Halfedge*> updated_he;
    updated_he.push_back(he_ki);
    updated_he.push_back(he_jk);
    return updated_he;
  }


  void legalize_edge(Vertex* pr, Halfedge* he){
    std::cout << "Halfedge to be legalized: " << he->origin->index << "->" << he->twin->origin->index << std::endl;
    Vertex* pi = he->origin;
    Vertex* pj = he->twin->origin;
    Vertex* pk = he->twin->next->next->origin;
    std::cout << "p_i: ";
    pi->print();
    std::cout << " p_j: ";
    pj->print();
    std::cout << " p_k: ";
    pk->print();
    std::cout << std::endl;
    if(is_illegal(pi, pj, pr, pk)){
        std::cout << "Halfedge is illegal." << std::endl;
        std::vector<Halfedge*> updated_he = flip(he);
        std::vector<Halfedge*>::iterator it = updated_he.begin();
        for(; it != updated_he.end(); it++){
          legalize_edge(pr, *it);
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
     * @brief Set the top point of the initial TriangleNode and return it.
     * 
     * @param top point
     * @return TriangleNode* root with point set
     */
    TriangleNode* update_root(Halfedge* he){
      root->depth = 0;
      return root;
    }

    /**
     * Search for triangle node that contain point and returns it,
     * if there isn't a triangle that contains, return a null pointer.
     *
     * @param tNode TriangleNode to search childs.
     * @param p Point to seach.
     * @return Pointer to TriangleNode or null pointer.
     */
    TriangleNode* search(TriangleNode* tNode, Point_2 p){
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
    TriangleNode* search(Point_2 p){
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
    void set_childs(TriangleNode* leafTri, std::vector<Halfedge*> newTri_he){
      //need to create the triangles from cgal
      //create the triangle nodes with this triangles
      //set the halfedges
      //and update leaftri childs
      std::vector<TriangleNode*> childs;
      std::vector<Halfedge*>::iterator it = newTri_he.begin();
      TriangleNode* curTri;
      for(; it != newTri_he.end(); it++){
        curTri = new TriangleNode(*it);
        curTri->depth = leafTri->depth + 1;
        childs.push_back(curTri);
      }
      leafTri->childs = childs;
      return; 
    }
};

class Delaunay{
  public:
    std::vector<Point_2> points;
    DCEL T;
    PointLocation D;

    void add(Point_2 p){
      points.push_back(p);
      std::cout << "Add point, vector length:" << points.size() << std::endl;
      return;
    }

    /**
     * Compute the point with the highest y-value, and if there is more than one,
     * the point with highest x-value.
     * 
     * @return std::vector<Point_2>::iterator iterator pointing to the point in the vector.
     */
    std::vector<Point_2>::iterator  rightmost_highest(){
      std::vector<Point_2>::iterator it_min = points.begin();
      std::vector<Point_2>::iterator it = points.begin();

      for(; it != points.end(); ++it){
        if((*it_min).y() < (*it).y()){
          it_min = it;
        }else if((*it_min).y() == (*it).y()){
          if((*it_min).x() < (*it).x()){
            it_min = it;
          }
        }
      }
      return it_min;
    }

    /**
     * @brief Start DCEL and PointLocation structures with the initial point.
     * 
     * @param p0 initial point.
     */
    void start_structures(Point_2 p0){
      //create a triangle in the dcel with two symbolic points
      //total of 3 halfedges
      //create a triangle node in the pointlocation (the root)
      
      //every halfedge with a link to this trianglenode
      //the trianglenode linked to one halfedge

      TriangleNode* tRoot = D.root;//D.start_structure();
      Halfedge* he = T.start_structure(p0, tRoot);
      D.update_root(he);
      T.print_edges();
      D.print();
      return;
    }

    void run(){
      int counter = 0;
      //Find the righmost highest point and remove it
      std::vector<Point_2>::iterator it_p0 = rightmost_highest();
      Point_2 p0 = *it_p0;
      std::cout << "The rightmost highest point is: " << p0 << std::endl;
      points.erase(it_p0);

      //start structures
      start_structures(p0);
      
      //Random order of the remaining points
      auto rng = std::default_random_engine {};
      std::shuffle(std::begin(points), std::end(points), rng);
      std::vector<Point_2>::iterator it = points.begin();
      for(; it != points.end(); ++it){
        std::cout << "Searching for the point: " << *it << std::endl;
        
        //find triangle in the point location that contains point
        TriangleNode* leafTri = D.search(*it);
        
        if(leafTri != nullptr){
          std::cout << "Triangle found." << std::endl;
          leafTri->print();
        }else{
          std::cout << "Triangle not found." << std::endl;
        }

        Vertex* pr = new Vertex(*it);
        //create new edges on the DCEL
        std::vector<Halfedge*> newTri_he = T.add_center_point(pr, leafTri->he);
        T.print_edges();

        //update links in the PointLocation
        D.set_childs(leafTri, newTri_he);
        D.print();

        std::vector<Halfedge*>::iterator it2 = newTri_he.begin();
        
        for(; it2 != newTri_he.end(); it2++){
          T.legalize_edge(pr, (*it2)->next);
          T.print_edges();
          D.print();
        }
        
        counter = counter + 1;
        if(counter >= 2){
          break;
        }
      }
      return;
    }
};

int main() {
  Delaunay delau;
  delau.add(Point_2 (1, 0));
  delau.add(Point_2 (2, 0.5));
  delau.add(Point_2 (0, 0.5));
  delau.add(Point_2 (1, 1));

  delau.run();
  return 0;
}
