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
  bool symbol_left = false;
  bool symbol_right = false;

  Vertex(){}

  Vertex(Point_2 new_p){
    p = new_p;
  }
};

struct Halfedge{
  Halfedge* twin;
  Halfedge* prev;
  Halfedge* next;
  Vertex* origin;

  Halfedge(){}
};

class DCEL{
  public:
  std::vector<Vertex*> vertices;
  std::vector<Halfedge*> halfedges;

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
    he1->origin = p1;
    he1->twin = he2;
    he2->origin = p2;
    he2->twin = he1;
    return he1;
  }

  /**
   * @brief Create the initial triangle of the structure with
   * the point p_left and p_right that are points with coordinates 
   * (-inf, max(points.y)), (inf, min(points.y)) and the rightmost highest
   * point between points.
   * 
   * @param p Rightmost highest Point_2.
   */
  void start_structure(Point_2 p){
    Vertex* v_p0 = new Vertex(p);
    Vertex* v_left = new Vertex();
    v_left->symbol_left = true;
    Vertex* v_right = new Vertex();
    v_right->symbol_right = true;

    vertices.push_back(v_left);
    vertices.push_back(v_p0);
    vertices.push_back(v_right);

    Halfedge* h1 = halfedges_between(v_left, v_p0);
    Halfedge* h2 = halfedges_between(v_p0, v_right);
    Halfedge* h3 = halfedges_between(v_right, v_left);

    Halfedge* h1_twin = h1->twin;
    Halfedge* h2_twin = h2->twin;
    Halfedge* h3_twin = h3->twin;

    h1->next = h2;
    h2->next = h3;
    h3->next = h1;

    h1->prev = h3;
    h2->prev = h1;
    h3->prev = h2;

    h3_twin->next = h2_twin;
    h2_twin->next = h1_twin;
    h1_twin->next = h3_twin;
    
    h3_twin->prev = h1_twin;
    h2_twin->prev = h3_twin;
    h1_twin->prev = h2_twin;

    halfedges.push_back(h1);
    halfedges.push_back(h2);
    halfedges.push_back(h3);
    halfedges.push_back(h1_twin);
    halfedges.push_back(h2_twin);
    halfedges.push_back(h3_twin);

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
    h1->next = h2;    h2->prev = h1;
    h2->next = h3;    h3->prev = h2;
    h3->next = h1;    h1->prev = h3;
    return;
  }


  /**
   * @brief dd a point inside the triangle that contains the halfedge he,
   * and for each vertex of this triangle, create a new halfedge to the center point.
   * 
   * @param p Point that will be linked to every other vertex inside this triangle
   * @param he1 One halfedge that is inside of the triangle
   * @return std::vector<TriangleNode*> to be appended on the PointLocation leaf
   */
  std::vector<TriangleNode*>  add_center_point(Point_2 p, Halfedge* he1){
    Vertex* v_p = new Vertex(p);
    Halfedge* he2 = he1->next;
    Halfedge* he3 = he2->next;

    Vertex* p1 = he1->origin;
    Vertex* p2 = he2->origin;
    Vertex* p3 = he3->origin;

    Halfedge* new_he1 = halfedges_between(v_p, p1);
    Halfedge* new_he2 = halfedges_between(v_p, p2);
    Halfedge* new_he3 = halfedges_between(v_p, p3);

    create_triangle(new_he1, he1, new_he2->twin);
    create_triangle(new_he2, he2, new_he3->twin);
    create_triangle(new_he3, he3, new_he1->twin);
    
    std::vector<TriangleNode*> childs;
    childs.push_back(new_he1);
    childs.push_back(new_he2);
    childs.push_back(new_he3);

    return childs;

  }
};

struct TriangleNode{
  std::vector<TriangleNode*> childs;
  Triangle_2 t;
  Halfedge* he;
  bool root_node;
  Point_2 p0;
  
  TriangleNode(){
    root_node = true;
  }

  void set_triangle(Triangle_2 new_triangle){
    t = new_triangle;
    root_node = false;
  }

  bool contains_point(Point_2 p){
    if(root_node){
      return p.y() <= p0.y();
    }else{
      return t.bounded_side(p) == CGAL::ON_BOUNDED_SIDE;
    }
  }
};

class PointLocation{
  public:
    TriangleNode* root = new TriangleNode();

    void start_structure(Point_2 p0, Halfedge* he){
      root->p0 = p0;
      root->he = he;
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

    void run(){
      //Find the righmost highest point and remove it
      std::vector<Point_2>::iterator it_p0 = rightmost_highest();
      Point_2 p0 = *it_p0;
      std::cout << "The rightmost highest point is: " << p0 << std::endl;
      points.erase(it_p0);

      //Start DCEL
      T.start_structure(p0);
      //Start PointLocation
      D.start_structure(p0, T.halfedges.front());
      

      //Random order of the remaining points
      auto rng = std::default_random_engine {};
      std::shuffle(std::begin(points), std::end(points), rng);
      std::vector<Point_2>::iterator it = points.begin();
      for(; it != points.end(); ++it){
        std::cout << "Searching for the point: " << *it << std::endl;
        
        
        //find triangle in the point location that contains point
        TriangleNode* leafTriang = D.search(*it);
        
        
        
        if(leafTriang != nullptr){
          std::cout << "Triangle found." << std::endl;
        }

        //create new edges on the DCEL
        if(leafTriang->t.bounded_side(*it) == CGAL::ON_BOUNDED_SIDE){
          std::cout << "Point inside a triangle." << std::endl;
          std::vector<TriangleNode*> childs = T.add_center_point(*it, leafTriang->he);
          D.insert(leafTriang, childs);

        }else if(leafTriang->t.bounded_side(*it) == CGAL::ON_BOUNDARY){
          std::cout << "Point on a boundary of a triangle." << std::endl;
        }
        
        

        //legalize wrong edges
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
