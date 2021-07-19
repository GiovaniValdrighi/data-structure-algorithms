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

struct Edge;

struct Vertex{
  Point_2 point;
  Edge* incident;

  Vertex(){}

  Vertex(Point_2 new_p){
    point = new_p;
  }
};

struct Edge{
  Edge* twin;
  Vertex* next;

  Edge(){}
};

class DCEL{
  public:
    std::vector<Vertex*> vertices;
    std::vector<Edge*> edges;
};

struct TriangleNode{
  std::vector<TriangleNode*> childs;
  Triangle_2 t;
  
  TriangleNode(){}

  TriangleNode(Triangle_2 new_t){
    t = new_t;
  }
};

class PointLocation{
  public:
    TriangleNode* root = new TriangleNode();

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
        if(tNode->t.bounded_side(p) == CGAL::ON_BOUNDED_SIDE){
          return tNode;
        }else{
          return nullptr;
        }
      }
      
      std::vector<TriangleNode*>::iterator it = tNode->childs.begin();
      for (; it != tNode->childs.end(); ++it){
        if((*it)->t.bounded_side(p) == CGAL::ON_BOUNDED_SIDE){
            return search(*it, p);
        }
      }
      return nullptr;
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
    DCEL* T;
    PointLocation* D;

    void insert(Point_2 p){
      points.push_back(p);
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
          it_min = &(*it);
        }else if((*it_min).y() == (*it).y()){
          if((*it_min).x() < (*it).x()){
            it_min = &(*it);
          }
        }
      }
      return it_min;
    }

    void run(){
      //Find the righmost highest point and remove it
      std::vector<Point_2>::iterator it_p0 = rightmost_highest();
      Point_2 p0 = *it;
      points.erase(it_p0);

      //Create root triangle

      //Random order of the remaining points
      auto rng = std::default_random_engine {};
      std::shuffle(std::begin(points), std::end(points), rng);
      std::vector<Point_2>::iterator it = points.begin();
      for(; it != points.end(); ++it){
        //find triangle in the point location that contains point

        //create new edges on the DCEL

        //legalize wrong edges
      }
      return;
    }
};

int main() {
  std::vector<int> myInts;
  for(int i = 0; i < 10; i++){
    myInts.push_back(2 * i);
  }

  

  for(std::vector<int>::iterator it = myInts.begin(); it != myInts.end(); ++it){
    std::cout << *it << std::endl;
  }
  
  
  return 0;
}
