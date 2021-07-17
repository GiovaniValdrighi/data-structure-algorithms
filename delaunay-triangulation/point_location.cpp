#include <iostream>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/HalfedgeDS_vertex_base.h>
#include <CGAL/HalfedgeDS_decorator.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2; 
typedef Kernel::Triangle_2 Triangle_2;
struct Traits { 
  typedef Kernel::Point_2 Point_2; 
};
typedef CGAL::HalfedgeDS_default<Traits> HDS;
typedef CGAL::HalfedgeDS_vertex_base<Traits> Vertex;
typedef CGAL::HalfedgeDS_decorator<HDS> Decorator;

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

int main() {
  PointLocation D;
  Point_2 p(0.0, 0.0), q(1.1, 1.1), r(1.0, 0);
  Triangle_2 test_triangle(p, q, r);
  std::vector<TriangleNode*> triangles_vector;
  triangles_vector.push_back(new TriangleNode(test_triangle));
  D.insert(triangles_vector);
  std::cout << D.root->childs.size() << std::endl;
  Point_2 m(0.5, 0.2);
  TriangleNode* result = D.search(D.root, m);
  std::cout << result->childs.size() << std::endl;
  
  return 0;
}