# Data Structure and Algorithms (Master)

Code for Data Structure and Algorithms discipline from Mathematical Modeling Master's course. Projects from the discipline are separated in folders.

## Delaunay Triangulation

Implementation of the Delaunay algorithm for the triangulation of points, i.e., with a set of points, indentify the optimal triangulation
that can be constructed with them. The implementation of the algorithm was done with C++ and the library CGAL, to demonstrate the results,
a application was made to triangulate an image and plot each triangle with the mean color of the vertex, this application was made with 
Python in Jupyter Notebooks and plotted in a webpage using ReGL.

### Getting started

#### Prerequisites

- CGAL version 5.3

> sudo apt-get install libcgal-dev

- Python opencv2 

#### Installation

- Clone the repository.

>

### Usage

1. Put the desired image in the folder _delaunay-triangulation/images/\<image-name\>.\<extension\>_ .
2. Run the functions from _processing_image.ipynb_, you can select if it will pick random points from the full image or only from the edges and the number of points.
3. Inside the folder delaunay-triangulation, run:

    > cgal_create_CMakeLists main
    >
    > cmake -DCMAKE_BUILD_TYPE=Release .
    > 
    > make 
    >
    > ./main images/\<image-name\>.txt

3. The result will be a json save in _delaunay-triangulation/images/\<image-name\>.json_, it will contain each triangle saved with coordinates and colors.
4. To visualize the result in the webpage, you need to update the links to the image in the script inside _docs/scripts.js_. 


