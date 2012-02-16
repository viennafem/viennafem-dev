//
// A simple hexahedral mesh generator for the unit cube
//


#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <cstdlib>

struct point
{
  double x;
  double y;
  double z;
};

class index_translator
{
  public:
    index_translator(std::size_t num) : points_per_coord_(num) {}

    std::size_t operator()(long i, long j, long k) const
    {
      return   i * points_per_coord_ * points_per_coord_
             + j * points_per_coord_
             + k
             + 1; //Netgen offset
    }
  private:
    std::size_t points_per_coord_;
};

int main()
{
  std::size_t points_per_coord = 8;
  std::vector< std::vector< std::vector<point> > > points(8, std::vector< std::vector<point> >(8, std::vector<point>(8)));

  //
  // setup points
  //
  for (std::size_t i=0; i<points_per_coord; ++i)
  {
    for (std::size_t j=0; j<points_per_coord; ++j)
    {
      for (std::size_t k=0; k<points_per_coord; ++k)
      {
        point p;
        p.x = i * 1.0 / (points_per_coord - 1.0);
        p.y = j * 1.0 / (points_per_coord - 1.0);
        p.z = k * 1.0 / (points_per_coord - 1.0);

        //add some jitter:
        if (p.x != 0.0 && p.x != 1.0)
          p.x += (rand() % 9 - 4.0) * 0.1 / (points_per_coord - 1.0);
        if (p.y != 0.0 && p.y != 1.0)
          p.y += (rand() % 9 - 4.0) * 0.1 / (points_per_coord - 1.0);
        if (p.z != 0.0 && p.z != 1.0)
          p.z += (rand() % 9 - 4.0) * 0.1 / (points_per_coord - 1.0);


        points[i][j][k] = p;
      }
    }
  }

  //
  // Write to file
  //

  // print points
  std::stringstream ss;
  ss << "cube" << (points_per_coord - 1) * (points_per_coord - 1) * (points_per_coord - 1) << "_hex.mesh";
  std::ofstream writer(ss.str().c_str());
  writer << (points_per_coord * points_per_coord * points_per_coord) << std::endl;
  
  for (std::size_t i=0; i<points_per_coord; ++i)
  {
    for (std::size_t j=0; j<points_per_coord; ++j)
    {
      for (std::size_t k=0; k<points_per_coord; ++k)
      {
        writer << points[i][j][k].x << " " 
               << points[i][j][k].y << " " 
               << points[i][j][k].z << std::endl;
      } 
    }
  }


  // print cells:
  writer << (points_per_coord-1) * (points_per_coord-1) * (points_per_coord-1) << std::endl;

  index_translator ijk_to_index(points_per_coord);
  for (std::size_t i=1; i<points_per_coord; ++i)
  {
    for (std::size_t j=1; j<points_per_coord; ++j)
    {
      for (std::size_t k=1; k<points_per_coord; ++k)
      {
         writer << "1 " 
                << ijk_to_index(i-1, j-1, k-1) << " "
                << ijk_to_index(i,   j-1, k-1) << " "
                << ijk_to_index(i-1, j, k-1) << " "
                << ijk_to_index(i,   j, k-1) << " "

                << ijk_to_index(i-1, j-1, k) << " "
                << ijk_to_index(i,   j-1, k) << " "
                << ijk_to_index(i-1, j, k) << " "
                << ijk_to_index(i,   j, k) << std::endl;
      } 
    }
  }

}
