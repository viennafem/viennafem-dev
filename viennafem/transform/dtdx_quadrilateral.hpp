/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_TRANSFORM_QUADRILATERAL_HPP
#define VIENNAFEM_TRANSFORM_QUADRILATERAL_HPP

#include <iostream>
#include <utility>
#include "viennagrid/topology/triangle.hpp"
#include "viennagrid/algorithm/spanned_volume.hpp"
#include "viennagrid/domain.hpp"
#include "viennafem/forwards.h"


namespace viennafem
{

  template <>
  struct dt_dx_handler<viennafem::unit_quadrilateral>
  {
    public:
      
      template <typename CellType>
      static void apply(CellType const & cell)
      {
        typedef typename CellType::config_type       Config;
        typedef typename viennagrid::result_of::point<Config>::type   PointType;
        
        PointType const & p0 = viennagrid::ncells<0>(cell)[0].point();
        PointType const & p1 = viennagrid::ncells<0>(cell)[1].point();
        PointType const & p2 = viennagrid::ncells<0>(cell)[2].point();
        PointType const & p3 = viennagrid::ncells<0>(cell)[3].point();
        
        //Step 1: store determinant:
        numeric_type x1_x0 = p1[0] - p0[0];
        numeric_type x2_x0 = p2[0] - p0[0];
        numeric_type y1_y0 = p1[1] - p0[1];
        numeric_type y2_y0 = p2[1] - p0[1];
        numeric_type coeff_x = p0[0] - p1[0] - p2[0] + p3[0];
        numeric_type coeff_y = p0[1] - p1[1] - p2[1] + p3[1];
        
        std::vector<numeric_type> coeff_J(3);
        coeff_J[0] =   x1_x0 * y2_y0   -   y1_y0 * x2_x0;  //constant coefficent
        coeff_J[1] =   x1_x0 * coeff_y -   y1_y0 * coeff_x;
        coeff_J[2] = coeff_x * y2_y0   - coeff_y * x2_x0;
        
        viennadata::access<det_dF_dt_key, std::vector<numeric_type> >()(cell) = coeff_J;
        
        //Step 2: store partial derivatives:
        typedef std::pair<numeric_type, numeric_type>   Pair;
        
        viennadata::access<dt_dx_key<0, 0>, Pair>()(cell) = std::make_pair(y2_y0, coeff_y);
        viennadata::access<dt_dx_key<0, 1>, Pair>()(cell) = std::make_pair(-x2_x0, -coeff_x);
        viennadata::access<dt_dx_key<1, 0>, Pair>()(cell) = std::make_pair(-y1_y0, -coeff_y);
        viennadata::access<dt_dx_key<1, 1>, Pair>()(cell) = std::make_pair(x1_x0, coeff_x);
        
      }

  };

  
} //namespace

#endif