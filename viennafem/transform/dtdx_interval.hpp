#ifndef VIENNAFEM_TRANSFORM_INTERVAL_HPP
#define VIENNAFEM_TRANSFORM_INTERVAL_HPP

/* =========================================================================
   Copyright (c) 2012, Institute for Microelectronics,
                       Institute for Analysis and Scientific Computing,
                       TU Wien.
                             -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                             -----------------

   Author:     Karl Rupp                          rupp@iue.tuwien.ac.at

   License:    MIT (X11), see file LICENSE in the ViennaFEM base directory
============================================================================ */

#include <iostream>
#include <utility>
#include "viennagrid/topology/triangle.hpp"
#include "viennagrid/algorithm/spanned_volume.hpp"
#include "viennagrid/domain.hpp"
#include "viennafem/forwards.h"

/** @file    dtdx_interval.hpp
    @brief   Provides the transformation coefficients of an arbitrary interval to the unit interval
*/

namespace viennafem
{

  template <>
  struct dt_dx_handler<viennafem::unit_interval>
  {
    public:
      
      template <typename CellType>
      static void apply(CellType const & cell)
      {
        typedef typename CellType::config_type       Config;
        typedef typename viennagrid::result_of::point<Config>::type   PointType;
        
        PointType const & p0 = viennagrid::ncells<0>(cell)[0].point();
        PointType const & p1 = viennagrid::ncells<0>(cell)[1].point();
        
        //Step 1: store determinant:
        numeric_type x1_x0 = p1[0] - p0[0];
        
        viennadata::access<det_dF_dt_key, numeric_type>()(cell) = x1_x0;
        viennadata::access<dt_dx_key<0, 0>, numeric_type>()(cell) = 1.0 / x1_x0;
      }

  };

  
} //namespace

#endif