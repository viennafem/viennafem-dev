#ifndef VIENNAFEM_BOUNDARY_HPP
#define VIENNAFEM_BOUNDARY_HPP

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

#include "viennafem/forwards.h"

#include "viennamath/forwards.h"
#include "viennamath/manipulation/substitute.hpp"
#include "viennamath/expression.hpp"
#include "viennadata/api.hpp"

/** @file  boundary.hpp
    @brief Provide convenience routines for setting boundary conditions
*/

namespace viennafem
{
  template <typename StorageType, typename VertexType>
  void set_dirichlet_boundary(StorageType& storage,
                              VertexType const & v,
                              numeric_type const & value,
                              std::size_t id = 0)
  {
    typedef viennafem::boundary_key      BoundaryKey;;

    //set flag:
    viennadata::access<BoundaryKey, bool>(storage, BoundaryKey(id), v) = true;

    //set data:
    viennadata::access<BoundaryKey, numeric_type>(storage, BoundaryKey(id), v) = value;
  }

  template <typename StorageType, typename VertexType, typename NumericT>
  void set_dirichlet_boundary(StorageType& storage,
                              VertexType const & v,
                              std::vector<NumericT> const & value,
                              std::size_t id = 0)
  {
    typedef viennafem::boundary_key      BoundaryKey;;

    //set flag:
    viennadata::access<BoundaryKey, bool>(storage, BoundaryKey(id), v) = true;

    //set data:
    viennadata::access<BoundaryKey, std::vector<NumericT> >(storage, BoundaryKey(id), v) = value;
  }

}
#endif
