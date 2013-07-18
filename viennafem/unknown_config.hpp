#ifndef VIENNAFEM_UNKNOWN_CONFIG_HPP
#define VIENNAFEM_UNKNOWN_CONFIG_HPP

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

/** @file   unknown_config.hpp
    @brief  Configures the finite element space for a particular unknown.
*/

namespace viennafem
{

  /** @brief A configuration class for a particular PDE. [SUBJECT TO CHANGE!]
   *
   * @tparam MatrixType         System matrix type
   * @tparam VectorType         Type of the load vector
   * @tparam BoundaryKeyType    The key type to use with ViennaData to query a Dirichlet boundary flag
   * @tparam MappingKeyType     The key type to use with ViennaData to access the unknown indices
   */
  template <typename MatrixType,
            typename VectorType,
            typename BoundaryKeyType = boundary_key,
            typename MappingKeyType = mapping_key>
  class unknown_config
  {
    public:
      typedef MatrixType         matrix_type;
      typedef VectorType         vector_type;
      typedef BoundaryKeyType    boundary_key_type;
      typedef MappingKeyType     mapping_key_type;

      unknown_config(MatrixType & matrix, VectorType & vector) : m(matrix),v(vector) {}

      BoundaryKeyType boundary_key() const { return BoundaryKeyType(0); }
      MappingKeyType mapping_key() const { return MappingKeyType(0); }

      MatrixType & system_matrix() { return m; }
      VectorType & load_vector() { return v; }

    private:
      MatrixType & m;
      VectorType & v;
  };

}
#endif
