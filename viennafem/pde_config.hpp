/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_PDE_CONFIG_HPP
#define VIENNAFEM_PDE_CONFIG_HPP

#include "viennafem/forwards.h"

namespace viennafem
{
  
  /** @brief A configuration class for a particular PDE.
   * 
   * @tparam BoundaryKeyType    The key type to use with ViennaData to query a Dirichlet boundary flag
   * @tparam MappingKeyType     The key type to use with ViennaData to access the unknown indices
   */
  template <typename BoundaryKeyType = boundary_key,
            typename MappingKeyType = mapping_key>
  struct pde_config
  {
    typedef BoundaryKeyType    boundary_key_type;
    typedef MappingKeyType     mapping_key_type;
    
    BoundaryKeyType boundary_key() const { return BoundaryKeyType(); } 
    MappingKeyType mapping_key() const { return MappingKeyType(); } 
  };  

}
#endif
