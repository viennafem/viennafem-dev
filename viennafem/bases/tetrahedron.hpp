/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_BASES_TETRAHEDRON_HPP
#define VIENNAFEM_BASES_TETRAHEDRON_HPP

#include <vector>
#include "viennagrid/topology/triangle.hpp"
#include "viennagrid/topology/tetrahedron.hpp"
#include "viennagrid/forwards.h"
#include "viennamath/expression.hpp"

namespace viennafem
{
  
  //TODO: Has to be generalized for:
  //       - different cell types
  //       - different degrees of basis functions
  //
  // For now, returns the linear basis:
  template <typename ExpressionType>
  std::vector<ExpressionType> get_basisfunctions(viennagrid::tetrahedron_tag)
  {
    typedef typename ExpressionType::interface_type    InterfaceType;
    
    //std::cout << "get_basisfunctions: entry" << std::endl;
    std::vector<ExpressionType> ret(4);
    
    //std::cout << "get_basisfunctions: Creating variables" << std::endl;
    viennamath::rt_variable<InterfaceType> x(0);
    viennamath::rt_variable<InterfaceType> y(1);
    viennamath::rt_variable<InterfaceType> z(2);
    
    //std::cout << "get_basisfunctions: filling 0" << std::endl;
    ret[0] = ExpressionType(1.0 - x - y - z);
    //std::cout << "get_basisfunctions: filling 1" << std::endl;
    ret[1] = ExpressionType(x);
    //std::cout << "get_basisfunctions: filling 2" << std::endl;
    ret[2] = ExpressionType(y);
    //std::cout << "get_basisfunctions: filling 3" << std::endl;
    ret[3] = ExpressionType(z);
    
    //std::cout << "get_basisfunctions: return" << std::endl;
    return ret;    
  }
  
  
  

}

#endif
