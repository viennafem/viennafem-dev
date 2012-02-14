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
#include "viennafem/forwards.h"

namespace viennafem
{
  
 
  //
  // Lagrange family on tetrahedra
  //
  
  
  // Vertex basis:
  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_tetrahedron,
                      0,   // vertex level
                      0>   // vertex (0,0,0)
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    
    typedef viennamath::ct_expr< viennamath::ct_expr< viennamath::ct_constant<1>,
                                                      viennamath::op_minus<viennafem::numeric_type>,
                                                      viennamath::ct_variable<0> >,
                                 viennamath::op_minus<viennafem::numeric_type>,
                                 viennamath::ct_expr< viennamath::ct_variable<1>,
                                                      viennamath::op_plus<viennafem::numeric_type>,
                                                      viennamath::ct_variable<2> >
                               >                 type;
                                 
    static std::vector<expression_type> get()
    {
      std::vector<expression_type> ret(1);
      ret[0] = expression_type(type());
      return ret;    
    }
  };

  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_tetrahedron,
                      0,   //vertex level
                      1>   //vertex (1,0,0)
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    
    typedef viennamath::ct_variable<0>           type;
    
    static std::vector<expression_type> get()
    {
      std::vector<expression_type> ret(1);
      ret[0] = expression_type(type());
      return ret;    
    }
  };

  
  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_tetrahedron,
                      0,   //vertex level
                      2>   //vertex (0,1,0)
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    
    typedef viennamath::ct_variable<1>           type;
    
    static std::vector<expression_type> get()
    {
      std::vector<expression_type> ret(1);
      ret[0] = expression_type(type());
      return ret;    
    }
  };
  
  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_tetrahedron,
                      0,   //vertex level
                      3>   //vertex (0,0,1)
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    
    typedef viennamath::ct_variable<2>           type;
    
    static std::vector<expression_type> get()
    {
      std::vector<expression_type> ret(1);
      ret[0] = expression_type(type());
      return ret;    
    }
  };
  
  
}

#endif
