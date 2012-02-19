/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_BASES_HEXAHEDRON_HPP
#define VIENNAFEM_BASES_HEXAHEDRON_HPP

#include <vector>
#include "viennagrid/topology/triangle.hpp"
#include "viennagrid/topology/tetrahedron.hpp"
#include "viennagrid/forwards.h"
#include "viennamath/expression.hpp"
#include "viennafem/forwards.h"
#include "viennafem/bases/quadrilateral.hpp"

namespace viennafem
{
  
 
  //
  // Lagrange family on hexahedra
  //
  
  
  // Vertex basis:
  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_hexahedron,
                      0,   // vertex level
                      0>   // vertex (0,0,0)
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    
    typedef viennamath::ct_binary_expr< typename local_basis<InterfaceType, viennafem::lagrange_tag<order>,
                                                             unit_quadrilateral, 0, 0>::type,
                                        viennamath::op_mult<viennafem::numeric_type>,
                                        viennamath::ct_binary_expr< viennamath::ct_constant<1>,
                                                                    viennamath::op_minus<viennafem::numeric_type>,
                                                                    viennamath::ct_variable<2> >
                                      >       type;                       
    
    static expression_type get() { return expression_type(type()); }
  };

  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_hexahedron,
                      0,   //vertex level
                      1>   //vertex (1,0,0)
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    
    typedef viennamath::ct_binary_expr< typename local_basis<InterfaceType, viennafem::lagrange_tag<order>,
                                                              unit_quadrilateral, 0, 1>::type,
                                        viennamath::op_mult<viennafem::numeric_type>,
                                        viennamath::ct_binary_expr< viennamath::ct_constant<1>,
                                                              viennamath::op_minus<viennafem::numeric_type>,
                                                              viennamath::ct_variable<2> >
                                      >       type;                       
    
    static expression_type get() { return expression_type(type()); }
  };

  
  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_hexahedron,
                      0,   //vertex level
                      2>   //vertex (0,1,0)
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    
    typedef viennamath::ct_binary_expr< typename local_basis<InterfaceType, viennafem::lagrange_tag<order>,
                                                              unit_quadrilateral, 0, 2>::type,
                                        viennamath::op_mult<viennafem::numeric_type>,
                                        viennamath::ct_binary_expr< viennamath::ct_constant<1>,
                                                              viennamath::op_minus<viennafem::numeric_type>,
                                                              viennamath::ct_variable<2> >
                                      >       type;                       
    
    static expression_type get() { return expression_type(type()); }
  };
  
  
  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_hexahedron,
                      0,   //vertex level
                      3>   //vertex (1,1,0)
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    
    typedef viennamath::ct_binary_expr< typename local_basis<InterfaceType, viennafem::lagrange_tag<order>,
                                                              unit_quadrilateral, 0, 3>::type,
                                        viennamath::op_mult<viennafem::numeric_type>,
                                        viennamath::ct_binary_expr< viennamath::ct_constant<1>,
                                                              viennamath::op_minus<viennafem::numeric_type>,
                                                              viennamath::ct_variable<2> >
                                      >       type;                       
    
    static expression_type get() { return expression_type(type()); }
  };
  

  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_hexahedron,
                      0,   // vertex level
                      4>   // vertex (0,0,1)
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    
    typedef viennamath::ct_binary_expr< typename local_basis<InterfaceType, viennafem::lagrange_tag<order>,
                                                              unit_quadrilateral, 0, 0>::type,
                                        viennamath::op_mult<viennafem::numeric_type>,
                                        viennamath::ct_variable<2>
                                      >       type;                       
    
    static expression_type get() { return expression_type(type()); }
  };

  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_hexahedron,
                      0,   //vertex level
                      5>   //vertex (1,0,1)
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    
    typedef viennamath::ct_binary_expr< typename local_basis<InterfaceType, viennafem::lagrange_tag<order>,
                                                              unit_quadrilateral, 0, 1>::type,
                                        viennamath::op_mult<viennafem::numeric_type>,
                                        viennamath::ct_variable<2>
                                      >       type;                       
    
    static expression_type get() { return expression_type(type()); }
  };

  
  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_hexahedron,
                      0,   //vertex level
                      6>   //vertex (0,1,1)
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    
    typedef viennamath::ct_binary_expr< typename local_basis<InterfaceType, viennafem::lagrange_tag<order>,
                                                              unit_quadrilateral, 0, 2>::type,
                                        viennamath::op_mult<viennafem::numeric_type>,
                                        viennamath::ct_variable<2>
                                      >       type;                       
    
    static expression_type get() { return expression_type(type()); }
  };
  
  
  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_hexahedron,
                      0,   //vertex level
                      7>   //vertex (1,1,1)
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    
    typedef viennamath::ct_binary_expr< typename local_basis<InterfaceType, viennafem::lagrange_tag<order>,
                                                              unit_quadrilateral, 0, 3>::type,
                                        viennamath::op_mult<viennafem::numeric_type>,
                                        viennamath::ct_variable<2>
                                      >       type;                       
    
    static expression_type get() { return expression_type(type()); }
  };
  
  
  
  //
  // quadratic:
  //
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_hexahedron,
                      1,   //edge level
                      0>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;
    
    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  0>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  1>::type       phi_1;
                         
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;
    
    static expression_type get() { return expression_type(type()); }
  };
  
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_hexahedron,
                      1,   //edge level
                      1>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;
    
    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  0>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  2>::type       phi_1;
                         
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;
    
    static expression_type get() { return expression_type(type()); }
  };
  
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_hexahedron,
                      1,   //edge level
                      2>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;
    
    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  0>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  4>::type       phi_1;
                         
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;
    
    static expression_type get() { return expression_type(type()); }
  };
  
  
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_hexahedron,
                      1,   //edge level
                      3>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;
    
    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  1>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  3>::type       phi_1;
                         
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;
    
    static expression_type get() { return expression_type(type()); }
  };
  
  
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_hexahedron,
                      1,   //edge level
                      4>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;
    
    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  1>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  5>::type       phi_1;
                         
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;
    
    static expression_type get() { return expression_type(type()); }
  };
  
  
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_hexahedron,
                      1,   //edge level
                      5>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;
    
    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  2>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  3>::type       phi_1;
                         
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;
    
    static expression_type get() { return expression_type(type()); }
  };
  
  
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_hexahedron,
                      1,   //edge level
                      6>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;
    
    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  2>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  6>::type       phi_1;
                         
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;
    
    static expression_type get() { return expression_type(type()); }
  };
  

  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_hexahedron,
                      1,   //edge level
                      7>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;
    
    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  3>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  7>::type       phi_1;
                         
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;
    
    static expression_type get() { return expression_type(type()); }
  };
  
  
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_hexahedron,
                      1,   //edge level
                      8>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;
    
    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  4>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  5>::type       phi_1;
                         
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;
    
    static expression_type get() { return expression_type(type()); }
  };
  

  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_hexahedron,
                      1,   //edge level
                      9>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;
    
    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  4>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  6>::type       phi_1;
                         
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;
    
    static expression_type get() { return expression_type(type()); }
  };
  
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_hexahedron,
                      1,   //edge level
                      10>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;
    
    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  5>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  7>::type       phi_1;
                         
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;
    
    static expression_type get() { return expression_type(type()); }
  };

  
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_hexahedron,
                      1,   //edge level
                      11>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;
    
    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  6>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_hexahedron,
                                  0,
                                  7>::type       phi_1;
                         
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                     > type;
    
    static expression_type get() { return expression_type(type()); }
  };
  
  
}

#endif
