#ifndef VIENNAFEM_BASES_LINE_HPP
#define VIENNAFEM_BASES_LINE_HPP

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

#include <vector>
#include "viennagrid/topology/triangle.hpp"
#include "viennagrid/topology/tetrahedron.hpp"
#include "viennagrid/forwards.h"
#include "viennamath/expression.hpp"
#include "viennafem/forwards.h"

/** @file   viennafem/bases/line.hpp
    @brief  Defines the various basis functions for lines.
*/

namespace viennafem
{
  
 
  //
  // Lagrange family on the unit interval
  //
  
  
  // Vertex basis:
  /** @brief Returns the left vertex basis function */
  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_interval,
                      0,   // vertex level
                      0>   // vertex at x=0
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    
    typedef viennamath::ct_binary_expr< viennamath::ct_constant<1>,
                                        viennamath::op_minus<viennafem::numeric_type>,
                                        viennamath::ct_variable<0> 
                                      >                 type;
                                 
    static expression_type get() { return expression_type(type()); }
  };

  /** @brief Returns the right vertex basis function */
  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_interval,
                      0,   //vertex level
                      1>   //vertex at x=1
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    
    typedef viennamath::ct_variable<0>           type;
    
    static expression_type get() { return expression_type(type()); }
  };

  
  //
  // quadratic
  //
  /** @brief Returns the quadratic basis function defined in the interior of the line */
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_interval,
                      1,   //edge level
                      0>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;
    
    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_interval,
                                  0,
                                  0>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_interval,
                                  0,
                                  1>::type       phi_1;
                         
    //x * (1-x)
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;
    
    static expression_type get() { return expression_type(type()); }
  };
  
  //
  // cubic
  //
  /** @brief Returns the first cubic basis function defined in the interior of the line */
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<3>,
                      unit_interval,
                      1,   //edge level
                      0>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;
    
    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<3>,
                                  unit_interval,
                                  0,
                                  0>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<3>,
                                  unit_interval,
                                  0,
                                  1>::type       phi_1;
                         
    //(1-x)^2 * x
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_0
                                      > phi_0_squared;

    typedef viennamath::ct_binary_expr<phi_0_squared,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;
                               
    static expression_type get() { return expression_type(type()); }
  };

  
  /** @brief Returns the second cubic basis function defined in the interior of the line */
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<3>,
                      unit_interval,
                      1,   //edge level
                      1>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;
    
    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<3>,
                                  unit_interval,
                                  0,
                                  0>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<3>,
                                  unit_interval,
                                  0,
                                  1>::type       phi_1;
    
    //(1-x) * x^2
    typedef viennamath::ct_binary_expr<phi_1,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > phi_1_squared;

    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1_squared
                                      > type;
    
    static expression_type get() { return expression_type(type()); }
  };
  
}

#endif
