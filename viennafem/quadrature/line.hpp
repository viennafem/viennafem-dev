#ifndef VIENNAFEM_QUADRATURE_LINE_HPP
#define VIENNAFEM_QUADRATURE_LINE_HPP

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
#include "viennafem/cell_quan.hpp"

#include "viennamath/forwards.h"
#include "viennamath/expression.hpp"
#include "viennamath/manipulation/substitute.hpp"
#include "viennamath/manipulation/diff.hpp"
#include "viennamath/manipulation/eval.hpp"
#include "viennamath/runtime/numerical_quadrature.hpp"

#include "viennagrid/topology/simplex.hpp"

/** @file   viennafem/quadrature/line.hpp
    @brief  Provides quadrature rules for lines
*/

namespace viennafem
{
  
  //
  //
  // Exact for polynomials up to order 1
  //
  //
  /** @brief Gaussian quadrature rule exact for polynomials up to order 1 */
  template <typename InterfaceType>
  class rt_gauss_quad_element <viennafem::unit_interval, 1, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_gauss_quad_element <viennafem::unit_interval, 1, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
      
    public:
      
      BaseType * clone() const { return new self_type(); }
      
      NumericT eval(viennamath::rt_interval<InterfaceType> const & /*interv*/,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & /*var*/) const
      {
        return viennamath::eval(e, 0.5);  //Note: reference line is [0,1] here.
      }
      
  };
  
  //
  //
  // Exact for polynomials up to degree 2
  //
  //

  

  //  
  //
  // Exact for polynomials up to degree 3:
  //
  //
  /** @brief Gaussian quadrature rule exact for polynomials up to order 3 */
  template <typename InterfaceType>
  class rt_gauss_quad_element <viennafem::unit_interval, 3, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_gauss_quad_element <viennafem::unit_interval, 3, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      
      BaseType * clone() const { return new self_type(); }
      
      NumericT eval(viennamath::rt_interval<InterfaceType> const & /*interv*/,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & /*var*/) const
      {
        return 0.5 * (  viennamath::eval(e, 0.7886751345948125) 
                      + viennamath::eval(e, 0.2113248654051875));  //Note: reference line is [0,1] here.
      }
      
  };
  
  
  //  
  //
  // Exact for polynomials up to degree 5:
  //
  //
  /** @brief Gaussian quadrature rule exact for polynomials up to order 5 */
  template <typename InterfaceType>
  class rt_gauss_quad_element <viennafem::unit_interval, 5, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_gauss_quad_element <viennafem::unit_interval, 5, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      
      BaseType * clone() const { return new self_type(); }
      
      NumericT eval(viennamath::rt_interval<InterfaceType> const & /*interv*/,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & /*var*/) const
      {
        return (  5.0 * viennamath::eval(e, 0.11270166537925829786) 
                + 8.0 * viennamath::eval(e, 0.5)
                + 5.0 * viennamath::eval(e, 0.88729833462074170214)) / 18.0;  //Note: reference line is [0,1] here.
      }
      
  };
  

}
#endif
