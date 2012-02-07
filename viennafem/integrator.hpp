/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_INTEGRATOR_HPP
#define VIENNAFEM_INTEGRATOR_HPP

#include "viennafem/forwards.h"
#include "viennafem/cell_quan.hpp"

#include "viennamath/forwards.h"
#include "viennamath/expression.hpp"
#include "viennamath/manipulation/substitute.hpp"
#include "viennamath/manipulation/diff.hpp"
#include "viennamath/runtime/numerical_quadrature.hpp"

#include "viennagrid/topology/triangle.hpp"
#include "viennagrid/topology/tetrahedron.hpp"

namespace viennafem
{

  template <typename ElementType, typename InterfaceType = viennamath::default_interface_type>
  struct rt_gauss_quad_1_element;
  
  template <typename InterfaceType>
  class rt_gauss_quad_1_element <viennagrid::triangle_tag, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
    public:
      explicit rt_gauss_quad_1_element() : p_(2)
      {
        p_[0] = 1.0/3.0;
        p_[1] = 1.0/3.0;
      }
      
      viennamath::rt_expr<InterfaceType> eval(viennamath::rt_interval<InterfaceType> const & interv,
                                              viennamath::rt_expr<InterfaceType> const & e,
                                              viennamath::rt_variable<InterfaceType> const & var) const
      {
        return viennamath::rt_expr<InterfaceType>(new viennamath::rt_constant<NumericT, InterfaceType>(0.5 * e.get()->eval(p_)));
      }
      
    private:
      std::vector<numeric_type> p_;
  };
  
  template <typename InterfaceType>
  class rt_gauss_quad_1_element <viennagrid::tetrahedron_tag, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
    public:
      explicit rt_gauss_quad_1_element() : p_(3)
      {
        p_[0] = 1.0/4.0;
        p_[1] = 1.0/4.0;
        p_[2] = 1.0/4.0;
      }
      
      viennamath::rt_expr<InterfaceType> eval(viennamath::rt_interval<InterfaceType> const & interv,
                                              viennamath::rt_expr<InterfaceType> const & e,
                                              viennamath::rt_variable<InterfaceType> const & var) const
      {
        return viennamath::rt_expr<InterfaceType>(new viennamath::rt_constant<NumericT, InterfaceType>(e.get()->eval(p_) / 6.0));
      }
      
    private:
      std::vector<numeric_type> p_;
  };
  

}
#endif
