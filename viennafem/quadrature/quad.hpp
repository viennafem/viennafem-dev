/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_QUADRATURE_QUAD_HPP
#define VIENNAFEM_QUADRATURE_QUAD_HPP

#include "viennafem/forwards.h"

#include "viennafem/quadrature/line.hpp"

#include "viennafem/quadrature/triangle.hpp"
#include "viennafem/quadrature/quadrilateral.hpp"


#include "viennafem/quadrature/tetrahedron.hpp"
#include "viennafem/quadrature/hexahedron.hpp"

namespace viennafem
{
  // prototype:
  template <typename ReferenceCellType, typename BasisTag>
  viennamath::numerical_quadrature make_quadrature_rule_impl(BasisTag)
  {
    std::cerr << "Quadrature rule not available!" << std::endl;
    throw "Quadrature rule not available!";
    return viennamath::numerical_quadrature(NULL);
  }

  //
  // Line:
  //
  template <>
  viennamath::numerical_quadrature make_quadrature_rule_impl<viennafem::unit_interval>(lagrange_tag<1>)
  {
    // conservative estimate: second-order accurate:
    return viennamath::numerical_quadrature(new rt_gauss_quad_element <viennafem::unit_interval, 3>());
  }

  template <>
  viennamath::numerical_quadrature make_quadrature_rule_impl<viennafem::unit_interval>(lagrange_tag<2>)
  {
    // conservative estimate: forth-order accurate:
    return viennamath::numerical_quadrature(new rt_gauss_quad_element <viennafem::unit_interval, 5>());
  }

  template <>
  viennamath::numerical_quadrature make_quadrature_rule_impl<viennafem::unit_interval>(lagrange_tag<3>)
  {
    // use fifth-order rule and hope the best (no sixth-order rule available so far)
    return viennamath::numerical_quadrature(new rt_gauss_quad_element <viennafem::unit_interval, 5>());
  }

  //
  // Quadrilateral:
  //
  template <>
  viennamath::numerical_quadrature make_quadrature_rule_impl<viennafem::unit_quadrilateral>(lagrange_tag<1>)
  {
    // conservative estimate: second-order accurate:
    return viennamath::numerical_quadrature(new rt_gauss_quad_element <viennafem::unit_quadrilateral, 3>());
  }

  template <>
  viennamath::numerical_quadrature make_quadrature_rule_impl<viennafem::unit_quadrilateral>(lagrange_tag<2>)
  {
    // conservative estimate: forth-order accuracy required:
    return viennamath::numerical_quadrature(new rt_gauss_quad_element <viennafem::unit_quadrilateral, 5>());
  }

  //
  // Triangle:
  //
  template <>
  viennamath::numerical_quadrature make_quadrature_rule_impl<viennafem::unit_triangle>(lagrange_tag<1>)
  {
    // conservative estimate: second-order accurate:
    return viennamath::numerical_quadrature(new rt_strang_quad_element <viennafem::unit_triangle, 2>());
  }

  template <>
  viennamath::numerical_quadrature make_quadrature_rule_impl<viennafem::unit_triangle>(lagrange_tag<2>)
  {
    // conservative estimate: forth-order accuracy required:
    return viennamath::numerical_quadrature(new rt_strang_quad_element <viennafem::unit_triangle, 4>());
  }


  //
  // Hexahedron:
  //
  template <>
  viennamath::numerical_quadrature make_quadrature_rule_impl<viennafem::unit_hexahedron>(lagrange_tag<1>)
  {
    // conservative estimate: second-order accurate:
    return viennamath::numerical_quadrature(new rt_gauss_quad_element <viennafem::unit_hexahedron, 3>());
  }

  template <>
  viennamath::numerical_quadrature make_quadrature_rule_impl<viennafem::unit_hexahedron>(lagrange_tag<2>)
  {
    // conservative estimate: forth-order accuracy required:
    return viennamath::numerical_quadrature(new rt_gauss_quad_element <viennafem::unit_hexahedron, 5>());
  }


  //
  // Tetrahedron:
  //
  template <>
  viennamath::numerical_quadrature make_quadrature_rule_impl<viennafem::unit_tetrahedron>(lagrange_tag<1>)
  {
    // conservative estimate: second-order accurate:
    return viennamath::numerical_quadrature(new rt_keast_quad_element <viennafem::unit_tetrahedron, 2>());
  }

  template <>
  viennamath::numerical_quadrature make_quadrature_rule_impl<viennafem::unit_tetrahedron>(lagrange_tag<2>)
  {
    // conservative estimate: forth-order accuracy required:
    return viennamath::numerical_quadrature(new rt_keast_quad_element <viennafem::unit_tetrahedron, 4>());
  }



  //
  // Public interface:
  //
  
  template <typename PDESystemType, typename DomainType>
  viennamath::numerical_quadrature make_quadrature_rule(PDESystemType const & pde_system, DomainType const & domain)
  {
    typedef typename DomainType::config_type              Config;
    typedef typename Config::cell_tag                     CellTag;
    
    // runtime switcher:
    switch (pde_system.option(0).trial_space_id())
    {
      case space_to_id< lagrange_tag<1> >::value:
          return make_quadrature_rule_impl<typename reference_cell_for_basis<CellTag, viennafem::lagrange_tag<1> >::type>(lagrange_tag<1>());
      case space_to_id< lagrange_tag<2> >::value:
          return make_quadrature_rule_impl<typename reference_cell_for_basis<CellTag, viennafem::lagrange_tag<2> >::type>(lagrange_tag<2>());
      case space_to_id< lagrange_tag<3> >::value:
          return make_quadrature_rule_impl<typename reference_cell_for_basis<CellTag, viennafem::lagrange_tag<3> >::type>(lagrange_tag<3>());
    }
    
    std::cerr << "make_quadrature_rule(): Cannot derive quadrature rule!" << std::endl;
    throw "Cannot derive quadrature rule!";
    
    return viennamath::numerical_quadrature(NULL);
  }
  
  
}


#endif
