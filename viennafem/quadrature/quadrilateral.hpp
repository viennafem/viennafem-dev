/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_QUADRATURE_QUADRILATERAL_HPP
#define VIENNAFEM_QUADRATURE_QUADRILATERAL_HPP

#include "viennafem/forwards.h"
#include "viennafem/cell_quan.hpp"

#include "viennamath/forwards.h"
#include "viennamath/expression.hpp"
#include "viennamath/manipulation/substitute.hpp"
#include "viennamath/manipulation/diff.hpp"
#include "viennamath/manipulation/eval.hpp"
#include "viennamath/runtime/numerical_quadrature.hpp"

#include "viennagrid/topology/triangle.hpp"
#include "viennagrid/topology/tetrahedron.hpp"

namespace viennafem
{
  //
  //
  //                 S E C T I O N    1
  //  
  //  Reference tetrahedron with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1)
  //
  //
  //
  
  
  
  
  //
  //
  // Exact for polynomials up to order 1
  //
  //
  template <typename InterfaceType>
  class rt_gauss_quad_element <viennafem::unit_quadrilateral, 1, InterfaceType> 
   : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_gauss_quad_element <viennafem::unit_quadrilateral, 1, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      enum { num_points = 1 };
      
      explicit rt_gauss_quad_element() : p_(2)
      {
        p_[0] = 0.5;
        p_[1] = 0.5;
      }
      
      BaseType * clone() const { return new self_type(); }
      
      NumericT eval(viennamath::rt_interval<InterfaceType> const & interv,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & var) const
      {
        return viennamath::eval(e, p_);
      }
      
    private:
      std::vector<NumericT> p_;
  };

  
  //
  //
  // Exact for polynomials up to order 3
  //
  //
  template <typename InterfaceType>
  class rt_gauss_quad_element <viennafem::unit_quadrilateral, 3, InterfaceType> 
   : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_gauss_quad_element <viennafem::unit_quadrilateral, 3, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      enum { num_points = 4 };
      
      explicit rt_gauss_quad_element() : abscissas_(num_points, std::vector<numeric_type>(2))
      {
        abscissas_[0][0] = 0.7886751345948125; abscissas_[0][1] = 0.7886751345948125;
        abscissas_[1][0] = 0.7886751345948125; abscissas_[1][1] = 0.2113248654051875;
        abscissas_[2][0] = 0.2113248654051875; abscissas_[2][1] = 0.7886751345948125;
        abscissas_[3][0] = 0.2113248654051875; abscissas_[3][1] = 0.2113248654051875;
      }
      
      BaseType * clone() const { return new self_type(); }
      
      NumericT eval(viennamath::rt_interval<InterfaceType> const & interv,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & var) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<num_points; ++i)
          result += viennamath::eval(e, abscissas_[i]);
        return 0.25 * result;
      }
      
    private:
      std::vector<std::vector<NumericT> > abscissas_;
  };
  

  //
  //
  // Exact for polynomials up to order 5
  //
  //
  template <typename InterfaceType>
  class rt_gauss_quad_element <viennafem::unit_quadrilateral, 5, InterfaceType> 
   : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_gauss_quad_element <viennafem::unit_quadrilateral, 5, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      enum { num_points = 9 };
      
      explicit rt_gauss_quad_element() : abscissas_(num_points, std::vector<numeric_type>(2)), weights_(num_points)
      {
        abscissas_[0][0] = 0.11270166537925829786; abscissas_[0][1] = 0.11270166537925829786;
        abscissas_[1][0] = 0.11270166537925829786; abscissas_[1][1] = 0.5;
        abscissas_[2][0] = 0.11270166537925829786; abscissas_[2][1] = 0.88729833462074170214;
        abscissas_[3][0] = 0.5;                    abscissas_[3][1] = 0.11270166537925829786;
        abscissas_[4][0] = 0.5;                    abscissas_[4][1] = 0.5;
        abscissas_[5][0] = 0.5;                    abscissas_[5][1] = 0.88729833462074170214;
        abscissas_[6][0] = 0.88729833462074170214; abscissas_[6][1] = 0.11270166537925829786;
        abscissas_[7][0] = 0.88729833462074170214; abscissas_[7][1] = 0.5;
        abscissas_[8][0] = 0.88729833462074170214; abscissas_[8][1] = 0.88729833462074170214;
        
        weights_[0] = (5.0 * 5.0) / (18.0 * 18.0);
        weights_[1] = (5.0 * 8.0) / (18.0 * 18.0);
        weights_[2] = (5.0 * 5.0) / (18.0 * 18.0);
        
        weights_[3] = (8.0 * 5.0) / (18.0 * 18.0);
        weights_[4] = (8.0 * 8.0) / (18.0 * 18.0);
        weights_[5] = (8.0 * 5.0) / (18.0 * 18.0);
        
        weights_[6] = (5.0 * 5.0) / (18.0 * 18.0);
        weights_[7] = (5.0 * 8.0) / (18.0 * 18.0);
        weights_[8] = (5.0 * 5.0) / (18.0 * 18.0);
      }
      
      BaseType * clone() const { return new self_type(); }
      
      NumericT eval(viennamath::rt_interval<InterfaceType> const & interv,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & var) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<num_points; ++i)
          result += weights_[i] * viennamath::eval(e, abscissas_[i]);
        return result;
      }
      
    private:
      std::vector<std::vector<NumericT> > abscissas_;
      std::vector<NumericT> weights_;
  };
  
  
  
  
  
  //
  // Quadrature rule generator
  //
  template <typename InterfaceType>
  viennamath::numerical_quadrature quadrature_for_reference_cell(viennafem::unit_quadrilateral const &,
                                                                 std::size_t order)
  {
    switch (order)
    {
      case 1: return viennamath::numerical_quadrature(new rt_gauss_quad_element <viennafem::unit_quadrilateral, 1, InterfaceType>());
      case 2:
      case 3: return viennamath::numerical_quadrature(new rt_gauss_quad_element <viennafem::unit_quadrilateral, 3, InterfaceType>());
      case 4:
      case 5: return viennamath::numerical_quadrature(new rt_gauss_quad_element <viennafem::unit_quadrilateral, 5, InterfaceType>());
    }
    
    std::cout << "Cannot find quadrature rule for order " << order << " - fallback to order 5." << std::endl;
    return viennamath::numerical_quadrature(new rt_gauss_quad_element <viennafem::unit_quadrilateral, 5, InterfaceType>());
  }
  
}
#endif
