/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_QUADRATURE_HEXAHEDRON_HPP
#define VIENNAFEM_QUADRATURE_HEXAHEDRON_HPP

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
  // Exact for polynomials up to order 1
  //
  //
  template <typename InterfaceType>
  class rt_gauss_quad_element <viennagrid::hexahedron_tag, 1, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
    public:
      explicit rt_gauss_quad_element() : p_(3)
      {
        p_[0] = 1.0/2.0;
        p_[1] = 1.0/2.0;
        p_[2] = 1.0/2.0;
      }
      
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
  // Exact for polynomials up to degree 2
  //
  //
  

  //  
  //
  // Exact for polynomials up to degree 3:
  //
  //
  
  template <typename InterfaceType>
  class rt_gauss_quad_element <viennagrid::hexahedron_tag, 3, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
    public:
      explicit rt_gauss_quad_element() : abscissas_(8, std::vector<numeric_type>(2))
      {
        abscissas_[0][0] = 0.7886751345948125; abscissas_[0][1] = 0.7886751345948125; abscissas_[0][2] = 0.7886751345948125;
        abscissas_[1][0] = 0.7886751345948125; abscissas_[1][1] = 0.7886751345948125; abscissas_[1][2] = 0.2113248654051875;
        abscissas_[2][0] = 0.7886751345948125; abscissas_[2][1] = 0.2113248654051875; abscissas_[2][2] = 0.7886751345948125;
        abscissas_[3][0] = 0.7886751345948125; abscissas_[3][1] = 0.2113248654051875; abscissas_[3][2] = 0.2113248654051875;
        abscissas_[4][0] = 0.2113248654051875; abscissas_[4][1] = 0.7886751345948125; abscissas_[4][2] = 0.7886751345948125;
        abscissas_[5][0] = 0.2113248654051875; abscissas_[5][1] = 0.7886751345948125; abscissas_[5][2] = 0.2113248654051875;
        abscissas_[6][0] = 0.2113248654051875; abscissas_[6][1] = 0.2113248654051875; abscissas_[6][2] = 0.7886751345948125;
        abscissas_[7][0] = 0.2113248654051875; abscissas_[7][1] = 0.2113248654051875; abscissas_[7][2] = 0.2113248654051875;
      }
      
      NumericT eval(viennamath::rt_interval<InterfaceType> const & interv,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & var) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<8; ++i)
          result += viennamath::eval(e, abscissas_[i]);
        return 0.25 * result;
      }
      
    private:
      std::vector<std::vector<NumericT> > abscissas_;
  };
  
  
  

}
#endif
