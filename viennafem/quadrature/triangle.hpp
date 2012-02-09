/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_QUADRATURE_TRIANGLE_HPP
#define VIENNAFEM_QUADRATURE_TRIANGLE_HPP

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
  class rt_gauss_quad_element <viennagrid::triangle_tag, 1, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
    public:
      explicit rt_gauss_quad_element() : p_(2)
      {
        p_[0] = 1.0/3.0;
        p_[1] = 1.0/3.0;
      }
      
      NumericT eval(viennamath::rt_interval<InterfaceType> const & interv,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & var) const
      {
        return 0.5 * viennamath::eval(e, p_);
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
  // Exact for polynomials up to degree 7:
  //
  //
  template <typename InterfaceType>
  class rt_gauss_quad_element <viennagrid::triangle_tag, 7, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
    public:
      explicit rt_gauss_quad_element() : abscissas_(16, std::vector<numeric_type>(2)), weights_(16)
      {
        abscissas_[0][0] = 0.0571041961; abscissas_[0][1] = 0.06546699455602246;
        abscissas_[1][0] = 0.2768430136; abscissas_[1][1] = 0.05021012321401679;
        abscissas_[2][0] = 0.5835904324; abscissas_[2][1] = 0.02891208422223085;
        abscissas_[3][0] = 0.8602401357; abscissas_[3][1] = 0.009703785123906346;
        abscissas_[4][0] = 0.0571041961; abscissas_[4][1] = 0.3111645522491480;
        abscissas_[5][0] = 0.2768430136; abscissas_[5][1] = 0.2386486597440242;
        abscissas_[6][0] = 0.5835904324; abscissas_[6][1] = 0.1374191041243166;
        abscissas_[7][0] = 0.8602401357; abscissas_[7][1] = 0.04612207989200404;
        abscissas_[8][0] = 0.0571041961; abscissas_[8][1] = 0.6317312516508520;
        abscissas_[9][0] = 0.2768430136; abscissas_[9][1] = 0.4845083266559759;
        abscissas_[10][0] = 0.5835904324; abscissas_[10][1] = 0.2789904634756834;
        abscissas_[11][0] = 0.8602401357; abscissas_[11][1] = 0.09363778440799593;
        abscissas_[12][0] = 0.0571041961; abscissas_[12][1] = 0.8774288093439775;
        abscissas_[13][0] = 0.2768430136; abscissas_[13][1] = 0.6729468631859832;
        abscissas_[14][0] = 0.5835904324; abscissas_[14][1] = 0.3874974833777692;
        abscissas_[15][0] = 0.8602401357; abscissas_[15][1] = 0.1300560791760936;
        
        weights_[0] = 0.04713673637581137;
        weights_[1] = 0.07077613579259895;
        weights_[2] = 0.04516809856187617;
        weights_[3] = 0.01084645180365496;
        weights_[4] = 0.08837017702418863;
        weights_[5] = 0.1326884322074010;
        weights_[6] = 0.08467944903812383;
        weights_[7] = 0.02033451909634504;
        weights_[8] = 0.08837017702418863;
        weights_[9] = 0.1326884322074010;
        weights_[10] = 0.08467944903812383;
        weights_[11] = 0.02033451909634504;
        weights_[12] = 0.04713673637581137;
        weights_[13] = 0.07077613579259895;
        weights_[14] = 0.04516809856187617;
        weights_[15] = 0.01084645180365496;
        
      }
      
      NumericT eval(viennamath::rt_interval<InterfaceType> const & interv,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & var) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<16; ++i)
          result += weights_[i] * viennamath::eval(e, abscissas_[i]);
        return 0.5 * result;
      }
      
    private:
      std::vector<std::vector<NumericT> > abscissas_;
      std::vector<NumericT> weights_;
  };
  

}
#endif
