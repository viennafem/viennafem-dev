#ifndef VIENNAFEM_QUADRATURE_HEXAHEDRON_HPP
#define VIENNAFEM_QUADRATURE_HEXAHEDRON_HPP

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

#include "viennagrid/topology/triangle.hpp"
#include "viennagrid/topology/tetrahedron.hpp"

/** @file   viennafem/quadrature/hexahedron.hpp
    @brief  Provides quadrature rules for hexahedra
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
  class rt_gauss_quad_element <viennafem::unit_cube, 1, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_gauss_quad_element <viennafem::unit_cube, 1, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      explicit rt_gauss_quad_element() : p_(3)
      {
        p_[0] = 1.0/2.0;
        p_[1] = 1.0/2.0;
        p_[2] = 1.0/2.0;
      }
      
      BaseType * clone() const { return new self_type(); }
      
      NumericT eval(viennamath::rt_interval<InterfaceType> const & /*interv*/,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & /*var*/) const
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
  
  /** @brief Gaussian quadrature rule exact for polynomials up to order 3 */
  template <typename InterfaceType>
  class rt_gauss_quad_element <viennafem::unit_cube, 3, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_gauss_quad_element <viennafem::unit_cube, 3, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      enum { num_points = 8 };
      
      explicit rt_gauss_quad_element() : abscissas_(num_points, std::vector<numeric_type>(3))
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
      
      BaseType * clone() const { return new self_type(); }
      
      NumericT eval(viennamath::rt_interval<InterfaceType> const & /*interv*/,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & /*var*/) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<num_points; ++i)
          result += viennamath::eval(e, abscissas_[i]);
        return result / 8.0;
      }
      
    private:
      std::vector<std::vector<NumericT> > abscissas_;
  };
  
  
  //
  //
  // Exact for polynomials up to order 5
  //
  //
  /** @brief Gaussian quadrature rule exact for polynomials up to order 5 */
  template <typename InterfaceType>
  class rt_gauss_quad_element <viennafem::unit_cube, 5, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_gauss_quad_element <viennafem::unit_cube, 5, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      enum { num_points = 27 };
      
      explicit rt_gauss_quad_element() : abscissas_(num_points, std::vector<numeric_type>(3)), weights_(num_points)
      {
        abscissas_[0][0] = 0.11270166537925829786; abscissas_[0][1] = 0.11270166537925829786; abscissas_[0][2] = 0.11270166537925829786;
        abscissas_[1][0] = 0.11270166537925829786; abscissas_[1][1] = 0.11270166537925829786; abscissas_[1][2] = 0.5;
        abscissas_[2][0] = 0.11270166537925829786; abscissas_[2][1] = 0.11270166537925829786; abscissas_[2][2] = 0.88729833462074170214;
        abscissas_[3][0] = 0.11270166537925829786; abscissas_[3][1] = 0.5;                    abscissas_[3][2] = 0.11270166537925829786;
        abscissas_[4][0] = 0.11270166537925829786; abscissas_[4][1] = 0.5;                    abscissas_[4][2] = 0.5;
        abscissas_[5][0] = 0.11270166537925829786; abscissas_[5][1] = 0.5;                    abscissas_[5][2] = 0.88729833462074170214;
        abscissas_[6][0] = 0.11270166537925829786; abscissas_[6][1] = 0.88729833462074170214; abscissas_[6][2] = 0.11270166537925829786;
        abscissas_[7][0] = 0.11270166537925829786; abscissas_[7][1] = 0.88729833462074170214; abscissas_[7][2] = 0.5;
        abscissas_[8][0] = 0.11270166537925829786; abscissas_[8][1] = 0.88729833462074170214; abscissas_[8][2] = 0.88729833462074170214;

        abscissas_[ 9][0] = 0.5;                   abscissas_[ 9][1] = 0.11270166537925829786; abscissas_[ 9][2] = 0.11270166537925829786;
        abscissas_[10][0] = 0.5;                   abscissas_[10][1] = 0.11270166537925829786; abscissas_[10][2] = 0.5;
        abscissas_[11][0] = 0.5;                   abscissas_[11][1] = 0.11270166537925829786; abscissas_[11][2] = 0.88729833462074170214;
        abscissas_[12][0] = 0.5;                   abscissas_[12][1] = 0.5;                    abscissas_[12][2] = 0.11270166537925829786;
        abscissas_[13][0] = 0.5;                   abscissas_[13][1] = 0.5;                    abscissas_[13][2] = 0.5;
        abscissas_[14][0] = 0.5;                   abscissas_[14][1] = 0.5;                    abscissas_[14][2] = 0.88729833462074170214;
        abscissas_[15][0] = 0.5;                   abscissas_[15][1] = 0.88729833462074170214; abscissas_[15][2] = 0.11270166537925829786;
        abscissas_[16][0] = 0.5;                   abscissas_[16][1] = 0.88729833462074170214; abscissas_[16][2] = 0.5;
        abscissas_[17][0] = 0.5;                   abscissas_[17][1] = 0.88729833462074170214; abscissas_[17][2] = 0.88729833462074170214;

        abscissas_[18][0] = 0.88729833462074170214; abscissas_[18][1] = 0.11270166537925829786; abscissas_[18][2] = 0.11270166537925829786;
        abscissas_[19][0] = 0.88729833462074170214; abscissas_[19][1] = 0.11270166537925829786; abscissas_[19][2] = 0.5;
        abscissas_[20][0] = 0.88729833462074170214; abscissas_[20][1] = 0.11270166537925829786; abscissas_[20][2] = 0.88729833462074170214;
        abscissas_[21][0] = 0.88729833462074170214; abscissas_[21][1] = 0.5;                    abscissas_[21][2] = 0.11270166537925829786;
        abscissas_[22][0] = 0.88729833462074170214; abscissas_[22][1] = 0.5;                    abscissas_[22][2] = 0.5;
        abscissas_[23][0] = 0.88729833462074170214; abscissas_[23][1] = 0.5;                    abscissas_[23][2] = 0.88729833462074170214;
        abscissas_[24][0] = 0.88729833462074170214; abscissas_[24][1] = 0.88729833462074170214; abscissas_[24][2] = 0.11270166537925829786;
        abscissas_[25][0] = 0.88729833462074170214; abscissas_[25][1] = 0.88729833462074170214; abscissas_[25][2] = 0.5;
        abscissas_[26][0] = 0.88729833462074170214; abscissas_[26][1] = 0.88729833462074170214; abscissas_[26][2] = 0.88729833462074170214;
        
        // weights:
        double denominator = 18.0 * 18.0 * 18.0;
        weights_[0] = (5.0 * 5.0 * 5.0) / denominator;
        weights_[1] = (5.0 * 5.0 * 8.0) / denominator;
        weights_[2] = (5.0 * 5.0 * 5.0) / denominator;
        weights_[3] = (5.0 * 8.0 * 5.0) / denominator;
        weights_[4] = (5.0 * 8.0 * 8.0) / denominator;
        weights_[5] = (5.0 * 8.0 * 5.0) / denominator;
        weights_[6] = (5.0 * 5.0 * 5.0) / denominator;
        weights_[7] = (5.0 * 5.0 * 8.0) / denominator;
        weights_[8] = (5.0 * 5.0 * 5.0) / denominator;
        
        weights_[ 9] = (8.0 * 5.0 * 5.0) / denominator;
        weights_[10] = (8.0 * 5.0 * 8.0) / denominator;
        weights_[11] = (8.0 * 5.0 * 5.0) / denominator;
        weights_[12] = (8.0 * 8.0 * 5.0) / denominator;
        weights_[13] = (8.0 * 8.0 * 8.0) / denominator;
        weights_[14] = (8.0 * 8.0 * 5.0) / denominator;
        weights_[15] = (8.0 * 5.0 * 5.0) / denominator;
        weights_[16] = (8.0 * 5.0 * 8.0) / denominator;
        weights_[17] = (8.0 * 5.0 * 5.0) / denominator;

        weights_[18] = (5.0 * 5.0 * 5.0) / denominator;
        weights_[19] = (5.0 * 5.0 * 8.0) / denominator;
        weights_[20] = (5.0 * 5.0 * 5.0) / denominator;
        weights_[21] = (5.0 * 8.0 * 5.0) / denominator;
        weights_[22] = (5.0 * 8.0 * 8.0) / denominator;
        weights_[23] = (5.0 * 8.0 * 5.0) / denominator;
        weights_[24] = (5.0 * 5.0 * 5.0) / denominator;
        weights_[25] = (5.0 * 5.0 * 8.0) / denominator;
        weights_[26] = (5.0 * 5.0 * 5.0) / denominator;        
      }
      
      BaseType * clone() const { return new self_type(); }
      
      NumericT eval(viennamath::rt_interval<InterfaceType> const & /*interv*/,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & /*var*/) const
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
  

}
#endif
