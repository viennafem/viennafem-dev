#ifndef VIENNAFEM_QUADRATURE_TETRAHEDRON_HPP
#define VIENNAFEM_QUADRATURE_TETRAHEDRON_HPP

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

/** @file   viennafem/quadrature/tetrahedron.hpp
    @brief  Provides quadrature rules for tetrahedra
*/

namespace viennafem
{

  //
  // Have a look at http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html
  //

  //
  //
  // Exact for polynomials up to order 1
  //
  //
  /** @brief Gaussian quadrature rule exact for polynomials up to order 1 */
  template <typename InterfaceType>
  class rt_gauss_quad_element <viennafem::unit_tetrahedron, 1, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_gauss_quad_element <viennafem::unit_tetrahedron, 1, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      explicit rt_gauss_quad_element() : p_(3)
      {
        p_[0] = 1.0/4.0;
        p_[1] = 1.0/4.0;
        p_[2] = 1.0/4.0;
      }

      BaseType * clone() const { return new self_type(); }

      NumericT eval(viennamath::rt_interval<InterfaceType> const & /*interv*/,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & /*var*/) const
      {
        return viennamath::eval(e, p_) / 6.0;
      }

    private:
      std::vector<NumericT> p_;
  };



  //
  //
  // Exact for polynomials up to degree 2
  //
  //
  /** @brief Keast rule, exact for polynomials up to degree 2.
   *
   *  Also have a look at the datasets by J. Burkardt at http://people.sc.fsu.edu/~jburkardt/datasets/datasets.html
   */
  template <typename InterfaceType>
  class rt_keast_quad_element <viennafem::unit_tetrahedron, 2, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_keast_quad_element <viennafem::unit_tetrahedron, 2, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      enum { num_points = 4 };

      explicit rt_keast_quad_element() : abscissas_(num_points, std::vector<numeric_type>(3))
      {
        abscissas_[0][0] = 0.5854101966249685; abscissas_[0][1] = 0.1381966011250105; abscissas_[0][2] = 0.1381966011250105;
        abscissas_[1][0] = 0.1381966011250105; abscissas_[1][1] = 0.1381966011250105; abscissas_[1][2] = 0.1381966011250105;
        abscissas_[2][0] = 0.1381966011250105; abscissas_[2][1] = 0.1381966011250105; abscissas_[2][2] = 0.5854101966249685;
        abscissas_[3][0] = 0.1381966011250105; abscissas_[3][1] = 0.5854101966249685; abscissas_[3][2] = 0.1381966011250105;
      }

      BaseType * clone() const { return new self_type(); }

      NumericT eval(viennamath::rt_interval<InterfaceType> const & /*interv*/,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & /*var*/) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<num_points; ++i)
          result += viennamath::eval(e, abscissas_[i]);
        return result / 24.0;
      }

    private:
      std::vector<std::vector<NumericT> > abscissas_;
  };



  //
  //
  // Exact for polynomials up to degree 3:
  //
  //

  /** @brief Keast rule, exact for polynomials up to degree 3. Uses a negative weight, thus be careful with numerical stability!
   *
   *  Also have a look at the datasets by J. Burkardt at http://people.sc.fsu.edu/~jburkardt/datasets/datasets.html
   */
  template <typename InterfaceType>
  class rt_keast_quad_element <viennafem::unit_tetrahedron, 3, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_keast_quad_element <viennafem::unit_tetrahedron, 3, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      enum { num_points = 5 };

      explicit rt_keast_quad_element() : abscissas_(num_points, std::vector<numeric_type>(3)),  weights_(num_points)
      {
        abscissas_[0][0] = 0.25;               abscissas_[0][1] = 0.25;               abscissas_[0][2] = 0.25;
        abscissas_[1][0] = 0.5;                abscissas_[1][1] = 0.1666666666666667; abscissas_[1][2] = 0.1666666666666667;
        abscissas_[2][0] = 0.1666666666666667; abscissas_[2][1] = 0.1666666666666667; abscissas_[2][2] = 0.1666666666666667;
        abscissas_[3][0] = 0.1666666666666667; abscissas_[3][1] = 0.1666666666666667; abscissas_[3][2] = 0.5;
        abscissas_[4][0] = 0.1666666666666667; abscissas_[4][1] = 0.5;                abscissas_[4][2] = 0.1666666666666667;

        weights_[0] = -0.8;
        weights_[1] = 0.45;
        weights_[2] = 0.45;
        weights_[3] = 0.45;
        weights_[4] = 0.45;
      }

      BaseType * clone() const { return new self_type(); }

      NumericT eval(viennamath::rt_interval<InterfaceType> const & /*interv*/,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & /*var*/) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<num_points; ++i)
          result += weights_[i] * viennamath::eval(e, abscissas_[i]);
        return result / 6.0;
      }

    private:
      std::vector<std::vector<NumericT> > abscissas_;
      std::vector<NumericT> weights_;
  };




  //
  //
  // Exact for polynomials up to degree 4:
  //
  //

  /** @brief Keast rule, exact for polynomials up to degree 4.
   *
   *  Also have a look at the datasets by J. Burkardt at http://people.sc.fsu.edu/~jburkardt/datasets/datasets.html
   */
  template <typename InterfaceType>
  class rt_keast_quad_element <viennafem::unit_tetrahedron, 4, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_keast_quad_element <viennafem::unit_tetrahedron, 4, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      enum { num_points = 14 };

      explicit rt_keast_quad_element() : abscissas_(num_points, std::vector<numeric_type>(3)),  weights_(num_points)
      {
        abscissas_[0][0] = 0.0; abscissas_[0][1] = 0.5; abscissas_[0][2] = 0.5;
        abscissas_[1][0] = 0.5; abscissas_[1][1] = 0.0; abscissas_[1][2] = 0.5;
        abscissas_[2][0] = 0.5; abscissas_[2][1] = 0.5; abscissas_[2][2] = 0.0;
        abscissas_[3][0] = 0.5; abscissas_[3][1] = 0.0; abscissas_[3][2] = 0.0;
        abscissas_[4][0] = 0.0; abscissas_[4][1] = 0.5; abscissas_[4][2] = 0.0;
        abscissas_[5][0] = 0.0; abscissas_[5][1] = 0.0; abscissas_[5][2] = 0.5;

        abscissas_[6][0] = 0.6984197043243866; abscissas_[6][1] = 0.1005267652252045; abscissas_[6][2] = 0.1005267652252045;
        abscissas_[7][0] = 0.1005267652252045; abscissas_[7][1] = 0.1005267652252045; abscissas_[7][2] = 0.1005267652252045;
        abscissas_[8][0] = 0.1005267652252045; abscissas_[8][1] = 0.1005267652252045; abscissas_[8][2] = 0.6984197043243866;
        abscissas_[9][0] = 0.1005267652252045; abscissas_[9][1] = 0.6984197043243866; abscissas_[9][2] = 0.1005267652252045;
        abscissas_[10][0] = 0.0568813795204234; abscissas_[10][1] = 0.3143728734931922; abscissas_[10][2] = 0.3143728734931922;
        abscissas_[11][0] = 0.3143728734931922; abscissas_[11][1] = 0.3143728734931922; abscissas_[11][2] = 0.3143728734931922;
        abscissas_[12][0] = 0.3143728734931922; abscissas_[12][1] = 0.3143728734931922; abscissas_[12][2] = 0.0568813795204234;
        abscissas_[13][0] = 0.3143728734931922; abscissas_[13][1] = 0.0568813795204234; abscissas_[13][2] = 0.3143728734931922;

        weights_[0] = 0.0190476190476190;
        weights_[1] = 0.0190476190476190;
        weights_[2] = 0.0190476190476190;
        weights_[3] = 0.0190476190476190;
        weights_[4] = 0.0190476190476190;
        weights_[5] = 0.0190476190476190;
        weights_[6] = 0.0885898247429807;
        weights_[7] = 0.0885898247429807;
        weights_[8] = 0.0885898247429807;
        weights_[9] = 0.0885898247429807;
        weights_[10] = 0.1328387466855907;
        weights_[11] = 0.1328387466855907;
        weights_[12] = 0.1328387466855907;
        weights_[13] = 0.1328387466855907;
      }

      BaseType * clone() const { return new self_type(); }

      NumericT eval(viennamath::rt_interval<InterfaceType> const & /*interv*/,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & /*var*/) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<num_points; ++i)
          result += weights_[i] * viennamath::eval(e, abscissas_[i]);
        return result / 6.0;
      }

    private:
      std::vector<std::vector<NumericT> > abscissas_;
      std::vector<NumericT> weights_;
  };





  //
  //
  // Exact for polynomials up to degree 5:
  //
  //

  /** @brief Keast rule, exact for polynomials up to degree 5.
   *
   *  Also have a look at the datasets by J. Burkardt at http://people.sc.fsu.edu/~jburkardt/datasets/datasets.html
   */
  template <typename InterfaceType>
  class rt_keast_quad_element <viennafem::unit_tetrahedron, 5, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_keast_quad_element <viennafem::unit_tetrahedron, 5, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      enum { num_points = 15 };

      explicit rt_keast_quad_element() : abscissas_(num_points, std::vector<numeric_type>(3)),  weights_(num_points)
      {
        abscissas_[0][0] = 0.25; abscissas_[0][1] = 0.25; abscissas_[0][2] = 0.25;

        abscissas_[1][0] = 0.0;     abscissas_[1][1] = 1.0/3.0; abscissas_[1][2] = 1.0/3.0;
        abscissas_[2][0] = 1.0/3.0; abscissas_[2][1] = 1.0/3.0; abscissas_[2][2] = 1.0/3.0;
        abscissas_[3][0] = 1.0/3.0; abscissas_[3][1] = 1.0/3.0; abscissas_[3][2] = 0.;
        abscissas_[4][0] = 1.0/3.0; abscissas_[4][1] = 0.;      abscissas_[4][2] = 1.0/3.0;

        abscissas_[5][0] = 0.7272727272727273; abscissas_[5][1] = 0.0909090909090909; abscissas_[5][2] = 0.0909090909090909;
        abscissas_[6][0] = 0.0909090909090909; abscissas_[6][1] = 0.0909090909090909; abscissas_[6][2] = 0.0909090909090909;
        abscissas_[7][0] = 0.0909090909090909; abscissas_[7][1] = 0.0909090909090909; abscissas_[7][2] = 0.7272727272727273;
        abscissas_[8][0] = 0.0909090909090909; abscissas_[8][1] = 0.7272727272727273; abscissas_[8][2] = 0.0909090909090909;

        abscissas_[ 9][0] = 0.4334498464263357; abscissas_[ 9][1] = 0.0665501535736643; abscissas_[ 9][2] = 0.0665501535736643;
        abscissas_[10][0] = 0.0665501535736643; abscissas_[10][1] = 0.4334498464263357; abscissas_[10][2] = 0.0665501535736643;
        abscissas_[11][0] = 0.0665501535736643; abscissas_[11][1] = 0.0665501535736643; abscissas_[11][2] = 0.4334498464263357;
        abscissas_[12][0] = 0.0665501535736643; abscissas_[12][1] = 0.4334498464263357; abscissas_[12][2] = 0.4334498464263357;
        abscissas_[13][0] = 0.4334498464263357; abscissas_[13][1] = 0.0665501535736643; abscissas_[13][2] = 0.4334498464263357;
        abscissas_[14][0] = 0.4334498464263357; abscissas_[14][1] = 0.4334498464263357; abscissas_[14][2] = 0.0665501535736643;

        weights_[0] = 0.1817020685825351;

        weights_[1] = 0.0361607142857143;
        weights_[2] = 0.0361607142857143;
        weights_[3] = 0.0361607142857143;
        weights_[4] = 0.0361607142857143;

        weights_[5] = 0.0698714945161738;
        weights_[6] = 0.0698714945161738;
        weights_[7] = 0.0698714945161738;
        weights_[8] = 0.0698714945161738;

        weights_[ 9] = 0.0656948493683187;
        weights_[10] = 0.0656948493683187;
        weights_[11] = 0.0656948493683187;
        weights_[12] = 0.0656948493683187;
        weights_[13] = 0.0656948493683187;
        weights_[14] = 0.0656948493683187;
      }

      BaseType * clone() const { return new self_type(); }

      NumericT eval(viennamath::rt_interval<InterfaceType> const & /*interv*/,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & /*var*/) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<num_points; ++i)
          result += weights_[i] * viennamath::eval(e, abscissas_[i]);
        return result / 6.0;
      }

    private:
      std::vector<std::vector<NumericT> > abscissas_;
      std::vector<NumericT> weights_;
  };

  //
  //
  // Exact for polynomials up to degree 6:
  //
  //

  /** @brief Keast rule, exact for polynomials up to degree 6.
   *
   *  Also have a look at the datasets by J. Burkardt at http://people.sc.fsu.edu/~jburkardt/datasets/datasets.html
   */
  template <typename InterfaceType>
  class rt_keast_quad_element <viennafem::unit_tetrahedron, 6, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_keast_quad_element <viennafem::unit_tetrahedron, 6, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      enum { num_points = 24 };

      explicit rt_keast_quad_element() : abscissas_(num_points, std::vector<numeric_type>(3)),  weights_(num_points)
      {
        abscissas_[ 0][0] = 0.3561913862225449; abscissas_[ 0][1] = 0.2146028712591517; abscissas_[ 0][2] = 0.2146028712591517;
        abscissas_[ 1][0] = 0.2146028712591517; abscissas_[ 1][1] = 0.2146028712591517; abscissas_[ 1][2] = 0.2146028712591517;
        abscissas_[ 2][0] = 0.2146028712591517; abscissas_[ 2][1] = 0.2146028712591517; abscissas_[ 2][2] = 0.3561913862225449;
        abscissas_[ 3][0] = 0.2146028712591517; abscissas_[ 3][1] = 0.3561913862225449; abscissas_[ 3][2] = 0.2146028712591517;

        abscissas_[ 4][0] = 0.8779781243961660; abscissas_[ 4][1] = 0.0406739585346113; abscissas_[ 4][2] = 0.0406739585346113;
        abscissas_[ 5][0] = 0.0406739585346113; abscissas_[ 5][1] = 0.0406739585346113; abscissas_[ 5][2] = 0.0406739585346113;
        abscissas_[ 6][0] = 0.0406739585346113; abscissas_[ 6][1] = 0.0406739585346113; abscissas_[ 6][2] = 0.8779781243961660;
        abscissas_[ 7][0] = 0.0406739585346113; abscissas_[ 7][1] = 0.8779781243961660; abscissas_[ 7][2] = 0.0406739585346113;

        abscissas_[ 8][0] = 0.0329863295731731; abscissas_[ 8][1] = 0.3223378901422757; abscissas_[ 8][2] = 0.3223378901422757;
        abscissas_[ 9][0] = 0.3223378901422757; abscissas_[ 9][1] = 0.3223378901422757; abscissas_[ 9][2] = 0.3223378901422757;
        abscissas_[10][0] = 0.3223378901422757; abscissas_[10][1] = 0.3223378901422757; abscissas_[10][2] = 0.0329863295731731;
        abscissas_[11][0] = 0.3223378901422757; abscissas_[11][1] = 0.0329863295731731; abscissas_[11][2] = 0.3223378901422757;

        abscissas_[12][0] = 0.2696723314583159; abscissas_[12][1] = 0.0636610018750175; abscissas_[12][2] = 0.0636610018750175;
        abscissas_[13][0] = 0.0636610018750175; abscissas_[13][1] = 0.2696723314583159; abscissas_[13][2] = 0.0636610018750175;
        abscissas_[14][0] = 0.0636610018750175; abscissas_[14][1] = 0.0636610018750175; abscissas_[14][2] = 0.2696723314583159;
        abscissas_[15][0] = 0.6030056647916491; abscissas_[15][1] = 0.0636610018750175; abscissas_[15][2] = 0.0636610018750175;

        abscissas_[16][0] = 0.0636610018750175; abscissas_[16][1] = 0.6030056647916491; abscissas_[16][2] = 0.0636610018750175;
        abscissas_[17][0] = 0.0636610018750175; abscissas_[17][1] = 0.0636610018750175; abscissas_[17][2] = 0.6030056647916491;
        abscissas_[18][0] = 0.0636610018750175; abscissas_[18][1] = 0.2696723314583159; abscissas_[18][2] = 0.6030056647916491;
        abscissas_[19][0] = 0.2696723314583159; abscissas_[19][1] = 0.6030056647916491; abscissas_[19][2] = 0.0636610018750175;

        abscissas_[20][0] = 0.6030056647916491; abscissas_[20][1] = 0.0636610018750175; abscissas_[20][2] = 0.2696723314583159;
        abscissas_[21][0] = 0.0636610018750175; abscissas_[21][1] = 0.6030056647916491; abscissas_[21][2] = 0.2696723314583159;
        abscissas_[22][0] = 0.2696723314583159; abscissas_[22][1] = 0.0636610018750175; abscissas_[22][2] = 0.6030056647916491;
        abscissas_[23][0] = 0.6030056647916491; abscissas_[23][1] = 0.2696723314583159; abscissas_[23][2] = 0.0636610018750175;

        weights_[ 0] = 0.0399227502581679;
        weights_[ 1] = 0.0399227502581679;
        weights_[ 2] = 0.0399227502581679;
        weights_[ 3] = 0.0399227502581679;

        weights_[ 4] = 0.0100772110553207;
        weights_[ 5] = 0.0100772110553207;
        weights_[ 6] = 0.0100772110553207;
        weights_[ 7] = 0.0100772110553207;

        weights_[ 8] = 0.0553571815436544;
        weights_[ 9] = 0.0553571815436544;
        weights_[10] = 0.0553571815436544;
        weights_[11] = 0.0553571815436544;

        weights_[12] = 0.0482142857142857;
        weights_[13] = 0.0482142857142857;
        weights_[14] = 0.0482142857142857;
        weights_[15] = 0.0482142857142857;

        weights_[16] = 0.0482142857142857;
        weights_[17] = 0.0482142857142857;
        weights_[18] = 0.0482142857142857;
        weights_[19] = 0.0482142857142857;

        weights_[20] = 0.0482142857142857;
        weights_[21] = 0.0482142857142857;
        weights_[22] = 0.0482142857142857;
        weights_[23] = 0.0482142857142857;
      }

      BaseType * clone() const { return new self_type(); }

      NumericT eval(viennamath::rt_interval<InterfaceType> const & /*interv*/,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & /*var*/) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<num_points; ++i)
          result += weights_[i] * viennamath::eval(e, abscissas_[i]);
        return result / 6.0;
      }

    private:
      std::vector<std::vector<NumericT> > abscissas_;
      std::vector<NumericT> weights_;
  };

}
#endif
