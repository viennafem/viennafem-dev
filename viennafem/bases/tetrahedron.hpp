#ifndef VIENNAFEM_BASES_TETRAHEDRON_HPP
#define VIENNAFEM_BASES_TETRAHEDRON_HPP

/* =========================================================================
   Copyright (c) 2012-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
                             -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                             -----------------

   Author:     Karl Rupp                          rupp@iue.tuwien.ac.at

   License:    MIT (X11), see file LICENSE in the ViennaFEM base directory
============================================================================ */

#include <vector>
#include "viennagrid/topology/simplex.hpp"
#include "viennamath/expression.hpp"
#include "viennafem/forwards.h"

/** @file   viennafem/bases/tetrahedron.hpp
    @brief  Defines the various basis functions for tetrahedra.
*/

namespace viennafem
{


  //
  // Lagrange family on tetrahedra
  //


  // Vertex basis:
  /** @brief Returns the first (linear) vertex basis function */
  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_tetrahedron,
                      0,   // vertex level
                      0>   // vertex (0,0,0)
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;

    typedef viennamath::ct_binary_expr< viennamath::ct_binary_expr< viennamath::ct_constant<1>,
                                                                    viennamath::op_minus<viennafem::numeric_type>,
                                                                    viennamath::ct_variable<0> >,
                                        viennamath::op_minus<viennafem::numeric_type>,
                                        viennamath::ct_binary_expr< viennamath::ct_variable<1>,
                                                                    viennamath::op_plus<viennafem::numeric_type>,
                                                                    viennamath::ct_variable<2> >
                                      >                 type;

    static expression_type get() { return expression_type(type()); }
  };

  /** @brief Returns the second (linear) vertex basis function */
  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_tetrahedron,
                      0,   //vertex level
                      1>   //vertex (1,0,0)
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;

    typedef viennamath::ct_variable<0>           type;

    static expression_type get() { return expression_type(type()); }
  };


  /** @brief Returns the third (linear) vertex basis function */
  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_tetrahedron,
                      0,   //vertex level
                      2>   //vertex (0,1,0)
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;

    typedef viennamath::ct_variable<1>           type;

    static expression_type get() { return expression_type(type()); }
  };

  /** @brief Returns the forth (linear) vertex basis function */
  template <typename InterfaceType, std::size_t order>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<order>,
                      unit_tetrahedron,
                      0,   //vertex level
                      3>   //vertex (0,0,1)
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;

    typedef viennamath::ct_variable<2>           type;

    static expression_type get() { return expression_type(type()); }
  };


  //
  // quadratic
  //
  /** @brief Returns the first (quadratic) edge basis function */
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_tetrahedron,
                      1,   //edge level
                      0>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_tetrahedron,
                                  0,
                                  0>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_tetrahedron,
                                  0,
                                  1>::type       phi_1;

    //x * (1-x)
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;

    static expression_type get() { return expression_type(type()); }
  };

  /** @brief Returns the second (quadratic) edge basis function */
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_tetrahedron,
                      1,   //edge level
                      1>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_tetrahedron,
                                  0,
                                  0>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_tetrahedron,
                                  0,
                                  2>::type       phi_1;

    //x * (1-x)
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;

    static expression_type get() { return expression_type(type()); }
  };

  /** @brief Returns the third (quadratic) edge basis function */
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_tetrahedron,
                      1,   //edge level
                      2>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_tetrahedron,
                                  0,
                                  0>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_tetrahedron,
                                  0,
                                  3>::type       phi_1;

    //x * (1-x)
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;

    static expression_type get() { return expression_type(type()); }
  };

  /** @brief Returns the forth (quadratic) edge basis function */
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_tetrahedron,
                      1,   //edge level
                      3>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_tetrahedron,
                                  0,
                                  1>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_tetrahedron,
                                  0,
                                  2>::type       phi_1;

    //x * (1-x)
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;

    static expression_type get() { return expression_type(type()); }
  };

  /** @brief Returns the fifth (quadratic) edge basis function */
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_tetrahedron,
                      1,   //edge level
                      4>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_tetrahedron,
                                  0,
                                  1>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_tetrahedron,
                                  0,
                                  3>::type       phi_1;

    //x * (1-x)
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;

    static expression_type get() { return expression_type(type()); }
  };

  /** @brief Returns the sixth (quadratic) edge basis function */
  template <typename InterfaceType>
  struct local_basis <InterfaceType,
                      viennafem::lagrange_tag<2>,
                      unit_tetrahedron,
                      1,   //edge level
                      5>
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef typename InterfaceType::numeric_type NumericT;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_tetrahedron,
                                  0,
                                  2>::type       phi_0;

    typedef typename local_basis <InterfaceType,
                                  viennafem::lagrange_tag<2>,
                                  unit_tetrahedron,
                                  0,
                                  3>::type       phi_1;

    //x * (1-x)
    typedef viennamath::ct_binary_expr<phi_0,
                                       viennamath::op_mult<NumericT>,
                                       phi_1
                                      > type;

    static expression_type get() { return expression_type(type()); }
  };

}

#endif
