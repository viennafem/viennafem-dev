/* =======================================================================
   Copyright (c) 2010, Institute for Microelectronics, TU Vienna.
   http://www.iue.tuwien.ac.at
                             -----------------
                 ViennaMath - Symbolic and Numeric Math in C++
                             -----------------

   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaMath base directory
======================================================================= */



#ifndef VIENNAMATH_WEAK_FORM_HPP
#define VIENNAMATH_WEAK_FORM_HPP

#include "viennamath/expression.hpp"
#include "viennamath/manipulation/substitute.hpp"
#include "viennamath/manipulation/integral.hpp"

namespace viennafem
{
  //TODO: Is this the right place for the derivation of a weak formulation? Shouldn't this go into ViennaFEM or ViennaFVM?
  
  
  //TODO: compile time derivation
  


  //tries to automatically derive the weak formulation from the strong formulation
  template <typename InterfaceType>
  viennamath::equation<InterfaceType> make_weak_form(viennamath::equation<InterfaceType> const & strong_formulation)
  {
    //TODO: More general derivations: Transform div(expr) to expr * grad(v)
    viennamath::expr<InterfaceType> new_lhs(viennamath::substitute( viennamath::laplace(viennamath::function_symbol<InterfaceType>(0, viennamath::unknown_tag<0>())),
                                            viennamath::constant<typename InterfaceType::numeric_type, InterfaceType>(-1) * (viennamath::grad(viennamath::function_symbol<InterfaceType>(0, viennamath::unknown_tag<0>())) * viennamath::grad(viennamath::function_symbol<InterfaceType>(0, viennamath::test_tag<0>()))),
                                            strong_formulation.lhs()
                                          )
                               );
    return viennamath::equation<InterfaceType>( viennamath::integral(viennamath::Omega(), new_lhs, viennamath::symbolic_tag()),
                                    viennamath::integral(viennamath::Omega(), strong_formulation.rhs() * viennamath::function_symbol<InterfaceType>(0, viennamath::test_tag<0>()), viennamath::symbolic_tag())
                                  );
  }

}

#endif