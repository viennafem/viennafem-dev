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
  
  //TODO: compile time derivation
  
  template <typename InterfaceType>
  struct weak_form_checker : public viennamath::traversal_interface<>
  {
    public:
      weak_form_checker() : is_weak_(false) {}
      
      void operator()(InterfaceType const * e) const 
      {
        if (viennamath::callback_if_castable< viennamath::unary_expr<InterfaceType> >::apply(e, *this))
          return;
      }
      
      void operator()(viennamath::unary_expr<InterfaceType> const & unary_expr) const
      {
        typedef viennamath::op_unary<viennamath::op_symbolic_integration<typename InterfaceType::numeric_type>, InterfaceType>  SymbolicIntegrationType;
        
        if (dynamic_cast<const SymbolicIntegrationType *>(unary_expr.op()) != NULL)
          is_weak_ = true;        
      }
      
      bool is_weak() const { return is_weak_; }
      
    private:
      mutable bool is_weak_; //TODO: Should not be necessary here...
  };
  
  
  template <typename InterfaceType>
  bool is_weak_form(viennamath::equation<InterfaceType> const & strong_formulation)
  {
    weak_form_checker<InterfaceType> * checker = new weak_form_checker<InterfaceType>();
    
    viennamath::traversal_wrapper<> wrapped_checker( checker ); //Note: checker pointer is auto-ptr'ed here, therefore no need for an explicit 'delete checker;' at the end of the function.

    strong_formulation.lhs().get()->recursive_traversal( wrapped_checker );
    bool result = checker->is_weak();
      
    return result;
  }

  //tries to automatically derive the weak formulation from the strong formulation
  template <typename InterfaceType>
  viennamath::equation<InterfaceType> make_weak_form(viennamath::equation<InterfaceType> const & strong_formulation)
  {
    if (is_weak_form(strong_formulation))
    {
      std::cout << "make_weak_form(): Nothing to do, problem already in weak form!" << std::endl;
      return strong_formulation;
    }
    
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