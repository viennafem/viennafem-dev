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
  //
  // Compile time transformation
  //
  
  //TODO



  //
  // Run time transformation
  //

  template <typename InterfaceType>
  struct weak_form_checker : public viennamath::rt_traversal_interface<>
  {
    public:
      weak_form_checker() : is_weak_(false) {}
      
      void operator()(InterfaceType const * e) const 
      {
        if (viennamath::callback_if_castable< viennamath::rt_unary_expr<InterfaceType> >::apply(e, *this))
          return;
      }
      
      void operator()(viennamath::rt_unary_expr<InterfaceType> const & unary_expr) const
      {
        typedef viennamath::op_unary<viennamath::op_rt_integral<InterfaceType>, InterfaceType>  SymbolicIntegrationType;
        
        if (dynamic_cast<const SymbolicIntegrationType *>(unary_expr.op()) != NULL)
          is_weak_ = true;        
      }
      
      bool is_weak() const { return is_weak_; }
      
    private:
      mutable bool is_weak_; //TODO: Should not be necessary here...
  };
  
  
  template <typename InterfaceType>
  bool is_weak_form(viennamath::rt_equation<InterfaceType> const & strong_formulation)
  {
    weak_form_checker<InterfaceType> * checker = new weak_form_checker<InterfaceType>();
    
    viennamath::rt_traversal_wrapper<> wrapped_checker( checker ); //Note: checker pointer is auto-ptr'ed here, therefore no need for an explicit 'delete checker;' at the end of the function.

    strong_formulation.lhs().get()->recursive_traversal( wrapped_checker );
    bool result = checker->is_weak();
      
    return result;
  }







  template <typename InterfaceType>
  struct weak_form_creator : public viennamath::rt_manipulation_interface<InterfaceType>
  {
    public:
      InterfaceType * operator()(InterfaceType const * e) const 
      {
        if (   !viennamath::callback_if_castable< viennamath::rt_unary_expr<InterfaceType> >::apply(e, *this)
            && !viennamath::callback_if_castable< viennamath::rt_binary_expr<InterfaceType> >::apply(e, *this))
        {
          //this if-body is only executed for trivial expressions such as 'u' (L^2-projection)
          
          //multiply with test function and integrate
          viennamath::rt_function_symbol<InterfaceType> test_func(0, viennamath::test_tag<0>());
          viennamath::rt_expr<InterfaceType> temp(e->clone());
          viennamath::rt_interval<InterfaceType> interv(true);
          integrated_expr = viennamath::integral(interv, temp * test_func, viennamath::ct_variable<0>());
        }
        
        return integrated_expr.get()->clone();
        
      }
      
      void operator()(viennamath::rt_unary_expr<InterfaceType> const & unary_expr) const
      {
        typedef typename InterfaceType::numeric_type   NumericType;
        typedef viennamath::op_unary<viennamath::op_divergence<NumericType>, InterfaceType>  DivergenceOperatorType;
        
        if (dynamic_cast<const DivergenceOperatorType *>(unary_expr.op()) != NULL) //this is something of the form div(expr) with some expression expr
        {
          //std::cout << "Detected divergence operator!" << std::endl;
          viennamath::rt_expr<InterfaceType> minus1 = -1;
          viennamath::rt_expr<InterfaceType> lhs(unary_expr.lhs()->clone());
          viennamath::rt_expr<InterfaceType> rhs = viennamath::grad(viennamath::rt_function_symbol<InterfaceType>(0, viennamath::test_tag<0>()));
          viennamath::rt_interval<InterfaceType> interv(true);
          integrated_expr = viennamath::integral(interv, minus1 * (lhs * rhs), viennamath::ct_variable<0>());
        }
        else
          throw "Cannot derive weak form!";
      }

      void operator()(viennamath::rt_binary_expr<InterfaceType> const & bin) const
      {
        typedef typename InterfaceType::numeric_type   NumericType;
        typedef viennamath::op_binary<viennamath::op_plus<NumericType>, InterfaceType>   PlusOperatorType;
        typedef viennamath::op_binary<viennamath::op_minus<NumericType>, InterfaceType>  MinusOperatorType;
        typedef viennamath::op_binary<viennamath::op_mult<NumericType>, InterfaceType>   ProductOperatorType;
        typedef viennamath::op_binary<viennamath::op_div<NumericType>, InterfaceType>  DivisionOperatorType;
        
        if (    dynamic_cast<const PlusOperatorType *>(bin.op()) != NULL
             || dynamic_cast<const MinusOperatorType *>(bin.op()) != NULL) //integration is additive :-)
        {
          viennamath::rt_manipulation_wrapper<InterfaceType> manipulator(new weak_form_creator<InterfaceType>());
          //Note: In the following, directly passing *this is not possible due to the need for a wrapper...
          integrated_expr = new viennamath::rt_binary_expr<InterfaceType>(bin.lhs()->recursive_manipulation(manipulator),
                                                                       bin.op()->clone(),
                                                                       bin.rhs()->recursive_manipulation(manipulator));
        }
        else if (dynamic_cast<const ProductOperatorType *>(bin.op()) != NULL
                  || dynamic_cast<const DivisionOperatorType *>(bin.op()) != NULL)
        {
          //handle the cases const * div(...) and div(...) * const
          
          if (bin.lhs()->is_constant())
          {
            viennamath::rt_manipulation_wrapper<InterfaceType> manipulator(new weak_form_creator<InterfaceType>());
            //Note: In the following, directly passing *this is not possible due to the need for a wrapper...
            integrated_expr = new viennamath::rt_binary_expr<InterfaceType>(bin.lhs()->clone(),
                                                                         bin.op()->clone(),
                                                                         bin.rhs()->recursive_manipulation(manipulator));
          }
          else if (bin.rhs()->is_constant())
          {
            viennamath::rt_manipulation_wrapper<InterfaceType> manipulator(new weak_form_creator<InterfaceType>());
            //Note: In the following, directly passing *this is not possible due to the need for a wrapper...
            integrated_expr = new viennamath::rt_binary_expr<InterfaceType>(bin.lhs()->recursive_manipulation(manipulator),
                                                                         bin.op()->clone(),
                                                                         bin.rhs()->clone());
          }
          else
          {
            throw "cannot derive weak form!";
          }
        }
        else //TODO: Add checks!
        {
          //multiply with test function and integrate
          viennamath::rt_function_symbol<InterfaceType> test_func(0, viennamath::test_tag<0>());
          viennamath::rt_interval<InterfaceType> interv(true);
          integrated_expr = viennamath::integral(interv, bin * test_func, viennamath::ct_variable<0>());
        }
      }

      bool modifies(InterfaceType const * e) const { return true; }
      
    private:
      mutable viennamath::rt_expr<InterfaceType> integrated_expr;
  };


  //tries to automatically derive the weak formulation from the strong formulation
  template <typename InterfaceType>
  viennamath::rt_equation<InterfaceType> make_weak_form(viennamath::rt_equation<InterfaceType> const & strong_formulation)
  {
    if (is_weak_form(strong_formulation))
    {
      std::cout << "make_weak_form(): Nothing to do, problem already in weak form!" << std::endl;
      return strong_formulation;
    }
    
    //TODO: More general derivations: Transform div(expr) to expr * grad(v)
    /*viennamath::expr<InterfaceType> new_lhs(viennamath::substitute( viennamath::laplace(viennamath::function_symbol<InterfaceType>(0, viennamath::unknown_tag<0>())),
                                            viennamath::constant<typename InterfaceType::numeric_type, InterfaceType>(-1) * (viennamath::grad(viennamath::function_symbol<InterfaceType>(0, viennamath::unknown_tag<0>())) * viennamath::grad(viennamath::function_symbol<InterfaceType>(0, viennamath::test_tag<0>()))),
                                            strong_formulation.lhs()
                                          )
                               );
    return viennamath::equation<InterfaceType>( viennamath::integral(viennamath::Omega(), new_lhs, viennamath::symbolic_tag()),
                                    viennamath::integral(viennamath::Omega(), strong_formulation.rhs() * viennamath::function_symbol<InterfaceType>(0, viennamath::test_tag<0>()), viennamath::symbolic_tag())
                                  );*/
    
    viennamath::rt_interval<InterfaceType> interv(true);
    viennamath::rt_manipulation_wrapper<InterfaceType> wrapped_checker( new weak_form_creator<InterfaceType>() );
    viennamath::rt_expr<InterfaceType> weak_lhs(strong_formulation.lhs().get()->recursive_manipulation( wrapped_checker ));
    viennamath::rt_expr<InterfaceType> weak_rhs = 
        viennamath::integral(interv,
                             strong_formulation.rhs() * viennamath::rt_function_symbol<InterfaceType>(0, viennamath::test_tag<0>()),
                             viennamath::ct_variable<0>());
    
    return viennamath::rt_equation<InterfaceType>( weak_lhs, weak_rhs);
    
  }

}

#endif