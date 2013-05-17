#ifndef VIENNAMATH_WEAK_FORM_HPP
#define VIENNAMATH_WEAK_FORM_HPP

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



#include "viennamath/expression.hpp"
#include "viennamath/manipulation/substitute.hpp"

/** @file   weak_form.hpp
    @brief  Derives the weak form of a PDE given in strong form.
*/

namespace viennafem
{
  //
  // Compile time transformation
  //
  
  //TODO: Not finished in time for the 1.0.0 release. Scheduled for ViennaFEM 1.1.0.



  //
  // Run time transformation
  //

  namespace detail
  {
    
    /** @brief A helper class which scans whether a ViennaMath equation is in a weak form already.
    * 
    * Currently only checks whether there is a symbolic integral somewhere. This should be sufficient as long as there are no integral equations passed, for which FEM is probably not useful at all...
    */
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
          typedef viennamath::op_unary<viennamath::op_rt_symbolic_integral<InterfaceType>, InterfaceType>  SymbolicIntegrationType;
          
          if (dynamic_cast<const SymbolicIntegrationType *>(unary_expr.op()) != NULL)
            is_weak_ = true;        
        }
        
        bool is_weak() const { return is_weak_; }
        
      private:
        mutable bool is_weak_; //TODO: Should not be necessary here...
    };
    
  } //namespace detail
  
  /** @brief Interface function for checking whether a certan equation is a weak formulation of a PDE. */
  template <typename InterfaceType>
  bool is_weak_form(viennamath::rt_equation<InterfaceType> const & strong_formulation)
  {
    detail::weak_form_checker<InterfaceType> * checker = new detail::weak_form_checker<InterfaceType>();
    
    viennamath::rt_traversal_wrapper<> wrapped_checker( checker ); //Note: checker pointer is auto-ptr'ed here, therefore no need for an explicit 'delete checker;' at the end of the function.

    strong_formulation.lhs().get()->recursive_traversal( wrapped_checker );
    return checker->is_weak();
  }






  namespace detail
  {
    /** @brief Transforms a strong formulation of an equation to a weak form, assuming homogeneous Neumann boundary conditions */
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
            integrated_expr = viennamath::integral(viennamath::symbolic_interval(), temp * test_func);
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
            integrated_expr = viennamath::integral(viennamath::symbolic_interval(),
                                                  minus1 * (lhs * rhs)            
                                                  );
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
            integrated_expr = viennamath::integral(viennamath::symbolic_interval(),
                                                  bin * test_func);
          }
        }

        bool modifies(InterfaceType const * /*e*/) const { return true; }
        
      private:
        mutable viennamath::rt_expr<InterfaceType> integrated_expr;
    };

  } //namespace detail

  //tries to automatically derive the weak formulation from the strong formulation
  /** @brief Transforms a ViennaMath equation from strong form to weak form. */
  template <typename InterfaceType>
  viennamath::rt_equation<InterfaceType> make_weak_form(viennamath::rt_equation<InterfaceType> const & strong_formulation)
  {
    if (is_weak_form(strong_formulation))
    {
      std::cout << "ViennaFEM: make_weak_form(): Nothing to do, problem already in weak form!" << std::endl;
      return strong_formulation;
    }
    
    viennamath::rt_manipulation_wrapper<InterfaceType> wrapped_checker( new detail::weak_form_creator<InterfaceType>() );
    viennamath::rt_expr<InterfaceType> weak_lhs(strong_formulation.lhs().get()->recursive_manipulation( wrapped_checker ));
    viennamath::rt_expr<InterfaceType> weak_rhs = 
        viennamath::integral(viennamath::symbolic_interval(),
                             strong_formulation.rhs() * viennamath::rt_function_symbol<InterfaceType>(0, viennamath::test_tag<0>())
                            );
    
    return viennamath::rt_equation<InterfaceType>(weak_lhs, weak_rhs);
    
  }

}

#endif
