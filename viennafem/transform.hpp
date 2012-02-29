#ifndef VIENNAFEM_TRANSFORM_HPP
#define VIENNAFEM_TRANSFORM_HPP

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

#include "viennamath/expression.hpp"
#include "viennamath/manipulation/substitute.hpp"
#include "viennamath/manipulation/diff.hpp"

#include "viennagrid/topology/triangle.hpp"
#include "viennagrid/topology/tetrahedron.hpp"


/** @file   transform.hpp
    @brief  Provides transformation routines for converting integrals on elements in physical space to elements in reference orientation.
*/

namespace viennafem
{

  namespace detail
  {
    /** @brief Wraps the cell quantity for the Jacobian on each cell depending on the reference element used. This is the overload for the unit interval. */
    template <typename CellType, typename InterfaceType>
    void wrap_jacobian(viennafem::cell_quan<CellType, InterfaceType> & det_dF_dt, viennafem::unit_interval)
    {
      det_dF_dt.wrap_constant( viennafem::det_dF_dt_key() );
    }

    /** @brief Wraps the cell quantity for the Jacobian on each cell depending on the reference element used. This is the overload for the unit triangle. */
    template <typename CellType, typename InterfaceType>
    void wrap_jacobian(viennafem::cell_quan<CellType, InterfaceType> & det_dF_dt, viennafem::unit_triangle)
    {
      det_dF_dt.wrap_constant( viennafem::det_dF_dt_key() );
    }

    /** @brief Wraps the cell quantity for the Jacobian on each cell depending on the reference element used. This is the overload for the unit tetrahedron. */
    template <typename CellType, typename InterfaceType>
    void wrap_jacobian(viennafem::cell_quan<CellType, InterfaceType> & det_dF_dt, viennafem::unit_tetrahedron)
    {
      det_dF_dt.wrap_constant( viennafem::det_dF_dt_key() );
    }

    /** @brief Wraps the cell quantity for the Jacobian on each cell depending on the reference element used. This is the overload for the unit square. */
    template <typename CellType, typename InterfaceType>
    void wrap_jacobian(viennafem::cell_quan<CellType, InterfaceType> & det_dF_dt, viennafem::unit_square)
    {
      det_dF_dt.wrap_expr( viennafem::det_dF_dt_key() );
    }

    /** @brief Wraps the cell quantity for the Jacobian on each cell depending on the reference element used. This is the overload for the unit cube. */
    template <typename CellType, typename InterfaceType>
    void wrap_jacobian(viennafem::cell_quan<CellType, InterfaceType> & det_dF_dt, viennafem::unit_cube)
    {
      det_dF_dt.wrap_expr( viennafem::det_dF_dt_key() );
    }
    
    /** @brief Multiplies a given expression with the Jacobian cell-quantity. If the expression is given by an integral, the integrand is multiplied with the Jacobian. */
    template <typename CellType, typename InterfaceType, typename ReferenceCellTag>
    struct jacobian_adder : public viennamath::rt_manipulation_interface<InterfaceType>
    {
      public:
        typedef typename InterfaceType::numeric_type        numeric_type;
        
        InterfaceType * operator()(InterfaceType const * e) const 
        {
          const viennamath::rt_unary_expr<InterfaceType> * integral_expression = dynamic_cast<const viennamath::rt_unary_expr<InterfaceType> *>(e);
          
          if (integral_expression != NULL)
          {
            typedef viennamath::op_unary<viennamath::op_rt_symbolic_integral<InterfaceType>, InterfaceType>    SymbolicIntegrationOperation;

            const SymbolicIntegrationOperation * op_symb_tmp = dynamic_cast<const SymbolicIntegrationOperation *>(integral_expression->op());

            viennafem::cell_quan<CellType, InterfaceType> det_dF_dt;
            detail::wrap_jacobian(det_dF_dt, ReferenceCellTag());
            
            if (op_symb_tmp != NULL)
            {
              return new viennamath::rt_unary_expr<InterfaceType>(
                          new viennamath::rt_binary_expr<InterfaceType>(integral_expression->lhs()->clone(),
                                                                        new viennamath::op_binary<viennamath::op_mult<numeric_type>, InterfaceType>(),
                                                                        det_dF_dt.clone()),
                          op_symb_tmp->clone());
            }
          }
          
          return e->clone();
        }
        
        bool modifies(InterfaceType const * e) const
        {
          const viennamath::rt_unary_expr<InterfaceType> * integral_expression = dynamic_cast<const viennamath::rt_unary_expr<InterfaceType> *>(e);
          
          if (integral_expression != NULL)
          {
            typedef viennamath::op_unary<viennamath::op_rt_symbolic_integral<InterfaceType>, InterfaceType>    SymbolicIntegrationOperation;

            const SymbolicIntegrationOperation * op_symb_tmp = dynamic_cast<const SymbolicIntegrationOperation *>(integral_expression->op());
            
            if (op_symb_tmp != NULL)
              return true;
          }
          return false;
        }
    };

    
    
    /** @brief Wraps the cell quantity for the partial derivatives of the reference mapping on each cell depending on the reference element used. This is the overload for the unit interval. */
    template <long i, long j, typename CellType, typename InterfaceType>
    void wrap_dt_dx(viennafem::cell_quan<CellType, InterfaceType> & dt_dx, viennafem::unit_interval)
    {
      dt_dx.wrap_constant( viennafem::dt_dx_key<i,j>() );
    }
    
    /** @brief Wraps the cell quantity for the partial derivatives of the reference mapping on each cell depending on the reference element used. This is the overload for the unit triangle. */
    template <long i, long j, typename CellType, typename InterfaceType>
    void wrap_dt_dx(viennafem::cell_quan<CellType, InterfaceType> & dt_dx, viennafem::unit_triangle)
    {
      dt_dx.wrap_constant( viennafem::dt_dx_key<i,j>() );
    }

    /** @brief Wraps the cell quantity for the partial derivatives of the reference mapping on each cell depending on the reference element used. This is the overload for the unit tetrahedron. */
    template <long i, long j, typename CellType, typename InterfaceType>
    void wrap_dt_dx(viennafem::cell_quan<CellType, InterfaceType> & dt_dx, viennafem::unit_tetrahedron)
    {
      dt_dx.wrap_constant( viennafem::dt_dx_key<i,j>() );
    }

    /** @brief Wraps the cell quantity for the partial derivatives of the reference mapping on each cell depending on the reference element used. This is the overload for the unit square. */
    template <long i, long j, typename CellType, typename InterfaceType>
    void wrap_dt_dx(viennafem::cell_quan<CellType, InterfaceType> & dt_dx, viennafem::unit_square)
    {
      dt_dx.wrap_expr( viennafem::dt_dx_key<i,j>() );
    }

    /** @brief Wraps the cell quantity for the partial derivatives of the reference mapping on each cell depending on the reference element used. This is the overload for the unit cube. */
    template <long i, long j, typename CellType, typename InterfaceType>
    void wrap_dt_dx(viennafem::cell_quan<CellType, InterfaceType> & dt_dx, viennafem::unit_cube)
    {
      dt_dx.wrap_expr( viennafem::dt_dx_key<i,j>() );
    }
    
    //////////////////////////         1d          //////////////////////////////
    
    //
    // Line:
    //
    /** @brief Transforms a weak formulation given in physical space to the reference cell. Also transforms derivatives where appropriate. Implementation for the 1d-case.*/
    template <typename CellType, typename EquationType, typename PDESystemType, typename ReferenceCellTag>
    EquationType transform_to_reference_cell_1d(EquationType const & weak_form,
                                                PDESystemType const & pde_system,
                                                ReferenceCellTag)
    {
      typedef typename EquationType::interface_type             InterfaceType;
      typedef typename InterfaceType::numeric_type              numeric_type;
      typedef viennamath::rt_function_symbol<InterfaceType>   FunctionSymbol;
      typedef viennamath::rt_expr<InterfaceType>              Expression;
      typedef viennamath::rt_unary_expr<InterfaceType>        UnaryExpression;
      typedef viennamath::rt_variable<InterfaceType>          Variable;
      typedef viennamath::op_unary<viennamath::op_partial_deriv<viennamath::default_numeric_type>, InterfaceType >   d_dx;
      
      
      
      FunctionSymbol  u(0, viennamath::unknown_tag<>());
      FunctionSymbol  v(0, viennamath::test_tag<>());
      Variable        r(0);
      
      //using local coordinates (r, s) and global coordinates (x, y)
      viennafem::cell_quan<CellType, InterfaceType>  dr_dx; detail::wrap_dt_dx<0,0>(dr_dx, ReferenceCellTag());
      
      Expression new_lhs = weak_form.lhs();
      Expression new_rhs = weak_form.rhs();
      
      //std::cout << "New lhs init: " << std::endl;
      //std::cout << new_lhs << std::endl;
      
      std::vector<Expression> search_types(2);
      search_types[0] = UnaryExpression(u.clone(), new d_dx(0));
      search_types[1] = UnaryExpression(v.clone(), new d_dx(0));
      
      std::vector<Expression> replace_types(2);
      replace_types[0] = viennamath::diff(u, r) * dr_dx;
      replace_types[1] = viennamath::diff(v, r) * dr_dx;
      
      //substitute expressions in lhs and rhs:
      new_lhs = viennamath::substitute(search_types, replace_types, weak_form.lhs());
      new_rhs = viennamath::substitute(search_types, replace_types, weak_form.rhs());
      
      viennamath::rt_manipulation_wrapper<InterfaceType> jacobian_manipulator( new detail::jacobian_adder<CellType, InterfaceType, ReferenceCellTag>() );
      Expression new_lhs_with_jacobian(new_lhs.get()->recursive_manipulation(jacobian_manipulator));
      Expression new_rhs_with_jacobian(new_rhs.get()->recursive_manipulation(jacobian_manipulator));
      
      return EquationType(new_lhs_with_jacobian, new_rhs_with_jacobian);
    }
    
    
    //////////////////////////         2d          //////////////////////////////
    
    
    //
    // Triangles and Quadrilaterals:
    //
    /** @brief Transforms a weak formulation given in physical space to the reference cell. Also transforms derivatives where appropriate.  Implementation for the 2d-case.*/
    template <typename CellType, typename EquationType, typename PDESystemType, typename ReferenceCellTag>
    EquationType transform_to_reference_cell_2d(EquationType const & weak_form,
                                                PDESystemType const & pde_system,
                                                ReferenceCellTag)
    {
      typedef typename EquationType::interface_type             InterfaceType;
      typedef typename InterfaceType::numeric_type              numeric_type;
      typedef viennamath::rt_function_symbol<InterfaceType>   FunctionSymbol;
      typedef viennamath::rt_expr<InterfaceType>              Expression;
      typedef viennamath::rt_unary_expr<InterfaceType>        UnaryExpression;
      typedef viennamath::rt_variable<InterfaceType>          Variable;
      typedef viennamath::op_unary<viennamath::op_partial_deriv<viennamath::default_numeric_type>, InterfaceType >   d_dx;
      
      
      
      FunctionSymbol  u(0, viennamath::unknown_tag<>());
      FunctionSymbol  v(0, viennamath::test_tag<>());
      Variable        r(0);
      Variable        s(1);
      
      //using local coordinates (r, s) and global coordinates (x, y)
      viennafem::cell_quan<CellType, InterfaceType>  dr_dx; detail::wrap_dt_dx<0,0>(dr_dx, ReferenceCellTag());
      viennafem::cell_quan<CellType, InterfaceType>  ds_dx; detail::wrap_dt_dx<1,0>(ds_dx, ReferenceCellTag());
      
      viennafem::cell_quan<CellType, InterfaceType>  dr_dy; detail::wrap_dt_dx<0,1>(dr_dy, ReferenceCellTag());
      viennafem::cell_quan<CellType, InterfaceType>  ds_dy; detail::wrap_dt_dx<1,1>(ds_dy, ReferenceCellTag());
      
      Expression new_lhs = weak_form.lhs();
      Expression new_rhs = weak_form.rhs();
      
      //std::cout << "New lhs init: " << std::endl;
      //std::cout << new_lhs << std::endl;
      
      std::vector<Expression> search_types(4);
      search_types[0] = UnaryExpression(u.clone(), new d_dx(0));
      search_types[1] = UnaryExpression(u.clone(), new d_dx(1));
      search_types[2] = UnaryExpression(v.clone(), new d_dx(0));
      search_types[3] = UnaryExpression(v.clone(), new d_dx(1));
      
      std::vector<Expression> replace_types(4);
      replace_types[0] = viennamath::diff(u, r) * dr_dx + viennamath::diff(u, s) * ds_dx;
      replace_types[1] = viennamath::diff(u, r) * dr_dy + viennamath::diff(u, s) * ds_dy;
      
      replace_types[2] = viennamath::diff(v, r) * dr_dx + viennamath::diff(v, s) * ds_dx;
      replace_types[3] = viennamath::diff(v, r) * dr_dy + viennamath::diff(v, s) * ds_dy;
      
      //substitute expressions in lhs and rhs:
      new_lhs = viennamath::substitute(search_types, replace_types, weak_form.lhs());
      new_rhs = viennamath::substitute(search_types, replace_types, weak_form.rhs());
      
      viennamath::rt_manipulation_wrapper<InterfaceType> jacobian_manipulator( new detail::jacobian_adder<CellType, InterfaceType, ReferenceCellTag>() );
      Expression new_lhs_with_jacobian(new_lhs.get()->recursive_manipulation(jacobian_manipulator));
      Expression new_rhs_with_jacobian(new_rhs.get()->recursive_manipulation(jacobian_manipulator));
      
      return EquationType(new_lhs_with_jacobian, new_rhs_with_jacobian);
    }




    //////////////////////////         3d          //////////////////////////////

    //
    // Tetrahedra and Hexahedra
    //
    /** @brief Transforms a weak formulation given in physical space to the reference cell. Also transforms derivatives where appropriate. Implementation for the 3d-case.*/
    template <typename CellType, typename EquationType, typename PDESystemType, typename ReferenceCellTag>
    EquationType transform_to_reference_cell_3d(EquationType const & weak_form,
                                                PDESystemType const & pde_system,
                                                ReferenceCellTag)
    {
      typedef typename EquationType::interface_type             InterfaceType;
      typedef viennamath::rt_function_symbol<InterfaceType>   FunctionSymbol;
      typedef viennamath::rt_expr<InterfaceType>              Expression;
      typedef viennamath::rt_unary_expr<InterfaceType>        UnaryExpression;
      typedef viennamath::rt_variable<InterfaceType>          Variable;
      typedef viennamath::op_unary<viennamath::op_partial_deriv<viennamath::default_numeric_type>, InterfaceType >   d_dx;

      
      //FunctionSymbol  u(0, viennamath::unknown_tag<>());
      //FunctionSymbol  v(0, viennamath::test_tag<>());
      
      std::size_t num_components = pde_system.unknown(0).size();
      std::cout << "Number of components: " << num_components << std::endl;
      std::vector< FunctionSymbol > u(num_components);
      for (std::size_t i=0; i<num_components; ++i)
        u[i] = FunctionSymbol(pde_system.unknown(0)[i].id(), viennamath::unknown_tag<>());

      std::vector< FunctionSymbol > v(num_components);
      for (std::size_t i=0; i<num_components; ++i)
        v[i] = FunctionSymbol(pde_system.unknown(0)[i].id(), viennamath::test_tag<>());
      
      Variable      r(0);
      Variable      s(1);
      Variable      t(2);
      
      //using local coordinates (r, s) and global coordinates (x, y)
      viennafem::cell_quan<CellType, InterfaceType>  dr_dx; detail::wrap_dt_dx<0,0>(dr_dx, ReferenceCellTag());
      viennafem::cell_quan<CellType, InterfaceType>  ds_dx; detail::wrap_dt_dx<1,0>(ds_dx, ReferenceCellTag());
      viennafem::cell_quan<CellType, InterfaceType>  dt_dx; detail::wrap_dt_dx<2,0>(dt_dx, ReferenceCellTag());
      
      viennafem::cell_quan<CellType, InterfaceType>  dr_dy; detail::wrap_dt_dx<0,1>(dr_dy, ReferenceCellTag());
      viennafem::cell_quan<CellType, InterfaceType>  ds_dy; detail::wrap_dt_dx<1,1>(ds_dy, ReferenceCellTag());
      viennafem::cell_quan<CellType, InterfaceType>  dt_dy; detail::wrap_dt_dx<2,1>(dt_dy, ReferenceCellTag());

      viennafem::cell_quan<CellType, InterfaceType>  dr_dz; detail::wrap_dt_dx<0,2>(dr_dz, ReferenceCellTag());
      viennafem::cell_quan<CellType, InterfaceType>  ds_dz; detail::wrap_dt_dx<1,2>(ds_dz, ReferenceCellTag());
      viennafem::cell_quan<CellType, InterfaceType>  dt_dz; detail::wrap_dt_dx<2,2>(dt_dz, ReferenceCellTag());
      
      Expression new_lhs = weak_form.lhs();
      Expression new_rhs = weak_form.rhs();

      
      std::vector<Expression> search_types(2 * 3 * num_components);
      std::size_t current_index = 0;
      for (std::size_t i=0; i<num_components; ++i)
      {
        search_types[current_index++] = UnaryExpression(u[i].clone(), new d_dx(0));
        search_types[current_index++] = UnaryExpression(u[i].clone(), new d_dx(1));
        search_types[current_index++] = UnaryExpression(u[i].clone(), new d_dx(2));
      
        search_types[current_index++] = UnaryExpression(v[i].clone(), new d_dx(0));
        search_types[current_index++] = UnaryExpression(v[i].clone(), new d_dx(1));
        search_types[current_index++] = UnaryExpression(v[i].clone(), new d_dx(2));
      }
      
      std::vector<Expression> replace_types(2 * 3 * num_components);
      current_index = 0;
      for (std::size_t i=0; i<num_components; ++i)
      {
        replace_types[current_index++] = viennamath::diff(u[i], r) * dr_dx + viennamath::diff(u[i], s) * ds_dx + viennamath::diff(u[i], t) * dt_dx;
        replace_types[current_index++] = viennamath::diff(u[i], r) * dr_dy + viennamath::diff(u[i], s) * ds_dy + viennamath::diff(u[i], t) * dt_dy;
        replace_types[current_index++] = viennamath::diff(u[i], r) * dr_dz + viennamath::diff(u[i], s) * ds_dz + viennamath::diff(u[i], t) * dt_dz;
        
        replace_types[current_index++] = viennamath::diff(v[i], r) * dr_dx + viennamath::diff(v[i], s) * ds_dx + viennamath::diff(v[i], t) * dt_dx;
        replace_types[current_index++] = viennamath::diff(v[i], r) * dr_dy + viennamath::diff(v[i], s) * ds_dy + viennamath::diff(v[i], t) * dt_dy;
        replace_types[current_index++] = viennamath::diff(v[i], r) * dr_dz + viennamath::diff(v[i], s) * ds_dz + viennamath::diff(v[i], t) * dt_dz;
      }
      
      //substitute expressions in lhs and rhs:
      new_lhs = viennamath::substitute(search_types, replace_types, weak_form.lhs());
      new_rhs = viennamath::substitute(search_types, replace_types, weak_form.rhs());
      
      viennamath::rt_manipulation_wrapper<InterfaceType> jacobian_manipulator( new detail::jacobian_adder<CellType,
                                                                                                          InterfaceType,
                                                                                                          ReferenceCellTag>() );
      Expression new_lhs_with_jacobian(new_lhs.get()->recursive_manipulation(jacobian_manipulator));
      Expression new_rhs_with_jacobian(new_rhs.get()->recursive_manipulation(jacobian_manipulator));
      
      return EquationType(new_lhs_with_jacobian, new_rhs_with_jacobian);
    }

  } //namespace detail
  

  ///////// Public Interface //////////////////
  
  // TODO: This can be more automatic/generative -> arbitrary dimensions!

  // 1d

  /** @brief Transforms a general weak form on any line to an integration over the unit interval */
  template <typename CellType, typename EquationType, typename PDESystemType>
  EquationType transform_to_reference_cell(EquationType const & weak_form,
                                           PDESystemType const & pde_system,
                                           viennafem::unit_interval)
  {
    return  detail::transform_to_reference_cell_1d<CellType>(weak_form, pde_system, viennafem::unit_interval());
  }

  // 2d

  /** @brief Transforms a general weak form on any quadrilateral to an integration over the unit square */
  template <typename CellType, typename EquationType, typename PDESystemType>
  EquationType transform_to_reference_cell(EquationType const & weak_form,
                                           PDESystemType const & pde_system,
                                           viennafem::unit_square)
  {
    return  detail::transform_to_reference_cell_2d<CellType>(weak_form, pde_system, viennafem::unit_square());
  }

  /** @brief Transforms a general weak form on any triangle to an integration over the unit triangle */
  template <typename CellType, typename EquationType, typename PDESystemType>
  EquationType transform_to_reference_cell(EquationType const & weak_form,
                                           PDESystemType const & pde_system,
                                           viennafem::unit_triangle)
  {
    return  detail::transform_to_reference_cell_2d<CellType>(weak_form, pde_system, viennafem::unit_triangle());
  }

  // 3d

  /** @brief Transforms a general weak form on any hexahedron to an integration over the unit cube */
  template <typename CellType, typename EquationType, typename PDESystemType>
  EquationType transform_to_reference_cell(EquationType const & weak_form,
                                           PDESystemType const & pde_system,
                                           viennafem::unit_cube)
  {
    return  detail::transform_to_reference_cell_3d<CellType>(weak_form, pde_system, viennafem::unit_cube());
  }

  /** @brief Transforms a general weak form on any tetrahedron to an integration over the unit tetrahedron */
  template <typename CellType, typename EquationType, typename PDESystemType>
  EquationType transform_to_reference_cell(EquationType const & weak_form,
                                           PDESystemType const & pde_system,
                                           viennafem::unit_tetrahedron)
  {
    return  detail::transform_to_reference_cell_3d<CellType>(weak_form, pde_system, viennafem::unit_tetrahedron());
  }


  // entry point:
  /** @brief Transforms a weak formulation given by an integration over an arbitrary cell to the reference element. */
  template <typename CellType, typename EquationType, typename PDESystemType>
  EquationType transform_to_reference_cell(EquationType const & weak_form,
                                           PDESystemType const & pde_system)
  {
    typedef typename CellType::tag        CellTag;
    
    return  transform_to_reference_cell<CellType>(weak_form, 
                                                  pde_system,
                                                  typename reference_cell_for_basis<CellTag, viennafem::lagrange_tag<1> >::type());
  }
  
  
  
  
  
  
  
  
  ///////// Basis function insertion ///////////////
  


  /** @brief Inserts vector-valued test and trial functions into the weak form.
   * 
   * Scalar-valued functions are a special case of vector-valued function (vectors of length 1) and are also handled here.
   *
   * @param weak_form     The weak form of the PDE
   * @param pde_system    Meta-information for the (system of) PDE(s) 
   * @param test_func     A vector of expressions. Typically something of the form (phi, 0, 0), (0, phi, 0) or (0, 0, phi) for a scalar-valued function phi
   * @param trial_func    A vector of expressions. Typically something of the form (phi, 0, 0), (0, phi, 0) or (0, 0, phi) for a scalar-valued function phi
   * @return The weak form with basis functions inserted for function symbols
   */
  template <typename EquationType, typename PDESystemType, typename ExpressionType>
  EquationType insert_test_and_trial_functions(EquationType weak_form,
                                               PDESystemType const & pde_system,
                                               std::vector<ExpressionType> test_func,
                                               std::vector<ExpressionType> trial_func
                                               )
  {
    typedef typename EquationType::interface_type             InterfaceType;
    typedef viennamath::rt_expr<InterfaceType>              Expression;
    typedef viennamath::rt_function_symbol<InterfaceType>   FunctionSymbol;
    typedef viennamath::rt_unary_expr<InterfaceType>        UnaryExpression;
    typedef viennamath::rt_variable<InterfaceType>          Variable;
    typedef viennamath::op_unary<viennamath::op_partial_deriv<viennamath::default_numeric_type>, InterfaceType >   d_dx;
    
    std::vector< FunctionSymbol > u(trial_func.size());
    for (std::size_t i=0; i<trial_func.size(); ++i)
      u[i] = FunctionSymbol(pde_system.unknown(0)[i].id(), viennamath::unknown_tag<>());

    std::vector< FunctionSymbol > v(test_func.size());
    for (std::size_t i=0; i<test_func.size(); ++i)
      v[i] = FunctionSymbol(pde_system.unknown(0)[i].id(), viennamath::test_tag<>());

    Variable r(0);
    Variable s(1);
    Variable t(2);
    
    std::vector<Expression> search_keys(4 * (test_func.size() + trial_func.size()));
    std::size_t current_index = 0;
    for (std::size_t i=0; i<trial_func.size(); ++i)
    {      
      search_keys[current_index++] = UnaryExpression(u[i].clone(), new d_dx(0));
      search_keys[current_index++] = UnaryExpression(u[i].clone(), new d_dx(1));
      search_keys[current_index++] = UnaryExpression(u[i].clone(), new d_dx(2));
      search_keys[current_index++] = u[i];
    }
    
    for (std::size_t i=0; i<test_func.size(); ++i)
    {      
      search_keys[current_index++] = UnaryExpression(v[i].clone(), new d_dx(0));
      search_keys[current_index++] = UnaryExpression(v[i].clone(), new d_dx(1));
      search_keys[current_index++] = UnaryExpression(v[i].clone(), new d_dx(2));
      search_keys[current_index++] = v[i];
    }
    
    
    std::vector<Expression> replacements(4 * (test_func.size() + trial_func.size()));
    current_index = 0;
    for (std::size_t i=0; i<trial_func.size(); ++i)
    {      
      replacements[current_index++] = viennamath::diff(trial_func[i], r);
      replacements[current_index++] = viennamath::diff(trial_func[i], s);
      replacements[current_index++] = viennamath::diff(trial_func[i], t);
      replacements[current_index++] = trial_func[i];
    }
    
    for (std::size_t i=0; i<test_func.size(); ++i)
    {      
      replacements[current_index++] = viennamath::diff(test_func[i], r);
      replacements[current_index++] = viennamath::diff(test_func[i], s);
      replacements[current_index++] = viennamath::diff(test_func[i], t);
      replacements[current_index++] = test_func[i];
    }
    
    ExpressionType new_lhs = viennamath::substitute(search_keys, replacements, weak_form.lhs());
    ExpressionType new_rhs = viennamath::substitute(search_keys, replacements, weak_form.rhs());
    
    //std::cout << "Old rhs: " << weak_form.rhs() << std::endl;
    //std::cout << "New rhs: " << new_rhs << std::endl;
    
    return EquationType(new_lhs, new_rhs);
  }
  


/* deprecated:
  template <typename EquationType, typename PDESystemType, typename ExpressionType>
  EquationType insert_test_and_trial_functions_vector(EquationType weak_form,
                                               PDESystemType const & pde_system,
                                               ExpressionType test_func,
                                               ExpressionType trial_func,
                                               std::size_t index_i,
                                               std::size_t index_j
                                               )
  {
    typedef typename EquationType::interface_type             InterfaceType;
    typedef viennamath::rt_expr<InterfaceType>              Expression;
    typedef viennamath::rt_function_symbol<InterfaceType>   FunctionSymbol;
    typedef viennamath::rt_unary_expr<InterfaceType>        UnaryExpression;
    typedef viennamath::rt_variable<InterfaceType>          Variable;
    typedef viennamath::op_unary<viennamath::op_partial_deriv<viennamath::default_numeric_type>, InterfaceType >   d_dx;
    
    std::vector< FunctionSymbol > u(3);
    for (std::size_t i=0; i<3; ++i)
    {
      u[i] = FunctionSymbol(pde_system.unknown(0)[i].id(), viennamath::unknown_tag<>());
      //std::cout << "ID of u at location " << i << ": " << pde_system.unknown(0)[i].id() << std::endl;
    }

    std::vector< FunctionSymbol > v(3);
    for (std::size_t i=0; i<3; ++i)
    {
      v[i] = FunctionSymbol(pde_system.unknown(0)[i].id(), viennamath::test_tag<>());
      //std::cout << "ID of v at location " << i << ": " << pde_system.unknown(0)[i].id() << std::endl;
    }

    Variable r(0);
    Variable s(1);
    Variable t(2);
    
    std::vector<Expression> search_keys(4 * (2 + 2));
    std::size_t current_index = 0;
    for (std::size_t i=0; i<3; ++i)
    {      
      if (i != index_j)
      {  
        search_keys[current_index++] = UnaryExpression(u[i].clone(), new d_dx(0));
        search_keys[current_index++] = UnaryExpression(u[i].clone(), new d_dx(1));
        search_keys[current_index++] = UnaryExpression(u[i].clone(), new d_dx(2));
        search_keys[current_index++] = u[i];
      }
    }
    
    for (std::size_t i=0; i<3; ++i)
    {      
      if (i != index_i)
      {  
        search_keys[current_index++] = UnaryExpression(v[i].clone(), new d_dx(0));
        search_keys[current_index++] = UnaryExpression(v[i].clone(), new d_dx(1));
        search_keys[current_index++] = UnaryExpression(v[i].clone(), new d_dx(2));
        search_keys[current_index++] = v[i];
      }
    }
    
    
    std::vector<Expression> replacements(4 * (2 + 2));
    current_index = 0;
    viennamath::rt_constant<double, InterfaceType> c0(0);
    for (std::size_t i=0; i<3; ++i)
    {      
      if (i != index_j)
      {
        replacements[current_index++] = c0;
        replacements[current_index++] = c0;
        replacements[current_index++] = c0;
        replacements[current_index++] = c0;
      }
    }
    
    for (std::size_t i=0; i<3; ++i)
    {      
      if (i != index_i)
      {
        replacements[current_index++] = c0;
        replacements[current_index++] = c0;
        replacements[current_index++] = c0;
        replacements[current_index++] = c0;
      }
    }
    
    ExpressionType new_lhs = viennamath::substitute(search_keys, replacements, weak_form.lhs());
    ExpressionType new_rhs = viennamath::substitute(search_keys, replacements, weak_form.rhs());
    
    return EquationType(new_lhs, new_rhs);
  }
  */

}
#endif
