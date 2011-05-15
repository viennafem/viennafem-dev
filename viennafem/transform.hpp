/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_TRANSFORM_HPP
#define VIENNAFEM_TRANSFORM_HPP

#include "viennafem/forwards.h"
#include "viennafem/cell_quan.hpp"

#include "viennamath/expression.hpp"
#include "viennamath/manipulation/substitute.hpp"
#include "viennamath/manipulation/diff.hpp"

#include "viennagrid/topology/triangle.hpp"
#include "viennagrid/topology/tetrahedron.hpp"

namespace viennafem
{

  
  template <typename CellType, typename EquationType>
  EquationType transform_to_reference_cell(EquationType const & weak_form, viennagrid::triangle_tag)
  {
    typedef typename EquationType::interface_type             InterfaceType;
    typedef typename InterfaceType::numeric_type              numeric_type;
    typedef viennamath::function_symbol<InterfaceType>   FunctionSymbol;
    typedef viennamath::expr<InterfaceType>              Expression;
    typedef viennamath::unary_expr<InterfaceType>        UnaryExpression;
    typedef viennamath::variable<InterfaceType>          Variable;
    typedef viennamath::op_unary<viennamath::op_partial_deriv<viennamath::default_numeric_type>, InterfaceType >   d_dx;
    
    
    
    FunctionSymbol  u(0, viennamath::unknown_tag<>());
    FunctionSymbol  v(0, viennamath::test_tag<>());
    Variable        r(0);
    Variable        s(1);
    
    //using local coordinates (r, s) and global coordinates (x, y)
    viennafem::cell_quan<CellType, InterfaceType>  dr_dx; dr_dx.wrap( viennafem::dt_dx_key<0,0>() );
    viennafem::cell_quan<CellType, InterfaceType>  ds_dx; ds_dx.wrap( viennafem::dt_dx_key<1,0>() );
    
    viennafem::cell_quan<CellType, InterfaceType>  dr_dy; dr_dy.wrap( viennafem::dt_dx_key<0,1>() );
    viennafem::cell_quan<CellType, InterfaceType>  ds_dy; ds_dy.wrap( viennafem::dt_dx_key<1,1>() );
    
    Expression new_lhs = weak_form.lhs();
    Expression new_rhs = weak_form.rhs();
    
    std::cout << "New lhs init: " << std::endl;
    std::cout << new_lhs << std::endl;
    
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
    
    //substitute expressions in lhs:
    new_lhs = viennamath::substitute(search_types, replace_types, weak_form.lhs());
    
    return EquationType(new_lhs, new_rhs);
  }



  template <typename CellType, typename EquationType>
  EquationType transform_to_reference_cell(EquationType const & weak_form, viennagrid::tetrahedron_tag)
  {
    typedef typename EquationType::interface_type             InterfaceType;
    typedef viennamath::function_symbol<InterfaceType>   FunctionSymbol;
    typedef viennamath::expr<InterfaceType>              Expression;
    typedef viennamath::unary_expr<InterfaceType>        UnaryExpression;
    typedef viennamath::variable<InterfaceType>          Variable;
    typedef viennamath::op_unary<viennamath::op_partial_deriv<viennamath::default_numeric_type>, InterfaceType >   d_dx;

    
    FunctionSymbol  u(0, viennamath::unknown_tag<>());
    FunctionSymbol  v(0, viennamath::test_tag<>());
    Variable      r(0);
    Variable      s(1);
    Variable      t(2);
    
    //using local coordinates (r, s) and global coordinates (x, y)
    viennafem::cell_quan<CellType, InterfaceType>  dr_dx; dr_dx.wrap( viennafem::dt_dx_key<0,0>() );
    viennafem::cell_quan<CellType, InterfaceType>  ds_dx; ds_dx.wrap( viennafem::dt_dx_key<1,0>() );
    viennafem::cell_quan<CellType, InterfaceType>  dt_dx; dt_dx.wrap( viennafem::dt_dx_key<2,0>() );
    
    viennafem::cell_quan<CellType, InterfaceType>  dr_dy; dr_dy.wrap( viennafem::dt_dx_key<0,1>() );
    viennafem::cell_quan<CellType, InterfaceType>  ds_dy; ds_dy.wrap( viennafem::dt_dx_key<1,1>() );
    viennafem::cell_quan<CellType, InterfaceType>  dt_dy; dt_dy.wrap( viennafem::dt_dx_key<2,1>() );

    viennafem::cell_quan<CellType, InterfaceType>  dr_dz; dr_dz.wrap( viennafem::dt_dx_key<0,2>() );
    viennafem::cell_quan<CellType, InterfaceType>  ds_dz; ds_dz.wrap( viennafem::dt_dx_key<1,2>() );
    viennafem::cell_quan<CellType, InterfaceType>  dt_dz; dt_dz.wrap( viennafem::dt_dx_key<2,2>() );
    
    Expression new_lhs = weak_form.lhs();
    Expression new_rhs = weak_form.rhs();

    
    std::vector<Expression> search_types(6);
    search_types[0] = UnaryExpression(u.clone(), new d_dx(0));
    search_types[1] = UnaryExpression(u.clone(), new d_dx(1));
    search_types[2] = UnaryExpression(u.clone(), new d_dx(2));
    
    search_types[3] = UnaryExpression(v.clone(), new d_dx(0));
    search_types[4] = UnaryExpression(v.clone(), new d_dx(1));
    search_types[5] = UnaryExpression(v.clone(), new d_dx(2));
    
    std::vector<Expression> replace_types(6);
    replace_types[0] = viennamath::diff(u, r) * dr_dx + viennamath::diff(u, s) * ds_dx + viennamath::diff(u, t) * dt_dx,
    replace_types[1] = viennamath::diff(u, r) * dr_dy + viennamath::diff(u, s) * ds_dy + viennamath::diff(u, t) * dt_dy,
    replace_types[2] = viennamath::diff(u, r) * dr_dz + viennamath::diff(u, s) * ds_dz + viennamath::diff(u, t) * dt_dz,
    
    replace_types[3] = viennamath::diff(v, r) * dr_dx + viennamath::diff(v, s) * ds_dx + viennamath::diff(v, t) * dt_dx,
    replace_types[4] = viennamath::diff(v, r) * dr_dy + viennamath::diff(v, s) * ds_dy + viennamath::diff(v, t) * dt_dy,
    replace_types[5] = viennamath::diff(v, r) * dr_dz + viennamath::diff(v, s) * ds_dz + viennamath::diff(v, t) * dt_dz,
    
    //substitute expressions in lhs:
    new_lhs = viennamath::substitute(search_types, replace_types, weak_form.lhs());
    
    return EquationType(new_lhs, new_rhs);
  }

  template <typename CellType, typename EquationType>
  EquationType transform_to_reference_cell(EquationType const & weak_form)
  {
    return  transform_to_reference_cell<CellType>(weak_form, typename CellType::element_tag());
  }
  



  template <typename EquationType, typename PDESystemType, typename ExpressionType>
  EquationType insert_test_and_trial_functions(EquationType weak_form,
                                               PDESystemType const & pde_system,
                                               std::vector<ExpressionType> test_func,
                                               std::vector<ExpressionType> trial_func
                                               )
  {
    typedef typename EquationType::interface_type             InterfaceType;
    typedef viennamath::expr<InterfaceType>              Expression;
    typedef viennamath::function_symbol<InterfaceType>   FunctionSymbol;
    typedef viennamath::unary_expr<InterfaceType>        UnaryExpression;
    typedef viennamath::variable<InterfaceType>          Variable;
    typedef viennamath::op_unary<viennamath::op_partial_deriv<viennamath::default_numeric_type>, InterfaceType >   d_dx;
    
    std::vector< FunctionSymbol > u(trial_func.size());
    for (size_t i=0; i<trial_func.size(); ++i)
      u[i] = FunctionSymbol(pde_system.unknown(0)[i].id(), viennamath::unknown_tag<>());

    std::vector< FunctionSymbol > v(test_func.size());
    for (size_t i=0; i<test_func.size(); ++i)
      v[i] = FunctionSymbol(pde_system.unknown(0)[i].id(), viennamath::test_tag<>());

    Variable r(0);
    Variable s(1);
    Variable t(2);
    
    std::vector<Expression> search_keys(4 * (test_func.size() + trial_func.size()));
    size_t current_index = 0;
    for (size_t i=0; i<trial_func.size(); ++i)
    {      
      search_keys[current_index++] = UnaryExpression(u[i].clone(), new d_dx(0));
      search_keys[current_index++] = UnaryExpression(u[i].clone(), new d_dx(1));
      search_keys[current_index++] = UnaryExpression(u[i].clone(), new d_dx(2));
      search_keys[current_index++] = u[i];
    }
    
    for (size_t i=0; i<test_func.size(); ++i)
    {      
      search_keys[current_index++] = UnaryExpression(v[i].clone(), new d_dx(0));
      search_keys[current_index++] = UnaryExpression(v[i].clone(), new d_dx(1));
      search_keys[current_index++] = UnaryExpression(v[i].clone(), new d_dx(2));
      search_keys[current_index++] = v[i];
    }
    
    
    std::vector<Expression> replacements(4 * (test_func.size() + trial_func.size()));
    current_index = 0;
    for (size_t i=0; i<trial_func.size(); ++i)
    {      
      replacements[current_index++] = viennamath::diff(trial_func[i], r);
      replacements[current_index++] = viennamath::diff(trial_func[i], s);
      replacements[current_index++] = viennamath::diff(trial_func[i], t);
      replacements[current_index++] = trial_func[i];
    }
    
    for (size_t i=0; i<test_func.size(); ++i)
    {      
      replacements[current_index++] = viennamath::diff(test_func[i], r);
      replacements[current_index++] = viennamath::diff(test_func[i], s);
      replacements[current_index++] = viennamath::diff(test_func[i], t);
      replacements[current_index++] = test_func[i];
    }
    
    ExpressionType new_lhs = viennamath::substitute(search_keys, replacements, weak_form.lhs());
    ExpressionType new_rhs = viennamath::substitute(search_keys, replacements, weak_form.rhs());
    
    return EquationType(new_lhs, new_rhs);
  }
  

}
#endif
