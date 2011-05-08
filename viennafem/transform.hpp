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
    Variable        s_temp(10);
    
    //using local coordinates (r, s) and global coordinates (x, y)
    viennafem::cell_quan<CellType, InterfaceType>  dr_dx; dr_dx.wrap( viennafem::dt_dx_key<0,0>() );
    viennafem::cell_quan<CellType, InterfaceType>  ds_dx; ds_dx.wrap( viennafem::dt_dx_key<1,0>() );
    
    viennafem::cell_quan<CellType, InterfaceType>  dr_dy; dr_dy.wrap( viennafem::dt_dx_key<0,1>() );
    viennafem::cell_quan<CellType, InterfaceType>  ds_dy; ds_dy.wrap( viennafem::dt_dx_key<1,1>() );
    
    Expression new_lhs = weak_form.lhs();
    Expression new_rhs = weak_form.rhs();
    
    std::cout << "New lhs init: " << std::endl;
    std::cout << new_lhs << std::endl;
    
    //transform du/dx:
    new_lhs = viennamath::substitute(UnaryExpression(u.clone(), new d_dx(0)),
                                     viennamath::diff(u, r) * dr_dx + viennamath::diff(u, s_temp) * ds_dx,
                                     new_lhs);
    std::cout << "New lhs: " << std::endl;
    std::cout << new_lhs << std::endl;
    //exit(0);

    //transform du/dy:
    new_lhs = viennamath::substitute(UnaryExpression(u.clone(), new d_dx(1) ),
                                     viennamath::diff(u, r) * dr_dy + viennamath::diff(u, s) * ds_dy,
                                     new_lhs);

    
    //transform dv/dx:
    new_lhs = viennamath::substitute(UnaryExpression(v.clone(), new d_dx(0) ),
                                     viennamath::diff(v, r) * dr_dx + viennamath::diff(v, s_temp) * ds_dx,
                                     new_lhs);

    //transform dv/dy:
    new_lhs = viennamath::substitute(UnaryExpression(v.clone(), new d_dx(1) ),
                                     viennamath::diff(v, r) * dr_dy + viennamath::diff(v, s) * ds_dy,
                                     new_lhs);
    
    //replace s_temp by s;
    new_lhs = viennamath::substitute(UnaryExpression(u.clone(), new d_dx(10) ),
                                     UnaryExpression(u.clone(), new d_dx(1) ),
                                     new_lhs);

    new_lhs = viennamath::substitute(UnaryExpression(v.clone(), new d_dx(10) ),
                                     UnaryExpression(v.clone(), new d_dx(1) ),
                                     new_lhs);
    
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
    Variable      s_temp(10);
    Variable      t_temp(11);
    
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

    
    ////// Step 1: Transform partial derivatives of u
    
    //transform du/dx:
    new_lhs = viennamath::substitute(UnaryExpression(u.clone(), new d_dx(0) ),
                                     viennamath::diff(u, r) * dr_dx + viennamath::diff(u, s_temp) * ds_dx + viennamath::diff(u, t_temp) * dt_dx,
                                     new_lhs);

    //transform du/dy:
    new_lhs = viennamath::substitute(UnaryExpression(u.clone(), new d_dx(1) ),
                                     viennamath::diff(u, r) * dr_dy + viennamath::diff(u, s_temp) * ds_dy + viennamath::diff(u, t_temp) * dt_dy,
                                     new_lhs);

    //transform du/dz:
    new_lhs = viennamath::substitute(UnaryExpression(u.clone(), new d_dx(2) ),
                                     viennamath::diff(u, r) * dr_dz + viennamath::diff(u, s) * ds_dz + viennamath::diff(u, t) * dt_dz,
                                     new_lhs);

    
    ////// Step 2: Transform partial derivatives of v
    
    //transform dv/dx:
    new_lhs = viennamath::substitute(UnaryExpression(v.clone(), new d_dx(0) ),
                                     viennamath::diff(v, r) * dr_dx + viennamath::diff(v, s_temp) * ds_dx + viennamath::diff(v, t_temp) * dt_dx,
                                     new_lhs);
    
    //transform dv/dy:
    new_lhs = viennamath::substitute(UnaryExpression(v.clone(), new d_dx(1) ),
                                     viennamath::diff(v, r) * dr_dy + viennamath::diff(v, s_temp) * ds_dy + viennamath::diff(v, t_temp) * dt_dy,
                                     new_lhs);

    //transform dv/dz:
    new_lhs = viennamath::substitute(UnaryExpression(v.clone(), new d_dx(2) ),
                                     viennamath::diff(v, r) * dr_dz + viennamath::diff(v, s) * ds_dz + viennamath::diff(v, t) * dt_dz,
                                     new_lhs);
    

    
    
    ///// Step 3: Replace s_temp and t_temp with s and t again
    
    
    //replace s_temp by s;
    new_lhs = viennamath::substitute(UnaryExpression(u.clone(), new d_dx(10) ),
                                     UnaryExpression(u.clone(), new d_dx(1) ),
                                     new_lhs);

    new_lhs = viennamath::substitute(UnaryExpression(v.clone(), new d_dx(10) ),
                                     UnaryExpression(v.clone(), new d_dx(1) ),
                                     new_lhs);

    //replace t_temp by t
    new_lhs = viennamath::substitute(UnaryExpression(u.clone(), new d_dx(11) ),
                                     UnaryExpression(u.clone(), new d_dx(2) ),
                                     new_lhs);

    new_lhs = viennamath::substitute(UnaryExpression(v.clone(), new d_dx(11) ),
                                     UnaryExpression(v.clone(), new d_dx(2) ),
                                     new_lhs);
    
    return EquationType(new_lhs, new_rhs);
  }

  template <typename CellType, typename EquationType>
  EquationType transform_to_reference_cell(EquationType const & weak_form)
  {
    return  transform_to_reference_cell<CellType>(weak_form, typename CellType::element_tag());
  }
  





  template <typename ExpressionType, typename EquationType>
  EquationType insert_test_and_trial_functions(ExpressionType test_func,
                                               ExpressionType trial_func,
                                               EquationType weak_form)
  {
    typedef typename EquationType::interface_type             InterfaceType;
    typedef viennamath::function_symbol<InterfaceType>   FunctionSymbol;
    typedef viennamath::unary_expr<InterfaceType>        UnaryExpression;
    typedef viennamath::variable<InterfaceType>          Variable;
    typedef viennamath::op_unary<viennamath::op_partial_deriv<viennamath::default_numeric_type>, InterfaceType >   dt_dx;
    
    
    FunctionSymbol  u(0, viennamath::unknown_tag<>());
    FunctionSymbol  v(0, viennamath::test_tag<>());
    Variable      r(0);
    Variable      s(1);
    Variable      t(2);
    ExpressionType new_lhs = weak_form.lhs();
    ExpressionType new_rhs = weak_form.rhs();
    
    //step 1: substitute derivatives: 
    new_lhs = viennamath::substitute(UnaryExpression(u.clone(), new dt_dx(0)),
                                     viennamath::diff(trial_func, r),
                                     new_lhs);
    new_lhs = viennamath::substitute(UnaryExpression(v.clone(), new dt_dx(0)),
                                     viennamath::diff(test_func, r),
                                     new_lhs);
    
    new_lhs = viennamath::substitute(UnaryExpression(u.clone(), new dt_dx(1)),
                                     viennamath::diff(trial_func, s),
                                     new_lhs);
    new_lhs = viennamath::substitute(UnaryExpression(v.clone(), new dt_dx(1)),
                                     viennamath::diff(test_func, s),
                                     new_lhs);
    
    new_lhs = viennamath::substitute(UnaryExpression(u.clone(), new dt_dx(2)),
                                     viennamath::diff(trial_func, t),
                                     new_lhs);
    new_lhs = viennamath::substitute(UnaryExpression(v.clone(), new dt_dx(2)),
                                     viennamath::diff(test_func, t),
                                     new_lhs);
    
    //step 2: substitute non-derivaties:
    //TODO  (not necessary for Poisson equation)
    
    new_rhs = viennamath::substitute(ExpressionType(v.clone()),
                                     test_func,
                                     new_rhs);
    
    
    //step 3: return result:
    return EquationType(new_lhs, new_rhs);
  }
  

}
#endif
