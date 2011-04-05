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

#include "viennagrid/celltags.hpp"

namespace viennafem
{

  
  template <typename CellType>
  viennamath::equation<> transform_to_reference_cell(viennamath::equation<> const & weak_form, viennagrid::triangle_tag)
  {
    typedef viennamath::op_unary<viennamath::op_partial_deriv<viennamath::default_numeric_type> >   d_dx;
    
    viennamath::function_symbol<>  u(0, viennamath::unknown_tag<>());
    viennamath::function_symbol<>  v(0, viennamath::test_tag<>());
    viennamath::variable<>      r(0);
    viennamath::variable<>      s(1);
    viennamath::variable<>     s_temp(10);
    
    //using local coordinates (r, s) and global coordinates (x, y)
    viennafem::cell_quan<CellType, viennafem::dt_dx_key<0,0> >  dr_dx;
    viennafem::cell_quan<CellType, viennafem::dt_dx_key<1,0> >  ds_dx;
    viennafem::cell_quan<CellType, viennafem::dt_dx_key<0,1> >  dr_dy;
    viennafem::cell_quan<CellType, viennafem::dt_dx_key<1,1> >  ds_dy;
    
    viennamath::expr<> new_lhs = weak_form.lhs();
    viennamath::expr<> new_rhs = weak_form.rhs();
    
    std::cout << "New lhs init: " << std::endl;
    std::cout << new_lhs << std::endl;
    
    //transform du/dx:
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(u.clone(), new d_dx(0)),
                                     viennamath::diff(u, r) * dr_dx + viennamath::diff(u, s_temp) * ds_dx,
                                     new_lhs);
    std::cout << "New lhs: " << std::endl;
    std::cout << new_lhs << std::endl;
    //exit(0);

    //transform du/dy:
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(u.clone(), new d_dx(1) ),
                                     viennamath::diff(u, r) * dr_dy + viennamath::diff(u, s) * ds_dy,
                                     new_lhs);

    
    //transform dv/dx:
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(v.clone(), new d_dx(0) ),
                                     viennamath::diff(v, r) * dr_dx + viennamath::diff(v, s_temp) * ds_dx,
                                     new_lhs);

    //transform dv/dy:
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(v.clone(), new d_dx(1) ),
                                     viennamath::diff(v, r) * dr_dy + viennamath::diff(v, s) * ds_dy,
                                     new_lhs);
    
    //replace s_temp by s;
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(u.clone(), new d_dx(10) ),
                                     viennamath::unary_expr<>(u.clone(), new d_dx(1) ),
                                     new_lhs);

    new_lhs = viennamath::substitute(viennamath::unary_expr<>(v.clone(), new d_dx(10) ),
                                     viennamath::unary_expr<>(v.clone(), new d_dx(1) ),
                                     new_lhs);
    
    return viennamath::equation<>(new_lhs, new_rhs);
  }



  template <typename CellType>
  viennamath::equation<> transform_to_reference_cell(viennamath::equation<> const & weak_form, viennagrid::tetrahedron_tag)
  {
    typedef viennamath::op_unary<viennamath::op_partial_deriv<viennamath::default_numeric_type> >   d_dx;
    
    viennamath::function_symbol<>  u(0, viennamath::unknown_tag<>());
    viennamath::function_symbol<>  v(0, viennamath::test_tag<>());
    viennamath::variable<>      r(0);
    viennamath::variable<>      s(1);
    viennamath::variable<>      t(2);
    viennamath::variable<>      s_temp(10);
    viennamath::variable<>      t_temp(11);
    
    //using local coordinates (r, s) and global coordinates (x, y)
    viennafem::cell_quan<CellType, viennafem::dt_dx_key<0,0> >  dr_dx;
    viennafem::cell_quan<CellType, viennafem::dt_dx_key<1,0> >  ds_dx;
    viennafem::cell_quan<CellType, viennafem::dt_dx_key<2,0> >  dt_dx;
    
    viennafem::cell_quan<CellType, viennafem::dt_dx_key<0,1> >  dr_dy;
    viennafem::cell_quan<CellType, viennafem::dt_dx_key<1,1> >  ds_dy;
    viennafem::cell_quan<CellType, viennafem::dt_dx_key<2,1> >  dt_dy;

    viennafem::cell_quan<CellType, viennafem::dt_dx_key<0,2> >  dr_dz;
    viennafem::cell_quan<CellType, viennafem::dt_dx_key<1,2> >  ds_dz;
    viennafem::cell_quan<CellType, viennafem::dt_dx_key<2,2> >  dt_dz;
    
    viennamath::expr<> new_lhs = weak_form.lhs();
    viennamath::expr<> new_rhs = weak_form.rhs();

    
    ////// Step 1: Transform partial derivatives of u
    
    //transform du/dx:
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(u.clone(), new d_dx(0) ),
                                     viennamath::diff(u, r) * dr_dx + viennamath::diff(u, s_temp) * ds_dx + viennamath::diff(u, t_temp) * dt_dx,
                                     new_lhs);

    //transform du/dy:
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(u.clone(), new d_dx(1) ),
                                     viennamath::diff(u, r) * dr_dy + viennamath::diff(u, s_temp) * ds_dy + viennamath::diff(u, t_temp) * dt_dy,
                                     new_lhs);

    //transform du/dz:
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(u.clone(), new d_dx(2) ),
                                     viennamath::diff(u, r) * dr_dz + viennamath::diff(u, s) * ds_dz + viennamath::diff(u, t) * dt_dz,
                                     new_lhs);

    
    ////// Step 2: Transform partial derivatives of v
    
    //transform dv/dx:
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(v.clone(), new d_dx(0) ),
                                     viennamath::diff(v, r) * dr_dx + viennamath::diff(v, s_temp) * ds_dx + viennamath::diff(v, t_temp) * dt_dx,
                                     new_lhs);
    
    //transform dv/dy:
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(v.clone(), new d_dx(1) ),
                                     viennamath::diff(v, r) * dr_dy + viennamath::diff(v, s_temp) * ds_dy + viennamath::diff(v, t_temp) * dt_dy,
                                     new_lhs);

    //transform dv/dz:
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(v.clone(), new d_dx(2) ),
                                     viennamath::diff(v, r) * dr_dz + viennamath::diff(v, s) * ds_dz + viennamath::diff(v, t) * dt_dz,
                                     new_lhs);
    

    
    
    ///// Step 3: Replace s_temp and t_temp with s and t again
    
    
    //replace s_temp by s;
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(u.clone(), new d_dx(10) ),
                                     viennamath::unary_expr<>(u.clone(), new d_dx(1) ),
                                     new_lhs);

    new_lhs = viennamath::substitute(viennamath::unary_expr<>(v.clone(), new d_dx(10) ),
                                     viennamath::unary_expr<>(v.clone(), new d_dx(1) ),
                                     new_lhs);

    //replace t_temp by t
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(u.clone(), new d_dx(11) ),
                                     viennamath::unary_expr<>(u.clone(), new d_dx(2) ),
                                     new_lhs);

    new_lhs = viennamath::substitute(viennamath::unary_expr<>(v.clone(), new d_dx(11) ),
                                     viennamath::unary_expr<>(v.clone(), new d_dx(2) ),
                                     new_lhs);
    
    return viennamath::equation<>(new_lhs, new_rhs);
  }

  template <typename CellType>
  viennamath::equation<> transform_to_reference_cell(viennamath::equation<> const & weak_form)
  {
    return  transform_to_reference_cell<CellType>(weak_form, typename CellType::element_tag());
  }
  





  
  viennamath::equation<> insert_test_and_trial_functions(viennamath::expr<> test_func,
                                                       viennamath::expr<> trial_func,
                                                       viennamath::equation<> weak_form)
  {
    typedef viennamath::op_unary<viennamath::op_partial_deriv<viennamath::default_numeric_type> >   dt_dx;
    
    viennamath::function_symbol<>  u(0, viennamath::unknown_tag<>());
    viennamath::function_symbol<>  v(0, viennamath::test_tag<>());
    viennamath::variable<>      r(0);
    viennamath::variable<>      s(1);
    viennamath::variable<>      t(2);
    viennamath::expr<> new_lhs = weak_form.lhs();
    viennamath::expr<> new_rhs = weak_form.rhs();
    
    //step 1: substitute derivatives: 
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(u.clone(), new dt_dx(0)),
                                     viennamath::diff(trial_func, r),
                                     new_lhs);
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(v.clone(), new dt_dx(0)),
                                     viennamath::diff(test_func, r),
                                     new_lhs);
    
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(u.clone(), new dt_dx(1)),
                                     viennamath::diff(trial_func, s),
                                     new_lhs);
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(v.clone(), new dt_dx(1)),
                                     viennamath::diff(test_func, s),
                                     new_lhs);
    
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(u.clone(), new dt_dx(2)),
                                     viennamath::diff(trial_func, t),
                                     new_lhs);
    new_lhs = viennamath::substitute(viennamath::unary_expr<>(v.clone(), new dt_dx(2)),
                                     viennamath::diff(test_func, t),
                                     new_lhs);
    
    //step 2: substitute non-derivaties:
    //TODO  (not necessary for Poisson equation)
    
    new_rhs = viennamath::substitute(viennamath::expr<>(v.clone()),
                                     test_func,
                                     new_rhs);
    
    
    //step 3: return result:
    return viennamath::equation<>(new_lhs, new_rhs);
  }
  

}
#endif
