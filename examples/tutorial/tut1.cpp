/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

// include necessary system headers
#include <iostream>


//#include "viennafem/afftrans.hpp"
//#include "viennafem/dtdx_tetrahedron.h"
#include "viennafem/dtdx_triangle.h"
#include "viennafem/typelist.h"
#include "viennafem/forwards.h"
#include "viennafem/BFStock.hpp"
//#include "viennafem/assembling.hpp"
//#include "viennafem/mapping.hpp"
#include "viennafem/cell_quan.hpp"

#include "viennagrid/domain.hpp"
#include "viennagrid/io/sgf_reader.hpp"
#include "viennagrid/io/vtk_writer.hpp"

#include "viennadata/interface.hpp"

#include "viennamath/expression.hpp"
#include "viennamath/substitute.hpp"
#include "viennamath/diff.hpp"
#include "viennamath/equation.hpp"
#include "viennamath/weak_form.hpp"
#include "viennamath/apply_coordinate_system.hpp"

struct TriangleConfig
{
  typedef double                                  numeric_type;
  typedef viennagrid::two_dimensions_tag          dimension_tag;
  typedef viennagrid::triangle_tag                cell_tag;

  //multigrid:
  //typedef viennagrid::full_multigrid_tag                       multigrid_tag;
  typedef viennagrid::no_multigrid_tag             multigrid_tag;
};


// define a key and configure viennadata to use a type-based dispatch:
struct example_key {};

std::ostream & operator<<(std::ostream & stream, example_key const & dummy)
{
  stream << "example_key";
  return stream;
}

namespace viennadata
{
  // tell ViennaData to use QuickKey:
  template <>
  struct dispatch_traits<example_key>
  {
    typedef type_key_dispatch_tag    tag;
  };
}


namespace viennafem
{
  template <typename CellType>
  viennamath::equation transform_to_reference_cell(viennamath::equation const & weak_form)
  {
    viennamath::unknown_func<0>  u;
    viennamath::test_func<0>     v;
    viennamath::variable<0>      r;
    viennamath::variable<1>      s;
    viennamath::variable<10>     s_temp;
    
    //using local coordinates (r, s) and global coordinates (x, y)
    viennafem::cell_quan<CellType, viennafem::dt_dx_key<0,0> >  dr_dx;
    viennafem::cell_quan<CellType, viennafem::dt_dx_key<1,0> >  ds_dx;
    viennafem::cell_quan<CellType, viennafem::dt_dx_key<0,1> >  dr_dy;
    viennafem::cell_quan<CellType, viennafem::dt_dx_key<1,1> >  ds_dy;
    
    viennamath::expr new_lhs = weak_form.lhs();
    viennamath::expr new_rhs = weak_form.rhs();
    
    //transform du/dx:
    new_lhs = viennamath::substitute(viennamath::unary_expr(u.clone(),
                                                            new viennamath::op_unary<viennamath::op_partial_deriv<0> >()
                                                           ),
                                     viennamath::diff(u, r) * dr_dx + viennamath::diff(u, s_temp) * ds_dx,
                                     new_lhs);

    //transform du/dy:
    new_lhs = viennamath::substitute(viennamath::unary_expr(u.clone(),
                                                            new viennamath::op_unary<viennamath::op_partial_deriv<1> >()
                                                           ),
                                     viennamath::diff(u, r) * dr_dy + viennamath::diff(u, s) * ds_dy,
                                     new_lhs);

    new_lhs = viennamath::substitute(viennamath::unary_expr(v.clone(),
                                                            new viennamath::op_unary<viennamath::op_partial_deriv<0> >()
                                                           ),
                                     viennamath::diff(v, r) * dr_dx + viennamath::diff(v, s_temp) * ds_dx,
                                     new_lhs);

    //transform du/dy:
    new_lhs = viennamath::substitute(viennamath::unary_expr(v.clone(),
                                                            new viennamath::op_unary<viennamath::op_partial_deriv<1> >()
                                                           ),
                                     viennamath::diff(v, r) * dr_dy + viennamath::diff(v, s) * ds_dy,
                                     new_lhs);
    
    //replace s_temp by s;
    new_lhs = viennamath::substitute(viennamath::unary_expr(u.clone(),
                                                            new viennamath::op_unary<viennamath::op_partial_deriv<10> >()
                                                           ),
                                     viennamath::unary_expr(u.clone(),
                                                            new viennamath::op_unary<viennamath::op_partial_deriv<1> >()
                                                           ),
                                     new_lhs);

    new_lhs = viennamath::substitute(viennamath::unary_expr(v.clone(),
                                                            new viennamath::op_unary<viennamath::op_partial_deriv<10> >()
                                                           ),
                                     viennamath::unary_expr(v.clone(),
                                                            new viennamath::op_unary<viennamath::op_partial_deriv<1> >()
                                                           ),
                                     new_lhs);
    
    return viennamath::equation(new_lhs, new_rhs);
  }
  
  template <typename CellType>
  viennamath::equation update_cell_quantities(CellType const & cell,
                                              viennamath::equation const & weak_form)
  {
    //step 1: update det_dF_dt
    viennamath::expr new_lhs =
       viennamath::substitute(viennamath::expr(new viennafem::cell_quan<CellType,
                                                                        viennafem::det_dF_dt_key>()),
                              viennamath::expr(new viennafem::cell_quan<CellType,
                                                                        viennafem::det_dF_dt_key>(cell)),
                              weak_form.lhs());
                              
    //step 2: update dt_dx<j, i>       
    new_lhs = viennamath::substitute(viennamath::expr(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<0, 0> >()),
                                     viennamath::expr(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<0, 0> >(cell)),
                                     new_lhs);

    new_lhs = viennamath::substitute(viennamath::expr(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<0, 1> >()),
                                     viennamath::expr(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<0, 1> >(cell)),
                                     new_lhs);

    new_lhs = viennamath::substitute(viennamath::expr(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<1, 0> >()),
                                     viennamath::expr(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<1, 0> >(cell)),
                                     new_lhs);

    new_lhs = viennamath::substitute(viennamath::expr(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<1, 1> >()),
                                     viennamath::expr(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<1, 1> >(cell)),
                                     new_lhs);
    
    return viennamath::equation(new_lhs,
                                weak_form.rhs());
  }
  
  
  viennamath::equation insert_test_and_trial_functions(viennamath::expr test_func,
                                                       viennamath::expr trial_func,
                                                       viennamath::equation weak_form)
  {
    viennamath::unknown_func<0>  u;
    viennamath::test_func<0>     v;
    viennamath::expr new_lhs = weak_form.lhs();
    viennamath::expr new_rhs = weak_form.rhs();
    
    //step 1: substitute derivatives: 
    new_lhs = viennamath::substitute(viennamath::unary_expr(u.clone(), new viennamath::op_unary<viennamath::op_partial_deriv<0> >()),
                                     viennamath::diff(trial_func, viennamath::variable<0>()),
                                     new_lhs);
    new_lhs = viennamath::substitute(viennamath::unary_expr(v.clone(), new viennamath::op_unary<viennamath::op_partial_deriv<0> >()),
                                     viennamath::diff(test_func, viennamath::variable<0>()),
                                     new_lhs);
    
    new_lhs = viennamath::substitute(viennamath::unary_expr(u.clone(), new viennamath::op_unary<viennamath::op_partial_deriv<1> >()),
                                     viennamath::diff(trial_func, viennamath::variable<1>()),
                                     new_lhs);
    new_lhs = viennamath::substitute(viennamath::unary_expr(v.clone(), new viennamath::op_unary<viennamath::op_partial_deriv<1> >()),
                                     viennamath::diff(test_func, viennamath::variable<1>()),
                                     new_lhs);
    
    //step 2: substitute non-derivaties:
    //TODO  (not necessary for Poisson equation)
    
    new_rhs = viennamath::substitute(viennamath::expr(v.clone()),
                                     test_func,
                                     new_rhs);
    
    
    //step 3: return result:
    return viennamath::equation(new_lhs, new_rhs);
  }
  
  //evaluates the weak form on a triangle using 1-point-rule
  //TODO generalize!!
  numeric_type eval_element_matrix_entry(viennamath::equation const & weak_form_substituted)
  {
    std::vector<numeric_type> p(2);
    p[0] = 1.0/3.0;
    p[1] = 1.0/3.0;
    
    return 0.5 * weak_form_substituted.lhs().get()->lhs()->eval(p); //TODO: this is pretty ugly...
  }
}



//does something:
template <typename DomainType>
void do_something(DomainType & domain)
{
  typedef typename DomainType::config_type              Config;
  typedef typename Config::cell_tag                     CellTag;
  
  typedef typename viennagrid::result_of::point_type<Config>::type                            PointType;
  typedef typename viennagrid::result_of::ncell_type<Config, CellTag::topology_level>::type   CellType;
  
  typedef typename viennagrid::result_of::ncell_container<DomainType, CellTag::topology_level>::type    CellContainer;
  typedef typename viennagrid::result_of::iterator<CellContainer>::type      CellIterator;

  typedef typename viennagrid::result_of::ncell_container<CellType, 0>::type    VertexOnCellContainer;
  typedef typename viennagrid::result_of::iterator<VertexOnCellContainer>::type      VertexOnCellIterator;
  
  viennamath::unknown_func<0> u;   //an unknown function
  
  viennamath::equation strong_form = viennamath::make_equation( viennamath::laplace(u), 1);
  viennamath::equation weak_form_general = viennamath::weak_form(strong_form);  
  viennamath::equation weak_form = viennamath::apply_coordinate_system(viennamath::cartesian<2>(),
                                                                       weak_form_general);
  
  std::cout << std::endl;
  std::cout << "Using weak form " << weak_form << std::endl;
  std::cout << std::endl;

  std::cout << "*" << std::endl;
  std::cout << "* Phase 1: Write dt_dx coefficients" << std::endl;
  std::cout << "*" << std::endl;
  //fill with cell quantities 
  CellContainer cells = viennagrid::ncells<CellTag::topology_level>(domain);
  for (CellIterator cell_iter = cells.begin();
       cell_iter != cells.end();
       ++cell_iter)
  {
    cell_iter->print_short();
    //viennadata::access<example_key, double>()(*cell_iter) = i; 
    viennafem::dt_dx_handler<CellTag>::apply(*cell_iter);
  }
  
  
  std::cout << "*" << std::endl;
  std::cout << "* Phase 2: Transform to reference element" << std::endl;
  std::cout << "*" << std::endl;
  
  viennamath::equation transformed_weak_form = viennafem::transform_to_reference_cell<CellType>(weak_form);
  
  std::cout << "Transformed weak form:" << std::endl;
  std::cout << transformed_weak_form << std::endl;
  std::cout << std::endl;

  std::cout << "*" << std::endl;
  std::cout << "* Phase 3: Assemble local element matrix" << std::endl;
  std::cout << "*" << std::endl;
  
  //transfer cell quantities to vtk:
  for (CellIterator cell_iter = cells.begin();
       cell_iter != cells.end();
       ++cell_iter)
  {
    //get basis:
    std::cout << "Getting basis..." << std::endl;
    std::vector<viennamath::expr> trial_functions = viennafem::get_basisfunctions(CellTag());
    std::vector<viennamath::expr> test_functions = viennafem::get_basisfunctions(CellTag());
    
    //set up element matrix:
    std::cout << "Creating element matrix..." << std::endl;
    std::vector<std::vector< viennafem::numeric_type > >  element_matrix(3);
    element_matrix[0].resize(3);
    element_matrix[1].resize(3);
    element_matrix[2].resize(3);
    
    //update cell_quantities:
    std::cout << "Updating cell quantities..." << std::endl;
    viennamath::equation cell_expr = viennafem::update_cell_quantities(*cell_iter, 
                                                                       transformed_weak_form);
    
    std::cout << "New cell_expr: " << cell_expr << std::endl;

    viennafem::cell_quan<CellType, viennafem::det_dF_dt_key> det_dF_dt(*cell_iter);
    
    //fill element_matrix:
    std::cout << "Filling element matrix..." << std::endl;
    for (size_t i = 0; i<test_functions.size(); ++i)
    {
      for (size_t j=0; j<trial_functions.size(); ++j)
      {
        viennamath::equation temp = viennafem::insert_test_and_trial_functions(test_functions[i],
                                                                               trial_functions[j],
                                                                               cell_expr);
        element_matrix[i][j] = viennafem::eval_element_matrix_entry(temp) * det_dF_dt.eval(1.0); 
      }
    }
    
    //print element matrix:
    for (size_t i = 0; i<test_functions.size(); ++i)
    {
      for (size_t j=0; j<trial_functions.size(); ++j)
      {
        std::cout << element_matrix[i][j] << " ";
      }
      std::cout << std::endl;
    }
    
    //write back to global matrix: TODO
  }


  viennagrid::io::Vtk_writer<DomainType> my_vtk_writer;
  my_vtk_writer.writeDomain(domain, "tut1.vtu");  

}



int main()
{
  viennagrid::domain<TriangleConfig> my_domain;
  
  try{
    viennagrid::io::sgf_reader my_sgf_reader;
    my_sgf_reader(my_domain, "../examples/data/square8.sgf");
  } catch (...){
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  
  do_something(my_domain);
  
  
  return EXIT_SUCCESS;
}
