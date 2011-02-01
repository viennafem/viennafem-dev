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


//#include "viennacl/matrix.hpp"
//#include "viennacl/vector.hpp"
//#include "viennacl/linalg/direct_solve.hpp"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>

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
struct boundary_key {};
struct mapping_key {};

namespace viennadata
{
  template <>
  struct dispatch_traits<boundary_key>
  {
    typedef type_key_dispatch_tag    tag;
  };
  
  template <>
  struct dispatch_traits<mapping_key>
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
  numeric_type eval_element_matrix_entry(viennamath::expr const & weak_form_lhs)
  {
    std::vector<numeric_type> p(2);
    p[0] = 1.0/3.0;
    p[1] = 1.0/3.0;
    
    return 0.5 * weak_form_lhs.get()->lhs()->eval(p); //TODO: this is pretty ugly...
  }
  
  //evaluates the weak form on a triangle using 1-point-rule
  //TODO generalize!!
  numeric_type eval_element_vector_entry(viennamath::expr const & weak_form_rhs)
  {
    std::vector<numeric_type> p(2);
    p[0] = 1.0/3.0;
    p[1] = 1.0/3.0;
    
    return 0.5 * weak_form_rhs.get()->lhs()->eval(p); //TODO: this is pretty ugly...
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

  typedef typename viennagrid::result_of::ncell_container<DomainType, 0>::type    VertexContainer;
  typedef typename viennagrid::result_of::iterator<VertexContainer>::type         VertexIterator;

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
  
  //setting some boundary flags:
  VertexContainer vertices = viennagrid::ncells<0>(domain);
  for (VertexIterator vit = vertices.begin();
       vit != vertices.end();
       ++vit)
  {
     if (vit->getPoint().get_x() == 0.0 || vit->getPoint().get_x() == 1.0 
       || vit->getPoint().get_y() == 0.0 || vit->getPoint().get_y() == 1.0 )
       viennadata::access<boundary_key, bool>()(*vit) = true;
     else
       viennadata::access<boundary_key, bool>()(*vit) = false;
  }
  

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
  std::cout << "* Phase 2: Create Mapping:" << std::endl;
  std::cout << "*" << std::endl;
  
  long map_index = 0;
  for (VertexIterator vit = vertices.begin();
       vit != vertices.end();
       ++vit)
  {
     if (viennadata::access<boundary_key, bool>()(*vit))
       viennadata::access<mapping_key, long>()(*vit) = -1;
     else
       viennadata::access<mapping_key, long>()(*vit) = map_index++;
  }
  std::cout << "Assigned degrees of freedom: " << map_index << std::endl;
  std::cout << std::endl;
  
  //build global system matrix and load vector:
  std::vector< std::vector<viennafem::numeric_type> > global_matrix(map_index);
  for (long i=0; i<map_index; ++i)
    global_matrix[i].resize(map_index);
  
  std::vector< viennafem::numeric_type > global_rhs(map_index);

  
  std::cout << "*" << std::endl;
  std::cout << "* Phase 3: Transform to reference element" << std::endl;
  std::cout << "*" << std::endl;
  
  viennamath::equation transformed_weak_form = viennafem::transform_to_reference_cell<CellType>(weak_form);
  
  std::cout << "Transformed weak form:" << std::endl;
  std::cout << transformed_weak_form << std::endl;
  std::cout << std::endl;

  std::cout << "*" << std::endl;
  std::cout << "* Phase 4: Assemble local element matrix" << std::endl;
  std::cout << "*" << std::endl;
  
  //transfer cell quantities to vtk:
  for (CellIterator cell_iter = cells.begin();
       cell_iter != cells.end();
       ++cell_iter)
  {
    //get basis:
    //std::cout << "Getting basis..." << std::endl;
    std::vector<viennamath::expr> trial_functions = viennafem::get_basisfunctions(CellTag());
    std::vector<viennamath::expr> test_functions = viennafem::get_basisfunctions(CellTag());
    
    //set up element matrix:
    //std::cout << "Creating element matrix..." << std::endl;
    std::vector<std::vector< viennafem::numeric_type > >  element_matrix(3);
    element_matrix[0].resize(3);
    element_matrix[1].resize(3);
    element_matrix[2].resize(3);
    
    std::vector<viennafem::numeric_type> element_vector(3);
    
    //update cell_quantities:
    //std::cout << "Updating cell quantities..." << std::endl;
    viennamath::equation cell_expr = viennafem::update_cell_quantities(*cell_iter, 
                                                                       transformed_weak_form);
    
    //std::cout << "New cell_expr: " << cell_expr << std::endl;

    viennafem::cell_quan<CellType, viennafem::det_dF_dt_key> det_dF_dt(*cell_iter);
    
    //fill element_matrix:
    //std::cout << "Filling element matrix..." << std::endl;
    for (size_t i = 0; i<test_functions.size(); ++i)
    {
      for (size_t j=0; j<trial_functions.size(); ++j)
      {
        viennamath::equation temp = viennafem::insert_test_and_trial_functions(test_functions[i],
                                                                               trial_functions[j],
                                                                               cell_expr);
        element_matrix[i][j] = viennafem::eval_element_matrix_entry(temp.lhs()) * det_dF_dt.eval(1.0); 
      }
      
      viennamath::equation temp = viennafem::insert_test_and_trial_functions(test_functions[i],
                                                                             trial_functions[0],
                                                                             cell_expr);
      element_vector[i] = viennafem::eval_element_vector_entry(temp.rhs()) * det_dF_dt.eval(1.0); 
    }
    
    //print element matrix:
    //for (size_t i = 0; i<test_functions.size(); ++i)
    //{
    //  for (size_t j=0; j<trial_functions.size(); ++j)
    //    std::cout << element_matrix[i][j] << " ";
    //  std::cout << " | " << element_vector[i] << std::endl;
    //}
    
    //write back to global matrix:
    VertexOnCellContainer vertices_on_cell = viennagrid::ncells<0>(*cell_iter);
    long global_index_i = 0;
    long global_index_j = 0;
    long local_index_i = 0;
    long local_index_j = 0;
    for (VertexOnCellIterator vocit_i = vertices_on_cell.begin();
         vocit_i != vertices_on_cell.end();
         ++vocit_i, ++local_index_i)
    {
      global_index_i = viennadata::access<mapping_key, long>()(*vocit_i);
      if (global_index_i == -1)
        continue;
      
      local_index_j = 0;
      for (VertexOnCellIterator vocit_j = vertices_on_cell.begin();
           vocit_j != vertices_on_cell.end();
           ++vocit_j, ++local_index_j)
      {
         global_index_j = viennadata::access<mapping_key, long>()(*vocit_j);
         
         if (global_index_j == -1)
           continue; //modify right-hand side here
         
         global_matrix[global_index_i][global_index_j] += element_matrix[local_index_i][local_index_j];
      }
      
      global_rhs[global_index_i] += element_vector[local_index_i];
    }
    
  }
  
  std::cout << std::endl;
  std::cout << "Done!" << std::endl;
  std::cout << std::endl;

  //print global matrix:
  //std::cout << "Global system: " << std::endl;
  //for (size_t i = 0; i<global_matrix.size(); ++i)
  //{
  // for (size_t j=0; j<global_matrix.size(); ++j)
  //    std::cout << global_matrix[i][j] << " ";
  //  std::cout << " | " << global_rhs[i] << std::endl;
  //}
    
    
  //TODO: Solve system
  /*
  viennacl::matrix<viennafem::numeric_type> vcl_matrix(map_index, map_index);
  viennacl::vector<viennafem::numeric_type> vcl_rhs(map_index);
  
  viennacl::copy(global_matrix, vcl_matrix);
  viennacl::copy(global_rhs, vcl_rhs);
  
  viennacl::linalg::lu_factorize(vcl_matrix);
  viennacl::linalg::lu_substitute(vcl_matrix, vcl_rhs);
  
  viennacl::copy(vcl_rhs, global_rhs);
  */
  
  std::cout << "*" << std::endl;
  std::cout << "* Phase 5: Solving linear system" << std::endl;
  std::cout << "*" << std::endl;
  
  boost::numeric::ublas::matrix<viennafem::numeric_type> ublas_matrix(map_index, map_index);
  boost::numeric::ublas::vector<viennafem::numeric_type> ublas_rhs(map_index);
  
  for (size_t i = 0; i<global_matrix.size(); ++i)
  {
    for (size_t j=0; j<global_matrix.size(); ++j)
      ublas_matrix(i, j) = global_matrix[i][j];
    ublas_rhs(i) = global_rhs[i];
  }
  
  boost::numeric::ublas::lu_factorize(ublas_matrix);
  boost::numeric::ublas::inplace_solve(ublas_matrix, ublas_rhs, boost::numeric::ublas::unit_lower_tag ());
  boost::numeric::ublas::inplace_solve(ublas_matrix, ublas_rhs, boost::numeric::ublas::upper_tag ());
  
  //print solution:
  std::cout << "Solution: ";
  for (size_t i=0; i<global_rhs.size(); ++i)
    std::cout << ublas_rhs(i) << " ";
  std::cout << std::endl;
  std::cout << std::endl;
  
  //Writing solution back to mesh:
  for (VertexIterator vit = vertices.begin();
       vit != vertices.end();
       ++vit)
  {
     long cur_index = viennadata::access<mapping_key, long>()(*vit);
     if (cur_index > -1)
       viennadata::access<std::string, double>("vtk_data")(*vit) = ublas_rhs(cur_index);
  }

  std::cout << "*" << std::endl;
  std::cout << "* Phase 6: Writing data to 'tut1.vtu' (can be viewed with e.g. Paraview)" << std::endl;
  std::cout << "*" << std::endl;

  viennagrid::io::Vtk_writer<DomainType> my_vtk_writer;
  my_vtk_writer.writeDomain(domain, "tut1.vtu");  

}



int main()
{
  viennagrid::domain<TriangleConfig> my_domain;
  
  try{
    viennagrid::io::sgf_reader my_sgf_reader;
    my_sgf_reader(my_domain, "../examples/data/square128.sgf");
  } catch (...){
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  
  do_something(my_domain);
  
  std::cout << "*****************************************" << std::endl;
  std::cout << "* Poisson solver finished successfully! *" << std::endl;
  std::cout << "*****************************************" << std::endl;
  return EXIT_SUCCESS;
}
