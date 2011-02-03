/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_PDE_SOLVER_HPP
#define VIENNAFEM_PDE_SOLVER_HPP

//ViennaFEM includes:
#include "viennafem/forwards.h"

//ViennaMath includes:
#include "viennamath/weak_form.hpp"
#include "viennamath/apply_coordinate_system.hpp"

//ViennaData includes:
#include "viennadata/interface.hpp"

//ViennaCL includes:
#ifndef VIENNACL_HAVE_UBLAS
 #define VIENNACL_HAVE_UBLAS
#endif
    
#ifdef USE_OPENCL
  #include "viennacl/matrix.hpp"
  #include "viennacl/vector.hpp"
#endif
#include "viennacl/linalg/cg.hpp"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>


namespace viennafem
{

  template <typename PDEConfig>
  class pde_solver
  {
    public:
      pde_solver(PDEConfig const & conf) : conf_(conf) {} 
      
      template <typename PDEDomain>
      std::vector<numeric_type> operator()(viennamath::equation const & strong_form,
                                           PDEDomain & domain) const;
      
    private:
      PDEConfig conf_;    
  };




  template <typename PDEConfig>  //template for class
  template <typename PDEDomain>  //template for operator()
  std::vector<numeric_type> pde_solver<PDEConfig>::operator()(viennamath::equation const & strong_form,
                                                              PDEDomain & domain) const
  {
    typedef PDEDomain                                     DomainType;
    
    typedef typename DomainType::config_type              Config;
    typedef typename Config::cell_tag                     CellTag;
    
    typedef typename viennagrid::result_of::point_type<Config>::type                            PointType;
    typedef typename viennagrid::result_of::ncell_type<Config, CellTag::topology_level>::type   CellType;

    typedef typename viennagrid::result_of::ncell_container<DomainType, 0>::type                VertexContainer;
    typedef typename viennagrid::result_of::iterator<VertexContainer>::type                     VertexIterator;

    typedef typename viennagrid::result_of::ncell_container<DomainType, CellTag::topology_level>::type    CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                                 CellIterator;

    typedef typename viennagrid::result_of::ncell_container<CellType, 0>::type                  VertexOnCellContainer;
    typedef typename viennagrid::result_of::iterator<VertexOnCellContainer>::type               VertexOnCellIterator;
    
    typedef typename PDEConfig::mapping_key_type          MappingKeyType;
    typedef typename PDEConfig::boundary_key_type         BoundaryKeyType;
    
    
    viennamath::equation weak_form_general = viennamath::weak_form(strong_form);  
    viennamath::equation weak_form = viennamath::apply_coordinate_system(viennamath::cartesian<Config::dimension_tag::value>(),
                                                                         weak_form_general);
    
    std::cout << "* pde_solver::operator(): Using weak form " << weak_form << std::endl;
    
    

    std::cout << "* pde_solver::operator(): Write dt_dx coefficients" << std::endl;
    //fill with cell quantities 
    CellContainer cells = viennagrid::ncells<CellTag::topology_level>(domain);
    for (CellIterator cell_iter = cells.begin();
        cell_iter != cells.end();
        ++cell_iter)
    {
      //cell_iter->print_short();
      //viennadata::access<example_key, double>()(*cell_iter) = i; 
      viennafem::dt_dx_handler<CellTag>::apply(*cell_iter);
    }

    std::cout << "* pde_solver::operator(): Create Mapping:" << std::endl;
   
    long map_index = 0;
    VertexContainer vertices = viennagrid::ncells<0>(domain);
    for (VertexIterator vit = vertices.begin();
        vit != vertices.end();
        ++vit)
    {
      if (viennadata::access<BoundaryKeyType, bool>(conf_.boundary_key())(*vit))
        viennadata::access<MappingKeyType, long>(conf_.mapping_key())(*vit) = -1;
      else
        viennadata::access<MappingKeyType, long>(conf_.mapping_key())(*vit) = map_index++;
    }
    std::cout << "* pde_solver::operator(): Assigned degrees of freedom: " << map_index << std::endl;
    
    //build global system matrix and load vector:
    boost::numeric::ublas::compressed_matrix<viennafem::numeric_type> ublas_matrix(map_index, map_index);
    boost::numeric::ublas::vector<viennafem::numeric_type> ublas_rhs(map_index);
    ublas_rhs.clear(); ublas_rhs.resize(map_index); //mind that ublas-vector may not be initialized to zero!!
    boost::numeric::ublas::vector<viennafem::numeric_type> ublas_result(map_index);
    
    std::cout << "* pde_solver::operator(): Transform to reference element" << std::endl;
    
    viennamath::equation transformed_weak_form = viennafem::transform_to_reference_cell<CellType>(weak_form);
    
    std::cout << "* pde_solver::operator(): Transformed weak form:" << std::endl;
    std::cout << transformed_weak_form << std::endl;
    std::cout << std::endl;

    std::cout << "* pde_solver::operator(): Assemble local element matrix" << std::endl;
    
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
        global_index_i = viennadata::access<MappingKeyType, long>(conf_.mapping_key())(*vocit_i);
        if (global_index_i == -1)
          continue;
        
        local_index_j = 0;
        for (VertexOnCellIterator vocit_j = vertices_on_cell.begin();
            vocit_j != vertices_on_cell.end();
            ++vocit_j, ++local_index_j)
        {
          global_index_j = viennadata::access<MappingKeyType, long>(conf_.mapping_key())(*vocit_j);
          
          if (global_index_j == -1)
            continue; //modify right-hand side here
          
          ublas_matrix(global_index_i, global_index_j) += element_matrix[local_index_i][local_index_j];
        }
        
        ublas_rhs(global_index_i) += element_vector[local_index_i];
      }
      
    }
    
    //      
    // Solve system of linear equations:
    //
    std::cout << "* pde_solver::operator(): Solving linear system" << std::endl;

  #ifdef USE_OPENCL
    viennacl::matrix<viennafem::numeric_type> vcl_matrix(map_index, map_index);
    viennacl::vector<viennafem::numeric_type> vcl_rhs(map_index);
    viennacl::vector<viennafem::numeric_type> vcl_result(map_index);
    
    viennacl::copy(global_matrix, vcl_matrix);
    viennacl::copy(global_rhs, vcl_rhs);
    
    vcl_result = viennacl::linalg::solve(vcl_matrix, vcl_rhs, viennacl::linalg::cg_tag());
    
    viennacl::copy(vcl_result, ublas_result);
  #else
    ublas_result = viennacl::linalg::solve(ublas_matrix, ublas_rhs, viennacl::linalg::cg_tag());
    std::cout << "Residual: " << norm_2(prod(ublas_matrix, ublas_result) - ublas_rhs) << std::endl;
  #endif
      
    std::cout << ublas_rhs << std::endl;
    
    //print solution:
    std::cout << "Solution: ";
    for (size_t i=0; i<ublas_result.size(); ++i)
      std::cout << ublas_result(i) << " ";
    std::cout << std::endl;
    std::cout << std::endl;
  
    std::vector<numeric_type> result(map_index);
    for (size_t i=0; i<ublas_result.size(); ++i)
      result[i] = ublas_result(i);
    
    return result;
  }
}
#endif
