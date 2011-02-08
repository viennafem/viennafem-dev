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
#include "viennafem/transform.hpp"
#include "viennafem/BFStock.hpp"
#include "viennafem/eval.hpp"
#include "viennafem/dtdx_triangle.h"

//ViennaMath includes:
#include "viennamath/weak_form.hpp"
#include "viennamath/apply_coordinate_system.hpp"

//ViennaData includes:
#include "viennadata/interface.hpp"

//ViennaGrid includes:
#include "viennagrid/domain.hpp"

namespace viennafem
{

  class pde_solver
  {
    public:
      
      template <typename EquationType, typename ConfigType, typename PDEDomain>
      void operator()(EquationType & strong_form,
                      ConfigType & conf,
                      PDEDomain & domain) const;
      
  };




  template <typename EquationType, typename ConfigType, typename DomainType>  //template for operator()
  void pde_solver::operator()(EquationType & strong_form,
                              ConfigType & config,
                              DomainType & domain) const
  {
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
    
    typedef typename ConfigType::matrix_type               MatrixType;
    typedef typename ConfigType::vector_type               VectorType;
    typedef typename ConfigType::mapping_key_type          MappingKeyType;
    typedef typename ConfigType::boundary_key_type         BoundaryKeyType;
    
    MatrixType & system_matrix = config.system_matrix();
    VectorType & load_vector = config.load_vector();
    
    
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
   
    size_t map_index = 0;
    VertexContainer vertices = viennagrid::ncells<0>(domain);
    for (VertexIterator vit = vertices.begin();
        vit != vertices.end();
        ++vit)
    {  
      if (viennadata::access<BoundaryKeyType, bool>(config.boundary_key())(*vit))
        viennadata::access<MappingKeyType, long>(config.mapping_key())(*vit) = -1;
      else
        viennadata::access<MappingKeyType, long>(config.mapping_key())(*vit) = map_index++;
    }
    std::cout << "* pde_solver::operator(): Assigned degrees of freedom: " << map_index << std::endl;
    
    //resize global system matrix and load vector if needed:
    if (map_index > system_matrix.size1())
    {
      system_matrix.resize(map_index, map_index, false);
      system_matrix.clear();
      system_matrix.resize(map_index, map_index, false);
    }
    
    if (map_index > load_vector.size())
    {
      load_vector.resize(map_index, false);
      load_vector.clear();
      load_vector.resize(map_index, false);
    }
    
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
        global_index_i = viennadata::access<MappingKeyType, long>(config.mapping_key())(*vocit_i);
        if (global_index_i == -1)
          continue;
        
        local_index_j = 0;
        for (VertexOnCellIterator vocit_j = vertices_on_cell.begin();
            vocit_j != vertices_on_cell.end();
            ++vocit_j, ++local_index_j)
        {
          global_index_j = viennadata::access<MappingKeyType, long>(config.mapping_key())(*vocit_j);
          
          if (global_index_j == -1)
            continue; //modify right-hand side here
          
          system_matrix(global_index_i, global_index_j) += element_matrix[local_index_i][local_index_j];
        }
        
        load_vector(global_index_i) += element_vector[local_index_i];
      }
      
    }
    
  }
}
#endif
