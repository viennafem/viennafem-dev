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
#include "viennafem/dtdx_tetrahedron.h"
#include "viennafem/weak_form.hpp"

//ViennaMath includes:
#include "viennamath/manipulation/apply_coordinate_system.hpp"

//ViennaData includes:
#include "viennadata/api.hpp"

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
    
    
    typedef typename EquationType::value_type      Expression;
    
    
    MatrixType & system_matrix = config.system_matrix();
    VectorType & load_vector = config.load_vector();
    
    
    EquationType weak_form_general = viennafem::make_weak_form(strong_form);  
    EquationType weak_form = viennamath::apply_coordinate_system(viennamath::cartesian<Config::dimension_tag::value>(), weak_form_general);
    
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
      {
        //std::cout << "boundary vertex" << std::endl;
        viennadata::access<MappingKeyType, long>(config.mapping_key())(*vit) = -1;
      }
      else
      {
        //std::cout << "interior vertex" << std::endl;
        viennadata::access<MappingKeyType, long>(config.mapping_key())(*vit) = map_index++;
      }
    }
    std::cout << "---------------------------" << std::endl;
    
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
    
    EquationType transformed_weak_form = viennafem::transform_to_reference_cell<CellType>(weak_form);
    
    std::cout << "* pde_solver::operator(): Transformed weak form:" << std::endl;
    std::cout << transformed_weak_form << std::endl;
    std::cout << std::endl;

    std::cout << "* pde_solver::operator(): Assemble local element matrix" << std::endl;

    //get basis:
    //std::cout << "Getting basis..." << std::endl;
    std::vector<Expression> trial_functions = viennafem::get_basisfunctions<Expression>(CellTag());
    std::vector<Expression> test_functions = viennafem::get_basisfunctions<Expression>(CellTag());
    
    
    std::vector<std::vector< EquationType > >  local_weak_form(viennagrid::subcell_traits<CellTag, 0>::num_elements);
    for (size_t i=0; i<viennagrid::subcell_traits<CellTag, 0>::num_elements; ++i)
      local_weak_form[i].resize(viennagrid::subcell_traits<CellTag, 0>::num_elements);
    
    for (size_t i = 0; i<test_functions.size(); ++i)
      for (size_t j=0; j<trial_functions.size(); ++j)
        local_weak_form[i][j] = viennafem::insert_test_and_trial_functions(test_functions[i],
                                                                           trial_functions[j],
                                                                           transformed_weak_form);
    
    
    //transfer cell quantities to vtk:
    for (CellIterator cell_iter = cells.begin();
        cell_iter != cells.end();
        ++cell_iter)
    {
      
      //set up element matrix:
      //std::cout << "Creating element matrix..." << std::endl;
      std::vector<std::vector< viennafem::numeric_type > >  element_matrix(viennagrid::subcell_traits<CellTag, 0>::num_elements);
      for (size_t i=0; i<viennagrid::subcell_traits<CellTag, 0>::num_elements; ++i)
        element_matrix[i].resize(viennagrid::subcell_traits<CellTag, 0>::num_elements);
      
      std::vector<viennafem::numeric_type> element_vector(viennagrid::subcell_traits<CellTag, 0>::num_elements);
      
      //update cell_quantities:
      //std::cout << "Updating cell quantities..." << std::endl;
      //EquationType cell_expr = viennafem::update_cell_quantities(*cell_iter, transformed_weak_form);
      transformed_weak_form.lhs().get()->update_cell(*cell_iter);
      transformed_weak_form.rhs().get()->update_cell(*cell_iter);
      //EquationType & cell_expr = transformed_weak_form;
      
      //std::cout << "New cell_expr: " << cell_expr << std::endl;

      viennafem::cell_quan<CellType, viennafem::det_dF_dt_key, typename EquationType::interface_type> det_dF_dt(*cell_iter);

      /*
      viennafem::cell_quan<CellType, viennafem::dt_dx_key<0,0> >  dr_dx(*cell_iter);
      viennafem::cell_quan<CellType, viennafem::dt_dx_key<1,0> >  ds_dx(*cell_iter);
      viennafem::cell_quan<CellType, viennafem::dt_dx_key<2,0> >  dt_dx(*cell_iter);
      
      viennafem::cell_quan<CellType, viennafem::dt_dx_key<0,1> >  dr_dy(*cell_iter);
      viennafem::cell_quan<CellType, viennafem::dt_dx_key<1,1> >  ds_dy(*cell_iter);
      viennafem::cell_quan<CellType, viennafem::dt_dx_key<2,1> >  dt_dy(*cell_iter);

      viennafem::cell_quan<CellType, viennafem::dt_dx_key<0,2> >  dr_dz(*cell_iter);
      viennafem::cell_quan<CellType, viennafem::dt_dx_key<1,2> >  ds_dz(*cell_iter);
      viennafem::cell_quan<CellType, viennafem::dt_dx_key<2,2> >  dt_dz(*cell_iter);
      
      std::cout << det_dF_dt.eval(1.0) << std::endl;
      std::cout << dr_dx.eval(0.0) << std::endl;
      std::cout << dr_dy.eval(0.0) << std::endl;
      std::cout << dr_dz.eval(0.0) << std::endl;
      
      std::cout << ds_dx.eval(0.0) << std::endl;
      std::cout << ds_dy.eval(0.0) << std::endl;
      std::cout << ds_dz.eval(0.0) << std::endl;
      
      std::cout << dt_dx.eval(0.0) << std::endl;
      std::cout << dt_dy.eval(0.0) << std::endl;
      std::cout << dt_dz.eval(0.0) << std::endl; */
      
      //fill element_matrix:
      //std::cout << "Filling element matrix..." << std::endl;
      //cell_iter->print_short();
      for (size_t i = 0; i<test_functions.size(); ++i)
      {
        for (size_t j=0; j<trial_functions.size(); ++j)
        {
          local_weak_form[i][j].lhs().get()->update_cell(*cell_iter);
          
          //EquationType temp = viennafem::insert_test_and_trial_functions(test_functions[i],
          //                                                                         trial_functions[j],
          //                                                                         cell_expr);
          element_matrix[i][j] = viennafem::eval_element_matrix_entry(local_weak_form[i][j].lhs(), CellTag()) * det_dF_dt.eval(1.0); 
        }
        
        //EquationType temp = viennafem::insert_test_and_trial_functions(test_functions[i],
        //                                                                         trial_functions[0],
        //                                                                         cell_expr);
        local_weak_form[i][0].rhs().get()->update_cell(*cell_iter);
        element_vector[i] = viennafem::eval_element_vector_entry(local_weak_form[i][0].rhs(), CellTag()) * det_dF_dt.eval(1.0); 
        
        //std::cout << std::endl;
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
        //std::cout << "glob_i: " << global_index_i << std::endl;
        if (global_index_i == -1)
          continue;
        
        local_index_j = 0;
        for (VertexOnCellIterator vocit_j = vertices_on_cell.begin();
            vocit_j != vertices_on_cell.end();
            ++vocit_j, ++local_index_j)
        {
          global_index_j = viennadata::access<MappingKeyType, long>(config.mapping_key())(*vocit_j);
          //std::cout << "glob_j: " << global_index_j << std::endl;
          if (global_index_j == -1)
            continue; //modify right-hand side here
          
          //std::cout << "incrementing sys matrix at " << global_index_i << " " << global_index_j << " by " << element_matrix[local_index_i][local_index_j] << std::endl;
          system_matrix(global_index_i, global_index_j) += element_matrix[local_index_i][local_index_j];
        }
        
        load_vector(global_index_i) += element_vector[local_index_i];
      }
      
    }
  }
}
#endif
