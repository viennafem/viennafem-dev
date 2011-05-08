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
#include "viennafem/cell_quan.hpp"
#include "viennafem/transform.hpp"
#include "viennafem/bases/tetrahedron.hpp"
#include "viennafem/bases/triangle.hpp"
#include "viennafem/eval.hpp"
#include "viennafem/transform/dtdx_triangle.hpp"
#include "viennafem/transform/dtdx_tetrahedron.hpp"
#include "viennafem/weak_form.hpp"
#include "viennafem/linear_pde_system.hpp"
#include "viennafem/linear_pde_options.hpp"


//ViennaMath includes:
#include "viennamath/manipulation/apply_coordinate_system.hpp"

//ViennaData includes:
#include "viennadata/api.hpp"

//ViennaGrid includes:
#include "viennagrid/domain.hpp"

namespace viennafem
{
  
  
  template <typename CellType, typename InterfaceType>
  struct cell_updater : public viennamath::traversal_interface<>
  {
    public:
      cell_updater(CellType const & cell) : cell_(cell) {}
      
      void operator()(InterfaceType const * e) const 
      {
        if (viennamath::callback_if_castable< viennafem::cell_quan<CellType, InterfaceType> >::apply(e, *this))
          return;
      }
      
      void operator()(viennafem::cell_quan<CellType, InterfaceType> const & cq) const
      {
        cq.update(cell_);
        //std::cout << "cell_quan updated!" << std::endl;
      }
      
    private:
      CellType const & cell_;
  };


  
  struct pde_assembler
  {
    
    template <typename EquationType, typename OptionType, typename DomainType, typename MatrixType, typename VectorType>
    void operator()(EquationType const & transformed_weak_form,
                    OptionType const & options,
                    DomainType & domain,
                    MatrixType & system_matrix,
                    VectorType & load_vector
                   ) const
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

      typedef viennafem::mapping_key                              MappingKeyType;
      
      typedef typename EquationType::value_type      Expression;
      
      MappingKeyType map_key(options.unknown_id());

      //get basis:
      //std::cout << "Getting basis..." << std::endl;
      std::vector<Expression> trial_functions = viennafem::get_basisfunctions<Expression>(CellTag());
      std::vector<Expression> test_functions = viennafem::get_basisfunctions<Expression>(CellTag());
      
      
      std::vector<std::vector< EquationType > >  local_weak_form(viennagrid::traits::subcell_desc<CellTag, 0>::num_elements);
      for (size_t i=0; i<viennagrid::traits::subcell_desc<CellTag, 0>::num_elements; ++i)
        local_weak_form[i].resize(viennagrid::traits::subcell_desc<CellTag, 0>::num_elements);
      
      for (size_t i = 0; i<test_functions.size(); ++i)
        for (size_t j=0; j<trial_functions.size(); ++j)
          local_weak_form[i][j] = viennafem::insert_test_and_trial_functions(test_functions[i],
                                                                            trial_functions[j],
                                                                            transformed_weak_form);
      
      viennafem::cell_quan<CellType, typename EquationType::interface_type> det_dF_dt;
      det_dF_dt.wrap( viennafem::det_dF_dt_key() );
      
      
      CellContainer cells = viennagrid::ncells<CellTag::topology_level>(domain);
      for (CellIterator cell_iter = cells.begin();
          cell_iter != cells.end();
          ++cell_iter)
      {
        //update cell_quantities:
        //std::cout << "Updating cell quantities..." << std::endl;
        //EquationType cell_expr = viennafem::update_cell_quantities(*cell_iter, transformed_weak_form);
        viennamath::traversal_wrapper<> updater( new cell_updater<CellType, typename EquationType::interface_type>(*cell_iter) );
        det_dF_dt.update(*cell_iter);
        
        
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
          global_index_i = viennadata::access<MappingKeyType, long>(map_key)(*vocit_i);
          //std::cout << "glob_i: " << global_index_i << std::endl;
          if (global_index_i == -1)
            continue;
          
          local_index_j = 0;
          for (VertexOnCellIterator vocit_j = vertices_on_cell.begin();
              vocit_j != vertices_on_cell.end();
              ++vocit_j, ++local_index_j)
          {
            global_index_j = viennadata::access<MappingKeyType, long>(map_key)(*vocit_j);
            //std::cout << "glob_j: " << global_index_j << std::endl;
            if (global_index_j == -1)
              continue; //modify right-hand side here
          
            local_weak_form[local_index_i][local_index_j].lhs().get()->recursive_traversal(updater);
          
            //std::cout << "incrementing sys matrix at " << global_index_i << " " << global_index_j << " by " << element_matrix[local_index_i][local_index_j] << std::endl;
            system_matrix(global_index_i, global_index_j) += viennafem::eval_element_matrix_entry(local_weak_form[local_index_i][local_index_j].lhs(), CellTag()) * det_dF_dt.eval(1.0); 
          }
          
          local_weak_form[local_index_i][0].rhs().get()->recursive_traversal(updater);
          load_vector(global_index_i) += viennafem::eval_element_vector_entry(local_weak_form[local_index_i][0].rhs(), CellTag()) * det_dF_dt.eval(1.0); 
        }
      }
      
    } //operator()
  };
  
  
  
  
  

  class pde_solver
  {
    public:
      
      template <typename SystemType, typename DomainType, typename MatrixType, typename VectorType>  //template for operator()
      void operator()(SystemType pde_system,
                      DomainType & domain,
                      MatrixType & system_matrix,
                      VectorType & load_vector
                     ) const;
      
  };




  template <typename SystemType, typename DomainType, typename MatrixType, typename VectorType>  //template for operator()
  void pde_solver::operator()(SystemType pde_system,
                              DomainType & domain,
                              MatrixType & system_matrix,
                              VectorType & load_vector
                             ) const
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
    
    typedef typename SystemType::equation_type                  EquationType;
    typedef typename SystemType::equation_type::value_type      Expression;
    
    typedef viennafem::boundary_key                             BoundaryKeyType;
    typedef viennafem::mapping_key                              MappingKeyType;
    
    
    EquationType weak_form_general = viennafem::make_weak_form(pde_system.pde(0));  
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
    BoundaryKeyType bnd_key(pde_system.option(0).unknown_id());
    MappingKeyType map_key(pde_system.option(0).unknown_id());
    
    VertexContainer vertices = viennagrid::ncells<0>(domain);
    for (VertexIterator vit = vertices.begin();
        vit != vertices.end();
        ++vit)
    {  
      if (viennadata::access<BoundaryKeyType, bool>(bnd_key)(*vit))
      {
        //std::cout << "boundary vertex" << std::endl;
        viennadata::access<MappingKeyType, long>(map_key)(*vit) = -1;
      }
      else
      {
        //std::cout << "interior vertex" << std::endl;
        viennadata::access<MappingKeyType, long>(map_key)(*vit) = map_index++;
      }
    }
    std::cout << "---------------------------" << std::endl;
    
    std::cout << "* pde_solver::operator(): Assigned degrees of freedom: " << map_index << std::endl;
    
    //resize global system matrix and load vector if needed:
    if (map_index > system_matrix.size1())
    {
      std::cout << "Resizing system matrix..." << std::endl;
      system_matrix.resize(map_index, map_index, false);
      system_matrix.clear();
      system_matrix.resize(map_index, map_index, false);
    }
    
    if (map_index > load_vector.size())
    {
      std::cout << "Resizing load vector..." << std::endl;
      load_vector.resize(map_index, false);
      load_vector.clear();
      load_vector.resize(map_index, false);
    }
    
    std::cout << "* pde_solver::operator(): Transform to reference element" << std::endl;
    
    EquationType transformed_weak_form = viennafem::transform_to_reference_cell<CellType>(weak_form);
    
    std::cout << "* pde_solver::operator(): Transformed weak form:" << std::endl;
    std::cout << transformed_weak_form << std::endl;
    std::cout << std::endl;

    std::cout << "* pde_solver::operator(): Assemble system" << std::endl;
    pde_assembler()(transformed_weak_form, pde_system.option(0), domain, system_matrix, load_vector);

  }
}
#endif
