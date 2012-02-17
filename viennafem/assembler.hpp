/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_ASSEMBLER_HPP
#define VIENNAFEM_ASSEMBLER_HPP

//ViennaFEM includes:
#include "viennafem/forwards.h"
#include "viennafem/cell_quan.hpp"
#include "viennafem/transform.hpp"
#include "viennafem/bases/all.hpp"
#include "viennafem/eval.hpp"
#include "viennafem/transform/dtdx_triangle.hpp"
#include "viennafem/transform/dtdx_tetrahedron.hpp"
#include "viennafem/weak_form.hpp"
#include "viennafem/linear_pde_system.hpp"
#include "viennafem/linear_pde_options.hpp"
#include "viennafem/mapping.hpp"
#include "viennafem/quadrature/quad.hpp"

//ViennaMath includes:
#include "viennamath/manipulation/apply_coordinate_system.hpp"
#include "viennamath/runtime/numerical_quadrature.hpp"
#include "viennamath/manipulation/eval.hpp"

//ViennaData includes:
#include "viennadata/api.hpp"

//ViennaGrid includes:
#include "viennagrid/domain.hpp"

namespace viennafem
{
  template <typename ExpressionType>
  std::vector<ExpressionType> make_full_function(ExpressionType const & func, size_t length, size_t index)
  {
    std::vector<ExpressionType> result(length);
    for (size_t i=0; i<length; ++i)
      result[i] = 0.0;
    
    result[index] = func;
    
    return result;
  }
  
  
  template <typename CellTag, typename EquationType, typename PDESystemType>
  std::vector< std::vector<EquationType> > make_local_weak_form(EquationType const & transformed_weak_form, PDESystemType const & pde_system)
  {
    typedef typename EquationType::value_type      Expression;
    typedef typename Expression::interface_type    InterfaceType;
    typedef typename reference_cell_for_basis<CellTag, viennafem::lagrange_tag<1> >::type    ReferenceCell;
    

    //test functions:
    std::vector<Expression> scalar_test_functions = viennafem::basis_factory<InterfaceType>::get(pde_system.option(0).test_space_id(), ReferenceCell());
    for (std::size_t i = 0; i<scalar_test_functions.size(); ++i)
      std::cout << "Test function " << i << ": " << scalar_test_functions[i] << std::endl;
    
    size_t local_size_i = pde_system.unknown(0).size() * scalar_test_functions.size();

    //std::cout << "Test functions: " << std::endl;
    std::vector< std::vector<Expression> > full_test_functions(local_size_i);
    size_t current_index = 0;
    for (size_t i=0; i<scalar_test_functions.size(); ++i)
    {
      for (size_t j=0; j<pde_system.unknown(0).size(); ++j)
      {  
        //std::cout << "No. " << current_index << ": ";
        full_test_functions[current_index++] = make_full_function(scalar_test_functions[i], pde_system.unknown(0).size(), j);
        //for (size_t k=0; k<pde_system.unknown(0).size(); ++k)
        //  std::cout << "[" << full_test_functions[current_index-1][k] << "]";
        //std::cout << std::endl;
      }
      //std::cout << std::endl;
    }

    //trial functions:
    std::vector<Expression> scalar_trial_functions = viennafem::basis_factory<InterfaceType>::get(pde_system.option(0).trial_space_id(), ReferenceCell());
    size_t local_size_j = pde_system.unknown(0).size() * scalar_trial_functions.size();
      
    //std::cout << "Trial functions: " << std::endl;
    std::vector< std::vector<Expression> > full_trial_functions(local_size_i);
    current_index = 0;
    for (size_t i=0; i<scalar_trial_functions.size(); ++i)
    {
      //std::cout << "No. i: ";
      for (size_t j=0; j<pde_system.unknown(0).size(); ++j)
      {  
        //std::cout << "No. " << current_index << ": ";
        full_trial_functions[current_index++] = make_full_function(scalar_trial_functions[i], pde_system.unknown(0).size(), j);
        //for (size_t k=0; k<pde_system.unknown(0).size(); ++k)
        //  std::cout << "[" << full_trial_functions[current_index-1][k] << "]";
        //std::cout << std::endl;
      }
      //std::cout << std::endl;
    }
    

    //plug functions into weak form to obtain local form:
    //local_size_i = 3;
    //local_size_j = 3;
    
    std::vector<std::vector< EquationType > >  local_weak_form(local_size_i);
    for (size_t i=0; i<local_size_i; ++i)
      local_weak_form[i].resize(local_size_j);
    
    for (size_t i = 0; i<local_size_i; ++i)
    {
      for (size_t j = 0; j<local_size_j; ++j)
      {
        local_weak_form[i][j] = viennafem::insert_test_and_trial_functions( transformed_weak_form,
                                                                            pde_system,
                                                                            full_test_functions[i],
                                                                            full_trial_functions[j]);
        /*local_weak_form[i][j] = viennafem::insert_test_and_trial_functions_vector( transformed_weak_form,
                                                                            pde_system,
                                                                            scalar_test_functions[0],
                                                                            scalar_trial_functions[0],
                                                                            i, j);*/
      }
    }
    
    return local_weak_form;    
  }
  
  
  template <typename CellType, typename InterfaceType>
  struct cell_updater : public viennamath::rt_traversal_interface<>
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


  template <typename T>
  bool is_uniform_basis(T const &) { return true; }   //TODO: Extend
  
  
  struct pde_assembler_internal
  {
    
    template <typename EquationType, typename PDESystemType, typename DomainType, typename LinearSystemT>
    void operator()(EquationType const & transformed_weak_form,
                    PDESystemType const & pde_system,
                    DomainType & domain,
                    LinearSystemT & linear_system
                   ) const
    {
      typedef typename DomainType::config_type              Config;
      typedef typename Config::cell_tag                     CellTag;
      
      typedef typename viennagrid::result_of::point<Config>::type                            PointType;
      typedef typename viennagrid::result_of::ncell<Config, CellTag::dim>::type              CellType;

      typedef typename viennagrid::result_of::ncell_range<DomainType, 0>::type               VertexContainer;
      typedef typename viennagrid::result_of::iterator<VertexContainer>::type                VertexIterator;

      typedef typename viennagrid::result_of::ncell_range<DomainType, CellTag::dim>::type    CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type                  CellIterator;

      typedef typename viennagrid::result_of::ncell_range<CellType, 0>::type                 VertexOnCellContainer;
      typedef typename viennagrid::result_of::iterator<VertexOnCellContainer>::type          VertexOnCellIterator;

      typedef typename PDESystemType::boundary_key_type  BoundaryKeyType;
      typedef std::vector<long>                          MappingContainer;
      
      typedef typename EquationType::value_type          Expression;
      
      BoundaryKeyType bnd_key(pde_system.option(0).data_id());
      
      std::cout << "Number of components: " << pde_system.unknown(0).size() << std::endl;

      
      //Set up element representations:
      std::vector<std::vector< EquationType > > local_weak_form;
      bool basis_is_uniform = is_uniform_basis(pde_system);
      
      if (basis_is_uniform)
      {
        //std::cout << "Using globally uniform basis";
        local_weak_form = make_local_weak_form<CellTag>(transformed_weak_form, pde_system);
      }
          
      /*
      for (size_t i=0; i<local_weak_form.size(); ++i)
      {
        for (size_t j=0; j<local_weak_form[i].size(); ++j)
        {
          std::cout << "(" << i <<  "," << j << "): " << local_weak_form[i][j] << std::endl;
        }
        std::cout << std::endl;
      }*/
      //exit(0);
      
      //Integrator setup:
      
      viennamath::numerical_quadrature integrator = viennafem::make_quadrature_rule(pde_system, domain);
      //viennamath::numerical_quadrature integrator(new viennafem::rt_gauss_quad_element<ReferenceCell, 3, typename EquationType::interface_type>());
      
      CellContainer cells = viennagrid::ncells<CellTag::dim>(domain);
      for (CellIterator cell_iter = cells.begin();
          cell_iter != cells.end();
          ++cell_iter)
      {
        //if (!basis_is_uniform)
        //  local_weak_form = make_local_weak_form(transformed_weak_form, pde_system, *cell_iter);
        
        //update cell_quantities:
        //std::cout << "Updating cell quantities..." << std::endl;
        viennamath::rt_traversal_wrapper<> updater( new cell_updater<CellType, typename EquationType::interface_type>(*cell_iter) );
        
        //write back to global matrix:
        long global_index_i = 0;
        long global_index_j = 0;
        long local_index_i = 0;
        long local_index_j = 0;

        MappingContainer map_indices_i = mapping_indices(pde_system, *cell_iter, 0);
        MappingContainer map_indices_j = mapping_indices(pde_system, *cell_iter, 0);
        
        for (typename MappingContainer::const_iterator map_iter_i = map_indices_i.begin();
             map_iter_i != map_indices_i.end();
             ++map_iter_i, ++local_index_i)
        {                  
          global_index_i = *map_iter_i;
          //std::cout << "glob_i: " << global_index_i << std::endl;
          
          if (global_index_i == -1) // This is a Dirichlet node -> skip
            continue;
          
          VertexOnCellContainer vertices_on_cell = viennagrid::ncells<0>(*cell_iter);
          VertexOnCellIterator vocit = vertices_on_cell.begin();
          local_index_j = 0;
          for (typename MappingContainer::const_iterator map_iter_j = map_indices_j.begin();
               map_iter_j != map_indices_j.end();
               ++map_iter_j, ++local_index_j)
          {
            if ( (local_index_j % pde_system.unknown(0).size()) == 0 && (local_index_j > 0) )
                ++vocit;
            
            local_weak_form[local_index_i][local_index_j].lhs().get()->recursive_traversal(updater);
            global_index_j = *map_iter_j;
            

            if (global_index_j == -1) // Dirichlet boundary
            {
              if (pde_system.unknown(0).size() == 1) //scalar valued unknowns
              {
                linear_system(global_index_i) -=
                   viennadata::access<BoundaryKeyType, double>(bnd_key)(*vocit) //TODO: Better encapsulate boundary value access
                   * integrator(local_weak_form[local_index_i][local_index_j].lhs());
              }
              else //vector valued unknowns
              {
                //TODO: Better encapsulate boundary value access
                std::vector<double> const & bnd_values = viennadata::access<BoundaryKeyType, std::vector<double> >(bnd_key)(*vocit);
                
                if (bnd_values.size() > 1) //allow homogeneous case without having the data vector initialized
                {
                  linear_system(global_index_i) -=
                     bnd_values[local_index_j % pde_system.unknown(0).size()] 
                     * integrator(local_weak_form[local_index_i][local_index_j].lhs());
                }
              }
            }
            else
            {
              //std::cout << " Evaluating LHS: " << local_weak_form[local_index_i][local_index_j].lhs() << std::endl;
              //std::cout << "Adding to (" << global_index_i << ", " << global_index_j << "): " << integrator(local_weak_form[local_index_i][local_index_j].lhs()) << std::endl;
               linear_system.add(
                  global_index_i, 
                  global_index_j, 
                  integrator(local_weak_form[local_index_i][local_index_j].lhs())
               );
            }
          }
          
          local_weak_form[local_index_i][0].rhs().get()->recursive_traversal(updater);
          
          //std::cout << "Evaluating RHS: " << local_weak_form[local_index_i][0].rhs() << std::endl;
          
          linear_system(global_index_i) += integrator(local_weak_form[local_index_i][0].rhs());
        }
      }
      
    } //operator()
    
  };
  
  
}
#endif
