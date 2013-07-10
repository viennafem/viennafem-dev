#ifndef VIENNAFEM_PDE_ASSEMBLER_HPP
#define VIENNAFEM_PDE_ASSEMBLER_HPP

/* =========================================================================
   Copyright (c) 2012, Institute for Microelectronics,
                       Institute for Analysis and Scientific Computing,
                       TU Wien.
                             -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                             -----------------

   Author:     Karl Rupp                          rupp@iue.tuwien.ac.at

   License:    MIT (X11), see file LICENSE in the ViennaFEM base directory
============================================================================ */

//ViennaFEM includes:
#include "viennafem/forwards.h"
#include "viennafem/cell_quan.hpp"
#include "viennafem/transform.hpp"
#include "viennafem/bases/tetrahedron.hpp"
#include "viennafem/bases/triangle.hpp"
#include "viennafem/transform/dtdx_interval.hpp"
#include "viennafem/transform/dtdx_triangle.hpp"
#include "viennafem/transform/dtdx_quadrilateral.hpp"
#include "viennafem/transform/dtdx_tetrahedron.hpp"
#include "viennafem/transform/dtdx_hexahedron.hpp"
#include "viennafem/weak_form.hpp"
#include "viennafem/linear_pde_system.hpp"
#include "viennafem/linear_pde_options.hpp"
#include "viennafem/detail/assembler.hpp"
#include "viennafem/log/api.hpp"


////ViennaMath includes:
#include "viennamath/manipulation/apply_coordinate_system.hpp"

////ViennaData includes:
#include "viennadata/api.hpp"

////ViennaGrid includes:
#include "viennagrid/config/default_configs.hpp"

//#define VIENNAFEM_DEBUG

/** @file   pde_assembler.hpp
    @brief  Defines the assembly function to be used by the library user.
*/

namespace viennafem
{

  namespace detail
  {
    /** @brief A simple wrapper class which abstracts a matrix and a vector into a linear equation */
    template<typename MatrixT, typename VectorT>
    struct equation_wrapper
    {
      equation_wrapper(MatrixT& matrix, VectorT& vector) : matrix(matrix), vector(vector) {}

      template<typename IndexT>
      typename VectorT::value_type& operator()(IndexT i)
      {
        return vector(i);
      }

      template<typename IndexT, typename NumericT>
      inline void add(IndexT col, IndexT row, NumericT value)
      {
        matrix(col,row) += value;
      }

      MatrixT& matrix; 
      VectorT& vector;
    };
    
  } //namespace detail

  /** @brief The main ViennaFEM assembler class */
  template<typename StorageType>
  class pde_assembler
  {
    public:
      
      pde_assembler(StorageType& storage) : storage(storage) {}
      
      /** @brief Functor interface for the assembly.
       * 
       * @param pde_system     The system of PDEs
       * @param domain         The ViennaGrid domain on which assembly is carried out
       * @param system_matrix  The system matrix. Can be any type supporting operator() access
       * @param load_vector    The load vector. Can be any type supporting operator() access
       */
      template <typename SystemType, typename DomainType, typename MatrixT, typename VectorT>  //template for operator()
      void operator()(SystemType pde_system,
                      DomainType & domain,
                      MatrixT    & system_matrix, 
                      VectorT    & load_vector
                     ) const
      {
        typedef typename viennagrid::result_of::cell_tag<DomainType>::type CellTag;
      
        typedef typename viennagrid::result_of::point<DomainType>::type                                   PointType;
        typedef typename viennagrid::result_of::element<DomainType, viennagrid::vertex_tag>::type         VertexType;
        typedef typename viennagrid::result_of::element<DomainType, CellTag>::type                        CellType;

        typedef typename viennagrid::result_of::element_range<DomainType, viennagrid::vertex_tag>::type   VertexContainer;
        typedef typename viennagrid::result_of::iterator<VertexContainer>::type                           VertexIterator;

        typedef typename viennagrid::result_of::element_range<DomainType, CellTag>::type                  CellContainer;
        typedef typename viennagrid::result_of::iterator<CellContainer>::type                             CellIterator;

        typedef typename viennagrid::result_of::element_range<CellType, VertexType>::type                 VertexOnCellContainer;
        typedef typename viennagrid::result_of::iterator<VertexOnCellContainer>::type                     VertexOnCellIterator;
        
        typedef typename SystemType::equation_type                  EquationType;
        typedef typename SystemType::equation_type::value_type      Expression;
        
     #ifdef VIENNAFEM_DEBUG
        std::cout << "Strong form: " << pde_system.pde(0) << std::endl;
     #endif
        log_strong_form(pde_system);
        EquationType weak_form_general = viennafem::make_weak_form(pde_system.pde(0));  
     #ifdef VIENNAFEM_DEBUG        
        std::cout << "* pde_solver::operator(): Using weak form general: " << weak_form_general << std::endl;
     #endif
        std::vector<EquationType> temp(1); temp[0] = weak_form_general;
        log_weak_form(temp, pde_system);
        EquationType weak_form = viennamath::apply_coordinate_system(viennamath::cartesian< PointType::dim >(), weak_form_general);
        //EquationType weak_form = viennamath::apply_coordinate_system(viennamath::cartesian<Config::coordinate_system_tag::dim>(), weak_form_general);
        temp[0] = weak_form;
        log_coordinated_weak_form(temp, pde_system);

     #ifdef VIENNAFEM_DEBUG        
        std::cout << "* pde_solver::operator(): Using weak form " << weak_form << std::endl;
        std::cout << "* pde_solver::operator(): Write dt_dx coefficients" << std::endl;
     #endif
     
        typedef typename reference_cell_for_basis<CellTag, viennafem::lagrange_tag<1> >::type    ReferenceCell;
        
        //
        // Create accessors for performance in the subsequent dt_dx_handler step
        //
        
        //viennafem::dtdx_assigner<DomainType, StorageType, ReferenceCell>::apply(domain, storage);
        
        viennafem::dt_dx_handler<DomainType, StorageType, ReferenceCell>  dt_dx_handler(domain, storage);
        
        //fill with cell quantities 
        CellContainer cells = viennagrid::elements<CellType>(domain);  
        for (CellIterator cell_iter = cells.begin();
            cell_iter != cells.end();
            ++cell_iter)
        {
          //cell_iter->print_short();
          //viennadata::access<example_key, double>()(*cell_iter) = i; 
          //viennafem::dt_dx_handler<ReferenceCell>::apply(storage, *cell_iter);
          dt_dx_handler(*cell_iter);
        }

     #ifdef VIENNAFEM_DEBUG        
        std::cout << "* pde_solver::operator(): Create Mapping:" << std::endl;
     #endif
        std::size_t map_index = create_mapping(storage, pde_system, domain);
        
     #ifdef VIENNAFEM_DEBUG                
        std::cout << "* pde_solver::operator(): Assigned degrees of freedom in domain so far: " << map_index << std::endl;
     #endif        
        // resize global system matrix and load vector if needed:
        // TODO: This can be a performance bottleneck for large numbers of segments! (lots of resize operations...)
        if (map_index > system_matrix.size1())
        {
          MatrixT temp = system_matrix;
          ////std::cout << "Resizing system matrix..." << std::endl;
          system_matrix.resize(map_index, map_index, false);
          system_matrix.clear();
          system_matrix.resize(map_index, map_index, false);
          for (typename MatrixT::iterator1 row_it = temp.begin1();
               row_it != temp.end1();
               ++row_it)
          {
            for (typename MatrixT::iterator2 col_it = row_it.begin();
                 col_it != row_it.end();
                 ++col_it)
                 system_matrix(col_it.index1(), col_it.index2()) = *col_it;
          }
        }
        if (map_index > load_vector.size())
        {
          VectorT temp = load_vector;
       #ifdef VIENNAFEM_DEBUG                          
          std::cout << "Resizing load vector..." << std::endl;
       #endif
          load_vector.resize(map_index, false);
          load_vector.clear();
          load_vector.resize(map_index, false);
          for (std::size_t i=0; i<temp.size(); ++i)
            load_vector(i) = temp(i);
        }

     #ifdef VIENNAFEM_DEBUG                        
        std::cout << "* pde_solver::operator(): Transform to reference element" << std::endl;
     #endif
        EquationType transformed_weak_form = viennafem::transform_to_reference_cell<CellType>(storage, weak_form, pde_system);
        temp[0] = transformed_weak_form;
        log_transformed_weak_form<CellType>(storage, temp, pde_system);
        
        std::cout << "* pde_solver::operator(): Transformed weak form:" << std::endl;
        std::cout << transformed_weak_form << std::endl;
        //std::cout << std::endl;

     #ifdef VIENNAFEM_DEBUG                        
        std::cout << "* pde_solver::operator(): Assemble system" << std::endl;
     #endif
     
        typedef detail::equation_wrapper<MatrixT, VectorT>    wrapper_type;
        wrapper_type wrapper(system_matrix, load_vector);
     
        detail::pde_assembler_internal()(storage, transformed_weak_form, pde_system, domain, wrapper);
//        pde_assembler_internal()(transformed_weak_form, pde_system, domain, system_matrix, load_vector);

      }
      
      // -----------------------------------------------------------------------
      //
      StorageType& storage;
  };

}
#endif
