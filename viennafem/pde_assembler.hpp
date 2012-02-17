/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_PDE_ASSEMBLER_HPP
#define VIENNAFEM_PDE_ASSEMBLER_HPP

//ViennaFEM includes:
#include "viennafem/forwards.h"
#include "viennafem/cell_quan.hpp"
#include "viennafem/transform.hpp"
#include "viennafem/bases/tetrahedron.hpp"
#include "viennafem/bases/triangle.hpp"
#include "viennafem/eval.hpp"
#include "viennafem/transform/dtdx_interval.hpp"
#include "viennafem/transform/dtdx_triangle.hpp"
#include "viennafem/transform/dtdx_tetrahedron.hpp"
#include "viennafem/transform/dtdx_quadrilateral.hpp"
#include "viennafem/transform/dtdx_hexahedron.hpp"
#include "viennafem/weak_form.hpp"
#include "viennafem/linear_pde_system.hpp"
#include "viennafem/linear_pde_options.hpp"
#include "viennafem/assembler.hpp"


////ViennaMath includes:
#include "viennamath/manipulation/apply_coordinate_system.hpp"

////ViennaData includes:
#include "viennadata/api.hpp"

////ViennaGrid includes:
#include "viennagrid/domain.hpp"

//#define VIENNAFEMDEBUG

namespace viennafem
{

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


  class pde_assembler
  {
    public:
      
      /// specialization for a linear solver object
      template <typename SystemType, typename DomainType, typename LinSolverT>  //template for operator()
      void operator()(SystemType pde_system,
                      DomainType & domain,
                      LinSolverT & linsolver
                     ) const
      {
        typedef typename DomainType::config_type              Config;
        typedef typename Config::cell_tag                     CellTag;
        
        typedef typename viennagrid::result_of::point<Config>::type                            PointType;
        typedef typename viennagrid::result_of::ncell<Config, CellTag::dim>::type   CellType;

        typedef typename viennagrid::result_of::ncell_range<DomainType, 0>::type                VertexContainer;
        typedef typename viennagrid::result_of::iterator<VertexContainer>::type                     VertexIterator;

        typedef typename viennagrid::result_of::ncell_range<DomainType, CellTag::dim>::type    CellContainer;
        typedef typename viennagrid::result_of::iterator<CellContainer>::type                                 CellIterator;

        typedef typename viennagrid::result_of::ncell_range<CellType, 0>::type                  VertexOnCellContainer;
        typedef typename viennagrid::result_of::iterator<VertexOnCellContainer>::type               VertexOnCellIterator;
        
        typedef typename SystemType::equation_type                  EquationType;
        typedef typename SystemType::equation_type::value_type      Expression;
        
     #ifdef VIENNAFEMDEBUG
        std::cout << "Strong form: " << pde_system.pde(0) << std::endl;
     #endif
        EquationType weak_form_general = viennafem::make_weak_form(pde_system.pde(0));  
     #ifdef VIENNAFEMDEBUG        
        std::cout << "* pde_solver::operator(): Using weak form general: " << weak_form_general << std::endl;
     #endif
        EquationType weak_form = viennamath::apply_coordinate_system(viennamath::cartesian<Config::coordinate_system_tag::dim>(), weak_form_general);

     #ifdef VIENNAFEMDEBUG        
        std::cout << "* pde_solver::operator(): Using weak form " << weak_form << std::endl;
        std::cout << "* pde_solver::operator(): Write dt_dx coefficients" << std::endl;
     #endif
        
        typedef typename reference_cell_for_basis<CellTag, viennafem::lagrange_tag<1> >::type    ReferenceCell;
     
        //fill with cell quantities 
        CellContainer cells = viennagrid::ncells<CellTag::dim>(domain);
        for (CellIterator cell_iter = cells.begin();
            cell_iter != cells.end();
            ++cell_iter)
        {
          //cell_iter->print_short();
          //viennadata::access<example_key, double>()(*cell_iter) = i; 
          viennafem::dt_dx_handler<ReferenceCell>::apply(*cell_iter);
        }

     #ifdef VIENNAFEMDEBUG        
        std::cout << "* pde_solver::operator(): Create Mapping:" << std::endl;
     #endif
        size_t map_index = create_mapping(pde_system, domain);
        
     #ifdef VIENNAFEMDEBUG                
        std::cout << "* pde_solver::operator(): Assigned degrees of freedom in domain so far: " << map_index << std::endl;
     #endif        
        // resize global system matrix and load vector if needed:
        // TODO: This can be a performance bottleneck for large numbers of segments! (lots of resize operations...)
////        if (map_index > system_matrix.size1())
////        {
////          MatrixType temp = system_matrix;
////          ////std::cout << "Resizing system matrix..." << std::endl;
////          system_matrix.resize(map_index, map_index, false);
////          system_matrix.clear();
////          system_matrix.resize(map_index, map_index, false);
////          for (typename MatrixType::iterator1 row_it = temp.begin1();
////               row_it != temp.end1();
////               ++row_it)
////          {
////            for (typename MatrixType::iterator2 col_it = row_it.begin();
////                 col_it != row_it.end();
////                 ++col_it)
////                 system_matrix(col_it.index1(), col_it.index2()) = *col_it;
////          }
////        }
////        if (map_index > load_vector.size())
////        {
////          VectorType temp = load_vector;
////       #ifdef VIENNAFEMDEBUG                          
////          std::cout << "Resizing load vector..." << std::endl;
////       #endif
////          load_vector.resize(map_index, false);
////          load_vector.clear();
////          load_vector.resize(map_index, false);
////          for (size_t i=0; i<temp.size(); ++i)
////            load_vector(i) = temp(i);
////        }

     #ifdef VIENNAFEMDEBUG                        
        std::cout << "* pde_solver::operator(): Transform to reference element" << std::endl;
     #endif
        EquationType transformed_weak_form = viennafem::transform_to_reference_cell<CellType>(weak_form, pde_system);
        
        //std::cout << "* pde_solver::operator(): Transformed weak form:" << std::endl;
        //std::cout << transformed_weak_form << std::endl;
        //std::cout << std::endl;

     #ifdef VIENNAFEMDEBUG                        
        std::cout << "* pde_solver::operator(): Assemble system" << std::endl;
     #endif
        pde_assembler_internal()(transformed_weak_form, pde_system, domain, linsolver);

      }
      
      
      /// specialization for matrix load-vector assembly
      template <typename SystemType, typename DomainType, typename MatrixT, typename VectorT>  //template for operator()
      void operator()(SystemType pde_system,
                      DomainType & domain,
                      MatrixT    & system_matrix, 
                      VectorT    & load_vector
                     ) const
      {
        typedef typename DomainType::config_type              Config;
        typedef typename Config::cell_tag                     CellTag;
        
        typedef typename viennagrid::result_of::point<Config>::type                            PointType;
        typedef typename viennagrid::result_of::ncell<Config, CellTag::dim>::type   CellType;

        typedef typename viennagrid::result_of::ncell_range<DomainType, 0>::type                VertexContainer;
        typedef typename viennagrid::result_of::iterator<VertexContainer>::type                     VertexIterator;

        typedef typename viennagrid::result_of::ncell_range<DomainType, CellTag::dim>::type    CellContainer;
        typedef typename viennagrid::result_of::iterator<CellContainer>::type                                 CellIterator;

        typedef typename viennagrid::result_of::ncell_range<CellType, 0>::type                  VertexOnCellContainer;
        typedef typename viennagrid::result_of::iterator<VertexOnCellContainer>::type               VertexOnCellIterator;
        
        typedef typename SystemType::equation_type                  EquationType;
        typedef typename SystemType::equation_type::value_type      Expression;
        
     #ifdef VIENNAFEMDEBUG
        std::cout << "Strong form: " << pde_system.pde(0) << std::endl;
     #endif
        EquationType weak_form_general = viennafem::make_weak_form(pde_system.pde(0));  
     #ifdef VIENNAFEMDEBUG        
        std::cout << "* pde_solver::operator(): Using weak form general: " << weak_form_general << std::endl;
     #endif
        EquationType weak_form = viennamath::apply_coordinate_system(viennamath::cartesian<Config::coordinate_system_tag::dim>(), weak_form_general);

     #ifdef VIENNAFEMDEBUG        
        std::cout << "* pde_solver::operator(): Using weak form " << weak_form << std::endl;
        std::cout << "* pde_solver::operator(): Write dt_dx coefficients" << std::endl;
     #endif
     
        typedef typename reference_cell_for_basis<CellTag, viennafem::lagrange_tag<1> >::type    ReferenceCell;
        
        //fill with cell quantities 
        CellContainer cells = viennagrid::ncells<CellTag::dim>(domain);
        for (CellIterator cell_iter = cells.begin();
            cell_iter != cells.end();
            ++cell_iter)
        {
          //cell_iter->print_short();
          //viennadata::access<example_key, double>()(*cell_iter) = i; 
          viennafem::dt_dx_handler<ReferenceCell>::apply(*cell_iter);
        }

     #ifdef VIENNAFEMDEBUG        
        std::cout << "* pde_solver::operator(): Create Mapping:" << std::endl;
     #endif
        size_t map_index = create_mapping(pde_system, domain);
        
     #ifdef VIENNAFEMDEBUG                
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
       #ifdef VIENNAFEMDEBUG                          
          std::cout << "Resizing load vector..." << std::endl;
       #endif
          load_vector.resize(map_index, false);
          load_vector.clear();
          load_vector.resize(map_index, false);
          for (size_t i=0; i<temp.size(); ++i)
            load_vector(i) = temp(i);
        }

     #ifdef VIENNAFEMDEBUG                        
        std::cout << "* pde_solver::operator(): Transform to reference element" << std::endl;
     #endif
        EquationType transformed_weak_form = viennafem::transform_to_reference_cell<CellType>(weak_form, pde_system);
        
        std::cout << "* pde_solver::operator(): Transformed weak form:" << std::endl;
        std::cout << transformed_weak_form << std::endl;
        //std::cout << std::endl;

     #ifdef VIENNAFEMDEBUG                        
        std::cout << "* pde_solver::operator(): Assemble system" << std::endl;
     #endif
     
        typedef equation_wrapper<MatrixT, VectorT>    wrapper_type;
        wrapper_type wrapper(system_matrix, load_vector);
     
        pde_assembler_internal()(transformed_weak_form, pde_system, domain, wrapper);
//        pde_assembler_internal()(transformed_weak_form, pde_system, domain, system_matrix, load_vector);

      }
      
  };

}
#endif
