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
#include "viennafem/transform/dtdx_triangle.hpp"
#include "viennafem/transform/dtdx_tetrahedron.hpp"
#include "viennafem/weak_form.hpp"
#include "viennafem/linear_pde_system.hpp"
#include "viennafem/linear_pde_options.hpp"
#include "viennafem/assembler.hpp"


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
      
      template <typename SystemType, typename DomainType, typename MatrixType, typename VectorType>  //template for operator()
      void operator()(SystemType pde_system,
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
        BoundaryKeyType bnd_key(pde_system.option(0).data_id());
        MappingKeyType map_key(pde_system.option(0).data_id());
        
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
            viennadata::access<MappingKeyType, long>(map_key)(*vit) = map_index;
            map_index += pde_system.unknown(0).size();
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
        
        EquationType transformed_weak_form = viennafem::transform_to_reference_cell<CellType>(weak_form, pde_system);
        
        std::cout << "* pde_solver::operator(): Transformed weak form:" << std::endl;
        std::cout << transformed_weak_form << std::endl;
        std::cout << std::endl;

        std::cout << "* pde_solver::operator(): Assemble system" << std::endl;
        pde_assembler()(transformed_weak_form, pde_system, domain, system_matrix, load_vector);

      }
      
  };

}
#endif
