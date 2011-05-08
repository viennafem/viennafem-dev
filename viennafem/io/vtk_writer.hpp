/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */


#ifndef VIENNAFEM_IO_VTKWRITER_HPP
#define VIENNAFEM_IO_VTKWRITER_HPP

// include necessary system headers
#include <iostream>

#include "viennafem/forwards.h"

// ViennaGrid includes:
#include "viennagrid/domain.hpp"
#include "viennagrid/io/sgf_reader.hpp"
#include "viennagrid/io/vtk_writer.hpp"

// ViennaData includes:
#include "viennadata/api.hpp"

namespace viennafem
{
  namespace io
  {

    template <typename VectorType,
              typename DomainType>
    void write_solution_to_VTK_file(VectorType const & result,
                                    std::string filename,
                                    DomainType const & domain,
                                    long id)
    {
      typedef typename viennagrid::result_of::const_ncell_container<DomainType, 0>::type    VertexContainer;
      typedef typename viennagrid::result_of::iterator<VertexContainer>::type         VertexIterator;
      
      typedef viennafem::mapping_key          MappingKeyType;
      typedef viennafem::boundary_key         BoundaryKeyType;
      
      MappingKeyType map_key(id);
      BoundaryKeyType bnd_key(id);
      
      
      std::cout << "* write_solution_to_VTK_file(): Writing result on mesh for later export" << std::endl;
      VertexContainer vertices = viennagrid::ncells<0>(domain);
      for (VertexIterator vit = vertices.begin();
          vit != vertices.end();
          ++vit)
      {
        long cur_index = viennadata::access<MappingKeyType, long>(map_key)(*vit);
        if (cur_index > -1)
          viennadata::access<std::string, double>("vtk_data")(*vit) = result[cur_index];
        else //use Dirichlet boundary data:
          viennadata::access<std::string, double>("vtk_data")(*vit) = 
            viennadata::access<BoundaryKeyType, double>(bnd_key)(*vit);
      }

      std::cout << "* write_solution_to_VTK_file(): Writing data to '"
                << filename
                << "' (can be viewed with e.g. Paraview)" << std::endl;

      viennagrid::io::vtk_writer<DomainType> my_vtk_writer;
      my_vtk_writer.writeDomain(domain, filename);  
    }


  }
}
#endif

