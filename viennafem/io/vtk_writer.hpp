#ifndef VIENNAFEM_IO_VTKWRITER_HPP
#define VIENNAFEM_IO_VTKWRITER_HPP

/* =========================================================================
   Copyright (c) 2012-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
                             -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                             -----------------

   Author:     Karl Rupp                          rupp@iue.tuwien.ac.at

   License:    MIT (X11), see file LICENSE in the ViennaFEM base directory
============================================================================ */


// include necessary system headers
#include <iostream>

#include "viennafem/forwards.h"

// ViennaGrid includes:
#include "viennagrid/io/vtk_writer.hpp"

// ViennaData includes:
#include "viennadata/api.hpp"

namespace viennafem
{
  namespace io
  {

    template <typename VectorType,
              typename DomainType,
              typename SegmentationType,
              typename StorageType>
    void write_solution_to_VTK_file(VectorType const & result,
                                    std::string filename,
                                    DomainType const & domain,
                                    SegmentationType const & segmentation,
                                    StorageType const & storage,
                                    std::vector<long> id_vector)
    {
      typedef typename viennagrid::result_of::element<DomainType, viennagrid::vertex_tag>::type               VertexType;
      typedef typename viennagrid::result_of::const_element_range<DomainType, viennagrid::vertex_tag>::type   VertexContainer;
      typedef typename viennagrid::result_of::iterator<VertexContainer>::type                                 VertexIterator;

      typedef viennafem::mapping_key          MappingKeyType;
      typedef viennafem::boundary_key         BoundaryKeyType;

      std::cout << "* write_solution_to_VTK_file(): Writing result on mesh for later export" << std::endl;
      viennagrid::io::vtk_writer<DomainType> my_vtk_writer;

      std::map< std::string, std::deque<double> > output_values;
      for (std::size_t i=0; i<id_vector.size(); ++i)
      {
        long id = id_vector[i];

        MappingKeyType  map_key(id);
        BoundaryKeyType bnd_key(id);

        typename viennadata::result_of::accessor<const StorageType, viennafem::mapping_key, long, VertexType>::type   mapping_accessor =
          viennadata::make_accessor<viennafem::mapping_key, long, VertexType>(storage, map_key);

        typename viennadata::result_of::accessor<const StorageType, BoundaryKeyType, double, VertexType>::type        boundary_accessor =
          viennadata::make_accessor<BoundaryKeyType, double, VertexType>(storage, bnd_key);

        std::stringstream ss;
        ss << "fem_result" << id;
        std::string result_string = ss.str(); // also used for passing staff over to VTK

        typename viennagrid::result_of::accessor< std::deque<double>, VertexType >::type output_value_accessor( output_values[result_string] );

        VertexContainer vertices = viennagrid::elements<VertexType>(domain);
        for (VertexIterator vit = vertices.begin();
                            vit != vertices.end();
                          ++vit)
        {
          long cur_index = mapping_accessor(*vit);
          if (cur_index > -1)
            output_value_accessor(*vit) = result[cur_index];
          else //use Dirichlet boundary data:
          {
            // TODO if Accessor concept takes care of that -> change!
            if (boundary_accessor.find(*vit))
              output_value_accessor(*vit) = boundary_accessor(*vit);
            else
              output_value_accessor(*vit) = false;
          }
          //output_value_accessor(*vit) = boundary_accessor(*vit);
        }

        my_vtk_writer.add_scalar_data_on_vertices( output_value_accessor, result_string );
      }

      std::cout << "* write_solution_to_VTK_file(): Writing data to '"
                << filename
                << "' (can be viewed with e.g. Paraview)" << std::endl;

      my_vtk_writer(domain, segmentation, filename);
    }

    template <typename VectorType,
              typename DomainType,
              typename SegmentationType,
              typename StorageType>
    void write_solution_to_VTK_file(VectorType const & result,
                                    std::string filename,
                                    DomainType const & domain,
                                    SegmentationType const & segmentation,
                                    StorageType const & storage,
                                    long id)
    {
      std::vector<long> id_vector(1);
      id_vector[0] = id;

      write_solution_to_VTK_file(result, filename, domain, segmentation, storage, id_vector);
    }
  }
}
#endif

