#ifndef VIENNAFEM_IO_VTKWRITER_HPP
#define VIENNAFEM_IO_VTKWRITER_HPP

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
//       typedef typename DomainType::config_type                                              ConfigType;
//       typedef typename ConfigType::cell_tag                                                 CellTag;

      typedef typename viennagrid::result_of::cell_tag<DomainType>::type CellTag;

      typedef typename viennagrid::result_of::element<DomainType, CellTag>::type         CellType;

      typedef typename viennagrid::result_of::const_element_range<DomainType, CellTag>::type   CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

      typedef viennafem::mapping_key          MappingKeyType;
      typedef viennafem::boundary_key         BoundaryKeyType;

      std::cout << "* write_solution_to_VTK_file(): Writing result on mesh for later export" << std::endl;
      viennagrid::io::vtk_writer<DomainType> my_vtk_writer;



      std::map< std::string, std::deque<double> > output_values;

      for (std::size_t i=0; i<id_vector.size(); ++i)
      {
        long id = id_vector[i];
        MappingKeyType map_key(id);
        BoundaryKeyType bnd_key(id);

        typename viennadata::result_of::accessor<const StorageType, viennafem::mapping_key, long, CellType>::type cell_mapping_accessor =
          viennadata::make_accessor<viennafem::mapping_key, long, CellType>(storage, map_key);
        
        typename viennadata::result_of::accessor<const StorageType, BoundaryKeyType, double, CellType>::type boundary_accessor =
          viennadata::make_accessor<BoundaryKeyType, double, CellType>(storage, bnd_key);

        std::stringstream ss;
        ss << "fvm_result" << id;
        std::string result_string = ss.str(); // also used for passing staff over to VTK

        typename viennagrid::result_of::accessor< std::deque<double>, CellType >::type output_value_accessor( output_values[result_string] );

        CellContainer cells = viennagrid::elements(domain);
        for (CellIterator cit = cells.begin();
                          cit != cells.end();
                        ++cit)
        {
          long cur_index = cell_mapping_accessor(*cit);
          if (cur_index > -1)
            output_value_accessor(*cit) = result[cur_index];
          else //use Dirichlet boundary data:
          {
            // TODO if Accessor concept takes care of that -> change!
            if (boundary_accessor.find(*cit))
              output_value_accessor(*cit) = boundary_accessor(*cit);
            else
              output_value_accessor(*cit) = false;
          }
//             output_value_accessor(*cit) = boundary_accessor(*cit);
        }

        my_vtk_writer.add_scalar_data_on_cells( output_value_accessor, result_string );
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

