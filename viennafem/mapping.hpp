#ifndef VIENNAFEM_MAPPING_HPP
#define VIENNAFEM_MAPPING_HPP

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


#include <vector>
#include "viennafem/forwards.h"
#include "viennagrid/config/domain_config.hpp"


/** @file   mapping.hpp
    @brief  Sets up and handles the local-to-global map (handles degrees of freedom).
*/

namespace viennafem
{
  
  namespace detail
  {
    /** @brief If EntityType is a ViennaGrid segment, returns the domain. If EntityType is already the domain, no changes.
    */
    template <typename EntityType>
    struct extract_domain
    {
      typedef EntityType  type;
      static EntityType & apply(EntityType & domain) { return domain; }
    };
    
    /** @brief Specialization of the domain extraction for a ViennaGrid segment */
    template <typename ConfigType>
    struct extract_domain<viennagrid::segment_t<ConfigType> >
    {
      typedef typename viennagrid::segment_t<ConfigType>::domain_type    type;
      static type & apply(viennagrid::segment_t<ConfigType> & seg) { return seg.segmentation().domain(); }
    };
  }
  
  
  /** @brief Distributes mapping indices over domain or segment.
   * 
   */
  template <typename StorageType, typename SystemType, typename DomainType>
  long create_mapping(StorageType& storage, 
                      SystemType & pde_system,
                      DomainType & domain)
  {
    typedef typename DomainType::config_type              Config;
    typedef typename Config::cell_tag                     CellTag;
    
    typedef typename viennagrid::result_of::point<Config>::type                                       PointType;
    typedef typename viennagrid::result_of::element<DomainType, viennagrid::vertex_tag>::type         VertexType;    
    typedef typename viennagrid::result_of::element<DomainType, CellTag>::type                        CellType;

    typedef typename viennagrid::result_of::element_range<DomainType, viennagrid::vertex_tag>::type   VertexRange;
    typedef typename viennagrid::result_of::iterator<VertexRange>::type                               VertexIterator;

    typedef typename SystemType::mapping_key_type   MappingKeyType;
    typedef typename SystemType::boundary_key_type  BoundaryKeyType;
    
    BoundaryKeyType bnd_key(pde_system.option(0).data_id());
    MappingKeyType  map_key(pde_system.option(0).data_id());
    
    long start_index = viennadata::access<MappingKeyType, long>(storage, map_key, detail::extract_domain<DomainType>::apply(domain));
    long map_index = start_index;
    bool init_done = viennadata::access<MappingKeyType, bool>(storage, map_key, detail::extract_domain<DomainType>::apply(domain));

    //eventually, map indices need to be set to invalid first:
    if (!init_done)
    {
      typedef typename detail::extract_domain<DomainType>::type                                         TrueDomainType;
      typedef typename viennagrid::result_of::element_range<DomainType, viennagrid::vertex_tag>::type   DomainVertexContainer;
      typedef typename viennagrid::result_of::iterator<DomainVertexContainer>::type                     DomainVertexIterator;
      
      DomainVertexContainer vertices = viennagrid::elements<VertexType>(detail::extract_domain<DomainType>::apply(domain));  
      for (DomainVertexIterator vit = vertices.begin();
          vit != vertices.end();
          ++vit)
      {  
        viennadata::access<MappingKeyType, long>(storage, map_key, *vit) = -1;
      }
      viennadata::access<MappingKeyType, bool>(storage, map_key, detail::extract_domain<DomainType>::apply(domain)) = true;
    }
    
    VertexRange vertices = viennagrid::elements<VertexType>(domain);  
    for (VertexIterator vit = vertices.begin();
        vit != vertices.end();
        ++vit)
    {  
      if (viennadata::access<BoundaryKeyType, bool>(storage, bnd_key, *vit))
      {
        //std::cout << "boundary vertex" << std::endl;
        viennadata::access<MappingKeyType, long>(storage, map_key, *vit) = -1;
      }
      else
      {
        //std::cout << "interior vertex" << std::endl;
        if (viennadata::access<MappingKeyType, long>(storage, map_key, *vit) < 0) //only assign if no dof assigned yet
        {
          viennadata::access<MappingKeyType, long>(storage, map_key, *vit) = map_index;
          map_index += pde_system.unknown(0).size();
        }
        //else
        //  std::cout << "Found vertex with DOF!" << std::endl;
      }
    }
    
    viennadata::access<MappingKeyType, long>(storage, map_key, detail::extract_domain<DomainType>::apply(domain)) = map_index;
    
    return map_index;
  }
  
  
  

  /** @brief Returns an array of mapping indices for the provided cell.
   */
  template <typename StorageType, typename SystemType, typename CellType>
  std::vector<long> mapping_indices(StorageType& storage, SystemType & pde_system, CellType const & cell, std::size_t pde_id = 0)
  {
    typedef typename CellType::config_type              Config;
    typedef typename Config::cell_tag                     CellTag;

    typedef typename viennagrid::result_of::element<Config, viennagrid::vertex_tag>::type            VertexType;
    typedef typename viennagrid::result_of::element_range<CellType, viennagrid::vertex_tag>::type    VertexOnCellContainer;
    typedef typename viennagrid::result_of::iterator<VertexOnCellContainer>::type                    VertexOnCellIterator;
    
//    typedef typename viennagrid::result_of::const_ncell_range<CellType, 0>::type                 VertexOnCellContainer;
//    typedef typename viennagrid::result_of::iterator<VertexOnCellContainer>::type               VertexOnCellIterator;

    typedef typename SystemType::mapping_key_type   MappingKeyType;
    typedef std::vector<long>                                   MappingContainer;
    
   
    MappingKeyType map_key(pde_system.option(pde_id).data_id());

    VertexOnCellContainer vertices_on_cell = viennagrid::elements<VertexType>(cell);  
    
    std::size_t unknown_components = pde_system.unknown(pde_id).size();
    //std::cout << "* mapping_indices() with space_id = " << space_id << " and unknown_components = " << unknown_components << std::endl;
    
    MappingContainer ret(viennagrid::boundary_elements<CellTag, viennagrid::vertex_tag>::num * unknown_components);
    
    
    long local_index = 0;
    for (VertexOnCellIterator vocit = vertices_on_cell.begin();
                              vocit != vertices_on_cell.end();
                              ++vocit)
    {
      long map_base = viennadata::access<MappingKeyType, long>(storage, map_key, *vocit);
      
      if (map_base < 0) //Dirichlet boundary
      {
        for (std::size_t i=0; i<unknown_components; ++i)
          ret[local_index++] = map_base;
      }
      else
      {
        for (std::size_t i=0; i<unknown_components; ++i)
          ret[local_index++] = map_base + i;
      }
    }
   
    return ret; //TODO: Avoid temporary
  }   
   

}

#endif
