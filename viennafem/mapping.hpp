/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */


#ifndef VIENNAFEM_MAPPING_HPP
#define VIENNAFEM_MAPPING_HPP

#include <vector>
#include "viennafem/forwards.h"
#include "viennafem/bases/BFStock.hpp"
#include "viennagrid/forwards.h"
#include "viennagrid/domain.hpp"
#include "viennagrid/segment.hpp"
#include "viennagrid/iterators.hpp"

namespace viennafem
{
  
  /** @brief If EntityType is a ViennaGrid segment, returns the domain. If EntityType is already the domain, no changes.
   */
  template <typename EntityType>
  struct extract_domain
  {
    typedef EntityType  type;
    static EntityType & apply(EntityType & domain) { return domain; }
  };
  
  template <typename ConfigType>
  struct extract_domain<viennagrid::segment_t<ConfigType> >
  {
    typedef viennagrid::domain<ConfigType> type;
    static type & apply(viennagrid::segment_t<ConfigType> & seg) { return seg.get_domain(); }
  };
  
  
  
  /** @brief Distributes mapping indices over domain or segment
   * 
   */
  template <typename SystemType, typename DomainType>
  long create_mapping(SystemType pde_system,
                      DomainType & domain)
  {
    typedef typename DomainType::config_type              Config;
    typedef typename Config::cell_tag                     CellTag;
    
    typedef typename viennagrid::result_of::point_type<Config>::type                            PointType;
    typedef typename viennagrid::result_of::ncell_type<Config, CellTag::topology_level>::type   CellType;

    typedef typename viennagrid::result_of::ncell_container<DomainType, 0>::type                VertexContainer;
    typedef typename viennagrid::result_of::iterator<VertexContainer>::type                     VertexIterator;

    typedef viennafem::boundary_key                             BoundaryKeyType;
    typedef viennafem::mapping_key                              MappingKeyType;
    
    BoundaryKeyType bnd_key(pde_system.option(0).data_id());
    MappingKeyType map_key(pde_system.option(0).data_id());
    
    long start_index = viennadata::access<MappingKeyType, long>(map_key)(extract_domain<DomainType>::apply(domain));
    long map_index = start_index;
    bool init_done = viennadata::access<MappingKeyType, bool>(map_key)(extract_domain<DomainType>::apply(domain));

    //eventually, map indices need to be set to invalid first:
    if (!init_done)
    {
      typedef typename extract_domain<DomainType>::type   TrueDomainType;
      typedef typename viennagrid::result_of::ncell_container<TrueDomainType, 0>::type            DomainVertexContainer;
      typedef typename viennagrid::result_of::iterator<DomainVertexContainer>::type               DomainVertexIterator;
      
      DomainVertexContainer vertices = viennagrid::ncells<0>(extract_domain<DomainType>::apply(domain));
      for (DomainVertexIterator vit = vertices.begin();
          vit != vertices.end();
          ++vit)
      {  
        viennadata::access<MappingKeyType, long>(map_key)(*vit) = -1;
      }
      viennadata::access<MappingKeyType, bool>(map_key)(extract_domain<DomainType>::apply(domain)) = true;
    }
    
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
        if (viennadata::access<MappingKeyType, long>(map_key)(*vit) < 0) //only assign if no dof assigned yet
          viennadata::access<MappingKeyType, long>(map_key)(*vit) = map_index;
        //else
        //  std::cout << "Found vertex with DOF!" << std::endl;
        map_index += pde_system.unknown(0).size();
      }
    }
    
    viennadata::access<MappingKeyType, long>(map_key)(extract_domain<DomainType>::apply(domain)) = map_index;
    
    return map_index;
  }
  
  
  

  /** @brief Returns an array of mapping indices 
   *  @param cell      The cell for which the mapping indices should be returned
   *  @param space_id  Accessor ID for retrieving the mapping key
   */
  template <typename CellType>
  std::vector<long> mapping_indices(CellType const & cell, long space_id, long unknown_components)
  {
    typedef typename CellType::config_type              Config;
    typedef typename Config::cell_tag                     CellTag;
    
    typedef typename viennagrid::result_of::const_ncell_container<CellType, 0>::type                  VertexOnCellContainer;
    typedef typename viennagrid::result_of::iterator<VertexOnCellContainer>::type               VertexOnCellIterator;

    typedef viennafem::mapping_key                              MappingKeyType;
    typedef std::vector<long>                                   MappingContainer;
    
    //std::cout << "* mapping_indices() with space_id = " << space_id << " and unknown_components = " << unknown_components << std::endl;
    
    MappingKeyType map_key(space_id);

    VertexOnCellContainer vertices_on_cell = viennagrid::ncells<0>(cell);
    
    MappingContainer ret(viennagrid::traits::subcell_desc<CellTag, 0>::num_elements * unknown_components);
    
    long local_index = 0;
    for (VertexOnCellIterator vocit = vertices_on_cell.begin();
                              vocit != vertices_on_cell.end();
                              ++vocit)
    {
      long map_base = viennadata::access<MappingKeyType, long>(map_key)(*vocit);
      
      if (map_base < 0) //Dirichlet boundary
      {
        for (long i=0; i<unknown_components; ++i)
          ret[local_index++] = map_base;
      }
      else
      {
        for (long i=0; i<unknown_components; ++i)
          ret[local_index++] = map_base + i;
      }
    }
   
    return ret; //TODO: Avoid temporary
  }
  
  
/*  
  
  //computeIndexFromBFKey: 
  template <typename BFKey, typename Permutator>
  long computeIndexFromBFKey_impl(BFKey const & bf, Permutator const & perm, long startindex, long degree, long tuple_len)
  {
    long ret = 0;
    long i_k = bf[perm(startindex)];
//     std::cout << "startindex: " << startindex << std::endl;
//     std::cout << "ubound: " << perm(startindex) << std::endl;
//     std::cout << "bf[ubound]: " << bf[perm(startindex)] << std::endl;
//     std::cout << "degree: " << degree << std::endl;
//     std::cout << "tuple_len: " << tuple_len << std::endl;

    for (long j=1; j<i_k; ++j)
      ret += PascalSimplex<0,0>::getInstance().getNum( degree - j, tuple_len - 1);

    //if index-tuple has still more than two entries, continue recursion.
    //if two entries, the final index was just computed, cause a one-tuple has only one possibility to distribute degrees.
    if (tuple_len > 2)
      ret += computeIndexFromBFKey_impl(bf, perm, startindex + 1,
                                        degree - i_k, tuple_len - 1);

    return ret;
  }

  template <typename BFKey, typename Permutator>
  long computeIndexFromBFKey(BFKey const & bfkey, Permutator const & perm)
  {
    return computeIndexFromBFKey_impl(bfkey, perm, 0, bfkey.getDegree(), bfkey.size());
  }

  //for Gamma<id>, Interface<id>, the assembly-cell-type is the domain facet-type
  template <typename DomainConfig, typename IntegrationDomain>
  struct ASSEMBLY_CELL_TYPE
  {
    typedef typename viennagrid::DomainTypes<DomainConfig>::FacetType    ResultType;
  };

  //for Omega, the assembly-cell-type is the domain cell-type
  template <typename DomainConfig>
  struct ASSEMBLY_CELL_TYPE<DomainConfig, Omega >
  {
    typedef typename viennagrid::DomainTypes<DomainConfig>::CellType    ResultType;
  };

  struct ElementOrientation {}; //to be removed again

  //Constructor helper for MappingIterator:
  template <typename IntegrationDomain>
  struct TwinHandler
  {
    //default case (Omega): pass arguments through:
    template <typename CellType, typename ElementType>
    static ElementType & getElement(CellType & cell, ElementType & elem)
    {
      return elem;
    }
    
    template <typename CellType, typename ElementType>
    static ElementOrientation const & getPermutation(CellType & cell, ElementType & elem, long orientation_index)
    {
      typedef typename ElementType::ElementTag                 ElementTag;
      return cell.template getLevelOrientation<ElementTag::TopoLevel>(orientation_index);
    }

    template <typename CellType, typename ElementType, typename MappingKey, typename SegmentType>
    static const std::vector<long> & getMappingVector(CellType & cell, ElementType & elem, MappingKey, SegmentType)
    {
      return elem.template retrieveQuantity< std::vector<long> >(MappingKey());
    }
  };


  //Iterate over mapping numbers on a given element:
  template <typename FEMConfig, typename ElementType, typename IntegrationDomain,
              bool map_u = false, //mapping for basefun<2>?
              long mapnum = PascalSimplex<FEMConfig::BasisfunctionTag::degree,
                                          viennagrid::subcell_traits<typename ElementType::ElementTag,
                                                        0>::ElementNum
                                         >::ReturnValue,
              bool cellMapping = true                          
              //                           
              bool cellMapping = (EQUALS<typename ASSEMBLY_CELL_TYPE
                                                     <typename ElementType::Configuration,
                                                      IntegrationDomain>::ResultType,
                                         ElementType>::ReturnValue > 0)
                                  && (mapnum != 1)
                                  && (NOT_OF_INTERFACE_TYPE<IntegrationDomain>::ReturnValue > 0 )//
             >
  class MappingIterator_impl
  {
      typedef typename ElementType::ElementTag                ElementTag;
      typedef typename ElementType::Configuration             DomainConfig;
      typedef typename ASSEMBLY_CELL_TYPE<DomainConfig, IntegrationDomain>::ResultType
                                                                 CellType;
      typedef typename viennagrid::DomainTypes<DomainConfig>::SegmentType SegmentType;

      typedef typename FEMConfig::MappingKey                  MappingKey;
      typedef typename FEMConfig::BoundaryData                BoundaryData;
      typedef typename FEMConfig::BoundaryData2               BoundaryData2;
      typedef typename FEMConfig::BoundaryData3               BoundaryData3;

      typedef typename viennagrid::IteratorTypes<CellType, ElementTag::TopoLevel>::ResultType    LevelIterator;

    public:
      //standard constructor:
      MappingIterator_impl(CellType & cell, ElementType & levelelement, long orientation_index) :
        curEl(TwinHandler<IntegrationDomain>::getElement(cell, levelelement)),
        currentVec(TwinHandler<IntegrationDomain>::getMappingVector(cell, curEl, MappingKey(), SegmentType())),
        retval(0),
        permutator( TwinHandler<IntegrationDomain>::getPermutation(cell, curEl, orientation_index)),
        bfkey(viennagrid::subcell_traits<ElementTag, 0>::ElementNum, FEMConfig::BasisfunctionTag::degree)
      {
        retval = currentVec[computeIndexFromBFKey(bfkey, permutator)];
      };

      long operator*() const
      {
        //return currentVec[computeIndexFromBFKey(bfkey, permutator)];
        return retval;
      }

      MappingIterator_impl & operator++()
      {
        //std::cout << "Incrementing edge-mapping" << std::endl;
        ++bfkey;
        retval = currentVec[computeIndexFromBFKey(bfkey, permutator)];
        return *this;
      }

      bool valid() const { return bfkey.valid(); }
//       bool operator==(bool b) { return Base::isValid == b; }
//       bool operator!=(bool b) { return Base::isValid != b; }

      double getBoundaryValue(int component)
      {
        //std::cout << "getting boundary value" << std::endl;
        if (component == 0)
          return (curEl.template retrieveQuantity< double >( BoundaryData() ));
        if (component == 1)
          return (curEl.template retrieveQuantity< double >( BoundaryData2() ));
        if (component == 2)
          return (curEl.template retrieveQuantity< double >( BoundaryData3() ));
        return 0.0;
      }

    private:
      ElementType & curEl;
      const std::vector<long> & currentVec;
      long retval;
      ElementOrientation const & permutator;
      BaseFunKey bfkey;
  };

  //Iterate over a single number per element (no interface)
  //performance gain: iterate over elements with direct map-index access
  //instead of iteration over elements plus iteration over MappingIterator_impl
  template <typename FEMConfig, typename ElementType, typename IntegrationDomain, bool map_u>
  class MappingIterator_impl<FEMConfig, ElementType, IntegrationDomain, map_u, 1, false>
  {
      typedef typename ElementType::ElementTag                ElementTag;
      typedef typename ElementType::Configuration             DomainConfig;
      typedef typename ASSEMBLY_CELL_TYPE<DomainConfig, IntegrationDomain>::ResultType
                                                                 CellType;

      typedef typename FEMConfig::MappingKey                  MappingKey;
      typedef typename FEMConfig::BoundaryData                BoundaryData;
      typedef typename FEMConfig::BoundaryData2               BoundaryData2;
      typedef typename FEMConfig::BoundaryData3               BoundaryData3;

      typedef typename viennagrid::IteratorTypes<CellType, ElementTag::TopoLevel>::ResultType    LevelIterator;

    public:
      MappingIterator_impl(CellType & cell, ElementType & levelelement, long orientation) : 
        lit(cell.template getLevelIteratorBegin<ElementTag::TopoLevel>()),
        litend(cell.template getLevelIteratorEnd<ElementTag::TopoLevel>()),
        isValid(true),
        retval(lit->template retrieveQuantity<long>(MappingKey()))
      {
        //std::cout << "No Interface, single " << std::endl;
      };

      long operator*() const
      {
        //return curEl.template retrieveQuantity<long>(MappingKey());
        return retval;
      };

      MappingIterator_impl & operator++()
      {
        if (++lit == litend)
          isValid = false;
        else
          retval = lit->template retrieveQuantity<long>(MappingKey());

        return *this;
      }

      bool valid() const { return isValid; }

      double getBoundaryValue(int component)
      {
        if (component == 0)
          return (lit->template retrieveQuantity< double >( BoundaryData() ));
        if (component == 1)
          return (lit->template retrieveQuantity< double >( BoundaryData2() ));
        if (component == 2)
          return (lit->template retrieveQuantity< double >( BoundaryData3() ));
        return 0.0;
      }

    private:
      LevelIterator lit;
      LevelIterator litend;
      bool isValid;
      long retval;
  };

  //Iterate over a single number per element on an interface
  //performance gain: iterate over elements with direct map-index access
  //instead of iteration over elements plus iteration over MappingIterator_impl
  
  // *
  template <typename FEMConfig, typename ElementType, long id>
  class MappingIterator_impl<FEMConfig, ElementType, Interface<id>, true, 1, false>
  {
      typedef typename ElementType::ElementTag                ElementTag;
      typedef typename ElementType::Configuration             DomainConfig;
      typedef typename ASSEMBLY_CELL_TYPE<DomainConfig, Interface<id> >::ResultType  CellType;
      typedef typename viennagrid::DomainTypes<DomainConfig>::SegmentType SegmentType;

      typedef typename FEMConfig::MappingKey                  MappingKey;
      typedef typename FEMConfig::BoundaryData                BoundaryData;
      typedef typename FEMConfig::BoundaryData2               BoundaryData2;
      typedef typename FEMConfig::BoundaryData3               BoundaryData3;

      typedef typename IteratorTypes<CellType, ElementTag::TopoLevel>::ResultType    LevelIterator;

    public:
      MappingIterator_impl(CellType & cell, ElementType & levelelement, long orientation) : 
        cell_(cell),
        lit(cell.template getLevelIteratorBegin<ElementTag::TopoLevel>()),
        litend(cell.template getLevelIteratorEnd<ElementTag::TopoLevel>()),
        isValid(true)
      {
        //std::cout << "Start: "; lit->print(); std::cout << std::endl;
        setRetval();
      };

      long operator*() const
      {
        //return curEl.template retrieveQuantity<long>(MappingKey());
        return retval;
      };

      MappingIterator_impl & operator++()
      {
        if (++lit == litend)
          isValid = false;
        else
          setRetval();

        return *this;
      }

      bool valid() const { return isValid; }

      //Dirichlet boundary is meaningless for interfaces:
      double getBoundaryValue(int component) { return 0.0; }

    private:

      void setRetval()
      {
        ElementType & twinElement = TwinHandler< Interface<id> >::getElement(cell_, *lit);

        //toggle active segment:
        CellType * twinCell = cell_.template retrieveQuantity<CellType *>(Interface<id>());
        SegmentType * twinSeg = cell_.template retrieveQuantity<SegmentType *>(Interface<id>());
        cell_.setCurrentSegment(*twinSeg);

        retval = twinElement.template retrieveQuantity<long>(MappingKey());

        //toggle active segment again:
        SegmentType * seg = twinCell->template retrieveQuantity<SegmentType *>(Interface<id>());
        twinCell->setCurrentSegment(*seg);
      }

      CellType & cell_;
      LevelIterator lit;
      LevelIterator litend;
      bool isValid;
      long retval;
  }; * //

  //iterate on cell only (bubble-functions, no permutator needed)
  template <typename FEMConfig, typename ElementType, typename IntegrationDomain, bool map_u, long mapnum>
  class MappingIterator_impl<FEMConfig, ElementType, IntegrationDomain, map_u, mapnum, true>
  {
      typedef typename ElementType::Configuration             DomainConfig;
      typedef typename ASSEMBLY_CELL_TYPE<DomainConfig, IntegrationDomain>::ResultType
                                                                 CellType;

      typedef typename ElementType::ElementTag                ElementTag;

      typedef typename FEMConfig::MappingKey                  MappingKey;
      typedef typename FEMConfig::BoundaryData                BoundaryData;
      typedef typename FEMConfig::BoundaryData2               BoundaryData2;
      typedef typename FEMConfig::BoundaryData3               BoundaryData3;

      typedef typename viennagrid::IteratorTypes<CellType, ElementTag::TopoLevel>::ResultType    LevelIterator;

    public:
      MappingIterator_impl(CellType & cell, ElementType & levelelement, long orientation) :
        currentVec(&(levelelement.template retrieveQuantity< std::vector<long> >(MappingKey()))),
        isValid(true),
        index_(0)
      {
        //std::cout << "Cell..." << std::endl;
      };

      long operator*() const
      {
        return (*currentVec)[index_];
      }

      MappingIterator_impl & operator++()
      {
        if (++index_ == currentVec->size())
          isValid = false;
        return *this;
      }

      bool valid() const { return isValid; }
//       bool operator==(bool b) { return Base::isValid == b; }
//       bool operator!=(bool b) { return Base::isValid != b; }

      //Dirichlet boundary is meaningless for interfaces:
      double getBoundaryValue(int component) { return 0.0; }

    private:
      const std::vector<long> * currentVec;
      bool isValid;
      unsigned long index_;
  };

  ////////////////// All-In-One Mapping Iterator ////////////////////
  template <typename FEMConfig, typename CellType, typename IntegrationDomain, bool map_u,
              long topolevel = 0,
              long mapnum = 3,
//              BasisFuncNum<typename CellType::ElementTag,
//                                         typename FEMConfig::BasisfunctionTag,
//                                         topolevel>::ReturnValue,
              bool topolevel_valid = (CellType::ElementTag::TopoLevel - topolevel >= 0)
              >
  struct MappingFiller
  {
    typedef typename CellType::Configuration             DomainConfig;
//     typedef typename ASSEMBLY_CELL_TYPE<DomainConfig, IntegrationDomain>::ResultType
//                                                                 CellType;
    typedef typename viennagrid::subcell_traits<typename CellType::ElementTag,
                                                topolevel>::ElementTag              ElementTag;
    typedef viennagrid::element<DomainConfig, ElementTag>                           ElementType;

    typedef typename viennagrid::IteratorTypes<CellType, topolevel>::ResultType  LevelIterator;

    //default case: iterate over elements at this level:
    static void apply(CellType & cell, long * mappingList, long mapindex)
    {
      //fill current topolevel:
      long i=0;
      long index = mapindex;
      for (LevelIterator lit = cell.template getLevelIteratorBegin<topolevel>();
            lit != cell.template getLevelIteratorEnd<topolevel>();
            ++lit, ++i)
      {
        MappingIterator_impl<FEMConfig, ElementType, IntegrationDomain, map_u> mapit(cell, *lit, i);
        while (mapit.valid())
        {
          mappingList[index] = (*mapit);
          ++index; ++mapit;
        }
      }

      //go to next topolevel:
      MappingFiller<FEMConfig, CellType, IntegrationDomain, map_u, topolevel + 1>::apply(cell, mappingList, index);
    }
  };


  //if only one mapping number per element, do not iterate over elements.
  template <typename FEMConfig, typename CellType, typename IntegrationDomain, bool map_u,
              long topolevel>
  struct MappingFiller<FEMConfig, CellType, IntegrationDomain, map_u, topolevel, 1, true>
  {
    typedef typename CellType::Configuration             DomainConfig;
    typedef typename viennagrid::subcell_traits<typename CellType::ElementTag,
                                                topolevel>::ElementTag              ElementTag;
    typedef viennagrid::element<DomainConfig, ElementTag>                           ElementType;

    static void apply(CellType & cell, long * mappingList, long mapindex)
    {
      long index = mapindex;
      //std::cout << "Filling..." << std::endl;

      MappingIterator_impl<FEMConfig, ElementType, IntegrationDomain, map_u>
        mapit(cell, *(cell.template getLevelIteratorBegin<topolevel>()), 0);
      while (mapit.valid())
      {
          mappingList[index] = (*mapit);
          ++index; ++mapit;
      }

      //go to next topolevel:
      MappingFiller<FEMConfig, CellType, IntegrationDomain,
                    map_u, topolevel + 1>::apply(cell, mappingList, index);
    }
  };

  //no bfs on this level, continue with next topolevel:
  template <typename FEMConfig, typename CellType, typename IntegrationDomain, bool map_u,
              long topolevel>
  struct MappingFiller<FEMConfig, CellType, IntegrationDomain, map_u, topolevel, 0, true>
  {
    typedef typename CellType::Configuration             DomainConfig;
    typedef typename viennagrid::subcell_traits<typename CellType::ElementTag,
                                      topolevel>::ElementTag              ElementTag;
    typedef viennagrid::element<DomainConfig, ElementTag>                           ElementType;

    static void apply(CellType & cell, long * mappingList, long mapindex)
    {
      //go to next topolevel:
      MappingFiller<FEMConfig, CellType, IntegrationDomain,
                    map_u, topolevel + 1>::apply(cell, mappingList, mapindex);
    }
  };

  //end of iteration:
  template <typename FEMConfig, typename ElementType, typename IntegrationDomain, bool map_u,
              long topolevel, long mapnum>
  struct MappingFiller<FEMConfig, ElementType, IntegrationDomain, map_u, topolevel, mapnum, false>
  {
    typedef typename ElementType::Configuration             DomainConfig;
    typedef typename ASSEMBLY_CELL_TYPE<DomainConfig, IntegrationDomain>::ResultType
                                                                CellType;

    static void apply(CellType & cell, long * mappingList, long mapindex) {}
  };

  //helper: determine mapping numbers of cell:
  template <typename CellTag, typename BasisfunctionTag, long topolevel = CellTag::TopoLevel>
  struct MappingNumForCell
  {
    enum{ ReturnValue = 1 // (BasisFuncNum< CellTag,
                          //              BasisfunctionTag,
                          //              topolevel>::ReturnValue *
                          //TopologyLevel<CellTag, topolevel>::ElementNum) +
                         //MappingNumForCell<CellTag, BasisfunctionTag, topolevel - 1>::ReturnValue 
        };
  };

  //stop iteration:
  template <typename CellTag, typename BasisfunctionTag>
  struct MappingNumForCell<CellTag, BasisfunctionTag, 0>
  {
    enum{ ReturnValue = 1 //BasisFuncNum< CellTag,
                          //             BasisfunctionTag,
                          //             0
                          //           >::ReturnValue * TopologyLevel<CellTag, 0>::ElementNum 
        };
  };


  //Final mapping iterator:
  template <typename FEMConfig, typename CellType, typename IntegrationDomain,
              bool map_u = false, typename BFTag = typename FEMConfig::BasisfunctionTag>
  class MappingIterator
  {
    typedef typename CellType::Configuration             DomainConfig;
// *    typedef typename ASSEMBLY_CELL_TYPE<DomainConfig, IntegrationDomain>::ResultType

                                                                CellType;* //
    typedef typename CellType::ElementTag                   CellTag;
    typedef typename viennagrid::DomainTypes<DomainConfig>::VertexType  VertexType;

    typedef typename FEMConfig::MappingKey                  MappingKey;
    typedef typename FEMConfig::BoundaryData                BoundaryData;
    typedef typename FEMConfig::BoundaryData2               BoundaryData2;
    typedef typename FEMConfig::BoundaryData3               BoundaryData3;

    enum{ mapnum = MappingNumForCell<CellTag, BFTag>::ReturnValue };


    public:
      MappingIterator(CellType & cell) :
        ppVertices(&(cell.vertices_[0])),
        curIndex(0)
      {
        //std::cout << "Setting up mapping iterator" << std::endl;
        MappingFiller<FEMConfig, CellType, IntegrationDomain, map_u>::apply(cell, mappingList, 0);
        //std::cout << "Finished" << std::endl;
      }

      //NOTE: Shallow copy of MappingIterator **DESIRED**, therefore no copy-constructor needed!

      long operator*() const { return mappingList[curIndex]; }
      MappingIterator & operator++() { ++curIndex; ++ppVertices; return *this; }
      bool valid() const { return curIndex < MappingNumForCell<CellTag, BFTag>::ReturnValue; }
      double getBoundaryValue(int component)
      {
        if (curIndex < viennagrid::subcell_traits<CellTag,0>::ElementNum)
        {
          if (component == 0)
            return ((*ppVertices)->template retrieveQuantity< double >( BoundaryData() ));
          if (component == 1)
            return ((*ppVertices)->template retrieveQuantity< double >( BoundaryData2() ));
          if (component == 2)
            return ((*ppVertices)->template retrieveQuantity< double >( BoundaryData3() ));
        }
        return 0.0;
      }


    private:
      VertexType **ppVertices;
      long mappingList[MappingNumForCell<CellTag, BFTag>::ReturnValue];
      unsigned long curIndex;
  };

  //for "simplest" finite elements, i.e. linear, do not use a mappingList:
  template <typename FEMConfig, typename CellType, typename IntegrationDomain>
  class MappingIterator<FEMConfig, CellType, IntegrationDomain, false, LinearBasisfunctionTag>
  {
      typedef typename CellType::Configuration             DomainConfig;
      typedef typename CellType::ElementTag                   CellTag;
      typedef typename viennagrid::DomainTypes<DomainConfig>::VertexType  VertexType;

      typedef typename FEMConfig::MappingKey                  MappingKey;
      typedef typename FEMConfig::BoundaryData                BoundaryData;
      typedef typename FEMConfig::BoundaryData2               BoundaryData2;
      typedef typename FEMConfig::BoundaryData3               BoundaryData3;

      typedef typename viennagrid::IteratorTypes<CellType, 0>::ResultType    LevelIterator;

    public:
      MappingIterator(CellType & cell) :
        ppVertices(cell.vertices_),
        index(0)
      {
        VertexType **temp = ppVertices;
        for (long j=0; j<viennagrid::subcell_traits<CellTag,0>::ElementNum; ++j, ++temp)
        {
          retval[j] = (*temp)->template retrieveQuantity<long>(MappingKey());
        };
      }

      long operator*() const
      {
        //return curEl.template retrieveQuantity<long>(MappingKey());
        return retval[index];
      };

      MappingIterator & operator++()
      {
        ++index; ++ppVertices;
        return *this;
      }

      bool valid() const { return index < viennagrid::subcell_traits<CellTag,0>::ElementNum; }

      double getBoundaryValue(int component)
      {
        if (component == 0)
          return ((*ppVertices)->template retrieveQuantity< double >( BoundaryData() ));
        if (component == 1)
          return ((*ppVertices)->template retrieveQuantity< double >( BoundaryData2() ));
        if (component == 2)
          return ((*ppVertices)->template retrieveQuantity< double >( BoundaryData3() ));
        return 0.0;
      }

    private:
      VertexType **ppVertices;
      long index;
      long retval[viennagrid::subcell_traits<CellTag,0>::ElementNum];
  };

  //degenerate 1D-boundary-integrals:
  template <typename FEMConfig, typename DomainConfig, typename IntegrationDomain, bool map_u, typename BFTag>
  class MappingIterator<FEMConfig, viennagrid::element<DomainConfig, viennagrid::point_tag>, IntegrationDomain, map_u, BFTag>
  {
      typedef typename viennagrid::DomainTypes<DomainConfig>::VertexType  VertexType;
      typedef viennagrid::element<DomainConfig, viennagrid::point_tag>                   CellType;

      typedef typename FEMConfig::MappingKey                  MappingKey;
      typedef typename FEMConfig::BoundaryData                BoundaryData;
      typedef typename FEMConfig::BoundaryData2               BoundaryData2;
      typedef typename FEMConfig::BoundaryData3               BoundaryData3;

      typedef typename viennagrid::IteratorTypes<CellType, 0>::ResultType    LevelIterator;

    public:
      MappingIterator(CellType & cell) :
        pVertices(&cell),
        isValid(true)
      { }

      long operator*() const
      {
        return pVertices->template retrieveQuantity<long>(MappingKey());
      };

      MappingIterator & operator++()
      {
        isValid = false;
        return *this;
      }

      bool valid() const { return isValid; }

      double getBoundaryValue(int component)
      {
        if (component == 0)
          return (pVertices->template retrieveQuantity< double >( BoundaryData() ));
        if (component == 1)
          return (pVertices->template retrieveQuantity< double >( BoundaryData2() ));
        if (component == 2)
          return (pVertices->template retrieveQuantity< double >( BoundaryData3() ));
        return 0.0;
      }

    private:
      VertexType *pVertices;
      bool isValid;
  };

  //degenerate 1D-boundary-integrals: map u = true
  // *
  template <typename FEMConfig, typename DomainConfig, long id, typename BFTag>
  class MappingIterator<FEMConfig, element<DomainConfig, PointTag>, Interface<id>, true, BFTag>
  {
      typedef typename viennagrid::DomainTypes<DomainConfig>::SegmentType SegmentType;
      typedef typename viennagrid::DomainTypes<DomainConfig>::VertexType  VertexType;
      typedef element<DomainConfig, PointTag>                   CellType;

      typedef typename FEMConfig::MappingKey                  MappingKey;
      typedef typename FEMConfig::BoundaryData                BoundaryData;
      typedef typename FEMConfig::BoundaryData2               BoundaryData2;
      typedef typename FEMConfig::BoundaryData3               BoundaryData3;

      typedef typename IteratorTypes<CellType, 0>::ResultType    LevelIterator;

    public:
      MappingIterator(CellType & cell) :
        isValid(true)
      {
        //toggle active segment:
        CellType * twinCell = cell.template retrieveQuantity<CellType *>(Interface<id>());
        SegmentType * twinSeg = cell.template retrieveQuantity<SegmentType *>(Interface<id>());
        cell.setCurrentSegment(*twinSeg);

        retval = twinCell->template retrieveQuantity<long>(MappingKey());

        //toggle active segment again:
        SegmentType * seg = twinCell->template retrieveQuantity<SegmentType *>(Interface<id>());
        twinCell->setCurrentSegment(*seg);
      }

      long operator*() const
      {
        return retval;
      };

      MappingIterator & operator++()
      {
        isValid = false;
        return *this;
      }

      bool valid() const { return isValid; }

      double getBoundaryValue(int component)
      {
        return 0.0;
      }

    private:
      long retval;
      bool isValid;
  }; * //

  //remove disambiguity:
  template <typename FEMConfig, typename DomainConfig, typename IntegrationDomain>
  class MappingIterator<FEMConfig, viennagrid::element<DomainConfig, viennagrid::point_tag>, IntegrationDomain, false, LinearBasisfunctionTag>
  {
      typedef typename viennagrid::DomainTypes<DomainConfig>::VertexType  VertexType;
      typedef viennagrid::element<DomainConfig, viennagrid::point_tag>                   CellType;

      typedef typename FEMConfig::MappingKey                  MappingKey;
      typedef typename FEMConfig::BoundaryData                BoundaryData;
      typedef typename FEMConfig::BoundaryData2               BoundaryData2;
      typedef typename FEMConfig::BoundaryData3               BoundaryData3;

      typedef typename viennagrid::IteratorTypes<CellType, 0>::ResultType    LevelIterator;

    public:
      MappingIterator(CellType & cell) :
        pVertices(&cell),
        isValid(true)
      { }

      long operator*() const
      {
        return pVertices->template retrieveQuantity<long>(MappingKey());
      };

      MappingIterator & operator++()
      {
        isValid = false;
        return *this;
      }

      bool valid() const { return isValid; }

      double getBoundaryValue(int component)
      {
        if (component == 0)
          return (pVertices->template retrieveQuantity< double >( BoundaryData() ));
        if (component == 1)
          return (pVertices->template retrieveQuantity< double >( BoundaryData2() ));
        if (component == 2)
          return (pVertices->template retrieveQuantity< double >( BoundaryData3() ));
        return 0.0;
      }

    private:
      VertexType *pVertices;
      bool isValid;
  };


  ////////////////// MAPPING SETUP /////////////////////////

  //recursive mapping: tag all elements at level 'topolevel'
  template <typename FEMConfig, long topolevel = 0, bool valid = true>
  struct MappingCreator
  {
    template <typename Segment>
    static long apply(Segment & seg, long startindex)
    {
      typedef typename Segment::Configuration                         DomainConfiguration;
      typedef typename DomainConfiguration::CellTag                   CellTag;
  
      typedef typename viennagrid::subcell_traits<CellTag, topolevel>::ElementTag    ElementTag;
      typedef viennagrid::element<DomainConfiguration, ElementTag>                    ElementType;
      typedef typename viennagrid::DomainTypes<DomainConfiguration>::SegmentType    SegmentType;

      typedef typename viennagrid::IteratorTypes<Segment, topolevel>::ResultType    ElementIterator;
  
      typedef typename FEMConfig::BasisfunctionTag                BasisFunctionTag;
      typedef typename FEMConfig::SegmentConnection               SegmentConnectionType;
      typedef typename FEMConfig::MappingTag                      MappingTag;
      typedef typename FEMConfig::MappingKey                      MappingKey;
      typedef typename FEMConfig::BoundaryKey                     BoundaryKey;
      typedef typename FEMConfig::BoundaryData                    BoundaryData;
      typedef typename FEMConfig::BoundaryData2                   BoundaryData2;
  
      long numIndex = startindex;
      MappingKey mkt;

      if (  1 ) //BasisFuncNum<ElementTag, BasisFunctionTag, topolevel>::ReturnValue == 1
      {
        //allocate needed memory:
        seg.template getLevelIteratorBegin<topolevel>()
          ->template reserveQuantity<long>(mkt);

        std::cout << "Tagging Level " << topolevel << std::endl;
        for (ElementIterator eit = seg.template getLevelIteratorBegin<topolevel>();
              eit != seg.template getLevelIteratorEnd<topolevel>();
              ++eit)
        {
          //std::cout << "Checking vertex: "; eit->print();
          if (MappingTag::apply(eit, BoundaryKey()))
          {
            (*eit).template storeQuantity< long >( mkt, -1 );
            //std::cout << "Boundary Key!!" << std::endl;
          }
          else
          {
            //is this segment coupled?
            if (eit->template hasQuantity<ElementType *>(SegmentConnectionType()))
            {
              ElementType * twinCell = 
                    eit->template retrieveQuantity<ElementType *>(SegmentConnectionType());
              //std::cout << twinCell << std::endl;
              SegmentType * twinSeg =
                    eit->template retrieveQuantity<SegmentType *>(SegmentConnectionType());

              //std::cout << twinSeg << std::endl;

              //element with lower address is the master:
              if (&(*eit) < twinCell)
              {
                //eit->print();
                eit->setCurrentSegment(*twinSeg);
                twinCell->template reserveQuantity<long>(mkt);
                twinCell->template storeQuantity< long >( mkt, numIndex );
                eit->setCurrentSegment(seg);

                (*eit).template storeQuantity< long >( mkt, numIndex );
                ++numIndex;
              }
            }
            else  //standard mapping index write:
            {
              (*eit).template storeQuantity< long >( mkt, numIndex );
              //eit->print();
              ++numIndex;
            }
          }
        };
      }
      else if ( 2 ) //BasisFuncNum<ElementTag, BasisFunctionTag, topolevel>::ReturnValue > 1)
      {
        seg.template getLevelIteratorBegin<topolevel>()
          ->template reserveQuantity< std::vector<long> >(mkt);
  
        std::cout << "Tagging Level " << topolevel << std::endl;
        for (ElementIterator eit = seg.template getLevelIteratorBegin<topolevel>();
              eit != seg.template getLevelIteratorEnd<topolevel>();
              ++eit)
        {
            std::vector<long> mapVec;
            mapVec.resize( 1 ); //BasisFuncNum<ElementTag, BasisFunctionTag, topolevel>::ReturnValue);
    
            if (MappingTag::apply(eit, BoundaryKey()))
            {
              for (long i=0;
                    i<1;//BasisFuncNum<ElementTag, BasisFunctionTag, topolevel>::ReturnValue;
                    ++i)
              {
                mapVec[i] = -1;
              }
            }
            else
            {
              for (long i=0;
                    i<1;//BasisFuncNum<ElementTag, BasisFunctionTag, topolevel>::ReturnValue;
                    ++i)
              {
                mapVec[i] = numIndex;
                ++numIndex;
              }
            }
    
            //is this segment coupled?
            if (eit->template hasQuantity<ElementType *>(SegmentConnectionType()))
            {
              std::cout << "WARNING: CONNECTED SEGMENTS NOT IMPLEMENTED FOR THIS BASIS FUNCTIONS!" << std::endl;
            }

            (*eit).template storeQuantity< std::vector<long> >( mkt, mapVec );
        };
      }

      //continue with next topology-level:
      return MappingCreator<FEMConfig,
                              topolevel + 1,
                              (CellTag::TopoLevel > topolevel)
                            >::apply(seg, numIndex);
    } //apply()
  }; //MappingCreator

  //end recursion when topolevel exceeds cell-level:
  template <typename FEMConfig, long topolevel>
  struct MappingCreator<FEMConfig, topolevel, false>
  {
    template <typename Segment>
    static long apply(Segment & seg, long startindex)
    {
      return startindex;
    }
  };

  //the interface:
  template <typename FEMConfig, typename Segment>
  long createMapping(Segment & seg, long startindex = 0)
  {
    seg.template getLevelIteratorBegin<0>()->setCurrentSegment(seg);

    return MappingCreator<FEMConfig>::apply(seg, startindex);
  } */

}

#endif
