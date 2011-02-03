/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_FORWARDS_H
#define VIENNAFEM_FORWARDS_H

#include "viennadata/interface.hpp"
#include "viennagrid/forwards.h"
#include "viennamath/forwards.h"

/** @file forwards.h
    @brief This file provides the forward declarations for the main types used within ViennaFEM
*/

namespace viennafem
{
  
  typedef double             numeric_type;

  //a tag for storing mapping indices on the grid
  template <long id>
  struct mapping_key_type {}; 

  
  //integration tags on cells:
  struct AnalyticIntegrationTag {};

  struct LinearIntegrationTag {};
  struct QuadraticIntegrationTag {};
  struct CubicIntegrationTag {};

  struct QuinticIntegrationTag {};

  struct OrderSevenIntegrationTag {};

  struct CellIntegrationTag {};


  
  //Basisfunction-Treatment:
  struct TypeListTag {};      //use typelists
  struct TypeErasureTag {};   //use type erasure

  //Tag for handling of functional determinant
  struct DtDxStoreAll{ };   //all elements of functional determinant are stored. Usually fastest, but high memory consumption.
  struct DtDxStoreDetOnly {}; //only determinant is stored, other elements computed on-the-fly. Needs less memory, but usually increases runtime a little bit
  struct DtDxOnAccess {};   //All values are (repeatedly) computed on-the-fly. Use only in situations where every bit of memory is important.
  struct DtDxStoreStatically {}; //computes the functional determinant at the point of first access to a cell in static variables. Must not be used together with multi-threading!
  
  
  //Continuity-Tag:
  struct C0_Tag {};   //ordinary H1-basisfunctions
  struct C1_Tag {};   //Hermite

  //Basis function tags
  struct LinearBasisfunctionTag
  {
    enum{ degree = 1 };
    //typedef C0_Tag                ContinuityTag;
  };

  struct QuadraticBasisfunctionTag
  {
    enum{ degree = 2 };
    //typedef C0_Tag                ContinuityTag;
  };

  struct CubicBasisfunctionTag
  {
    enum{ degree = 3 };
  };

  struct QuarticBasisfunctionTag
  {
    enum{ degree = 4 };
  };

  struct QuinticBasisfunctionTag
  {
    enum{ degree = 5 };
  };

  // A tag for scalar-valued PDEs
  struct ScalarTag {
    enum { dim = 1 };
  };

  // A tag for vector-valued PDEs
  template <typename DimensionTag>
  struct VectorTag {
    enum { dim = DimensionTag::dim };
  };

  
  //Integration domain:
  struct Omega {};
  
  template <long id>
  struct Gamma {};
  
  
  
  template <unsigned long local_index,
            unsigned long global_index>
  struct dt_dx_key {};
  
  template <unsigned long local_index,
            unsigned long global_index>
  std::ostream & operator<<(std::ostream & stream,
                            dt_dx_key<local_index, global_index> const & dummy)
  {
    stream << "dt_dx_key<" << local_index << "," << global_index << ">";
    return stream;
  }


  struct det_dF_dt_key {};
  
  std::ostream & operator<<(std::ostream & stream, det_dF_dt_key const & dummy)
  {
    stream << "det_dF_dt_key";
    return stream;
  }
  
  
  template <typename CellTag>
  struct dt_dx_handler;
  
  //Mapping Tags:
  
  /** @brief A tag indicating that degrees of freedoms are assigned regardless of whether a (Dirichlet) boundary condition is specified there. */
  struct FullMappingTag             //Create the 'full' mapping, including boundary vertices. This means that 
  {
    template <typename Iterator, typename BoundaryKey>
    static bool apply(Iterator & it, BoundaryKey const & bk) { return false; }
  };

  /** @brief A tag that specifies that no degrees of freedom are assigned at (Dirichlet) boundaries */
  struct NoBoundaryMappingTag       //exclude boundary vertices from mapping.
  {
    template <typename Iterator, typename BoundaryKey>
    static bool apply(Iterator & it, BoundaryKey const & bk)
    {
      return viennadata::access<BoundaryKey, bool>(bk)(*it);
    }
  };
  
  // define a key and configure viennadata to use a type-based dispatch:
  struct boundary_key {};
  struct mapping_key {};

  

}

//configure vienndata for type-based dispatch on boundary_key, mapping_key, dt_dx and det_F:
namespace viennadata
{
  template <>
  struct dispatch_traits<viennafem::boundary_key>
  {
    typedef type_key_dispatch_tag    tag;
  };
  
  template <>
  struct dispatch_traits<viennafem::mapping_key>
  {
    typedef type_key_dispatch_tag    tag;
  };
  
  template <unsigned long local_index,
            unsigned long global_index>
  struct dispatch_traits<viennafem::dt_dx_key<local_index,
                                              global_index>
                        >
  {
    typedef type_key_dispatch_tag    tag;
  };
  
  template <>
  struct dispatch_traits<viennafem::det_dF_dt_key>
  {
    typedef type_key_dispatch_tag    tag;
  };
  
}


#endif