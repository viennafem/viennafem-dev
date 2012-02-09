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

#include "viennadata/api.hpp"
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

  
  //
  // Integration
  //
  
  //prototype for Gauss quadrature rules
  template <typename ElementType, unsigned long order, typename InterfaceType = viennamath::default_interface_type>
  struct rt_gauss_quad_element;
  
  //prototype for Strang quadrature rules on triangles
  template <typename ElementType, unsigned long order, typename InterfaceType = viennamath::default_interface_type>
  struct rt_strang_quad_element;

  
  //prototype for Strang quadrature rules on tetrahedra
  template <typename ElementType, unsigned long order, typename InterfaceType = viennamath::default_interface_type>
  struct rt_keast_quad_element;

  
  //
  //
  //
  
  
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
  
  template <typename T>
  struct space_to_id;
  
  template <>
  struct space_to_id<LinearBasisfunctionTag>
  {
    enum { value = 42 };
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
  class boundary_key 
  {
    public:
      boundary_key(long id) : id_(id) {}
      
      bool operator<(boundary_key const & other) const { return id_ < other.id_; }
    
    private:
      long id_;
  };
  
  class mapping_key
  {
    public:
      mapping_key(long id) : id_(id) {}
      
      bool operator<(mapping_key const & other) const { return id_ < other.id_; }
    
    private:
      long id_;
  };

  

}

//configure vienndata for type-based dispatch on boundary_key, mapping_key, dt_dx and det_F:
namespace viennadata
{
  namespace config
  {
    template <unsigned long local_index,
              unsigned long global_index>
    struct key_dispatch<viennafem::dt_dx_key<local_index,
                                             global_index>
                       >
    {
      typedef type_key_dispatch_tag    tag;
    };
    
    template <>
    struct key_dispatch<viennafem::det_dF_dt_key>
    {
      typedef type_key_dispatch_tag    tag;
    };
    
    //
    // tell ViennaData to use the get_id() member for vertices as identification mechanism
    //
    template <typename ConfigType>
    struct object_identifier<viennagrid::element_t<ConfigType, viennagrid::point_tag> >
    {
      typedef object_provided_id    tag;
      typedef size_t                id_type;

      static size_t get(viennagrid::element_t<ConfigType, viennagrid::point_tag> const & obj) { return obj.id(); }
    };

    template <typename ConfigType>
    struct object_identifier<viennagrid::element_t<ConfigType, viennagrid::tetrahedron_tag> >
    {
      typedef object_provided_id    tag;
      typedef size_t                id_type;

      static size_t get(viennagrid::element_t<ConfigType, viennagrid::tetrahedron_tag> const & obj) { return obj.id(); }
    };
    
    //
    // store data densely, no matter which key type is used:
    //
    template <typename KeyType, typename ValueType, typename ConfigType>
    struct storage<KeyType, ValueType, viennagrid::element_t<ConfigType, viennagrid::point_tag> >
    {
      typedef dense_data_tag    tag;
    };

    template <typename KeyType, typename ValueType, typename ConfigType>
    struct storage<KeyType, ValueType, viennagrid::element_t<ConfigType, viennagrid::tetrahedron_tag> >
    {
      typedef dense_data_tag    tag;
    };
    
  }
}


#endif