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
  // Reference elements
  //
  
  /** @brief A tag representing the interval [-1,1]. Particularly useful for higher-order FEM. */
  struct symmetric_interval {};

  /** @brief A tag for the unit interval [0,1].
   *
   * Some families of basis functions can be defined in a more compact manner here than on [-1,1]. */
  struct unit_interval {};

  
  /** @brief A tag for the triangle with vertices (0,0), (1,0), (0,1).
   *
   * Home of form functions.
   */
  struct unit_triangle {};
  
  /** @brief A tag for the triangle with vertices at (-1,-1), (1,-1), (-1,1).
   *
   * Commonly used for the construction of high-order FEM bases.
   */
  struct symmetric_triangle {};

  
  /** @brief A tag for the quadrilateral with vertices at (0,0), (1,0), (0,1), (1,1)
   * 
   * Often used in the engineering community for first and second order FEM.
   */
  struct unit_quadrilateral {};
  
  /** @brief A tag for the quadrilateral with vertices at (-1,-1), (1, -1), (-1, 1), (1, 1)
   * 
   * Allows for a construction of high-order FEM in a very natural way out of [-1,1].
   */
  struct symmetric_quadrilateral {};
  
  
  /** @brief A tag for the tetrahedron with vertices at (0,0,0), (1,0,0), (0,1,0), (0,0,1)
   * 
   * Commonly used for simple form functions of first or second order. 
   */
  struct unit_tetrahedron {};
  
  /** @brief A tag for the tetrahedron with vertices at (-1,-1,-1), (1,-1,-1), (-1,1,-1), (-1,-1,1)
   * 
   * Handy for high-order FEM, since hierarchical bases can be constructed out of the interval [-1,1]
   */
  struct symmetric_tetrahedron {};
  
  
  /** @brief A tag for the hexahedron with vertices at (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1), (1,0,1), (0,1,1), (1,1,1)
   * 
   * Commonly used for low-order form functions
   */
  struct unit_hexahedron {};
  
  /** @brief A tag for the hexahedron with vertices at (-1,-1,-1), (1,-1,-1), (-1,1,-1), (1,1,-1), (-1,-1,1), (1,-1,1), (-1,1,1), (1,1,1)
   * 
   * Used for high-order FEM, since basis functions can be constructed more easily out of [-1,1].
   */
  struct symmetric_hexahedron {};
  
  
  //
  // Basis functions
  //
  
  struct LinearBasisfunctionTag {};
  
  template <typename T>
  struct space_to_id {};
  
  template <>
  struct space_to_id<LinearBasisfunctionTag>
  {
    enum { value = 1 };
  };
  
  template <typename Cell, typename T>
  struct reference_cell_for_basis {};
  
  template <>
  struct reference_cell_for_basis < viennagrid::line_tag, LinearBasisfunctionTag>
  {
    typedef unit_interval   type;
  };

  template <>
  struct reference_cell_for_basis < viennagrid::triangle_tag, LinearBasisfunctionTag>
  {
    typedef unit_triangle   type;
  };

  template <>
  struct reference_cell_for_basis < viennagrid::tetrahedron_tag, LinearBasisfunctionTag>
  {
    typedef unit_tetrahedron   type;
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