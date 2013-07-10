#ifndef VIENNAFEM_FORWARDS_H
#define VIENNAFEM_FORWARDS_H

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

/**
 @mainpage Source Code Documentation for ViennaFEM 1.0.0

 Thank you for downloading ViennaFEM!
 These pages provide the Doxygen-generated source code documentation with the many details on the functions implemented in ViennaFEM.
 
 For a general overview over the types and functionality provided by ViennaFEM, please refer to the file doc/viennafem.pdf. 
 Since ViennaFEM depends on the libraries ViennaCL, ViennaData, ViennaGrid and ViennaMath, you might also want to check out the PDF manuals
 doc/viennacl.pdf, doc/viennadata.pdf, doc/viennagrid.pdf and doc/viennamath.pdf.
*/


#include "viennadata/api.hpp"
#include "viennagrid/forwards.hpp"
#include "viennamath/forwards.h"

/** @file forwards.h
    @brief This file provides the forward declarations for the main types used within ViennaFEM
*/

/** @brief The main namespace for ViennaFEM. All user functionality resides in this namespace. */
namespace viennafem
{
  
  /** @brief The default floating point type to be used in ViennaFEM.
   * 
   *  Feel free to change this typedef to a high-precision type if required. 
   *  Keep in mind that only float and double types can be used for GPU acceleration.
   */
  typedef double             numeric_type;

  
  /** @brief A tag for storing mapping indices on the grid  
   * 
   *  @tparam id    An ID for the simulation to allow for several simulations on the same grid.
   */
  template <long id>
  struct mapping_key_type {}; 

  
  struct none;
  
  
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
  
  /** @brief A tag for the unit interval [0,1].
   *
   * Some families of basis functions can be defined in a more compact manner here than on [-1,1]. */
  struct unit_interval {};

  /** @brief A tag representing the interval [-1,1]. Particularly useful for higher-order FEM. */
  struct symmetric_interval {};


  
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

  
  /** @brief A tag for the square with vertices at (0,0), (1,0), (0,1), (1,1)
   * 
   * Often used in the engineering community for first and second order FEM.
   */
  struct unit_square {};
  
  /** @brief A tag for the square with vertices at (-1,-1), (1, -1), (-1, 1), (1, 1)
   * 
   * Allows for a construction of high-order FEM in a very natural way out of [-1,1].
   */
  struct symmetric_square {};
  
  
    
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
  
  
  /** @brief A tag for the cube with vertices at (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1), (1,0,1), (0,1,1), (1,1,1)
   * 
   * Commonly used for low-order form functions
   */
  struct unit_cube {};
  
  /** @brief A tag for the cube with vertices at (-1,-1,-1), (1,-1,-1), (-1,1,-1), (1,1,-1), (-1,-1,1), (1,-1,1), (-1,1,1), (1,1,1)
   * 
   * Used for high-order FEM, since basis functions can be constructed more easily out of [-1,1].
   */
  struct symmetric_cube {};
  
  
  //
  // Basis functions
  //
  
  /** @brief A tag identifying the family of Lagrange basis functions
   * 
   * @tparam order     Polynomial degree of the basis along edges (degree might be higher or lower in the cell interior)
   */
  template <std::size_t order>
  struct lagrange_tag {};
  
  /** @brief Provides a unique ID from a basis function tag (compiletime-runtime translation)  */
  template <typename T>
  struct space_to_id {};
  
  /** @brief Specialization of the unique ID facility for the Lagrange family */
  template <std::size_t order>
  struct space_to_id< lagrange_tag<order> >
  {
    enum { value = order };
  };
  
  
  /** @brief Metafunction for returning the reference cell for a cell type and a basis function type */
  template <typename Cell, typename T>
  struct reference_cell_for_basis {};
  
  /** @brief Defining the unit interval to be used for the Lagrange basis on lines */
  template <std::size_t order>
  struct reference_cell_for_basis < viennagrid::line_tag, lagrange_tag<order> >
  {
    typedef unit_interval   type;
  };

  /** @brief Defining the unit interval to be used for the Lagrange basis on lines */
  template <std::size_t order>
  struct reference_cell_for_basis < viennagrid::triangle_tag, lagrange_tag<order> >
  {
    typedef unit_triangle   type;
  };

  /** @brief Defining the unit interval to be used for the Lagrange basis on lines */
  template <std::size_t order>
  struct reference_cell_for_basis < viennagrid::quadrilateral_tag, lagrange_tag<order> >
  {
    typedef unit_square   type;
  };

  
  /** @brief Defining the unit tetrahedron to be used for the Lagrange basis on lines */
  template <std::size_t order>
  struct reference_cell_for_basis < viennagrid::tetrahedron_tag, lagrange_tag<order> >
  {
    typedef unit_tetrahedron   type;
  };
  
  /** @brief Defining the unit cube to be used for the Lagrange basis on lines */
  template <std::size_t order>
  struct reference_cell_for_basis < viennagrid::hexahedron_tag, lagrange_tag<order> >
  {
    typedef unit_cube   type;
  };
  

  
  // Basis function retrieval:
  /** @brief Facility for retrieving the local basis of a cell. Used internally only.
   *
   *  @tparam InterfaceType     Runtime interface type of ViennaMath expressions
   *  @tparam BasisTag          A tag identifying the basis
   *  @tparam ReferenceCell     The reference element used for the basis
   *  @tparam TopologyDim       Topological dimension of the the element for which the basis function should be retrieved
   *  @tparam ElementID         The ID of the topological element on the reference cell
   */
  template <typename InterfaceType,
            typename BasisTag,
            typename ReferenceCell,
            std::size_t TopologyDim,  //topology dimension of the element for which the basis is to be returned
            std::size_t ElementID>    //element ID (0, ... num-1) of the element at dimension 'topology_dim'
  struct local_basis
  {
    typedef typename BasisTag::ERROR_BASIS_NOT_DEFINED_FOR_THE_PROVIDED_PARAMETER_SET   error_type;
  };                    
  
  
  
  /** @brief A tag class used for storing and accessing the partial derivatives of the element mappings */
  template <unsigned long local_index,
            unsigned long global_index>
  struct dt_dx_key 
  {
    // When using a trivial key operator< MUST return false
    bool operator<(dt_dx_key<0,0> const &)
    {
        return false;
    }
  };
  
  /** @brief Convenience overload for converting a det_dF_dt_key to a string and streaming it to an output-stream */
  template <unsigned long local_index,
            unsigned long global_index>
  std::ostream & operator<<(std::ostream & stream,
                            dt_dx_key<local_index, global_index> const & dummy)
  {
    stream << "dt_dx_key<" << local_index << "," << global_index << ">";
    return stream;
  }


  /** @brief A tag class used for storing and accessing the Jacobian determinant via ViennaData on ViennaGrid objects */
  struct det_dF_dt_key {
    // When using a trivial key operator< MUST return false
    bool operator<(det_dF_dt_key const &)
    {
        return false;
    }  
  };
  
  /** @brief Convenience overload for converting a det_dF_dt_key to a string and streaming it to an output-stream */
  inline std::ostream & operator<<(std::ostream & stream, det_dF_dt_key const & /*dummy*/)
  {
    stream << "det_dF_dt_key";
    return stream;
  }
  
  /** @brief The main facility for computing transformations from the reference element to the physical element   */
  template <typename DomainType, typename StorageType, typename CellTag>
  struct dt_dx_handler;
  

  /** @brief A tag class used for storing and accessing boundary data via ViennaData from ViennaGrid objects */
  class boundary_key 
  {
    public:
      boundary_key(long id) : id_(id) {}
      
      bool operator<(boundary_key const & other) const { return id_ < other.id_; }
    
    private:
      long id_;
  };

  
  /** @brief A tag class used for storing and accessing unknown indices (global basisfunction numbers) via ViennaData from ViennaGrid objects */
  class mapping_key
  {
    public:
      mapping_key(long id) : id_(id) {}
      
      bool operator<(mapping_key const & other) const { return id_ < other.id_; }
    
    private:
      long id_;
  };
  
  /** @brief Contains implementation details of ViennaFEM not intended to be used by the library user directly. */
  namespace detail {}
  
  /** @brief Provides input and output routines for Finite Element related quantities. Refer to viennagrid::io for details on entirely mesh-related input-output functionality */
  namespace io {}

} // namespace viennafem







#endif
