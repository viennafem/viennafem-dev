/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_BASES_ALL_HPP
#define VIENNAFEM_BASES_ALL_HPP

#include "viennafem/forwards.h"

#include "viennafem/bases/line.hpp"

#include "viennafem/bases/quadrilateral.hpp"
#include "viennafem/bases/triangle.hpp"

#include "viennafem/bases/hexahedron.hpp"
#include "viennafem/bases/tetrahedron.hpp"

namespace viennafem
{

  template <typename InterfaceType,
            typename BasisTag = viennafem::none>
  struct basis_factory 
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    
    template <typename T>
    static std::vector<expression_type> get(std::size_t space_id, T const & reference_cell)
    {
      switch (space_id)
      {
        case space_to_id< lagrange_tag<1> >::value:   return basis_factory<InterfaceType, lagrange_tag<1> >::get(reference_cell);
      }
      
      std::cerr << "Space ID unknown!" << std::endl;
      throw "unknown";
      
      return std::vector<expression_type>();
    }
    
  };
  
  
  template <typename InterfaceType>
  struct basis_factory <InterfaceType, lagrange_tag<1> >
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef lagrange_tag<1>                      basis_tag;
    
    static std::vector<expression_type> get(viennafem::unit_interval)
    {
      std::vector<expression_type> ret;
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_interval, 0, 0>::get()[0] );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_interval, 0, 1>::get()[0] );
      
      return ret;
    }
    
    static std::vector<expression_type> get(viennafem::unit_triangle)
    {
      std::vector<expression_type> ret;
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_triangle, 0, 0>::get()[0] );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_triangle, 0, 1>::get()[0] );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_triangle, 0, 2>::get()[0] );
      
      return ret;
    }

    static std::vector<expression_type> get(viennafem::unit_tetrahedron)
    {
      std::vector<expression_type> ret;
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 0, 0>::get()[0] );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 0, 1>::get()[0] );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 0, 2>::get()[0] );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 0, 3>::get()[0] );
      
      return ret;
    }

  };
  
  
  
}


#endif
