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
        case space_to_id< lagrange_tag<2> >::value:   return basis_factory<InterfaceType, lagrange_tag<2> >::get(reference_cell);
        //case space_to_id< lagrange_tag<3> >::value:   return basis_factory<InterfaceType, lagrange_tag<3> >::get(reference_cell);
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
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_interval, 0, 0>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_interval, 0, 1>::get() );
      
      return ret;
    }
    
    static std::vector<expression_type> get(viennafem::unit_quadrilateral)
    {
      std::vector<expression_type> ret;
      ret.reserve(4);
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_quadrilateral, 0, 0>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_quadrilateral, 0, 1>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_quadrilateral, 0, 2>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_quadrilateral, 0, 3>::get() );
      
      return ret;
    }

    static std::vector<expression_type> get(viennafem::unit_triangle)
    {
      std::vector<expression_type> ret;
      ret.reserve(3);
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_triangle, 0, 0>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_triangle, 0, 1>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_triangle, 0, 2>::get() );
      
      return ret;
    }

    static std::vector<expression_type> get(viennafem::unit_hexahedron)
    {
      std::vector<expression_type> ret;
      ret.reserve(8);
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 0, 0>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 0, 1>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 0, 2>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 0, 3>::get() );

      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 0, 4>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 0, 5>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 0, 6>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 0, 7>::get() );
      
      return ret;
    }

    static std::vector<expression_type> get(viennafem::unit_tetrahedron)
    {
      std::vector<expression_type> ret;
      ret.reserve(4);
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 0, 0>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 0, 1>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 0, 2>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 0, 3>::get() );
      
      return ret;
    }

  };
  

  //
  // quadratic
  //
  template <typename InterfaceType>
  struct basis_factory <InterfaceType, lagrange_tag<2> >
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef lagrange_tag<2>                      basis_tag;
    
    static std::vector<expression_type> get(viennafem::unit_interval)
    {
      std::vector<expression_type> ret;
      ret.reserve(3);
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_interval, 0, 0>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_interval, 0, 1>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_interval, 1, 0>::get() );
      
      return ret;
    }
    
    static std::vector<expression_type> get(viennafem::unit_quadrilateral)
    {
      std::vector<expression_type> ret;
      ret.reserve(8);
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_quadrilateral, 0, 0>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_quadrilateral, 0, 1>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_quadrilateral, 0, 2>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_quadrilateral, 0, 3>::get() );

      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_quadrilateral, 1, 0>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_quadrilateral, 1, 1>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_quadrilateral, 1, 2>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_quadrilateral, 1, 3>::get() );
      
      return ret;
    }

    static std::vector<expression_type> get(viennafem::unit_triangle)
    {
      std::vector<expression_type> ret;
      ret.reserve(6);
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_triangle, 0, 0>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_triangle, 0, 1>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_triangle, 0, 2>::get() );

      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_triangle, 1, 0>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_triangle, 1, 1>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_triangle, 1, 2>::get() );
      
      return ret;
    }

    static std::vector<expression_type> get(viennafem::unit_hexahedron)
    {
      std::vector<expression_type> ret;
      ret.reserve(20);
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 0, 0>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 0, 1>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 0, 2>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 0, 3>::get() );

      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 0, 4>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 0, 5>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 0, 6>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 0, 7>::get() );

      //edges:
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 1, 0>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 1, 1>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 1, 2>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 1, 3>::get() );

      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 1, 4>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 1, 5>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 1, 6>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 1, 7>::get() );

      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 1,  8>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 1,  9>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 1, 10>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_hexahedron, 1, 11>::get() );
      
      return ret;
    }

    static std::vector<expression_type> get(viennafem::unit_tetrahedron)
    {
      std::vector<expression_type> ret;
      ret.reserve(10);
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 0, 0>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 0, 1>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 0, 2>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 0, 3>::get() );

      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 1, 0>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 1, 1>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 1, 2>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 1, 3>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 1, 4>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_tetrahedron, 1, 5>::get() );
      
      return ret;
    }

  };

  
  //
  // cubic
  //
  template <typename InterfaceType>
  struct basis_factory <InterfaceType, lagrange_tag<3> >
  {
    typedef viennamath::rt_expr<InterfaceType>   expression_type;
    typedef lagrange_tag<3>                      basis_tag;
    
    static std::vector<expression_type> get(viennafem::unit_interval)
    {
      std::vector<expression_type> ret;
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_interval, 0, 0>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_interval, 0, 1>::get() );
      ret.push_back( local_basis<InterfaceType, basis_tag, viennafem::unit_interval, 1, 0>::get() );
      
      return ret;
    }
    
    // generic compile time error handler for not yet implemented functionality:
    template <typename T>
    static std::vector<expression_type> get(T)
    {
       typedef typename T::ERROR_BASIS_FUNCTION_NOT_IMPLEMENTED       error_type;
       return std::vector<expression_type>(1);
    }
  };  
}


#endif
