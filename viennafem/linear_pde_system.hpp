/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_LINEAR_PDE_SYSTEM_HPP
#define VIENNAFEM_LINEAR_PDE_SYSTEM_HPP

#include <vector>

//ViennaFEM includes:
#include "viennafem/forwards.h"
#include "viennafem/cell_quan.hpp"
#include "viennafem/transform.hpp"
#include "viennafem/bases/tetrahedron.hpp"
#include "viennafem/bases/triangle.hpp"
#include "viennafem/eval.hpp"
#include "viennafem/transform/dtdx_triangle.hpp"
#include "viennafem/transform/dtdx_tetrahedron.hpp"
#include "viennafem/weak_form.hpp"
#include "viennafem/linear_pde_options.hpp"

//ViennaMath includes:
#include "viennamath/manipulation/apply_coordinate_system.hpp"

//ViennaData includes:
#include "viennadata/api.hpp"

//ViennaGrid includes:
#include "viennagrid/domain.hpp"

namespace viennafem
{
  template <typename InterfaceType = viennamath::rt_expression_interface<viennamath::default_numeric_type>, 
            typename MappingKeyType  = viennafem::mapping_key, 
            typename BoundaryKeyType = viennafem::boundary_key >
  class linear_pde_system
  {
    public:
      typedef InterfaceType                                 interface_type;
      typedef MappingKeyType                                mapping_key_type;
      typedef BoundaryKeyType                               boundary_key_type;    
      typedef viennamath::rt_equation<InterfaceType>           equation_type;
      typedef viennamath::rt_function_symbol<InterfaceType>    unknown_type;
      typedef viennafem::linear_pde_options                 option_type;
      
      
      void add_pde(viennamath::rt_equation<InterfaceType> const & pde,
                   std::vector< viennamath::rt_function_symbol<InterfaceType> > const & unknown,
                   viennafem::linear_pde_options const & option)
      {
        pdes_.push_back(pde); 
        unknowns_.push_back(unknown);
        options_.push_back(option);
      }
      
      viennamath::rt_equation<InterfaceType> pde(size_t index) const { return pdes_[index]; }
      std::vector<viennamath::rt_function_symbol<InterfaceType> > const & unknown(size_t index) const { return unknowns_[index]; }
      viennafem::linear_pde_options option(size_t index) const { return options_[index]; }
      
      size_t size() const { return pdes_.size(); }
      
    private:
      std::vector<viennamath::rt_equation<InterfaceType> >                       pdes_;
      std::vector<std::vector<viennamath::rt_function_symbol<InterfaceType> > >  unknowns_;
      std::vector<viennafem::linear_pde_options>                              options_;
  };
  
  
  template <typename InterfaceType>
  linear_pde_system<InterfaceType> make_linear_pde_system(viennamath::rt_equation<InterfaceType> equ_1,
                                                          viennamath::rt_function_symbol<InterfaceType> unknown_1,
                                                          viennafem::linear_pde_options options_1)
  {
    linear_pde_system<InterfaceType> ret;
    std::vector<viennamath::rt_function_symbol<InterfaceType> > unknown_vec_1(1);
    unknown_vec_1[0] = unknown_1;
    ret.add_pde(equ_1, unknown_vec_1, options_1);
    return ret;
  }
  
  template <typename InterfaceType>
  linear_pde_system<InterfaceType> make_linear_pde_system(viennamath::rt_equation<InterfaceType> equ_1,
                                                          std::vector<viennamath::rt_function_symbol<InterfaceType> > unknowns_1,
                                                          viennafem::linear_pde_options options_1)
  {
    linear_pde_system<InterfaceType> ret;
    ret.add_pde(equ_1, unknowns_1, options_1);
    return ret;
  }
  
}
#endif
