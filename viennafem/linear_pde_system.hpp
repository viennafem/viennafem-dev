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
  template <typename InterfaceType>
  class linear_pde_system
  {
    public:
      typedef viennamath::equation<InterfaceType>    equation_type;
      typedef viennamath::function_symbol<InterfaceType>    unknown_type;
      typedef viennafem::linear_pde_options          option_type;
      
      
      void add_pde(viennamath::equation<InterfaceType> const & pde,
                   viennamath::function_symbol<InterfaceType> const & unknown,
                   viennafem::linear_pde_options const & option)
      {
        pdes_.push_back(pde); 
        unknowns_.push_back(unknown);
        options_.push_back(option);
      }
      
      viennamath::equation<InterfaceType> pde(size_t index) const { return pdes_[index]; }
      viennamath::variable<InterfaceType> unknown(size_t index) const { return unknowns_[index]; }
      viennafem::linear_pde_options option(size_t index) const { return options_[index]; }
      
      size_t size() const { return pdes_.size(); }
      
    private:
      std::vector<viennamath::equation<InterfaceType> >  pdes_;
      std::vector<viennamath::function_symbol<InterfaceType> >  unknowns_;
      std::vector<viennafem::linear_pde_options>         options_;
  };
  
  
  template <typename InterfaceType>
  linear_pde_system<InterfaceType> make_linear_pde_system(viennamath::equation<InterfaceType> equ_1,
                                                          viennamath::function_symbol<InterfaceType> unknown_1,
                                                          viennafem::linear_pde_options options_1)
  {
    linear_pde_system<InterfaceType> ret;
    ret.add_pde(equ_1, unknown_1, options_1);
    return ret;
  }
}
#endif
