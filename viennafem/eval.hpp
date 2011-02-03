/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_EVAL_HPP
#define VIENNAFEM_EVAL_HPP

#include "viennafem/forwards.h"
#include "viennafem/cell_quan.hpp"

#include "viennamath/equation.hpp"
#include "viennamath/op_tags.hpp"
#include "viennamath/substitute.hpp"
#include "viennamath/diff.hpp"

namespace viennafem
{
  
  //evaluates the weak form on a triangle using 1-point-rule
  //TODO generalize!!
  numeric_type eval_element_matrix_entry(viennamath::expr const & weak_form_lhs)
  {
    std::vector<numeric_type> p(2);
    p[0] = 1.0/3.0;
    p[1] = 1.0/3.0;
    
    return 0.5 * weak_form_lhs.get()->lhs()->eval(p); //TODO: this is pretty ugly...
  }
  
  //evaluates the weak form on a triangle using 1-point-rule
  //TODO generalize!!
  numeric_type eval_element_vector_entry(viennamath::expr const & weak_form_rhs)
  {
    std::vector<numeric_type> p(2);
    p[0] = 1.0/3.0;
    p[1] = 1.0/3.0;
    
    return 0.5 * weak_form_rhs.get()->lhs()->eval(p); //TODO: this is pretty ugly...
  }

}
#endif
