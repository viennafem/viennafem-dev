/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_LINEAR_PDE_OPTIONS_HPP
#define VIENNAFEM_LINEAR_PDE_OPTIONS_HPP

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

//ViennaMath includes:
#include "viennamath/manipulation/apply_coordinate_system.hpp"

//ViennaData includes:
#include "viennadata/api.hpp"

//ViennaGrid includes:
#include "viennagrid/domain.hpp"

namespace viennafem
{
  
  class linear_pde_options
  {
    public:      
      long data_id() const { return data_id_; }
      void data_id(long new_id) { data_id_ = new_id; }
      
      long trial_space_id() const { return trial_space_id_; }
      void trial_space_id(long new_id) { trial_space_id_ = new_id; }
      
      long test_space_id() const { return test_space_id_; }
      void test_space_id(long new_id) { test_space_id_ = new_id; }
      
    private:
      long data_id_;
      long trial_space_id_;
      long test_space_id_;
  };
  
  template <typename TrialSpaceTag, typename TestSpaceTag>
  linear_pde_options make_linear_pde_options(long data_id, TrialSpaceTag, TestSpaceTag)
  {
    linear_pde_options options;
    options.data_id(data_id);
    options.trial_space_id(viennafem::space_to_id<TrialSpaceTag>::value);
    options.test_space_id(viennafem::space_to_id<TestSpaceTag>::value);
    return options;
  }
  
}
#endif
