#ifndef VIENNAFEM_LOG_API_HPP
#define VIENNAFEM_LOG_API_HPP

/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */


#include <vector>
#include "viennafem/forwards.h"
#include "viennafem/bases/BFStock.hpp"
#include "viennagrid/forwards.h"
#include "viennagrid/domain.hpp"
#include "viennagrid/segment.hpp"
#include "viennagrid/iterators.hpp"

namespace viennafem
{
  
  template <typename PDESystem>
  void log_strong_form(PDESystem const & pde_system)
  {
    pde_system.logger().write_strong_form(pde_system.pdes());
  }

  template <typename EquationArray, typename PDESystem>
  void log_weak_form(EquationArray const & weak_form, PDESystem const & pde_system)
  {
    pde_system.logger().write_weak_form(weak_form);
  }

  template <typename EquationArray, typename PDESystem>
  void log_coordinated_weak_form(EquationArray const & weak_form, PDESystem const & pde_system)
  {
    pde_system.logger().write_coordinated_weak_form(weak_form);
  }


  template <typename EquationArray, typename PDESystem>
  void log_transformed_weak_form(EquationArray const & weak_form, PDESystem const & pde_system)
  {
    pde_system.logger().write_transformed_weak_form(weak_form);
  }

  template <typename PDESystem>
  void log_test_and_trial_space(PDESystem const & pde_system)
  {
    pde_system.logger().write_test_and_trial_space(pde_system);
  }
  
  
  /*void write_linear_solver_stats(latex_logger & log)
  {
    
  }*/

}

#endif
