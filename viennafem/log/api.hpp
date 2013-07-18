#ifndef VIENNAFEM_LOG_API_HPP
#define VIENNAFEM_LOG_API_HPP

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


#include <vector>
#include "viennafem/forwards.h"
#include "viennagrid/config/domain_config.hpp"

#include "viennafem/log/latex.hpp"

/** @file   viennafem/log/api.hpp
    @brief  Provides functions for logging
*/

namespace viennafem
{
  /** @brief Writes the strong form of the PDE to the logger */
  template <typename PDESystem>
  void log_strong_form(PDESystem const & pde_system)
  {
    pde_system.logger().write_strong_form(pde_system.pdes());
  }

  /** @brief Writes the weak form of the PDE to the logger */
  template <typename EquationArray, typename PDESystem>
  void log_weak_form(EquationArray const & weak_form, PDESystem const & pde_system)
  {
    pde_system.logger().write_weak_form(weak_form);
  }

  /** @brief Writes the weak form (with the coordinate system applied) to the logger */
  template <typename EquationArray, typename PDESystem>
  void log_coordinated_weak_form(EquationArray const & weak_form, PDESystem const & pde_system)
  {
    pde_system.logger().write_coordinated_weak_form(weak_form);
  }

  /** @brief Writes the weak form transformed to the reference cell to the logger */
  template <typename CellType, typename EquationArray, typename PDESystem>
  void log_transformed_weak_form(EquationArray const & weak_form, PDESystem const & pde_system)
  {
    // make cell type known to logger
    // TODO: Better encapsulation
    pde_system.logger().translator().add_processor(new rt_latex_dt_dx_processor<CellType, viennamath::default_interface_type>());

    pde_system.logger().write_transformed_weak_form(weak_form);
  }

  /** @brief Writes the test and trial functions to the logger */
  template <typename EquationArray, typename PDESystem>
  void log_test_and_trial_space(EquationArray const & test_space,
                                EquationArray const & trial_space,
                                PDESystem const & pde_system)
  {
    pde_system.logger().write_test_and_trial_space(test_space, trial_space);
  }


  /*void write_linear_solver_stats(latex_logger & log)
  {

  }*/

}

#endif
