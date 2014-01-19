#ifndef VIENNAFEM_LINEAR_PDE_OPTIONS_HPP
#define VIENNAFEM_LINEAR_PDE_OPTIONS_HPP

/* =========================================================================
   Copyright (c) 2012-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
                             -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                             -----------------

   Author:     Karl Rupp                          rupp@iue.tuwien.ac.at

   License:    MIT (X11), see file LICENSE in the ViennaFEM base directory
============================================================================ */

#include <vector>

//ViennaFEM includes:
#include "viennafem/forwards.h"

/** @file   linear_pde_options.hpp
    @brief  Provides a set of options (i.e. metainformation) for each PDE specified via ViennaMath
*/

namespace viennafem
{
  class linear_pde_options
  {
    public:
      /** @brief Returns the ID used for storing mapping and boundary information */
      long data_id() const { return data_id_; }
      /** @brief Sets the ID used for storing mapping and boundary information */
      void data_id(long new_id) { data_id_ = new_id; }

      /** @brief Returns the ID of the trial space to be used */
      long trial_space_id() const { return trial_space_id_; }
      /** @brief Sets the ID of the trial space to be used */
      void trial_space_id(long new_id) { trial_space_id_ = new_id; }

      /** @brief Returns the ID of the test space to be used */
      long test_space_id() const { return test_space_id_; }
      /** @brief Sets the ID of the test space to be used */
      void test_space_id(long new_id) { test_space_id_ = new_id; }

      /** @brief Specifies whether an existing mapping from a previous FEM run should be reused */
      bool check_existing_mapping() const { return check_mapping_; }
      /** @brief Specifies whether an existing mapping from a previous FEM run should be reused */
      void check_existing_mapping(bool b) { check_mapping_ = b; }

    private:
      long data_id_;
      long trial_space_id_;
      long test_space_id_;
      bool check_mapping_;
  };

  /** @brief Convenience function for creating a PDE config object. */
  template <typename TrialSpaceTag, typename TestSpaceTag>
  linear_pde_options make_linear_pde_options(long data_id, TrialSpaceTag, TestSpaceTag, bool existing_mapping = false)
  {
    linear_pde_options options;
    options.data_id(data_id);
    options.trial_space_id(viennafem::space_to_id<TrialSpaceTag>::value);
    options.test_space_id(viennafem::space_to_id<TestSpaceTag>::value);
    options.check_existing_mapping(existing_mapping);
    return options;
  }

}
#endif
