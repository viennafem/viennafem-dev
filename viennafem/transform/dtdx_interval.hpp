#ifndef VIENNAFEM_TRANSFORM_INTERVAL_HPP
#define VIENNAFEM_TRANSFORM_INTERVAL_HPP

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

#include <iostream>
#include <utility>
#include "viennagrid/topology/simplex.hpp"
#include "viennagrid/algorithm/spanned_volume.hpp"
#include "viennagrid/config/default_configs.hpp"
#include "viennafem/forwards.h"

/** @file    dtdx_interval.hpp
    @brief   Provides the transformation coefficients of an arbitrary interval to the unit interval
*/

namespace viennafem
{

  template <typename DomainType, typename StorageType>
  struct dt_dx_handler <DomainType, StorageType, viennafem::unit_interval>
  {
  public:
    typedef typename viennagrid::result_of::cell_tag<DomainType>::type                    CellTag;
    typedef typename viennagrid::result_of::element<DomainType, CellTag>::type            CellType;
    typedef typename viennagrid::result_of::point<DomainType>::type                       PointType;
    typedef typename viennagrid::result_of::default_point_accessor<DomainType>::type      PointAccessorType;

    typedef typename viennadata::result_of::accessor<StorageType, det_dF_dt_key,   viennafem::numeric_type, CellType>::type   det_dF_dt_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<0, 0>, viennafem::numeric_type, CellType>::type   dt_dx_key_00_AccessorType;

    dt_dx_handler(DomainType& domain, StorageType& storage) : pnt_acc(viennagrid::default_point_accessor(domain))
    {
      det_dF_dt_acc    = viennadata::make_accessor(storage, det_dF_dt_key());
      dt_dx_key_00_acc = viennadata::make_accessor(storage, dt_dx_key<0, 0>());
    }

    template <typename CellType>
    void operator()(CellType const & cell)
    {
      PointType const & p0 = pnt_acc( viennagrid::vertices(cell)[0] );
      PointType const & p1 = pnt_acc( viennagrid::vertices(cell)[1] );

      //Step 1: store determinant:
      numeric_type x1_x0 = p1[0] - p0[0];

      det_dF_dt_acc(cell)    = x1_x0;
      dt_dx_key_00_acc(cell) = 1.0 / x1_x0;
    }

    PointAccessorType           pnt_acc;
    det_dF_dt_AccessorType      det_dF_dt_acc;
    dt_dx_key_00_AccessorType   dt_dx_key_00_acc;
  };

} //namespace

#endif
