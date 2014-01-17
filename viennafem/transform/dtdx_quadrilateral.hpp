#ifndef VIENNAFEM_TRANSFORM_QUADRILATERAL_HPP
#define VIENNAFEM_TRANSFORM_QUADRILATERAL_HPP

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

#include <iostream>
#include <utility>
#include "viennagrid/topology/simplex.hpp"
#include "viennagrid/algorithm/spanned_volume.hpp"
#include "viennagrid/config/mesh_config.hpp"
#include "viennafem/forwards.h"
#include "viennamath/manipulation/simplify.hpp"

/** @file    dtdx_quadrilateral.hpp
    @brief   Provides the transformation coefficients of an arbitrary quadrilateral to the unit square
*/

namespace viennafem
{
  template <typename DomainType, typename StorageType>
  struct dt_dx_handler <DomainType, StorageType, viennafem::unit_square>
  {
  public:
    typedef typename viennagrid::result_of::cell_tag<DomainType>::type                    CellTag;
    typedef typename viennagrid::result_of::element<DomainType, CellTag>::type            CellType;
    typedef typename viennagrid::result_of::point<DomainType>::type                       PointType;
    typedef typename viennagrid::result_of::default_point_accessor<DomainType>::type      PointAccessorType;

    typedef typename viennadata::result_of::accessor<StorageType, det_dF_dt_key,   viennamath::expr, CellType>::type   det_dF_dt_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<0, 0>, viennamath::expr, CellType>::type   dt_dx_key_00_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<0, 1>, viennamath::expr, CellType>::type   dt_dx_key_01_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<1, 0>, viennamath::expr, CellType>::type   dt_dx_key_10_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<1, 1>, viennamath::expr, CellType>::type   dt_dx_key_11_AccessorType;

    dt_dx_handler(DomainType& domain, StorageType& storage) : pnt_acc(viennagrid::default_point_accessor(domain))
    {
      det_dF_dt_acc    = viennadata::make_accessor(storage, det_dF_dt_key());
      dt_dx_key_00_acc = viennadata::make_accessor(storage, dt_dx_key<0, 0>());
      dt_dx_key_01_acc = viennadata::make_accessor(storage, dt_dx_key<0, 1>());
      dt_dx_key_10_acc = viennadata::make_accessor(storage, dt_dx_key<1, 0>());
      dt_dx_key_11_acc = viennadata::make_accessor(storage, dt_dx_key<1, 1>());
    }

    template <typename CellType>
    void operator()(CellType const & cell)
    {
      PointType const& p0 = pnt_acc( viennagrid::vertices(cell)[0] );
      PointType const& p1 = pnt_acc( viennagrid::vertices(cell)[1] );
      PointType const& p2 = pnt_acc( viennagrid::vertices(cell)[2] );
      PointType const& p3 = pnt_acc( viennagrid::vertices(cell)[3] );

      //Step 1: store determinant:
      numeric_type x1_x0 = p1[0] - p0[0];
      numeric_type x2_x0 = p2[0] - p0[0];
      numeric_type y1_y0 = p1[1] - p0[1];
      numeric_type y2_y0 = p2[1] - p0[1];
      numeric_type coeff_x = p0[0] - p1[0] - p2[0] + p3[0];
      numeric_type coeff_y = p0[1] - p1[1] - p2[1] + p3[1];

      viennamath::variable  xi(0);
      viennamath::variable eta(1);

      viennamath::expr det_J = (x1_x0 + eta * coeff_x) * (y2_y0 + xi * coeff_y) - (y1_y0 + eta * coeff_y) * (x2_x0 + xi * coeff_x);
      viennamath::inplace_simplify(det_J);

      det_dF_dt_acc(cell) = det_J;

      //Step 2: store partial derivatives:
      dt_dx_key_00_acc(cell) = ( y2_y0 +  xi * coeff_y) / det_J;
      dt_dx_key_01_acc(cell) = (-x2_x0 -  xi * coeff_x) / det_J;
      dt_dx_key_10_acc(cell) = (-y1_y0 - eta * coeff_y) / det_J;
      dt_dx_key_11_acc(cell) = ( x1_x0 + eta * coeff_x) / det_J;
    }

    PointAccessorType           pnt_acc;
    det_dF_dt_AccessorType      det_dF_dt_acc;
    dt_dx_key_00_AccessorType   dt_dx_key_00_acc;
    dt_dx_key_01_AccessorType   dt_dx_key_01_acc;
    dt_dx_key_10_AccessorType   dt_dx_key_10_acc;
    dt_dx_key_11_AccessorType   dt_dx_key_11_acc;

  };


} //namespace

#endif
