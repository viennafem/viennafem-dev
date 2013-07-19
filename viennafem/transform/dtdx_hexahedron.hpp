#ifndef VIENNAFEM_TRANSFORM_HEXAHEDRON_HPP
#define VIENNAFEM_TRANSFORM_HEXAHEDRON_HPP

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
#include "viennagrid/config/domain_config.hpp"
#include "viennafem/forwards.h"

#include "viennamath/manipulation/simplify.hpp"

/** @file    dtdx_hexahedron.hpp
    @brief   Provides the transformation coefficients of an arbitrary hexahedron to the unit hexahedron
*/

namespace viennafem
{

  template <typename DomainType, typename StorageType>
  struct dt_dx_handler <DomainType, StorageType, viennafem::unit_cube>
  {
    public:

    typedef typename viennagrid::result_of::cell_tag<DomainType>::type                    CellTag;
    typedef typename viennagrid::result_of::element<DomainType, CellTag>::type            CellType;
    typedef typename viennagrid::result_of::point<DomainType>::type                       PointType;
    typedef typename viennagrid::result_of::default_point_accessor<DomainType>::type      PointAccessorType;

    typedef typename viennadata::result_of::accessor<StorageType, det_dF_dt_key,   viennamath::expr, CellType>::type   det_dF_dt_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<0, 0>, viennamath::expr, CellType>::type   dt_dx_key_00_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<0, 1>, viennamath::expr, CellType>::type   dt_dx_key_01_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<0, 2>, viennamath::expr, CellType>::type   dt_dx_key_02_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<1, 0>, viennamath::expr, CellType>::type   dt_dx_key_10_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<1, 1>, viennamath::expr, CellType>::type   dt_dx_key_11_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<1, 2>, viennamath::expr, CellType>::type   dt_dx_key_12_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<2, 0>, viennamath::expr, CellType>::type   dt_dx_key_20_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<2, 1>, viennamath::expr, CellType>::type   dt_dx_key_21_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<2, 2>, viennamath::expr, CellType>::type   dt_dx_key_22_AccessorType;

    dt_dx_handler(DomainType& domain, StorageType& storage) : pnt_acc(viennagrid::default_point_accessor(domain))
    {

      det_dF_dt_acc    = viennadata::make_accessor(storage, det_dF_dt_key());
      dt_dx_key_00_acc = viennadata::make_accessor(storage, dt_dx_key<0, 0>());
      dt_dx_key_01_acc = viennadata::make_accessor(storage, dt_dx_key<0, 1>());
      dt_dx_key_02_acc = viennadata::make_accessor(storage, dt_dx_key<0, 2>());
      dt_dx_key_10_acc = viennadata::make_accessor(storage, dt_dx_key<1, 0>());
      dt_dx_key_11_acc = viennadata::make_accessor(storage, dt_dx_key<1, 1>());
      dt_dx_key_12_acc = viennadata::make_accessor(storage, dt_dx_key<1, 2>());
      dt_dx_key_20_acc = viennadata::make_accessor(storage, dt_dx_key<2, 0>());
      dt_dx_key_21_acc = viennadata::make_accessor(storage, dt_dx_key<2, 1>());
      dt_dx_key_22_acc = viennadata::make_accessor(storage, dt_dx_key<2, 2>());
    }

    template <typename CellType>
    void operator()(CellType const & cell)
    {
      PointType const& p0 = pnt_acc( viennagrid::vertices(cell)[0] );
      PointType const& p1 = pnt_acc( viennagrid::vertices(cell)[1] );
      PointType const& p2 = pnt_acc( viennagrid::vertices(cell)[2] );
      PointType const& p3 = pnt_acc( viennagrid::vertices(cell)[3] );
      PointType const& p4 = pnt_acc( viennagrid::vertices(cell)[4] );
      PointType const& p5 = pnt_acc( viennagrid::vertices(cell)[5] );
      PointType const& p6 = pnt_acc( viennagrid::vertices(cell)[6] );
      PointType const& p7 = pnt_acc( viennagrid::vertices(cell)[7] );

      // Write mapping from local coordinates (xi, eta, nu) to global coordinates (x,y,z) in the form
      //
      // ( x )
      // ( y ) = a + xi * a_0 + eta * a_1 + nu * a_2 + xi * eta * a_01 + xi * nu * a_02 + eta * nu * a_12 + xi * eta * nu * a_123
      // ( z )
      //
      // with vectors a, a_0, etc.
      PointType a_0   = p1 - p0;
      PointType a_1   = p2 - p0;
      PointType a_2   = p4 - p0;
      PointType a_01  = p0 - p1 - p2 + p3;
      PointType a_02  = p0 - p1 - p4 + p5;
      PointType a_12  = p0 - p2 - p4 + p6;
      PointType a_012 = p1 - p0 + p2 - p3 + p4 - p5 - p6 + p7;

      viennamath::variable  xi(0);
      viennamath::variable eta(1);
      viennamath::variable  nu(2);

      viennamath::expr J_00 = a_0[0] + a_01[0] * eta + a_02[0] *  nu + a_012[0] * eta *  nu; viennamath::inplace_simplify(J_00);
      viennamath::expr J_10 = a_0[1] + a_01[1] * eta + a_02[1] *  nu + a_012[1] * eta *  nu; viennamath::inplace_simplify(J_10);
      viennamath::expr J_20 = a_0[2] + a_01[2] * eta + a_02[2] *  nu + a_012[2] * eta *  nu; viennamath::inplace_simplify(J_20);

      viennamath::expr J_01 = a_1[0] + a_01[0] *  xi + a_12[0] *  nu + a_012[0] *  xi *  nu; viennamath::inplace_simplify(J_01);
      viennamath::expr J_11 = a_1[1] + a_01[1] *  xi + a_12[1] *  nu + a_012[1] *  xi *  nu; viennamath::inplace_simplify(J_11);
      viennamath::expr J_21 = a_1[2] + a_01[2] *  xi + a_12[2] *  nu + a_012[2] *  xi *  nu; viennamath::inplace_simplify(J_21);

      viennamath::expr J_02 = a_2[0] + a_02[0] *  xi + a_12[0] * eta + a_012[0] *  xi * eta; viennamath::inplace_simplify(J_02);
      viennamath::expr J_12 = a_2[1] + a_02[1] *  xi + a_12[1] * eta + a_012[1] *  xi * eta; viennamath::inplace_simplify(J_12);
      viennamath::expr J_22 = a_2[2] + a_02[2] *  xi + a_12[2] * eta + a_012[2] *  xi * eta; viennamath::inplace_simplify(J_22);

      // determinant:
      viennamath::expr det_J =   J_00 * J_11 * J_22 + J_01 * J_12 * J_20 + J_02 * J_10 * J_21
                               - J_20 * J_11 * J_02 - J_21 * J_12 * J_00 - J_22 * J_10 * J_01;

      viennamath::inplace_simplify(det_J);

      det_dF_dt_acc(cell) = det_J;

      //Step 2: store partial derivatives:
      typedef viennamath::expr   ValueType;

      dt_dx_key_00_acc(cell) = viennamath::simplify(J_11*J_22 - J_21*J_12) / det_J;
      dt_dx_key_01_acc(cell) = viennamath::simplify(J_21*J_02 - J_01*J_22) / det_J;
      dt_dx_key_02_acc(cell) = viennamath::simplify(J_01*J_12 - J_11*J_02) / det_J;

      dt_dx_key_10_acc(cell) = viennamath::simplify(J_20*J_12 - J_10*J_22) / det_J;
      dt_dx_key_11_acc(cell) = viennamath::simplify(J_00*J_22 - J_20*J_02) / det_J;
      dt_dx_key_12_acc(cell) = viennamath::simplify(J_10*J_02 - J_00*J_12) / det_J;

      dt_dx_key_20_acc(cell) = viennamath::simplify(J_10*J_21 - J_20*J_11) / det_J;
      dt_dx_key_21_acc(cell) = viennamath::simplify(J_20*J_01 - J_00*J_21) / det_J;
      dt_dx_key_22_acc(cell) = viennamath::simplify(J_00*J_11 - J_10*J_01) / det_J;
    }

    PointAccessorType           pnt_acc;
    det_dF_dt_AccessorType      det_dF_dt_acc;
    dt_dx_key_00_AccessorType   dt_dx_key_00_acc;
    dt_dx_key_01_AccessorType   dt_dx_key_01_acc;
    dt_dx_key_02_AccessorType   dt_dx_key_02_acc;
    dt_dx_key_10_AccessorType   dt_dx_key_10_acc;
    dt_dx_key_11_AccessorType   dt_dx_key_11_acc;
    dt_dx_key_12_AccessorType   dt_dx_key_12_acc;
    dt_dx_key_20_AccessorType   dt_dx_key_20_acc;
    dt_dx_key_21_AccessorType   dt_dx_key_21_acc;
    dt_dx_key_22_AccessorType   dt_dx_key_22_acc;
  };


} //namespace

#endif
