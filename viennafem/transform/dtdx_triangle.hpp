#ifndef VIENNAFEM_TRANSFORM_TRIANGLE_HPP
#define VIENNAFEM_TRANSFORM_TRIANGLE_HPP

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
#include "viennafem/forwards.h"
#include "viennagrid/topology/simplex.hpp"
#include "viennagrid/algorithm/spanned_volume.hpp"
#include "viennagrid/config/domain_config.hpp"
#include "viennagrid/point.hpp"


/** @file    dtdx_triangle.hpp
    @brief   Provides the transformation coefficients of an arbitrary triangle to the unit triangle
*/

namespace viennafem
{

  //memory-intensive: Compute them once and store the computed values until next update
  template <typename DomainType, typename StorageType>
  struct dt_dx_handler <DomainType, StorageType, viennafem::unit_triangle>
  {
    typedef typename viennagrid::result_of::cell_tag<DomainType>::type                    CellTag;
    typedef typename viennagrid::result_of::element<DomainType, CellTag>::type            CellType;
    typedef typename viennagrid::result_of::point<DomainType>::type                       PointType;
    typedef typename viennagrid::result_of::default_point_accessor<DomainType>::type      PointAccessorType;

    typedef typename viennadata::result_of::accessor<StorageType, det_dF_dt_key,   viennafem::numeric_type, CellType>::type   det_dF_dt_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<0, 0>, viennafem::numeric_type, CellType>::type   dt_dx_key_00_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<0, 1>, viennafem::numeric_type, CellType>::type   dt_dx_key_01_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<1, 0>, viennafem::numeric_type, CellType>::type   dt_dx_key_10_AccessorType;
    typedef typename viennadata::result_of::accessor<StorageType, dt_dx_key<1, 1>, viennafem::numeric_type, CellType>::type   dt_dx_key_11_AccessorType;

    dt_dx_handler(DomainType& domain, StorageType& storage) : pnt_acc(viennagrid::default_point_accessor(domain))
    {
      det_dF_dt_acc    = viennadata::make_accessor<det_dF_dt_key,   viennafem::numeric_type, CellType>(storage, det_dF_dt_key());
      dt_dx_key_00_acc = viennadata::make_accessor<dt_dx_key<0, 0>, viennafem::numeric_type, CellType>(storage, dt_dx_key<0, 0>());
      dt_dx_key_01_acc = viennadata::make_accessor<dt_dx_key<0, 1>, viennafem::numeric_type, CellType>(storage, dt_dx_key<0, 1>());
      dt_dx_key_10_acc = viennadata::make_accessor<dt_dx_key<1, 0>, viennafem::numeric_type, CellType>(storage, dt_dx_key<1, 0>());
      dt_dx_key_11_acc = viennadata::make_accessor<dt_dx_key<1, 1>, viennafem::numeric_type, CellType>(storage, dt_dx_key<1, 1>());
    }

    template <typename CellType>
    void operator()(CellType const & cell)
    {
      PointType const& p0 = pnt_acc( viennagrid::vertices(cell)[0] );
      PointType const& p1 = pnt_acc( viennagrid::vertices(cell)[1] );
      PointType const& p2 = pnt_acc( viennagrid::vertices(cell)[2] );

      //Step 1: store determinant:
      double det_dF_dt = 2.0 * viennagrid::spanned_volume(p0, p1, p2);

      assert(det_dF_dt > 0);
      //std::cout << "dt_dx triangle!" << std::endl;

      det_dF_dt_acc(cell) = det_dF_dt;

      //Step 2: store partial derivatives:
      dt_dx_key_00_acc(cell) =   (p2[1] - p0[1]) / det_dF_dt;
      dt_dx_key_01_acc(cell) = - (p2[0] - p0[0]) / det_dF_dt;
      dt_dx_key_10_acc(cell) = - (p1[1] - p0[1]) / det_dF_dt;
      dt_dx_key_11_acc(cell) =   (p1[0] - p0[0]) / det_dF_dt;
    }

    PointAccessorType           pnt_acc;
    det_dF_dt_AccessorType      det_dF_dt_acc;
    dt_dx_key_00_AccessorType   dt_dx_key_00_acc;
    dt_dx_key_01_AccessorType   dt_dx_key_01_acc;
    dt_dx_key_10_AccessorType   dt_dx_key_10_acc;
    dt_dx_key_11_AccessorType   dt_dx_key_11_acc;

  };

  /*
  template <typename T_Configuration>
  struct dt_dx_handler<T_Configuration, TriangleTag, DtDxStoreDetOnly>
  {
    typedef typename T_Configuration::CoordType                     ScalarType;
    typedef typename DomainTypes<T_Configuration>::PointType        PointType;
    typedef typename DomainTypes<T_Configuration>::VertexType       VertexType;

    public:
      dt_dx_handler() : Base() {};
      dt_dx_handler( const dt_dx_handler & ddh ) : Base(ddh) {};

      //returns the element dt_i/dx_j of the functional determinant induced by the mapping to the reference element. i and j start at 0.
      ScalarType get_dt_dx(int i, int j) const
      {
        PointType & p0 = Base::vertices_[0]->getPoint();
        double ret = 0.0;

        if (i == 0)
        {
          PointType & p2 = Base::vertices_[2]->getPoint();
          if (j == 0)
            ret = (p2.get_y() - p0.get_y()) / det_dF_dt;
          else if (j == 1)
            ret = - (p2.get_x() - p0.get_x()) / det_dF_dt;
          else
            std::cerr << "ERROR: Accessing invalid elements of functional determinant!!";
        }
        else if (i == 1)
        {
          PointType & p1 = Base::vertices_[1]->getPoint();
          if (j == 0)
            ret = - (p1.get_y() - p0.get_y()) / det_dF_dt;
          else if (j == 1)
            ret = (p1.get_x() - p0.get_x()) / det_dF_dt;
          else
            std::cerr << "ERROR: Accessing invalid elements of functional determinant!!";
        }

        return ret; //dt_dx[2*i + j];
      }

      ScalarType get_det_dF_dt() const
      {
        return det_dF_dt;
      }

      void update_dt_dx()
      {
        checkOrientation();
      };

      void init_dt_dx() {}

    protected:

      void checkOrientation()
      {
        PointType & p0 = Base::vertices_[0]->getPoint();
        PointType & p1 = Base::vertices_[1]->getPoint();
        PointType & p2 = Base::vertices_[2]->getPoint();

        //volume:
        det_dF_dt = spannedVolume(p1 - p0, p2 - p0);

        if (det_dF_dt < 0)
        {
          det_dF_dt = - det_dF_dt;
          VertexType *temp = Base::vertices_[0];
          Base::vertices_[0] = Base::vertices_[1];
          Base::vertices_[1] = temp;
        }
        else if (det_dF_dt == 0.0)
          std::cout << "ERROR: detected degenerated element!" << std::endl;
      }

      void print(long indent) const
      {
        for (long i = 0; i<indent; ++i)
          std::cout << "   ";
        std::cout << "* dt-dx-Handler: StoreDetOnly" << std::endl;
        Base::print(indent);
      }

    private:
      ScalarType det_dF_dt;                         //determinant of Jacobian matrix
  };

  //save as much memory as possible: compute all values on access
  //however, in general the DtDxStoreDetOnly is much faster while having only moderate additional memory requirements.
  template <typename T_Configuration>
  struct dt_dx_handler<T_Configuration, TriangleTag, DtDxOnAccess>
    : public lower_level_holder<T_Configuration, TriangleTag, TriangleTag::TopoLevel - 1>
  {
    typedef lower_level_holder<T_Configuration, TriangleTag, TriangleTag::TopoLevel - 1>   Base;
    typedef typename T_Configuration::CoordType                     ScalarType;
    typedef typename DomainTypes<T_Configuration>::PointType        PointType;
    typedef typename DomainTypes<T_Configuration>::VertexType       VertexType;

    public:
      dt_dx_handler() : Base() {};
      dt_dx_handler( const dt_dx_handler & ddh ) : Base(ddh) {};

      //returns the element dt_i/dx_j of the functional determinant induced by the mapping to the reference element. i and j start at 0.
      ScalarType get_dt_dx(int i, int j) const
      {
        PointType & p0 = Base::vertices_[0]->getPoint();
        PointType & p1 = Base::vertices_[1]->getPoint();
        PointType & p2 = Base::vertices_[2]->getPoint();
        double ret = 0.0;
        double det_dF_dt = spannedVolume(p1 - p0, p2 - p0);

        if (i == 0)
        {
          if (j == 0)
            ret = (p2.get_y() - p0.get_y()) / det_dF_dt;
          else if (j == 1)
            ret = - (p2.get_x() - p0.get_x()) / det_dF_dt;
          else
            std::cerr << "ERROR: Accessing invalid elements of functional determinant!!";
        }
        else if (i == 1)
        {
          if (j == 0)
            ret = - (p1.get_y() - p0.get_y()) / det_dF_dt;
          else if (j == 1)
            ret = (p1.get_x() - p0.get_x()) / det_dF_dt;
          else
            std::cerr << "ERROR: Accessing invalid elements of functional determinant!!";
        }

        return ret; //dt_dx[2*i + j];
      }

      ScalarType get_det_dF_dt() const
      {
        PointType & p0 = Base::vertices_[0]->getPoint();
        PointType & p1 = Base::vertices_[1]->getPoint();
        PointType & p2 = Base::vertices_[2]->getPoint();
        return spannedVolume(p1 - p0, p2 - p0);
      }

      void update_dt_dx()
      {
        checkOrientation();
      };

      void init_dt_dx() {}

    protected:

      void checkOrientation()
      {
        PointType & p0 = Base::vertices_[0]->getPoint();
        PointType & p1 = Base::vertices_[1]->getPoint();
        PointType & p2 = Base::vertices_[2]->getPoint();

        //volume:
        double det_dF_dt = spannedVolume(p1 - p0, p2 - p0);

        if (det_dF_dt < 0)
        {
          det_dF_dt = - det_dF_dt;
          VertexType *temp = Base::vertices_[0];
          Base::vertices_[0] = Base::vertices_[1];
          Base::vertices_[1] = temp;
        }
      }

      void print(long indent = 0) const
      {
        for (long i = 0; i<indent; ++i)
          std::cout << "   ";
        std::cout << "* dt-dx-Handler: OnAccess" << std::endl;
        Base::print(indent);
      }

  };

  //fast and memory-saving: Compute an element's Jacobian as soon as the first entry is accessed. Use static memory, so that only one Jacobian set is stored at the same time. Must not be used in multi-threaded applications!!
  template <typename T_Configuration>
  struct dt_dx_handler<T_Configuration, TriangleTag, DtDxStoreStatically>
    : public lower_level_holder<T_Configuration, TriangleTag, TriangleTag::TopoLevel - 1>
  {
    typedef lower_level_holder<T_Configuration, TriangleTag, TriangleTag::TopoLevel - 1>   Base;
    typedef typename T_Configuration::CoordType                     ScalarType;
    typedef typename DomainTypes<T_Configuration>::PointType        PointType;
    typedef typename DomainTypes<T_Configuration>::VertexType       VertexType;

    public:
      dt_dx_handler() : Base() {};
      dt_dx_handler( const dt_dx_handler & ddh ) : Base(ddh) {};

      //returns the element dt_i/dx_j of the functional determinant induced by the mapping to the reference element. i and j start at 0.
      ScalarType get_dt_dx(int i, int j) const
      {
        return dt_dx[2*i + j];
      }

      ScalarType get_det_dF_dt() const
      {
        return det_dF_dt;
      }

      void update_dt_dx()
      {
        checkOrientation();
      };

      void init_dt_dx()
      {
        computeCellJacobian();
      }

    protected:

      void computeCellJacobian()
      {
        PointType & p0 = Base::vertices_[0]->getPoint();
        PointType & p1 = Base::vertices_[1]->getPoint();
        PointType & p2 = Base::vertices_[2]->getPoint();

        det_dF_dt = spannedVolume(p1 - p0, p2 - p0);

        if (det_dF_dt != 0)
        {
          //dt_1/dx_1
          dt_dx[0] = ( p2.get_y() - p0.get_y()) / det_dF_dt;
          //dt_1/dx_2
          dt_dx[1] = - (p2.get_x() - p0.get_x()) / det_dF_dt;

          //dt_2/dx_1
          dt_dx[2] = - (p1.get_y() - p0.get_y()) / det_dF_dt;
          //dt_2/dx_2
          dt_dx[3] = ( p1.get_x() - p0.get_x()) / det_dF_dt;
        }

      } //computeCellJacobian()

      void checkOrientation()
      {
        PointType & p0 = Base::vertices_[0]->getPoint();
        PointType & p1 = Base::vertices_[1]->getPoint();
        PointType & p2 = Base::vertices_[2]->getPoint();

        //volume:
        ScalarType detdFdt = spannedVolume(p1 - p0, p2 - p0);

        if (detdFdt < 0)
        {
          VertexType *temp = Base::vertices_[0];
          Base::vertices_[0] = Base::vertices_[1];
          Base::vertices_[1] = temp;
        }
        else if (detdFdt == 0.0)
        {
          std::cout << "ERROR: detected degenerated element!" << std::endl;
          Base::vertices_[0]->print();
          Base::vertices_[1]->print();
          Base::vertices_[2]->print();
          std::cout << "Enter char to continue" << std::endl;
          char dummy;
          std::cin >> dummy;
        }
      }

      void print(long indent) const
      {
        for (long i = 0; i<indent; ++i)
          std::cout << "   ";
        std::cout << "* dt-dx-Handler: StoreStatically" << std::endl;
        Base::print(indent);
      }

    private:
      static ScalarType dt_dx[4];                          //inverse of Jacobian matrix from mapping
      static ScalarType det_dF_dt;                         //determinant of Jacobian matrix
  };

  //initialize static variables:
  template <typename T_Configuration>
  typename T_Configuration::CoordType
  dt_dx_handler<T_Configuration, TriangleTag, DtDxStoreStatically>::det_dF_dt = 0;

  template <typename T_Configuration>
  typename T_Configuration::CoordType
  dt_dx_handler<T_Configuration, TriangleTag, DtDxStoreStatically>::dt_dx[] = {0, 0, 0, 0};
 */

} //namespace

#endif
