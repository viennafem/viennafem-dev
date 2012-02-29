#ifndef VIENNAFEM_TRANSFORM_TETRAHEDRON_HPP
#define VIENNAFEM_TRANSFORM_TETRAHEDRON_HPP

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
#include "viennagrid/topology/tetrahedron.hpp"
#include "viennagrid/algorithm/spanned_volume.hpp"
#include "viennagrid/domain.hpp"
#include "viennafem/forwards.h"

/** @file    dtdx_tetrahedron.hpp
    @brief   Provides the transformation coefficients of an arbitrary triangle to the unit tetrahedron
*/

namespace viennafem
{
  
  //memory-intensive: Compute them once and store the computed values until next update
  template <>
  struct dt_dx_handler<viennafem::unit_tetrahedron>
  {
    public:
      
      template <typename CellType>
      static void apply(CellType const & cell)
      {
        typedef typename CellType::config_type       Config;
        typedef typename viennagrid::result_of::point<Config>::type   PointType;
        
        PointType p0 = viennagrid::ncells<0>(cell)[0].point();
        PointType p1 = viennagrid::ncells<0>(cell)[1].point() - p0;
        PointType p2 = viennagrid::ncells<0>(cell)[2].point() - p0;
        PointType p3 = viennagrid::ncells<0>(cell)[3].point() - p0;
        
        //Step 1: store determinant:
        numeric_type det_dF_dt = 6.0 * viennagrid::spanned_volume(viennagrid::ncells<0>(cell)[0].point(),
                                                                  viennagrid::ncells<0>(cell)[1].point(),
                                                                  viennagrid::ncells<0>(cell)[2].point(),
                                                                  viennagrid::ncells<0>(cell)[3].point());
        
        viennadata::access<det_dF_dt_key, numeric_type>()(cell) = det_dF_dt;
        
        //Step 2: store partial derivatives:
        viennadata::access<dt_dx_key<0, 0>, numeric_type>()(cell) = (  + p2[1] * p3[2] - p2[2] * p3[1] ) / det_dF_dt;
        viennadata::access<dt_dx_key<0, 1>, numeric_type>()(cell) = (  - p2[0] * p3[2] + p2[2] * p3[0] ) / det_dF_dt;
        viennadata::access<dt_dx_key<0, 2>, numeric_type>()(cell) = (  + p2[0] * p3[1] - p2[1] * p3[0] ) / det_dF_dt;
        
        viennadata::access<dt_dx_key<1, 0>, numeric_type>()(cell) = (  - p1[1] * p3[2] + p1[2] * p3[1] ) / det_dF_dt;
        viennadata::access<dt_dx_key<1, 1>, numeric_type>()(cell) = (  + p1[0] * p3[2] - p1[2] * p3[0] ) / det_dF_dt;
        viennadata::access<dt_dx_key<1, 2>, numeric_type>()(cell) = (  - p1[0] * p3[1] + p1[1] * p3[0] ) / det_dF_dt;
        
        viennadata::access<dt_dx_key<2, 0>, numeric_type>()(cell) = (  + p1[1] * p2[2] - p1[2] * p2[1] ) / det_dF_dt;
        viennadata::access<dt_dx_key<2, 1>, numeric_type>()(cell) = (  - p1[0] * p2[2] + p1[2] * p2[0] ) / det_dF_dt;
        viennadata::access<dt_dx_key<2, 2>, numeric_type>()(cell) = (  + p1[0] * p2[1] - p1[1] * p2[0] ) / det_dF_dt;
      }

  };
  
  
  
/*  Old code to follow...
 

  //memory-intensive: Compute them once and store the computed values until next update
  template <typename T_Configuration>
  struct dt_dx_handler<T_Configuration, viennagrid::tetrahedron_tag, DtDxStoreAll>
  {
    typedef typename T_Configuration::CoordType                     ScalarType;
    typedef typename viennagrid::DomainTypes<T_Configuration>::point_type        PointType;
    typedef typename viennagrid::DomainTypes<T_Configuration>::vertex_type       VertexType;

    public:

      //returns the element dt_i/dx_j of the functional determinant induced by the mapping to the reference element. i and j start at 0.
      ScalarType get_dt_dx(int i, int j) const
      {
        return dt_dx[3*i + j];
      }

      ScalarType get_det_dF_dt() const
      {
        return det_dF_dt;
      }

      template <typename CellType>
      void update_dt_dx(CellType const & cell)
      {
        //checkOrientation();

        PointType p0 = cell.getPoint(0);
        PointType p1 = cell.getPoint(1) - p0;
        PointType p2 = cell.getPoint(2) - p0;
        PointType p3 = cell.getPoint(3) - p0;

        //dt_1/dx_1
        dt_dx[0] = (  + p2[1] * p3[2] - p2[2] * p3[1] ) / det_dF_dt;
        //dt_1/dx_2
        dt_dx[1] = (  - p2[0] * p3[2] + p2[2] * p3[0] ) / det_dF_dt;
        //dt_1/dx_3
        dt_dx[2] = (  + p2[0] * p3[1] - p2[1] * p3[0] ) / det_dF_dt;

        //dt_2/dx_1
        dt_dx[3] = (  - p1[1] * p3[2] + p1[2] * p3[1] ) / det_dF_dt;
        //dt_2/dx_2
        dt_dx[4] = (  + p1[0] * p3[2] - p1[2] * p3[0] ) / det_dF_dt;
        //dt_2/dx_3
        dt_dx[5] = (  - p1[0] * p3[1] + p1[1] * p3[0] ) / det_dF_dt;

        //dt_3/dx_1
        dt_dx[6] = (  + p1[1] * p2[2] - p1[2] * p2[1] ) / det_dF_dt;
        //dt_3/dx_2
        dt_dx[7] = (  - p1[0] * p2[2] + p1[2] * p2[0] ) / det_dF_dt;
        //dt_3/dx_3
        dt_dx[8] = (  + p1[0] * p2[1] - p1[1] * p2[0] ) / det_dF_dt;

      }

      void init_dt_dx() {}

    protected:

      / *
      void checkOrientation()
      {
        PointType & p0 = cell.getPoint(0);
        PointType & p1 = cell.getPoint(1) - p0;
        PointType & p2 = cell.getPoint(2) - p0;
        PointType & p3 = cell.getPoint(3) - p0;

        //volume:
        det_dF_dt = spannedVolume(p1 - p0, p2 - p0, p3 - p0);

        if (det_dF_dt < 0)
        {
          det_dF_dt = - det_dF_dt;
          VertexType *temp = Base::vertices_[0];
          Base::vertices_[0] = Base::vertices_[1];
          Base::vertices_[1] = temp;
        }
        else if (det_dF_dt == 0.0)
          std::cout << "ERROR: detected degenerated element!" << std::endl;
      } * /

      void print(long indent) const
      {
        for (long i = 0; i<indent; ++i)
          std::cout << "   ";
        std::cout << "* dt-dx-Handler: StoreAll" << std::endl;
      }

    private:
      ScalarType dt_dx[9];                          //inverse of Jacobian matrix from mapping
      ScalarType det_dF_dt;                         //determinant of Jacobian matrix
  };

  //memory-computation-tradeoff: Store value of Jacobian only
  template <typename T_Configuration>
  struct dt_dx_handler<T_Configuration, viennagrid::tetrahedron_tag, DtDxStoreDetOnly>
  {
    typedef typename T_Configuration::numeric_type                  ScalarType;
    typedef typename viennagrid::DomainTypes<T_Configuration>::point_type        PointType;
    typedef typename viennagrid::DomainTypes<T_Configuration>::vertex_type       VertexType;

    public:

      //returns the element dt_i/dx_j of the functional determinant induced by the mapping to the reference element. i and j start at 0.
      template <typename CellType>
      ScalarType get_dt_dx(CellType const & cell, int i, int j) const
      {
        PointType p0 = cell.getPoint(0);
        PointType p1 = cell.getPoint(1) - p0;
        PointType p2 = cell.getPoint(2) - p0;
        PointType p3 = cell.getPoint(3) - p0;

        double ret = 0.0;

        //dt_1/dx_1
        if (i==0)
        {
          if (j==0)
            ret = (  + p2[1] * p3[2] - p2[2] * p3[1] ) / det_dF_dt;
          else if (j==1)
            ret = (  - p2[0] * p3[2] + p2[2] * p3[0] ) / det_dF_dt;
          else if (j==2)
            ret = (  + p2[0] * p3[1] - p2[1] * p3[0] ) / det_dF_dt;
          else
            std::cerr << "ERROR: Accessing invalid elements of functional determinant!!";
        }
        else if (i==1)
        {
          if (j==0)
            ret = (  - p1[1] * p3[2] + p1[2] * p3[1] ) / det_dF_dt;
          else if (j==1)
            ret = (  + p1[0] * p3[2] - p1[2] * p3[0] ) / det_dF_dt;
          else if (j==2)
            ret = (  - p1[0] * p3[1] + p1[1] * p3[0] ) / det_dF_dt;
          else
            std::cerr << "ERROR: Accessing invalid elements of functional determinant!!";
        }
        else if (i==2)
        {
          if (j==0)
            ret = (  + p1[1] * p2[2] - p1[2] * p2[1] ) / det_dF_dt;
          else if (j==1)
            ret = (  - p1[0] * p2[2] + p1[2] * p2[0] ) / det_dF_dt;
          else if (j==2)
            ret = (  + p1[0] * p2[1] - p1[1] * p2[0] ) / det_dF_dt;
          else
            std::cerr << "ERROR: Accessing invalid elements of functional determinant!!";
        }
        else
            std::cerr << "ERROR: Accessing invalid elements of functional determinant!!";

        return ret;
      }

      ScalarType get_det_dF_dt() const
      {
        return det_dF_dt;
      }

      void update_dt_dx()
      {
        //checkOrientation();
      }

      void init_dt_dx() {}

    protected:

      / *
      void checkOrientation()
      {
        PointType & p0 = Base::vertices_[0]->getPoint();
        PointType & p1 = Base::vertices_[1]->getPoint();
        PointType & p2 = Base::vertices_[2]->getPoint();
        PointType & p3 = Base::vertices_[3]->getPoint();

        //volume:
        det_dF_dt = spannedVolume(p1 - p0, p2 - p0, p3 - p0);

        if (det_dF_dt < 0)
        {
          det_dF_dt = - det_dF_dt;
          VertexType *temp = Base::vertices_[0];
          Base::vertices_[0] = Base::vertices_[1];
          Base::vertices_[1] = temp;
        }
        else if (det_dF_dt == 0.0)
          std::cout << "ERROR: detected degenerated element!" << std::endl;
      } * /

      void print(long indent) const
      {
        for (long i = 0; i<indent; ++i)
          std::cout << "   ";
        std::cout << "* dt-dx-Handler: StoreDetOnly"<< std::endl;
      }

    private:
      ScalarType det_dF_dt;                         //determinant of Jacobian matrix
  };

  //save as much memory as possible: compute all values on access
  //however, in general the DtDxStoreDetOnly is much faster while having only moderate additional memory requirements.
  template <typename T_Configuration>
  struct dt_dx_handler<T_Configuration, viennagrid::tetrahedron_tag, DtDxOnAccess>
  {
    typedef typename T_Configuration::numeric_type                     ScalarType;
    typedef typename viennagrid::DomainTypes<T_Configuration>::point_type        PointType;
    typedef typename viennagrid::DomainTypes<T_Configuration>::vertex_type       VertexType;

    public:

      //returns the element dt_i/dx_j of the functional determinant induced by the mapping to the reference element. i and j start at 0.
      template <typename CellType>
      ScalarType get_dt_dx(CellType const & cell, int i, int j) const
      {
        PointType p0 = cell.getPoint(0);
        PointType p1 = cell.getPoint(1) - p0;
        PointType p2 = cell.getPoint(2) - p0;
        PointType p3 = cell.getPoint(3) - p0;

        double ret = 0.0;
        double det_dF_dt = get_det_dF_dt();

        //dt_1/dx_1
        if (i==0)
        {
          if (j==0)
            ret = (  + p2[1] * p3[2] - p2[2] * p3[1] ) / det_dF_dt;
          else if (j==1)
            ret = (  - p2[0] * p3[2] + p2[2] * p3[0] ) / det_dF_dt;
          else if (j==2)
            ret = (  + p2[0] * p3[1] - p2[1] * p3[0] ) / det_dF_dt;
          else
            std::cerr << "ERROR: Accessing invalid elements of functional determinant!!";
        }
        else if (i==1)
        {
          if (j==0)
            ret = (  - p1[1] * p3[2] + p1[2] * p3[1] ) / det_dF_dt;
          else if (j==1)
            ret = (  + p1[0] * p3[2] - p1[2] * p3[0] ) / det_dF_dt;
          else if (j==2)
            ret = (  - p1[0] * p3[1] + p1[1] * p3[0] ) / det_dF_dt;
          else
            std::cerr << "ERROR: Accessing invalid elements of functional determinant!!";
        }
        else if (i==2)
        {
          if (j==0)
            ret = (  + p1[1] * p2[2] - p1[2] * p2[1] ) / det_dF_dt;
          else if (j==1)
            ret = (  - p1[0] * p2[2] + p1[2] * p2[0] ) / det_dF_dt;
          else if (j==2)
            ret = (  + p1[0] * p2[1] - p1[1] * p2[0] ) / det_dF_dt;
          else
            std::cerr << "ERROR: Accessing invalid elements of functional determinant!!";
        }
        else
            std::cerr << "ERROR: Accessing invalid elements of functional determinant!!";

        return ret;
      }

      template <typename CellType> 
      ScalarType get_det_dF_dt(CellType const & cell) const
      {
        PointType p0 = cell.getPoint(0);
        PointType p1 = cell.getPoint(1);
        PointType p2 = cell.getPoint(2);
        PointType p3 = cell.getPoint(3);

        //volume:
        return spannedVolume(p1 - p0, p2 - p0, p3 - p0);
      }

    protected:

      / *
      void checkOrientation()
      {
        PointType & p0 = Base::vertices_[0]->getPoint();
        PointType & p1 = Base::vertices_[1]->getPoint();
        PointType & p2 = Base::vertices_[2]->getPoint();
        PointType & p3 = Base::vertices_[3]->getPoint();

        //volume:
        double det_dF_dt = spannedVolume(p1 - p0, p2 - p0, p3 - p0);

        if (det_dF_dt < 0)
        {
          det_dF_dt = - det_dF_dt;
          VertexType *temp = Base::vertices_[0];
          Base::vertices_[0] = Base::vertices_[1];
          Base::vertices_[1] = temp;
        }
        else if (det_dF_dt == 0.0)
          std::cout << "ERROR: detected degenerated element!" << std::endl;
      } * /

      void print(long indent) const
      {
        for (long i = 0; i<indent; ++i)
          std::cout << "   ";
        std::cout << "* dt-dx-Handler: OnAccess" << std::endl;
      }

  };

/ *  //fast and memory-saving: Compute an element's Jacobian as soon as the first entry is accessed. Use static memory, so that only one Jacobian set is stored at the same time. Must not be used in multi-threaded applications!!
  template <typename T_Configuration>
  struct dt_dx_handler<T_Configuration, TetrahedronTag, DtDxStoreStatically>
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
        return dt_dx[3*i + j];
      }

      ScalarType get_det_dF_dt() const
      {
        return det_dF_dt;
      }

      void update_dt_dx()
      {
        checkOrientation();
      }

      void print(long indent = 0) const
      {
        for (long i = 0; i<indent; ++i)
          std::cout << "   ";
        std::cout << "* dt-dx-Handler: StoreStatically" << std::endl;
        Base::print(indent);
      }

      void init_dt_dx() { computeCellJacobian(); }

    protected:

      void computeCellJacobian() const
      {
        PointType & p0 = Base::vertices_[0]->getPoint();
        PointType p1 = Base::vertices_[1]->getPoint() - p0;
        PointType p2 = Base::vertices_[2]->getPoint() - p0;
        PointType p3 = Base::vertices_[3]->getPoint() - p0;

        det_dF_dt = spannedVolume(p1, p2, p3);

        //dt_1/dx_1
        dt_dx[0] = (  + p2[1] * p3[2] - p2[2] * p3[1] ) / det_dF_dt;
        //dt_1/dx_2
        dt_dx[1] = (  - p2[0] * p3[2] + p2[2] * p3[0] ) / det_dF_dt;
        //dt_1/dx_3
        dt_dx[2] = (  + p2[0] * p3[1] - p2[1] * p3[0] ) / det_dF_dt;

        //dt_2/dx_1
        dt_dx[3] = (  - p1[1] * p3[2] + p1[2] * p3[1] ) / det_dF_dt;
        //dt_2/dx_2
        dt_dx[4] = (  + p1[0] * p3[2] - p1[2] * p3[0] ) / det_dF_dt;
        //dt_2/dx_3
        dt_dx[5] = (  - p1[0] * p3[1] + p1[1] * p3[0] ) / det_dF_dt;

        //dt_3/dx_1
        dt_dx[6] = (  + p1[1] * p2[2] - p1[2] * p2[1] ) / det_dF_dt;
        //dt_3/dx_2
        dt_dx[7] = (  - p1[0] * p2[2] + p1[2] * p2[0] ) / det_dF_dt;
        //dt_3/dx_3
        dt_dx[8] = (  + p1[0] * p2[1] - p1[1] * p2[0] ) / det_dF_dt;

      } //computeCellJacobian()

      void checkOrientation()
      {
        PointType & p0 = Base::vertices_[0]->getPoint();
        PointType & p1 = Base::vertices_[1]->getPoint();
        PointType & p2 = Base::vertices_[2]->getPoint();
        PointType & p3 = Base::vertices_[3]->getPoint();

        //volume:
        ScalarType detdFdt = spannedVolume(p1 - p0, p2 - p0, p3 - p0);

        if (detdFdt < 0)
        {
          VertexType *temp = Base::vertices_[0];
          Base::vertices_[0] = Base::vertices_[1];
          Base::vertices_[1] = temp;
        }
        else if (detdFdt == 0.0)
          std::cout << "ERROR: detected degenerated element!" << std::endl;
      }

    private:
      static ScalarType dt_dx[9];                          //inverse of Jacobian matrix from mapping
      static ScalarType det_dF_dt;                         //determinant of Jacobian matrix
  };

  //initialize static variables:
  template <typename T_Configuration>
  typename T_Configuration::CoordType
  dt_dx_handler<T_Configuration, TetrahedronTag, DtDxStoreStatically>::det_dF_dt = 0;

  template <typename T_Configuration>
  typename T_Configuration::CoordType
  dt_dx_handler<T_Configuration, TetrahedronTag, DtDxStoreStatically>::dt_dx[] = {0, 0, 0,
                                   0, 0, 0,
                                   0, 0, 0};

                                   
                                   */
}
#endif