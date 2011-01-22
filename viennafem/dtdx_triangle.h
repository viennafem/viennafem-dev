/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_DTDX_TRIANGLE_HPP
#define VIENNAFEM_DTDX_TRIANGLE_HPP

#include <iostream>
#include "viennagrid/celltags.h"
#include "viennagrid/domain.hpp"
#include "viennafem/forwards.h"


namespace viennafem
{

  //memory-intensive: Compute them once and store the computed values until next update
  template <typename T_Configuration>
  struct dt_dx_handler<T_Configuration, viennagrid::triangle_tag, DtDxStoreAll>
  {
    typedef typename T_Configuration::numeric_type                              ScalarType;
    typedef typename viennagrid::DomainTypes<T_Configuration>::point_type                 PointType;
    typedef typename viennagrid::DomainTypes<T_Configuration>::vertex_type                VertexType;

    public:

      //returns the element dt_i/dx_j of the functional determinant induced by the mapping to the reference element. i and j start at 0.
      ScalarType get_dt_dx(int i, int j) const
      {
        return dt_dx[2*i + j];
      }

      ScalarType get_det_dF_dt() const
      {
        return det_dF_dt;
      }

      template <typename CellType>
      void update_dt_dx(CellType const & cell)
      {
        PointType & p0 = cell.getPoint(0);
        PointType & p1 = cell.getPoint(0);
        PointType & p2 = cell.getPoint(0);

        checkOrientation(p0, p1, p2);

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
      }

      void init_dt_dx() {}

    protected:      
      void checkOrientation(PointType & p0, PointType & p1, PointType & p2)
      {
        //volume:
        det_dF_dt = spannedVolume(p1 - p0, p2 - p0);

        if (det_dF_dt < 0)
        {
          //swap two points:
          //det_dF_dt = - det_dF_dt;
          //VertexType *temp = Base::vertices_[0];
          //Base::vertices_[0] = Base::vertices_[1];
          //Base::vertices_[1] = temp;
        }
        else if (det_dF_dt == 0.0)
          std::cout << "ERROR: detected degenerated element!" << std::endl;
      }

      void print(long indent = 0) const
      {
        for (long i = 0; i<indent; ++i)
          std::cout << "   ";
        std::cout << "* dt-dx-Handler: StoreAll" << std::endl;
      }

    private:
      ScalarType dt_dx[4];                          //inverse of Jacobian matrix from mapping
      ScalarType det_dF_dt;                         //determinant of Jacobian matrix
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