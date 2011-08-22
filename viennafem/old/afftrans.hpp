/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

/******
*        AFFINE TRANSFORMATION MODULE
*
*   Provides tranformations from the reference element to the element in space
*
*
************/

#ifndef VIENNAFEM_AFFTRANS_GUARD
#define VIENNAFEM_AFFTRANS_GUARD

//#include "typelist.h"

#include "viennagrid/celltags.hpp"
#include "viennagrid/domain.hpp"

namespace viennafem
{

  /******************* mapToElement: *****************************/
  //Transforms a point 'refpoint' given within the reference cell to the corresponding point in the cell.

  //first the overloads:
  //namespace {

    //cell is a triangle:
    template <typename CellType>
    typename viennagrid::result_of::point<typename CellType::config_type>::type 
    mapToElement_impl(CellType const & cell,
                      typename viennagrid::result_of::point<typename CellType::config_type>::type const & refpoint,
                      viennagrid::line_tag)
    {
      //typename viennagrid::DomainTypes<typename CellType::Configuration>::PointType & point1 = cell.getPoint(0);
      //typename viennagrid::DomainTypes<typename CellType::Configuration>::PointType & point2 = cell.getPoint(1);
      //double scalar = 2.0;
      //return point1 + point2;
      return cell.getPoint(0) * (1.0 - refpoint.get_x())
             + cell.getPoint(1) * refpoint.get_x();
    }

    //cell is a triangle:
    template <typename CellType>
    typename viennagrid::DomainTypes<typename CellType::Configuration>::point_type
    mapToElement_impl(CellType const & cell,
                      typename viennagrid::DomainTypes<typename CellType::config_type>::point_type const & refpoint,
                      viennagrid::triangle_tag)
    {
      //typename viennagrid::DomainTypes<typename CellType::Configuration>::PointType & point1 = cell.getPoint(0);
      //typename viennagrid::DomainTypes<typename CellType::Configuration>::PointType & point2 = cell.getPoint(1);
      //double scalar = 2.0;
      //return point1 + point2;
      return cell.getPoint(0) * (1.0 - refpoint.get_x() - refpoint.get_y())
             + cell.getPoint(1) * refpoint.get_x()
             + cell.getPoint(2) * refpoint.get_y();
    }

    /*
    template <typename CellType, typename T, typename U>
    typename viennagrid::DomainTypes<typename CellType::Configuration>::PointType
    mapToElement_impl(CellType const & cell,
                      TypeList<T, U> const & refpoint,
                      TriangleTag)
    {
      //typename viennagrid::DomainTypes<typename CellType::Configuration>::PointType & point1 = cell.getPoint(0);
      //typename viennagrid::DomainTypes<typename CellType::Configuration>::PointType & point2 = cell.getPoint(1);
      //double scalar = 2.0;
      //return point1 + point2;
      return cell.getPoint(0) * (1.0 - typename ElementAt<TypeList<T,U>, 0>::ResultType()() 
                                      - typename ElementAt<TypeList<T,U>, 1>::ResultType()())
             + cell.getPoint(1) * typename ElementAt<TypeList<T,U>, 0>::ResultType()()
             + cell.getPoint(2) * typename ElementAt<TypeList<T,U>, 1>::ResultType()();
    }*/

    //cell is a tetrahedron:
    template <typename CellType>
    typename viennagrid::DomainTypes<typename CellType::config_type>::point_type
    mapToElement_impl(CellType const & cell,
                      typename viennagrid::DomainTypes<typename CellType::config_type>::point_type const & refpoint,
                      viennagrid::tetrahedron_tag)
    {
      return cell.getPoint(0) * (1.0 - refpoint.get_x() - refpoint.get_y() - refpoint.get_z())
            + cell.getPoint(1) * refpoint.get_x()
            + cell.getPoint(2) * refpoint.get_y()
            + cell.getPoint(3) * refpoint.get_z();
    }

    /*
    template <typename CellType, typename T, typename U>
    typename viennagrid::DomainTypes<typename CellType::Configuration>::PointType
    mapToElement_impl(CellType const & cell,
                      TypeList<T, U> const & refpoint,
                      TetrahedronTag)
    {
      return cell.getPoint(0) * (1.0 - typename ElementAt<TypeList<T,U>, 0>::ResultType()() 
                                      - typename ElementAt<TypeList<T,U>, 1>::ResultType()() 
                                      - typename ElementAt<TypeList<T,U>, 2>::ResultType()())
            + cell.getPoint(1) * typename ElementAt<TypeList<T,U>, 0>::ResultType()()
            + cell.getPoint(2) * typename ElementAt<TypeList<T,U>, 1>::ResultType()()
            + cell.getPoint(3) * typename ElementAt<TypeList<T,U>, 2>::ResultType()();
    } */

  //}

  //the interface:
  template <typename CellType>
  typename viennagrid::DomainTypes<typename CellType::config_type>::point_type
  mapToElement( CellType const & cell,
                typename viennagrid::DomainTypes<typename CellType::config_type>::point_type const & refpoint)
  {
    //overload with respect to the cell_tag:
    return mapToElement_impl<CellType>(cell, refpoint, typename CellType::element_tag());
  }

 
  //template <typename CellType, typename T, typename U>
  //typename viennagrid::DomainTypes<typename CellType::Configuration>::PointType mapToElement(CellType const & cell, TypeList<T,U> const & refpoint)
  //{
  //  //overload with respect to the cell_tag:
  //  return mapToElement_impl<CellType>(cell, refpoint, typename CellType::element_tag());
  //}


  //a functor:
  template <typename CellType>
  class MappingFunctor
  {
    public:
      MappingFunctor(CellType & cell) : cell_(cell) {};

      typename viennagrid::DomainTypes<typename CellType::Configuration>::PointType operator()(typename viennagrid::DomainTypes<typename CellType::Configuration>::PointType const & refpoint)
      {
        return mapToElement(cell_, refpoint);
      }

      CellType & getCell() const { return cell_; }

    private:
      CellType & cell_;
  };

}
#endif
