/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_ASSEMBLING_HPP
#define VIENNAFEM_ASSEMBLING_HPP

#include "viennafem/typelist.h"
#include "viennafem/forwards.h"
#include "viennafem/mapping.hpp"

#include "viennagrid/domain.hpp"
#include "viennagrid/iterators.hpp"

namespace viennafem
{
  
  
  /////////////////// OLD CODE TO FOLLOW - IGNORE ////////////////////////
  
  

  struct AssembleMatrixTag{};
  struct AssembleRHSTag{};

  //helper: for topolevel=0, returns 1
  //        for topolevel>0, returns basisfunction-degree
  template <long topolevel, typename BasisFunctionTag>  
  struct BasisFunctionDegreeOnTopoLevel
  {
    //enum{ ReturnValue = TopologyLevel<typename ElementType::ElementTag, 0>::ElementNum };
    enum{ ReturnValue = BasisFunctionTag::degree };
  };

  template <typename BasisFunctionTag>  
  struct BasisFunctionDegreeOnTopoLevel<0, BasisFunctionTag>
  {
    enum{ ReturnValue = 1 };
  };

  //Meta-Class for writing Matrix-Entries:
  template <typename TListIter>
  struct EntryWriter
  {
    //for matrix:
    template <typename MatrixType,
                typename AssemblyCellType,
                typename BasisFunction1,
                typename BasisFunction2,
                typename VectorType>
    static void apply(MatrixType & matrix, long rowindex, long colindex, AssemblyCellType & cell, BasisFunction1 const & bf1, BasisFunction2 const & bf2, VectorType & prev_result1, VectorType & prev_result2)
    {
      typedef typename TListIter::ResultType      ExpressionType;
/*      std::cout << "Writing at " << rowindex << "," << colindex << std::endl;
      bf1.print(); std::cout << ","; bf2.print(); std::cout << std::endl;*/
      matrix(rowindex + TListIter::rowindex, colindex + TListIter::colindex)
          += ExpressionType::evaluate(prev_result1, prev_result2, cell, bf1, bf2);

      EntryWriter<typename TListIter::IncrementType>::apply(matrix, rowindex, colindex, cell, bf1, bf2, prev_result1, prev_result2);
    }

    //for rhs: 
    template <typename VectorType,
                typename AssemblyCellType,
                typename BasisFunction1>
    static void apply(VectorType & rhs, long rowindex, AssemblyCellType & cell, BasisFunction1 const & bf1, VectorType & prev_result1, VectorType & prev_result2)
    {
      typedef typename TListIter::ResultType      ExpressionType;

/*      std::cout << "Adding at rhs(" << rowindex << "+" << TListIter::rowindex << "): ";
      ExpressionType().print(); //  ::evaluate(prev_result1, prev_result2, cell, bf1)
      std::cout << " with basisfun "; bf1.print(); std::cout << std::endl;*/
//      cell.print();
      rhs(rowindex + TListIter::rowindex)
          += ExpressionType::evaluate(prev_result1, prev_result2, cell, bf1);
       //std::cout << "Writing RHS at " << rowindex << ", " << rhs(rowindex + TListIter::rowindex) << std::endl;
//       bf1.print(); std::cout << std::endl;

      EntryWriter<typename TListIter::IncrementType>::apply(rhs, rowindex, cell, bf1, prev_result1, prev_result2);
    }


  };

  //end iteration:
  template <>
  struct EntryWriter<TypeListEnd>
  {
    template <typename MatrixType,
                typename AssemblyCellType,
                typename BasisFunction1,
                typename BasisFunction2,
                typename VectorType>
    static void apply(MatrixType & matrix, long rowindex, long colindex, AssemblyCellType & cell, BasisFunction1 const & bf1, BasisFunction2 const & bf2, VectorType & prev_result1, VectorType & prev_result2)
    {
      //do nothing
    }

    //for rhs:
    template <typename VectorType,
                typename AssemblyCellType,
                typename BasisFunction1>
    static void apply(VectorType & rhs, long rowindex, AssemblyCellType & cell, BasisFunction1 const & bf1, VectorType & prev_result1, VectorType & prev_result2)
    {
      //do nothing
    }

  };

  //Dirichlet boundary data:
  //Part 1: Collect boundary values: (i holds the number of boundary values missing)
  template <typename TListIter, long i = 0, bool allfetched = false>
  struct BoundaryWriter
  {
    template <typename VectorType,
                typename AssemblyCellType,
                typename MapIterator,
                typename BasisFunction1,
                typename BasisFunction2>
    static void apply(  VectorType & rhs,
                        long index,
                        AssemblyCellType & cell,
                        MapIterator & mapit,
                        BasisFunction1 const & bf1,
                        BasisFunction2 const & bf2,
                        VectorType & prev_result1,
                        VectorType & prev_result2  )
    {
      //first fetch: set up double array:
      double bnd[TListIter::colnum];
      bnd[0] = mapit.getBoundaryValue(0);

      BoundaryWriter<TListIter, i+1, TListIter::colnum == (i+1)>::apply(rhs, cell, mapit, index, bf1, bf2, prev_result1, prev_result2, bnd);
    }

    template <typename VectorType,
                typename AssemblyCellType,
                typename MapIterator,
                typename BasisFunction1,
                typename BasisFunction2>
    static void apply(  VectorType & rhs,
                        AssemblyCellType & cell,
                        MapIterator & mapit,
                        long index,
                        BasisFunction1 const & bf1,
                        BasisFunction2 const & bf2,
                        VectorType & prev_result1,
                        VectorType & prev_result2,
                        double* bnd )
    {
      bnd[i] = mapit.getBoundaryValue(i);

      BoundaryWriter<TListIter, i+1, TListIter::colnum == (i+1)>::apply(rhs, cell, mapit, index, bf1, bf2, prev_result1, prev_result2, bnd);
    }
  };

  //Part 2: write to rhs:
  template <typename TListIter, long i>
  struct BoundaryWriter<TListIter, i, true>
  {
    template <typename VectorType,
                typename AssemblyCellType,
                typename MapIterator,
                typename BasisFunction1,
                typename BasisFunction2>
    static void apply(  VectorType & rhs,
                        AssemblyCellType & cell,
                        MapIterator & mapit,
                        long index,
                        BasisFunction1 const & bf1,
                        BasisFunction2 const & bf2,
                        VectorType & prev_result1,
                        VectorType & prev_result2,
                        double* bnd)
    {
      //iterate over typelist:
      typedef typename TListIter::ResultType                  ExpressionType;

      rhs(index + TListIter::rowindex) -= bnd[TListIter::colindex] * ExpressionType::evaluate(prev_result1, prev_result2, cell, bf1, bf2);

      BoundaryWriter<typename TListIter::IncrementType, i, true>::apply(rhs, cell, mapit, index, bf1, bf2, prev_result1, prev_result2, bnd);
    }

  };

  //Part 3: end of recursion:
  template <long i>
  struct BoundaryWriter<TypeListEnd, i, true>
  {
    template <typename VectorType,
                typename AssemblyCellType,
                typename MapIterator,
                typename BasisFunction1,
                typename BasisFunction2>
    static void apply(  VectorType & rhs,
                        AssemblyCellType & cell,
                        MapIterator & mapit,
                        long index,
                        BasisFunction1 const & bf1,
                        BasisFunction2 const & bf2,
                        VectorType & prev_result1,
                        VectorType & prev_result2,
                        double* bnd) {}
  };

  //Meta class for iteration over all basisfunction-pairs and writing the entry:
  template <typename FEMConfig,
              typename AssemblyCellType,
              typename Equation,
              typename BFIterator_v,
              typename BFIterator_u
             >
  struct MatrixIterator
  {
    typedef typename AssemblyCellType::ElementTag                 AssemblyCellTag;

    template <typename MatrixType, typename VectorType, typename MappingIterator_v, typename MappingIterator_u>
    static void assembleRow(MatrixType & matrix, VectorType & rhs, AssemblyCellType & cell, VectorType & prev_result1, VectorType & prev_result2, MappingIterator_v const & mit1, MappingIterator_u & mit2)
    {
      typedef typename FEMConfig::ResultDimension         ResultDimension;
      typedef typename BFIterator_v::ResultType           BFType_v;
      typedef typename BFIterator_u::ResultType           BFType_u;

/*      std::cout << "MatrixIterator at (" << element_id_v << "," << bf_id_v << ") and ("
                                         <<element_id_u << "," << bf_id_u << ")" << std::endl;*/
      //This function will not be called for *mit1 == -1. Only check *mit2:
      if (*mit2 != -1)
      {
        //write:
//         std::cout << "Writing at mit1=" << *mit1 << " and mit2=" << *mit2 << std::endl;
//         std::cout << "Bf u: "; BFType_v().print(); std::cout << std::endl;
//         std::cout << "Bf v: "; BFType_u().print(); std::cout << std::endl;

        EntryWriter < TypeListIterator< typename Equation::LHSType,
                                        ResultDimension::dim,
                                        ResultDimension::dim
                                      >
                      >::apply( matrix,
                                ResultDimension::dim * (*mit1),
                                ResultDimension::dim * (*mit2),
                                cell,
                                BFType_v(),
                                BFType_u(),
                                prev_result1, prev_result2);
      }
      else
      {
        //std::cout << "Writing to boundary..." << std::endl;
        BoundaryWriter  < TypeListIterator< typename Equation::LHSType,
                                        ResultDimension::dim,
                                        ResultDimension::dim
                                      >
                        >::apply( rhs,
                                  ResultDimension::dim * (*mit1),
                                  cell,
                                  mit2,
                                  BFType_v(),
                                  BFType_u(),
                                  prev_result1, prev_result2
                                );
      }

      //proceed in row to next col:
      MatrixIterator< FEMConfig, AssemblyCellType, Equation,
                      BFIterator_v,
                      typename BFIterator_u::NextType
                      >::assembleRow(matrix, rhs, cell, prev_result1, prev_result2, mit1, ++mit2);
    }

    //iterate over basisfunction-tuples: in row-direction
    //Call row-assembly-routine (assembly with increasing column-index) for each iteration
    template <typename MatrixType, typename VectorType, typename MappingIterator_v, typename MappingIterator_u>
    static void assemble(MatrixType & matrix, VectorType & rhs, AssemblyCellType & cell, VectorType & prev_result1, VectorType & prev_result2, MappingIterator_v & mit1, MappingIterator_u & mit2)
    {

//       std::cout << "Basisfunction set 1: "; BFList1().print(); 
//       std::cout << "Basisfunction set 2: "; BFList2().print(); 
      //std::cout << "MatrixIterator at (" << element_id_v << "," << bf_id_v << ")" << std::endl;

      if (*mit1 != -1)
      {
        //std::cout << "Starting row-wise assembly on element with mapit1=" << *mit1 << std::endl;
        MappingIterator_u temp(mit2);
        assembleRow(matrix, rhs, cell, prev_result1, prev_result2, mit1, temp);
      }

      //proceed with next basisfunction on this level (col-index):
      MatrixIterator<FEMConfig, AssemblyCellType, Equation,
                      typename BFIterator_v::NextType,
                      BFIterator_u
                    >::assemble(matrix, rhs, cell, prev_result1, prev_result2, ++mit1, mit2);
    }

    //type erasure iterator:
    template <typename MatrixType, typename VectorType, typename BasisFunction, typename MappingIterator_v, typename MappingIterator_u>
    static void assemble(MatrixType & matrix, VectorType & rhs, AssemblyCellType & cell, VectorType & prev_result1, VectorType & prev_result2, BasisFunction * basisfuns_v, MappingIterator_v & mapit1, BasisFunction * basisfuns_u, MappingIterator_u & mit2)
    {
      typedef typename FEMConfig::ResultDimension                     ResultDimension;

      //Iterate over all basisfunctions and mapping indices simultaneously:
      long i=0;
      long j=0;

      for (; mapit1.valid(); ++i, ++mapit1)
      {
        if (*mapit1 == -1)    //nothing to do because of Dirichlet boundary conditions
          continue;

        j=0;
        MappingIterator_u mapit2(mit2);
        for (; mapit2.valid(); ++j, ++mapit2)
        {
          if (*mapit2 != -1)
          {

            //std::cout << "Accessing " << i << "," << j << std::endl;

            //write to matrix (do not forget about multiple contributions from boundary cells!)
            EntryWriter < TypeListIterator< typename Equation::LHSType,
                                             ResultDimension::dim,
                                             ResultDimension::dim
                                           >
                          >::apply( matrix,
                                  ResultDimension::dim * (*mapit1),
                                  ResultDimension::dim * (*mapit2),
                                  cell,
                                  basisfuns_v[i],
                                  basisfuns_u[j],
                                  prev_result1, prev_result2 );
          }
          else
          {
            //write Dirichlet boundary conditions directly to rhs-vector:
            // (again, do not forget possible multiple contributions from boundary cells!)
            BoundaryWriter  < TypeListIterator< typename Equation::LHSType,
                                                ResultDimension::dim,
                                                ResultDimension::dim
                                              >
                            >::apply( rhs,
                                      ResultDimension::dim * (*mapit1),
                                      cell,
                                      mapit2,
                                      basisfuns_v[i],
                                      basisfuns_u[j],
                                      prev_result1, prev_result2
                                    );
          }
        } //for j

        //write RHS:
        EntryWriter < TypeListIterator< typename Equation::RHSType,
                                        ResultDimension::dim
                                      >
                    >::apply( rhs,
                              ResultDimension::dim * (*mapit1),
                              cell,
                              basisfuns_v[i],
                              prev_result1, prev_result2
                            );

      } //for i

    }

  };

  //end of column:
  template <typename FEMConfig,
              typename AssemblyCellType,
              typename Equation,
              typename BFIterator_u
             >
  struct MatrixIterator<FEMConfig, AssemblyCellType, Equation,
                         IteratorStop, BFIterator_u>
  {
    template <typename MatrixType, typename VectorType, typename MappingIterator_v, typename MappingIterator_u>
    static void assemble(MatrixType & matrix, VectorType & rhs, AssemblyCellType & cell, VectorType & prev_result1, VectorType & prev_result2, MappingIterator_v & mit1, MappingIterator_u & mit2)
    {
      //std::cout << "Matrix-Iterator: Column finished!" << std::endl;
    }
  };

  //for RHS-entries:
  template <typename FEMConfig,
              typename AssemblyCellType,
              typename Equation,
              typename BFIterator_v
             >
  struct MatrixIterator<FEMConfig, AssemblyCellType, Equation,
                         BFIterator_v, IteratorStop>
  {
    //end of row:
    template <typename MatrixType, typename VectorType, typename MappingIterator_v, typename MappingIterator_u>
    static void assembleRow(MatrixType & matrix, VectorType & rhs, AssemblyCellType & cell, VectorType & prev_result1, VectorType & prev_result2, MappingIterator_v const & mit1, MappingIterator_u & mit2)
    {
      //std::cout << "Matrix-Iterator: Row finished!" << std::endl;
      typedef typename FEMConfig::ResultDimension         ResultDimension;
      typedef typename BFIterator_v::ResultType           BFType_v;

      // Right hand side entries:
      EntryWriter < TypeListIterator< typename Equation::RHSType,
                                      ResultDimension::dim
                                    >
                  >::apply( rhs,
                            ResultDimension::dim * (*mit1),
                            cell,
                            BFType_v(),
                            prev_result1, prev_result2
                          );
    }
  };

  //helper for next function:
  template <typename IntDomain>
  struct IntegrationDomainChecker
  {
    //check for element on Gamma<id>
    template <typename AssemblyCellType>
    static bool apply(AssemblyCellType & cell)
    {
      return cell.template hasQuantity<bool>(IntDomain());
    }
  };

  /*template <long id>
  struct IntegrationDomainChecker< Interface<id> >
  {
    //check for element on Interface<id>
    template <typename AssemblyCellType>
    static bool apply(AssemblyCellType & cell)
    {
      return cell.template hasQuantity<AssemblyCellType *>(Interface<id>());
    }
  };*/

  template <>
  struct IntegrationDomainChecker < Omega >
  {
    //for Omega, every element (cell) is valid:
    template <typename AssemblyCellType>
    static bool apply(AssemblyCellType & cell) { return true; }
  };

  // unifies assembleMatrix() and assembleRHS() into a single assembly instance.
  template <typename FEMConfig, typename Segment, typename MatrixType, typename VectorType, typename EquationArray>
  void assemble_impl(Segment & seg, MatrixType & matrix, VectorType & rhs, EquationArray const & eqnarray, VectorType & prev_result1, VectorType & prev_result2, TypeListTag)
  {
    typedef typename FEMConfig::ResultDimension      ResultDimension;
    typedef typename FEMConfig::BaseFunTreatmentTag  BaseFunTreatmentTag;
    typedef typename FEMConfig::BasisfunctionTag     BasisfunctionTag;

    typedef typename Segment::Configuration                      DomainConfiguration;
    typedef typename viennagrid::DomainTypes<DomainConfiguration>::VertexType       VertexType;
    typedef typename viennagrid::DomainTypes<DomainConfiguration>::CellType         Cell;
    typedef typename Cell::ElementTag                                   CellTag;
    typedef typename viennagrid::DomainTypes<DomainConfiguration>::CellIterator     CellIterator;

    // Rearrange equations such that matrix-contributions are on the lhs, rhs-contributions on the rhs: //
    //typedef typename ASSEMBLY_REARRANGER<EquationArray>::ResultType     Equation;
    typedef EquationArray     Equation;   //debug only

    // Set up domain iterator (Omega, Gamma1, ...)
    //typedef typename ASSEMBLY_GET_FIRST_INT_DOMAIN<Equation>::ResultType    IntDomain;
    typedef Omega   IntDomain;

    //Iteration has to be done over all "cells" of the integration domain
    //i.e. for Omega: CellIterator, for Gamma: FacetIterator
    typedef typename viennagrid::IteratorTypes<Segment, 1
                                      //GET_INTEGRATION_TOPO_LEVEL<IntDomain, CellTag>::ReturnValue
                                     >::ResultType                          AssemblyCellIterator;

    //set up the AssemblyCellType:
    typedef typename AssemblyCellIterator::value_type                    AssemblyCellType;
    typedef typename AssemblyCellType::ElementTag                        AssemblyCellTag;

    //basis functions:
    typedef BasisFunctionIterator<AssemblyCellType, BasisfunctionTag>     BFIterator;

    //mapping iterator:
    typedef MappingIterator<FEMConfig, AssemblyCellType, IntDomain, false>       MapIterator_v;
    typedef MappingIterator<FEMConfig, AssemblyCellType, IntDomain, true>       MapIterator_u;

    /************** Fill the Matrix *******************/
    std::cout << "Assembling "; Equation().print(); std::cout << std::endl;

    //Iterate over all AssemblyCells:
    for (AssemblyCellIterator acit = seg.template getLevelIteratorBegin<AssemblyCellTag::TopoLevel>();
        acit != seg.template getLevelIteratorEnd<AssemblyCellTag::TopoLevel>();
        ++acit)
    {
      //iterate over elements in the integration domain only, skip others:
      if ( IntegrationDomainChecker<IntDomain>::apply(*acit) == true)
      {
        MapIterator_v mapit_v(*acit);
        MapIterator_u mapit_u(*acit);

        //acit->print_short();
        acit->init_dt_dx();

        //let MatrixIterator do the job now:
        MatrixIterator<FEMConfig, AssemblyCellType, Equation,
                        BFIterator, BFIterator
                      >::assemble(matrix, rhs, *acit,
                                  prev_result1, prev_result2, mapit_v, mapit_u);
      }
    } //for
  }

  // type erasure implementation:
  template <typename FEMConfig, typename Segment, typename MatrixType, typename VectorType, typename EquationArray>
  void assemble_impl(Segment & seg, MatrixType & matrix, VectorType & rhs, EquationArray const & eqnarray, VectorType & prev_result1, VectorType & prev_result2, TypeErasureTag)
  {
    typedef typename FEMConfig::ResultDimension      ResultDimension;
    typedef typename FEMConfig::BaseFunTreatmentTag  BaseFunTreatmentTag;
    typedef typename FEMConfig::BasisfunctionTag     BasisfunctionTag;

    typedef typename Segment::Configuration                      DomainConfiguration;
    typedef typename viennagrid::DomainTypes<DomainConfiguration>::PointType        PointType;
    typedef typename viennagrid::DomainTypes<DomainConfiguration>::VertexType       VertexType;
    typedef typename viennagrid::DomainTypes<DomainConfiguration>::CellType         Cell;
    typedef typename Cell::ElementTag                                   CellTag;
    typedef typename viennagrid::DomainTypes<DomainConfiguration>::CellIterator     CellIterator;
    //typedef basisfunction<PointType>                                      BasisFunction;
    typedef double BasisFunction;

    // Rearrange equations such that matrix-contributions are on the lhs, rhs-contributions on the rhs: //
    //typedef typename ASSEMBLY_REARRANGER<EquationArray>::ResultType     Equation;
    //typedef EquationArray     Equation;   //debug only
    typedef float Equation;

    // Set up domain iterator (Omega, Gamma1, ...)
    //typedef typename ASSEMBLY_GET_FIRST_INT_DOMAIN<Equation>::ResultType    IntDomain;
    typedef Omega   IntDomain;

    //Iteration has to be done over all "cells" of the integration domain
    //i.e. for Omega: CellIterator, for Gamma: FacetIterator
    typedef typename viennagrid::IteratorTypes<Segment, 1
                                      //GET_INTEGRATION_TOPO_LEVEL<IntDomain, CellTag>::ReturnValue
                                     >::ResultType                          AssemblyCellIterator;

    //set up the AssemblyCellType:
    typedef typename AssemblyCellIterator::value_type                    AssemblyCellType;
    typedef typename AssemblyCellType::ElementTag                        AssemblyCellTag;

    //basis functions:
    typedef BasisFunctionIterator<AssemblyCellType, BasisfunctionTag>     BFIterator;

    //mapping iterator:
    typedef MappingIterator<FEMConfig, AssemblyCellType, IntDomain, false>       MapIterator_v;
    typedef MappingIterator<FEMConfig, AssemblyCellType, IntDomain, true>       MapIterator_u;

    /************** Fill the Matrix *******************/
    //std::cout << "Assembling "; Equation().print(); std::cout << std::endl;

    //set up type erasure basisfunctions here (if needed)
    long bfnum = MappingNumForCell<CellTag, BasisfunctionTag>::ReturnValue;
    BasisFunction basisfuns_v[ bfnum ];
    BFFiller<FEMConfig, AssemblyCellType, BFIterator>::apply(basisfuns_v);

    BasisFunction basisfuns_u[ bfnum ];
    BFFiller<FEMConfig, AssemblyCellType, BFIterator>::apply(basisfuns_u);

    //Iterate over all AssemblyCells:
    for (AssemblyCellIterator acit = seg.template getLevelIteratorBegin<AssemblyCellTag::TopoLevel>();
        acit != seg.template getLevelIteratorEnd<AssemblyCellTag::TopoLevel>();
        ++acit)
    {
      //iterate over elements in the integration domain only, skip others:
      if ( IntegrationDomainChecker<IntDomain>::apply(*acit) == true)
      {
        MapIterator_v mapit_v(*acit);
        MapIterator_u mapit_u(*acit);

        //acit->print();

        acit->init_dt_dx();

        //let MatrixIterator do the job now:
        MatrixIterator<FEMConfig, AssemblyCellType, Equation,
                        BFIterator, BFIterator
                      >::assemble(matrix, rhs, *acit,
                                  prev_result1, prev_result2,
                                  basisfuns_v, mapit_v,
                                  basisfuns_u, mapit_u);
      }
    } //for
  };

  // unifies old assembleMatrix() and assembleRHS() into a single assembly instance.

  // Terminates assembly recursion as soon as all integrals are assembled.
  //template <typename FEMConfig, typename Segment, typename MatrixType, typename VectorType>
  //void assemble(Segment & seg, MatrixType & matrix, VectorType & rhs, EmptyEquation eetype, VectorType & prev_result1, VectorType & prev_result2)
  //{ } 

  template <typename FEMConfig, typename Segment, typename MatrixType, typename VectorType, typename EquationArray>
  void assemble(Segment & seg, MatrixType & matrix, VectorType & rhs, EquationArray const & eqnarray, VectorType & prev_result1, VectorType & prev_result2)
  {
    typedef typename Segment::Configuration                      DomainConfiguration;
    typedef typename viennagrid::DomainTypes<DomainConfiguration>::CellType::ElementTag  CellTag;

    typedef typename FEMConfig::BaseFunTreatmentTag       BaseFunTreatmentTag;
    typedef typename FEMConfig::BasisfunctionTag          BasisFunctionTag;

    seg.template getLevelIteratorBegin<0>()->setCurrentSegment(seg);

    // First step: Collect all contributions from first integration domain and assemble according to the assembly tag:
    //typedef typename ASSEMBLY_GET_FIRST_INT_DOMAIN<EquationArray>::ResultType  IntDomain;
    typedef Omega   IntDomain;
    
    //typedef typename ASSEMBLY_EXTRACT_DOMAIN<EquationArray, IntDomain>::ResultType   CurrentEquation;
    typedef double CurrentEquation;

//     typedef typename TopologyLevel<CellTag,
//                                       GET_INTEGRATION_TOPO_LEVEL<IntDomain, CellTag>::ReturnValue
//                                      >::ElementTag                             AssemblyCellTag;

    assemble_impl<FEMConfig>(seg, matrix, rhs, CurrentEquation(), prev_result1, prev_result2, BaseFunTreatmentTag());

    // Second step: call this function again with already assembled contributions removed:
    //assemble<FEMConfig>(seg, matrix, rhs, typename ASSEMBLY_REMOVE_DOMAIN<EquationArray, IntDomain>::ResultType(), prev_result1, prev_result2);

  } //assemble()

  //some more public interfaces:
  template <typename FEMConfig, typename Segment, typename MatrixType, typename VectorType, typename EquationArray>
  void assemble(Segment & seg, MatrixType & matrix, VectorType & rhs, EquationArray const & eqnarray)
  {
    VectorType dummy(1);
    assemble<FEMConfig>(seg, matrix, rhs, eqnarray, dummy, dummy);
  }

  template <typename FEMConfig, typename Segment, typename MatrixType, typename VectorType, typename EquationArray>
  void assemble(Segment & seg, MatrixType & matrix, VectorType & rhs, EquationArray const & eqnarray, VectorType & prev_result1)
  {
    VectorType dummy(1);
    assemble<FEMConfig>(seg, matrix, rhs, eqnarray, prev_result1, dummy);
  }



  template <typename Matrix, typename Vector, typename VertexBoundaryIterator, typename BoundaryList>
  struct BoundaryConditionsWriter
  {
    static void apply(Matrix & matrix, Vector & rhs, VertexBoundaryIterator & vbit, long startIndex)
    {
      typedef typename BoundaryList::ResultType         BoundaryData;
      BoundaryConditionsWriter<Matrix, Vector, VertexBoundaryIterator,
                               typename BoundaryList::Tail
      >::apply(matrix, rhs, vbit, startIndex,
               (*vbit).template retrieveQuantity<double>( BoundaryData() ));
    }

    static void apply(Matrix & matrix, Vector & rhs, VertexBoundaryIterator & vbit, long startIndex, double bnd1)
    {
      typedef typename BoundaryList::ResultType         BoundaryData;
      BoundaryConditionsWriter<Matrix, Vector, VertexBoundaryIterator,
                               typename BoundaryList::Tail
      >::apply(matrix, rhs, vbit, startIndex, bnd1, (*vbit).template retrieveQuantity<double>( BoundaryData() ) );
    }

    static void apply(Matrix & matrix, Vector & rhs, VertexBoundaryIterator & vbit, long startIndex, double bnd1, double bnd2)
    {
      typedef typename BoundaryList::ResultType         BoundaryData;
      BoundaryConditionsWriter<Matrix, Vector, VertexBoundaryIterator,
                               typename BoundaryList::Tail
      >::apply(matrix, rhs, vbit, startIndex, bnd1, bnd2, (*vbit).template retrieveQuantity<double>( BoundaryData() ) );
    }
  };

  template <typename Matrix, typename Vector, typename VertexBoundaryIterator>
  struct BoundaryConditionsWriter<Matrix, Vector, VertexBoundaryIterator, TypeListEnd>
  {
    static void apply(Matrix & matrix, Vector & rhs, VertexBoundaryIterator & vbit, long startIndex, double bnd1, double bnd2, double bnd3)
    {

      for (long i=0; i<matrix.getRows(); ++i)
      {
        rhs(i) -= bnd1 * matrix(i, startIndex);
        rhs(i) -= bnd2 * matrix(i, startIndex + 1);
        rhs(i) -= bnd3 * matrix(i, startIndex + 2);
      }
      matrix.clearCol(startIndex);
      matrix.clearCol(startIndex + 1);
      matrix.clearCol(startIndex + 2);
  
      matrix.clearRow(startIndex);
      matrix.clearRow(startIndex + 1);
      matrix.clearRow(startIndex + 2);
  
      matrix(startIndex,     startIndex)     = 1.0;
      matrix(startIndex + 1, startIndex + 1) = 1.0;
      matrix(startIndex + 2, startIndex + 2) = 1.0;
      rhs(startIndex)     = bnd1;
      rhs(startIndex + 1) = bnd2;
      rhs(startIndex + 2) = bnd3;
    }

    static void apply(Matrix & matrix, Vector & rhs, VertexBoundaryIterator & vbit, long startIndex, double bnd1, double bnd2)
    {
      for (long i=0; i<matrix.getRows(); ++i)
      {
        rhs(i) -= bnd1 * matrix(i, startIndex);
        rhs(i) -= bnd2 * matrix(i, startIndex + 1);
      }
      matrix.clearCol(startIndex);
      matrix.clearCol(startIndex + 1);

      matrix.clearRow(startIndex);
      matrix.clearRow(startIndex + 1);

      matrix(startIndex,     startIndex)     = 1.0;
      matrix(startIndex + 1, startIndex + 1) = 1.0;
      rhs(startIndex)     = bnd1;
      rhs(startIndex + 1) = bnd2;
    }

    static void apply(Matrix & matrix, Vector & rhs, VertexBoundaryIterator & vbit, long startIndex, double bnd1)
    {
      //std::cout << "Applying boundary values with bnd1=" << bnd1 << std::endl;
      for (long i=0; i<matrix.getRows(); ++i)
      {
        rhs(i) -= bnd1 * matrix(i, startIndex);
        //matrix(i, startIndex) = 0.0;
      }

      matrix.clearCol(startIndex);
      matrix.clearRow(startIndex);
  
      matrix(startIndex,     startIndex)     = 1.0;
      rhs(startIndex)     = bnd1;
    }

  };
  
  //creates the typelist for boundary conditions in dependence of the problem dimension:
  template <typename FEMConfig, long dim>
  struct BoundaryDataListFactory
  { };

  template <typename FEMConfig>
  struct BoundaryDataListFactory<FEMConfig, 3>
  {
    typedef TypeList< typename FEMConfig::BoundaryData,
                        TypeList< typename FEMConfig::BoundaryData2,
                                  TypeList< typename FEMConfig::BoundaryData3,
                                            TypeListEnd >
                                >
                      >                                 ResultType;
  };

  template <typename FEMConfig>
  struct BoundaryDataListFactory<FEMConfig, 2>
  {
    typedef TypeList< typename FEMConfig::BoundaryData,
                        TypeList< typename FEMConfig::BoundaryData2,
                                  TypeListEnd
                                >
                      >                                 ResultType;
  };

  template <typename FEMConfig>
  struct BoundaryDataListFactory<FEMConfig, 1>
  {
    typedef TypeList< typename FEMConfig::BoundaryData, TypeListEnd>      ResultType;
  };


  template <typename FEMConfig, typename Segment, typename Matrix, typename Vector>
  long applyBoundaryConditions_impl(Segment & seg, Matrix & matrix, Vector & rhs, NoBoundaryMappingTag)
  {
    return matrix.getRowNum();
  }

  template <typename FEMConfig, typename Segment, typename Matrix, typename Vector>
  long applyBoundaryConditions_impl(Segment & seg, Matrix & matrix, Vector & rhs, FullMappingTag)
  {
    typedef typename FEMConfig::ResultDimension      ResultDimension;
    typedef typename FEMConfig::MappingKey           MappingKey;
    typedef typename FEMConfig::BoundaryKey          BoundaryKey;
    typedef typename BoundaryDataListFactory<FEMConfig, ResultDimension::dim>::ResultType    BoundaryList;

    typedef typename Segment::Configuration                      DomainConfiguration;
    typedef typename viennagrid::DomainTypes<DomainConfiguration>::PointType            PointType;
    typedef typename viennagrid::DomainTypes<DomainConfiguration>::VertexIterator       VertexIterator;

    //Iterate over boundary:
    std::cout << "Boundary: " << std::endl;

    long dim = matrix.getRows();

    for (VertexIterator vit = seg.getVertexBegin();
          vit != seg.getVertexEnd();
          ++vit)
    {
      //vbit->print(); std::cout << std::endl;
      if (vit->template hasQuantity<bool>(BoundaryKey()) == false)
        continue;

      long vertexID = vit->getID();
      
      BoundaryConditionsWriter< Matrix,
                                Vector,
                                VertexIterator,
                                BoundaryList
                              >::apply(matrix, rhs, vit, ResultDimension::dim * vertexID); 

      dim -= ResultDimension::dim;
    }

    return dim;
  }

  //interface:
  template <typename FEMConfig, typename Segment, typename Matrix, typename Vector>
  long applyBoundaryConditions(Segment & seg, Matrix & matrix, Vector & rhs)
  {
    typedef typename FEMConfig::MappingTag             MappingTag;
    return applyBoundaryConditions_impl<FEMConfig>(seg, matrix, rhs, MappingTag());
  }


}
#endif
