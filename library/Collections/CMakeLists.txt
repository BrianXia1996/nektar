SET(COLLECTIONS_SOURCES
  CoalescedGeomData.cpp
  Collection.cpp
  CollectionOptimisation.cpp
  Operator.cpp
  BwdTrans.cpp
  Helmholtz.cpp
  IProductWRTBase.cpp
  IProductWRTDerivBase.cpp
  IProduct.cpp
  PhysDeriv.cpp
  )

SET(COLLECTIONS_HEADERS
  CoalescedGeomData.h
  Collection.h
  CollectionOptimisation.h
  IProduct.h
  MatrixFreeBase.h
  Operator.h      
)

ADD_NEKTAR_LIBRARY(Collections
    SOURCES ${COLLECTIONS_SOURCES}
    HEADERS ${COLLECTIONS_HEADERS}
    DEPENDS LocalRegions MatrixFreeOps
    SUMMARY "Nektar++ Collections library"
    DESCRIPTION "This library provides efficient multi-expansion implementations of operators on the various supported regions.")


INSTALL(DIRECTORY ./
        DESTINATION ${NEKTAR_INCLUDE_DIR}/Collections
        COMPONENT dev
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
