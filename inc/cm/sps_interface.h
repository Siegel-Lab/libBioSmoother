#include "cm/annotation_index.h"
#include "sps/default.h"
#include "sps/index.h"
#include <type_traits>

#pragma once

namespace cm
{

static const bool BIN_SEARCH_SPARSE = false;
#define MAX_ANNO_FILTERS 3

template <bool CACHED> class SpsInterface;

template <template <size_t, size_t, bool> typename storage_t, size_t D, size_t O, bool CACHED>
std::shared_ptr<sps::Index<storage_t<D, O, BIN_SEARCH_SPARSE>>> getIndexHelper( SpsInterface<CACHED>* );

// @todo some of these are not actually needed
#define RELEVANT_COMBINATIONS( MACRO )                                                                                 \
    MACRO( 1, 0 )                                                                                                      \
    MACRO( 2, 1 )                                                                                                      \
                                                                                                                       \
    MACRO( 2, 0 )                                                                                                      \
    MACRO( 4, 0 )                                                                                                      \
    MACRO( 4, 1 )                                                                                                      \
    MACRO( 3, 1 )                                                                                                      \
    MACRO( 4, 2 )                                                                                                      \
    MACRO( 5, 1 )                                                                                                      \
                                                                                                                       \
    MACRO( 3, 0 )                                                                                                      \
    MACRO( 5, 0 )                                                                                                      \
    MACRO( 5, 2 )                                                                                                      \
    MACRO( 6, 2 )                                                                                                      \
    MACRO( 7, 2 )


#define DECLARE_INDEX_TYPE( D, O ) using I##D##O##_t = sps::Index<storage_t<D, O>>;

#define DECLARE_INDEX_OBJECT( D, O ) std::shared_ptr<I##D##O##_t> pIndex##D##O = nullptr;


#define INDEX_LOADED( D, O )                                                                                           \
    if( pIndex##D##O == nullptr )                                                                                      \
        return false;

#define INIT_INDEX_OBJECT( D, O )                                                                                      \
    this->pIndex##D##O =                                                                                               \
        std::make_shared<I##D##O##_t>( sFilePrefix + "/" + std::to_string( D ) + "." + std::to_string( O ), bWrite );


#define D_INC 8

#define COMBINE( D, O ) ( D << D_INC ) | O

#define GENERATE( D, O )                                                                                               \
    case COMBINE( D, O ):                                                                                              \
        uiId = pIndex##D##O->generate( fFac, uiVerbosity );                                                            \
        break;

#define COUNT( D, O )                                                                                                  \
    case COMBINE( D, O ): {                                                                                            \
        std::array<size_t, D - O> vF##D##O;                                                                            \
        std::array<size_t, D - O> vT##D##O;                                                                            \
        size_t uiI##D##O = 0;                                                                                          \
                                                                                                                       \
        for( bool bHasAnno : vHasAnno )                                                                                \
            if( bHasAnno )                                                                                             \
            {                                                                                                          \
                vF##D##O[ uiI##D##O ] = vFromAnnoFilter[ uiI##D##O ];                                                  \
                vT##D##O[ uiI##D##O ] = vToAnnoFilter[ uiI##D##O ];                                                    \
                uiI##D##O++;                                                                                           \
            }                                                                                                          \
        vF##D##O[ uiI##D##O ] = vFrom[ 0 ];                                                                            \
        vT##D##O[ uiI##D##O ] = vTo[ 0 ];                                                                              \
        uiI##D##O++;                                                                                                   \
        if constexpr( b2DPoints )                                                                                      \
        {                                                                                                              \
            vF##D##O[ uiI##D##O ] = vFrom[ 1 ];                                                                        \
            vT##D##O[ uiI##D##O ] = vTo[ 1 ];                                                                          \
            uiI##D##O++;                                                                                               \
        }                                                                                                              \
        if( bHasMapQ )                                                                                                 \
        {                                                                                                              \
            vF##D##O[ uiI##D##O ] = uiMapQMax;                                                                         \
            vT##D##O[ uiI##D##O ] = uiMapQMin;                                                                         \
        }                                                                                                              \
        uiVal = pIndex##D##O->count( iDataSetId, vF##D##O, vT##D##O, xIntersect, uiVerbosity );                        \
    }                                                                                                                  \
    break;

#define GRID_COUNT( D, O )                                                                                             \
    case COMBINE( D, O ): {                                                                                            \
        std::array<size_t, D - O> vF##D##O;                                                                            \
        std::array<size_t, D - O> vS##D##O;                                                                            \
        std::array<size_t, D - O> vN##D##O;                                                                            \
        size_t uiI##D##O = 0;                                                                                          \
                                                                                                                       \
        for( bool bHasAnno : vHasAnno )                                                                                \
            if( bHasAnno && false )                                                                                    \
            {                                                                                                          \
                vF##D##O[ uiI##D##O ] = vFromAnnoFilter[ uiI##D##O ];                                                  \
                vS##D##O[ uiI##D##O ] = 1+vToAnnoFilter[ uiI##D##O ] - vFromAnnoFilter[ uiI##D##O ];                     \
                vN##D##O[ uiI##D##O ] = 1;                                                                             \
                uiI##D##O++;                                                                                           \
            }                                                                                                          \
        vF##D##O[ uiI##D##O ] = vFrom[ 0 ];                                                                            \
        vS##D##O[ uiI##D##O ] = vSize[ 0 ];                                                                            \
        vN##D##O[ uiI##D##O ] = vNum[ 0 ];                                                                             \
        uiI##D##O++;                                                                                                   \
        if constexpr( b2DPoints )                                                                                      \
        {                                                                                                              \
            vF##D##O[ uiI##D##O ] = vFrom[ 1 ];                                                                        \
            vS##D##O[ uiI##D##O ] = vSize[ 1 ];                                                                        \
            vN##D##O[ uiI##D##O ] = vNum[ 1 ];                                                                         \
            uiI##D##O++;                                                                                               \
        }                                                                                                              \
        if( bHasMapQ )                                                                                                 \
        {                                                                                                              \
            vF##D##O[ uiI##D##O ] = uiMapQMax;                                                                         \
            vS##D##O[ uiI##D##O ] = 1+uiMapQMin - uiMapQMax;                                                           \
            vN##D##O[ uiI##D##O ] = 1;                                                                                 \
        }                                                                                                              \
        std::cout << "vF##D##O " << vF##D##O << "vS##D##O " << vS##D##O << "vN##D##O " << vN##D##O << std::endl; \
        xRet = pIndex##D##O->gridCount( iDataSetId, vF##D##O, vS##D##O, vN##D##O, xIntersect, uiVerbosity );       \
    }                                                                                                                  \
    break;

#define INSERT_COUNT( D, O )                                                                                           \
    case COMBINE( D, O ):                                                                                              \
        if( vStart.size( ) != D - O )                                                                                  \
            throw std::logic_error( "start position of count must have dimensionality equal to uiD - uiO" );           \
        if( vEnd.size( ) != D - O )                                                                                    \
            throw std::logic_error( "end position of count must have dimensionality equal to uiD - uiO" );             \
                                                                                                                       \
        if constexpr( O == 0 )                                                                                         \
        {                                                                                                              \
            std::array<uint64_t, D - O> aStart##D##O;                                                                  \
            std::copy_n( vStart.begin( ), D - O, aStart##D##O.begin( ) );                                              \
            pIndex##D##O->addPoint( aStart##D##O );                                                                    \
        }                                                                                                              \
        else                                                                                                           \
        {                                                                                                              \
            std::array<uint64_t, D - O> aStart##D##O;                                                                  \
            std::copy_n( vStart.begin( ), D - O, aStart##D##O.begin( ) );                                              \
            std::array<uint64_t, D - O> aEnd##D##O;                                                                    \
            std::copy_n( vEnd.begin( ), D - O, aEnd##D##O.begin( ) );                                                  \
            pIndex##D##O->addPoint( aStart##D##O, aEnd##D##O );                                                        \
        }                                                                                                              \
        break;

#define CLEAR_P_D( D, O ) pIndex##D##O->clearCorners( );

template <bool CACHED> class SpsInterface
{
  private:
    template <size_t D, size_t orthope>
    using storage_t = typename std::conditional<CACHED, CachedTypeDef<D, orthope, BIN_SEARCH_SPARSE>,
                                                DiskTypeDef<D, orthope, BIN_SEARCH_SPARSE>>::type;

    RELEVANT_COMBINATIONS( DECLARE_INDEX_TYPE )

  public:
    RELEVANT_COMBINATIONS( DECLARE_INDEX_OBJECT )


    using anno_t = AnnotationDescIndex<DiskVecGenerator>;

    anno_t vAnno;


    SpsInterface( std::string sFilePrefix, bool bWrite )
        : vAnno( sFilePrefix + "/anno", bWrite ) //
    { //
        RELEVANT_COMBINATIONS( INIT_INDEX_OBJECT ) //
    }

    SpsInterface( )
    {}

    bool loaded( )
    {
        RELEVANT_COMBINATIONS( INDEX_LOADED )

        return true;
    }

    template <size_t D, size_t O> std::shared_ptr<sps::Index<storage_t<D + O, O>>> getIndex( )
    {
        if constexpr( CACHED )
            return getIndexHelper<CachedTypeDef, D + O, O, true>( this );
        else
            return getIndexHelper<DiskTypeDef, D + O, O, false>( this );
    }

  private:
    constexpr uint16_t combine( uint16_t uiD, uint16_t uiO )
    {
        if( uiO >= 1 << D_INC )
            throw std::logic_error( "number of orthotope dimensions too high" );
        if( uiD >= 1 << D_INC )
            throw std::logic_error( "number of dimensions too high" );
        return COMBINE( uiD, uiO );
    }

  public:
    void insert( size_t uiD, size_t uiO, std::vector<uint64_t> vStart, std::vector<uint64_t> vEnd )
    {
        switch( combine( uiD, uiO ) )
        {
            RELEVANT_COMBINATIONS( INSERT_COUNT )

            default:
                throw std::logic_error( "invalid combination of dimensions and orthotope dimensions" );
                break;
        }
    }

    std::vector<unsigned int> gridCount( size_t uiD, size_t uiO, size_t iDataSetId, std::array<size_t, 2> vFrom,
                                   std::array<size_t, 2> vSize, std::array<size_t, 2> vNum, bool bHasMapQ,
                                   size_t uiMapQMin, size_t uiMapQMax, std::array<bool, MAX_ANNO_FILTERS> vHasAnno,
                                   std::array<size_t, MAX_ANNO_FILTERS>& vFromAnnoFilter,
                                   std::array<size_t, MAX_ANNO_FILTERS>& vToAnnoFilter,
                                   sps::IntersectionType xIntersect, size_t uiVerbosity = 1 )
    {
        constexpr bool b2DPoints = true;
        std::vector<unsigned int> xRet;
        switch( combine( uiD, uiO ) )
        {
            RELEVANT_COMBINATIONS( GRID_COUNT )

            default:
                throw std::logic_error( "invalid combination of dimensions and orthotope dimensions" );
                break;
        }
        return xRet;
    }

    size_t count( size_t uiD, size_t uiO, size_t iDataSetId, std::array<size_t, 2> vFrom, std::array<size_t, 2> vTo,
                  bool bHasMapQ, size_t uiMapQMin, size_t uiMapQMax, std::array<bool, MAX_ANNO_FILTERS> vHasAnno,
                  std::array<size_t, MAX_ANNO_FILTERS>& vFromAnnoFilter,
                  std::array<size_t, MAX_ANNO_FILTERS>& vToAnnoFilter, sps::IntersectionType xIntersect,
                  size_t uiVerbosity = 1 )
    {
        constexpr bool b2DPoints = true;
        size_t uiVal;
        switch( combine( uiD, uiO ) )
        {
            RELEVANT_COMBINATIONS( COUNT )

            default:
                throw std::logic_error( "invalid combination of dimensions and orthotope dimensions" );
                break;
        }
        return uiVal;
    }

    size_t count( size_t uiD, size_t uiO, size_t iDataSetId, std::array<size_t, 1> vFrom, std::array<size_t, 1> vTo,
                  bool bHasMapQ, size_t uiMapQMin, size_t uiMapQMax, std::array<bool, MAX_ANNO_FILTERS> vHasAnno,
                  std::array<size_t, MAX_ANNO_FILTERS>& vFromAnnoFilter,
                  std::array<size_t, MAX_ANNO_FILTERS>& vToAnnoFilter, sps::IntersectionType xIntersect,
                  size_t uiVerbosity = 1 )
    {
        constexpr bool b2DPoints = false;
        size_t uiVal;
        switch( combine( uiD, uiO ) )
        {
            RELEVANT_COMBINATIONS( COUNT )

            default:
                throw std::logic_error( "invalid combination of dimensions and orthotope dimensions" );
                break;
        }
        return uiVal;
    }

    size_t generate( size_t uiD, size_t uiO, double fFac = -1, size_t uiVerbosity = 1 )
    {
        size_t uiId;
        switch( combine( uiD, uiO ) )
        {
            RELEVANT_COMBINATIONS( GENERATE )

            default:
                throw std::logic_error( "invalid combination of dimensions and orthotope dimensions" );
                break;
        }
        return uiId;
    }

    void clearPointsAndDesc( )
    {
        RELEVANT_COMBINATIONS( CLEAR_P_D )
    }
};


#define IMPLEMENT_GET_INDEX_HELPER( D, O )                                                                             \
    template <>                                                                                                        \
    std::shared_ptr<sps::Index<CachedTypeDef<D, O, BIN_SEARCH_SPARSE>>> getIndexHelper<CachedTypeDef, D, O, true>(     \
        SpsInterface<true> * pInterface )                                                                              \
    {                                                                                                                  \
        return pInterface->pIndex##D##O;                                                                               \
    }                                                                                                                  \
    template <>                                                                                                        \
    std::shared_ptr<sps::Index<DiskTypeDef<D, O, BIN_SEARCH_SPARSE>>> getIndexHelper<DiskTypeDef, D, O, false>(        \
        SpsInterface<false> * pInterface )                                                                             \
    {                                                                                                                  \
        return pInterface->pIndex##D##O;                                                                               \
    }

RELEVANT_COMBINATIONS( IMPLEMENT_GET_INDEX_HELPER )

} // namespace cm