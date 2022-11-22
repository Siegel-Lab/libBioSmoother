#include "cm/annotation_index.h"
#include "sps/default.h"
#include "sps/index.h"
#include <type_traits>

#pragma once

namespace cm
{

template <bool CACHED> class SpsInterface;

template <template <size_t D, size_t orthope> typename storage_t, size_t D, size_t O, bool CACHED>
std::shared_ptr<sps::Index<storage_t<D, O>>> getIndexHelper( SpsInterface<CACHED>* );


#define RELEVANT_COMBINATIONS( MACRO )                                                                                 \
    MACRO( 1, 0 )                                                                                                      \
    MACRO( 2, 1 )                                                                                                      \
                                                                                                                       \
    MACRO( 2, 0 )                                                                                                      \
    MACRO( 3, 1 )                                                                                                      \
    MACRO( 4, 2 )                                                                                                      \
                                                                                                                       \
    MACRO( 3, 0 )                                                                                                      \
    MACRO( 5, 2 )


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

#define INSERT_COUNT( D, O )                                                                                           \
    case COMBINE( D, O ):                                                                                              \
        if( vStart.size( ) != D - O )                                                                                  \
            throw std::logic_error( "start position of count must have dimensionality equal to uiD - uiO" );           \
        if( vEnd.size( ) != D - O )                                                                                    \
            throw std::logic_error( "end position of count must have dimensionality equal to uiD - uiO" );             \
                                                                                                                       \
        std::array<uint64_t, D - O> aStart##D##O;                                                                      \
        std::copy_n( vStart.begin( ), D - O, aStart##D##O.begin( ) );                                                  \
                                                                                                                       \
        if constexpr( O == 0 )                                                                                         \
            pIndex##D##O->addPoint( aStart##D##O, 1 );                                                                 \
        else                                                                                                           \
        {                                                                                                              \
            std::array<uint64_t, D - O> aEnd##D##O;                                                                    \
            std::copy_n( vEnd.begin( ), D - O, aEnd##D##O.begin( ) );                                                  \
                                                                                                                       \
            pIndex##D##O->addPoint( aStart##D##O, aEnd##D##O, 1 );                                                     \
        }                                                                                                              \
        break;

#define CLEAR_P_D( D, O ) pIndex##D##O->clearCorners( );

template <bool CACHED> class SpsInterface
{
  private:
    template <size_t D, size_t orthope>
    using storage_t = typename std::conditional<CACHED, CachedTypeDef<D, orthope>, DiskTypeDef<D, orthope>>::type;

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
    std::shared_ptr<sps::Index<CachedTypeDef<D, O>>> getIndexHelper<CachedTypeDef, D, O, true>( SpsInterface<true> *   \
                                                                                                pInterface )           \
    {                                                                                                                  \
        return pInterface->pIndex##D##O;                                                                               \
    }                                                                                                                  \
    template <>                                                                                                        \
    std::shared_ptr<sps::Index<DiskTypeDef<D, O>>> getIndexHelper<DiskTypeDef, D, O, false>( SpsInterface<false> *     \
                                                                                             pInterface )              \
    {                                                                                                                  \
        return pInterface->pIndex##D##O;                                                                               \
    }

RELEVANT_COMBINATIONS( IMPLEMENT_GET_INDEX_HELPER )

} // namespace cm