#include "sps/default.h"
#include "sps/index.h"
#include <type_traits>

namespace cm
{

template <bool CACHED> class SpsInterface;

template <template <size_t D, bool dependant_dim, bool uniform_overlay_grid, size_t orthope> typename storage_t,
          size_t D, size_t O, bool CACHED>
sps::Index<storage_t<D, false, true, O>>& getIndexHelper( SpsInterface<CACHED>* );

template <bool CACHED> class SpsInterface
{
  private:
    template <size_t D, bool dependant_dim, bool uniform_overlay_grid, size_t orthope>
    using storage_t = typename std::conditional<CACHED, CachedTypeDef<D, dependant_dim, uniform_overlay_grid, orthope>,
                                                DiskTypeDef<D, dependant_dim, uniform_overlay_grid, orthope>>::type;

    using I10_t = sps::Index<storage_t<1, false, true, 0>>;
    using I21_t = sps::Index<storage_t<2, false, true, 1>>;

    using I20_t = sps::Index<storage_t<2, false, true, 0>>;
    using I31_t = sps::Index<storage_t<3, false, true, 1>>;
    using I42_t = sps::Index<storage_t<4, false, true, 2>>;

    using I30_t = sps::Index<storage_t<3, false, true, 0>>;
    using I52_t = sps::Index<storage_t<5, false, true, 2>>;

  public:
    I10_t xIndex10;
    I21_t xIndex21;
    I20_t xIndex20;
    I31_t xIndex31;
    I42_t xIndex42;
    I30_t xIndex30;
    I52_t xIndex52;


    SpsInterface( std::string sFilePrefix )
        : xIndex10( sFilePrefix + ".1.0" ),
          xIndex21( sFilePrefix + ".2.1" ),
          xIndex20( sFilePrefix + ".2.0" ),
          xIndex31( sFilePrefix + ".3.1" ),
          xIndex42( sFilePrefix + ".4.2" ),
          xIndex30( sFilePrefix + ".3.0" ),
          xIndex52( sFilePrefix + ".5.2" )
    {}

    template <size_t D, size_t O> sps::Index<storage_t<D + O, false, true, O>>& getIndex( )
    {
        return getIndexHelper<storage_t, D + O, O>( this );
    }
};


#define IMPLEMENT_GET_INDEX_HELPER( D, O )                                                                             \
    template <>                                                                                                        \
    sps::Index<CachedTypeDef<D, false, true, O>>& getIndexHelper<CachedTypeDef, D, O, true>( SpsInterface<true> *      \
                                                                                             pInterface )              \
    {                                                                                                                  \
        return pInterface->xIndex##D##O;                                                                               \
    }                                                                                                                  \
    template <>                                                                                                        \
    sps::Index<DiskTypeDef<D, false, true, O>>& getIndexHelper<DiskTypeDef, D, O, false>( SpsInterface<false> *        \
                                                                                          pInterface )                 \
    {                                                                                                                  \
        return pInterface->xIndex##D##O;                                                                               \
    }

IMPLEMENT_GET_INDEX_HELPER( 1, 0 )
IMPLEMENT_GET_INDEX_HELPER( 2, 1 )

IMPLEMENT_GET_INDEX_HELPER( 2, 0 )
IMPLEMENT_GET_INDEX_HELPER( 3, 1 )
IMPLEMENT_GET_INDEX_HELPER( 4, 2 )

IMPLEMENT_GET_INDEX_HELPER( 3, 0 )
IMPLEMENT_GET_INDEX_HELPER( 5, 2 )

} // namespace cm