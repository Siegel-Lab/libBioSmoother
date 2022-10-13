#include "sps/default.h"
#include "sps/index.h"
#include <type_traits>

namespace cm
{
template <bool CACHED> class SpsInterface
{
  private:
    template <size_t D, bool dependant_dim, bool uniform_overlay_grid, size_t orthope>
    using storage_t = typename std::conditional<CACHED,
                                                CachedTypeDef<D, dependant_dim, uniform_overlay_grid, orthope>,
                                                DiskTypeDef<D, dependant_dim, uniform_overlay_grid, orthope>>::type;

    using I10_t = sps::Index<storage_t<1, false, true, 0>>;
    using I21_t = sps::Index<storage_t<2, false, true, 1>>;

    using I20_t = sps::Index<storage_t<2, false, true, 0>>;
    using I31_t = sps::Index<storage_t<3, false, true, 1>>;
    using I42_t = sps::Index<storage_t<4, false, true, 2>>;

    using I30_t = sps::Index<storage_t<3, false, true, 0>>;
    using I52_t = sps::Index<storage_t<5, false, true, 2>>;

    I10_t xIndex10;
    I21_t xIndex21;
    I20_t xIndex20;
    I31_t xIndex31;
    I42_t xIndex42;
    I30_t xIndex30;
    I52_t xIndex52;


  public:
    SpsInterface( std::string sFilePrefix )
        : xIndex10( sFilePrefix + ".1.0" ),
            xIndex21( sFilePrefix + ".2.1" ),
            xIndex20( sFilePrefix + ".2.0" ),
            xIndex31( sFilePrefix + ".3.1" ),
            xIndex42( sFilePrefix + ".4.2" ),
            xIndex30( sFilePrefix + ".3.0" ),
            xIndex52( sFilePrefix + ".5.2" )
    {}

    // @todo resolve this via overloading the function @continue_here
    template <size_t D, size_t O> constexpr sps::Index<storage_t<D + O, false, true, O>>& getIndex( )
    {
        switch (D)
        {
            case 1:
                switch (O)
                {
                case 0:
                    return xIndex10;
                case 1:
                    return xIndex21;
                default:
                    throw std::runtime_error("invalid value for O");
                }
            case 2:
                switch (O)
                {
                case 0:
                    return xIndex20;
                case 1:
                    return xIndex31;
                case 2:
                    return xIndex42;
                default:
                    throw std::runtime_error("invalid value for O");
                }
            case 3:
                switch (O)
                {
                case 0:
                    return xIndex30;
                case 2:
                    return xIndex52;
                default:
                    throw std::runtime_error("invalid value for O");
                }
            
            default:
                throw std::runtime_error("invalid value for D");
        }
    }
};
} // namespace cm