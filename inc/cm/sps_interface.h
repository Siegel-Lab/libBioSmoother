#include "sps/index.h"
#include <type_traits>

namespace cm
{
template <bool CACHED> class SpsInterface
{
  private:
    template <size_t D, bool dependant_dim, bool uniform_overlay_grid, size_t orthope>
    using storage_t = typename std::conditional < CACHED,
          CachedTypeDef < D, dependant_dim, uniform_overlay_grid, orthope,
          DiskTypeDef<D, dependant_dim, uniform_overlay_grid, orthope>::type;

    std::tuple < std::tuple < sps::Index < storage_t<1, false, true, 0>, sps::Index<storage_t<1, false, true, 1>>,
        std::tuple < sps::Index < storage_t<2, false, true, 0>, sps::Index < storage_t<2, false, true, 1>,
        sps::Index<storage_t<2, false, true, 2>>,
        std::tuple < sps::Index<storage_t<3, false, true, 0>, sps::Index<storage_t<3, false, true, 2>>, > xIndices;

    SpsInterface( std::string sFilePrefix )
        : xIndices{ { { sFilePrefix + ".1.0" }, { sFilePrefix + ".1.1" } },
                    { { sFilePrefix + ".2.0" }, { sFilePrefix + ".2.1" }, { sFilePrefix + ".2.2" } },
                    { { sFilePrefix + ".3.0" }, { sFilePrefix + ".3.2" } } }
    {}

    template <size_t D, size_t O> sps::Index < storage_t<D, false, true, O> getIndex( )
    {
        constexpr std::array<std::array<size_t, 3>, 3> vConv = {
            { 0, 1, -1 },
            { 0, 1, 2 },
            { 0, -1, 1 }
        } return std::get<vConv[ D - 1 ][ O ]>( std::get<D - 1>( xIndices ) );
    }
};
} // namespace cm