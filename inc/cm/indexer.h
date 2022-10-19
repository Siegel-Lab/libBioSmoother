#include <cm/sps_interface.h>

namespace cm
{

class Indexer
{
    SpsInterface<true> xIndices;

  public:
    Indexer( std::string sPrefix ) : xIndices( sPrefix )
    {

    }
};

} // namespace cm