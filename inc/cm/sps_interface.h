#include "cm/annotation_index.h"
#include "sps/default.h"
#include "sps/index.h"
#include <nlohmann/json.hpp>
#include <pybind11_json/pybind11_json.hpp>
#include <type_traits>

#pragma once

using json = nlohmann::json;

namespace cm
{

class HasSession
{
  protected:
    std::string sFilePrefix;
    json xSession;

    std::string replace(const std::string& sStr, const char cFrom, std::string sTo) {
        std::string sRet = "";
        sRet.reserve( sStr.size( ) * 2 );
        for( const char c : sStr )
        {
            if (c == cFrom)
                sRet += sTo;
            else
                sRet += std::string( 1, c );
        }
        return sRet;
    }

    std::string escape(const std::string& sKey)
    {
        return replace(replace(sKey, '~', "~0"), '/', "~1");
    }

    std::string toString( std::vector<std::string>& vKeys )
    {
        std::string sPtr = "";
        for( const std::string& sKey : vKeys )
            sPtr += "/" + escape(sKey);
        return sPtr;
    }

    nlohmann::json::json_pointer toPointer( std::vector<std::string>& vKeys )
    {
        return nlohmann::json::json_pointer( toString( vKeys ) );
    }

  public:
    HasSession( )
    {}

    HasSession( const std::string sFilePrefix )
        : sFilePrefix( sFilePrefix ), xSession( json::parse( std::ifstream( sFilePrefix + "/session.json" ) ) )
    {}

    HasSession( const std::shared_ptr<HasSession> rOther )
        : sFilePrefix( rOther != nullptr ? rOther->sFilePrefix : std::string( "" ) ),
          xSession( rOther != nullptr ? rOther->copySession( ) : json::parse( "{}" ) )
    {}

    template <typename T> T getSessionValue( std::vector<std::string> vKeys )
    {
        try
        {
            return this->xSession[ toPointer( vKeys ) ].get<T>( );
        }
        catch( ... )
        {
            throw std::runtime_error( std::string( "error querying session value: " ) + toString( vKeys ) );
        }
    }

    template <typename T> void setSessionValue( std::vector<std::string> vKeys, T xVal )
    {
        this->xSession[ toPointer( vKeys ) ] = xVal;
    }

    json copySession( ) const
    {
        return json::parse( this->xSession.dump( ) );
    }

    void saveSession( )
    {
        std::ofstream o( sFilePrefix + "/session.json" );
        o << xSession << std::endl;
    }

    HasSession& operator=( const HasSession& rOther )
    {
        sFilePrefix = rOther.sFilePrefix;
        xSession = rOther.copySession( );
        return *this;
    }
};

template <bool CACHED> class SpsInterface : public HasSession
{
  public:
    using anno_t = AnnotationDescIndex<DiskVecGenerator>;
    anno_t vAnno;

  private:
    static const size_t D = 8;
    static const size_t O = 2;
    static const size_t uiBiasAccuracy = 1000000;
    static const size_t uiOrthotopePower = 10;
#ifdef WITH_STXXL
    using index_t = sps::Index<typename std::conditional<CACHED, PowerCachedTypeDef<D, O, uiOrthotopePower>,
                                                         PowerDiskTypeDef<D, O, uiOrthotopePower>>::type>;
#else
    using index_t = sps::Index<PowerDiskTypeDef<D, O, uiOrthotopePower>>;
#endif
    std::shared_ptr<index_t> pIndex;
    const std::string sFilePrefix;

  public:
    using coordinate_t = typename index_t::coordinate_t;

    SpsInterface( std::string sFilePrefix, bool bWrite )
        : HasSession( sFilePrefix ), //
          vAnno( sFilePrefix + "/anno", bWrite ), //
          pIndex( std::make_shared<index_t>( sFilePrefix + "/sps", bWrite ) ),
          sFilePrefix( sFilePrefix )
    {}

    SpsInterface( std::string sFilePrefix ) : SpsInterface( sFilePrefix, false )
    {}

    const std::string getFilePrefix( ) const
    {
        return sFilePrefix;
    }

    bool loaded( )
    {
        return pIndex != nullptr;
    }

    void insert( std::array<uint64_t, D - O> vStart, std::array<uint64_t, D - O> vEnd, int iValue )
    {
        pIndex->addPoint( vStart, vEnd, static_cast<size_t>( iValue ) );
    }

    void insert( std::vector<uint64_t> vStart, std::vector<uint64_t> vEnd, int iValue )
    {
        assert( vStart.size( ) == D - O );
        assert( vEnd.size( ) == D - O );

        std::array<uint64_t, D - O> aStart;
        std::copy_n( vStart.begin( ), D - O, aStart.begin( ) );
        std::array<uint64_t, D - O> aEnd;
        std::copy_n( vEnd.begin( ), D - O, aEnd.begin( ) );

        insert( aStart, aEnd, iValue );
    }


    std::vector<uint32_t> gridCount( size_t iDataSetId, std::array<std::vector<coordinate_t>, D - O> vGrid,
                                     sps::IntersectionType xIntersect, size_t uiVerbosity = 1 )
    {
        return pIndex->gridCount( iDataSetId, vGrid, xIntersect, uiVerbosity );
    }

    size_t count( size_t iDataSetId, std::array<coordinate_t, D - O> vFrom, std::array<coordinate_t, D - O> vTo,
                  sps::IntersectionType xIntersect, bool bNoPoints, size_t uiVerbosity = 1 )
    {
        return pIndex->count( iDataSetId, vFrom, vTo, xIntersect, bNoPoints, uiVerbosity );
    }

    size_t generate( double fFac = -1, size_t uiVerbosity = 1 )
    {
        return pIndex->generate( fFac, uiVerbosity );
    }

    void clearPointsAndDesc( )
    {
        pIndex->clearCorners( );
    }

    size_t getNumOverlays( size_t iDataSetId )
    {
        return pIndex->getNumOverlays( iDataSetId );
    }
};


} // namespace cm