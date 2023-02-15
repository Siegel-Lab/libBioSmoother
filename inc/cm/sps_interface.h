#include "cm/annotation_index.h"
#include "sps/default.h"
#include "sps/index.h"
#include <type_traits>
#include <nlohmann/json.hpp>
#include <pybind11_json/pybind11_json.hpp>

#pragma once

using json = nlohmann::json;

namespace cm
{

class HasSession{
protected:
    std::string sFilePrefix;
    json xSession;
    
    std::string toString( std::vector<std::string>& vKeys )
    {
        std::string sPtr = "";
        for( std::string sKey : vKeys )
            sPtr += "/" + sKey;
        return sPtr;
    }

    nlohmann::json::json_pointer toPointer( std::vector<std::string>& vKeys )
    {
        return nlohmann::json::json_pointer( toString( vKeys ) );
    }

public:
    HasSession() {}

    HasSession(const std::string sFilePrefix) : 
        sFilePrefix(sFilePrefix), 
        xSession(json::parse( std::ifstream( sFilePrefix + "/session.json" ) )) 
    {}

    HasSession(const HasSession& rOther) : sFilePrefix(rOther.sFilePrefix), xSession(rOther.copySession()) {}

    template <typename T> T getSessionValue( std::vector<std::string> vKeys )
    {
        return this->xSession[ toPointer( vKeys ) ].get<T>( );
    }

    template <typename T> void setSessionValue( std::vector<std::string> vKeys, T xVal )
    {
        this->xSession[ toPointer( vKeys ) ] = xVal;
    }

    json copySession() const
    {
        return json::parse( this->xSession.dump( ) );
    }

    void saveSession( )
    {
        std::ofstream o( sFilePrefix + +"/session.json" );
        o << xSession << std::endl;
    }

    HasSession& operator=(const HasSession& rOther)
    {
        sFilePrefix = rOther.sFilePrefix;
        xSession = rOther.copySession();
        return *this;
    }
};

template <bool CACHED> class SpsInterface: public HasSession
{
  public:
    using anno_t = AnnotationDescIndex<DiskVecGenerator>;
    anno_t vAnno;

  private:
    static const bool BIN_SEARCH_SPARSE = false;
    static const size_t D = 6;
    static const size_t O = 2;

    using index_t = sps::Index<typename std::conditional<CACHED, CachedTypeDef<D, O, BIN_SEARCH_SPARSE>,
                                                         DiskTypeDef<D, O, BIN_SEARCH_SPARSE>>::type>;
    std::shared_ptr<index_t> pIndex;

  public:
    SpsInterface( std::string sFilePrefix, bool bWrite )
        : HasSession( sFilePrefix ), //
          vAnno( sFilePrefix + "/anno", bWrite ), //
          pIndex(
              std::make_shared<index_t>( sFilePrefix + "/" + std::to_string( D ) + "." + std::to_string( O ), bWrite ) )
    {}

    SpsInterface( std::string sFilePrefix ) :
        SpsInterface( sFilePrefix, false )
    {}

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

    std::vector<unsigned int> gridCount( size_t iDataSetId, std::array<std::vector<size_t>, D - O> vGrid,
                                         sps::IntersectionType xIntersect, size_t uiVerbosity = 1 )
    {
        std::vector<unsigned int> xRet;
        return xRet;
    }

    size_t count( size_t iDataSetId, std::array<size_t, D - O> vFrom, std::array<size_t, D - O> vTo,
                  sps::IntersectionType xIntersect, size_t uiVerbosity = 1 )
    {
        return pIndex->count( iDataSetId, vFrom, vTo, xIntersect, uiVerbosity );
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