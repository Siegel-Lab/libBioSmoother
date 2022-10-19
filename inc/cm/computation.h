#include "cm/sps_interface.h"
#include <nlohmann/json.hpp>

#pragma once

using json = nlohmann::json;

namespace cm
{

struct ChromDesc
{
    std::string sName;
    size_t uiLength;
};

struct AxisCoord
{
    std::string sChromosome;
    size_t uiScreenPos;
    size_t uiIndexPos;
    size_t uiSize;
};

struct BinCoord
{
    std::string sChromosomeX, sChromosomeY;
    size_t uiScreenX, uiScreenY;
    size_t uiIndexX, uiIndexY;
    size_t uiW, uiH;
};


class ContactMapping
{
  public:
    size_t uiVerbosity = 1;

  private:
    enum NodeNames
    {
        BinSize = 0,
        RenderArea,
        ActiveChrom,
        AxisCoords,
        Symmetry,
        BinCoords,
        IntersectionType,
        ActiveReplicates,
        ActiveCoverage,
        CoverageValues,
        FlatCoverageValues,
        BinValues,
        FlatValues,
        InGroup,
        Normalized,
        Colors,
        BetweenGroup,
        Combined,
        Colored,
        SIZE
    };
    struct ComputeNode
    {
        std::string sNodeName;
        void ( ContactMapping::*fFunc )( void );
        std::vector<NodeNames> vIncomingFunctions;
        std::vector<std::vector<std::string>> vIncomingSession;
        size_t uiLastUpdated;
    };

    size_t uiCurrTime;
    std::vector<ComputeNode> vGraph;

    json xSession;
    std::map<std::vector<std::string>, size_t> xSessionTime;


    bool updateSettings( json& xSettings, json& xOldSettings, std::vector<std::string> vPrefix = { } )
    {
        bool bUpdated = false;
        for( auto& xEle : xSettings.items( ) )
        {
            if( xEle.key( ) != "previous" && xEle.key( ) != "next" )
            {
                std::vector<std::string> vPrefix2( vPrefix );
                vPrefix2.push_back( xEle.key( ) );

                if( xEle.value( ).is_object( ) )
                    bUpdated = updateSettings( xEle.value( ), xOldSettings[ xEle.key( ) ], vPrefix2 ) || bUpdated;
                else if( xOldSettings[ xEle.key( ) ] != xEle.value( ) )
                {
                    xSessionTime[ vPrefix2 ] = uiCurrTime;
                    bUpdated = true;
                }
            }
        }
        if( bUpdated )
            xSessionTime[ vPrefix ] = uiCurrTime;
        return bUpdated;
    }

    void registerNode( NodeNames xName, ComputeNode xNode )
    {
        vGraph[ xName ] = xNode;
    }

    std::string vec2str( std::vector<std::string> vVec )
    {
        std::string sRet;
        for( std::string sStr : vVec )
        {
            if( sRet.size( ) > 0 )
                sRet += ".";
            sRet += sStr;
        }
        return sRet;
    }

    bool update( NodeNames xNodeName )
    {
        ComputeNode& xNode = vGraph[ xNodeName ];
        if( xNode.uiLastUpdated == uiCurrTime )
            return true;

        std::vector<std::string> vOutOfDateReason;
        auto p1 = std::chrono::high_resolution_clock::now( );
        for( NodeNames xPred : xNode.vIncomingFunctions )
            if( update( xPred ) )
                vOutOfDateReason.push_back( "change in previous node " + vGraph[ xPred ].sNodeName );
        auto p2 = std::chrono::high_resolution_clock::now( );
        auto pms_int = std::chrono::duration_cast<std::chrono::milliseconds>( p2 - p1 );

        for( std::vector<std::string> vSetting : xNode.vIncomingSession )
            if( xSessionTime.count( vSetting ) == 0 || xSessionTime[ vSetting ] > xNode.uiLastUpdated )
                vOutOfDateReason.push_back( "change in variable " + vec2str( vSetting ) );

        if( vOutOfDateReason.size( ) > 0 )
        {
            if( uiVerbosity >= 3 )
            {
                std::cout << "running node " << xNode.sNodeName << " due to: " << std::endl;
                for( std::string& sStr : vOutOfDateReason )
                    std::cout << "\t- " << sStr << std::endl;
            }
            auto t1 = std::chrono::high_resolution_clock::now( );

            ( this->*xNode.fFunc )( );

            auto t2 = std::chrono::high_resolution_clock::now( );
            auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 );
            if( uiVerbosity >= 1 )
                std::cout << ms_int.count( ) << "ms (" << xNode.sNodeName << ")" << std::endl;
            if( uiVerbosity >= 2 )
                std::cout << pms_int.count( ) << "ms (predecessors)" << std::endl << std::endl;
        }

        xNode.uiLastUpdated = uiCurrTime;
        return vOutOfDateReason.size( ) > 0;
    }


  public:
    void setSession( json& xSession )
    {
        ++uiCurrTime;
        updateSettings( xSession, this->xSession );
        this->xSession = xSession;
    }

    const json& getSession( )
    {
        return this->xSession;
    }

    template <typename T> const T& getSession( std::vector<std::string> vKeys )
    {
        json& rCurr = this->xSession;
        for( std::string sKey : vKeys )
            rCurr = rCurr[ sKey ];
        return rCurr.get<T>( );
    }

    template <typename T> void changeSession( std::vector<std::string> vKeys, T xVal )
    {
        ++uiCurrTime;

        // remove old redo
        this->xSession[ "next" ] = nullptr;

        // deep copy json
        json xNew = json::parse( this->xSession.dump( ) ); // is this necessary or is a normal copy enough
        xNew[ "previous" ] = this->xSession; // back up previous

        // change setting
        json& rCurr = this->xSession;
        std::vector<std::string> vCurrKey;
        for( std::string sKey : vKeys )
        {
            xSessionTime[ vCurrKey ] = uiCurrTime;
            rCurr = rCurr[ sKey ];
            vCurrKey.push_back( sKey );
        }
        xSessionTime[ vCurrKey ] = uiCurrTime;
        rCurr = xVal;

        // make sure undo's and redo's do not pile up infinitely
        int uiMaxRedo = this->xSession[ "settings" ][ "interface" ][ "max_undo" ][ "val" ].get<int>( );
        json& rToNull = this->xSession;
        while( --uiMaxRedo >= 0 && !rToNull[ "previous" ].is_null( ) )
            rToNull = rToNull[ "previous" ];
        rToNull[ "previous" ] = nullptr;
    }

    bool hasUndo( )
    {
        return !this->xSession[ "previous" ].is_null( );
    }

    void undo( )
    {
        if( hasUndo( ) )
        {
            json xNow = json::parse( this->xSession.dump( ) ); // is this necessary of is a normal copy enough
            setSession( this->xSession[ "previous" ] );
            this->xSession[ "next" ] = xNow;
        }
    }

    bool hasRedo( )
    {
        return !this->xSession[ "next" ].is_null( );
    }

    void redo( )
    {
        if( hasRedo( ) )
            setSession( this->xSession[ "next" ] );
    }

  private:
    SpsInterface<false> xIndices;

    size_t uiBinWidth, uiBinHeight;
    int64_t iStartX, iStartY, iEndX, iEndY;

    std::array<std::vector<ChromDesc>, 2> vActiveChromosomes;
    std::array<std::vector<AxisCoord>, 2> vAxisCords;

    std::vector<std::string> vActiveReplicates;

    sps::IntersectionType xIntersect;

    std::vector<std::array<BinCoord, 2>> vBinCoords;

    std::array<std::vector<std::pair<std::string, bool>>, 2> vActiveCoverage;
    size_t uiSymmetry;
    size_t iInGroupSetting, iBetweenGroupSetting;

    std::vector<std::vector<size_t>> vvBinValues;
    std::vector<std::array<size_t, 2>> vvFlatValues;
    std::vector<std::array<double, 2>> vvNormalized;
    std::vector<double> vCombined;
    std::vector<std::string> vColored;

    std::array<std::vector<std::vector<size_t>>, 2> vvCoverageValues;
    // outer array: column / row
    std::array<std::vector<double>, 2> vvFlatCoverageValues;

    std::vector<std::string> vColorPalette;

    // bin_size.h
    size_t nextEvenNumber( double fX );


    // bin_size.h
    void setBinSize( );
    // bin_size.h
    void setRenderArea( );

    // bin_size.h
    void regBinSize( );

    // coords.h
    void setActiveChrom( );
    // coords.h
    void setAxisCoords( );

    // coords.h
    void setSymmetry( );
    // coords.h
    void setBinCoords( );

    // coords.h
    void regCoords( );

    // replicates.h
    void setIntersectionType( );

    // replicates.h
    void setActiveReplicates( );

    // replicates.h
    void regReplicates( );

    // coverage.h
    void setActiveCoverage( );

    // coverage.h
    void setCoverageValues( );

    // coverage.h
    void setFlatCoverageValues( );

    // coverage.h
    void regCoverage( );

    // replicates.h
    size_t symmetry( size_t uiA, size_t uiB );

    // replicates.h
    void setBinValues( );

    // replicates.h
    size_t getFlatValue( std::vector<size_t> vCollected );
    // replicates.h
    double getMixedValue( double uiA, double uiB );

    // replicates.h
    void setFlatValues( );

    // replicates.h
    void setInGroup( );

    // replicates.h
    void setBetweenGroup( );

    // normalization.h
    void doNotNormalize( );
    // normalization.h
    void normalizeTracks( );
    // normalization.h
    void normalizeSize( size_t uiSize );
    // normalization.h
    void normalizeBinominalTest( );
    // normalization.h
    void normalizeIC( );
    // normalization.h
    void setNormalized( );

    // normalization.h
    void regNormalization( );

    // colors.h
    void setColors( );
    // colors.h
    void setCombined( );
    // colors.h
    void setColored( );
    // colors.h
    double logScale( double );
    // colors.h
    size_t colorRange( double fX );

    // colors.h
    void regColors( );

  protected:
    virtual std::vector<std::array<double, 2>> normalizeBinominalTestTrampoline( std::vector<std::array<size_t, 2>>&,
                                                                                 size_t, size_t, double )
    {
        throw std::logic_error( "Function not implemented" );
    }

    virtual std::vector<std::string> colorPalette( std::string )
    {
        throw std::logic_error( "Function not implemented" );
    }

  public:
    ContactMapping( std::string sPrefix ) : uiCurrTime( 0 ), vGraph( NodeNames::SIZE ), xIndices( sPrefix )
    {
        regBinSize( );
        regCoords( );
        regReplicates( );
        regCoverage( );
        regNormalization( );
        regColors( );
    }

    virtual ~ContactMapping( )
    {}
    
    // colors.h
    std::vector<std::string> getColors();
};


} // namespace cm
