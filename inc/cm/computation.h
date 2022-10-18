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
        SIZE
    };
    struct ComputeNode
    {
        std::string sNodeName;
        void ( ContactMapping::*fFunc )( void );
        std::vector<NodeNames> vIncomingFunctions;
        std::vector<std::vector<std::string>> vIncomingSession;
        std::vector<std::vector<std::string>> vIncomingRender;
        size_t uiLastUpdated;
    };

    size_t uiCurrTime;
    std::vector<ComputeNode> vGraph;

    json xRenderSettings, xSession;
    std::map<std::vector<std::string>, size_t> xRenderTime;
    std::map<std::vector<std::string>, size_t> xSessionTime;


    bool updateSettings( std::map<std::vector<std::string>, size_t>& xMap, json& xSettings, json& xOldSettings,
                         std::vector<std::string> vPrefix = { } )
    {
        bool bUpdated = false;
        for( auto& xEle : xSettings.items( ) )
        {
            std::vector<std::string> vPrefix2( vPrefix );
            vPrefix2.push_back( xEle.key( ) );

            if( xEle.value( ).is_object( ) )
                bUpdated = updateSettings( xMap, xEle.value( ), xOldSettings[ xEle.key( ) ], vPrefix2 ) || bUpdated;
            else if( xOldSettings[ xEle.key( ) ] != xEle.value( ) )
            {
                xMap[ vPrefix2 ] = uiCurrTime;
                bUpdated = true;
            }
        }
        if( bUpdated )
            xMap[ vPrefix ] = uiCurrTime;
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
                vOutOfDateReason.push_back( "change in session value " + vec2str( vSetting ) );

        for( std::vector<std::string> vSetting : xNode.vIncomingRender )
            if( xRenderTime.count( vSetting ) == 0 || xRenderTime[ vSetting ] > xNode.uiLastUpdated )
                vOutOfDateReason.push_back( "change in setting " + vec2str( vSetting ) );

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
    void setRenderSettings( json& xRenderSettings )
    {
        ++uiCurrTime;
        updateSettings( xRenderTime, xRenderSettings, this->xRenderSettings );
        this->xRenderSettings = xRenderSettings;
    }

    void setSession( json& xSession )
    {
        ++uiCurrTime;
        updateSettings( xSessionTime, xSession, this->xSession );
        this->xSession = xSession;
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
    size_t iInGroupSetting;

    std::vector<std::vector<size_t>> vvBinValues;
    std::vector<std::array<size_t, 2>> vvFlatValues;
    std::vector<std::array<double, 2>> vvNormalized;

    std::array<std::vector<std::vector<size_t>>, 2> vvCoverageValues;
    // outer array: column / row
    std::array<std::vector<size_t>, 2> vvFlatCoverageValues;

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
    void setFlatValues( );

    // replicates.h
    void setInGroup( );

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

  protected:
    virtual std::vector<std::array<double, 2>> normalizeBinominalTestTrampoline( std::vector<std::array<size_t, 2>>&,
                                                                                 size_t, size_t, double )
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
    }

    virtual ~ContactMapping( )
    {}

    void computeAll( json xRenderSettings, json xSession )
    {
        this->xRenderSettings = xRenderSettings;
        this->xSession = xSession;
        setBinSize( );
        setRenderArea( );

        setActiveChrom( );
        setAxisCoords( );

        setBinCoords( );

        setActiveReplicates( );
        setIntersectionType( );

        setSymmetry( );
        setInGroup( );

        setActiveCoverage( );
        setCoverageValues( );
        setFlatCoverageValues( );

        setBinValues( );
        setFlatValues( );

        setNormalized( );
    }
};


} // namespace cm
