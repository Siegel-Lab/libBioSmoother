#include "cm/sps_interface.h"
#include <cctype>
#include <mutex>
#include <nlohmann/json.hpp>
#include <pybind11_json/pybind11_json.hpp>

#pragma once

using json = nlohmann::json;

#define CANCEL_RETURN                                                                                                  \
    if( this->bCancel )                                                                                                \
    return false
#define END_RETURN return true

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
    size_t uiScreenSize;
    size_t uiIndexSize;
    bool bFiltered = false;
};

struct DecayCoord
{
    std::string sChromosome;
    int64_t iFrom;
    int64_t iTo;

    friend bool operator<( const DecayCoord& l, const DecayCoord& r )
    {
        if( l.sChromosome != r.sChromosome )
            return l.sChromosome < r.sChromosome;
        if( l.iFrom != r.iFrom )
            return l.iFrom < r.iFrom;
        if( l.iTo != r.iTo )
            return l.iTo < r.iTo;
        return false;
    }
};

struct AxisRegion : AxisCoord
{
    size_t uiCoordStartIdx;
    size_t uiNumCoords;
};

struct BinCoord
{
    std::string sChromosomeX, sChromosomeY;
    size_t uiScreenX, uiScreenY;
    size_t uiIndexX, uiIndexY;
    size_t uiScreenW, uiScreenH;
    size_t uiIndexW, uiIndexH;
    size_t uiDecayCoordIndex = std::numeric_limits<size_t>::max( );
};

struct Annotation
{
    std::string sInfo;
    bool bForw;
    double fScreenX, fScreenY;
    size_t uiIndexX, uiIndexY;
    size_t uiRow;
    std::string sChromosome;
};


class PartialQuarry
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
        FilteredCoords,
        Symmetry,
        BinCoords,
        DecayCoords,
        FlatDecay,
        IntersectionType,
        ActiveReplicates,
        ActiveCoverage,
        CoverageValues,
        FlatCoverageValues,
        NormalizeCoverageValues,
        CombinedCoverageValues,
        BinValues,
        DecayValues,
        FlatValues,
        InGroup,
        Normalized,
        DistDepDecayRemoved,
        DecayCDS,
        Colors,
        BetweenGroup,
        Combined,
        Colored,
        ActivateAnnotation,
        AnnotationValues,
        AnnotationCDS,
        AnnotationColors,
        ActivateAnnotationCDS,
        HeatmapCDS,
        HeatmapExport,
        TrackExport,
        Scaled,
        Ticks,
        Tracks,
        Divided,
        Palette,
        LCS,
        RankedSlicesCDS,
        CanvasSize,
        SIZE
    };
    struct ComputeNode
    {
        std::string sNodeName;
        bool ( PartialQuarry::*fFunc )( void );
        std::vector<NodeNames> vIncomingFunctions;
        std::vector<std::vector<std::string>> vIncomingSession;
        std::vector<std::vector<std::string>> vSessionsIncomingInPrevious;
    };

    struct ComputeNodeData
    {
        size_t uiLastUpdated = 0;
        std::string sError = "";
        bool bRegistered = false;
        bool bTerminal = true;
    };

    std::string sPrefix;
    size_t uiCurrTime;
    std::vector<ComputeNode> vGraph;
    std::vector<ComputeNodeData> vGraphData;
    NodeNames xCurrNodeName = NodeNames::SIZE;

    json xSession;
    std::map<std::vector<std::string>, size_t> xSessionTime;

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

    bool updateSettings( const json& xNewSettings, std::vector<std::string> vPrefix = { } )
    {
        bool bUpdated = false;
        for( auto& xEle : xNewSettings[ toPointer( vPrefix ) ].items( ) )
        {
            if( xEle.key( ) != "previous" && xEle.key( ) != "next" )
            {
                std::vector<std::string> vPrefix2( vPrefix );
                vPrefix2.push_back( xEle.key( ) );

                if( xNewSettings[ toPointer( vPrefix2 ) ].is_object( ) )
                    bUpdated = updateSettings( xNewSettings, vPrefix2 ) || bUpdated;
                else if( this->xSession[ toPointer( vPrefix2 ) ] != xNewSettings[ toPointer( vPrefix2 ) ] )
                {
                    if( xSessionTime.count( vPrefix2 ) > 0 )
                        xSessionTime[ vPrefix2 ] = uiCurrTime;
                    bUpdated = true;
                }
            }
        }
        if( bUpdated && xSessionTime.count( vPrefix ) > 0 )
            xSessionTime[ vPrefix ] = uiCurrTime;
        return bUpdated;
    }

    void registerNode( NodeNames xName, ComputeNode xNode )
    {
        assert( !vGraphData[ xName ].bRegistered );

        for( std::vector<std::string>& rKey : xNode.vIncomingSession )
            xSessionTime[ rKey ] = uiCurrTime + 1;
        for( NodeNames& rIncoming : xNode.vIncomingFunctions )
            vGraphData[ rIncoming ].bTerminal = false;
        vGraph[ xName ] = xNode;
        vGraphData[ xName ].uiLastUpdated = uiCurrTime;
        vGraphData[ xName ].bRegistered = true;
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

    void setError( std::string sError )
    {
        assert( xCurrNodeName != NodeNames::SIZE );
        vGraphData[ xCurrNodeName ].sError = sError;
    }

    size_t update_helper( NodeNames xNodeName )
    {
        ComputeNode& xNode = vGraph[ xNodeName ];
        ComputeNodeData& xNodeData = vGraphData[ xNodeName ];
        if( xNodeData.uiLastUpdated < uiCurrTime )
        {
            std::vector<std::string> vOutOfDateReason{ };
            auto p1 = std::chrono::high_resolution_clock::now( );
            for( NodeNames xPred : xNode.vIncomingFunctions )
            {
                size_t uiUpdateTime = update_helper( xPred );
                if( uiUpdateTime == 0 )
                {
                    if( uiVerbosity >= 1 )
                        std::cout << xNode.sNodeName << " was cancelled" << std::endl;
                    return 0;
                }
                else if( uiUpdateTime > xNodeData.uiLastUpdated )
                    vOutOfDateReason.push_back( "change in previous node " + vGraph[ xPred ].sNodeName );
            }
            auto p2 = std::chrono::high_resolution_clock::now( );
            auto pms_int = std::chrono::duration_cast<std::chrono::milliseconds>( p2 - p1 );

            for( std::vector<std::string> vSetting : xNode.vIncomingSession )
                if( xSessionTime[ vSetting ] > xNodeData.uiLastUpdated )
                    vOutOfDateReason.push_back( "change in variable " + vec2str( vSetting ) +
                                                " (last set: " + std::to_string( xSessionTime[ vSetting ] ) +
                                                "; now: " + std::to_string( uiCurrTime ) + ")" );

            if( vOutOfDateReason.size( ) > 0 )
            {
                if( uiVerbosity >= 3 )
                {
                    std::cout << "running node " << xNode.sNodeName << " due to: " << std::endl;
                    for( std::string& sStr : vOutOfDateReason )
                        std::cout << "\t- " << sStr << std::endl;
                    std::cout << "last updated: " << xNodeData.uiLastUpdated << "; now: " << uiCurrTime << std::endl;
                }

                xCurrNodeName = xNodeName;
                setError( "" );

                auto t1 = std::chrono::high_resolution_clock::now( );
                bool bUpdateDone = ( this->*xNode.fFunc )( );
                auto t2 = std::chrono::high_resolution_clock::now( );

                xCurrNodeName = NodeNames::SIZE;
                if( !bUpdateDone )
                {
                    if( uiVerbosity >= 1 )
                        std::cout << xNode.sNodeName << " was cancelled" << std::endl;
                    // update was cancelled -> we do not know how messed up the result form the last call was
                    xNodeData.uiLastUpdated = 0;
                    return 0;
                }

                auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 );
                if( uiVerbosity >= 1 )
                    print( std::to_string( ms_int.count( ) ) + " ms (" + xNode.sNodeName + ")" );
                if( uiVerbosity >= 2 )
                    std::cout << pms_int.count( ) << " ms (predecessors)" << std::endl << std::endl;

                xNodeData.uiLastUpdated = uiCurrTime;
            }
            else if( uiVerbosity >= 3 )
                std::cout << "no reason found to run node " << xNode.sNodeName << std::endl;
        }
        else if( uiVerbosity >= 3 )
            std::cout << "node up to date from previous check " << xNode.sNodeName << std::endl;

        size_t uiRet = xNodeData.uiLastUpdated;
        // node must be up to date now
        xNodeData.uiLastUpdated = uiCurrTime;

        return uiRet;
    }


    bool update_no_throw( NodeNames xNodeName )
    {
        // allow python multithreading during this call
        pybind11::gil_scoped_release release;

        // make sure that no two threads can enter here (basically allow multithreading other than this)
        std::lock_guard<std::mutex> xGuard( xUpdateMutex );

        bCancel = false;
        size_t uiUpdateTime = update_helper( xNodeName );

        return uiUpdateTime != 0;
    }

    void update( NodeNames xNodeName )
    {
        if( !update_no_throw( xNodeName ) )
            throw std::runtime_error( "update was cancelled" );
    }


    void checkUnnecessaryDependency( NodeNames xNodeCurr, std::map<NodeNames, std::string>& rDependents,
                                     std::map<std::vector<std::string>, std::string>& rDependentSettings )
    {
        for( std::vector<std::string>& rSetting : vGraph[ xNodeCurr ].vIncomingSession )
            rDependentSettings[ rSetting ] = vGraph[ xNodeCurr ].sNodeName;

        for( NodeNames& rPrev : vGraph[ xNodeCurr ].vIncomingFunctions )
        {
            rDependents[ rPrev ] = vGraph[ xNodeCurr ].sNodeName;
            checkUnnecessaryDependency( rPrev, rDependents, rDependentSettings );
        }
    }

    void checkUnnecessaryDependencies( )
    {
        std::set<std::string> xNodeNames;
        for( size_t uiNodeName = 0; uiNodeName < NodeNames::SIZE; uiNodeName++ )
        {
            if( xNodeNames.count( vGraph[ uiNodeName ].sNodeName ) > 0 )
                std::cerr << "two nodes have been registered under the same name: " << vGraph[ uiNodeName ].sNodeName
                          << std::endl;
            xNodeNames.insert( vGraph[ uiNodeName ].sNodeName );

            if( !vGraphData[ uiNodeName ].bRegistered )
                std::cerr << "Node has not been registered. id: " << uiNodeName << std::endl;
            else
            {
                std::map<NodeNames, std::string> xDependents;
                std::map<std::vector<std::string>, std::string> xDependentSettings;
                for( NodeNames& rPrev : vGraph[ uiNodeName ].vIncomingFunctions )
                    checkUnnecessaryDependency( rPrev, xDependents, xDependentSettings );

                for( std::vector<std::string>& rSetting : vGraph[ uiNodeName ].vIncomingSession )
                    if( xDependentSettings.count( rSetting ) > 0 )
                        std::cerr << "unnecessary session dependency warning: " << vGraph[ uiNodeName ].sNodeName
                                  << " depends on " << toString( rSetting ) << ", but the previous node "
                                  << xDependentSettings[ rSetting ] << " shares this dependency." << std::endl;

                for( std::vector<std::string>& rSetting : vGraph[ uiNodeName ].vSessionsIncomingInPrevious )
                    if( xDependentSettings.count( rSetting ) == 0 )
                        std::cerr << "missing session dependency warning: " << vGraph[ uiNodeName ].sNodeName
                                  << " declares " << toString( rSetting )
                                  << " a dependency of a previous node, but no such previous node could be found."
                                  << std::endl;

                for( NodeNames& rPrev : vGraph[ uiNodeName ].vIncomingFunctions )
                    if( xDependents.count( rPrev ) > 0 )
                        std::cerr << "unnecessary session dependency warning: " << vGraph[ uiNodeName ].sNodeName
                                  << " depends on " << vGraph[ rPrev ].sNodeName << ", but the previous node "
                                  << xDependents[ rPrev ] << " shares this dependency" << std::endl;
            }
        }
    }

  public:
    std::string getError( )
    {
        std::string sRet = "";
        for( size_t uiNodeName = 0; uiNodeName < NodeNames::SIZE; uiNodeName++ )
            if( vGraphData[ uiNodeName ].sError.size( ) > 0 )
                sRet += vGraph[ uiNodeName ].sNodeName + ": " + vGraphData[ uiNodeName ].sError + "\n";
        if( sRet.size( ) > 0 )
            sRet = sRet.substr( 0, sRet.size( ) - 1 );
        return sRet;
    }

    void cancel( )
    {
        bCancel = true;
    }

    void updateCDS( )
    {
        bool bContinue;
        do
        {
            bContinue = false;
            for( auto& rNode : { HeatmapCDS, Tracks, AnnotationCDS, ActivateAnnotationCDS, Ticks, Tracks, Palette,
                                 DecayCDS, RankedSlicesCDS } )
                if( !update_no_throw( rNode ) )
                {
                    bContinue = true;
                    break;
                }
        } while( bContinue );
    }

    void setSession( const json& xSession )
    {
        ++uiCurrTime;
        updateSettings( xSession );
        this->xSession = xSession;
    }

    const json& getSession( )
    {
        return this->xSession;
    }

    template <typename T> T getValue( std::vector<std::string> vKeys )
    {
#ifndef NDEBUG
        if( xCurrNodeName != NodeNames::SIZE )
        {
            // currently computing a node
            bool bFound = false;

            for( std::vector<std::vector<std::string>> vSettings :
                 { vGraph[ xCurrNodeName ].vIncomingSession, vGraph[ xCurrNodeName ].vSessionsIncomingInPrevious } )
                if( !bFound )
                    for( std::vector<std::string>& rSetting : vSettings )
                        if( vKeys.size( ) >= rSetting.size( ) )
                        {
                            bool bMatches = true;
                            for( size_t uiI = 0; uiI < rSetting.size( ) && bMatches; uiI++ )
                                bMatches = rSetting[ uiI ] == vKeys[ uiI ];
                            if( bMatches )
                            {
                                bFound = true;
                                break;
                            }
                        }

            if( !bFound )
                std::cout << "undeclared session dependency warning: " << vGraph[ xCurrNodeName ].sNodeName
                          << " accesses " << toPointer( vKeys ) << " but does not declare this." << std::endl;
        }
#endif
        return this->xSession[ toPointer( vKeys ) ].get<T>( );
    }

    template <typename T> void setValue( std::vector<std::string> vKeys, T xVal )
    {
        if( this->xSession[ toPointer( vKeys ) ].is_null( ) || this->xSession[ toPointer( vKeys ) ] != json( xVal ) )
        {
            if( this->xSessionTime.count( vKeys ) > 0 )
            {
                ++uiCurrTime;

                // remove old redo
                this->xSession[ "next" ] = nullptr;

                // deep copy json
                this->xSession[ "previous" ] = json::parse( this->xSession.dump( ) ); // back up previous

                // change setting
                this->xSession[ toPointer( vKeys ) ] = xVal;

                // update time
                std::vector<std::string> vCurrKey;
                if( xSessionTime.count( vCurrKey ) > 0 )
                    xSessionTime[ vCurrKey ] = uiCurrTime;
                for( std::string sKey : vKeys )
                {
                    vCurrKey.push_back( sKey );
                    if( xSessionTime.count( vCurrKey ) > 0 )
                        xSessionTime[ vCurrKey ] = uiCurrTime;
                }

                // make sure undo's and redo's do not pile up infinitely
                if( this->xSession[ "settings" ].is_object( ) )
                {
                    int uiMaxRedo = this->xSession[ "settings" ][ "interface" ][ "max_undo" ][ "val" ].get<int>( );
                    std::vector<std::string> vPrevKey{ "previous" };
                    while( --uiMaxRedo >= 0 && !this->xSession[ toPointer( vPrevKey ) ].is_null( ) )
                        vPrevKey.push_back( "previous" );
                    this->xSession[ toPointer( vPrevKey ) ] = nullptr;
                }
            }
            else
                this->xSession[ toPointer( vKeys ) ] = xVal;
        }
    }

    bool hasUndo( )
    {
        return !this->xSession[ "previous" ].is_null( );
    }

    void undo( )
    {
        if( hasUndo( ) )
        {
            json xNow = json::parse( this->xSession.dump( ) );
            json xPrev = json::parse( this->xSession[ "previous" ].dump( ) );
            setSession( xPrev );
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
    std::vector<std::array<DecayCoord, 2>> vDistDepDecCoords;
    std::array<std::vector<AxisRegion>, 2> vAxisRegions;

    std::vector<std::string> vActiveReplicates;
    std::vector<size_t> vActiveReplicatesTotal;
    std::array<std::vector<size_t>, 2> vInGroup;

    sps::IntersectionType xIntersect;

    std::vector<std::array<BinCoord, 2>> vBinCoords;

    std::array<std::vector<std::pair<std::string, bool>>, 2> vActiveCoverage;
    std::array<std::vector<size_t>, 2> vActiveCoverageTotal;

    size_t uiSymmetry;
    size_t iInGroupSetting, iBetweenGroupSetting;

    std::vector<std::vector<size_t>> vvBinValues;
    std::vector<std::vector<double>> vvDecayValues;
    std::vector<std::array<size_t, 2>> vvFlatValues;
    std::vector<std::array<double, 2>> vvFlatDecay;
    std::array<size_t, 2> vvFlatTotal;
    std::vector<std::array<double, 2>> vvNormalized;
    std::vector<double> vCombined;
    std::vector<double> vDivided;
    std::vector<double> vScaled;
    std::vector<std::string> vColored;

    std::string sBackgroundColor;

    std::array<std::vector<std::vector<size_t>>, 2> vvCoverageValues;
    std::array<std::array<std::vector<size_t>, 3>, 2> vInGroupCoverage;
    std::array<std::array<std::vector<size_t>, 2>, 2> vNormCoverage;
    // outer array: column / row
    std::array<std::vector<std::array<size_t, 2>>, 2> vvFlatCoverageValues;
    std::array<std::array<size_t, 2>, 2> vFlatCoverageTotal;
    std::array<std::vector<std::array<double, 2>>, 2> vvNormalizedCoverageValues;
    std::array<std::vector<double>, 2> vvCombinedCoverageValues;
    std::array<std::vector<std::array<size_t, 2>>, 2> vFlatNormValues;

    std::vector<std::string> vColorPalette;
    pybind11::list vRenderedPalette;
    std::vector<std::string> vColorPaletteAnnotation;
    std::vector<std::string> vColorPaletteAnnotationDark;
    std::array<std::vector<std::string>, 2> vActiveAnnotation;
    std::array<std::vector<std::string>, 2> vFilterAnnotation;
    std::array<std::vector<std::pair<std::vector<size_t>, std::vector<Annotation>>>, 2> vAnnotationValues;
    std::array<size_t, 2> vMaxAnnoRows;
    std::array<std::vector<size_t>, 2> vMaxRowsPerAnno;
    std::array<pybind11::dict, 2> vAnnotationCDS;
    std::array<pybind11::list, 2> vActiveAnnotationCDS;
    pybind11::dict xHeatmapCDS;
    pybind11::dict xDistDepDecCDS;
    std::array<pybind11::dict, 2> vRankedSliceCDS;
    std::vector<std::tuple<std::string, size_t, size_t, std::string, size_t, size_t, double>> vHeatmapExport;
    std::array<std::vector<std::string>, 2> vTrackExportNames;
    std::array<std::vector<std::tuple<std::string, size_t, size_t, std::vector<double>>>, 2> vTrackExport;
    std::array<pybind11::dict, 2> xTicksCDS;
    std::array<pybind11::dict, 2> xTracksCDS;
    std::array<pybind11::list, 2> vTickLists;
    std::array<size_t, 2> vCanvasSize;
    std::array<std::array<double, 2>, 2> vvMinMaxTracks;
    double fMax, fMin;
    double fDataMax, fDataMin;
    size_t uiLogestCommonSuffix;

    bool bCancel = false;
    std::mutex xUpdateMutex{ };


    // bin_size.h
    size_t nextEvenNumber( double );


    // bin_size.h
    bool setBinSize( );
    // bin_size.h
    bool setRenderArea( );
    // bin_size.h
    std::string readableBp( size_t );

    // bin_size.h
    void regBinSize( );

    // coords.h
    bool setActiveChrom( );
    // coords.h
    bool setAxisCoords( );
    // coords.h
    bool setFilteredCoords( );

    // coords.h
    bool setSymmetry( );
    // coords.h
    bool setBinCoords( );
    // coords.h
    bool setDecayCoords( );
    // coords.h
    bool setLCS( );
    // coords.h
    bool setTicks( );
    // coords.h
    bool setCanvasSize( );

    // coords.h
    void regCoords( );

    // replicates.h
    bool setIntersectionType( );

    // replicates.h
    bool setActiveReplicates( );

    // replicates.h
    void regReplicates( );

    // coverage.h
    bool setActiveCoverage( );

    // coverage.h
    bool setCoverageValues( );
    // coverage.h
    bool setNormalizeCoverageValues( );
    // coverage.h
    bool setCombinedCoverageValues( );
    // coverage.h
    bool setTracks( );
    // coverage.h
    bool setTrackExport( );
    // coverage.h
    size_t getCoverageFromRepl( bool,
                                bool,
                                bool,
                                size_t,
                                size_t,
                                const AxisCoord&,
                                const json&,
                                bool,
                                size_t,
                                const std::vector<ChromDesc>&,
                                bool,
                                bool );
    // coverage.h
    std::tuple<size_t, int64_t, size_t, size_t> makeHeapTuple( bool, bool, size_t, size_t, bool, bool, const AxisCoord&,
                                                               int64_t, size_t, size_t );

    // coverage.h
    bool setFlatCoverageValues( );
    bool setRankedSlicesCDS( );

    // coverage.h
    void regCoverage( );

    // replicates.h
    template <typename v_t> v_t symmetry( v_t, v_t );

    // replicates.h
    bool setBinValues( );
    // replicates.h
    bool setDecayValues( );

    // replicates.h
    template <typename v_t> v_t getFlatValue( std::vector<v_t> );
    // replicates.h
    double getMixedValue( double, double );

    // replicates.h
    bool setFlatValues( );
    // replicates.h
    bool setFlatDecay( );
    // replicates.h
    bool setDecayCDS( );

    // replicates.h
    bool setInGroup( );

    // replicates.h
    bool setBetweenGroup( );

    // normalization.h
    bool doNotNormalize( );
    // normalization.h
    bool setDivided( );
    // normalization.h
    bool normalizeSize( size_t );
    // normalization.h
    bool normalizeBinominalTest( );
    // normalization.h
    bool normalizeIC( );
    // normalization.h
    bool normalizeCoolIC( );
    // normalization.h
    bool setNormalized( );
    // normalization.h
    bool setDistDepDecayRemoved( );

    struct IceData
    {
        std::array<std::vector<double>, 2> vSliceBias;
        std::array<std::vector<double>, 2> vSliceMargin;
        std::vector<double> vBiases;
    };

    // normalization.h
    size_t iceGetCount( IceData&, size_t, size_t, bool );
    // normalization.h
    void icePreFilter( IceData&, bool, size_t, size_t, bool );
    // normalization.h
    void iceFilter( IceData&, size_t, size_t );
    // normalization.h
    void iceTimesOuterProduct( IceData&, bool, size_t, size_t );
    // normalization.h
    void iceMarginalize( IceData&, bool, size_t, size_t );
    // normalization.h
    double iceNonZeroMarginVariance( IceData&, bool, double );
    // normalization.h
    double iceNonZeroMarginMean( IceData&, bool );
    // normalization.h
    void iceDivByMargin( IceData&, bool, double, size_t, size_t );
    // normalization.h
    void iceApplyBias( IceData&, bool, size_t, size_t );
    // normalization.h
    double iceMaxBias( IceData&, bool );
    // normalization.h
    void regNormalization( );

    // colors.h
    bool setColors( );
    // colors.h
    bool setAnnotationColors( );
    // colors.h
    bool setCombined( );
    // colors.h
    bool setScaled( );
    // colors.h
    bool setColored( );
    // colors.h
    double logScale( double, double );
    // colors.h
    double colorRange( double, double, double );
    // colors.h
    size_t colorIndex( double );
    // colors.h
    bool setHeatmapCDS( );
    // colors.h
    bool setHeatmapExport( );
    // colors.h
    bool setPalette( );

    // colors.h
    void regColors( );


    // annotation.h
    bool setActivateAnnotation( );
    // annotation.h
    bool setAnnotationValues( );
    // annotation.h
    bool setAnnotationCDS( );
    // annotation.h
    bool setActivateAnnotationCDS( );
    // annotation.h
    void regAnnotation( );

    void registerAll( )
    {
        regBinSize( );
        regCoords( );
        regReplicates( );
        regCoverage( );
        regNormalization( );
        regColors( );
        regAnnotation( );

#ifndef NDEBUG
        checkUnnecessaryDependencies( );
#endif
    }

  protected:
    virtual std::vector<std::array<double, 2>> normalizeBinominalTestTrampoline( std::vector<std::array<size_t, 2>>&,
                                                                                 std::vector<std::array<size_t, 2>>&,
                                                                                 size_t, double )
    {
        throw std::logic_error( "Function not implemented" );
    }

    virtual std::vector<double> normalizeCoolerTrampoline( std::vector<size_t>&, size_t )
    {
        throw std::logic_error( "Function not implemented" );
    }

    virtual std::vector<std::string> colorPalette( std::string, std::string, std::string )
    {
        throw std::logic_error( "Function not implemented" );
    }

    virtual void print( std::string )
    {
        throw std::logic_error( "Function not implemented" );
    }

  public:
    PartialQuarry( )
        : sPrefix( "" ),
          uiCurrTime( 1 ),
          vGraph( NodeNames::SIZE ),
          vGraphData( NodeNames::SIZE ),
          xSession( ),
          xIndices( )
    {
        registerAll( );
        ++uiCurrTime;
    }

    PartialQuarry( std::string sPrefix )
        : sPrefix( sPrefix ),
          uiCurrTime( 1 ),
          vGraph( NodeNames::SIZE ),
          vGraphData( NodeNames::SIZE ),
          xSession( json::parse( std::ifstream( sPrefix + "/session.json" ) ) ),
          xIndices( sPrefix, false )
    {
        registerAll( );
        ++uiCurrTime;
    }

    virtual ~PartialQuarry( )
    {}

    // colors.h
    const pybind11::list getPalette( );

    // colors.h
    const std::string& getBackgroundColor( );

    // coords.h
    const std::vector<std::array<BinCoord, 2>>& getBinCoords( );

    // coords.h
    const pybind11::list getAnnotationList( bool );

    // colors.h
    const pybind11::dict getTicks( bool );

    // colors.h
    const pybind11::list getTickList( bool );

    // coords.h
    const std::vector<AxisCoord>& getAxisCoords( bool );

    // annotation.h
    const pybind11::dict getAnnotation( bool );

    // annotation.h
    const pybind11::list getDisplayedAnnos( bool );

    // colors.h
    const pybind11::dict getHeatmap( );

    // colors.h
    const std::array<double, 4> getPaletteTicks( );

    // colors.h
    const decltype( vHeatmapExport ) getHeatmapExport( );

    // bin_size.h
    const std::array<int64_t, 4> getDrawingArea( );

    // bin_size.h
    const std::array<size_t, 2> getBinSize( );

    // coords.h
    const std::array<size_t, 2> getCanvasSize( );

    // coverage.h
    const pybind11::dict getTracks( bool );

    // coverage.h
    const pybind11::dict getRankedSlices( bool );

    // replicates.h
    const pybind11::dict getDecayCDS( );

    // coverage.h
    const decltype( vTrackExport[ 0 ] ) getTrackExport( bool );

    // coverage.h
    const std::vector<std::string> getTrackExportNames( bool );

    // coverage.h
    const std::array<double, 2> getMinMaxTracks( bool bAxis );

  private:
    size_t stringCompare( std::string sQuery, std::string sRef )
    {
        if( sQuery.size( ) < 3 )
            return sQuery == sRef ? sQuery.size( ) : 0;

        size_t uiRet = 0;
        size_t uiPos = 0;

        for( char c : sQuery )
        {
            size_t uiCurr = uiPos;
            while( uiCurr < sRef.size( ) && std::tolower( c ) != std::tolower( sRef[ uiCurr ] ) )
                ++uiCurr;
            if( uiCurr < sRef.size( ) )
            {
                uiPos = uiCurr;
                ++uiRet;
            }
        }

        if( uiRet < sQuery.size( ) / 3 )
            return 0;

        return uiRet;
    }

    std::string substringChr( std::string sChr )
    {
        if( sChr.size( ) > uiLogestCommonSuffix + 2 )
            return sChr.substr( 0, sChr.size( ) - uiLogestCommonSuffix );
        else
            return sChr;
    }

  public:
    const pybind11::int_ interpretName( std::string sName, std::vector<bool> vbXAxis, bool bBottom )
    {
        size_t uiMaxEquality = 0;
        size_t uiMinRefSize = 0;
        size_t uiMaxPos = 0;
        const bool bSqueeze =
            this->xSession[ "settings" ][ "filters" ][ "anno_in_multiple_bins" ].get<std::string>( ) == "squeeze";
        for( bool bX : vbXAxis )
        {
            std::string sCoords =
                this->xSession[ "contigs" ][ bX ? "column_coordinates" : "row_coordinates" ].get<std::string>( );
            const bool bFullGenome = sCoords == "full_genome";
            size_t uiRunningPos = 0;
            for( auto xChr : vActiveChromosomes[ bX ? 0 : 1 ] )
            {
                std::string sSearch = "chromosome=" + xChr.sName;
                size_t uiEq = stringCompare( sName, sSearch );
                size_t uiLen;
                if( bFullGenome )
                    uiLen = xChr.uiLength;
                else
                {
                    int64_t iDataSetId =
                        this->xSession[ "annotation" ][ "by_name" ][ sCoords ][ xChr.sName ].get<int64_t>( );
                    if( bSqueeze )
                        uiLen = xIndices.vAnno.numIntervals( iDataSetId );
                    else
                        uiLen = xIndices.vAnno.totalIntervalSize( iDataSetId );
                }
                if( uiEq > uiMaxEquality || ( uiEq == uiMaxEquality && xChr.sName.size( ) < uiMinRefSize ) )
                {
                    uiMaxEquality = uiEq;
                    if( bBottom )
                        uiMaxPos = uiRunningPos;
                    else
                        uiMaxPos = uiRunningPos + uiLen;
                    uiMinRefSize = xChr.sName.size( );
                }
                uiRunningPos += uiLen;
            }

            if( ( bBottom && sName == "*" ) || sName.size( ) == 0 || sName == "start" )
                return pybind11::int_( 0 );
            if( ( !bBottom && sName == "*" ) || sName.size( ) == 0 || sName == "end" )
                return pybind11::int_( uiRunningPos );

            for( std::string sAnno : vActiveAnnotation[ bX ? 0 : 1 ] )
            {
                auto& rJson = this->xSession[ "annotation" ][ "by_name" ][ sAnno ];
                for( auto xChr : vActiveChromosomes[ bX ? 0 : 1 ] )
                {
                    int64_t iDataSetId = rJson[ xChr.sName ].get<int64_t>( );
                    xIndices.vAnno.iterate(
                        iDataSetId,
                        [ & ]( std::tuple<size_t, size_t, std::string, bool> xTup ) {
                            std::string sSearch = sAnno + "=" + std::get<2>( xTup );
                            size_t uiEq = stringCompare( sName, sSearch );
                            if( uiEq > uiMaxEquality ||
                                ( uiEq == uiMaxEquality && std::get<2>( xTup ).size( ) < uiMinRefSize ) )
                            {
                                uiMaxEquality = uiEq;
                                uiMinRefSize = std::get<2>( xTup ).size( );
                                if( bBottom )
                                    uiMaxPos = std::get<0>( xTup );
                                else
                                    uiMaxPos = std::get<1>( xTup );
                            }

                            return true;
                        },
                        !bFullGenome && !bSqueeze, !bFullGenome && bSqueeze );
                }
            }
        }
        return pybind11::int_( uiMaxPos );
    }

    void printSizes( )
    {
        std::cout << "vBinCoords " << vBinCoords.size( ) << std::endl;
        std::cout << "vvBinValues " << vvBinValues.size( ) << std::endl;
        std::cout << "vvFlatValues " << vvFlatValues.size( ) << std::endl;
        std::cout << "vvNormalized " << vvNormalized.size( ) << std::endl;
        std::cout << "vCombined " << vCombined.size( ) << std::endl;
        std::cout << "vColored " << vColored.size( ) << std::endl;
    }

    void printNumReadsNumOverlays()
    {
        std::cout << "replicate\tcontig_a\tcontig_b\treads\toverlays" << std::endl; 
        for(const json& rJson : this->xSession["replicates"]["list"])
        {
            std::string sRep = rJson.get<std::string>();
            for(const json& rJson : this->xSession["contigs"]["list"])
            {
                std::string sContigA = rJson.get<std::string>();
                size_t uiSizeA = this->xSession["contigs"]["lengths"][sContigA].get<size_t>();
                for(const json& rJson : this->xSession["contigs"]["list"])
                {
                    std::string sContigB = rJson.get<std::string>();
                    size_t uiSizeB = this->xSession["contigs"]["lengths"][sContigB].get<size_t>();
                    std::cout << sRep << "\t" << sContigA << "\t" << sContigB << "\t";

                    size_t iDataSetId = this->xSession["replicates"]["by_name"][sRep]["ids"][sContigA][sContigB].get<size_t>();
                    bool bHasMapQ = getValue<bool>( { "replicates", "by_name", sRep, "has_map_q" } );
                    bool bHasMultiMap = getValue<bool>( { "replicates", "by_name", sRep, "has_multimapping" } );

                    if( bHasMapQ && bHasMultiMap )
                    {
                        std::cout << xIndices.getIndex<3, 2>( )->count(
                            iDataSetId,
                            { 0, 0, 0 },
                            { uiSizeA, uiSizeB, 256 },
                            sps::IntersectionType::overlaps,
                            0 ) << "\t" << xIndices.getIndex<3, 2>( )->getNumOverlays(iDataSetId);
                    }
                    else if( !bHasMapQ && bHasMultiMap )
                    {
                        std::cout << xIndices.getIndex<2, 2>( )->count(
                            iDataSetId,
                            { 0, 0 },
                            { uiSizeA, uiSizeB },
                            sps::IntersectionType::overlaps,
                            0 ) << "\t" << xIndices.getIndex<2, 2>( )->getNumOverlays(iDataSetId);
                    }
                    else if( bHasMapQ && !bHasMultiMap )
                    {
                        std::cout << xIndices.getIndex<3, 0>( )->count(
                            iDataSetId,
                            { 0, 0, 0 },
                            { uiSizeA, uiSizeB, 256 },
                            sps::IntersectionType::overlaps,
                            0 ) << "\t" << xIndices.getIndex<3, 0>( )->getNumOverlays(iDataSetId);
                    }
                    else // if(!bHasMapQ && !bHasMultiMap)
                    {
                        std::cout << xIndices.getIndex<2, 0>( )->count(
                            iDataSetId,
                            { 0, 0 },
                            { uiSizeA, uiSizeB },
                            sps::IntersectionType::overlaps,
                            0 ) << "\t" << xIndices.getIndex<2, 0>( )->getNumOverlays(iDataSetId);
                    }

                    std::cout << std::endl;
                }
            }
        }
    }

    std::string getDOT( )
    {
        std::string sRet = "digraph libContactMappingFlowDiagram {\n";
        for( size_t uiI = 0; uiI < NodeNames::SIZE; uiI++ )
        {
            const ComputeNode& rNode = vGraph[ uiI ];
            const ComputeNodeData& rNodeData = vGraphData[ uiI ];
            std::string sJoined = "";
            for( std::vector<std::string> vParam : rNode.vIncomingSession )
            {
                for( std::string rS : vParam )
                {
                    if( sJoined.size( ) > 0 && sJoined.back( ) != '>' )
                        sJoined += ".";
                    sJoined += rS;
                }
                sJoined += "<br/>";
            }
            if( sJoined.size( ) > 0 )
            {
                sRet += "\t" + rNode.sNodeName + "_in [shape=box, label=<" + sJoined + ">];\n";
                sRet += "\t" + rNode.sNodeName + "_in -> " + rNode.sNodeName + ";\n";
            }
            if( rNodeData.bTerminal )
                sRet += "\t" + rNode.sNodeName + " [shape=hexagon];\n";
            for( const NodeNames& rIncoming : rNode.vIncomingFunctions )
                sRet += "\t" + vGraph[ rIncoming ].sNodeName + " -> " + rNode.sNodeName + ";\n";
        }
        return sRet + "}";
    }

    void saveSession( )
    {
        std::ofstream o( sPrefix + +"/session.json" );
        o << this->getSession( ) << std::endl;
    }
};

template <> pybind11::object PartialQuarry::getValue( std::vector<std::string> vKeys )
{
    return pyjson::from_json( getValue<json>( vKeys ) );
}

} // namespace cm
