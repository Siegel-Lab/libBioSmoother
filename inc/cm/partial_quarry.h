#include "cm/sps_interface.h"
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
    size_t uiSize;
};

struct BinCoord
{
    std::string sChromosomeX, sChromosomeY;
    size_t uiScreenX, uiScreenY;
    size_t uiIndexX, uiIndexY;
    size_t uiW, uiH;
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
        ActivateAnnotation,
        AnnotationValues,
        AnnotationCDS,
        AnnotationColors,
        ActivateAnnotationCDS,
        HeatmapCDS,
        Scaled,
        Ticks,
        Tracks,
        Divided,
        Palette,
        SIZE
    };
    struct ComputeNode
    {
        std::string sNodeName;
        bool ( PartialQuarry::*fFunc )( void );
        std::vector<NodeNames> vIncomingFunctions;
        std::vector<std::vector<std::string>> vIncomingSession;
        size_t uiLastUpdated;
    };

    std::string sPrefix;
    size_t uiCurrTime;
    std::vector<ComputeNode> vGraph;

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
        for( std::vector<std::string>& rKey : xNode.vIncomingSession )
            xSessionTime[ rKey ] = uiCurrTime + 1;
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


    size_t update_helper( NodeNames xNodeName )
    {
        ComputeNode& xNode = vGraph[ xNodeName ];
        if( xNode.uiLastUpdated < uiCurrTime )
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
                else if( uiUpdateTime > xNode.uiLastUpdated )
                    vOutOfDateReason.push_back( "change in previous node " + vGraph[ xPred ].sNodeName );
            }
            auto p2 = std::chrono::high_resolution_clock::now( );
            auto pms_int = std::chrono::duration_cast<std::chrono::milliseconds>( p2 - p1 );

            for( std::vector<std::string> vSetting : xNode.vIncomingSession )
                if( xSessionTime[ vSetting ] > xNode.uiLastUpdated )
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
                    std::cout << "last updated: " << xNode.uiLastUpdated << "; now: " << uiCurrTime << std::endl;
                }
                auto t1 = std::chrono::high_resolution_clock::now( );

                bool bUpdateDone = ( this->*xNode.fFunc )( );
                if( !bUpdateDone )
                {
                    if( uiVerbosity >= 1 )
                        std::cout << xNode.sNodeName << " was cancelled" << std::endl;
                    // update was cancelled -> we do not know how messed up the result form the last call was
                    xNode.uiLastUpdated = 0;
                    return 0;
                }

                auto t2 = std::chrono::high_resolution_clock::now( );
                auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 );
                if( uiVerbosity >= 1 )
                    std::cout << ms_int.count( ) << " ms (" << xNode.sNodeName << ")" << std::endl;
                if( uiVerbosity >= 2 )
                    std::cout << pms_int.count( ) << " ms (predecessors)" << std::endl << std::endl;

                xNode.uiLastUpdated = uiCurrTime;
            }
            else if( uiVerbosity >= 3 )
                std::cout << "no reason found to run node " << xNode.sNodeName << std::endl;
        }
        else if( uiVerbosity >= 3 )
            std::cout << "node up to date from previous check " << xNode.sNodeName << std::endl;

        size_t uiRet = xNode.uiLastUpdated;
        // node must be up to date now
        xNode.uiLastUpdated = uiCurrTime;

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
        for( size_t uiNodeName = 0; uiNodeName < NodeNames::SIZE; uiNodeName++ )
        {
            std::map<NodeNames, std::string> xDependents;
            std::map<std::vector<std::string>, std::string> xDependentSettings;
            for( NodeNames& rPrev : vGraph[ uiNodeName ].vIncomingFunctions )
                checkUnnecessaryDependency( rPrev, xDependents, xDependentSettings );

            for( std::vector<std::string>& rSetting : vGraph[ uiNodeName ].vIncomingSession )
                if( xDependentSettings.count( rSetting ) > 0 )
                    std::cerr << "unnecessary session dependency warning: " << vGraph[ uiNodeName ].sNodeName
                              << " depends on " << toString( rSetting ) << ", but the previous node "
                              << xDependentSettings[ rSetting ] << " shares this dependency" << std::endl;

            for( NodeNames& rPrev : vGraph[ uiNodeName ].vIncomingFunctions )
                if( xDependents.count( rPrev ) > 0 )
                    std::cerr << "unnecessary session dependency warning: " << vGraph[ uiNodeName ].sNodeName
                              << " depends on " << vGraph[ rPrev ].sNodeName << ", but the previous node "
                              << xDependents[ rPrev ] << " shares this dependency" << std::endl;
        }
    }

  public:
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
            for( auto& rNode : { HeatmapCDS, Tracks, AnnotationCDS, ActivateAnnotationCDS, Ticks, Tracks, Palette } )
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
        // @todo the session should be guarded by a mutex since there can now be multithreading in here
        // also changing a value while rendering should automatically call the cancel function?
        auto xRet = this->xSession[ toPointer( vKeys ) ].get<T>( );

        return xRet;
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

    std::vector<std::string> vActiveReplicates;
    std::array<std::vector<size_t>, 2> vInGroup;

    sps::IntersectionType xIntersect;

    std::vector<std::array<BinCoord, 2>> vBinCoords;

    std::array<std::vector<std::pair<std::string, bool>>, 2> vActiveCoverage;

    size_t uiSymmetry;
    size_t iInGroupSetting, iBetweenGroupSetting;

    std::vector<std::vector<size_t>> vvBinValues;
    std::vector<std::array<size_t, 2>> vvFlatValues;
    std::vector<std::array<double, 2>> vvNormalized;
    std::vector<double> vCombined;
    std::vector<double> vDivided;
    std::vector<double> vScaled;
    std::vector<std::string> vColored;

    std::string sBackgroundColor;

    std::array<std::vector<std::vector<size_t>>, 2> vvCoverageValues;
    std::array<std::array<std::vector<size_t>, 3>, 2> vInGroupCoverage;
    std::array<std::vector<size_t>, 2> vNormCoverage;
    // outer array: column / row
    std::array<std::vector<double>, 2> vvFlatCoverageValues;
    std::vector<std::array<size_t, 2>> vFlatNormValues;

    std::vector<std::string> vColorPalette;
    pybind11::list vRenderedPalette;
    std::vector<std::string> vColorPaletteAnnotation;
    std::array<std::vector<std::string>, 2> vActiveAnnotation;
    std::array<std::vector<std::vector<size_t>>, 2> vAnnotationValues;
    std::array<pybind11::dict, 2> vAnnotationCDS;
    std::array<pybind11::list, 2> vActiveAnnotationCDS;
    pybind11::dict xHeatmapCDS;
    std::array<pybind11::dict, 2> xTicksCDS;
    std::array<pybind11::dict, 2> xTracksCDS;
    std::array<pybind11::list, 2> vTickLists;
    std::array<size_t, 2> vCanvasSize;
    std::array<std::array<int64_t, 2>, 2> vvMinMaxTracks;
    double fMax, fMin;

    bool bCancel = false;
    std::mutex xUpdateMutex{ };


    // bin_size.h
    size_t nextEvenNumber( double );


    // bin_size.h
    bool setBinSize( );
    // bin_size.h
    bool setRenderArea( );

    // bin_size.h
    void regBinSize( );

    // coords.h
    bool setActiveChrom( );
    // coords.h
    bool setAxisCoords( );

    // coords.h
    bool setSymmetry( );
    // coords.h
    bool setBinCoords( );
    // coords.h
    bool setTicks( );

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
    bool setTracks( );

    // coverage.h
    bool setFlatCoverageValues( );

    // coverage.h
    void regCoverage( );

    // replicates.h
    size_t symmetry( size_t, size_t );

    // replicates.h
    bool setBinValues( );

    // replicates.h
    size_t getFlatValue( std::vector<size_t> );
    // replicates.h
    double getMixedValue( double, double );

    // replicates.h
    bool setFlatValues( );

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
    bool setNormalized( );

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

    virtual std::vector<std::string> colorPalette( std::string, std::string, std::string )
    {
        throw std::logic_error( "Function not implemented" );
    }

  public:
    PartialQuarry( ) : sPrefix( "" ), uiCurrTime( 1 ), vGraph( NodeNames::SIZE ), xSession( ), xIndices( )
    {
        registerAll( );
        ++uiCurrTime;
    }

    PartialQuarry( std::string sPrefix )
        : sPrefix( sPrefix ),
          uiCurrTime( 1 ),
          vGraph( NodeNames::SIZE ),
          xSession( json::parse( std::ifstream( sPrefix + "/session.json" ) ) ),
          xIndices( sPrefix )
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

    // colors.h
    const pybind11::dict getTicks( bool bXAxis );

    // colors.h
    const pybind11::list getTickList( bool bXAxis );

    // coords.h
    const std::vector<AxisCoord>& getAxisCoords( bool bXAxis );

    // annotation.h
    const pybind11::dict getAnnotation( bool bXAxis );

    // annotation.h
    const pybind11::list getDisplayedAnnos( bool bXAxis );

    // colors.h
    const pybind11::dict getHeatmap( );

    // bin_size.h
    const std::array<int64_t, 4> getDrawingArea( );

    // coords.h
    const std::array<size_t, 2> getCanvasSize( );

    // coverage.h
    const pybind11::dict getTracks( bool bXAxis );

    // coverage.h
    const std::array<int64_t, 2> getMinMaxTracks( bool bAxis );

    void printSizes( )
    {
        std::cout << "vBinCoords " << vBinCoords.size( ) << std::endl;
        std::cout << "vvBinValues " << vvBinValues.size( ) << std::endl;
        std::cout << "vvFlatValues " << vvFlatValues.size( ) << std::endl;
        std::cout << "vvNormalized " << vvNormalized.size( ) << std::endl;
        std::cout << "vCombined " << vCombined.size( ) << std::endl;
        std::cout << "vColored " << vColored.size( ) << std::endl;
    }

    std::string getDOT( )
    {
        std::string sRet = "digraph libContactMappingFlowDiagram {\n";
        for( const ComputeNode& rNode : vGraph )
        {
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
