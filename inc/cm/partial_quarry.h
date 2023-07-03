#include "cm/sps_interface.h"
#include <cctype>
#include <memory>
#include <mutex>

#pragma once


#define USE_GRID_QUERIES 0

#define CANCEL_RETURN                                                                                                  \
    if( this->bCancel )                                                                                                \
        return false;
#if 0 // this causes segmentation faults :(
    if( this->bAllowCtrlCCancel && PyErr_CheckSignals( ) != 0 ) /* allow Ctrl-C cancel */                                                         \
        throw pybind11::error_already_set( )
#endif

#define END_RETURN return true

namespace cm
{

struct ChromDesc
{
    std::string sName;
    size_t uiUnadjustedLength;
    size_t uiLength;
    size_t uiId;
};

struct IndexCoord
{
    size_t uiChromosome;
    size_t uiIndexPos;
    size_t uiIndexSize;
};

struct AnnoCoord : IndexCoord
{
    size_t uiAnnoId;
};

struct AxisCoord : IndexCoord
{
    size_t uiScreenPos;
    size_t uiScreenSize;
    size_t uiRegionIdx;
    size_t uiIdx;
};

struct DecayCoord
{
    size_t uiChromosomeX;
    size_t uiChromosomeY;
    int64_t iFrom;
    int64_t iTo;

    friend bool operator<( const DecayCoord& l, const DecayCoord& r )
    {
        if( l.uiChromosomeX != r.uiChromosomeX )
            return l.uiChromosomeX < r.uiChromosomeX;
        if( l.uiChromosomeY != r.uiChromosomeY )
            return l.uiChromosomeY < r.uiChromosomeY;
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

struct BinCoordBase
{
    size_t uiChromosomeX = std::numeric_limits<size_t>::max( );
    size_t uiChromosomeY = std::numeric_limits<size_t>::max( );
    size_t uiScreenX, uiScreenY;
    size_t uiIndexX, uiIndexY;
    size_t uiScreenW, uiScreenH;
    size_t uiIndexW, uiIndexH;
    size_t uiXAxisIdx, uiYAxisIdx;
};

struct BinCoord : BinCoordBase
{
    size_t uiDecayCoordIndex = std::numeric_limits<size_t>::max( );
};

struct BinCoordRegion : BinCoordBase
{
    size_t uiCoordStartIdxX = 0, uiCoordStartIdxY = 0;
    size_t uiNumCoordsX = 0, uiNumCoordsY = 0;
};

struct Annotation
{
    std::string sInfo;
    bool bForw;
    double fScreenX, fScreenY;
    size_t uiIndexX, uiIndexY;
    size_t uiRow;
    size_t uiChromosome;
};


/*
 * iterates over the indices between uiFrom and uiTo in the following order:
 * every index is assigned a priority:
 *   - 1st: the maximal power of two the index is evenly devidable by
 *   - 2nd: the index intself
 * indices are iterates in the order of their 1st priority (highest to lowest)
 * if two indices have the same 1st priority they are iterated in order of their 2nd priority (lowest to highest)
 */
void iterateEvenlyDividableByMaxPowerOfTwo( size_t uiFrom, size_t uiTo, std::function<bool( size_t )> fDo )
{
    if( uiFrom + 1 >= uiTo )
    {
        fDo( uiFrom );
        return;
    }

    size_t uiN = std::log2( uiTo - uiFrom ) + 2;
    while( uiN > 0 )
    {
        --uiN;
        size_t uiX = 1 << uiN;
        if( uiFrom % uiX == 0 && uiFrom % ( uiX * 2 ) != 0 )
        {
            if( !fDo( uiFrom ) )
                return;
        }
        size_t uiY = uiX * ( uiFrom / uiX ) + uiX;
        while( uiY < uiTo )
        {
            assert( uiY >= uiFrom );
            if( uiY % ( uiX * 2 ) != 0 )
            {
                if( !fDo( uiY ) )
                    return;
                uiY += uiX * 2;
            }
            else
                uiY += uiX;
        }
        assert( uiY + uiX >= uiTo );
    }
    if( uiFrom == 0 )
        fDo( 0 );
}

size_t getEvenlyDividableByMaxTwoPowNIn( size_t uiFrom, size_t uiTo )
{
    size_t uiN;
    iterateEvenlyDividableByMaxPowerOfTwo( uiFrom, uiTo, [ & ]( size_t uiX ) {
        uiN = uiX;
        return false;
    } );
    return uiN;
}


class PartialQuarry : public HasSession
{
  public:
    size_t uiVerbosity = 1;
    bool bAllowCtrlCCancel = false;

  private:
    enum NodeNames
    {
        BinSize = 0,
        RenderArea,
        ActiveChrom,
        AxisCoords,
        Symmetry,
        BinCoords,
        DecayCoords,
        FlatDecay,
        IntersectionType,
        ActiveReplicates,
        ActiveCoverage,
        CoverageValues,
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
        Palette,
        AnnoFilters,
        LCS,
        CanvasSize,
        MappingQuality,
        Directionality,
        RankedSlicesCDS,
        GridSeqCoverage,
        RadiclSeqCoverage,
        RnaAssociatedGenesFilter,
        RnaAssociatedBackground,
        GridSeqSamples,
        RadiclSeqSamples,
        DatasetIdPerRepl,
        ActiveChromLength,
        V4cCoords,
        Flat4C,
        IceCoords,
        SIZE
    };
    struct ComputeNode
    {
        std::string sNodeName;
        bool ( PartialQuarry::*fFunc )( void );
        std::vector<NodeNames> vIncomingFunctions;
        std::vector<std::vector<std::string>> vIncomingSession;
        std::vector<std::vector<std::string>> vSessionsIncomingInPrevious;
        bool bHidden = false;
    };

    struct ComputeNodeData
    {
        size_t uiLastUpdated = 0;
        size_t uiLastChecked = 0;
        std::string sError = "";
        bool bRegistered = false;
        bool bTerminal = true;
        int64_t uiLastRuntime = 0;
    };

    size_t uiCurrTime;
    std::vector<ComputeNode> vGraph;
    std::vector<ComputeNodeData> vGraphData;
    NodeNames xCurrNodeName = NodeNames::SIZE;

    std::map<std::vector<std::string>, size_t> xSessionTime;


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

    size_t update_helper( NodeNames xNodeName, const std::function<void( const std::string& )>& fPyPrint )
    {
        ComputeNode& xNode = vGraph[ xNodeName ];
        ComputeNodeData& xNodeData = vGraphData[ xNodeName ];
        if( xNodeData.uiLastChecked < uiCurrTime )
        {
            std::vector<std::string> vOutOfDateReason{ };
            auto p1 = std::chrono::high_resolution_clock::now( );
            for( NodeNames xPred : xNode.vIncomingFunctions )
            {
                size_t uiUpdateTime = update_helper( xPred, fPyPrint );
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
                    xNodeData.uiLastChecked = 0;
                    return 0;
                }

                auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 );
                if( uiVerbosity >= 1 )
                    fPyPrint( std::to_string( ms_int.count( ) ) + " ms (" + xNode.sNodeName + ")" );
                if( uiVerbosity >= 2 )
                    std::cout << pms_int.count( ) << " ms (predecessors)" << std::endl << std::endl;

                xNodeData.uiLastUpdated = uiCurrTime;
                auto mics_int = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 );
                xNodeData.uiLastRuntime = mics_int.count( );
            }
            else if( uiVerbosity >= 3 )
                std::cout << "no reason found to run node " << xNode.sNodeName << std::endl;
        }
        else if( uiVerbosity >= 3 )
            std::cout << "node up to date from previous check " << xNode.sNodeName << std::endl;

        // node must be up to date now
        xNodeData.uiLastChecked = uiCurrTime;

        if( uiVerbosity >= 4 )
            std::cout << "check time of node " << xNode.sNodeName << " to " << uiCurrTime << std::endl;

        return xNodeData.uiLastUpdated;
    }


    bool update_no_throw( NodeNames xNodeName, const std::function<void( const std::string& )>& fPyPrint )
    {
        // allow python multithreading during this call
        pybind11::gil_scoped_release release;

        // make sure that no two threads can enter here (basically allow multithreading other than this)
        std::lock_guard<std::mutex> xGuard( xUpdateMutex.xX );

        bCancel = false;
        size_t uiUpdateTime = update_helper( xNodeName, fPyPrint );

        return uiUpdateTime != 0;
    }

    void update( NodeNames xNodeName, const std::function<void( const std::string& )>& fPyPrint )
    {
        if( !update_no_throw( xNodeName, fPyPrint ) )
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
            if( vGraphData[ uiNodeName ].bRegistered && xNodeNames.count( vGraph[ uiNodeName ].sNodeName ) > 0 )
                std::cerr << "two nodes have been registered under the same name: " << vGraph[ uiNodeName ].sNodeName
                          << std::endl;
            if( vGraphData[ uiNodeName ].bRegistered )
                xNodeNames.insert( vGraph[ uiNodeName ].sNodeName );

            if( !vGraphData[ uiNodeName ].bRegistered )
            {
                std::cerr << "Node has not been registered. id: " << uiNodeName;
                if( uiNodeName > 0 )
                    std::cerr << " (the previous node is " << vGraph[ uiNodeName - 1 ].sNodeName << ")";
                std::cerr << std::endl;
            }
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

    void updateCDS( const std::function<void( const std::string& )>& fPyPrint )
    {
        bool bContinue;
        do
        {
            bContinue = false;
            for( auto& rNode : { HeatmapCDS, Tracks, AnnotationCDS, ActivateAnnotationCDS, Ticks, Palette, DecayCDS,
                                 RankedSlicesCDS } )
                if( !update_no_throw( rNode, fPyPrint ) )
                {
                    bContinue = true;
                    break;
                }
        } while( bContinue );
    }

    void clearCache( )
    {
        ++uiCurrTime;
        for( auto& rX : xSessionTime )
            rX.second = uiCurrTime;
    }

    void updateAll( const std::function<void( const std::string& )>& fPyPrint )
    {
        for( size_t uiNodeName = 0; uiNodeName < NodeNames::SIZE; uiNodeName++ )
            update_helper( (NodeNames)uiNodeName, fPyPrint );
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
        try
        {
            return this->xSession[ toPointer( vKeys ) ].get<T>( );
        }
        catch( ... )
        {
            throw std::runtime_error( std::string( "error querying json value: " ) + toString( vKeys ) );
        }
    }

    bool hasValue( std::vector<std::string> vKeys )
    {
        std::string sKey = vKeys.back( );
        vKeys.pop_back( );
        return getValue<json>( vKeys ).contains( sKey );
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
#ifndef NDEBUG
            std::cout << "session updated: " << toString( vKeys ) << std::endl;
#endif
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

    using interface_t = SpsInterface<false>;
    std::shared_ptr<interface_t> pIndices;

    using coordinate_t = interface_t::coordinate_t;

    const std::string getFilePrefix( )
    {
        return pIndices->getFilePrefix( );
    }

  private:
    /**
     * @brief Dependency Declared Access
     *
     * will throw warnings in DEBUG mode if accessed without declaration in the comp graph.
     *
     * @tparam T
     */
    template <typename T> class DepDec
    {
      private:
        T xData;
#ifndef NDEBUG
        PartialQuarry* rQ;
        std::string sName;
        std::set<NodeNames> vWriteAccess;
        std::set<NodeNames> vReadAccess;
#endif

      public:
        DepDec( [[maybe_unused]] PartialQuarry* rQ, [[maybe_unused]] const std::string& sName )
#ifndef NDEBUG
            : rQ( rQ ), sName( sName )
#endif
        {}


        const T& r( ) const
        {
#ifndef NDEBUG
            if( rQ->xCurrNodeName != NodeNames::SIZE && vReadAccess.count( rQ->xCurrNodeName ) == 0 &&
                vWriteAccess.count( rQ->xCurrNodeName ) == 0 )
                std::cout << "undeclared data dependency warning: " << rQ->vGraph[ rQ->xCurrNodeName ].sNodeName
                          << " reads from " << sName << " but does not declare this." << std::endl;
#endif
            return xData;
        }

        T& w( )
        {
#ifndef NDEBUG
            if( rQ->xCurrNodeName != NodeNames::SIZE && vWriteAccess.count( rQ->xCurrNodeName ) == 0 )
                std::cout << "undeclared data dependency warning: " << rQ->vGraph[ rQ->xCurrNodeName ].sNodeName
                          << " writes to " << sName << " but does not declare this." << std::endl;
#endif
            return xData;
        }

        void registerWrite( [[maybe_unused]] const NodeNames& xNode )
        {
#ifndef NDEBUG
            vWriteAccess.insert( xNode );
#endif
        }

        void registerRead( [[maybe_unused]] const NodeNames& xNode )
        {
#ifndef NDEBUG
            vReadAccess.insert( xNode );
#endif
        }
    };

    DepDec<size_t> uiBinWidth, uiBinHeight;
    int64_t iStartX, iStartY, iEndX, iEndY;

    /* Coord systems are:
     * - normal window
     * - virtual 4 C on x axis
     * - virtual 4 C on y axis
     *
     */
    const static size_t NUM_COORD_SYSTEMS = 3;

    std::array<std::vector<ChromDesc>, 2> vActiveChromosomes;
    std::array<std::vector<AxisCoord>, 2> vAxisCords;
    std::array<std::vector<AxisCoord>, 2> vIceAxisCoords;
    std::array<std::vector<AxisCoord>, 2> vV4cCoords;
    std::array<std::vector<AxisCoord>, 2> vIceV4cCoords;
    std::array<std::vector<std::array<DecayCoord, 2>>, NUM_COORD_SYSTEMS> vDistDepDecCoords;
    std::array<std::vector<AxisRegion>, 2> vAxisRegions;
    std::array<std::vector<std::vector<size_t>>, 2> vvGridSeqCoverageValues;

    std::vector<std::string> vActiveReplicates;
    std::vector<size_t> vActiveReplicatesTotal;
    std::array<std::vector<size_t>, 2> vInGroup;

    sps::IntersectionType xIntersect;

    std::array<std::vector<std::array<BinCoord, 2>>, NUM_COORD_SYSTEMS> vBinCoordsIce;
    std::array<std::vector<std::array<BinCoord, 2>>, NUM_COORD_SYSTEMS> vBinCoords;

    std::array<std::vector<std::string>, 2> vActiveCoverage;

    size_t uiSymmetry;
    size_t iInGroupSetting, iBetweenGroupSetting;

    std::array<std::vector<std::vector<size_t>>, NUM_COORD_SYSTEMS> vvBinValues;
    std::array<std::vector<std::vector<double>>, NUM_COORD_SYSTEMS> vvDecayValues;
    std::array<std::vector<std::array<size_t, 2>>, NUM_COORD_SYSTEMS> vvFlatValues;
    std::array<std::vector<std::array<double, 2>>, NUM_COORD_SYSTEMS> vvFlatDecay;
    std::array<std::array<size_t, 2>, NUM_COORD_SYSTEMS> vvFlatTotal;
    std::array<std::vector<std::array<double, 2>>, NUM_COORD_SYSTEMS> vvNormalized;
    std::array<std::vector<std::array<double, 2>>, NUM_COORD_SYSTEMS> vvNormalizedDDD;
    std::array<std::vector<double>, 3> vCombined;
    std::array<std::vector<double>, 2> vFlat4C;
    std::vector<double> vDivided;
    std::vector<double> vScaled;
    std::vector<double> vRanged;
    std::vector<std::string> vColored;

    std::string sBackgroundColor;

    size_t uiFromAnnoFilter;
    size_t uiToAnnoFilter;
    size_t uiFromSameStrandFilter;
    size_t uiToSameStrandFilter;
    size_t ui1DFromStrandFilter;
    size_t ui1DToStrandFilter;
    size_t uiFromYStrandFilter;
    size_t uiToYStrandFilter;
    size_t uiMapQMin;
    size_t uiMapQMax;

    std::array<std::vector<std::vector<size_t>>, 2> vvCoverageValues;
    std::array<std::vector<size_t>, 2> vNormCoverage;

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
    std::array<pybind11::dict, 2> xContigTicksCDS;
    std::array<pybind11::dict, 2> xTracksCDS;
    std::array<pybind11::list, 2> vTickLists;
    std::array<pybind11::list, 2> vContigStartList;
    std::array<size_t, 2> vCanvasSize;
    std::array<std::array<double, 2>, 2> vvMinMaxTracks;
    double fMax, fMin;
    double fDataMax, fDataMin;
    size_t uiLogestCommonSuffix;

    std::vector<std::array<size_t, 2>> vGridSeqAnnoCoverage;
    std::vector<std::array<bool, 2>> vGridSeqFiltered;
    std::vector<std::pair<size_t, size_t>> vChromIdForAnnoIdx;
    std::vector<size_t> vBackgroundGridSeq;

    std::vector<AnnoCoord> vGridSeqSamples;
    std::vector<AnnoCoord> vRadiclSeqSamples;
    std::vector<std::array<size_t, 2>> vRadiclSeqCoverage;
    std::vector<std::array<size_t, 2>> vRadiclSeqNumNonEmptyBins;

    std::vector<size_t> vvDatasetIdsPerReplAndChr;
    size_t uiContigListSize;

    size_t uiIceFilterIgnoreDiags;
    std::array<std::array<std::array<std::vector<double>, 2>, 2>, NUM_COORD_SYSTEMS> vIceSliceBias;

    size_t getDatasetIdfromReplAndChr( size_t uiRepl, size_t uiChrX, size_t uiChrY )
    {
        return vvDatasetIdsPerReplAndChr[ uiRepl ] + uiChrY + uiChrX * uiContigListSize;
    }

    bool bCancel = false;

    template <typename T> class NoCopy
    {
      public:
        T xX;

        NoCopy& operator=( const NoCopy& )
        {
            // do not copy xX!!!!
            return *this;
        }
    };

    NoCopy<std::mutex> xUpdateMutex{ };


    // bin_size.h
    size_t nextEvenNumber( double );

    // coords.h
    size_t getAnnoListId( std::string );

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
    bool setActiveChromLength( );
    // coords.h
    bool setAnnoFilters( );
    // coords.h
    std::pair<std::vector<AxisCoord>, std::vector<AxisRegion>> setAxisCoordsHelper( bool );
    bool setAxisCoords( );

    // coords.h
    template <typename out_t, typename in_t> std::array<out_t, 2> binObjFromCoords( const in_t&, const in_t& );
    // coords.h
    std::array<BinCoordRegion, 2> binRegionObjFromCoords( const AxisRegion&, const AxisRegion& );
    // coords.h
    bool setSymmetry( );
    // coords.h
    bool setMappingQuality( );
    // coords.h
    bool setDirectionality( );
    // coords.h
    const std::vector<AxisCoord>& pickXCoords( const size_t );
    // coords.h
    const std::vector<AxisCoord>& pickYCoords( const size_t );
    // coords.h
    bool setBinCoords( );
    // coords.h
    bool setDecayCoords( );
    // coords.h
    bool sampleAndMerge( size_t, size_t, const std::vector<AxisCoord>&, std::vector<AxisCoord>& );
    // coords.h
    bool setIceCoords( );
    // coords.h
    bool setLCS( );
    // coords.h
    bool setTicks( );
    // coords.h
    bool setCanvasSize( );
    // coords.h
    bool setV4cCoords( );

    // coords.h
    void regCoords( );

    // replicates.h
    bool setIntersectionType( );

    // replicates.h
    bool setActiveReplicates( );

    // replicates.h
    bool setDatasetIdPerRepl( );

    // replicates.h
    void regReplicates( );

    // coverage.h
    bool setActiveCoverage( );

    // coverage.h
    bool setCoverageValues( );
    // coverage.h
    bool setTracks( );
    // coverage.h
    bool setRankedSlicesCDS( );
    // coverage.h
    bool setTrackExport( );
    // coverage.h
    size_t getMaxCoverageFromRepl( const AxisCoord&, const size_t, size_t, bool, bool );
    // coverage.h
    size_t getMaxCoverageFromRepl( const size_t, const size_t, const size_t, const size_t, size_t, bool, bool );
    // coverage.h
    size_t getCoverageFromRepl( const AxisCoord&, const size_t, bool, bool );
    // coverage.h
    size_t getCoverageFromRepl( const size_t, const size_t, const size_t, const size_t, bool, bool );
    // coverage.h
    std::tuple<size_t, int64_t, size_t, size_t> makeHeapTuple( bool, bool, const size_t, const size_t, int64_t, size_t,
                                                               size_t );

    // coverage.h
    void regCoverage( );

    // replicates.h
    bool hasChr( size_t, const std::string& );
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
    enum GetSamplesMode
    {
        OneAnnotation,
        Bins
        //, BinnedAnno
    };
    bool getSamples( const GetSamplesMode&, const size_t, const std::string&, const size_t, const bool, const bool,
                     std::vector<AnnoCoord>& );
    // normalization.h
    size_t getChromIdxForAnnoIdx( size_t );
    // normalization.h
    bool setRadiclSeqSamples( );
    // normalization.h
    bool setICESamples( );
    // normalization.h
    bool setGridSeqSamples( );
    // normalization.h
    bool setGridSeqCoverage( );
    // normalization.h
    bool setRnaAssociatedGenesFilter( );
    // normalization.h
    bool setRnaAssociatedBackground( );
    // normalization.h
    bool doNotNormalize( );
    // normalization.h
    bool normalizeSize( size_t );
    // normalization.h
    bool normalizeBinominalTest( );
    // normalization.h
    bool setRadiclSeqCoverage( );
    // normalization.h
    bool normalizeGridSeq( );
    // normalization.h
    bool normalizeIC( );
    bool normalizeSymmIC( );
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
    struct SymmIceData
    {
        std::vector<double> vSliceBias;
        std::vector<double> vSliceMargin;
        std::vector<double> vBiases;
    };

    // normalization.h
    size_t iceGetCount( IceData&, size_t, size_t, size_t, bool );
    size_t iceGetCount( SymmIceData&, size_t, size_t, size_t, bool );
    // normalization.h
    void icePreFilter( IceData&, bool, size_t, size_t, size_t, bool );
    // normalization.h
    void iceFilter( IceData&, size_t, size_t );
    void iceFilter( SymmIceData&, size_t, size_t );
    // normalization.h
    void iceTimesOuterProduct( IceData&, bool, size_t, size_t, size_t );
    void iceTimesOuterProduct( SymmIceData&, bool, size_t, size_t, size_t );
    // normalization.h
    void iceMarginalize( IceData&, bool, size_t, size_t );
    void iceMarginalize( SymmIceData&, size_t, size_t );
    // normalization.h
    double iceNonZeroMarginVariance( IceData&, bool, double );
    double iceNonZeroMarginVariance( SymmIceData&, double );
    // normalization.h
    double iceNonZeroMarginMean( IceData&, bool );
    double iceNonZeroMarginMean( SymmIceData& );
    // normalization.h
    void iceDivByMargin( IceData&, bool, double, size_t, size_t );
    void iceDivByMargin( SymmIceData&, double, size_t, size_t );
    // normalization.h
    void rescaleBias( IceData&, bool, double, size_t, size_t );
    void rescaleBias( SymmIceData&, double, size_t, size_t );
    // normalization.h
    double iceMaxBias( IceData&, bool );
    double iceMaxBias( SymmIceData& );
    // normalization.h
    void regNormalization( );

    // colors.h
    bool setColors( );
    // colors.h
    bool setAnnotationColors( );
    // colors.h
    bool setCombined( );
    // colors.h
    bool setFlat4C( );
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
    virtual std::vector<std::array<double, 2>>
    normalizeBinominalTestTrampoline( const std::vector<std::array<size_t, 2>>&,
                                      const std::vector<std::array<size_t, 2>>&,
                                      const std::vector<std::array<size_t, 2>>&, size_t, size_t, double, bool, size_t )
    {
        throw std::logic_error( "Function not implemented" );
    }

    virtual std::vector<double> normalizeCoolerTrampoline( const std::vector<size_t>&, size_t )
    {
        throw std::logic_error( "Function not implemented" );
    }

    virtual std::vector<std::string> colorPalette( std::string, std::string, std::string )
    {
        throw std::logic_error( "Function not implemented" );
    }

  public:
    PartialQuarry( std::shared_ptr<interface_t> pIndex )
        : HasSession( pIndex ),
          uiCurrTime( 1 ),
          vGraph( NodeNames::SIZE ),
          vGraphData( NodeNames::SIZE ),
          pIndices( pIndex ),
          uiBinWidth( this, "uiBinWidth" ),
          uiBinHeight( this, "uiBinHeight" )
    {
        registerAll( );
        ++uiCurrTime;

        // @todo-low-prio these should be declared while registering the nodes
        // also every member variable should be wrapped with DepDec.
        uiBinWidth.registerWrite( NodeNames::BinSize );
        uiBinWidth.registerRead( NodeNames::ActiveChromLength );
        uiBinWidth.registerRead( NodeNames::RenderArea );
        uiBinWidth.registerRead( NodeNames::AxisCoords );
        uiBinWidth.registerRead( NodeNames::IceCoords );
        uiBinWidth.registerRead( NodeNames::RadiclSeqSamples );
        uiBinWidth.registerRead( NodeNames::Normalized );

        uiBinHeight.registerWrite( NodeNames::BinSize );
        uiBinHeight.registerRead( NodeNames::ActiveChromLength );
        uiBinHeight.registerRead( NodeNames::RenderArea );
        uiBinHeight.registerRead( NodeNames::AxisCoords );
        uiBinHeight.registerRead( NodeNames::IceCoords );
        uiBinHeight.registerRead( NodeNames::RadiclSeqSamples );
        uiBinHeight.registerRead( NodeNames::Normalized );
    }

    PartialQuarry( ) : PartialQuarry( std::shared_ptr<interface_t>{ } )
    {}

    PartialQuarry( std::string sPrefix ) : PartialQuarry( std::make_shared<interface_t>( sPrefix ) )
    {}

    void copyFrom( const PartialQuarry& rOther )
    {
        *this = rOther;
    }

    virtual ~PartialQuarry( )
    {}

    // colors.h
    const pybind11::list getPalette( const std::function<void( const std::string& )>& );

    // colors.h
    const std::string& getBackgroundColor( const std::function<void( const std::string& )>& );

    // coords.h
    const std::vector<std::array<BinCoord, 2>>& getBinCoords( const std::function<void( const std::string& )>& );

    // coords.h
    const pybind11::list getAnnotationList( bool, const std::function<void( const std::string& )>& );

    // colors.h
    const pybind11::dict getTicks( bool, const std::function<void( const std::string& )>& );
    const pybind11::dict getContigTicks( bool, const std::function<void( const std::string& )>& );

    // colors.h
    const pybind11::list getTickList( bool, const std::function<void( const std::string& )>& );
    // colors.h
    const pybind11::list getContigStartList( bool, const std::function<void( const std::string& )>& );

    // coords.h
    const std::vector<AxisCoord>& getAxisCoords( bool, const std::function<void( const std::string& )>& );

    // coords.h
    size_t getAxisSize( bool, const std::function<void( const std::string& )>& );

    // annotation.h
    const pybind11::dict getAnnotation( bool, const std::function<void( const std::string& )>& );

    // annotation.h
    const pybind11::list getDisplayedAnnos( bool, const std::function<void( const std::string& )>& );

    // colors.h
    const pybind11::dict getHeatmap( const std::function<void( const std::string& )>& );

    // colors.h
    const std::array<double, 4> getPaletteTicks( const std::function<void( const std::string& )>& );

    // colors.h
    const decltype( vHeatmapExport ) getHeatmapExport( const std::function<void( const std::string& )>& );

    // normalization.h
    const decltype( vDivided ) getDivided( const std::function<void( const std::string& )>& );

    // normalization.h
    const decltype( vScaled ) getScaled( const std::function<void( const std::string& )>& );

    // colors.h
    const std::vector<double> getCombined( const std::function<void( const std::string& )>& );

    // bin_size.h
    const std::array<int64_t, 4> getDrawingArea( const std::function<void( const std::string& )>& );

    // bin_size.h
    const std::array<size_t, 2> getBinSize( const std::function<void( const std::string& )>& );

    // coords.h
    const std::array<size_t, 2> getCanvasSize( const std::function<void( const std::string& )>& );

    // coverage.h
    const pybind11::dict getTracks( bool, const std::function<void( const std::string& )>& );

    // replicates.h
    const pybind11::dict getDecayCDS( const std::function<void( const std::string& )>& );

    // coverage.h
    const decltype( vTrackExport[ 0 ] ) getTrackExport( bool, const std::function<void( const std::string& )>& );

    // coverage.h
    const std::vector<std::string> getTrackExportNames( bool, const std::function<void( const std::string& )>& );

    // coverage.h
    const std::array<double, 2> getMinMaxTracks( bool, const std::function<void( const std::string& )>& );

    // coverage.h
    const pybind11::dict getRankedSlices( bool, const std::function<void( const std::string& )>& );

    size_t getLongestCommonSuffix( const std::function<void( const std::string& )>& fPyPrint )
    {
        update( NodeNames::LCS, fPyPrint );
        return uiLogestCommonSuffix;
    }

  private:
    size_t stringCompare( std::string sQuery, std::string sRef )
    {
        if( sQuery.size( ) < 3 )
            return sQuery == sRef ? sQuery.size( ) : 0;

        size_t uiPos = 0;

        for( char c : sRef )
            if( uiPos < sQuery.size( ) && std::tolower( c ) == std::tolower( sQuery[ uiPos ] ) )
                ++uiPos;

        if( uiPos < sQuery.size( ) / 3 )
            return 0;

        return uiPos;
    }

    std::string substringChr( std::string sChr )
    {
        if( sChr.size( ) > uiLogestCommonSuffix + 2 )
            return sChr.substr( 0, sChr.size( ) - uiLogestCommonSuffix );
        else
            return sChr;
    }

    std::string to_lower( std::string s ) const
    {
        std::string ret;
        for( const char c : s )
            ret += std::tolower( c );
        return ret;
    }

  public:
    const pybind11::int_ interpretName( std::string sName, bool bXAxis, bool bBottom, bool bGenomicCoords )
    {
        sName = to_lower( sName );
        if( ( bBottom && sName == "*" ) || sName.size( ) == 0 || sName == "start" )
            return pybind11::int_( 0 );
        if( !bGenomicCoords && ( ( !bBottom && sName == "*" ) || sName.size( ) == 0 || sName == "end" ) )
            return pybind11::int_( vCanvasSize[ bXAxis ? 0 : 1 ] );
        else if( bGenomicCoords && ( ( !bBottom && sName == "*" ) || sName.size( ) == 0 || sName == "end" ) )
            return pybind11::int_( getValue<size_t>( { "contigs", "genome_size" } ) );

        size_t uiMaxPos = 0;
        const bool bSqueeze =
            this->xSession[ "settings" ][ "filters" ][ "anno_in_multiple_bins" ].get<std::string>( ) == "squeeze";
        std::string sCoords = getValue<std::string>( { "contigs", "annotation_coordinates" } );
        const bool bFullGenome =
            bGenomicCoords ||
            !getValue<bool>( { "settings", "filters", bXAxis ? "anno_coords_col" : "anno_coords_row" } );
        size_t uiRunningPos = 0;
        for( auto xChr : vActiveChromosomes[ bXAxis ? 0 : 1 ] )
        {
            bool bMatch = to_lower( xChr.sName ).find( sName ) != std::string::npos;
            size_t uiLen = 0;
            if( bFullGenome )
                uiLen = xChr.uiLength;
            else if( this->xSession[ "annotation" ][ "by_name" ][ sCoords ].contains( xChr.sName ) )
            {
                int64_t iDataSetId =
                    this->xSession[ "annotation" ][ "by_name" ][ sCoords ][ xChr.sName ].get<int64_t>( );
                if( bSqueeze )
                    uiLen = pIndices->vAnno.numIntervals( iDataSetId );
                else
                    uiLen = pIndices->vAnno.totalIntervalSize( iDataSetId );
            }
            if( bMatch )
            {
                if( bBottom )
                {
                    uiMaxPos = uiRunningPos;
                    break;
                }
                else
                    uiMaxPos = uiRunningPos + uiLen;
            }
            uiRunningPos += uiLen;
        }

        return pybind11::int_( uiMaxPos );
    }

    void printSizes( )
    {
        std::cout << "vBinCoords " << vBinCoords[ 0 ].size( ) << std::endl;
        std::cout << "vvBinValues " << vvBinValues[ 0 ].size( ) << std::endl;
        std::cout << "vvFlatValues " << vvFlatValues[ 0 ].size( ) << std::endl;
        std::cout << "vvNormalized " << vvNormalized[ 0 ].size( ) << std::endl;
        std::cout << "vCombined " << vCombined.size( ) << std::endl;
        std::cout << "vColored " << vColored.size( ) << std::endl;
    }

    void printNumReadsNumOverlays( )
    {
        std::cout << "replicate\tcontig_a\tcontig_b\treads\toverlays" << std::endl;
        for( const json& rJson : this->xSession[ "replicates" ][ "list" ] )
        {
            std::string sRep = rJson.get<std::string>( );
            for( const json& rJson : this->xSession[ "contigs" ][ "list" ] )
            {
                std::string sContigA = rJson.get<std::string>( );
                size_t uiSizeA = this->xSession[ "contigs" ][ "lengths" ][ sContigA ].get<size_t>( );
                for( const json& rJson : this->xSession[ "contigs" ][ "list" ] )
                {
                    std::string sContigB = rJson.get<std::string>( );
                    size_t uiSizeB = this->xSession[ "contigs" ][ "lengths" ][ sContigB ].get<size_t>( );
                    std::cout << sRep << "\t" << sContigA << "\t" << sContigB << "\t";
                    size_t uiNumAnno = this->xSession[ "annotation" ][ "list" ].size( );

                    size_t iDataSetId =
                        this->xSession[ "replicates" ][ "by_name" ][ sRep ][ "ids" ][ sContigA ][ sContigB ]
                            .get<size_t>( );

                    std::cout << pIndices->count( iDataSetId,
                                                  { 0, 0, 0, 0 },
                                                  { uiSizeA, uiSizeB, 256, uiNumAnno * 3 + 2 },
                                                  sps::IntersectionType::overlaps,
                                                  0 )
                              << "\t" << pIndices->getNumOverlays( iDataSetId );

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
            for( const NodeNames& rIncoming : rNode.vIncomingFunctions )
                if( vGraph[ rIncoming ].bHidden )
                    sJoined += "<i>[" + vGraph[ rIncoming ].sNodeName + "]</i><br/>";
            if( sJoined.size( ) > 0 )
            {
                sRet += "\t" + rNode.sNodeName + "_in [shape=box, label=<" + sJoined + ">];\n";
                sRet += "\t" + rNode.sNodeName + "_in -> " + rNode.sNodeName + ";\n";
            }
            if( rNodeData.bTerminal )
                sRet += "\t" + rNode.sNodeName + " [shape=hexagon];\n";
            for( const NodeNames& rIncoming : rNode.vIncomingFunctions )
                if( !vGraph[ rIncoming ].bHidden )
                    sRet += "\t" + vGraph[ rIncoming ].sNodeName + " -> " + rNode.sNodeName + ";\n";
        }
        return sRet + "}";
    }

    std::vector<std::pair<std::string, size_t>> getRuntimes( )
    {
        std::vector<std::pair<std::string, size_t>> vRet;
        vRet.reserve( NodeNames::SIZE );
        for( size_t uiNodeName = 0; uiNodeName < NodeNames::SIZE; uiNodeName++ )
            vRet.push_back( std::make_pair( vGraph[ uiNodeName ].sNodeName, vGraphData[ uiNodeName ].uiLastRuntime ) );

        return vRet;
    }
};

template <> pybind11::object PartialQuarry::getValue( std::vector<std::string> vKeys )
{
    return pyjson::from_json( getValue<json>( vKeys ) );
}

} // namespace cm
