#include "cm/sps_interface.h"
#include <nlohmann/json.hpp>
#include <pybind11_json/pybind11_json.hpp>

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


class PartialQuarry
{
  public:
    size_t uiVerbosity = 3;

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
        SIZE
    };
    struct ComputeNode
    {
        std::string sNodeName;
        void ( PartialQuarry::*fFunc )( void );
        std::vector<NodeNames> vIncomingFunctions;
        std::vector<std::vector<std::string>> vIncomingSession;
        size_t uiLastUpdated;
    };

    size_t uiCurrTime;
    std::vector<ComputeNode> vGraph;

    json xSession;
    std::map<std::vector<std::string>, size_t> xSessionTime;


    nlohmann::json::json_pointer toPointer( std::vector<std::string>& vKeys )
    {
        std::string sPtr = "";
        for( std::string sKey : vKeys )
            sPtr += "/" + sKey;
        return nlohmann::json::json_pointer( sPtr );
    }

    bool updateSettings( const json& xNewSettings, std::vector<std::string> vPrefix = { } )
    {
        bool bUpdated = false;
        for( auto& xEle : xNewSettings.items( ) )
        {
            if( xEle.key( ) != "previous" && xEle.key( ) != "next" )
            {
                std::vector<std::string> vPrefix2( vPrefix );
                vPrefix2.push_back( xEle.key( ) );

                if( xNewSettings[ toPointer( vPrefix2 ) ].is_object( ) )
                    bUpdated = updateSettings( xNewSettings, vPrefix2 ) || bUpdated;
                else if( this->xSession[ toPointer( vPrefix2 ) ] != xNewSettings[ toPointer( vPrefix2 ) ] )
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

    size_t update( NodeNames xNodeName )
    {
        ComputeNode& xNode = vGraph[ xNodeName ];
        if( xNode.uiLastUpdated < uiCurrTime )
        {
            std::vector<std::string> vOutOfDateReason{ };
            auto p1 = std::chrono::high_resolution_clock::now( );
            for( NodeNames xPred : xNode.vIncomingFunctions )
                if( update( xPred ) > xNode.uiLastUpdated )
                    vOutOfDateReason.push_back( "change in previous node " + vGraph[ xPred ].sNodeName );
            auto p2 = std::chrono::high_resolution_clock::now( );
            auto pms_int = std::chrono::duration_cast<std::chrono::milliseconds>( p2 - p1 );

            for( std::vector<std::string> vSetting : xNode.vIncomingSession )
                if( xSessionTime.count( vSetting ) == 0 || xSessionTime[ vSetting ] > xNode.uiLastUpdated )
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

                ( this->*xNode.fFunc )( );

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


  public:
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
        auto xRet = this->xSession[ toPointer( vKeys ) ].get<T>( );

        return xRet;
    }

    template <typename T> void setValue( std::vector<std::string> vKeys, T xVal )
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
        xSessionTime[ vCurrKey ] = uiCurrTime;
        for( std::string sKey : vKeys )
        {
            vCurrKey.push_back( sKey );
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
    std::vector<double> vScaled;
    std::vector<std::string> vColored;

    std::string sBackgroundColor;

    std::array<std::vector<std::vector<size_t>>, 2> vvCoverageValues;
    std::array<std::array<std::vector<size_t>, 2>, 2> vInGroupCoverage;
    // outer array: column / row
    std::array<std::vector<double>, 2> vvFlatCoverageValues;

    std::vector<std::string> vColorPalette;
    std::vector<std::string> vColorPaletteAnnotation;
    std::array<std::vector<std::string>, 2> vActiveAnnotation;
    std::array<std::vector<std::vector<size_t>>, 2> vAnnotationValues;
    std::array<pybind11::dict, 2> vAnnotationCDS;
    std::array<pybind11::list, 2> vActiveAnnotationCDS;
    pybind11::dict xHeatmapCDS;
    std::array<pybind11::dict, 2> xTicksCDS;
    std::array<pybind11::list, 2> vTickLists;
    std::array<size_t, 2> vCanvasSize;


    // bin_size.h
    size_t nextEvenNumber( double );


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
    void setTicks( );

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
    size_t symmetry( size_t, size_t );

    // replicates.h
    void setBinValues( );

    // replicates.h
    size_t getFlatValue( std::vector<size_t> );
    // replicates.h
    double getMixedValue( double, double );

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
    void normalizeSize( size_t );
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
    void setAnnotationColors( );
    // colors.h
    void setCombined( );
    // colors.h
    void setScaled( );
    // colors.h
    void setColored( );
    // colors.h
    double logScale( double );
    // colors.h
    size_t colorRange( double );
    // colors.h
    void setHeatmapCDS( );

    // colors.h
    void regColors( );


    // annotation.h
    void setActivateAnnotation( );
    // annotation.h
    void setAnnotationValues( );
    // annotation.h
    void setAnnotationCDS( );
    // annotation.h
    void setActivateAnnotationCDS( );
    // annotation.h
    void regAnnotation( );

  protected:
    virtual std::vector<std::array<double, 2>> normalizeBinominalTestTrampoline( std::vector<std::array<size_t, 2>>&,
                                                                                 size_t, size_t, double )
    {
        throw std::logic_error( "Function not implemented" );
    }

    virtual std::vector<std::string> colorPalette( std::string, std::string, std::string )
    {
        throw std::logic_error( "Function not implemented" );
    }

  public:
    PartialQuarry( ) : uiCurrTime( 0 ), vGraph( NodeNames::SIZE ), xSession( ), xIndices( )
    {
        regBinSize( );
        regCoords( );
        regReplicates( );
        regCoverage( );
        regNormalization( );
        regColors( );
        regAnnotation( );
    }

    PartialQuarry( std::string sPrefix )
        : uiCurrTime( 0 ),
          vGraph( NodeNames::SIZE ),
          xSession( json::parse( std::ifstream( sPrefix + "/default_session.json" ) ) ),
          xIndices( sPrefix )
    {
        regBinSize( );
        regCoords( );
        regReplicates( );
        regCoverage( );
        regNormalization( );
        regColors( );
        regAnnotation( );
    }

    virtual ~PartialQuarry( )
    {}

    // colors.h
    const std::vector<std::string>& getColors( );

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
    const std::array<size_t, 2> getCanvasSize();

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
};

template <> pybind11::object PartialQuarry::getValue( std::vector<std::string> vKeys )
{
    return pyjson::from_json( getValue<json>( vKeys ) );
}

} // namespace cm
