#include "cm/partial_quarry.h"
#include "sps/desc.h"

#pragma once

namespace cm
{

bool PartialQuarry::setActivateAnnotation( )
{
    for( size_t uiX : { 0, 1 } )
    {
        vActiveAnnotation[ uiX ].clear( );
        auto& rList = this->xSession[ "annotation" ][ uiX == 0 ? "visible_x" : "visible_y" ];
        vActiveAnnotation[ uiX ].reserve( rList.size( ) );
        for( auto& rRep : rList )
        {
            CANCEL_RETURN;
            vActiveAnnotation[ uiX ].push_back( rRep.get<std::string>( ) );
        }
    }
    END_RETURN;
}

bool PartialQuarry::setActivateAnnotationCDS( )
{
    pybind11::gil_scoped_acquire acquire;
    for( size_t uiX : { 0, 1 } )
    {
        pybind11::list xL;
        for( std::string sAnno : vActiveAnnotation[ uiX ] )
        {
            CANCEL_RETURN;
            xL.append( sAnno );
        }
        vActiveAnnotationCDS[ uiX ] = xL;
    }
    END_RETURN;
}

bool PartialQuarry::setAnnotationValues( )
{
    size_t uiMaxDetailedDisplay =
        this->xSession[ "settings" ][ "interface" ][ "max_detailed_anno_display" ].get<size_t>( );
    for( size_t uiX : { 0, 1 } )
    {
        size_t uiMinAnnoDist = 2;
        vAnnotationValues[ uiX ].clear( );
        vAnnotationValues[ uiX ].reserve( vActiveAnnotation[ uiX ].size( ) );
        vMaxRowsPerAnno[ uiX ].clear( );
        vMaxRowsPerAnno[ uiX ].reserve( vActiveAnnotation[ uiX ].size( ) );

        size_t uiTotalCount = 0;
        vMaxAnnoRows[ uiX ] = 1;
        for( std::string sCurrAnno : vActiveAnnotation[ uiX ] )
        {
            auto& rJson = this->xSession[ "annotation" ][ "by_name" ][ sCurrAnno ];
            for( AxisRegion& xRegion : vAxisRegions[ uiX ] )
            {
                CANCEL_RETURN;

                int64_t iDataSetId = rJson[ xRegion.sChromosome ].get<int64_t>( );

                uiTotalCount +=
                    xIndices.vAnno.count( iDataSetId, xRegion.uiIndexPos, xRegion.uiIndexPos + xRegion.uiIndexSize );
            }
        }
        for( std::string sCurrAnno : vActiveAnnotation[ uiX ] )
        {
            auto& rJson = this->xSession[ "annotation" ][ "by_name" ][ sCurrAnno ];
            vAnnotationValues[ uiX ].emplace_back( );
            if( uiTotalCount < uiMaxDetailedDisplay )
            {
                vAnnotationValues[ uiX ].back( ).second.reserve( uiTotalCount );
                std::vector<size_t> vEndPos;
                for( AxisRegion& xRegion : vAxisRegions[ uiX ] )
                {
                    int64_t iDataSetId = rJson[ xRegion.sChromosome ].get<int64_t>( );
                    for( auto& xAnno : xIndices.vAnno.query(
                             iDataSetId, xRegion.uiIndexPos, xRegion.uiIndexPos + xRegion.uiIndexSize ) )
                    {
                        size_t uiStartPos =
                            ( std::get<0>( xAnno ) > xRegion.uiIndexPos ? std::get<0>( xAnno ) - xRegion.uiIndexPos
                                                                        : 0 ) +
                            xRegion.uiScreenPos;
                        size_t uiEndPos = std::min( ( std::get<1>( xAnno ) > xRegion.uiIndexPos
                                                          ? std::get<1>( xAnno ) - xRegion.uiIndexPos
                                                          : 0 ) +
                                                        xRegion.uiScreenPos,
                                                    xRegion.uiScreenPos + xRegion.uiScreenSize );
                        if( uiEndPos > uiStartPos )
                        {
                            size_t uiRow = 0;
                            while( uiRow < vEndPos.size( ) && uiStartPos < vEndPos[ uiRow ] + uiMinAnnoDist )
                                ++uiRow;
                            if( uiRow == vEndPos.size( ) )
                                vEndPos.emplace_back( 0 );
                            vEndPos[ uiRow ] = uiEndPos;
                            vAnnotationValues[ uiX ].back( ).second.push_back(
                                Annotation{ .sInfo = std::get<2>( xAnno ),
                                            .bForw = std::get<3>( xAnno ),
                                            .uiScreenX = uiStartPos,
                                            .uiScreenY = uiEndPos,
                                            .uiIndexX = std::get<0>( xAnno ),
                                            .uiIndexY = std::get<1>( xAnno ),
                                            .uiRow = uiRow,
                                            .sChromosome = xRegion.sChromosome } );
                        }
                    }
                }
                vMaxRowsPerAnno[ uiX ].push_back( vEndPos.size( ) );
                vMaxAnnoRows[ uiX ] = std::max( vMaxAnnoRows[ uiX ], vEndPos.size( ) );
            }
            else
            {
                vAnnotationValues[ uiX ].back( ).first.reserve( vAxisCords[ uiX ].size( ) );
                for( AxisCoord& xCoords : vAxisCords[ uiX ] )
                {
                    int64_t iDataSetId = rJson[ xCoords.sChromosome ].get<int64_t>( );
                    vAnnotationValues[ uiX ].back( ).first.push_back( xIndices.vAnno.count(
                        iDataSetId, xCoords.uiIndexPos, xCoords.uiIndexPos + xCoords.uiIndexSize ) );
                }
            }
        }
    }
    END_RETURN;
}

/// http://www.martinbroadhurst.com/how-to-split-a-string-in-c.html
template <class Container_t> Container_t splitString( const std::string& sString, char sDelim = ' ' )
{
    std::stringstream xStream( sString );
    std::string sCurrToken;
    Container_t xCont;
    while( std::getline( xStream, sCurrToken, sDelim ) )
        xCont.push_back( sCurrToken );
    return xCont;
} // method

bool PartialQuarry::setAnnotationCDS( )
{
    using namespace pybind11::literals;
    pybind11::gil_scoped_acquire acquire;

    size_t uiDividend = this->xSession[ "dividend" ].get<size_t>( );

    double fMinAnnoDist = this->xSession[ "settings" ][ "interface" ][ "min_anno_dist" ].get<double>( );

    for( size_t uiX : { 0, 1 } )
    {
        pybind11::list vChr;
        pybind11::list vIndexStart;
        pybind11::list vIndexEnd;
        pybind11::list vAnnoName;
        pybind11::list vScreenStart;
        pybind11::list vScreenEnd;
        pybind11::list vColor;
        pybind11::list vNumAnno;
        pybind11::list vInfo;
        pybind11::list vSize;
        pybind11::list vStrand;
        pybind11::list vID;
        pybind11::list vDesc;

        double fAnnoHeight = ( 1.0 - fMinAnnoDist ) / vMaxAnnoRows[ uiX ];

        for( size_t uiN = 0; uiN < vActiveAnnotation[ uiX ].size( ); uiN++ )
        {
            std::string& rAnnoName = vActiveAnnotation[ uiX ][ uiN ];
            for( size_t uiA = 0; uiA < vAnnotationValues[ uiX ][ uiN ].first.size( ); uiA++ )
            {
                CANCEL_RETURN;

                size_t uiVal = vAnnotationValues[ uiX ][ uiN ].first[ uiA ];
                auto& rCoords = vAxisCords[ uiX ][ uiA ];

                if( uiVal > 0 )
                {
                    vChr.append( rCoords.sChromosome );
                    vIndexStart.append( readableBp( rCoords.uiIndexPos * uiDividend ) );
                    vIndexEnd.append( readableBp( ( rCoords.uiIndexPos + rCoords.uiIndexSize ) * uiDividend ) );
                    vAnnoName.append( py::make_tuple( rAnnoName, 0 /*anno overlap*/ ) );
                    vScreenStart.append( rCoords.uiScreenPos );
                    vScreenEnd.append( rCoords.uiScreenPos + rCoords.uiScreenSize );
                    vColor.append( vColorPaletteAnnotation[ uiN % vColorPaletteAnnotation.size( ) ] );
                    vNumAnno.append( uiVal );
                    vInfo.append( "n/a" );
                    vSize.append( fAnnoHeight * ( 1.0 - fMinAnnoDist ) );
                    vStrand.append( "n/a" );
                    vID.append( "n/a" );
                    vDesc.append( "n/a" );
                }
            }
            for( Annotation& rA : vAnnotationValues[ uiX ][ uiN ].second )
            {
                CANCEL_RETURN;

                vChr.append( rA.sChromosome );
                vIndexStart.append( readableBp( rA.uiIndexX * uiDividend ) );
                vIndexEnd.append( readableBp( rA.uiIndexY * uiDividend ) );
                vAnnoName.append( py::make_tuple(
                    rAnnoName,
                    ( rA.uiRow % 2 == 0 ? 1.0 : -1.0 ) * fAnnoHeight * ( ( rA.uiRow + 1 ) / 2 ) +
                        ( vMaxRowsPerAnno[ uiX ][ uiN ] % 2 == 0 ? 0.5 * fAnnoHeight : 0 ) ) ); // uiRow - 1
                vScreenStart.append( rA.uiScreenX );
                vScreenEnd.append( rA.uiScreenY );
                vColor.append( rA.bForw ? vColorPaletteAnnotation[ uiN % vColorPaletteAnnotation.size( ) ]
                                        : vColorPaletteAnnotationDark[ uiN % vColorPaletteAnnotationDark.size( ) ] );
                vNumAnno.append( 1 );
                vSize.append( fAnnoHeight * ( 1.0 - fMinAnnoDist ) );
                vStrand.append( rA.bForw ? "+" : "-" );
                std::vector<std::string> vSplit = splitString<std::vector<std::string>>( rA.sInfo, '\n' );
                std::string sId = "n/a";
                std::string sDesc = "n/a";
                std::string sInfo = "";
                const std::string csID = "ID=";
                const std::string csDesc = "description=";
                for( std::string sX : vSplit )
                {
                    if( sX.substr( 0, csID.size( ) ) == csID )
                        sId = sX.substr( csID.size( ) );
                    else if( sX.substr( 0, csDesc.size( ) ) == csDesc )
                        sDesc = sX.substr( csDesc.size( ) );
                    else
                        sInfo += sX + "\n";
                }
                vInfo.append( sInfo.substr( 0, sInfo.size( ) - 1 ) );
                vID.append( sId );
                vDesc.append( sDesc );
            }
        }

        vAnnotationCDS[ uiX ] = pybind11::dict( "chr"_a = vChr,
                                                "index_start"_a = vIndexStart,
                                                "index_end"_a = vIndexEnd,
                                                "anno_name"_a = vAnnoName,
                                                "screen_start"_a = vScreenStart,
                                                "screen_end"_a = vScreenEnd,
                                                "color"_a = vColor,
                                                "num_anno"_a = vNumAnno,
                                                "info"_a = vInfo,
                                                "size"_a = vSize,
                                                "strand"_a = vStrand,
                                                "id"_a = vID,
                                                "desc"_a = vDesc );
    }
    END_RETURN;
}

const pybind11::dict PartialQuarry::getAnnotation( bool bXAxis )
{
    update( NodeNames::AnnotationCDS );
    return vAnnotationCDS[ bXAxis ? 0 : 1 ];
}

const pybind11::list PartialQuarry::getDisplayedAnnos( bool bXAxis )
{
    update( NodeNames::ActivateAnnotationCDS );
    return vActiveAnnotationCDS[ bXAxis ? 0 : 1 ];
}

void PartialQuarry::regAnnotation( )
{
    registerNode( NodeNames::ActivateAnnotation,
                  ComputeNode{ .sNodeName = "active_annotation",
                               .fFunc = &PartialQuarry::setActivateAnnotation,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "annotation", "visible_x" }, { "annotation", "visible_y" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::AnnotationValues,
                  ComputeNode{ .sNodeName = "annotation_values",
                               .fFunc = &PartialQuarry::setAnnotationValues,
                               .vIncomingFunctions = { NodeNames::ActivateAnnotation, NodeNames::AxisCoords },
                               .vIncomingSession = { { "annotation", "by_name" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::AnnotationCDS,
                  ComputeNode{ .sNodeName = "annotation_cds",
                               .fFunc = &PartialQuarry::setAnnotationCDS,
                               .vIncomingFunctions = { NodeNames::AnnotationValues, NodeNames::AnnotationColors },
                               .vIncomingSession = { },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::ActivateAnnotationCDS,
                  ComputeNode{ .sNodeName = "active_annotation_cds",
                               .fFunc = &PartialQuarry::setActivateAnnotationCDS,
                               .vIncomingFunctions = { NodeNames::ActivateAnnotation },
                               .vIncomingSession = { },
                               .uiLastUpdated = uiCurrTime } );
}

} // namespace cm
