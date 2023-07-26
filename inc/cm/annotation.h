#include "cm/partial_quarry.h"
#include "sps/desc.h"
#include <algorithm>

#pragma once

namespace cm
{


std::vector<std::string> shortAnnoNames( const std::vector<std::string>& rAnnos, size_t uiMaxChars )
{
    uiMaxChars = std::max( uiMaxChars, size_t( 7 ) );
    std::vector<std::string> vShortAnnos;
    std::map<std::string, size_t> xAnnoCount;

    vShortAnnos.reserve( rAnnos.size( ) );
    for( auto& sAnno : rAnnos )
    {
        std::string sCurrName;
        if( sAnno.size( ) <= uiMaxChars )
            sCurrName = sAnno;
        else
            sCurrName = sAnno.substr( 0, uiMaxChars - 3 ) + "...";

        if( xAnnoCount.count( sCurrName ) == 0 )
            xAnnoCount[ sCurrName ] = 0;

        if( xAnnoCount[ sCurrName ] == 0 )
            vShortAnnos.push_back( sCurrName );
        else
        {
            if( sAnno.size( ) + 4 <= uiMaxChars )
                vShortAnnos.push_back( sAnno + " (" + std::to_string( xAnnoCount[ sCurrName ] ) + ")" );
            else
                vShortAnnos.push_back( sAnno.substr( 0, uiMaxChars - 7 ) + "... (" +
                                       std::to_string( xAnnoCount[ sCurrName ] ) + ")" );
        }
        xAnnoCount[ sCurrName ] += 1;
    }

    return vShortAnnos;
}

bool PartialQuarry::setActivateAnnotation( )
{
    const size_t uiMaxChar = getValue<size_t>( { "settings", "interface", "axis_label_max_char", "val" } );
    for( size_t uiX : { 0, 1 } )
    {
        vActiveAnnotation[ uiX ].clear( );
        auto rList = getValue<json>( { "annotation", uiX == 0 ? "visible_x" : "visible_y" } );
        vActiveAnnotation[ uiX ].reserve( rList.size( ) );
        for( auto& rRep : rList )
        {
            CANCEL_RETURN;
            std::string sAnno = rRep.get<std::string>( );
            vActiveAnnotation[ uiX ].push_back( sAnno );
        }
        vShortAnnotation[ uiX ] = shortAnnoNames( vActiveAnnotation[ uiX ], uiMaxChar );
    }
    END_RETURN;
}

bool PartialQuarry::setActivateAnnotationCDS( )
{
    pybind11::gil_scoped_acquire acquire;
    for( size_t uiX : { 0, 1 } )
    {
        pybind11::list xL;
        for( std::string sAnno : vShortAnnotation[ uiX ] )
        {
            CANCEL_RETURN;
            xL.append( sAnno );
        }
        if( getValue<bool>( { "settings", "interface", "v4c", uiX == 0 ? "do_row" : "do_col" } ) )
            xL.append( "V4C" );
        vActiveAnnotationCDS[ uiX ] = xL;
    }
    END_RETURN;
}


size_t multiple_bins( std::string sVal )
{
    if( sVal == "separate" ) // make several bins for the annotation
        return 0;
    if( sVal == "stretch" ) // stretch one bin over entire annotation
        return 1;
    if( sVal == "squeeze" ) // squeeze annotation in a bin of size 1
        return 2;
    throw std::logic_error( "unknown anno_in_multiple_bins value" );
}

bool PartialQuarry::setAnnotationValues( )
{
    size_t uiMaxDetailedDisplay = getValue<size_t>( { "settings", "interface", "max_detailed_anno_display" } );
    const bool bSqueeze = getValue<std::string>( { "settings", "filters", "anno_in_multiple_bins" } ) == "squeeze";
    const uint32_t uiDividend = getValue<uint32_t>( { "dividend" } );
    for( size_t uiX : { 0, 1 } )
    {
        const bool bIsGenomeCoords =
            !getValue<bool>( { "settings", "filters", ( uiX == 0 ) ? "anno_coords_col" : "anno_coords_row" } );
        size_t uiMinAnnoDist = 2;
        vAnnotationValues[ uiX ].clear( );
        vAnnotationValues[ uiX ].reserve( vActiveAnnotation[ uiX ].size( ) );
        vMaxRowsPerAnno[ uiX ].clear( );
        vMaxRowsPerAnno[ uiX ].reserve( vActiveAnnotation[ uiX ].size( ) );

        size_t uiTotalCount = 0;
        vMaxAnnoRows[ uiX ] = 1;
        for( std::string sCurrAnno : vActiveAnnotation[ uiX ] )
        {
            auto uiFistAnnoIdx = getValue<size_t>( { "annotation", "by_name", sCurrAnno } );
            for( AxisRegion& xRegion : vAxisRegions[ uiX ] )
            {
                CANCEL_RETURN;

                size_t iDataSetId = uiFistAnnoIdx + vActiveChromosomes[ uiX ][ xRegion.uiChromosome ].uiActualContigId;

                uiTotalCount += pIndices->vAnno.count( iDataSetId, xRegion.uiIndexPos,
                                                       ( xRegion.uiIndexPos + xRegion.uiIndexSize ), false, false );
            }
        }
        for( std::string sCurrAnno : vActiveAnnotation[ uiX ] )
        {
            auto uiFistAnnoIdx = getValue<size_t>( { "annotation", "by_name", sCurrAnno } );
            vAnnotationValues[ uiX ].emplace_back( );
            if( uiTotalCount < uiMaxDetailedDisplay )
            {
                vAnnotationValues[ uiX ].back( ).second.reserve( uiTotalCount );
                std::vector<size_t> vEndPos;
                for( AxisRegion& xRegion : vAxisRegions[ uiX ] )
                {
                    size_t iDataSetId =
                        uiFistAnnoIdx + vActiveChromosomes[ uiX ][ xRegion.uiChromosome ].uiActualContigId;
                    for( auto& xAnno : pIndices->vAnno.query( iDataSetId, xRegion.uiIndexPos,
                                                              xRegion.uiIndexPos + xRegion.uiIndexSize, false, false ) )
                    {
                        double uiStartPos;
                        double uiEndPos;
                        if( bSqueeze && !bIsGenomeCoords )
                        {
                            uiStartPos = xRegion.uiScreenPos;
                            uiEndPos = ( xRegion.uiScreenPos + xRegion.uiScreenSize );
                        }
                        else
                        {
                            double fAnnoF = (double)std::get<0>( xAnno );
                            double fAnnoT = (double)std::get<1>( xAnno );
                            if( !bSqueeze || bIsGenomeCoords )
                            {
                                fAnnoF /= (double)uiDividend;
                                fAnnoT /= (double)uiDividend;
                            }
                            double fRegF = xRegion.uiIndexPos;
                            double fRegF2 = xRegion.uiScreenPos;
                            double fRegT = xRegion.uiScreenPos + xRegion.uiScreenSize;
                            uiStartPos = ( fAnnoF > fRegF ? fAnnoF - fRegF : 0.0 ) + fRegF2;
                            uiEndPos = std::min( ( fAnnoT > fRegF ? fAnnoT - fRegF : 0.0 ) + fRegF2, fRegT );
                        }

                        if( uiEndPos >= uiStartPos )
                        {
                            size_t uiRow = 0;
                            while( uiRow < vEndPos.size( ) && uiStartPos < vEndPos[ uiRow ] + uiMinAnnoDist )
                                ++uiRow;
                            if( uiRow == vEndPos.size( ) )
                                vEndPos.emplace_back( 0 );
                            vEndPos[ uiRow ] = uiEndPos;
                            vAnnotationValues[ uiX ].back( ).second.push_back(
                                Annotation{ /*.sInfo =*/std::get<2>( xAnno ),
                                            /*.bForw =*/std::get<3>( xAnno ),
                                            /*.fScreenX =*/uiStartPos,
                                            /*.fScreenY =*/uiEndPos,
                                            /*.uiIndexX =*/std::get<0>( xAnno ),
                                            /*.uiIndexY =*/std::get<1>( xAnno ),
                                            /*.uiRow =*/uiRow,
                                            /*.uiChromosome =*/xRegion.uiChromosome } );
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
                    const size_t uiChrIdx = vActiveChromosomes[ uiX ][ xCoords.uiChromosome ].uiActualContigId;
                    int64_t iDataSetId = uiFistAnnoIdx + uiChrIdx;
                    vAnnotationValues[ uiX ].back( ).first.push_back( pIndices->vAnno.count(
                        iDataSetId, xCoords.uiIndexPos, ( xCoords.uiIndexPos + xCoords.uiIndexSize ), false, false ) );
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

    size_t uiDividend = getValue<size_t>( { "dividend" } );

    double fMinAnnoDist = getValue<double>( { "settings", "interface", "min_anno_dist" } );

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
            std::string rAnnoName = vShortAnnotation[ uiX ][ uiN ];
            for( size_t uiA = 0; uiA < vAnnotationValues[ uiX ][ uiN ].first.size( ); uiA++ )
            {
                CANCEL_RETURN;

                size_t uiVal = vAnnotationValues[ uiX ][ uiN ].first[ uiA ];
                auto& rCoords = vAxisCords[ uiX ][ uiA ];

                if( uiVal > 0 )
                {
                    std::string sChromName = vActiveChromosomes[ uiX ][ rCoords.uiChromosome ].sName;
                    vChr.append( sChromName );
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

                std::string sChromName = vActiveChromosomes[ uiX ][ rA.uiChromosome ].sName;
                vChr.append( sChromName );
                vIndexStart.append( readableBp( rA.uiIndexX ) );
                vIndexEnd.append( readableBp( rA.uiIndexY ) );
                vAnnoName.append( py::make_tuple(
                    rAnnoName,
                    ( rA.uiRow % 2 == 0 ? 1.0 : -1.0 ) * fAnnoHeight * ( ( rA.uiRow + 1 ) / 2 ) +
                        ( vMaxRowsPerAnno[ uiX ][ uiN ] % 2 == 0 ? 0.5 * fAnnoHeight : 0 ) ) ); // uiRow - 1
                vScreenStart.append( rA.fScreenX );
                vScreenEnd.append( rA.fScreenY );
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
        for( const AxisCoord& rCoords : vV4cCoords[ uiX ] )
        {
            std::string sChromName = vActiveChromosomes[ uiX ][ rCoords.uiChromosome ].sName;
            vChr.append( sChromName );
            vIndexStart.append( readableBp( rCoords.uiIndexPos * uiDividend ) );
            vIndexEnd.append( readableBp( ( rCoords.uiIndexPos + rCoords.uiIndexSize ) * uiDividend ) );
            vAnnoName.append( py::make_tuple( "V4C", 0 /*anno overlap*/ ) );
            vScreenStart.append( rCoords.uiScreenPos );
            vScreenEnd.append( rCoords.uiScreenPos + rCoords.uiScreenSize );
            vColor.append(
                vColorPaletteAnnotation[ vActiveAnnotation[ uiX ].size( ) % vColorPaletteAnnotation.size( ) ] );
            vNumAnno.append( "n/a" );
            vInfo.append( "n/a" );
            vSize.append( fAnnoHeight * ( 1.0 - fMinAnnoDist ) );
            vStrand.append( "n/a" );
            vID.append( "n/a" );
            vDesc.append( "The region that ist selected for Virtual4C" );
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

const pybind11::dict PartialQuarry::getAnnotation( bool bXAxis,
                                                   const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::AnnotationCDS, fPyPrint );
    return vAnnotationCDS[ bXAxis ? 0 : 1 ];
}

const pybind11::list PartialQuarry::getDisplayedAnnos( bool bXAxis,
                                                       const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::ActivateAnnotationCDS, fPyPrint );
    return vActiveAnnotationCDS[ bXAxis ? 0 : 1 ];
}

void PartialQuarry::regAnnotation( )
{
    registerNode( NodeNames::ActivateAnnotation,
                  ComputeNode{ /*.sNodeName =*/"active_annotation",
                               /*.fFunc =*/&PartialQuarry::setActivateAnnotation,
                               /*.vIncomingFunctions =*/{ },
                               /*.vIncomingSession =*/
                               { { "annotation", "visible_x" },
                                 { "annotation", "visible_y" },
                                 { "settings", "interface", "axis_label_max_char", "val" } },
                               /*.vSessionsIncomingInPrevious =*/{ },
                               /*bHidden =*/false } );

    registerNode( NodeNames::AnnotationValues,
                  ComputeNode{ /*.sNodeName =*/"annotation_values",
                               /*.fFunc =*/&PartialQuarry::setAnnotationValues,
                               /*.vIncomingFunctions =*/{ NodeNames::ActivateAnnotation, NodeNames::AxisCoords },
                               /*.vIncomingSession =*/{ { "settings", "interface", "max_detailed_anno_display" } },
                               /*.vSessionsIncomingInPrevious =*/
                               { { "settings", "filters", "anno_in_multiple_bins" },
                                 { "dividend" },
                                 { "annotation", "by_name" },
                                 { "settings", "filters", "anno_coords_row" },
                                 { "settings", "filters", "anno_coords_col" } },
                               /*bHidden =*/false } );

    registerNode(
        NodeNames::AnnotationCDS,
        ComputeNode{
            /*.sNodeName =*/"annotation_cds",
            /*.fFunc =*/&PartialQuarry::setAnnotationCDS,
            /*.vIncomingFunctions =*/{ NodeNames::AnnotationValues, NodeNames::AnnotationColors, NodeNames::V4cCoords },
            /*.vIncomingSession =*/
            { { "settings", "interface", "min_anno_dist" } },
            /*.vSessionsIncomingInPrevious =*/{ { "dividend" } },
            /*bHidden =*/false } );

    registerNode( NodeNames::ActivateAnnotationCDS,
                  ComputeNode{ /*.sNodeName =*/"active_annotation_cds",
                               /*.fFunc =*/&PartialQuarry::setActivateAnnotationCDS,
                               /*.vIncomingFunctions =*/{ NodeNames::ActivateAnnotation },
                               /*.vIncomingSession =*/
                               {
                                   { "settings", "interface", "v4c", "do_col" },
                                   { "settings", "interface", "v4c", "do_row" },
                               },
                               /*.vSessionsIncomingInPrevious =*/{ },
                               /*bHidden =*/false } );
}

} // namespace cm
