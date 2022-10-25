#include "cm/partial_quarry.h"

#pragma once

namespace cm
{

void PartialQuarry::setActivateAnnotation( )
{
    for( size_t uiX : { 0, 1 } )
    {
        vActiveAnnotation[ uiX ].clear( );
        auto& rList = this->xSession[ "annotation" ][ uiX == 0 ? "visible_x" : "visible_y" ];
        vActiveAnnotation[ uiX ].resize( rList.size( ) );
        for( auto& rRep : rList )
            vActiveAnnotation[ uiX ].push_back( rRep.get<std::string>( ) );
    }
}

void PartialQuarry::setAnnotationValues( )
{
    for( size_t uiX : { 0, 1 } )
    {
        vAnnotationValues[ uiX ].clear( );
        vAnnotationValues[ uiX ].resize( vActiveAnnotation[ uiX ].size( ) );

        for( std::string sCurrAnno : vActiveAnnotation[ uiX ] )
        {
            vAnnotationValues[ uiX ].emplace_back( );
            vAnnotationValues[ uiX ].back( ).reserve( vAxisCords[ uiX ].size( ) );
            auto& rJson = this->xSession[ "annotation" ][ "by_name" ][ sCurrAnno ];
            for( AxisCoord& xCoords : vAxisCords[ uiX ] )
            {
                if( xCoords.sChromosome != "" )
                {
                    int64_t iDataSetId = rJson[ xCoords.sChromosome ].get<int64_t>( );
                    vAnnotationValues[ uiX ].back( ).push_back(
                        xIndices.getIndex<1, 1>( )->count( iDataSetId,
                                                           { xCoords.uiIndexPos },
                                                           { xCoords.uiIndexPos + xCoords.uiSize },
                                                           sps::IntersectionType::overlaps,
                                                           0 ) );
                }
                else
                    vAnnotationValues[ uiX ].back( ).push_back( false );
            }
        }
    }
}

void PartialQuarry::setAnnotationCDS( )
{
    using namespace pybind11::literals;
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

        for( size_t uiN = 0; uiN < vActiveAnnotation[ uiX ].size( ); uiN++ )
        {
            std::string& rAnnoName = vActiveAnnotation[ uiX ][ uiN ];
            for( size_t uiA = 0; uiA < vAnnotationValues[ uiX ][ uiN ].size( ); uiA++ )
            {
                size_t uiVal = vAnnotationValues[ uiX ][ uiN ][ uiA ];
                auto& rCoords = vAxisCords[ uiX ][ uiA ];

                if( uiVal > 0 )
                {
                    vChr.append( rCoords.sChromosome );
                    vIndexStart.append( rCoords.uiIndexPos );
                    vIndexEnd.append( rCoords.uiIndexPos + rCoords.uiSize );
                    vAnnoName.append( rAnnoName );
                    vScreenStart.append( rCoords.uiScreenPos );
                    vScreenEnd.append( rCoords.uiScreenPos + rCoords.uiSize );
                    vColor.append( vColorPaletteAnnotation[ uiN % vColorPaletteAnnotation.size( ) ] );
                    vNumAnno.append( uiVal );
                    vInfo.append( "" );
                }
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
                                                "info"_a = vInfo );
    }
}

const pybind11::dict& PartialQuarry::getAnnotation( bool bXAxis )
{
    update( NodeNames::AnnotationCDS );
    return vAnnotationCDS[ bXAxis ? 0 : 1 ];
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
}

} // namespace cm
