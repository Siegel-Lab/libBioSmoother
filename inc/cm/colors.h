#include "cm/partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

void PartialQuarry::setColors( )
{
    vColorPalette = colorPalette( this->xSession[ "settings" ][ "interface" ][ "color_palette" ].get<std::string>( ),
                                  this->xSession[ "settings" ][ "interface" ][ "color_low" ].get<std::string>( ),
                                  this->xSession[ "settings" ][ "interface" ][ "color_high" ].get<std::string>( ) );
}

void PartialQuarry::setAnnotationColors( )
{
    vColorPaletteAnnotation.clear( );
    auto& rList = this->xSession[ "settings" ][ "interface" ][ "annotation_color_palette" ];
    vColorPaletteAnnotation.reserve( rList.size( ) );
    for( auto& rC : rList )
        vColorPaletteAnnotation.push_back( rC.get<std::string>( ) );
}

void PartialQuarry::setCombined( )
{
    vCombined.clear( );
    vCombined.reserve( vvNormalized.size( ) );

    for( auto& vArr : vvNormalized )
        vCombined.push_back( getMixedValue( vArr[ 0 ], vArr[ 1 ] ) );
}

double PartialQuarry::logScale( double fX )
{
    bool bLargerZero = fX > 0;
    if( iBetweenGroupSetting == 2 ) // between_groups == sub
        fX = std::abs( fX );

    double fBase = this->xSession[ "settings" ][ "normalization" ][ "log_base" ][ "val" ].get<double>( );
    if( fBase == 0 )
        return fX;
    double fA = std::pow( 2.0, fBase );
    double fRet = std::log( fA * ( std::min( 1.0, fX ) * ( 1.0 - ( 1.0 / fA ) ) + ( 1.0 / fA ) ) ) / std::log( fA );

    if( iBetweenGroupSetting == 2 ) // between_groups == sub
        fRet = fRet * ( bLargerZero ? 1 : -1 ) / 2.0 + 0.5;

    return fRet;
}

size_t PartialQuarry::colorRange( double fX )
{
    double fMin = this->xSession[ "settings" ][ "normalization" ][ "color_range" ][ "val_min" ].get<double>( );
    double fMax = this->xSession[ "settings" ][ "normalization" ][ "color_range" ][ "val_max" ].get<double>( );
    if( fMin == fMax )
        return 0;
    fX = ( fX - fMin ) / ( fMax / fMin );
    return std::min( vColorPalette.size( ) - 1,
                     std::max( (size_t)0, (size_t)( (double)( vColorPalette.size( ) - 1 ) * fX ) ) );
}

void PartialQuarry::setColored( )
{
    vColored.clear( );
    vColored.reserve( vCombined.size( ) );

    for( double fC : vCombined )
        vColored.push_back( vColorPalette[ colorRange( logScale( fC ) ) ] );
}


const std::vector<std::string>& PartialQuarry::getColors( )
{
    update( NodeNames::Colored );
    return vColored;
}

void PartialQuarry::regColors( )
{
    registerNode( NodeNames::Colors,
                  ComputeNode{ .sNodeName = "color_palette",
                               .fFunc = &PartialQuarry::setColors,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "settings", "interface", "color_palette" },
                                                     { "settings", "interface", "color_low" },
                                                     { "settings", "interface", "color_high" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::AnnotationColors,
                  ComputeNode{ .sNodeName = "annotation_color_palette",
                               .fFunc = &PartialQuarry::setAnnotationColors,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "settings", "interface", "annotation_color_palette" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::Combined,
                  ComputeNode{ .sNodeName = "combined_bins",
                               .fFunc = &PartialQuarry::setCombined,
                               .vIncomingFunctions = { NodeNames::Normalized, NodeNames::BetweenGroup },
                               .vIncomingSession = { },
                               .uiLastUpdated = uiCurrTime } );

    registerNode(
        NodeNames::Colored,
        ComputeNode{ .sNodeName = "colored_bins",
                     .fFunc = &PartialQuarry::setColored,
                     .vIncomingFunctions = { NodeNames::Colors, NodeNames::Combined, NodeNames::BetweenGroup },
                     .vIncomingSession = { { "settings", "normalization", "log_base", "val" },
                                           { "settings", "normalization", "color_range", "val_max" },
                                           { "settings", "normalization", "color_range", "val_min" } },
                     .uiLastUpdated = uiCurrTime } );
}

} // namespace cm