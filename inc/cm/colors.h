#include "cm/computation.h"
#include <cmath>

#pragma once

namespace cm
{

void ContactMapping::setColors( )
{
    vColorPalette = colorPalette( this->xSession[ "settings" ][ "interface" ][ "color_palette" ].get<std::string>( ) );
}

void ContactMapping::setCombined( )
{
    vCombined.reserve( vvNormalized.size( ) );
    vCombined.clear( );

    for( auto& vArr : vvNormalized )
        vCombined.push_back( getMixedValue( vArr[ 0 ], vArr[ 1 ] ) );
}

double ContactMapping::logScale( double fX )
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

size_t ContactMapping::colorRange( double fX )
{
    double fMin = this->xSession[ "settings" ][ "normalization" ][ "color_range" ][ "val_min" ].get<double>( );
    double fMax = this->xSession[ "settings" ][ "normalization" ][ "color_range" ][ "val_max" ].get<double>( );
    if( fMin == fMax )
        return 0;
    fX = ( fX - fMin ) / ( fMax / fMin );
    return std::min( vColorPalette.size( ) - 1,
                     std::max( (size_t)0, (size_t)( (double)( vColorPalette.size( ) - 1 ) * fX ) ) );
}

void ContactMapping::setColored( )
{
    vColored.reserve( vCombined.size( ) );
    vColored.clear( );

    for( double fC : vCombined )
        vColored.push_back( vColorPalette[ colorRange( logScale( fC ) ) ] );
}


std::vector<std::string> ContactMapping::getColors()
{
    update(NodeNames::Colored);
    return vColored;
}

void ContactMapping::regColors( )
{
    registerNode( NodeNames::Colors,
                  ComputeNode{ .sNodeName = "color_palette",
                               .fFunc = &ContactMapping::setColors,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "settings", "interface", "color_palette" },
                                                     { "settings", "interface", "color_low" },
                                                     { "settings", "interface", "color_high" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::Combined,
                  ComputeNode{ .sNodeName = "combined",
                               .fFunc = &ContactMapping::setCombined,
                               .vIncomingFunctions = { NodeNames::Normalized, NodeNames::BetweenGroup },
                               .vIncomingSession = { },
                               .uiLastUpdated = uiCurrTime } );

    registerNode(
        NodeNames::Colored,
        ComputeNode{ .sNodeName = "colored",
                     .fFunc = &ContactMapping::setColored,
                     .vIncomingFunctions = { NodeNames::Colors, NodeNames::Combined, NodeNames::BetweenGroup },
                     .vIncomingSession = { { "settings", "normalization", "log_base", "val" },
                                           { "settings", "normalization", "color_range", "val_max" },
                                           { "settings", "normalization", "color_range", "val_min" } },
                     .uiLastUpdated = uiCurrTime } );
}

} // namespace cm