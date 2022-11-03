#include "cm/partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

bool PartialQuarry::setColors( )
{
    vColorPalette.clear( );
    vColorPalette = colorPalette( this->xSession[ "settings" ][ "interface" ][ "color_palette" ].get<std::string>( ),
                                  this->xSession[ "settings" ][ "interface" ][ "color_low" ].get<std::string>( ),
                                  this->xSession[ "settings" ][ "interface" ][ "color_high" ].get<std::string>( ) );
    if( iBetweenGroupSetting == 2 )
        sBackgroundColor = vColorPalette[ vColorPalette.size( ) / 2 ];
    else
        sBackgroundColor = vColorPalette.front( );
    END_RETURN;
}

bool PartialQuarry::setAnnotationColors( )
{
    vColorPaletteAnnotation.clear( );
    auto& rList = this->xSession[ "settings" ][ "interface" ][ "annotation_color_palette" ];
    vColorPaletteAnnotation.reserve( rList.size( ) );
    for( auto& rC : rList )
    {
        CANCEL_RETURN;
        vColorPaletteAnnotation.push_back( rC.get<std::string>( ) );
    }
    END_RETURN;
}

bool PartialQuarry::setCombined( )
{
    vCombined.clear( );
    vCombined.reserve( vvNormalized.size( ) );

    for( auto& vArr : vvNormalized )
    {
        CANCEL_RETURN;
        vCombined.push_back( getMixedValue( vArr[ 0 ], vArr[ 1 ] ) );
    }

    END_RETURN;
}

bool PartialQuarry::setScaled( )
{
    vScaled.clear( );
    vScaled.reserve( vDivided.size( ) );

    fMax = std::numeric_limits<double>::min( );
    fMin = std::numeric_limits<double>::max( );
    for( double fVal : vDivided )
    {
        CANCEL_RETURN;
        if( !std::isnan( fVal ) )
        {
            fMax = std::max( fMax, fVal );
            fMin = std::min( fMin, fVal );
        }
    }

    if( this->xSession[ "settings" ][ "normalization" ][ "scale" ].get<std::string>( ) == "minmax" )
    {
        for( double fVal : vDivided )
        {
            CANCEL_RETURN;
            if( std::isnan( fVal ) )
                vScaled.push_back( fVal );
            else
                vScaled.push_back( ( fVal + fMin ) / ( fMax - fMin ) );
        }
        fMin = 0;
        fMax = 1;
    }
    else if( this->xSession[ "settings" ][ "normalization" ][ "scale" ].get<std::string>( ) == "max" )
    {
        for( double fVal : vDivided )
        {
            CANCEL_RETURN;
            if( std::isnan( fVal ) )
                vScaled.push_back( fVal );
            else
                vScaled.push_back( fVal / fMax );
        }
        fMax = 1;
    }
    else if( this->xSession[ "settings" ][ "normalization" ][ "scale" ].get<std::string>( ) == "dont" )
    {
        for( double fVal : vDivided )
        {
            CANCEL_RETURN;
            vScaled.push_back( fVal );
        }
    }
    else
        throw std::logic_error( "invlaid scale value" );
    END_RETURN;
}

double PartialQuarry::logScale( double fX, double fBase )
{
    bool bLargerZero = fX > 0;
    if( iBetweenGroupSetting == 2 ) // between_groups == sub
        fX = std::abs( fX );

    double fRet;
    if( fBase == 0 )
        fRet = fX;
    else
    {
        double fA = std::pow( 2.0, fBase );
        fRet = std::log( fA * ( std::min( 1.0, fX ) * ( 1.0 - ( 1.0 / fA ) ) + ( 1.0 / fA ) ) ) / std::log( fA );
    }

    if( iBetweenGroupSetting == 2 ) // between_groups == sub
        fRet = fRet * ( bLargerZero ? 1 : -1 ) / 2.0 + 0.5;

    return fRet;
}

double PartialQuarry::colorRange( double fMin, double fMax, double fX )
{
    if( fMin == fMax )
        return 0;
    return ( fX - fMin ) / ( fMax - fMin );
}

size_t PartialQuarry::colorIndex( double fX )
{
    return std::min( vColorPalette.size( ) - 1,
                     (size_t)std::max( 0.0, ( (double)( vColorPalette.size( ) - 1 ) * fX ) ) );
}

bool PartialQuarry::setColored( )
{
    vColored.clear( );
    vColored.reserve( vScaled.size( ) );
    double fMin = this->xSession[ "settings" ][ "normalization" ][ "color_range" ][ "val_min" ].get<double>( );
    double fMax = this->xSession[ "settings" ][ "normalization" ][ "color_range" ][ "val_max" ].get<double>( );
    double fBase = this->xSession[ "settings" ][ "normalization" ][ "log_base" ][ "val" ].get<double>( );

    for( double fC : vScaled )
    {
        CANCEL_RETURN;
        if( std::isnan( fC ) )
            vColored.push_back( sBackgroundColor );
        else
            vColored.push_back( vColorPalette[ colorIndex( logScale( colorRange( fMin, fMax, fC ), fBase ) ) ] );
    }
    END_RETURN;
}

bool PartialQuarry::setPalette( )
{
    pybind11::gil_scoped_acquire acquire;
    pybind11::list vNew;
    double fMinR = this->xSession[ "settings" ][ "normalization" ][ "color_range" ][ "val_min" ].get<double>( );
    double fMaxR = this->xSession[ "settings" ][ "normalization" ][ "color_range" ][ "val_max" ].get<double>( );
    double fBase = this->xSession[ "settings" ][ "normalization" ][ "log_base" ][ "val" ].get<double>( );

    for( double fC = std::min( fMin, fMinR ); fC <= std::max( fMax, fMaxR );
         fC += ( std::max( fMax, fMaxR ) - std::min( fMin, fMinR ) ) / 255.0 )
    {
        CANCEL_RETURN;
        vNew.append( vColorPalette[ colorIndex( logScale( colorRange( fMinR, fMaxR, fC ), fBase ) ) ] );
    }

    vRenderedPalette = vNew;
    END_RETURN;
}

const pybind11::list PartialQuarry::getPalette( )
{
    update( NodeNames::Palette );
    return vRenderedPalette;
}

bool PartialQuarry::setHeatmapCDS( )
{
    using namespace pybind11::literals;
    pybind11::gil_scoped_acquire acquire;

    pybind11::list vScreenBottom;
    pybind11::list vScreenLeft;
    pybind11::list vScreenTop;
    pybind11::list vScreenRight;

    pybind11::list vColor;

    pybind11::list vChrX;
    pybind11::list vChrY;
    pybind11::list vIndexLeft;
    pybind11::list vIndexRight;
    pybind11::list vIndexBottom;
    pybind11::list vIndexTop;

    pybind11::list vChrXSym;
    pybind11::list vChrYSym;
    pybind11::list vIndexSymLeft;
    pybind11::list vIndexSymRight;
    pybind11::list vIndexSymBottom;
    pybind11::list vIndexSymTop;

    pybind11::list vScoreTotal;
    pybind11::list vScoreA;
    pybind11::list vScoreB;

    for( size_t uiI = 0; uiI < vColored.size( ); uiI++ )
    {
        CANCEL_RETURN;

        if( vColored[ uiI ] != sBackgroundColor )
        {
            vScreenBottom.append( vBinCoords[ uiI ][ 0 ].uiScreenY );
            vScreenLeft.append( vBinCoords[ uiI ][ 0 ].uiScreenX );
            vScreenTop.append( vBinCoords[ uiI ][ 0 ].uiScreenY + vBinCoords[ uiI ][ 0 ].uiH );
            vScreenRight.append( vBinCoords[ uiI ][ 0 ].uiScreenX + vBinCoords[ uiI ][ 0 ].uiW );

            vColor.append( vColored[ uiI ] );

            vChrX.append( vBinCoords[ uiI ][ 0 ].sChromosomeX );
            vChrY.append( vBinCoords[ uiI ][ 0 ].sChromosomeY );
            vIndexLeft.append( vBinCoords[ uiI ][ 0 ].uiIndexX * this->xSession[ "dividend" ].get<size_t>( ) );
            vIndexRight.append( ( vBinCoords[ uiI ][ 0 ].uiIndexX + vBinCoords[ uiI ][ 0 ].uiW ) *
                                this->xSession[ "dividend" ].get<size_t>( ) );
            vIndexBottom.append( vBinCoords[ uiI ][ 0 ].uiIndexY * this->xSession[ "dividend" ].get<size_t>( ) );
            vIndexTop.append( ( vBinCoords[ uiI ][ 0 ].uiIndexY + vBinCoords[ uiI ][ 0 ].uiH ) *
                              this->xSession[ "dividend" ].get<size_t>( ) );

            if( vBinCoords[ uiI ][ 1 ].sChromosomeX != "" )
            {
                vChrXSym.append( vBinCoords[ uiI ][ 1 ].sChromosomeX );
                vChrYSym.append( vBinCoords[ uiI ][ 1 ].sChromosomeY );
                vIndexSymLeft.append( vBinCoords[ uiI ][ 1 ].uiIndexX * this->xSession[ "dividend" ].get<size_t>( ) );
                vIndexSymRight.append( ( vBinCoords[ uiI ][ 1 ].uiIndexX + vBinCoords[ uiI ][ 1 ].uiW ) *
                                       this->xSession[ "dividend" ].get<size_t>( ) );
                vIndexSymBottom.append( vBinCoords[ uiI ][ 1 ].uiIndexY * this->xSession[ "dividend" ].get<size_t>( ) );
                vIndexSymTop.append( ( vBinCoords[ uiI ][ 1 ].uiIndexY + vBinCoords[ uiI ][ 1 ].uiH ) *
                                     this->xSession[ "dividend" ].get<size_t>( ) );
            }
            else
            {
                vChrXSym.append( nullptr );
                vChrYSym.append( nullptr );
                vIndexSymLeft.append( nullptr );
                vIndexSymRight.append( nullptr );
                vIndexSymBottom.append( nullptr );
                vIndexSymTop.append( nullptr );
            }

            vScoreTotal.append( vCombined[ uiI ] );
            vScoreA.append( vvNormalized[ uiI ][ 0 ] );
            vScoreB.append( vvNormalized[ uiI ][ 1 ] );
        }
    }


    xHeatmapCDS = pybind11::dict( "screen_bottom"_a = vScreenBottom,
                                  "screen_left"_a = vScreenLeft,
                                  "screen_top"_a = vScreenTop,
                                  "screen_right"_a = vScreenRight,

                                  "color"_a = vColor,

                                  "chr_x"_a = vChrX,
                                  "chr_y"_a = vChrY,
                                  "index_left"_a = vIndexLeft,
                                  "index_right"_a = vIndexRight,
                                  "index_bottom"_a = vIndexBottom,
                                  "index_top"_a = vIndexTop,

                                  "chr_x_symmetry"_a = vChrXSym,
                                  "chr_y_symmetry"_a = vChrYSym,
                                  "index_symmetry_left"_a = vIndexSymLeft,
                                  "index_symmetry_right"_a = vIndexSymRight,
                                  "index_symmetry_bottom"_a = vIndexSymBottom,
                                  "index_symmetry_top"_a = vIndexSymTop,

                                  "score_total"_a = vScoreTotal,
                                  "score_a"_a = vScoreA,
                                  "score_b"_a = vScoreB );
    END_RETURN;
}

const pybind11::dict PartialQuarry::getHeatmap( )
{
    update( NodeNames::HeatmapCDS );
    return xHeatmapCDS;
}

const std::string& PartialQuarry::getBackgroundColor( )
{
    update( NodeNames::Colors );
    return sBackgroundColor;
}

void PartialQuarry::regColors( )
{
    registerNode( NodeNames::Colors,
                  ComputeNode{ .sNodeName = "color_palette",
                               .fFunc = &PartialQuarry::setColors,
                               .vIncomingFunctions = { NodeNames::BetweenGroup },
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
                               .vIncomingFunctions = { NodeNames::Normalized },
                               .vIncomingSession = { },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::Scaled,
                  ComputeNode{ .sNodeName = "colored_bins",
                               .fFunc = &PartialQuarry::setScaled,
                               .vIncomingFunctions = { NodeNames::Divided },
                               .vIncomingSession = { { "settings", "normalization", "scale" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::Colored,
                  ComputeNode{ .sNodeName = "colored_bins",
                               .fFunc = &PartialQuarry::setColored,
                               .vIncomingFunctions = { NodeNames::Colors, NodeNames::Scaled },
                               .vIncomingSession = { { "settings", "normalization", "log_base", "val" },
                                                     { "settings", "normalization", "color_range", "val_max" },
                                                     { "settings", "normalization", "color_range", "val_min" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::HeatmapCDS,
                  ComputeNode{ .sNodeName = "heatmap_cds",
                               .fFunc = &PartialQuarry::setHeatmapCDS,
                               .vIncomingFunctions = { NodeNames::Colored },
                               .vIncomingSession = { },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::Palette,
                  ComputeNode{ .sNodeName = "color_palette",
                               .fFunc = &PartialQuarry::setPalette,
                               .vIncomingFunctions = { NodeNames::Scaled },
                               .vIncomingSession = { },
                               .uiLastUpdated = uiCurrTime } );
}

} // namespace cm