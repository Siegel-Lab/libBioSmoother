#include "cm/partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

bool PartialQuarry::setColors( )
{
    vColorPalette.clear( );
    vColorPalette = colorPalette( getValue<std::string>( { "settings", "interface", "color_palette" } ),
                                  getValue<std::string>( { "settings", "interface", "color_low" } ),
                                  getValue<std::string>( { "settings", "interface", "color_high" } ) );

    double fBase = getValue<double>( { "settings", "normalization", "log_base", "val" } );

    if( iBetweenGroupSetting == 2 )
        sBackgroundColor = vColorPalette[ colorIndex( logScale( 0.5, fBase ) ) ];
    else
        sBackgroundColor = vColorPalette.front( );
    END_RETURN;
}

bool PartialQuarry::setAnnotationColors( )
{
    vColorPaletteAnnotation.clear( );
    auto rList = getValue<json>( { "settings", "interface", "annotation_color_palette" } );
    vColorPaletteAnnotation.reserve( rList.size( ) );
    for( auto& rC : rList )
    {
        CANCEL_RETURN;
        vColorPaletteAnnotation.push_back( rC.get<std::string>( ) );
    }
    vColorPaletteAnnotationDark.clear( );
    auto rList2 = getValue<json>( { "settings", "interface", "annotation_color_palette_dark" } );
    vColorPaletteAnnotationDark.reserve( rList2.size( ) );
    for( auto& rC : rList2 )
    {
        CANCEL_RETURN;
        vColorPaletteAnnotationDark.push_back( rC.get<std::string>( ) );
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

    if( getValue<std::string>( { "settings", "normalization", "scale" } ) == "minmax" )
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
    else if( getValue<std::string>( { "settings", "normalization", "scale" } ) == "abs" )
    {
        double fAbs = std::max( std::abs( fMax ), std::abs( fMin ) );
        for( double fVal : vDivided )
        {
            CANCEL_RETURN;
            if( std::isnan( fVal ) )
                vScaled.push_back( fVal );
            else
                vScaled.push_back( fVal / fAbs );
        }
        fMin = std::max( -1.0, fMin );
        fMax = std::min( 1.0, fMin );
    }
    else if( getValue<std::string>( { "settings", "normalization", "scale" } ) == "max" )
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
    else if( getValue<std::string>( { "settings", "normalization", "scale" } ) == "dont" )
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
    bool bLargerZero = fX > 0.5;
    if( iBetweenGroupSetting == 2 ) // between_groups == sub
        fX = std::abs( fX - 0.5 );

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
    vRanged.clear( );
    vRanged.reserve( vScaled.size( ) );

    double fMin = getValue<double>( { "settings", "normalization", "color_range", "val_min" } );
    double fMax = getValue<double>( { "settings", "normalization", "color_range", "val_max" } );
    double fBase = getValue<double>( { "settings", "normalization", "log_base", "val" } );

    fDataMin = std::numeric_limits<double>::max( );
    fDataMax = std::numeric_limits<double>::min( );

    for( double fC : vScaled )
    {
        CANCEL_RETURN;
        if( std::isnan( fC ) )
            vColored.push_back( sBackgroundColor );
        else
        {
            fDataMin = std::min( fDataMin, fC );
            fDataMax = std::max( fDataMax, fC );

            vRanged.push_back( colorRange( fMin, fMax, fC ) );
            vColored.push_back( vColorPalette[ colorIndex( logScale( vRanged.back( ), fBase ) ) ] );
        }
    }
    END_RETURN;
}

bool PartialQuarry::setPalette( )
{
    pybind11::gil_scoped_acquire acquire;
    pybind11::list vNew;
    double fMinR = getValue<double>( { "settings", "normalization", "color_range", "val_min" } );
    double fMaxR = getValue<double>( { "settings", "normalization", "color_range", "val_max" } );
    double fBase = getValue<double>( { "settings", "normalization", "log_base", "val" } );

    for( double fC = fMinR; fC <= fMaxR; fC += std::max( 1.0, fMaxR - fMinR ) / 255.0 )
    {
        CANCEL_RETURN;
        vNew.append( vColorPalette[ colorIndex( logScale( colorRange( fMinR, fMaxR, fC ), fBase ) ) ] );
    }

    vRenderedPalette = vNew;
    END_RETURN;
}

const pybind11::list PartialQuarry::getPalette( const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::Palette, fPyPrint );
    return vRenderedPalette;
}

// min_data max_data min_range max_range
const std::array<double, 4> PartialQuarry::getPaletteTicks( const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::Palette, fPyPrint );
    double fMinR = getValue<double>( { "settings", "normalization", "color_range", "val_min" } );
    double fMaxR = getValue<double>( { "settings", "normalization", "color_range", "val_max" } );

    double fCRMin = std::min( fDataMin, fMinR );
    double fCRMax = std::max( fDataMax, fMaxR );

    return std::array<double, 4>{ ( fDataMin - fCRMin ) / ( fCRMax - fCRMin ),
                                  ( fDataMax - fCRMin ) / ( fCRMax - fCRMin ), ( fMinR - fCRMin ) / ( fCRMax - fCRMin ),
                                  ( fMaxR - fCRMin ) / ( fCRMax - fCRMin ) };
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
    pybind11::list vRangedOut;
    pybind11::list vScoreA;
    pybind11::list vScoreB;
    pybind11::list vZero;

    size_t uiDividend = getValue<size_t>( { "dividend" } );

    for( size_t uiI = 0; uiI < vColored.size( ); uiI++ )
    {
        CANCEL_RETURN;

        if( vColored[ uiI ] != sBackgroundColor || true )
        {
            vScreenBottom.append( vBinCoords[ uiI ][ 0 ].uiScreenY );
            vScreenLeft.append( vBinCoords[ uiI ][ 0 ].uiScreenX );
            vScreenTop.append( vBinCoords[ uiI ][ 0 ].uiScreenY + vBinCoords[ uiI ][ 0 ].uiScreenH );
            vScreenRight.append( vBinCoords[ uiI ][ 0 ].uiScreenX + vBinCoords[ uiI ][ 0 ].uiScreenW );

            vColor.append( vColored[ uiI ] );


            std::string sChromNameX =
                substringChr( vActiveChromosomes[ 0 ][ vBinCoords[ uiI ][ 0 ].uiChromosomeX ].sName );
            std::string sChromNameY =
                substringChr( vActiveChromosomes[ 1 ][ vBinCoords[ uiI ][ 0 ].uiChromosomeY ].sName );
            vChrX.append( sChromNameX );
            vChrY.append( sChromNameY );
            vIndexLeft.append( readableBp( vBinCoords[ uiI ][ 0 ].uiIndexX * uiDividend ) );
            vIndexRight.append(
                readableBp( ( vBinCoords[ uiI ][ 0 ].uiIndexX + vBinCoords[ uiI ][ 0 ].uiIndexW ) * uiDividend ) );
            vIndexBottom.append( readableBp( vBinCoords[ uiI ][ 0 ].uiIndexY * uiDividend ) );
            vIndexTop.append(
                readableBp( ( vBinCoords[ uiI ][ 0 ].uiIndexY + vBinCoords[ uiI ][ 0 ].uiIndexH ) * uiDividend ) );

            if( vBinCoords[ uiI ][ 1 ].uiChromosomeX != std::numeric_limits<size_t>::max( ) )
            {
                std::string sChromNameX =
                    substringChr( vActiveChromosomes[ 0 ][ vBinCoords[ uiI ][ 1 ].uiChromosomeX ].sName );
                std::string sChromNameY =
                    substringChr( vActiveChromosomes[ 1 ][ vBinCoords[ uiI ][ 1 ].uiChromosomeY ].sName );
                vChrXSym.append( sChromNameX );
                vChrYSym.append( sChromNameY );
                vIndexSymLeft.append( readableBp( vBinCoords[ uiI ][ 1 ].uiIndexX * uiDividend ) );
                vIndexSymRight.append(
                    readableBp( ( vBinCoords[ uiI ][ 1 ].uiIndexX + vBinCoords[ uiI ][ 1 ].uiIndexW ) * uiDividend ) );
                vIndexSymBottom.append( readableBp( vBinCoords[ uiI ][ 1 ].uiIndexY * uiDividend ) );
                vIndexSymTop.append(
                    readableBp( ( vBinCoords[ uiI ][ 1 ].uiIndexY + vBinCoords[ uiI ][ 1 ].uiIndexH ) * uiDividend ) );
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

            vScoreTotal.append( vScaled[ uiI ] );
            vRangedOut.append( vRanged[ uiI ] );
            vScoreA.append( vvNormalized[ uiI ][ 0 ] );
            vScoreB.append( vvNormalized[ uiI ][ 1 ] );
            vZero.append( 0 );
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
                                  "ranged_score"_a = vRangedOut,
                                  "score_a"_a = vScoreA,
                                  "score_b"_a = vScoreB,
                                  "0"_a = vZero );
    END_RETURN;
}

bool PartialQuarry::setHeatmapExport( )
{
    vHeatmapExport.clear( );
    vHeatmapExport.reserve( vScaled.size( ) );
    size_t uiDividend = getValue<size_t>( { "dividend" } );
    for( size_t uiI = 0; uiI < vScaled.size( ); uiI++ )
    {
        CANCEL_RETURN;
        std::string sChromNameX = vActiveChromosomes[ 0 ][ vBinCoords[ uiI ][ 0 ].uiChromosomeX ].sName;
        std::string sChromNameY = vActiveChromosomes[ 1 ][ vBinCoords[ uiI ][ 0 ].uiChromosomeY ].sName;
        vHeatmapExport.emplace_back( sChromNameX,
                                     vBinCoords[ uiI ][ 0 ].uiIndexX * uiDividend,
                                     ( vBinCoords[ uiI ][ 0 ].uiIndexX + vBinCoords[ uiI ][ 0 ].uiIndexW ) * uiDividend,
                                     sChromNameY,
                                     vBinCoords[ uiI ][ 0 ].uiIndexY * uiDividend,
                                     ( vBinCoords[ uiI ][ 0 ].uiIndexY + vBinCoords[ uiI ][ 0 ].uiIndexH ) * uiDividend,
                                     vScaled[ uiI ] );
    }
    END_RETURN;
}

const pybind11::dict PartialQuarry::getHeatmap( const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::HeatmapCDS, fPyPrint );
    return xHeatmapCDS;
}

const decltype( PartialQuarry::vHeatmapExport )
PartialQuarry::getHeatmapExport( const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::HeatmapExport, fPyPrint );
    return vHeatmapExport;
}

const std::string& PartialQuarry::getBackgroundColor( const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::Colors, fPyPrint );
    return sBackgroundColor;
}

void PartialQuarry::regColors( )
{
    registerNode( NodeNames::Colors,
                  ComputeNode{ /*.sNodeName =*/"color_palette",
                               /*.fFunc =*/&PartialQuarry::setColors,
                               /*.vIncomingFunctions =*/{ NodeNames::BetweenGroup },
                               /*.vIncomingSession =*/
                               { { "settings", "interface", "color_palette" },
                                 { "settings", "interface", "color_low" },
                                 { "settings", "interface", "color_high" },
                                 { "settings", "normalization", "log_base", "val" } },
                               /*.vSessionsIncomingInPrevious =*/{} } );

    registerNode( NodeNames::AnnotationColors,
                  ComputeNode{ /*.sNodeName =*/"annotation_color_palette",
                               /*.fFunc =*/&PartialQuarry::setAnnotationColors,
                               /*.vIncomingFunctions =*/{ },
                               /*.vIncomingSession =*/
                               { { "settings", "interface", "annotation_color_palette" },
                                 { "settings", "interface", "annotation_color_palette_dark" } },
                               /*.vSessionsIncomingInPrevious =*/{} } );

    registerNode( NodeNames::Combined,
                  ComputeNode{ /*.sNodeName =*/"combined_bins",
                               /*.fFunc =*/&PartialQuarry::setCombined,
                               /*.vIncomingFunctions =*/{ NodeNames::DistDepDecayRemoved },
                               /*.vIncomingSession =*/{ },
                               /*.vSessionsIncomingInPrevious =*/{} } );

    registerNode( NodeNames::Scaled,
                  ComputeNode{ /*.sNodeName =*/"scaled_bins",
                               /*.fFunc =*/&PartialQuarry::setScaled,
                               /*.vIncomingFunctions =*/{ NodeNames::Divided },
                               /*.vIncomingSession =*/{ { "settings", "normalization", "scale" } },
                               /*.vSessionsIncomingInPrevious =*/{} } );

    registerNode(
        NodeNames::Colored,
        ComputeNode{ /*.sNodeName =*/"colored_bins",
                     /*.fFunc =*/&PartialQuarry::setColored,
                     /*.vIncomingFunctions =*/{ NodeNames::Colors, NodeNames::Scaled },
                     /*.vIncomingSession =*/
                     { { "settings", "normalization", "color_range", "val_max" },
                       { "settings", "normalization", "color_range", "val_min" } },
                     /*.vSessionsIncomingInPrevious =*/{ { "settings", "normalization", "log_base", "val" } } } );

    registerNode( NodeNames::HeatmapCDS,
                  ComputeNode{ /*.sNodeName =*/"heatmap_cds",
                               /*.fFunc =*/&PartialQuarry::setHeatmapCDS,
                               /*.vIncomingFunctions =*/{ NodeNames::Colored },
                               /*.vIncomingSession =*/{ },
                               /*.vSessionsIncomingInPrevious =*/{ { "dividend" } } } );

    registerNode( NodeNames::HeatmapExport,
                  ComputeNode{ /*.sNodeName =*/"heatmap_export",
                               /*.fFunc =*/&PartialQuarry::setHeatmapExport,
                               /*.vIncomingFunctions =*/{ NodeNames::Scaled },
                               /*.vIncomingSession =*/{ },
                               /*.vSessionsIncomingInPrevious =*/{ { "dividend" } } } );

    registerNode(
        NodeNames::Palette,
        ComputeNode{ /*.sNodeName =*/"rendered_palette",
                     /*.fFunc =*/&PartialQuarry::setPalette,
                     /*.vIncomingFunctions =*/{ NodeNames::Scaled, NodeNames::Colors },
                     /*.vIncomingSession =*/
                     { { "settings", "normalization", "color_range", "val_max" },
                       { "settings", "normalization", "color_range", "val_min" } },
                     /*.vSessionsIncomingInPrevious =*/{ { "settings", "normalization", "log_base", "val" } } } );
}

} // namespace cm