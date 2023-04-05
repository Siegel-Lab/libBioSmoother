#include "cm/partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

size_t PartialQuarry::nextEvenNumber( double fX )
{
    if( !getValue<bool>( { "settings", "interface", "snap_bin_size" } ) )
        return std::ceil( fX );

    size_t uiN = 0;
    while( true )
    {
        for( size_t uiI : getValue<json>( { "settings", "interface", "snap_factors" } ) )
            if( uiI * std::pow( 10, uiN ) > fX )
                return uiI * std::pow( 10, uiN );
        uiN += 1;
    }
}


bool PartialQuarry::setBinSize( )
{
    size_t uiT = getValue<size_t>( { "settings", "interface", "min_bin_size", "val" } );
    size_t uiMinBinSize =
        std::max( (size_t)1,
                  (size_t)std::ceil( ( 1 + uiT % 9 ) * std::pow( 10, uiT / 9 ) ) / getValue<size_t>( { "dividend" } ) );
    size_t uiMaxNumBins = getValue<size_t>( { "settings", "interface", "max_num_bins", "val" } ) *
                          getValue<size_t>( { "settings", "interface", "max_num_bins_factor" } );
    if( getValue<bool>( { "settings", "interface", "fixed_number_of_bins" } ) )
    {
        size_t uiNumX = getValue<size_t>( { "settings", "interface", "fixed_num_bins_x", "val" } );
        size_t uiNumY = getValue<size_t>( { "settings", "interface", "fixed_num_bins_y", "val" } );
        uiBinHeight.w( ) = std::max(
            (size_t)1,
            (size_t)( ( getValue<double>( { "area", "y_end" } ) - getValue<double>( { "area", "y_start" } ) ) /
                      static_cast<double>( uiNumY ) ) );
        uiBinWidth.w( ) = std::max(
            (size_t)1,
            (size_t)( ( getValue<double>( { "area", "x_end" } ) - getValue<double>( { "area", "x_start" } ) ) /
                      static_cast<double>( uiNumX ) ) );
    }
    else if( getValue<bool>( { "settings", "interface", "fixed_bin_size" } ) )
    {
        size_t uiNumX = getValue<size_t>( { "settings", "interface", "fixed_bin_size_x", "val" } );
        size_t uiNumY = getValue<size_t>( { "settings", "interface", "fixed_bin_size_y", "val" } );
        uiBinHeight.w( ) = std::max( (size_t)1, (size_t)uiNumY );
        uiBinWidth.w( ) = std::max( (size_t)1, (size_t)uiNumX );
    }
    else if( !getValue<bool>( { "settings", "interface", "squared_bins" } ) )
    {
        uiBinHeight.w( ) =
            nextEvenNumber( ( getValue<double>( { "area", "y_end" } ) - getValue<double>( { "area", "y_start" } ) ) /
                            std::sqrt( uiMaxNumBins ) );
        uiBinHeight.w( ) = std::max( uiBinHeight.r( ), uiMinBinSize );

        uiBinWidth.w( ) =
            nextEvenNumber( ( getValue<double>( { "area", "x_end" } ) - getValue<double>( { "area", "x_start" } ) ) /
                            std::sqrt( uiMaxNumBins ) );
        uiBinWidth.w( ) = std::max( uiBinWidth.r( ), uiMinBinSize );
    }
    else
    {
        size_t uiArea = ( getValue<size_t>( { "area", "x_end" } ) - getValue<size_t>( { "area", "x_start" } ) ) *
                        ( getValue<size_t>( { "area", "y_end" } ) - getValue<size_t>( { "area", "y_start" } ) ) /
                        uiMaxNumBins;
        uiBinHeight.w( ) = nextEvenNumber( std::sqrt( uiArea ) );
        uiBinHeight.w( ) = std::max( uiBinHeight.r( ), uiMinBinSize );

        uiBinWidth.w( ) = uiBinHeight.r( );
    }
    END_RETURN;
}

bool PartialQuarry::setRenderArea( )
{
    if( getValue<bool>( { "settings", "export", "do_export_full" } ) )
    {
        iStartX = 0;
        iStartY = 0;

        iEndX = getValue<size_t>( { "contigs", "genome_size" } );
        iEndY = getValue<size_t>( { "contigs", "genome_size" } );
    }
    else
    {
        double fScale = getValue<double>( { "settings", "interface", "add_draw_area", "val" } ) / 100.0;

        int64_t uiL = getValue<size_t>( { "area", "x_start" } );
        int64_t uiR = getValue<size_t>( { "area", "x_end" } );

        int64_t uiW = uiR - uiL;

        int64_t uiB = getValue<size_t>( { "area", "y_start" } );
        int64_t uiT = getValue<size_t>( { "area", "y_end" } );

        int64_t uiH = uiT - uiB;

        uiL -= uiW * fScale;
        uiR += uiW * fScale;

        uiB -= uiH * fScale;
        uiT += uiH * fScale;

        iStartX = uiL - ( uiL % uiBinWidth.r( ) );
        iStartY = uiB - ( uiB % uiBinHeight.r( ) );

        iEndX = uiR + uiBinWidth.r( ) - ( uiR % uiBinWidth.r( ) );
        iEndY = uiT + uiBinHeight.r( ) - ( uiT % uiBinHeight.r( ) );
    }
    END_RETURN;
}

const std::array<int64_t, 4> PartialQuarry::getDrawingArea( const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::RenderArea, fPyPrint );
    return std::array<int64_t, 4>{ { iStartX, iStartY, iEndX, iEndY } };
}


const std::array<size_t, 2> PartialQuarry::getBinSize( const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::BinSize, fPyPrint );
    size_t uiD = getValue<size_t>( { "dividend" } );
    return std::array<size_t, 2>{ uiBinWidth.r( ) * uiD, uiBinHeight.r( ) * uiD };
}


template <typename CharT> struct Sep : public std::numpunct<CharT>
{
    virtual std::string do_grouping( ) const
    {
        return "\003";
    }
};

std::string putCommas( size_t uiBp )
{
    std::string sNum = std::to_string( uiBp );
    std::string sRet = "";
    for( size_t uiI = 0; uiI < sNum.size( ); uiI++ )
    {
        if( uiI % 3 == 0 && uiI > 0 )
            sRet += ",";
        sRet += sNum[ sNum.size( ) - uiI - 1 ];
    }
    std::reverse( sRet.begin( ), sRet.end( ) );
    return sRet;
}

std::string PartialQuarry::readableBp( size_t uiBp )
{
    if( uiBp == 0 )
        return "0 bp";
    else if( uiBp % 1000000 == 0 )
        return putCommas( uiBp / 1000000 ) + " mbp";
    else if( uiBp % 1000 == 0 )
        return putCommas( uiBp / 1000 ) + " kbp";
    else
        return putCommas( uiBp ) + " bp";
}

void PartialQuarry::regBinSize( )
{
    registerNode( NodeNames::BinSize,
                  ComputeNode{ /*.sNodeName =*/"bin_size",
                               /*.fFunc =*/&PartialQuarry::setBinSize,
                               /*.vIncomingFunctions =*/{ },
                               /*.vIncomingSession =*/
                               { { "settings", "interface", "snap_bin_size" },
                                 { "settings", "interface", "snap_factors" },
                                 { "settings", "interface", "min_bin_size", "val" },
                                 { "settings", "interface", "max_num_bins", "val" },
                                 { "settings", "interface", "max_num_bins_factor" },
                                 { "settings", "interface", "squared_bins" },
                                 { "settings", "interface", "fixed_number_of_bins" },
                                 { "settings", "interface", "fixed_num_bins_x", "val" },
                                 { "settings", "interface", "fixed_num_bins_y", "val" },
                                 { "settings", "interface", "fixed_bin_size" },
                                 { "settings", "interface", "fixed_bin_size_x", "val" },
                                 { "settings", "interface", "fixed_bin_size_y", "val" },
                                 { "dividend" },
                                 { "area" } },
                               /*.vSessionsIncomingInPrevious =*/{ },
                               /*bHidden =*/false } );

    registerNode( NodeNames::RenderArea,
                  ComputeNode{ /*.sNodeName =*/"render_area",
                               /*.fFunc =*/&PartialQuarry::setRenderArea,
                               /*.vIncomingFunctions =*/{ NodeNames::BinSize },
                               /*.vIncomingSession =*/
                               {
                                   { "settings", "export", "do_export_full" },
                                   { "settings", "contigs", "genome_size" },
                                   { "settings", "interface", "add_draw_area", "val" },
                                   {"contigs", "genome_size"},
                               },
                               /*.vSessionsIncomingInPrevious =*/{ { "area" } },
                               /*bHidden =*/true } );
}


} // namespace cm