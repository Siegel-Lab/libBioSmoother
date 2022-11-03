#include "cm/partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

size_t PartialQuarry::nextEvenNumber( double fX )
{
    if( !this->xSession[ "settings" ][ "interface" ][ "snap_bin_size" ].get<bool>( ) )
        return std::ceil( fX );

    size_t uiN = 0;
    while( true )
    {
        for( size_t uiI : this->xSession[ "settings" ][ "interface" ][ "snap_factors" ] )
            if( uiI * std::pow( 10, uiN ) > fX )
                return uiI * std::pow( 10, uiN );
        uiN += 1;
    }
}


bool PartialQuarry::setBinSize( )
{
    size_t uiT = this->xSession[ "settings" ][ "interface" ][ "min_bin_size" ][ "val" ].get<size_t>( );
    size_t uiMinBinSize = std::max( (size_t)1,
                                    (size_t)std::ceil( ( 1 + uiT % 9 ) * std::pow( 10, uiT / 9 ) ) /
                                        this->xSession[ "dividend" ].get<size_t>( ) );
    size_t uiMaxNumBins = this->xSession[ "settings" ][ "interface" ][ "max_num_bins" ][ "val" ].get<size_t>( ) *
                          this->xSession[ "settings" ][ "interface" ][ "max_num_bins_factor" ].get<size_t>( );
    if( !this->xSession[ "settings" ][ "interface" ][ "squared_bins" ].get<bool>( ) )
    {
        uiBinHeight = nextEvenNumber( ( this->xSession[ "area" ][ "y_end" ].get<double>( ) -
                                        this->xSession[ "area" ][ "y_start" ].get<double>( ) ) /
                                      std::sqrt( uiMaxNumBins ) );
        uiBinHeight = std::max( uiBinHeight, uiMinBinSize );

        uiBinWidth = nextEvenNumber( ( this->xSession[ "area" ][ "x_end" ].get<double>( ) -
                                       this->xSession[ "area" ][ "x_start" ].get<double>( ) ) /
                                     std::sqrt( uiMaxNumBins ) );
        uiBinWidth = std::max( uiBinWidth, uiMinBinSize );
    }
    else
    {
        size_t uiArea = ( this->xSession[ "area" ][ "x_end" ].get<size_t>( ) -
                          this->xSession[ "area" ][ "x_start" ].get<size_t>( ) ) *
                        ( this->xSession[ "area" ][ "y_end" ].get<size_t>( ) -
                          this->xSession[ "area" ][ "y_start" ].get<size_t>( ) ) /
                        uiMaxNumBins;
        uiBinHeight = nextEvenNumber( std::sqrt( uiArea ) );
        uiBinHeight = std::max( uiBinHeight, uiMinBinSize );

        uiBinWidth = uiBinHeight;
    }
    END_RETURN;
}

bool PartialQuarry::setRenderArea( )
{
    if( this->xSession[ "settings" ][ "export" ][ "do_export_full" ] )
    {
        iStartX = 0;
        iStartY = 0;

        iEndX = this->xSession[ "contigs" ][ "genome_size" ].get<int64_t>( );
        iEndY = this->xSession[ "contigs" ][ "genome_size" ].get<int64_t>( );
    }
    else
    {
        double fScale = this->xSession[ "settings" ][ "interface" ][ "add_draw_area" ][ "val" ].get<double>( ) / 100.0;

        int64_t uiL = this->xSession[ "area" ][ "x_start" ].get<int64_t>( );
        int64_t uiR = this->xSession[ "area" ][ "x_end" ].get<int64_t>( );

        int64_t uiW = uiR - uiL;

        int64_t uiB = this->xSession[ "area" ][ "y_start" ].get<int64_t>( );
        int64_t uiT = this->xSession[ "area" ][ "y_end" ].get<int64_t>( );

        int64_t uiH = uiT - uiB;

        uiL -= uiW * fScale;
        uiR += uiW * fScale;

        uiB -= uiH * fScale;
        uiT += uiH * fScale;

        iStartX = uiL - ( uiL % uiBinWidth );
        iStartY = uiB - ( uiB % uiBinHeight );

        iEndX = uiR + uiBinWidth - ( uiR % uiBinWidth );
        iEndY = uiT + uiBinHeight - ( uiT % uiBinHeight );
    }
    END_RETURN;
}

const std::array<int64_t, 4> PartialQuarry::getDrawingArea( )
{
    update( NodeNames::RenderArea );
    return std::array<int64_t, 4>{ { iStartX, iStartY, iEndX, iEndY } };
}


const std::array<size_t, 2> PartialQuarry::getBinSize( )
{
    update( NodeNames::BinSize );
    size_t uiD = this->xSession[ "dividend" ].get<size_t>( );
    return std::array<size_t, 2>{ uiBinWidth * uiD, uiBinHeight * uiD };
}


template<typename CharT>
struct Sep : public std::numpunct<CharT>
{
    virtual std::string do_grouping() const {return "\003";}
};

std::string putCommas(size_t uiBp)
{
    std::stringstream ss;
    ss.imbue(std::locale(std::cout.getloc(), new Sep <char>()));
    ss << uiBp;
    return ss.str();
}

std::string PartialQuarry::readableBp( size_t uiBp )
{
    
    if( uiBp % 1000000 == 0 )
        return putCommas( uiBp / 1000000 ) + "mbp";
    else if( uiBp % 1000 == 0 )
        return putCommas( uiBp / 1000 ) + "kbp";
    else
        return putCommas(uiBp) + "bp";
}

void PartialQuarry::regBinSize( )
{
    registerNode( NodeNames::BinSize,
                  ComputeNode{ .sNodeName = "bin_size",
                               .fFunc = &PartialQuarry::setBinSize,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "settings", "interface", "snap_bin_size" },
                                                     { "settings", "interface", "snap_factors" },
                                                     { "settings", "interface", "min_bin_size", "val" },
                                                     { "settings", "interface", "max_num_bins", "val" },
                                                     { "settings", "interface", "max_num_bins_factor" },
                                                     { "settings", "interface", "squared_bins" },
                                                     { "dividend" },
                                                     { "area" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::RenderArea,
                  ComputeNode{ .sNodeName = "render_area",
                               .fFunc = &PartialQuarry::setRenderArea,
                               .vIncomingFunctions = { NodeNames::BinSize },
                               .vIncomingSession =
                                   {
                                       { "settings", "export", "do_export_full" },
                                       { "settings", "contigs", "genome_size" },
                                       { "settings", "interface", "add_draw_area", "val" },
                                   },
                               .uiLastUpdated = uiCurrTime } );
}


} // namespace cm