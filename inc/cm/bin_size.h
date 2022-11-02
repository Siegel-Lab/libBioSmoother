#include "cm/partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

size_t PartialQuarry::nextEvenNumber( double fX )
{
    if( this->xSession[ "settings" ][ "interface" ][ "snap_bin_size" ] == "no" )
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
    if( this->xSession[ "settings" ][ "interface" ][ "bin_aspect_ratio" ] == "view" )
    {
        uiBinHeight = nextEvenNumber( ( this->xSession[ "area" ][ "y_end" ].get<double>( ) -
                                        this->xSession[ "area" ][ "y_start" ].get<double>( ) ) /
                                      std::sqrt( uiMaxNumBins ) );
        uiBinHeight = std::max( uiBinHeight, uiMinBinSize );

        uiBinWidth = nextEvenNumber( ( this->xSession[ "area" ][ "x_end" ].get<double>( ) -
                                       this->xSession[ "area" ][ "x_start" ].get<double>( ) ) /
                                     std::sqrt( uiMaxNumBins ) );
        uiBinWidth = std::max( uiBinWidth, uiMinBinSize );

        std::cout << "uiBinHeight: " << uiBinHeight << " uiBinWidth: " << uiBinWidth << std::endl;
    }
    else if( this->xSession[ "settings" ][ "interface" ][ "bin_aspect_ratio" ].get<std::string>( ) == "coord" )
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
    else
        throw std::logic_error( "invlaid bin_aspect_ratio value" );
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
        iStartX = this->xSession[ "area" ][ "x_start" ].get<int64_t>( ) -
                  ( this->xSession[ "area" ][ "x_start" ].get<int64_t>( ) % uiBinWidth );
        iStartY = this->xSession[ "area" ][ "y_start" ].get<int64_t>( ) -
                  ( this->xSession[ "area" ][ "y_start" ].get<int64_t>( ) % uiBinHeight );

        iEndX = this->xSession[ "area" ][ "x_end" ].get<int64_t>( ) + uiBinWidth -
                ( this->xSession[ "area" ][ "x_end" ].get<int64_t>( ) % uiBinWidth );
        iEndY = this->xSession[ "area" ][ "y_end" ].get<int64_t>( ) + uiBinHeight -
                ( this->xSession[ "area" ][ "y_end" ].get<int64_t>( ) % uiBinHeight );

    }
    END_RETURN;
}

const std::array<int64_t, 4> PartialQuarry::getDrawingArea( )
{
    update( NodeNames::RenderArea );
    return std::array<int64_t, 4>{ { iStartX, iStartY, iEndX, iEndY } };
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
                                                     { "settings", "interface", "bin_aspect_ratio" },
                                                     { "dividend" },
                                                     { "area" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::RenderArea,
                  ComputeNode{ .sNodeName = "render_area",
                               .fFunc = &PartialQuarry::setRenderArea,
                               .vIncomingFunctions = { NodeNames::BinSize },
                               .vIncomingSession = { { "settings", "export", "do_export_full" },
                                                     { "settings", "contigs", "genome_size" } },
                               .uiLastUpdated = uiCurrTime } );
}


} // namespace cm