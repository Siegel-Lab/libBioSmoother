#include "cm/computation.h"
#include <cmath>

#pragma once

namespace cm
{

size_t ContactMapping::nextEvenNumber( double fX )
{
    if( this->xRenderSettings[ "interface" ][ "snap_bin_size" ] == "no" )
        return std::ceil( fX );

    size_t uiN = 0;
    while( true )
    {
        for( size_t uiI : this->xRenderSettings[ "interface" ][ "snap_factors" ] )
            if( std::pow( uiI * 10, uiN ) > fX )
                return std::pow( uiI * 10, uiN );
        uiN += 1;
    }
}


void ContactMapping::setBinSize( )
{
    size_t uiT = this->xRenderSettings[ "interface" ][ "min_bin_size" ][ "val" ].get<size_t>( );
    size_t uiMinBinSize = std::max( (size_t)1,
                                    (size_t)std::ceil( ( 1 + uiT % 9 ) * std::pow( 10, uiT / 9 ) ) /
                                        this->xSession[ "dividend" ].get<size_t>( ) );
    size_t uiMaxNumBins = this->xRenderSettings[ "interface" ][ "max_num_bins" ][ "val" ].get<size_t>( ) *
                          this->xRenderSettings[ "interface" ][ "max_num_bins_factor" ].get<size_t>( );

    if( this->xRenderSettings[ "interface" ][ "bin_aspect_ratio" ] == "view" )
    {
        uiBinHeight = nextEvenNumber( ( this->xSession[ "area" ][ "y_end" ].get<size_t>( ) -
                                        this->xSession[ "area" ][ "y_start" ].get<size_t>( ) ) /
                                      uiMaxNumBins );
        uiBinHeight = std::max( uiBinHeight, uiMinBinSize );

        uiBinWidth = nextEvenNumber( ( this->xSession[ "area" ][ "x_end" ].get<size_t>( ) -
                                       this->xSession[ "area" ][ "x_start" ].get<size_t>( ) ) /
                                     uiMaxNumBins );
        uiBinWidth = std::max( uiBinWidth, uiMinBinSize );
    }
    else if( this->xRenderSettings[ "interface" ][ "bin_aspect_ratio" ].get<std::string>( ) == "coord" )
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
}

void ContactMapping::setRenderArea( )
{
    if( this->xRenderSettings[ "export" ][ "do_export_full" ] )
    {
        iStartX = 0;
        iStartY = 0;

        iEndX = this->xSession[ "contigs" ][ "genome_size" ].get<size_t>( );
        iEndY = this->xSession[ "contigs" ][ "genome_size" ].get<size_t>( );
    }
    else
    {
        iStartX = this->xSession[ "area" ][ "x_start" ].get<size_t>( ) -
                  ( this->xSession[ "area" ][ "x_start" ].get<size_t>( ) % uiBinWidth );
        iStartY = this->xSession[ "area" ][ "y_start" ].get<size_t>( ) -
                  ( this->xSession[ "area" ][ "y_start" ].get<size_t>( ) % uiBinHeight );

        iEndX = this->xSession[ "area" ][ "x_end" ].get<size_t>( ) + uiBinWidth -
                ( this->xSession[ "area" ][ "x_end" ].get<size_t>( ) % uiBinWidth );
        iEndY = this->xSession[ "area" ][ "y_end" ].get<size_t>( ) + uiBinHeight -
                ( this->xSession[ "area" ][ "y_end" ].get<size_t>( ) % uiBinHeight );
    }
}

void ContactMapping::regBinSize( )
{
    registerNode( NodeNames::BinSize,
                  ComputeNode{ .sNodeName = "bin_size",
                               .fFunc = &ContactMapping::setBinSize,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "interface", "snap_bin_size" },
                                                     { "interface", "snap_factors" },
                                                     { "interface", "min_bin_size", "val" },
                                                     { "interface", "max_num_bins", "val" },
                                                     { "interface", "max_num_bins_factor" },
                                                     { "interface", "bin_aspect_ratio" } },
                               .vIncomingRender = { { "dividend" },
                                                    { "area", "y_end" },
                                                    { "area", "y_start" },
                                                    { "area", "x_end" },
                                                    { "area", "x_start" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::RenderArea,
                  ComputeNode{ .sNodeName = "render_area",
                               .fFunc = &ContactMapping::setRenderArea,
                               .vIncomingFunctions = { NodeNames::BinSize },
                               .vIncomingSession = { { "export", "do_export_full" } },
                               .vIncomingRender = { { "contigs", "genome_size" },
                                                    { "area", "y_end" },
                                                    { "area", "y_start" },
                                                    { "area", "x_end" },
                                                    { "area", "x_start" } },
                               .uiLastUpdated = uiCurrTime } );
}


} // namespace cm