#include "cm/computation.h"
#include <cmath>

#pragma once

namespace cm
{

size_t Computation::nextEvenNumber( double fX )
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


void Computation::setBinSize( )
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
        throw std::runtime_error( "invlaid bin_aspect_ratio value" );
}

void Computation::setRenderArea( )
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


} // namespace cm