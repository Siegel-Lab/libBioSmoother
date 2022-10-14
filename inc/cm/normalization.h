#include "cm/computation.h"
#include <cmath>

#pragma once

namespace cm
{

void ContactMapping::normalizeSize( size_t uiSize )
{
    std::array<size_t, 2> vTotalReads = { 0, 0 };
    for( std::string& rRep : this->vActiveReplicates )
    {
        if( this->xSession[ "replicates" ][ "by_name" ][ rRep ][ "in_group_a" ].get<bool>( ) )
            vTotalReads[ 0 ] += this->xSession[ "replicates" ][ "by_name" ][ rRep ][ "total_reads" ].get<size_t>( );
        if( this->xSession[ "replicates" ][ "by_name" ][ rRep ][ "in_group_b" ].get<bool>( ) )
            vTotalReads[ 1 ] += this->xSession[ "replicates" ][ "by_name" ][ rRep ][ "total_reads" ].get<size_t>( );
    }
    for( auto& vVal : vvFlatValues )
        vvNormalized.push_back( std::array<double, 2>{ uiSize * vVal[ 0 ] / (double)vTotalReads[ 0 ],
                                                       uiSize * vVal[ 1 ] / (double)vTotalReads[ 1 ] } );
}

void ContactMapping::doNotNormalize( )
{
    for( auto& vVal : vvFlatValues )
        vvNormalized.push_back( std::array<double, 2>{ (double)vVal[ 0 ], (double)vVal[ 1 ] } );
}

/*

            elif self.settings['normalization']['normalize_by'] == "column":
                ret.append([x / max(ns[idx][idx_2 // len(rows)], 1) for idx_2, x in enumerate(bins)])
            elif self.settings['normalization']['normalize_by'] == "row":
                ret.append([x / max(ns[idx][idx_2 % len(rows)], 1) for idx_2, x in enumerate(bins)])

*/

void ContactMapping::normalizeTracks( )
{
    for( size_t uiI = 0; uiI < vvFlatValues.size(); uiI++ )
        vvNormalized.push_back( std::array<double, 2>{ (double)vvFlatValues[uiI][ 0 ], 
                                                       (double)vvFlatValues[uiI][ 1 ] } );
}

void ContactMapping::normalize( )
{
    vvNormalized.reserve( vvFlatValues.size( ) );
    if( this->xRenderSettings[ "normalization" ][ "normalize_by" ].get<std::string>( ) == "dont" )
        doNotNormalize( );
    else if( this->xRenderSettings[ "normalization" ][ "normalize_by" ].get<std::string>( ) == "tracks" )
        normalizeTracks();
    else if( this->xRenderSettings[ "normalization" ][ "normalize_by" ].get<std::string>( ) == "radicl-seq" )
        ;
    else if( this->xRenderSettings[ "normalization" ][ "normalize_by" ].get<std::string>( ) == "rpm" )
        normalizeSize( 1000000 );
    else if( this->xRenderSettings[ "normalization" ][ "normalize_by" ].get<std::string>( ) == "rpk" )
        normalizeSize( 1000 );
    else if( this->xRenderSettings[ "normalization" ][ "normalize_by" ].get<std::string>( ) == "hi-c" )
        ;
    else
        throw std::runtime_error( "invalid value for normalize_by" );
}

} // namespace cm