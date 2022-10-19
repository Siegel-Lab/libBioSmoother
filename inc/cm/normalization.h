#include "cm/partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

void PartialQuarry::normalizeSize( size_t uiSize )
{
    std::array<size_t, 2> vTotalReads = { 0, 0 };
    for( std::string& rRep : this->vActiveReplicates )
    {
        if( this->xSession[ "replicates" ][ "in_group" ][ rRep ][ "in_group_a" ].get<bool>( ) )
            vTotalReads[ 0 ] += this->xSession[ "replicates" ][ "by_name" ][ rRep ][ "total_reads" ].get<size_t>( );
        if( this->xSession[ "replicates" ][ "in_group" ][ rRep ][ "in_group_b" ].get<bool>( ) )
            vTotalReads[ 1 ] += this->xSession[ "replicates" ][ "by_name" ][ rRep ][ "total_reads" ].get<size_t>( );
    }
    for( auto& vVal : vvFlatValues )
        vvNormalized.push_back( std::array<double, 2>{ uiSize * vVal[ 0 ] / (double)vTotalReads[ 0 ],
                                                       uiSize * vVal[ 1 ] / (double)vTotalReads[ 1 ] } );
}

void PartialQuarry::doNotNormalize( )
{
    for( auto& vVal : vvFlatValues )
        vvNormalized.push_back( std::array<double, 2>{ (double)vVal[ 0 ], (double)vVal[ 1 ] } );
}


void PartialQuarry::normalizeTracks( )
{
    for( size_t uiI = 0; uiI < vvFlatValues.size( ); uiI++ )
        vvNormalized.push_back(
            std::array<double, 2>{ (double)vvFlatValues[ uiI ][ 0 ] /
                                       std::max( 1.0, vvFlatCoverageValues[ 0 ][ uiI / vAxisCords[ 0 ].size( ) ] ),
                                   (double)vvFlatValues[ uiI ][ 1 ] /
                                       std::max( 1.0, vvFlatCoverageValues[ 1 ][ uiI % vAxisCords[ 0 ].size( ) ] ) } );
}

void PartialQuarry::normalizeBinominalTest( )
{
    size_t uiNumBinsInRowTotal = ( this->xSession[ "contigs" ][ "genome_size" ].get<size_t>( ) - 1 ) / uiBinWidth + 1;
    vvNormalized = normalizeBinominalTestTrampoline(
        vvFlatValues, '?', uiNumBinsInRowTotal,
        this->xSession[ "settings" ][ "normalization" ][ "p_accept" ][ "val" ].get<double>( ) );
}

void PartialQuarry::normalizeIC( )
{
    throw std::logic_error( "Function not implemented" );
}

void PartialQuarry::setNormalized( )
{
    vvNormalized.reserve( vvFlatValues.size( ) );
    vvNormalized.clear( );

    if( this->xSession[ "settings" ][ "normalization" ][ "normalize_by" ].get<std::string>( ) == "dont" )
        doNotNormalize( );
    else if( this->xSession[ "settings" ][ "normalization" ][ "normalize_by" ].get<std::string>( ) == "tracks" )
        normalizeTracks( );
    else if( this->xSession[ "settings" ][ "normalization" ][ "normalize_by" ].get<std::string>( ) == "radicl-seq" )
        normalizeBinominalTest( );
    else if( this->xSession[ "settings" ][ "normalization" ][ "normalize_by" ].get<std::string>( ) == "rpm" )
        normalizeSize( 1000000 );
    else if( this->xSession[ "settings" ][ "normalization" ][ "normalize_by" ].get<std::string>( ) == "rpk" )
        normalizeSize( 1000 );
    else if( this->xSession[ "settings" ][ "normalization" ][ "normalize_by" ].get<std::string>( ) == "hi-c" )
        normalizeIC( );
    else
        throw std::logic_error( "invalid value for normalize_by" );
}

void PartialQuarry::regNormalization( )
{
    registerNode( NodeNames::Normalized,
                  ComputeNode{ .sNodeName = "normalization",
                               .fFunc = &PartialQuarry::setNormalized,
                               .vIncomingFunctions = { NodeNames::FlatValues, NodeNames::FlatCoverageValues },
                               .vIncomingSession = { { "normalization", "normalize_by" },
                                                     { "replicates", "in_group" },
                                                     { "settings", "normalization", "p_accept", "val" } },
                               .uiLastUpdated = uiCurrTime } );
}

} // namespace cm