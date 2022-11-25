#include "cm/partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

bool PartialQuarry::normalizeSize( size_t uiSize )
{
    for( auto& vVal : vvFlatValues )
    {
        CANCEL_RETURN;
        vvNormalized.push_back( std::array<double, 2>{
            vvFlatTotal[ 0 ] == 0 ? 0 : (double)uiSize * (double)vVal[ 0 ] / (double)vvFlatTotal[ 0 ],
            vvFlatTotal[ 1 ] == 0 ? 0 : (double)uiSize * (double)vVal[ 1 ] / (double)vvFlatTotal[ 1 ] } );
    }
    END_RETURN;
}

bool PartialQuarry::doNotNormalize( )
{
    for( auto& vVal : vvFlatValues )
    {
        CANCEL_RETURN;
        vvNormalized.push_back( std::array<double, 2>{ (double)vVal[ 0 ], (double)vVal[ 1 ] } );
    }
    END_RETURN;
}


bool PartialQuarry::normalizeBinominalTest( )
{
    size_t uiNumBinsInRowTotal = ( getValue<size_t>( { "contigs", "genome_size" } ) - 1 ) / uiBinWidth + 1;
    vvNormalized =
        normalizeBinominalTestTrampoline( vvFlatValues, vFlatNormValues[ 1 ], uiNumBinsInRowTotal,
                                          getValue<double>( { "settings", "normalization", "p_accept", "val" } ) );
    CANCEL_RETURN;
    END_RETURN;
}

/*
original ICing implementation:
@ https://github.com/open2c/cooler/blob/3d284485df070255f4a904102178b186c14c1cee/cooler/balance.py
@ https://github.com/open2c/cooler/blob/3d284485df070255f4a904102178b186c14c1cee/cooler/tools.py

*/


size_t PartialQuarry::iceGetCount( IceData& rIceData, size_t uiX, size_t uiY, bool bA )
{
    assert( uiX < rIceData.vSliceBias[ 0 ].size( ) );
    assert( uiY < rIceData.vSliceBias[ 1 ].size( ) );
    size_t uiIdx = uiY + uiX * ( rIceData.vSliceBias[ 1 ].size( ) );
    assert( uiIdx < vvFlatValues.size( ) );
    return vvFlatValues[ uiIdx ][ bA ? 0 : 1 ];
}

void PartialQuarry::iceFilter( IceData& /*rIceData*/, size_t /*uiFrom*/, size_t /*uiTo*/ )
{}

void PartialQuarry::iceTimesOuterProduct( IceData& rIceData, bool bA, size_t uiFrom, size_t uiTo )
{
    const size_t uiH = rIceData.vSliceBias[ 1 ].size( );
    for( size_t uiI = uiFrom; uiI < uiTo; uiI++ )
    {
        assert( uiI / uiH < rIceData.vSliceBias[ 0 ].size( ) );
        assert( uiI % uiH < rIceData.vSliceBias[ 1 ].size( ) );
        assert( uiI < rIceData.vBiases.size( ) );
        rIceData.vBiases[ uiI ] = rIceData.vSliceBias[ 0 ][ uiI / uiH ] * rIceData.vSliceBias[ 1 ][ uiI % uiH ] *
                                  iceGetCount( rIceData, uiI / uiH, uiI % uiH, bA );
    }
}

void PartialQuarry::iceMarginalize( IceData& rIceData, bool bCol, size_t uiFrom, size_t uiTo )
{
    const size_t uiH = rIceData.vSliceBias[ 1 ].size( );
    const size_t uiWSlice = rIceData.vSliceBias[ bCol ? 1 : 0 ].size( );
    for( size_t uiI = uiFrom; uiI < uiTo; uiI++ )
    {
        assert( uiI < rIceData.vSliceMargin[ bCol ? 0 : 1 ].size( ) );
        rIceData.vSliceMargin[ bCol ? 0 : 1 ][ uiI ] = 0;
        for( size_t uiJ = 0; uiJ < uiWSlice; uiJ++ )
        {
            size_t uiK = bCol ? uiH * uiI + uiJ : uiI + uiH * uiJ;
            assert( uiK < rIceData.vBiases.size( ) );

            rIceData.vSliceMargin[ bCol ? 0 : 1 ][ uiI ] += rIceData.vBiases[ uiK ];
        }
    }
}

double PartialQuarry::iceNonZeroMarginMean( IceData& rIceData, bool bCol )
{
    double fMean = 0;
    size_t uiNum = 0;
    for( double fVal : rIceData.vSliceMargin[ bCol ? 0 : 1 ] )
        if( fVal != 0.0 )
        {
            fMean += fVal;
            ++uiNum;
        }
    if( uiNum == 0 )
        return 0;
    return fMean / (double)uiNum;
}

double PartialQuarry::iceNonZeroMarginVariance( IceData& rIceData, bool bCol, double fMean )
{
    double fSum = 0;
    size_t uiNum = 0;
    for( double fVal : rIceData.vSliceMargin[ bCol ? 0 : 1 ] )
        if( fVal != 0.0 )
        {
            fSum += ( fVal - fMean ) * ( fVal - fMean );
            ++uiNum;
        }
    if( uiNum == 0 )
        return 0;
    return fSum / (double)uiNum;
}

void PartialQuarry::iceDivByMargin( IceData& rIceData, bool bCol, double fMean, size_t uiFrom, size_t uiTo )
{
    for( size_t uiI = uiFrom; uiI < uiTo; uiI++ )
    {
        assert( uiI < rIceData.vSliceMargin[ bCol ? 0 : 1 ].size( ) );

        double fVal = rIceData.vSliceMargin[ bCol ? 0 : 1 ][ uiI ];
        fVal /= fMean;
        if( fVal == 0 )
            fVal = 1;

        assert( uiI < rIceData.vSliceBias[ bCol ? 0 : 1 ].size( ) );

        rIceData.vSliceBias[ bCol ? 0 : 1 ][ uiI ] /= fVal;
    }
}


void PartialQuarry::iceApplyBias( IceData& rIceData, bool bA, size_t uiFrom, size_t uiTo )
{
    const size_t uiH = rIceData.vSliceBias[ 1 ].size( );
    for( size_t uiI = uiFrom; uiI < uiTo; uiI++ )
    {
        assert( uiI < vvNormalized.size( ) );

        vvNormalized[ uiI ][ bA ? 0 : 1 ] = rIceData.vSliceBias[ 0 ][ uiI / uiH ] //
                                            * rIceData.vSliceBias[ 1 ][ uiI % uiH ] //
                                            * iceGetCount( rIceData, uiI / uiH, uiI % uiH, bA );
    }
}

double PartialQuarry::iceMaxBias( IceData& rIceData, bool bCol )
{
    double bMaxB = 0.0;
    for( size_t uiI = 0; uiI < rIceData.vSliceBias[ bCol ? 0 : 1 ].size( ); uiI++ )
        bMaxB = std::max( bMaxB, rIceData.vSliceBias[ bCol ? 0 : 1 ][ uiI ] );
    return bMaxB;
}

bool PartialQuarry::normalizeCoolIC( )
{
    if( vAxisCords[ 0 ].size( ) != vAxisCords[ 1 ].size( ) )
        doNotNormalize( );
    else
    {
        vvNormalized.resize( vvFlatValues.size( ) );
        for( size_t uiI = 0; uiI < 2; uiI++ )
        {
            std::vector<size_t> vCnt;
            for( auto& rFlat : vvFlatValues )
                vCnt.push_back( rFlat[ uiI ] );
            auto vRet = normalizeCoolerTrampoline( vCnt, vAxisCords[ 0 ].size( ) );
            for( size_t uiX = 0; uiX < vRet.size( ); uiX++ )
                vvNormalized[ uiX ][ uiI ] = vRet[ uiX ];
        }
    }
    CANCEL_RETURN;
    END_RETURN;
}

bool PartialQuarry::normalizeIC( )
{
    vvNormalized.resize( vvFlatValues.size( ) );

    size_t uiW = vAxisCords[ 0 ].size( );
    size_t uiH = vAxisCords[ 1 ].size( );

    const size_t uiMaxIters = 200;
    const double fTol = 1e-5;
    for( size_t uiI = 0; uiI < 2; uiI++ )
    {
        CANCEL_RETURN;
        IceData xData = { .vSliceBias =
                              std::array<std::vector<double>, 2>{
                                  std::vector<double>( uiW, 1.0 ),
                                  std::vector<double>( uiH, 1.0 ),
                              },
                          .vSliceMargin =
                              std::array<std::vector<double>, 2>{
                                  std::vector<double>( uiW, 0.0 ),
                                  std::vector<double>( uiH, 0.0 ),
                              },
                          .vBiases = std::vector<double>( uiW * uiH, 1.0 ) };
        std::array<double, 2> vVar{ 0, 0 };
        std::array<double, 2> vMean{ 0, 0 };
        for( size_t uiItr = 0; uiItr < uiMaxIters; uiItr++ )
        {
            CANCEL_RETURN;
            iceFilter( xData, 0, xData.vBiases.size( ) );
            iceTimesOuterProduct( xData, uiI == 0, 0, xData.vBiases.size( ) );
            for( bool bCol : { true, false } )
            {
                CANCEL_RETURN;
                iceMarginalize( xData, bCol, 0, xData.vSliceBias[ bCol ? 0 : 1 ].size( ) );
                double fMean = iceNonZeroMarginMean( xData, bCol );
                iceDivByMargin( xData, bCol, fMean, 0, xData.vSliceBias[ bCol ? 0 : 1 ].size( ) );
                vVar[ bCol ? 0 : 1 ] = iceNonZeroMarginVariance( xData, bCol, fMean );
                vMean[ bCol ? 0 : 1 ] = fMean;
            }

            if( vVar[ 0 ] < fTol && vVar[ 1 ] < fTol )
                break;
        }
        CANCEL_RETURN;
        size_t uiMaxFlat = 0;
        for( size_t uiJ = 0; uiJ < vvFlatValues.size( ); uiJ++ )
        {
            CANCEL_RETURN;
            uiMaxFlat = std::max( uiMaxFlat, vvFlatValues[ uiJ ][ uiI ] );
        }
        if( uiMaxFlat > 0 )
        {
            if( vVar[ 0 ] >= fTol || vVar[ 1 ] >= fTol )
            {
                setError( "iterative correction did not converge (var=" + std::to_string( vVar[ 0 ] ) + ", " +
                          std::to_string( vVar[ 1 ] ) + " mean=" + std::to_string( vMean[ 0 ] ) + ", " +
                          std::to_string( vMean[ 1 ] ) + "), showing data anyways" );
                iceApplyBias( xData, uiI == 0, 0, vvNormalized.size( ) );
            }
            else if( iceMaxBias( xData, true ) == 0 || iceMaxBias( xData, false ) == 0 )
            {
                setError( "iterative correction converged to zero, showing un-normalized data" );
                for( size_t uiJ = 0; uiJ < vvFlatValues.size( ); uiJ++ )
                {
                    CANCEL_RETURN;
                    vvNormalized[ uiJ ][ uiI ] = vvFlatValues[ uiJ ][ uiI ];
                }
            }
            else
                iceApplyBias( xData, uiI == 0, 0, vvNormalized.size( ) );
        }
    }
    END_RETURN;
}

bool PartialQuarry::setNormalized( )
{
    vvNormalized.clear( );
    vvNormalized.reserve( vvFlatValues.size( ) );

    if( getValue<std::string>( { "settings", "normalization", "normalize_by" } ) == "dont" )
        return doNotNormalize( );
    else if( getValue<std::string>( { "settings", "normalization", "normalize_by" } ) == "radicl-seq" )
        return normalizeBinominalTest( );
    else if( getValue<std::string>( { "settings", "normalization", "normalize_by" } ) == "rpm" )
        return normalizeSize( 1000000 );
    else if( getValue<std::string>( { "settings", "normalization", "normalize_by" } ) == "rpk" )
        return normalizeSize( 1000 );
    else if( getValue<std::string>( { "settings", "normalization", "normalize_by" } ) == "hi-c" )
        return normalizeIC( );
    else if( getValue<std::string>( { "settings", "normalization", "normalize_by" } ) == "cool-hi-c" )
        return normalizeCoolIC( );
    else
        throw std::logic_error( "invalid value for normalize_by" );
}

bool PartialQuarry::setDistDepDecayRemoved( )
{
    if( getValue<bool>( { "settings", "normalization", "ddd" } ) )
        for( size_t uiI = 0; uiI < vvNormalized.size( ); uiI++ )
            for( size_t uiJ = 0; uiJ < 2; uiJ++ )
                if( vBinCoords[ uiI ][ uiJ ].uiDecayCoordIndex != std::numeric_limits<size_t>::max( ) )
                {
                    CANCEL_RETURN;
                    if( vvFlatDecay[ vBinCoords[ uiI ][ uiJ ].uiDecayCoordIndex ][ uiJ ] > 0 )
                        vvNormalized[ uiI ][ uiJ ] /= vvFlatDecay[ vBinCoords[ uiI ][ uiJ ].uiDecayCoordIndex ][ uiJ ];
                    else
                        vvNormalized[ uiI ][ uiJ ] = 0;
                }
    END_RETURN;
}

bool PartialQuarry::setDivided( )
{
    vDivided.clear( );
    vDivided.reserve( vCombined.size( ) );

    const bool bByCol = getValue<bool>( { "settings", "normalization", "divide_by_column_coverage" } ) &&
                        vvCombinedCoverageValues[ 0 ].size( ) != 0;
    const bool bByRow = getValue<bool>( { "settings", "normalization", "divide_by_row_coverage" } ) &&
                        vvCombinedCoverageValues[ 1 ].size( ) != 0;

    for( size_t uiI = 0; uiI < vCombined.size( ); uiI++ )
    {
        CANCEL_RETURN;

        double fVal = vCombined[ uiI ];

        if( bByCol )
        {
            if( vvCombinedCoverageValues[ 0 ][ uiI / vAxisCords[ 0 ].size( ) ] == 0 )
                fVal = std::numeric_limits<double>::quiet_NaN( );
            else
                fVal /= vvCombinedCoverageValues[ 0 ][ uiI / vAxisCords[ 0 ].size( ) ];
        }
        if( bByRow )
        {
            if( vvCombinedCoverageValues[ 1 ][ uiI % vAxisCords[ 0 ].size( ) ] == 0 || std::isnan( fVal ) )
                fVal = std::numeric_limits<double>::quiet_NaN( );
            else
                fVal /= vvCombinedCoverageValues[ 1 ][ uiI % vAxisCords[ 0 ].size( ) ];
        }

        vDivided.push_back( fVal );
    }

    END_RETURN;
}

void PartialQuarry::regNormalization( )
{
    registerNode( NodeNames::Normalized,
                  ComputeNode{ .sNodeName = "normalized_bins",
                               .fFunc = &PartialQuarry::setNormalized,
                               .vIncomingFunctions = { NodeNames::FlatValues },
                               .vIncomingSession = { { "settings", "normalization", "p_accept", "val" } },
                               .vSessionsIncomingInPrevious = { { "replicates", "by_name" },
                                                                { "settings", "normalization", "normalize_by" } } } );

    registerNode( NodeNames::DistDepDecayRemoved,
                  ComputeNode{ .sNodeName = "dist_dep_dec_normalized_bins",
                               .fFunc = &PartialQuarry::setDistDepDecayRemoved,
                               .vIncomingFunctions = { NodeNames::Normalized, NodeNames::FlatDecay },
                               .vIncomingSession = { { "settings", "normalization", "ddd" } },
                               .vSessionsIncomingInPrevious = {} } );

    registerNode( NodeNames::Divided,
                  ComputeNode{ .sNodeName = "divided_by_tracks",
                               .fFunc = &PartialQuarry::setDivided,
                               .vIncomingFunctions = { NodeNames::Combined },
                               .vIncomingSession = { { "settings", "normalization", "divide_by_column_coverage" },
                                                     { "settings", "normalization", "divide_by_row_coverage" } },
                               .vSessionsIncomingInPrevious = {} } );
}

} // namespace cm