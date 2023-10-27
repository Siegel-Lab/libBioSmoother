#include "cm/partial_quarry.h"
#include "partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

bool PartialQuarry::normalizeSize( size_t uiSize )
{
    for( size_t uiY = 0; uiY < 3; uiY++ )
        for( auto& vVal : vvNormalizedDDD[ uiY ] )
        {
            CANCEL_RETURN;
            vvNormalized[ uiY ].push_back( std::array<double, 2>{
                vvFlatTotal[ uiY ][ 0 ] == 0 ? 0 : (double)uiSize * vVal[ 0 ] / (double)vvFlatTotal[ uiY ][ 0 ],
                vvFlatTotal[ uiY ][ 1 ] == 0 ? 0 : (double)uiSize * vVal[ 1 ] / (double)vvFlatTotal[ uiY ][ 1 ] } );
        }
    END_RETURN;
}

bool PartialQuarry::doNotNormalize( )
{
    for( size_t uiY = 0; uiY < 3; uiY++ )
        for( auto& vVal : vvNormalizedDDD[ uiY ] )
        {
            CANCEL_RETURN;
            vvNormalized[ uiY ].push_back( std::array<double, 2>{ vVal[ 0 ], vVal[ 1 ] } );
        }
    END_RETURN;
}


bool PartialQuarry::normalizeBinominalTest( )
{
    const bool bAxisIsCol = getValue<bool>( { "settings", "normalization", "radicl_seq_axis_is_column" } );

    for( size_t uiY = 0; uiY < 3; uiY++ )
    {
        vRadiclSeqCoverage[ uiY ].clear( );
        vRadiclSeqNumNonEmptyBins[ uiY ].clear( );

        const auto& rYCoords = pickYCoords( uiY );
        const auto& rXCoords = pickXCoords( uiY );

        vRadiclSeqCoverage[ uiY ].resize( bAxisIsCol ? rXCoords.size( ) : rYCoords.size( ) );
        vRadiclSeqNumNonEmptyBins[ uiY ].resize( bAxisIsCol ? rXCoords.size( ) : rYCoords.size( ) );

        for( size_t uiI = 0; uiI < vvPloidyValues[ uiY ].size( ); uiI++ )
        {
            CANCEL_RETURN;
            for( size_t uiJ = 0; uiJ < 2; uiJ++ )
            {
                size_t uiIdx = bAxisIsCol ? uiI / rYCoords.size( ) : uiI % rYCoords.size( );

                vRadiclSeqCoverage[ uiY ][ uiIdx ][ uiJ ] += vvPloidyValues[ uiY ][ uiI ][ uiJ ];
                if( vvPloidyValues[ uiY ][ uiI ][ uiJ ] > 0 )
                    vRadiclSeqNumNonEmptyBins[ uiY ][ uiIdx ][ uiJ ] += 1;
            }
        }

        std::vector<std::array<double, 2>> vIntermediate = normalizeBinominalTestTrampoline(
            vvPloidyValues[ uiY ], vRadiclSeqCoverage[ uiY ], vRadiclSeqNumNonEmptyBins[ uiY ],
            getValue<double>( { "settings", "normalization", "p_accept", "val" } ), bAxisIsCol, rYCoords.size( ) );
        CANCEL_RETURN;

        size_t uiH2;
        if( uiY == 2 )
            uiH2 = vV4cCoords[ 1 ].size( );
        else
            uiH2 = vAxisCords[ 1 ].size( );
        vvNormalized[ uiY ].resize( vBinCoords[ uiY ].size( ) );
        for( size_t uiX = 0; uiX < vIntermediate.size( ); uiX++ )
            for( size_t uiZ = 0; uiZ < 2; uiZ++ )
            {
                CANCEL_RETURN;
                const size_t uiXAxisId = vBinCoordsSampled[ uiY ][ uiX ][ uiZ ].uiXAxisIdx;
                const size_t uiYAxisId = vBinCoordsSampled[ uiY ][ uiX ][ uiZ ].uiYAxisIdx;
                if( uiXAxisId != SAMPLE_COORD && uiYAxisId != SAMPLE_COORD )
                {
                    const size_t uiX2 = uiYAxisId + uiH2 * uiXAxisId;
                    assert( uiX2 < vvNormalized[ uiY ].size( ) );
                    vvNormalized[ uiY ][ uiX2 ][ uiZ ] = vIntermediate[ uiX ][ uiZ ];
                }
            }
    }
    CANCEL_RETURN;
    END_RETURN;
}

/*
original ICing implementation:
@ https://github.com/open2c/cooler/blob/3d284485df070255f4a904102178b186c14c1cee/cooler/balance.py
@ https://github.com/open2c/cooler/blob/3d284485df070255f4a904102178b186c14c1cee/cooler/tools.py

*/


double PartialQuarry::iceGetCount( IceData& rIceData, size_t uiX, size_t uiY, size_t uiY_, bool bA,
                                   size_t uiIgnoreDiags )
{
    assert( uiX < rIceData.vSliceBias[ 0 ].size( ) );
    assert( uiY < rIceData.vSliceBias[ 1 ].size( ) );
    size_t uiIdx = uiY + uiX * ( rIceData.vSliceBias[ 1 ].size( ) );
    if( vBinCoordsSampled[ uiY_ ][ uiIdx ][ bA ? 0 : 1 ].uiActualContigIdX ==
        vBinCoordsSampled[ uiY_ ][ uiIdx ][ bA ? 0 : 1 ].uiActualContigIdY )
    {
        size_t uiXPos = vBinCoordsSampled[ uiY_ ][ uiIdx ][ bA ? 0 : 1 ].uiIndexX / uiBinWidth.r( );
        size_t uiYPos = vBinCoordsSampled[ uiY_ ][ uiIdx ][ bA ? 0 : 1 ].uiIndexY / uiBinHeight.r( );
        if( ( uiXPos > uiYPos ? uiXPos - uiYPos : uiYPos - uiXPos ) < uiIgnoreDiags )
            return 0;
    }
    assert( uiIdx < vvNormalizedDDD[ uiY_ ].size( ) );
    return vvNormalizedDDD[ uiY_ ][ uiIdx ][ bA ? 0 : 1 ];
}

void PartialQuarry::iceFilter( IceData& /*rIceData*/, size_t /*uiFrom*/, size_t /*uiTo*/ )
{}

void median_ready_sort( std::vector<double>& vX )
{
    // place the middle element at the center
    auto xItMiddle = vX.begin( ) + vX.size( ) / 2;
    std::nth_element( vX.begin( ), xItMiddle, vX.end( ) );
    if( vX.size( ) > 0 && vX.size( ) % 2 == 0 )
        // if the number of elements is even, the median is the average of the two middle elements
        // hence we place the element before the middle at the correct position as well
        // for this we use max_element, since its runtime complexity is O(n) instead of the O(n log n) of nth_element
        // swap then places the found element at the correct spot.
        std::swap( *std::max_element( vX.begin( ), xItMiddle ), *( xItMiddle - 1 ) );
}
double median( const std::vector<double>& vX )
{
    if( vX.size( ) == 0 )
        return 0;
    else if( vX.size( ) % 2 == 0 )
        return ( vX[ vX.size( ) / 2 ] + vX[ vX.size( ) / 2 - 1 ] ) / 2;
    else
        return vX[ vX.size( ) / 2 ];
}

void PartialQuarry::icePreFilter( IceData& rIceData, bool bCol, size_t uiFrom, size_t uiTo, size_t uiY, bool bA,
                                  size_t uiIgnoreDiags )
{
    const size_t uiMinNZ = getValue<size_t>( { "settings", "normalization", "ice_min_nz", "val" } );
    if( uiMinNZ > 0 )
    {
        const size_t uiWSlice = rIceData.vSliceBias[ bCol ? 1 : 0 ].size( );
        for( size_t uiI = uiFrom; uiI < uiTo; uiI++ )
        {
            size_t uiCnt = 0;
            for( size_t uiJ = 0; uiJ < uiWSlice; uiJ++ )
                if( iceGetCount( rIceData, bCol ? uiI : uiJ, bCol ? uiJ : uiI, uiY, bA, uiIgnoreDiags ) > 0 )
                    ++uiCnt;
            if( uiCnt < uiMinNZ )
                rIceData.vSliceBias[ bCol ? 0 : 1 ][ uiI ] = 0;
        }
    }


    // filter out rows and columns that have less than 1/4 of their cells filled
    double fFilter = getValue<double>( { "settings", "normalization", "ice_sparse_slice_filter", "val" } ) / 100.0;
    if( fFilter > 0 )
    {
        const size_t uiWSlice = rIceData.vSliceBias[ bCol ? 1 : 0 ].size( );
        for( size_t uiI = uiFrom; uiI < uiTo; uiI++ )
        {
            size_t uiCnt = 0;
            for( size_t uiJ = 0; uiJ < uiWSlice; uiJ++ )
                if( iceGetCount( rIceData, bCol ? uiI : uiJ, bCol ? uiJ : uiI, uiY, bA, uiIgnoreDiags ) > 0 )
                    ++uiCnt;
            if( (double)uiCnt <= uiWSlice * fFilter )
                rIceData.vSliceBias[ bCol ? 0 : 1 ][ uiI ] = 0;
        }
    }

    const size_t uiMaxMax = getValue<size_t>( { "settings", "normalization", "ice_mad_max", "val" } );
    if( uiMaxMax > 0 )
    {
        iceInit( rIceData, bA, 0, rIceData.vBiases.size( ), uiY, uiIgnoreDiags );
        iceFilter( rIceData, 0, rIceData.vBiases.size( ) );
        for( bool bCol : { true, false } )
        {
            iceMarginalize( rIceData, bCol, 0, rIceData.vSliceBias[ bCol ? 0 : 1 ].size( ) );

            std::vector<double> vMargins( rIceData.vSliceMargin[ bCol ? 0 : 1 ] );

            std::vector<double> vLogMargins( vMargins );
            vLogMargins.erase(
                std::remove_if( vLogMargins.begin( ), vLogMargins.end( ), []( double fVal ) { return fVal == 0; } ),
                vLogMargins.end( ) );
            median_ready_sort( vLogMargins );

            double fMedian = median( vLogMargins );
            for( double& fVal : vMargins )
                fVal /= fMedian;

            for( double& fVal : vLogMargins )
                fVal = std::log( fVal / fMedian );

            double fLogMedian = median( vLogMargins );
            for( double& fVal : vLogMargins )
                fVal = std::abs( fVal - fLogMedian );
            median_ready_sort( vLogMargins );
            double fDevLogMedian = median( vLogMargins );

            double fCutOff = std::exp( fLogMedian - fDevLogMedian * ( (double)uiMaxMax ) );

            for( size_t uiI = 0; uiI < rIceData.vSliceBias[ bCol ? 0 : 1 ].size( ); uiI++ )
                if( vMargins[ uiI ] < fCutOff )
                    rIceData.vSliceBias[ bCol ? 0 : 1 ][ uiI ] = 0;
        }
    }
}


void PartialQuarry::iceInit( IceData& rIceData, bool bA, size_t uiFrom, size_t uiTo, size_t uiY, size_t uiIgnoreDiags )
{
    const size_t uiH = rIceData.vSliceBias[ 1 ].size( );
    for( size_t uiI = uiFrom; uiI < uiTo; uiI++ )
    {
        assert( uiI / uiH < rIceData.vSliceBias[ 0 ].size( ) );
        assert( uiI % uiH < rIceData.vSliceBias[ 1 ].size( ) );
        assert( uiI < rIceData.vBiases.size( ) );
        rIceData.vBiases[ uiI ] = iceGetCount( rIceData, uiI / uiH, uiI % uiH, uiY, bA, uiIgnoreDiags );
    }
}

void PartialQuarry::iceTimesOuterProduct( IceData& rIceData, size_t uiFrom, size_t uiTo )
{
    const size_t uiH = rIceData.vSliceBias[ 1 ].size( );
    for( size_t uiI = uiFrom; uiI < uiTo; uiI++ )
    {
        assert( uiI / uiH < rIceData.vSliceBias[ 0 ].size( ) );
        assert( uiI % uiH < rIceData.vSliceBias[ 1 ].size( ) );
        assert( uiI < rIceData.vBiases.size( ) );
        rIceData.vBiases[ uiI ] *= rIceData.vSliceBias[ 0 ][ uiI / uiH ] * rIceData.vSliceBias[ 1 ][ uiI % uiH ];
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

            rIceData.vSliceMargin[ bCol ? 0 : 1 ][ uiI ] += 2 * rIceData.vBiases[ uiK ];
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

void PartialQuarry::rescaleBias( IceData& rIceData, bool bCol, double fMean, size_t uiFrom, size_t uiTo )
{
    double fMeanSqrt = std::sqrt( fMean );
    for( size_t uiI = uiFrom; uiI < uiTo; uiI++ )
    {
        assert( uiI < rIceData.vSliceMargin[ bCol ? 0 : 1 ].size( ) );

        rIceData.vSliceBias[ bCol ? 0 : 1 ][ uiI ] /= fMeanSqrt;
    }
}


bool PartialQuarry::hasChr( size_t uiI, const std::string& sName )
{
    for( size_t uiJ = 0; uiJ < this->vActiveChromosomes[ uiI ].size( ); uiJ++ )
        if( this->vActiveChromosomes[ uiI ][ uiJ ].sName == sName )
            return true;
    return false;
}

size_t PartialQuarry::getChromIdxForAnnoIdx( size_t uiAnnoIdx )
{
    const auto& rIt = std::lower_bound( vChromIdForAnnoIdx.begin( ),
                                        vChromIdForAnnoIdx.end( ),
                                        std::make_pair( uiAnnoIdx, this->vActiveChromosomes[ 0 ].size( ) + 1 ) );
    assert( rIt != vChromIdForAnnoIdx.begin( ) );
    return ( rIt - 1 )->second;
}

// @todo-low-prio ICe samples could also be computed using this...
bool PartialQuarry::getSamples( const GetSamplesMode& rMode, const size_t uiNumSamples, const std::string& sAnno,
                                const size_t uiBinSize, const bool bAxisIsCol, const bool bOnlyChromsOnBothAxes,
                                std::vector<AnnoCoord>& rOut )
{
    const uint32_t uiDividend = getValue<uint32_t>( { "dividend" } );

    const size_t uiColumn = bAxisIsCol ? 0 : 1;

    std::vector<std::pair<size_t, size_t>> vChromIdForSampleIdx;
    vChromIdForSampleIdx.reserve( this->vActiveChromosomes[ uiColumn ].size( ) );


    size_t uiFistAnnoIdx;
    if( rMode != GetSamplesMode::Bins )
        uiFistAnnoIdx = getValue<size_t>( { "annotation", "by_name", sAnno } );

    size_t uiMaxSamples = 0;
    for( size_t uiI = 0; uiI < this->vActiveChromosomes[ uiColumn ].size( ); uiI++ )
    {
        const ChromDesc& rActiveChrom = this->vActiveChromosomes[ uiColumn ][ uiI ];
        if( ( !bOnlyChromsOnBothAxes || hasChr( 1 - uiColumn, rActiveChrom.sName ) ) )
        {
            CANCEL_RETURN;
            switch( rMode )
            {
                case GetSamplesMode::OneAnnotation:
                    vChromIdForSampleIdx.emplace_back( uiMaxSamples, uiI );
                    uiMaxSamples += pIndices->vAnno.numIntervals( uiFistAnnoIdx + rActiveChrom.uiActualContigId );
                    break;
                case GetSamplesMode::Bins:
                    vChromIdForSampleIdx.emplace_back( uiMaxSamples, uiI );
                    uiMaxSamples += rActiveChrom.uiLength / uiBinSize;
                    break;
                    // case GetSamplesMode::BinnedAnno:
                    //     vChromIdForSampleIdx.emplace_back( uiMaxSamples, uiI );
                    //     uiMaxSamples +=
                    //         pIndices->vAnno.totalIntervalSize( rAnnoJson[ rActiveChrom.sName ].get<int64_t>( ) ) /
                    //         uiBinSize;
                    //     break;

                default:
                    assert( false );
                    break;
            }
        }
    }

    rOut.clear( );
    rOut.reserve( std::min( uiNumSamples, uiMaxSamples ) );
    std::set<size_t> vSeenDescIds;
    iterateEvenlyDividableByMaxPowerOfTwo( 0, uiMaxSamples, [ & ]( size_t uiVal ) {
        // the CANCEL_RETURN works since this also returns false
        // and the lambda function breaks the loop by returning false
        CANCEL_RETURN;

        const auto& rIt = std::lower_bound( vChromIdForSampleIdx.begin( ),
                                            vChromIdForSampleIdx.end( ),
                                            std::make_pair( uiVal, this->vActiveChromosomes[ 0 ].size( ) + 1 ) );
        assert( rIt != vChromIdForSampleIdx.begin( ) );
        const size_t uiChrom = ( rIt - 1 )->second;


        assert( uiChrom < this->vActiveChromosomes[ uiColumn ].size( ) );
        const size_t uiChr = this->vActiveChromosomes[ uiColumn ][ uiChrom ].uiActualContigId;

        switch( rMode )
        {
            case GetSamplesMode::OneAnnotation: {
                const auto rIntervalIt = pIndices->vAnno.get( uiFistAnnoIdx + uiChr, uiVal );
                const size_t uiDescId = rIntervalIt->uiDescId;
                if( vSeenDescIds.count( uiDescId ) == 0 )
                {
                    const size_t uiAnnoStart = rIntervalIt->uiAnnoStart / uiDividend;
                    const size_t uiAnnoEnd = std::max( uiAnnoStart + 1, rIntervalIt->uiAnnoEnd / uiDividend );
                    vSeenDescIds.insert( uiDescId );
                    rOut.push_back( AnnoCoord{
                        /*{*/
                        /*.uiChromosome =*/uiChrom, /*.uiIndexPos =*/uiAnnoStart,
                        /*.uiIndexSize =*/uiAnnoEnd - uiAnnoStart,
                        /* .uiActualContigId =*/this->vActiveChromosomes[ uiColumn ][ uiChrom ].uiActualContigId,
                        /* .uiChromId =*/this->vActiveChromosomes[ uiColumn ][ uiChrom ].uiCorrectedContigId,
                        /*}*/
                        /*.uiAnnoId =*/uiVal } );
                }
                break;
            }
            case GetSamplesMode::Bins:
                rOut.push_back(
                    AnnoCoord{ /*{*/
                               /*.uiChromosome =*/uiChrom,
                               /*.uiIndexPos =*/( uiVal - ( rIt - 1 )->first ) * uiBinSize,
                               /*.uiIndexSize =*/uiBinSize,
                               /* .uiActualContigId =*/this->vActiveChromosomes[ uiColumn ][ uiChrom ].uiActualContigId,
                               /* .uiChromId =*/this->vActiveChromosomes[ uiColumn ][ uiChrom ].uiCorrectedContigId,
                               /*}*/
                               /*.uiAnnoId =*/0 } );
                break;

                // case GetSamplesMode::BinnedAnno:
                //     break;

            default:
                assert( false );
                break;
        }


        return rOut.size( ) < uiNumSamples;
    } );
    // then this cancel_return is necessary to catch the cancelled iterateEvenlyDividableByMaxPowerOfTwo loop
    CANCEL_RETURN;
    END_RETURN;
}


bool PartialQuarry::setGridSeqSamples( )
{
    const std::string sNorm = getValue<std::string>( { "settings", "normalization", "normalize_by" } );

    vGridSeqSamples.clear( );

    if( sNorm != "grid-seq" )
        END_RETURN;


    const std::string sAnno = getValue<std::string>( { "settings", "normalization", "grid_seq_annotation" } );
    const bool bAxisIsCol = getValue<bool>( { "settings", "normalization", "grid_seq_axis_is_column" } );
    const bool bAllSamples = getValue<bool>( { "settings", "normalization", "grid_seq_global" } );
    size_t uiGridSeqSamples = getValue<size_t>( { "settings", "normalization", "grid_seq_samples", "val" } );
    if( bAllSamples )
        uiGridSeqSamples = std::numeric_limits<size_t>::max( );

    getSamples( GetSamplesMode::OneAnnotation, uiGridSeqSamples, sAnno, 0, bAxisIsCol, true, vGridSeqSamples );
    // then this cancel_return is necessary to catch the cancelled getSamples
    CANCEL_RETURN;


    END_RETURN;
}

bool PartialQuarry::setGridSeqCoverage( )
{
    const std::string sNorm = getValue<std::string>( { "settings", "normalization", "normalize_by" } );

    vGridSeqAnnoCoverage.clear( );

    if( sNorm != "grid-seq" )
        END_RETURN;

    const size_t uiMinuend = getValue<size_t>( { "settings", "normalization", "min_interactions", "val" } );
    const bool bAxisIsCol = getValue<bool>( { "settings", "normalization", "grid_seq_axis_is_column" } );
    const bool bFilterIntersection = getValue<bool>( { "settings", "normalization", "grid_seq_filter_intersection" } );
    const size_t uiCoverageMaxBinSize =
        getValue<size_t>( { "settings", "normalization", "grid_seq_max_bin_size", "val" } );

    CANCEL_RETURN;

    vGridSeqAnnoCoverage.reserve( vGridSeqSamples.size( ) );

    // collect & combine replicate data
    for( const IndexCoord& rSample : vGridSeqSamples )
    {
        if( bFilterIntersection )
            vGridSeqAnnoCoverage.push_back(
                { std::numeric_limits<size_t>::max( ), std::numeric_limits<size_t>::max( ) } );
        else
            vGridSeqAnnoCoverage.push_back( { 0, 0 } );

        for( size_t uiI = 0; uiI < vActiveReplicates.size( ); uiI++ )
        {
            CANCEL_RETURN;

            const size_t uiAnnoStart = rSample.uiIndexPos;
            const size_t uiAnnoEnd = rSample.uiIndexPos + rSample.uiIndexSize;

            for( size_t uiK = 0; uiK < 2; uiK++ )
            {
                std::array<size_t, 2> vVals;
                for( size_t uiJ = 0; uiJ < 2; uiJ++ )
                {
                    if( uiK == 0 )
                        vVals[ uiJ ] = getCoverageFromRepl( rSample.uiActualContigId, uiAnnoStart, uiAnnoEnd, uiI,
                                                            bAxisIsCol, uiJ != 0 );
                    else
                        vVals[ uiJ ] = getMaxCoverageFromRepl( rSample.uiActualContigId, uiAnnoStart, uiAnnoEnd, uiI,
                                                               uiCoverageMaxBinSize, !bAxisIsCol, uiJ != 0 );
                    vVals[ uiJ ] = vVals[ uiJ ] > uiMinuend ? vVals[ uiJ ] - uiMinuend : 0;
                }

                const size_t uiVal = symmetry( vVals[ 0 ], vVals[ 1 ] );

                if( bFilterIntersection )
                    vGridSeqAnnoCoverage.back( )[ uiK ] = std::min( vGridSeqAnnoCoverage.back( )[ uiK ], uiVal );
                else
                    vGridSeqAnnoCoverage.back( )[ uiK ] = std::max( vGridSeqAnnoCoverage.back( )[ uiK ], uiVal );
            }
        }
    }

    END_RETURN;
}

bool PartialQuarry::setRnaAssociatedGenesFilter( )
{

    vGridSeqFiltered.clear( );

    if( getValue<std::string>( { "settings", "normalization", "normalize_by" } ) != "grid-seq" )
        END_RETURN;

    const size_t uiRnaMin = getValue<size_t>( { "settings", "normalization", "grid_seq_rna_filter", "val_min" } );
    const size_t uiRnaMax = getValue<size_t>( { "settings", "normalization", "grid_seq_rna_filter", "val_max" } );
    const size_t uiDnaMin = getValue<size_t>( { "settings", "normalization", "grid_seq_dna_filter", "val_min" } );
    const size_t uiDnaMax = getValue<size_t>( { "settings", "normalization", "grid_seq_dna_filter", "val_max" } );

    vGridSeqFiltered.reserve( vGridSeqAnnoCoverage.size( ) );

    for( const auto& vArr : vGridSeqAnnoCoverage )
    {
        CANCEL_RETURN;
        vGridSeqFiltered.push_back(
            { vArr[ 0 ] >= uiRnaMin && vArr[ 0 ] < uiRnaMax, vArr[ 1 ] >= uiDnaMin && vArr[ 1 ] < uiDnaMax } );
    }

    END_RETURN;
}

// @todo-low-prio think about blacking out columns that are not chromatin associated
bool PartialQuarry::setRnaAssociatedBackground( )
{

    vBackgroundGridSeq.clear( );

    if( getValue<std::string>( { "settings", "normalization", "normalize_by" } ) != "grid-seq" )
        END_RETURN;

    const bool bAxisIsCol = getValue<bool>( { "settings", "normalization", "grid_seq_axis_is_column" } );
    const bool bIgnoreCis = getValue<bool>( { "settings", "normalization", "grid_seq_ignore_cis" } );
    bool bAtLeastOneTransContig = false;

    vBackgroundGridSeq.reserve( vAxisCords[ bAxisIsCol ? 0 : 1 ].size( ) );
    // for each eligible element
    for( const AxisCoord& rAxis : vAxisCords[ bAxisIsCol ? 0 : 1 ] )
    {
        vBackgroundGridSeq.push_back( 0 );
        for( size_t uiI = 0; uiI < vGridSeqAnnoCoverage.size( ); uiI++ )
            if( vGridSeqFiltered[ uiI ][ 0 ] && vGridSeqFiltered[ uiI ][ 1 ] ) // element is not filtered out
            {
                const IndexCoord& rSample = vGridSeqSamples[ uiI ];
                if( rSample.uiChromosome != rAxis.uiChromosome || !bIgnoreCis ) // ignore cis-contigs
                {
                    bAtLeastOneTransContig = true;
                    const size_t uiAnnoStart = rSample.uiIndexPos;
                    const size_t uiAnnoEnd = rSample.uiIndexPos + rSample.uiIndexSize;
                    for( size_t uiRepl = 0; uiRepl < vActiveReplicates.size( ); uiRepl++ )
                    {
                        size_t uiStartX;
                        size_t uiStartY;
                        size_t uiEndX;
                        size_t uiEndY;
                        size_t uiChrX;
                        size_t uiChrY;
                        if( bAxisIsCol )
                        {
                            uiStartX = rAxis.uiIndexPos;
                            uiStartY = uiAnnoStart;
                            uiEndX = rAxis.uiIndexPos + rAxis.uiIndexSize;
                            uiEndY = uiAnnoEnd;
                            uiChrX = rAxis.uiChromosome;
                            uiChrY = rSample.uiChromosome;
                        }
                        else
                        {
                            uiStartX = uiAnnoStart;
                            uiStartY = rAxis.uiIndexPos;
                            uiEndX = uiAnnoEnd;
                            uiEndY = rAxis.uiIndexPos + rAxis.uiIndexSize;
                            uiChrX = rSample.uiChromosome;
                            uiChrY = rAxis.uiChromosome;
                        }
                        size_t iDataSetId = getDatasetIdfromReplAndChr(
                            uiRepl,
                            vActiveChromosomes[ bAxisIsCol ? 0 : 1 ][ uiChrX ].uiActualContigId,
                            vActiveChromosomes[ bAxisIsCol ? 1 : 0 ][ uiChrY ].uiActualContigId );

                        vBackgroundGridSeq.back( ) +=
                            indexCount( uiRepl, iDataSetId, uiStartY, uiStartX, uiEndY, uiEndX, true );
                    }
                }
            }
    }

    if( !bAtLeastOneTransContig && bIgnoreCis )
        setError( "You are only cis interactions for the assoc. slices normalization; however there was not a single "
                  "trans interaction in the index. If your genome has only a single contig, you will have to include "
                  "cis interactions in the normalization." );

    END_RETURN;
}

bool PartialQuarry::normalizeGridSeq( )
{
    const bool bAxisIsCol = getValue<bool>( { "settings", "normalization", "grid_seq_axis_is_column" } );

    for( size_t uiY = 0; uiY < 3; uiY++ )
    {
        const auto& rYCoords = uiY == 2 ? vV4cCoords[ 1 ] : vAxisCords[ 1 ];
        for( size_t uiI = 0; uiI < vvNormalizedDDD[ uiY ].size( ); uiI++ )
        {
            CANCEL_RETURN;
            std::array<double, 2> vdVal = { vvNormalizedDDD[ uiY ][ uiI ][ 0 ], vvNormalizedDDD[ uiY ][ uiI ][ 1 ] };

            for( size_t uiK = 0; uiK < 2; uiK++ )
                vdVal[ uiK ] /= (double)std::max(
                    (size_t)1, vBackgroundGridSeq[ bAxisIsCol ? uiI / rYCoords.size( ) : uiI % rYCoords.size( ) ] );

            vvNormalized[ uiY ].push_back( vdVal );
        }
    }
    END_RETURN;
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
        vvNormalized[ 0 ].resize( vvNormalizedDDD[ 0 ].size( ) );
        for( size_t uiI = 0; uiI < 2; uiI++ )
        {
            std::vector<size_t> vCnt;
            for( auto& rFlat : vvNormalizedDDD[ 0 ] )
                vCnt.push_back( rFlat[ uiI ] );
            auto vRet = normalizeCoolerTrampoline( vCnt, vAxisCords[ 0 ].size( ) );
            for( size_t uiX = 0; uiX < vRet.size( ); uiX++ )
                vvNormalized[ 0 ][ uiX ][ uiI ] = vRet[ uiX ];
        }

        for( size_t uiY = 1; uiY < 3; uiY++ )
        {
            vvNormalized[ uiY ].resize( vvNormalizedDDD[ uiY ].size( ) );
            for( auto& vVal : vvNormalizedDDD[ uiY ] )
            {
                CANCEL_RETURN;
                vvNormalized[ uiY ].push_back( std::array<double, 2>{ vVal[ 0 ], vVal[ 1 ] } );
            }
        }
    }
    CANCEL_RETURN;
    END_RETURN;
}


bool PartialQuarry::normalizeIC( )
{
    const size_t uiIgnoreDiags = getValue<size_t>( { "settings", "normalization", "ice_ignore_n_diags", "val" } );
    for( size_t uiY = 0; uiY < NUM_COORD_SYSTEMS; uiY++ )
    {
        const auto& rXCoords = pickXCoords( uiY );
        const auto& rYCoords = pickYCoords( uiY );

        size_t uiW = rXCoords.size( );
        size_t uiH = rYCoords.size( );


        const size_t uiMaxIters = 200;
        const double fTol = 1e-5;
        for( size_t uiI = 0; uiI < 2; uiI++ )
        {
            CANCEL_RETURN;
            IceData xData = { /*.vSliceBias =*/
                              std::array<std::vector<double>, 2>{
                                  std::vector<double>( uiW, 1.0 ),
                                  std::vector<double>( uiH, 1.0 ),
                              },
                              /*.vSliceMargin =*/
                              std::array<std::vector<double>, 2>{
                                  std::vector<double>( uiW, 0.0 ),
                                  std::vector<double>( uiH, 0.0 ),
                              },
                              /*.vBiases =*/std::vector<double>( uiW * uiH, 1.0 ) };
            std::array<double, 2> vVar{ 0, 0 };
            std::array<double, 2> vMean{ 0, 0 };
            for( bool bCol : { true, false } )
                icePreFilter( xData, bCol, 0, xData.vSliceBias[ bCol ? 0 : 1 ].size( ), uiY, uiI == 0, uiIgnoreDiags );
            for( size_t uiItr = 0; uiItr < uiMaxIters; uiItr++ )
            {
                CANCEL_RETURN;
                iceInit( xData, uiI == 0, 0, xData.vBiases.size( ), uiY, uiIgnoreDiags );
                iceFilter( xData, 0, xData.vBiases.size( ) );
                iceTimesOuterProduct( xData, 0, xData.vBiases.size( ) );
                for( bool bCol : { true, false } )
                {
                    CANCEL_RETURN;
                    iceMarginalize( xData, bCol, 0, xData.vSliceBias[ bCol ? 0 : 1 ].size( ) );
                    double fMean = iceNonZeroMarginMean( xData, bCol );
                    iceDivByMargin( xData, bCol, fMean, 0, xData.vSliceBias[ bCol ? 0 : 1 ].size( ) );
                    vVar[ bCol ? 0 : 1 ] = iceNonZeroMarginVariance( xData, bCol, fMean );
                    vMean[ bCol ? 0 : 1 ] = fMean;
                }

                if( ( vVar[ 0 ] + vVar[ 1 ] ) / 2.0 < fTol )
                    break;
            }
            for( bool bCol : { true, false } )
                rescaleBias( xData, bCol, vMean[ bCol ? 0 : 1 ], 0, xData.vSliceBias[ bCol ? 0 : 1 ].size( ) );
            CANCEL_RETURN;
            double fMaxFlat = 0;
            for( size_t uiJ = 0; uiJ < vvNormalizedDDD[ uiY ].size( ); uiJ++ )
            {
                CANCEL_RETURN;
                fMaxFlat = std::max( fMaxFlat, vvNormalizedDDD[ uiY ][ uiJ ][ uiI ] );
            }
            if( fMaxFlat > 0 )
            {
                std::string sErrorWhere =
                    std::array<std::string, 3>{ "heatmap", "V4C on rows", "V4C on columns" }[ uiY ];
                if( vVar[ 0 ] >= fTol || vVar[ 1 ] >= fTol )
                {
                    setError( "iterative correction did not converge for the " + sErrorWhere +
                              " (var=" + std::to_string( vVar[ 0 ] ) + ", " + std::to_string( vVar[ 1 ] ) +
                              " mean=" + std::to_string( vMean[ 0 ] ) + ", " + std::to_string( vMean[ 1 ] ) +
                              "), showing data anyways" );
                }
                else if( iceMaxBias( xData, true ) == 0 || iceMaxBias( xData, false ) == 0 )
                {
                    setError( "iterative correction converged to zero for the " + sErrorWhere +
                              ", showing un-normalized data" );
                    for( size_t uiY = 0; uiY < 2; uiY++ )
                        for( size_t uiJ = 0; uiJ < xData.vSliceBias[ uiY ].size( ); uiJ++ )
                        {
                            CANCEL_RETURN;
                            xData.vSliceBias[ uiY ][ uiJ ] = 1;
                        }
                }
            }
            else
            {
                for( size_t uiY = 0; uiY < 2; uiY++ )
                    for( size_t uiJ = 0; uiJ < xData.vSliceBias[ uiY ].size( ); uiJ++ )
                    {
                        CANCEL_RETURN;
                        xData.vSliceBias[ uiY ][ uiJ ] = 1;
                    }
            }

            vIceSliceBias[ uiY ][ 0 ][ uiI ].swap( xData.vSliceBias[ 0 ] );
            vIceSliceBias[ uiY ][ 1 ][ uiI ].swap( xData.vSliceBias[ 1 ] );
        }
        vvNormalized[ uiY ].resize( vBinCoords[ uiY ].size( ) );
        size_t uiH2;
        if( uiY == 2 )
            uiH2 = vV4cCoords[ 1 ].size( );
        else
            uiH2 = vAxisCords[ 1 ].size( );
        for( size_t uiX = 0; uiX < vvNormalizedDDD[ uiY ].size( ); uiX++ )
            for( size_t uiZ = 0; uiZ < 2; uiZ++ )
            {
                const size_t uiXAxisId = vBinCoordsSampled[ uiY ][ uiX ][ uiZ ].uiXAxisIdx;
                const size_t uiYAxisId = vBinCoordsSampled[ uiY ][ uiX ][ uiZ ].uiYAxisIdx;
                if( uiXAxisId != SAMPLE_COORD && uiYAxisId != SAMPLE_COORD )
                {
                    const size_t uiX2 = uiYAxisId + uiH2 * uiXAxisId;
                    assert( uiX2 < vvNormalized[ uiY ].size( ) );
                    vvNormalized[ uiY ][ uiX2 ][ uiZ ] = vIceSliceBias[ uiY ][ 0 ][ uiZ ][ uiX / uiH ] //
                                                         * vIceSliceBias[ uiY ][ 1 ][ uiZ ][ uiX % uiH ] //
                                                         * vvNormalizedDDD[ uiY ][ uiX ][ uiZ ];
                }
            }
    }

    END_RETURN;
}

bool PartialQuarry::setNormalized( )
{
    for( size_t uiY = 0; uiY < NUM_COORD_SYSTEMS; uiY++ )
    {
        vvNormalized[ uiY ].clear( );
        vvNormalized[ uiY ].reserve( vvNormalizedDDD[ uiY ].size( ) );
        for( size_t uiX = 0; uiX < 2; uiX++ )
            for( size_t uiI = 0; uiI < 2; uiI++ )
                vIceSliceBias[ uiY ][ uiX ][ uiI ].clear( );
    }

    const std::string sNorm = getValue<std::string>( { "settings", "normalization", "normalize_by" } );

    if( sNorm == "dont" )
        return doNotNormalize( );
    else if( sNorm == "radicl-seq" )
        return normalizeBinominalTest( );
    else if( sNorm == "rpm" )
        return normalizeSize( 1000000 );
    else if( sNorm == "rpk" )
        return normalizeSize( 1000 );
    else if( sNorm == "ice" )
        return normalizeIC( );
    else if( sNorm == "cool-ice" )
        return normalizeCoolIC( );
    else if( sNorm == "grid-seq" )
        return normalizeGridSeq( );
    else
        throw std::logic_error( "invalid value for normalize_by" );
}

bool PartialQuarry::setDistDepDecayRemoved( )
{
    if( getValue<bool>( { "settings", "normalization", "ddd" } ) )
    {
        for( size_t uiY = 0; uiY < NUM_COORD_SYSTEMS; uiY++ )
        {
            vvNormalizedDDD[ uiY ].resize( vvPloidyValues[ uiY ].size( ) );
            assert( vBinCoordsSampled[ uiY ].size( ) <= vvPloidyValues[ uiY ].size( ) );
            for( size_t uiI = 0; uiI < vBinCoordsSampled[ uiY ].size( ); uiI++ )
                for( size_t uiJ = 0; uiJ < 2; uiJ++ )
                    if( vBinCoordsSampled[ uiY ][ uiI ][ uiJ ].uiDecayCoordIndex !=
                        std::numeric_limits<size_t>::max( ) )
                    {
                        CANCEL_RETURN;
                        if( vvFlatDecay[ uiY ][ vBinCoordsSampled[ uiY ][ uiI ][ uiJ ].uiDecayCoordIndex ][ uiJ ] > 0 )
                            vvNormalizedDDD[ uiY ][ uiI ][ uiJ ] =
                                (double)vvPloidyValues[ uiY ][ uiI ][ uiJ ] /
                                (double)vvFlatDecay[ uiY ][ vBinCoordsSampled[ uiY ][ uiI ][ uiJ ].uiDecayCoordIndex ]
                                                   [ uiJ ];
                        else
                            vvNormalizedDDD[ uiY ][ uiI ][ uiJ ] = 0;
                    }
        }
    }
    else
        for( size_t uiY = 0; uiY < 3; uiY++ )
        {
            vvNormalizedDDD[ uiY ].resize( vvPloidyValues[ uiY ].size( ) );
            for( size_t uiI = 0; uiI < vBinCoordsSampled[ uiY ].size( ); uiI++ )
                for( size_t uiJ = 0; uiJ < 2; uiJ++ )
                    vvNormalizedDDD[ uiY ][ uiI ][ uiJ ] = (double)vvPloidyValues[ uiY ][ uiI ][ uiJ ];
        }
    END_RETURN;
}


const decltype( PartialQuarry::vScaled )
PartialQuarry::getScaled( const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::Scaled, fPyPrint );
    return vScaled;
}

void PartialQuarry::regNormalization( )
{
    registerNode( NodeNames::Normalized,
                  ComputeNode{ /*.sNodeName =*/"normalized_bins",
                               /*.sNodeDesc =*/"Normalize the heatmap bins.",
                               /*.fFunc =*/&PartialQuarry::setNormalized,
                               /*.vIncomingFunctions =*/
                               { NodeNames::DistDepDecayRemoved, NodeNames::RnaAssociatedBackground },
                               /*.vIncomingSession =*/
                               { { "settings", "normalization", "p_accept", "val" },
                                 { "settings", "normalization", "ice_symmetric_bias" },
                                 { "settings", "normalization", "ice_min_nz", "val" },
                                 { "settings", "normalization", "ice_ignore_n_diags", "val" },
                                 { "settings", "normalization", "ice_mad_max", "val" },
                                 { "settings", "normalization", "ice_sparse_slice_filter", "val" } },
                               /*.vSessionsIncomingInPrevious =*/
                               { { "contigs", "genome_size" },
                                 { "settings", "normalization", "normalize_by" },
                                 { "settings", "normalization", "grid_seq_axis_is_column" },
                                 { "settings", "normalization", "ice_local" },
                                 { "settings", "normalization", "radicl_seq_axis_is_column" } },
                               /*bHidden =*/false } );

    registerNode( NodeNames::GridSeqSamples,
                  ComputeNode{ /*.sNodeName =*/"grid_seq_samples",
                               /*.sNodeDesc =*/"Compute the sample bin coordinates for assoc. slices normalizaition.",
                               /*.fFunc =*/&PartialQuarry::setGridSeqSamples,
                               /*.vIncomingFunctions =*/{ NodeNames::AxisCoords },
                               /*.vIncomingSession =*/
                               { { "settings", "normalization", "normalize_by" },
                                 { "settings", "normalization", "grid_seq_annotation" },
                                 { "settings", "normalization", "grid_seq_samples", "val" },
                                 { "settings", "normalization", "grid_seq_axis_is_column" },
                                 { "settings", "normalization", "grid_seq_global" } },
                               /*.vSessionsIncomingInPrevious =*/{ { "annotation", "by_name" }, { "dividend" } },
                               /*bHidden =*/false } );

    registerNode( NodeNames::GridSeqCoverage,
                  ComputeNode{ /*.sNodeName =*/"grid_seq_coverage",
                               /*.sNodeDesc =*/"Compute the coverage of the assoc. slices normalization.",
                               /*.fFunc =*/&PartialQuarry::setGridSeqCoverage,
                               /*.vIncomingFunctions =*/
                               { NodeNames::ActiveReplicates, NodeNames::MappingQuality, NodeNames::Directionality,
                                 NodeNames::GridSeqSamples },
                               /*.vIncomingSession =*/
                               { { "settings", "normalization", "grid_seq_filter_intersection" },
                                 { "settings", "normalization", "grid_seq_max_bin_size", "val" },
                                 { "replicates", "by_name" },
                                 { "settings", "normalization", "min_interactions", "val" } },
                               /*.vSessionsIncomingInPrevious =*/
                               { { "settings", "normalization", "grid_seq_axis_is_column" },
                                 { "settings", "normalization", "normalize_by" },
                                 { "settings", "normalization", "grid_seq_annotation" },
                                 { "settings", "normalization", "grid_seq_samples", "val" },
                                 { "annotation", "by_name" },
                                 { "dividend" } },
                               /*bHidden =*/false } );

    registerNode( NodeNames::RnaAssociatedGenesFilter,
                  ComputeNode{ /*.sNodeName =*/"rna_associated_genes_filter",
                               /*.sNodeDesc =*/"Compute the associated annotations of the assoc. slices normalization.",
                               /*.fFunc =*/&PartialQuarry::setRnaAssociatedGenesFilter,
                               /*.vIncomingFunctions =*/{ NodeNames::GridSeqCoverage },
                               /*.vIncomingSession =*/
                               { { "settings", "normalization", "grid_seq_rna_filter", "val_min" },
                                 { "settings", "normalization", "grid_seq_rna_filter", "val_max" },
                                 { "settings", "normalization", "grid_seq_dna_filter", "val_min" },
                                 { "settings", "normalization", "grid_seq_dna_filter", "val_max" } },
                               /*.vSessionsIncomingInPrevious =*/{ { "settings", "normalization", "normalize_by" } },
                               /*bHidden =*/false } );

    registerNode(
        NodeNames::RnaAssociatedBackground,
        ComputeNode{ /*.sNodeName =*/"rna_associated_background",
                     /*.sNodeDesc =*/"Compute the RNA associated background of the assoc. slices normalization.",
                     /*.fFunc =*/&PartialQuarry::setRnaAssociatedBackground,
                     /*.vIncomingFunctions =*/{ NodeNames::RnaAssociatedGenesFilter },
                     /*.vIncomingSession =*/{ { "settings", "normalization", "grid_seq_ignore_cis" } },
                     /*.vSessionsIncomingInPrevious =*/
                     { { "settings", "normalization", "normalize_by" },
                       { "settings", "normalization", "grid_seq_axis_is_column" },
                       { "replicates", "by_name" },
                       { "dividend" } },
                     /*bHidden =*/false } );

    registerNode( NodeNames::DistDepDecayRemoved,
                  ComputeNode{ /*.sNodeName =*/"dist_dep_dec_normalized_bins",
                               /*.sNodeDesc =*/"Remove the distance dependent decay from the normalized bins.",
                               /*.fFunc =*/&PartialQuarry::setDistDepDecayRemoved,
                               /*.vIncomingFunctions =*/{ NodeNames::PloidyValues, NodeNames::FlatDecay },
                               /*.vIncomingSession =*/{ },
                               /*.vSessionsIncomingInPrevious =*/{ { "settings", "normalization", "ddd" } },
                               /*bHidden =*/false } );
}

} // namespace cm