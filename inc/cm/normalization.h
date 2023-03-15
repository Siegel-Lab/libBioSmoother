#include "cm/partial_quarry.h"
#include "partial_quarry.h"
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
    const bool bIsCol = getValue<bool>( { "settings", "normalization", "radicl_seq_axis_is_column" } );
    const size_t uiNumBinsInRowTotal =
        ( vCanvasSize[ bIsCol ? 1 : 0 ] - 1 ) / ( bIsCol ? uiBinHeight : uiBinWidth ) + 1;
    const size_t uiNumSaples = getValue<size_t>( { "settings", "normalization", "radicl_seq_samples", "val" } );
    vvNormalized = normalizeBinominalTestTrampoline(
        vvFlatValues, vRadiclSeqCoverage, vRadiclSeqNumNonEmptyBins, uiNumSaples, uiNumBinsInRowTotal,
        getValue<double>( { "settings", "normalization", "p_accept", "val" } ), bIsCol, vAxisCords[ 1 ].size( ) );
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
    assert( uiIdx < vvFlatValues.size( ) ); // @todo this assert triggered with max coverage per bin columns active
    return vvFlatValues[ uiIdx ][ bA ? 0 : 1 ];
}

void PartialQuarry::iceFilter( IceData& /*rIceData*/, size_t /*uiFrom*/, size_t /*uiTo*/ )
{}

void PartialQuarry::icePreFilter( IceData& rIceData, bool bCol, size_t uiFrom, size_t uiTo, bool bA )
{
    // filter out rows and columns that have less than 1/4 of their cells filled
    double fFilter = getValue<double>( { "settings", "normalization", "ice_sparse_slice_filter", "val" } ) / 100.0;
    if( fFilter > 0 )
    {
        const size_t uiWSlice = rIceData.vSliceBias[ bCol ? 1 : 0 ].size( );
        for( size_t uiI = uiFrom; uiI < uiTo; uiI++ )
        {
            size_t uiCnt = 0;
            for( size_t uiJ = 0; uiJ < uiWSlice; uiJ++ )
                if( iceGetCount( rIceData, bCol ? uiI : uiJ, bCol ? uiJ : uiI, bA ) > 0 )
                    ++uiCnt;
            if( (double)uiCnt <= uiWSlice * fFilter )
                rIceData.vSliceBias[ bCol ? 0 : 1 ][ uiI ] = 0;
        }
    }
}


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

bool PartialQuarry::getSamples( const GetSamplesMode& rMode, const size_t uiNumSamples, const std::string& sAnno,
                                const size_t uiBinSize, const bool bAxisIsCol, const bool bOnlyChromsOnBothAxes,
                                std::vector<AnnoCoord>& rOut )
{
    const uint32_t uiDividend = getValue<uint32_t>( { "dividend" } );

    const size_t uiColumn = bAxisIsCol ? 0 : 1;

    std::vector<std::pair<size_t, size_t>> vChromIdForSampleIdx;
    vChromIdForSampleIdx.reserve( this->vActiveChromosomes[ uiColumn ].size( ) );

    json rAnnoJson;
    if( rMode != GetSamplesMode::Bins )
        rAnnoJson = getValue<json>( { "annotation", "by_name", sAnno } );

    size_t uiMaxSamples = 0;
    for( size_t uiI = 0; uiI < this->vActiveChromosomes[ uiColumn ].size( ); uiI++ )
    {
        const ChromDesc& rActiveChrom = this->vActiveChromosomes[ uiColumn ][ uiI ];
        if( ( !bOnlyChromsOnBothAxes || hasChr( 1 - uiColumn, rActiveChrom.sName ) ) &&
            ( rMode == GetSamplesMode::Bins || rAnnoJson.contains( rActiveChrom.sName ) ) )
        {
            CANCEL_RETURN;
            switch( rMode )
            {
                case GetSamplesMode::OneAnnotation:
                    vChromIdForSampleIdx.emplace_back( uiMaxSamples, uiI );
                    uiMaxSamples += pIndices->vAnno.numIntervals( rAnnoJson[ rActiveChrom.sName ].get<int64_t>( ) );
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
    rOut.reserve( uiNumSamples );
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

        const std::string& rChr = this->vActiveChromosomes[ uiColumn ][ uiChrom ].sName;

        switch( rMode )
        {
            case GetSamplesMode::OneAnnotation: {
                const auto rIntervalIt = pIndices->vAnno.get( rAnnoJson[ rChr ].get<int64_t>( ), uiVal );
                const size_t uiDescId = rIntervalIt->uiDescId;
                if( vSeenDescIds.count( uiDescId ) == 0 )
                {
                    const size_t uiAnnoStart = rIntervalIt->uiAnnoStart / uiDividend;
                    const size_t uiAnnoEnd = std::max( uiAnnoStart + 1, rIntervalIt->uiAnnoEnd / uiDividend );
                    vSeenDescIds.insert( uiDescId );
                    rOut.push_back( AnnoCoord{ /*{*/ /*.uiChromosome =*/uiChrom, /*.uiIndexPos =*/uiAnnoStart,
                                               /*.uiIndexSize =*/uiAnnoEnd - uiAnnoStart /*}*/,
                                               /*.uiAnnoId =*/uiVal } );
                }
                break;
            }
            case GetSamplesMode::Bins:
                rOut.push_back( AnnoCoord{ /*{*/ /*.uiChromosome =*/uiChrom,
                                           /*.uiIndexPos =*/( uiVal - ( rIt - 1 )->first ) * uiBinSize,
                                           /*.uiIndexSize =*/uiBinSize /*}*/,
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
    if( sNorm != "grid-seq" )
        END_RETURN;


    const std::string sAnno = getValue<std::string>( { "settings", "normalization", "grid_seq_annotation" } );
    const bool bAxisIsCol = getValue<bool>( { "settings", "normalization", "grid_seq_axis_is_column" } );
    const size_t uiGridSeqSamples = getValue<size_t>( { "settings", "normalization", "grid_seq_samples", "val" } );

    getSamples( GetSamplesMode::OneAnnotation, uiGridSeqSamples, sAnno, 0, bAxisIsCol, true, vGridSeqSamples );
    // then this cancel_return is necessary to catch the cancelled getSamples
    CANCEL_RETURN;


    END_RETURN;
}

bool PartialQuarry::setRadiclSeqSamples( )
{
    const std::string sNorm = getValue<std::string>( { "settings", "normalization", "normalize_by" } );
    if( sNorm != "radicl-seq" )
        END_RETURN;


    const bool bAxisIsCol = getValue<bool>( { "settings", "normalization", "radicl_seq_axis_is_column" } );
    const size_t uiRadiclSeqSamples = getValue<size_t>( { "settings", "normalization", "radicl_seq_samples", "val" } );

    getSamples( GetSamplesMode::Bins, uiRadiclSeqSamples, "", bAxisIsCol ? uiBinHeight : uiBinWidth, bAxisIsCol, false,
                vRadiclSeqSamples );
    // then this cancel_return is necessary to catch the cancelled getSamples
    CANCEL_RETURN;


    END_RETURN;
}

bool PartialQuarry::setICESamples( )
{
    const std::string sNorm = getValue<std::string>( { "settings", "normalization", "normalize_by" } );
    if( sNorm != "hi-c" )
        END_RETURN;


    const size_t uiICESamples = getValue<size_t>( { "settings", "normalization", "hi_c_samples", "val" } );

    for( size_t uiI = 0; uiI < 2; uiI++ )
    {
        getSamples( GetSamplesMode::Bins, uiICESamples, "", uiI == 0 ? uiBinHeight : uiBinWidth, uiI == 0, false,
                    vICESamples[ uiI ] );
        // then this cancel_return is necessary to catch the cancelled getSamples
        CANCEL_RETURN;
    }

    END_RETURN;
}

bool PartialQuarry::setRadiclSeqCoverage( )
{
    const std::string sNorm = getValue<std::string>( { "settings", "normalization", "normalize_by" } );
    if( sNorm != "radicl-seq" )
        END_RETURN;

    const bool bAxisIsCol = getValue<bool>( { "settings", "normalization", "radicl_seq_axis_is_column" } );
    const size_t uiMinuend = getValue<size_t>( { "settings", "normalization", "min_interactions", "val" } );
    std::vector<size_t> vRawCoverage;
    vRawCoverage.reserve( vActiveReplicates.size( ) );

    vRadiclSeqCoverage.clear( );
    vRadiclSeqCoverage.reserve( vAxisCords[ bAxisIsCol ? 0 : 1 ].size( ) );
    vRadiclSeqNumNonEmptyBins.clear( );
    vRadiclSeqNumNonEmptyBins.reserve( vAxisCords[ bAxisIsCol ? 0 : 1 ].size( ) );

    for( const AxisCoord& rAxis : vAxisCords[ bAxisIsCol ? 0 : 1 ] )
    {
        CANCEL_RETURN;
        vRawCoverage.clear( );

        // collect raw values
        for( size_t uiI = 0; uiI < vActiveReplicates.size( ); uiI++ )
        {
            CANCEL_RETURN;

            std::array<size_t, 2> vVals;
            for( size_t uiJ = 0; uiJ < 2; uiJ++ )
            {
                CANCEL_RETURN;
                vVals[ uiJ ] = getCoverageFromRepl( rAxis, uiI, bAxisIsCol, uiJ != 0 );
                vVals[ uiJ ] = vVals[ uiJ ] > uiMinuend ? vVals[ uiJ ] - uiMinuend : 0;
            }

            vRawCoverage.push_back( symmetry( vVals[ 0 ], vVals[ 1 ] ) );
        }

        // combine raw values
        std::array<size_t, 2> vVals;
        for( size_t uiJ = 0; uiJ < 2; uiJ++ )
        {
            std::vector<size_t> vCollected;
            vCollected.reserve( vInGroup[ uiJ ].size( ) );
            for( size_t uiX : vInGroup[ uiJ ] )
            {
                CANCEL_RETURN;
                vCollected.push_back( vRawCoverage[ uiX ] );
            }

            vVals[ uiJ ] = getFlatValue( vCollected );
        }
        vRadiclSeqCoverage.push_back( vVals );

        vRadiclSeqNumNonEmptyBins.push_back( { 0, 0 } );
        for( const AnnoCoord& rSample : vRadiclSeqSamples )
        {
            CANCEL_RETURN;
            vRawCoverage.clear( );

            // collect raw values
            for( size_t uiRepl = 0; uiRepl < vActiveReplicates.size( ); uiRepl++ )
            {
                CANCEL_RETURN;

                std::array<size_t, 2> vVals;
                for( size_t uiJ = 0; uiJ < 2; uiJ++ )
                {
                    CANCEL_RETURN;
                    const bool bSymPart = uiJ != 0;

                    const size_t uiChrX = bAxisIsCol != bSymPart ? rAxis.uiChromosome : rSample.uiChromosome;
                    const size_t uiChrY = bAxisIsCol != bSymPart ? rSample.uiChromosome : rAxis.uiChromosome;

                    size_t uiDataSetId = getDatasetIdfromReplAndChr( uiRepl, uiChrX, uiChrY );
                    if( uiDataSetId != std::numeric_limits<size_t>::max( ) )
                    {
                        const size_t uiXMin = bAxisIsCol != bSymPart ? rAxis.uiIndexPos : rSample.uiIndexPos;
                        const size_t uiXMax = bAxisIsCol != bSymPart ? rAxis.uiIndexPos + rAxis.uiIndexSize
                                                                     : rSample.uiIndexPos + rSample.uiIndexSize;
                        const size_t uiYMin = bAxisIsCol != bSymPart ? rSample.uiIndexPos : rAxis.uiIndexPos;
                        const size_t uiYMax = bAxisIsCol != bSymPart ? rSample.uiIndexPos + rSample.uiIndexSize
                                                                     : rAxis.uiIndexPos + rAxis.uiIndexSize;

                        vVals[ uiJ ] = pIndices->count(
                            uiDataSetId,
                            { uiYMin, uiXMin, uiMapQMin, uiFromAnnoFilter, uiFromSameStrandFilter,
                              uiFromYStrandFilter },
                            { uiYMax, uiXMax, uiMapQMax, uiToAnnoFilter, uiToSameStrandFilter, uiToYStrandFilter },
                            xIntersect, 0 );
                        vVals[ uiJ ] = vVals[ uiJ ] > uiMinuend ? vVals[ uiJ ] - uiMinuend : 0;
                    }
                    else
                        vVals[ uiJ ] = 0;
                }

                vRawCoverage.push_back( symmetry( vVals[ 0 ], vVals[ 1 ] ) );
            }
            // @todo-low-prio code duplication...
            for( size_t uiJ = 0; uiJ < 2; uiJ++ )
            {
                std::vector<size_t> vCollected;
                vCollected.reserve( vInGroup[ uiJ ].size( ) );
                for( size_t uiX : vInGroup[ uiJ ] )
                {
                    CANCEL_RETURN;
                    vCollected.push_back( vRawCoverage[ uiX ] );
                }

                const size_t uiVal = getFlatValue( vCollected );
                if( uiVal > 0 )
                    ++vRadiclSeqNumNonEmptyBins.back( )[ uiJ ];
            }
        }
    }

    END_RETURN;
}

bool PartialQuarry::setGridSeqCoverage( )
{
    const std::string sNorm = getValue<std::string>( { "settings", "normalization", "normalize_by" } );
    if( sNorm != "grid-seq" )
        END_RETURN;

    const size_t uiMinuend = getValue<size_t>( { "settings", "normalization", "min_interactions", "val" } );
    const bool bAxisIsCol = getValue<bool>( { "settings", "normalization", "grid_seq_axis_is_column" } );
    const bool bFilterIntersection = getValue<bool>( { "settings", "normalization", "grid_seq_filter_intersection" } );
    const size_t uiCoverageMaxBinSize =
        getValue<size_t>( { "settings", "normalization", "grid_seq_max_bin_size", "val" } );

    CANCEL_RETURN;

    vGridSeqAnnoCoverage.clear( );
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
                        vVals[ uiJ ] = getCoverageFromRepl( rSample.uiChromosome, uiAnnoStart, uiAnnoEnd, uiI,
                                                            bAxisIsCol, uiJ != 0 );
                    else
                        vVals[ uiJ ] = getMaxCoverageFromRepl( rSample.uiChromosome, uiAnnoStart, uiAnnoEnd, uiI,
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
    if( getValue<std::string>( { "settings", "normalization", "normalize_by" } ) != "grid-seq" )
        END_RETURN;

    const size_t uiRnaMin = getValue<size_t>( { "settings", "normalization", "grid_seq_rna_filter", "val_min" } );
    const size_t uiRnaMax = getValue<size_t>( { "settings", "normalization", "grid_seq_rna_filter", "val_max" } );
    const size_t uiDnaMin = getValue<size_t>( { "settings", "normalization", "grid_seq_dna_filter", "val_min" } );
    const size_t uiDnaMax = getValue<size_t>( { "settings", "normalization", "grid_seq_dna_filter", "val_max" } );

    vGridSeqFiltered.clear( );
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
    if( getValue<std::string>( { "settings", "normalization", "normalize_by" } ) != "grid-seq" )
        END_RETURN;

    const bool bAxisIsCol = getValue<bool>( { "settings", "normalization", "grid_seq_axis_is_column" } );
    const std::string sAnno = getValue<std::string>( { "settings", "normalization", "grid_seq_annotation" } );
    const auto rJson = getValue<json>( { "annotation", "by_name", sAnno } );

    vBackgroundGridSeq.clear( );
    vBackgroundGridSeq.reserve( vAxisCords[ bAxisIsCol ? 0 : 1 ].size( ) );
    // for each eligible element
    for( const AxisCoord& rAxis : vAxisCords[ bAxisIsCol ? 0 : 1 ] )
    {
        vBackgroundGridSeq.push_back( 0 );
        for( size_t uiI = 0; uiI < vGridSeqAnnoCoverage.size( ); uiI++ )
            if( vGridSeqFiltered[ uiI ][ 0 ] && vGridSeqFiltered[ uiI ][ 1 ] ) // element is not filtered out
            {
                const IndexCoord& rSample = vGridSeqSamples[ uiI ];
                if( rSample.uiChromosome != rAxis.uiChromosome ) // ignore cis-contigs
                {
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
                        size_t iDataSetId = getDatasetIdfromReplAndChr( uiRepl, uiChrX, uiChrY );

                        if( iDataSetId != std::numeric_limits<size_t>::max( ) )
                            vBackgroundGridSeq.back( ) += pIndices->count(
                                iDataSetId,
                                { uiStartY, uiStartX, uiMapQMin, uiFromAnnoFilter, uiFromSameStrandFilter,
                                  uiFromYStrandFilter },
                                { uiEndY, uiEndX, uiMapQMax, uiToAnnoFilter, uiToSameStrandFilter, uiToYStrandFilter },
                                xIntersect,
                                0 );
                    }
                }
            }
    }

    END_RETURN;
}

bool PartialQuarry::normalizeGridSeq( )
{
    const bool bAxisIsCol = getValue<bool>( { "settings", "normalization", "grid_seq_axis_is_column" } );

    for( size_t uiI = 0; uiI < vvFlatValues.size( ); uiI++ )
    {
        CANCEL_RETURN;
        std::array<double, 2> vdVal = { (double)vvFlatValues[ uiI ][ 0 ], (double)vvFlatValues[ uiI ][ 1 ] };
        for( size_t uiK = 0; uiK < 2; uiK++ )
            vdVal[ uiK ] /= (double)std::max(
                (size_t)1,
                vBackgroundGridSeq[ bAxisIsCol ? uiI / vAxisCords[ 1 ].size( ) : uiI % vAxisCords[ 1 ].size( ) ] );

        vvNormalized.push_back( vdVal );
    }
    END_RETURN;
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
            icePreFilter( xData, bCol, 0, xData.vSliceBias[ bCol ? 0 : 1 ].size( ), uiI == 0 );
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

    const std::string sNorm = getValue<std::string>( { "settings", "normalization", "normalize_by" } );

    if( sNorm == "dont" )
        return doNotNormalize( );
    else if( sNorm == "radicl-seq" )
        return normalizeBinominalTest( );
    else if( sNorm == "rpm" )
        return normalizeSize( 1000000 );
    else if( sNorm == "rpk" )
        return normalizeSize( 1000 );
    else if( sNorm == "hi-c" )
        return normalizeIC( );
    else if( sNorm == "cool-hi-c" )
        return normalizeCoolIC( );
    else if( sNorm == "grid-seq" )
        return normalizeGridSeq( );
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

    const json& rJson = getValue<json>( { "coverage", "list" } );
    const std::string sByCol = getValue<std::string>( { "settings", "normalization", "divide_by_column_coverage" } );
    const size_t uiByCol =
        ( sByCol == "dont" ? std::numeric_limits<size_t>::max( ) : ( rJson.find( sByCol ) - rJson.begin( ) ) );
    const std::string sByRow = getValue<std::string>( { "settings", "normalization", "divide_by_row_coverage" } );
    const size_t uiByRow =
        ( sByRow == "dont" ? std::numeric_limits<size_t>::max( ) : ( rJson.find( sByRow ) - rJson.begin( ) ) );


    for( size_t uiI = 0; uiI < vCombined.size( ); uiI++ )
    {
        CANCEL_RETURN;

        double fVal = vCombined[ uiI ];

        if( uiByCol != std::numeric_limits<size_t>::max( ) )
        {
            if( vvCoverageValues[ 0 ][ uiByCol ][ uiI / vAxisCords[ 1 ].size( ) ] == 0 )
                fVal = std::numeric_limits<double>::quiet_NaN( );
            else
                fVal /= vvCoverageValues[ 0 ][ uiByCol ][ uiI / vAxisCords[ 1 ].size( ) ];
        }
        if( uiByRow != std::numeric_limits<size_t>::max( ) )
        {
            if( vvCoverageValues[ 1 ][ uiByRow ][ uiI % vAxisCords[ 1 ].size( ) ] == 0 || std::isnan( fVal ) )
                fVal = std::numeric_limits<double>::quiet_NaN( );
            else
                fVal /= vvCoverageValues[ 1 ][ uiByRow ][ uiI % vAxisCords[ 1 ].size( ) ];
        }

        vDivided.push_back( fVal );
    }

    END_RETURN;
}

const decltype( PartialQuarry::vDivided )
PartialQuarry::getDivided( const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::Divided, fPyPrint );
    return vDivided;
}

const decltype( PartialQuarry::vScaled )
PartialQuarry::getScaled( const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::Scaled, fPyPrint );
    return vScaled;
}

void PartialQuarry::regNormalization( )
{
    registerNode(
        NodeNames::Normalized,
        ComputeNode{ /*.sNodeName =*/"normalized_bins",
                     /*.fFunc =*/&PartialQuarry::setNormalized,
                     /*.vIncomingFunctions =*/
                     { NodeNames::FlatValues, NodeNames::RnaAssociatedBackground, NodeNames::RadiclSeqCoverage },
                     /*.vIncomingSession =*/
                     { { "settings", "normalization", "p_accept", "val" },
                       { "settings", "normalization", "ice_sparse_slice_filter", "val" },
                       { "contigs", "genome_size" } },
                     /*.vSessionsIncomingInPrevious =*/
                     { { "settings", "normalization", "normalize_by" },
                       { "settings", "normalization", "grid_seq_axis_is_column" },
                       { "settings", "normalization", "radicl_seq_samples", "val" },
                       { "settings", "normalization", "radicl_seq_axis_is_column" } } } );

    registerNode( NodeNames::GridSeqSamples,
                  ComputeNode{ /*.sNodeName =*/"grid_seq_samples",
                               /*.fFunc =*/&PartialQuarry::setGridSeqSamples,
                               /*.vIncomingFunctions =*/{ NodeNames::ActiveChromLength },
                               /*.vIncomingSession =*/
                               { { "annotation", "by_name" },
                                 { "dividend" },
                                 { "settings", "normalization", "normalize_by" },
                                 { "settings", "normalization", "grid_seq_annotation" },
                                 { "settings", "normalization", "grid_seq_samples", "val" },
                                 { "settings", "normalization", "grid_seq_axis_is_column" } },
                               /*.vSessionsIncomingInPrevious =*/{} } );

    registerNode( NodeNames::RadiclSeqSamples,
                  ComputeNode{ /*.sNodeName =*/"radicl_seq_samples",
                               /*.fFunc =*/&PartialQuarry::setRadiclSeqSamples,
                               /*.vIncomingFunctions =*/{ NodeNames::ActiveChromLength },
                               /*.vIncomingSession =*/
                               { { "dividend" },
                                 { "settings", "normalization", "normalize_by" },
                                 { "settings", "normalization", "radicl_seq_samples", "val" },
                                 { "settings", "normalization", "radicl_seq_axis_is_column" } },
                               /*.vSessionsIncomingInPrevious =*/{} } );
    // registerNode( NodeNames::ICESamples,
    //               ComputeNode{ /*.sNodeName =*/ "ice_samples",
    //                            /*.fFunc=*/ &PartialQuarry::setICESamples,
    //                            /*.vIncomingFunctions =*/ { NodeNames::ActiveChromLength },
    //                            /*.vIncomingSession =*/ { { "dividend" },
    //                                                  { "settings", "normalization", "normalize_by" },
    //                                                  { "settings", "normalization", "hi_c_samples", "val" } },
    //                            /*.vSessionsIncomingInPrevious =*/ {} } );

    registerNode( NodeNames::GridSeqCoverage,
                  ComputeNode{ /*.sNodeName =*/"grid_seq_coverage",
                               /*.fFunc =*/&PartialQuarry::setGridSeqCoverage,
                               /*.vIncomingFunctions =*/
                               { NodeNames::ActiveReplicates, NodeNames::MappingQuality, NodeNames::Directionality,
                                 NodeNames::GridSeqSamples, NodeNames::AxisCoords },
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
                                 { "dividend" } } } );
    registerNode(
        NodeNames::RadiclSeqCoverage,
        ComputeNode{
            /*.sNodeName =*/"radicl_coverage",
            /*.fFunc =*/&PartialQuarry::setRadiclSeqCoverage,
            /*.vIncomingFunctions =*/
            { NodeNames::ActiveReplicates, NodeNames::MappingQuality, NodeNames::AxisCoords,
              NodeNames::RadiclSeqSamples, NodeNames::IntersectionType, NodeNames::Directionality },
            /*.vIncomingSession =*/
            { { "replicates", "by_name" },
              { "settings", "normalization", "normalize_by" },
              { "settings", "normalization", "min_interactions", "val" } },
            /*.vSessionsIncomingInPrevious =*/{ { "settings", "normalization", "radicl_seq_axis_is_column" } } } );

    registerNode(
        NodeNames::RnaAssociatedGenesFilter,
        ComputeNode{ /*.sNodeName =*/"rna_associated_genes_filter",
                     /*.fFunc =*/&PartialQuarry::setRnaAssociatedGenesFilter,
                     /*.vIncomingFunctions =*/{ NodeNames::GridSeqCoverage },
                     /*.vIncomingSession =*/
                     { { "settings", "normalization", "grid_seq_rna_filter", "val_min" },
                       { "settings", "normalization", "grid_seq_rna_filter", "val_max" },
                       { "settings", "normalization", "grid_seq_dna_filter", "val_min" },
                       { "settings", "normalization", "grid_seq_dna_filter", "val_max" } },
                     /*.vSessionsIncomingInPrevious =*/{ { "settings", "normalization", "normalize_by" } } } );

    registerNode( NodeNames::RnaAssociatedBackground,
                  ComputeNode{ /*.sNodeName =*/"rna_associated_background",
                               /*.fFunc =*/&PartialQuarry::setRnaAssociatedBackground,
                               /*.vIncomingFunctions =*/{ NodeNames::RnaAssociatedGenesFilter },
                               /*.vIncomingSession =*/{ },
                               /*.vSessionsIncomingInPrevious =*/
                               { { "settings", "normalization", "normalize_by" },
                                 { "settings", "normalization", "grid_seq_axis_is_column" },
                                 { "settings", "normalization", "grid_seq_annotation" },
                                 { "annotation", "by_name" },
                                 { "replicates", "by_name" },
                                 { "dividend" } } } );

    registerNode( NodeNames::DistDepDecayRemoved,
                  ComputeNode{ /*.sNodeName =*/"dist_dep_dec_normalized_bins",
                               /*.fFunc =*/&PartialQuarry::setDistDepDecayRemoved,
                               /*.vIncomingFunctions =*/{ NodeNames::Normalized, NodeNames::FlatDecay },
                               /*.vIncomingSession =*/{ },
                               /*.vSessionsIncomingInPrevious =*/{ { "settings", "normalization", "ddd" } } } );

    registerNode( NodeNames::Divided,
                  ComputeNode{ /*.sNodeName =*/"divided_by_tracks",
                               /*.fFunc =*/&PartialQuarry::setDivided,
                               /*.vIncomingFunctions =*/{ NodeNames::Combined },
                               /*.vIncomingSession =*/
                               { { "settings", "normalization", "divide_by_column_coverage" },
                                 { "settings", "normalization", "divide_by_row_coverage" },
                                 { "coverage", "list" } },
                               /*.vSessionsIncomingInPrevious =*/{} } );
}

} // namespace cm