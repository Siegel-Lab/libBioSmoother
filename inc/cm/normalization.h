#include "cm/partial_quarry.h"
#include "partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

bool PartialQuarry::normalizeSize( size_t uiSize )
{
    for( size_t uiY = 0; uiY < 3; uiY++ )
        for( auto& vVal : vvFlatValues[ uiY ] )
        {
            CANCEL_RETURN;
            vvNormalized[ uiY ].push_back( std::array<double, 2>{
                vvFlatTotal[ uiY ][ 0 ] == 0 ? 0 : (double)uiSize * (double)vVal[ 0 ] / (double)vvFlatTotal[ uiY ][ 0 ],
                vvFlatTotal[ uiY ][ 1 ] == 0 ? 0
                                             : (double)uiSize * (double)vVal[ 1 ] / (double)vvFlatTotal[ uiY ][ 1 ] } );
        }
    END_RETURN;
}

bool PartialQuarry::doNotNormalize( )
{
    for( size_t uiY = 0; uiY < 3; uiY++ )
        for( auto& vVal : vvFlatValues[ uiY ] )
        {
            CANCEL_RETURN;
            vvNormalized[ uiY ].push_back( std::array<double, 2>{ (double)vVal[ 0 ], (double)vVal[ 1 ] } );
        }
    END_RETURN;
}


bool PartialQuarry::normalizeBinominalTest( )
{
    const bool bIsCol = getValue<bool>( { "settings", "normalization", "radicl_seq_axis_is_column" } );
    const size_t uiNumBinsInRowTotal =
        ( vCanvasSize[ bIsCol ? 1 : 0 ] - 1 ) / ( bIsCol ? uiBinHeight.r( ) : uiBinWidth.r( ) ) + 1;
    const size_t uiNumSaples = getValue<size_t>( { "settings", "normalization", "radicl_seq_samples", "val" } );
    for( size_t uiY = 0; uiY < 3; uiY++ )
        vvNormalized[ uiY ] = normalizeBinominalTestTrampoline(
            vvFlatValues[ uiY ], vRadiclSeqCoverage, vRadiclSeqNumNonEmptyBins, uiNumSaples, uiNumBinsInRowTotal,
            getValue<double>( { "settings", "normalization", "p_accept", "val" } ), bIsCol,
            uiY == 2 ? vV4cCoords[ 1 ].size( ) : vAxisCords[ 1 ].size( ) );
    CANCEL_RETURN;
    END_RETURN;
}

/*
original ICing implementation:
@ https://github.com/open2c/cooler/blob/3d284485df070255f4a904102178b186c14c1cee/cooler/balance.py
@ https://github.com/open2c/cooler/blob/3d284485df070255f4a904102178b186c14c1cee/cooler/tools.py

*/


size_t PartialQuarry::iceGetCount( IceData& rIceData, size_t uiX, size_t uiY, size_t uiY_, bool bA )
{
    assert( uiX < rIceData.vSliceBias[ 0 ].size( ) );
    assert( uiY < rIceData.vSliceBias[ 1 ].size( ) );
    size_t uiIdx = uiY + uiX * ( rIceData.vSliceBias[ 1 ].size( ) );
    assert( uiIdx <
            vvFlatValues[ uiY_ ].size( ) ); // @todo this assert triggered with max coverage per bin columns active
    return vvFlatValues[ uiY_ ][ uiIdx ][ bA ? 0 : 1 ];
}

void PartialQuarry::iceFilter( IceData& /*rIceData*/, size_t /*uiFrom*/, size_t /*uiTo*/ )
{}

void PartialQuarry::icePreFilter( IceData& rIceData, bool bCol, size_t uiFrom, size_t uiTo, size_t uiY, bool bA )
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
                if( iceGetCount( rIceData, bCol ? uiI : uiJ, bCol ? uiJ : uiI, uiY, bA ) > 0 )
                    ++uiCnt;
            if( (double)uiCnt <= uiWSlice * fFilter )
                rIceData.vSliceBias[ bCol ? 0 : 1 ][ uiI ] = 0;
        }
    }
}


void PartialQuarry::iceTimesOuterProduct( IceData& rIceData, bool bA, size_t uiFrom, size_t uiTo, size_t uiY )
{
    const size_t uiH = rIceData.vSliceBias[ 1 ].size( );
    for( size_t uiI = uiFrom; uiI < uiTo; uiI++ )
    {
        assert( uiI / uiH < rIceData.vSliceBias[ 0 ].size( ) );
        assert( uiI % uiH < rIceData.vSliceBias[ 1 ].size( ) );
        assert( uiI < rIceData.vBiases.size( ) );
        rIceData.vBiases[ uiI ] = rIceData.vSliceBias[ 0 ][ uiI / uiH ] * rIceData.vSliceBias[ 1 ][ uiI % uiH ] *
                                  iceGetCount( rIceData, uiI / uiH, uiI % uiH, uiY, bA );
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

// @todo-low-prio ICe samples could also be computed using this...
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

    getSamples( GetSamplesMode::Bins, uiRadiclSeqSamples, "", bAxisIsCol ? uiBinHeight.r( ) : uiBinWidth.r( ),
                bAxisIsCol, false, vRadiclSeqSamples );
    // then this cancel_return is necessary to catch the cancelled getSamples
    CANCEL_RETURN;


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

    for( size_t uiY = 0; uiY < 3; uiY++ )
    {
        const auto& rYCoords = uiY == 2 ? vV4cCoords[ 1 ] : vAxisCords[ 1 ];
        for( size_t uiI = 0; uiI < vvFlatValues[ uiY ].size( ); uiI++ )
        {
            CANCEL_RETURN;
            std::array<double, 2> vdVal = { (double)vvFlatValues[ uiY ][ uiI ][ 0 ],
                                            (double)vvFlatValues[ uiY ][ uiI ][ 1 ] };

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
        vvNormalized[ 0 ].resize( vvFlatValues[ 0 ].size( ) );
        for( size_t uiI = 0; uiI < 2; uiI++ )
        {
            std::vector<size_t> vCnt;
            for( auto& rFlat : vvFlatValues[ 0 ] )
                vCnt.push_back( rFlat[ uiI ] );
            auto vRet = normalizeCoolerTrampoline( vCnt, vAxisCords[ 0 ].size( ) );
            for( size_t uiX = 0; uiX < vRet.size( ); uiX++ )
                vvNormalized[ 0 ][ uiX ][ uiI ] = vRet[ uiX ];
        }

        for( size_t uiY = 1; uiY < 3; uiY++ )
        {
            vvNormalized[ uiY ].resize( vvFlatValues[ uiY ].size( ) );
            for( auto& vVal : vvFlatValues[ uiY ] )
            {
                CANCEL_RETURN;
                vvNormalized[ uiY ].push_back( std::array<double, 2>{ (double)vVal[ 0 ], (double)vVal[ 1 ] } );
            }
        }
    }
    CANCEL_RETURN;
    END_RETURN;
}

bool PartialQuarry::normalizeIC( )
{
    for( size_t uiY = 3; uiY < NUM_COORD_SYSTEMS; uiY++ )
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
                icePreFilter( xData, bCol, 0, xData.vSliceBias[ bCol ? 0 : 1 ].size( ), uiY, uiI == 0 );
            for( size_t uiItr = 0; uiItr < uiMaxIters; uiItr++ )
            {
                CANCEL_RETURN;
                iceFilter( xData, 0, xData.vBiases.size( ) );
                iceTimesOuterProduct( xData, uiI == 0, 0, xData.vBiases.size( ), uiY );
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
            for( size_t uiJ = 0; uiJ < vvFlatValues[ uiY ].size( ); uiJ++ )
            {
                CANCEL_RETURN;
                uiMaxFlat = std::max( uiMaxFlat, vvFlatValues[ uiY ][ uiJ ][ uiI ] );
            }
            if( uiMaxFlat > 0 )
            {
                if( vVar[ 0 ] >= fTol || vVar[ 1 ] >= fTol )
                {
                    setError( "iterative correction did not converge (var=" + std::to_string( vVar[ 0 ] ) + ", " +
                              std::to_string( vVar[ 1 ] ) + " mean=" + std::to_string( vMean[ 0 ] ) + ", " +
                              std::to_string( vMean[ 1 ] ) + "), showing data anyways" );
                }
                else if( iceMaxBias( xData, true ) == 0 || iceMaxBias( xData, false ) == 0 )
                {
                    setError( "iterative correction converged to zero, showing un-normalized data" );
                    for( size_t uiJ = 0; uiJ < xData.vSliceBias[ 4 - uiY ].size( ); uiJ++ )
                    {
                        CANCEL_RETURN;
                        xData.vSliceBias[ 4 - uiY ][ uiJ ] = 1;
                    }
                }
            }
            else
            {
                for( size_t uiJ = 0; uiJ < xData.vSliceBias[ 4 - uiY ].size( ); uiJ++ )
                {
                    CANCEL_RETURN;
                    xData.vSliceBias[ 4 - uiY ][ uiJ ] = 1;
                }
            }
            vIceSliceBias[ 4 - uiY ][ uiI ].swap( xData.vSliceBias[ 4 - uiY ] );
        }
    }

    {
        const auto& rYCoords = pickYCoords( 0 );
        size_t uiH = rYCoords.size( );
        vvNormalized[ 0 ].resize( vvFlatValues[ 0 ].size( ) );
        for( size_t uiA : { 0, 1 } )
            for( size_t uiI = 0; uiI < vvFlatValues[ 0 ].size( ); uiI++ )
                vvNormalized[ 0 ][ uiI ][ uiA ] = vIceSliceBias[ 0 ][ uiA ][ uiI / uiH ] //
                                                  * vIceSliceBias[ 1 ][ uiA ][ uiI % uiH ] //
                                                  * (double)vvFlatValues[ 0 ][ uiI ][ uiA ];
    }

    {
        const auto& rYCoords = pickYCoords( 0 );
        size_t uiH = rYCoords.size( );
        vvNormalized[ 1 ].resize( vvFlatValues[ 1 ].size( ) );
        for( size_t uiA : { 0, 1 } )
            for( size_t uiI = 0; uiI < vvFlatValues[ 1 ].size( ); uiI++ )
                vvNormalized[ 1 ][ uiI ][ uiA ] = vIceSliceBias[ 0 ][ uiA ][ uiI / uiH ] //
                                                  * (double)vvFlatValues[ 1 ][ uiI ][ uiA ];
    }

    {
        const auto& rYCoords = pickYCoords( 0 );
        size_t uiH = rYCoords.size( );
        vvNormalized[ 2 ].resize( vvFlatValues[ 2 ].size( ) );
        for( size_t uiA : { 0, 1 } )
            for( size_t uiI = 0; uiI < vvFlatValues[ 2 ].size( ); uiI++ )
                vvNormalized[ 1 ][ uiI ][ uiA ] = vIceSliceBias[ 1 ][ uiA ][ uiI % uiH ] //
                                                  * (double)vvFlatValues[ 2 ][ uiI ][ uiA ];
    }
    END_RETURN;
}

bool PartialQuarry::setNormalized( )
{
    for( size_t uiY = 0; uiY < 3; uiY++ )
    {
        vvNormalized[ uiY ].clear( );
        vvNormalized[ uiY ].reserve( vvFlatValues[ uiY ].size( ) );
    }
    for( size_t uiY = 0; uiY < 2; uiY++ )
        for( size_t uiI = 0; uiI < 2; uiI++ )
            vIceSliceBias[ uiY ][ uiI ].clear( );

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
        for( size_t uiY = 0; uiY < 5; uiY++ )
            for( size_t uiI = 0; uiI < vvNormalized[ uiY ].size( ); uiI++ )
                for( size_t uiJ = 0; uiJ < 2; uiJ++ )
                    if( vBinCoords[ uiY ][ uiI ][ uiJ ].uiDecayCoordIndex != std::numeric_limits<size_t>::max( ) )
                    {
                        CANCEL_RETURN;
                        if( vvFlatDecay[ uiY ][ vBinCoords[ 0 ][ uiI ][ uiJ ].uiDecayCoordIndex ][ uiJ ] > 0 )
                            vvNormalized[ uiY ][ uiI ][ uiJ ] /=
                                vvFlatDecay[ uiY ][ vBinCoords[ 0 ][ uiI ][ uiJ ].uiDecayCoordIndex ][ uiJ ];
                        else
                            vvNormalized[ uiY ][ uiI ][ uiJ ] = 0;
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
    registerNode(
        NodeNames::Normalized,
        ComputeNode{ /*.sNodeName =*/"normalized_bins",
                     /*.fFunc =*/&PartialQuarry::setNormalized,
                     /*.vIncomingFunctions =*/
                     { NodeNames::FlatValues, NodeNames::RnaAssociatedBackground, NodeNames::RadiclSeqCoverage },
                     /*.vIncomingSession =*/
                     { { "settings", "normalization", "p_accept", "val" },
                       { "settings", "normalization", "ice_sparse_slice_filter", "val" } },
                     /*.vSessionsIncomingInPrevious =*/
                     { { "contigs", "genome_size" },
                       { "settings", "normalization", "normalize_by" },
                       { "settings", "normalization", "grid_seq_axis_is_column" },
                       { "settings", "normalization", "ice_local" },
                       { "settings", "normalization", "radicl_seq_samples", "val" },
                       { "settings", "normalization", "radicl_seq_axis_is_column" } },
                     /*bHidden =*/false } );

    registerNode( NodeNames::GridSeqSamples,
                  ComputeNode{ /*.sNodeName =*/"grid_seq_samples",
                               /*.fFunc =*/&PartialQuarry::setGridSeqSamples,
                               /*.vIncomingFunctions =*/{ NodeNames::ActiveChromLength },
                               /*.vIncomingSession =*/
                               { { "annotation", "by_name" },
                                 { "settings", "normalization", "normalize_by" },
                                 { "settings", "normalization", "grid_seq_annotation" },
                                 { "settings", "normalization", "grid_seq_samples", "val" },
                                 { "settings", "normalization", "grid_seq_axis_is_column" } },
                               /*.vSessionsIncomingInPrevious =*/{ { "dividend" } },
                               /*bHidden =*/false } );

    registerNode( NodeNames::RadiclSeqSamples,
                  ComputeNode{ /*.sNodeName =*/"radicl_seq_samples",
                               /*.fFunc =*/&PartialQuarry::setRadiclSeqSamples,
                               /*.vIncomingFunctions =*/{ NodeNames::ActiveChromLength },
                               /*.vIncomingSession =*/
                               { { "settings", "normalization", "normalize_by" },
                                 { "settings", "normalization", "radicl_seq_samples", "val" },
                                 { "settings", "normalization", "radicl_seq_axis_is_column" } },
                               /*.vSessionsIncomingInPrevious =*/{ { "dividend" } },
                               /*bHidden =*/false } );
    // registerNode( NodeNames::ICESamples,
    //               ComputeNode{ /*.sNodeName =*/ "ice_samples",
    //                            /*.fFunc=*/ &PartialQuarry::setICESamples,
    //                            /*.vIncomingFunctions =*/ { NodeNames::ActiveChromLength },
    //                            /*.vIncomingSession =*/ { { "dividend" },
    //                                                  { "settings", "normalization", "normalize_by" },
    //                                                  { "settings", "normalization", "hi_c_samples", "val" } },
    //                            /*.vSessionsIncomingInPrevious =*/ {},
    //                           /*bHidden =*/false } );

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
                                 { "dividend" } },
                               /*bHidden =*/false } );
    registerNode( NodeNames::RadiclSeqCoverage,
                  ComputeNode{ /*.sNodeName =*/"radicl_coverage",
                               /*.fFunc =*/&PartialQuarry::setRadiclSeqCoverage,
                               /*.vIncomingFunctions =*/
                               { NodeNames::ActiveReplicates, NodeNames::MappingQuality, NodeNames::AxisCoords,
                                 NodeNames::RadiclSeqSamples, NodeNames::IntersectionType, NodeNames::Directionality },
                               /*.vIncomingSession =*/
                               { { "replicates", "by_name" },

                                 { "settings", "normalization", "min_interactions", "val" } },
                               /*.vSessionsIncomingInPrevious =*/
                               { { "settings", "normalization", "normalize_by" },
                                 { "settings", "normalization", "radicl_seq_axis_is_column" } },
                               /*bHidden =*/false } );

    registerNode( NodeNames::RnaAssociatedGenesFilter,
                  ComputeNode{ /*.sNodeName =*/"rna_associated_genes_filter",
                               /*.fFunc =*/&PartialQuarry::setRnaAssociatedGenesFilter,
                               /*.vIncomingFunctions =*/{ NodeNames::GridSeqCoverage },
                               /*.vIncomingSession =*/
                               { { "settings", "normalization", "grid_seq_rna_filter", "val_min" },
                                 { "settings", "normalization", "grid_seq_rna_filter", "val_max" },
                                 { "settings", "normalization", "grid_seq_dna_filter", "val_min" },
                                 { "settings", "normalization", "grid_seq_dna_filter", "val_max" } },
                               /*.vSessionsIncomingInPrevious =*/{ { "settings", "normalization", "normalize_by" } },
                               /*bHidden =*/false } );

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
                                 { "dividend" } },
                               /*bHidden =*/false } );

    registerNode( NodeNames::DistDepDecayRemoved,
                  ComputeNode{ /*.sNodeName =*/"dist_dep_dec_normalized_bins",
                               /*.fFunc =*/&PartialQuarry::setDistDepDecayRemoved,
                               /*.vIncomingFunctions =*/{ NodeNames::Normalized, NodeNames::FlatDecay },
                               /*.vIncomingSession =*/{ },
                               /*.vSessionsIncomingInPrevious =*/{ { "settings", "normalization", "ddd" } },
                               /*bHidden =*/false } );
}

} // namespace cm