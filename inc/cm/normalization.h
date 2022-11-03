#include "cm/partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

bool PartialQuarry::normalizeSize( size_t uiSize )
{
    std::array<size_t, 2> vTotalReads = { 0, 0 };

    for( size_t uiJ = 0; uiJ < 2; uiJ++ )
        for( size_t uiX : vInGroup[ uiJ ] )
        {
            CANCEL_RETURN;
            vTotalReads[ uiJ ] +=
                this->xSession[ "replicates" ][ "by_name" ][ this->vActiveReplicates[ uiX ] ][ "total_reads" ]
                    .get<size_t>( );
        }
    for( auto& vVal : vvFlatValues )
    {
        CANCEL_RETURN;
        vvNormalized.push_back( std::array<double, 2>{
            vTotalReads[ 0 ] == 0 ? 0 : (double)uiSize * (double)vVal[ 0 ] / (double)vTotalReads[ 0 ],
            vTotalReads[ 1 ] == 0 ? 0 : (double)uiSize * (double)vVal[ 1 ] / (double)vTotalReads[ 1 ] } );
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
    size_t uiNumBinsInRowTotal = ( this->xSession[ "contigs" ][ "genome_size" ].get<size_t>( ) - 1 ) / uiBinWidth + 1;
    vvNormalized = normalizeBinominalTestTrampoline(
        vvFlatValues, vFlatNormValues, 
        uiNumBinsInRowTotal,
        this->xSession[ "settings" ][ "normalization" ][ "p_accept" ][ "val" ].get<double>( ) );
    CANCEL_RETURN;
    END_RETURN;
}

bool PartialQuarry::normalizeIC( )
{
    /*

    def hi_c_normalization(self, bins, cols, rows):
        # max 50 iterations
        for num_itr in range(50):
            assert len(bins) == len(cols) * len(rows)
            # compute coverage & width
            cov_rows = [(sum(bins[y + len(cols) * x] for y in range(len(cols))), v[1]) for x, v in enumerate(rows)]
            cov_cols = [(sum(bins[y + len(cols) * x] for x in range(len(rows))), v[1]) for y, v in enumerate(cols)]
            def unit_mean(l):
                total_width = sum(w for v, w in l if v > 0)
                cnt = 0
                m = 1
                for m, w in sorted([(v, w) for v, w in l if v > 0]):
                    if cnt + w > total_width / 2:
                        break
                    cnt += w
                return [r if r != 0 else 1 for r in [x / m for x, _ in l]]

            cov_cols = unit_mean(cov_cols)
            cov_rows = unit_mean(cov_rows)
            assert len(cov_cols) == len(cols)
            assert len(cov_rows) == len(rows)
            assert len(bins) == len(cov_cols) * len(cov_rows)

            max_bias_delta = 0
            for idx in range(len(bins)):
                bias_delta = cov_cols[idx % len(cols)] * cov_rows[idx // len(cols)]
                bins[idx] = bins[idx] / bias_delta
                max_bias_delta = max(max_bias_delta, abs(1-bias_delta))

            if max_bias_delta < 0.01:
                print("stopped at iteration", num_itr, "since max_bias_delta is", max_bias_delta)
                break

        n = max(bins + [1])
        return [x/n for x in bins]

    */
    throw std::logic_error( "Function not implemented" );
}

bool PartialQuarry::setNormalized( )
{
    vvNormalized.clear( );
    vvNormalized.reserve( vvFlatValues.size( ) );

    if( this->xSession[ "settings" ][ "normalization" ][ "normalize_by" ].get<std::string>( ) == "dont" )
        return doNotNormalize( );
    else if( this->xSession[ "settings" ][ "normalization" ][ "normalize_by" ].get<std::string>( ) == "radicl-seq" )
        return normalizeBinominalTest( );
    else if( this->xSession[ "settings" ][ "normalization" ][ "normalize_by" ].get<std::string>( ) == "rpm" )
        return normalizeSize( 1000000 );
    else if( this->xSession[ "settings" ][ "normalization" ][ "normalize_by" ].get<std::string>( ) == "rpk" )
        return normalizeSize( 1000 );
    else if( this->xSession[ "settings" ][ "normalization" ][ "normalize_by" ].get<std::string>( ) == "hi-c" )
        return normalizeIC( );
    else
        throw std::logic_error( "invalid value for normalize_by" );
}

bool PartialQuarry::setDivided( )
{
    vDivided.clear( );
    vDivided.reserve( vCombined.size( ) );

    const bool bByCol = this->xSession[ "settings" ][ "normalization" ][ "divide_by_column_coverage" ].get<bool>( ) &&
                        vvFlatCoverageValues[ 0 ].size( ) != 0;
    const bool bByRow = this->xSession[ "settings" ][ "normalization" ][ "divide_by_row_coverage" ].get<bool>( ) &&
                        vvFlatCoverageValues[ 1 ].size( ) != 0;

    for( size_t uiI = 0; uiI < vCombined.size( ); uiI++ )
    {
        CANCEL_RETURN;

        double fVal = vCombined[ uiI ];

        if( bByCol )
        {
            if( vvFlatCoverageValues[ 0 ][ uiI / vAxisCords[ 0 ].size( ) ] == 0 )
                fVal = std::numeric_limits<double>::quiet_NaN( );
            else
                fVal /= vvFlatCoverageValues[ 0 ][ uiI / vAxisCords[ 0 ].size( ) ];
        }
        if( bByRow )
        {
            if( vvFlatCoverageValues[ 1 ][ uiI % vAxisCords[ 0 ].size( ) ] == 0 || std::isnan( fVal ) )
                fVal = std::numeric_limits<double>::quiet_NaN( );
            else
                fVal /= vvFlatCoverageValues[ 1 ][ uiI % vAxisCords[ 0 ].size( ) ];
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
                               .vIncomingFunctions = { NodeNames::FlatValues, NodeNames::FlatCoverageValues },
                               .vIncomingSession = { { "settings", "normalization", "normalize_by" },
                                                     { "replicates", "by_name" },
                                                     { "settings", "normalization", "p_accept", "val" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::Divided,
                  ComputeNode{ .sNodeName = "divided_by_tracks",
                               .fFunc = &PartialQuarry::setDivided,
                               .vIncomingFunctions = { NodeNames::Combined, NodeNames::FlatCoverageValues },
                               .vIncomingSession = { { "settings", "normalization", "divide_by_column_coverage" },
                                                     { "settings", "normalization", "divide_by_row_coverage" } },
                               .uiLastUpdated = uiCurrTime } );
}

} // namespace cm