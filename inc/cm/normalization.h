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
        vvFlatValues, vFlatNormValues[1], uiNumBinsInRowTotal,
        this->xSession[ "settings" ][ "normalization" ][ "p_accept" ][ "val" ].get<double>( ) );
    CANCEL_RETURN;
    END_RETURN;
}

/*

PLAN:

    generate remainders of each column and row
    generate corner remainder
    repeat x times
        generate marginals (i.e. row & cols sums)
        divide bias by marginal

@ https://github.com/open2c/cooler/blob/3d284485df070255f4a904102178b186c14c1cee/cooler/balance.py
@ https://github.com/open2c/cooler/blob/3d284485df070255f4a904102178b186c14c1cee/cooler/tools.py


*/

struct ICBias
{
    double fVal;
    size_t uiX, uiY;
    size_t uiW, uiH;
};

struct ICMargin
{
    double fVal;
    size_t uiP;
    size_t uiS;
};

bool PartialQuarry::normalizeIC( )
{
    if( vAxisCords[0].size() != vAxisCords[1].size() )
        doNotNormalize();
    else
    {
        vvNormalized.resize(vvFlatValues.size());
        for(size_t uiI = 0; uiI < 2; uiI++)
        {
            std::vector<size_t> vCnt;
            for(auto& rFlat : vvFlatValues)
                vCnt.push_back(rFlat[uiI]);
            auto vRet = normalizeCoolerTrampoline( vCnt, vAxisCords[0].size() );
            for(size_t uiX = 0; uiX < vRet.size(); uiX++)
                vvNormalized[uiX][uiI] = vRet[uiX];
        }
    }
    CANCEL_RETURN;
    END_RETURN;
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
                               .vIncomingSession = { { "replicates", "by_name" },
                                                     { "settings", "normalization", "p_accept", "val" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::Divided,
                  ComputeNode{ .sNodeName = "divided_by_tracks",
                               .fFunc = &PartialQuarry::setDivided,
                               .vIncomingFunctions = { NodeNames::Combined },
                               .vIncomingSession = { { "settings", "normalization", "divide_by_column_coverage" },
                                                     { "settings", "normalization", "divide_by_row_coverage" } },
                               .uiLastUpdated = uiCurrTime } );
}

} // namespace cm