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

@ https://github.com/mirnylab/hiclib-legacy/
@ https://github.com/mirnylab/hiclib-legacy/blob/518546e41987dca8a40f45ddc63601a5aaf46bfa/src/hiclib/highResBinnedData.py#L565

def _marginalError(self, marginals=None, percentile=99.9):
    """Checks after each pass of IC if marginals are close enough to 1.
    The error is calculated as the specified percentile of deviation
    from the mean marginal.
    """
    if marginals is None:
        if not hasattr(self, "marginals"):
            return 99999
        marginals = self.marginals
    marginals = np.concatenate(marginals)
    marginals = marginals[marginals != 0] # removes all zero values
    error = np.percentile(np.abs(marginals - marginals.mean()), percentile)
    return error / marginals.mean()

def getMarginals(self, normalizeForIC=False):
    """
    Returns a sum over each row/column, and saves them to self.marginals
    normalizeForIc=True will normalize them to mean 1.
    """
    self._hasData()
    marginals = [np.zeros(i, float) for i in self.genome.chrmLensBin]
    for chr1, chr2 in self.data:
        m2, m1 = self.data[(chr1, chr2)].getSums()
        marginals[chr1] += m1
        if chr1 != chr2:
            marginals[chr2] += m2

    self.marginals = marginals
    return marginals

def divideByVectors(self, vectors):
    """
    Divides each row and column by correspoinding
        value from vectors[0] and vectors[1]
    Does it without calling getData twice!
    """

    vecX = vectors[0]
    vecY = vectors[1]
    vecX[vecX == 0] = 1
    vecY[vecY == 0] = 1
    data = self.getData()
    assert data.shape[1] == len(vecX)
    assert data.shape[0] == len(vecY)
    data /= vecX[None, :] # puts entire array into new array
    data /= vecY[:, None] # puts each element of the array into a new array
    self.setData(data)

def divideByMarginals(self, marginals=None):
    """Divides matrix by a vector.
    If the vector is not provided, will divide it by a
    marginals calculated previously.
    """
    self._hasData()
    if marginals is None:
        marginals = self.marginals

    for chr1, chr2 in self.data:
        m2 = marginals[chr1]
        m1 = marginals[chr2]
        self.data[(chr1, chr2)].divideByVectors((np.sqrt(m1), np.sqrt(m2)))

def iterativeCorrection(self, tolerance=1e-2):
    """
    Performs iterative correction in place.
    Tolerance is the maximum allowed relative deviation of the marginals.
    Tries to set biases to self.biases after IC.
    """
    self._hasData()
    curPass = 0
    marginals = np.ones(self.genome.numBins, float)
    while self._marginalError() > tolerance:
        m = self.getMarginals(normalizeForIC=True)
        marginals *= np.concatenate(m)
        self.divideByMarginals()
        print("Pass = %d, Error = %lf" % (curPass, self._marginalError()))
        curPass += 1
    self.biases = marginals
    return self.biases

*/
bool PartialQuarry::normalizeIC( )
{
    for( auto& vVal : vvFlatValues )
    {
        CANCEL_RETURN;
        vvNormalized.push_back( std::array<double, 2>{ (double)vVal[ 0 ], (double)vVal[ 1 ] } );
    }
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