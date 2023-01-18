#include "cm/partial_quarry.h"
#include <cmath>

#pragma once

namespace cm
{

/**
 * @brief
 *
 * iSmallerBins
 * - 0: add smaller bin
 * - 1: add smaller bin only if there is at least one full sized bin
 * - 2: make last bin larger
 * - 3: skip bin
 * - 4: fit chromosome size to bin size (smaller)
 * - 5: fit chromosome size to bin size (larger)
 *
 * @param uiBinSize
 * @param uiScreenStartPos
 * @param uiScreenEndPos
 * @param bEmptyBinForChromBorder
 * @param iSmallerBins
 * @param vChromosomes
 * @return std::vector<AxisCoord>
 */
std::pair<std::vector<AxisCoord>, std::vector<AxisRegion>> axisCoordsHelper( size_t uiBinSize, size_t uiScreenStartPos,
                                                                             size_t uiScreenEndPos, size_t iSmallerBins,
                                                                             std::vector<ChromDesc> vChromosomes,
                                                                             bool& bCancel )
{
    std::vector<AxisCoord> vRet;
    std::vector<AxisRegion> vRet2;
    vRet.reserve( 2 * ( uiScreenEndPos - uiScreenStartPos ) / uiBinSize );
    size_t uiCurrScreenPos = uiScreenStartPos;
    size_t uiChromosomeStartPos = 0;
    size_t uiChr = 0;
    for( size_t uiI = 0; uiI < vChromosomes.size( ); uiI++ )
    {
        size_t uiChromosomeEndPos = uiChromosomeStartPos + vChromosomes[ uiI ].uiLength;
        size_t uiItrEndPos = std::min( uiScreenEndPos, uiChromosomeEndPos );
        size_t uiStartIdx = vRet.size( );
        size_t uiStartScreenPos = uiCurrScreenPos;
        size_t uiStartChromPos = uiCurrScreenPos >= uiChromosomeStartPos ? uiCurrScreenPos - uiChromosomeStartPos : 0;
        while( uiCurrScreenPos >= uiChromosomeStartPos && uiCurrScreenPos < uiItrEndPos )
        {
            if( bCancel )
                return std::make_pair( vRet, vRet2 );
            size_t uiIndexPos = uiCurrScreenPos - uiChromosomeStartPos;
            size_t uiCurrBinSize;
            switch( iSmallerBins )
            {
                case 0:
                case 1:
                    assert( vChromosomes[ uiI ].uiLength >= uiIndexPos );
                    uiCurrBinSize = std::min( uiBinSize, vChromosomes[ uiI ].uiLength - uiIndexPos );
                    break;
                case 2:
                    assert( uiIndexPos + 2 * uiBinSize <= vChromosomes[ uiI ].uiLength ||
                            vChromosomes[ uiI ].uiLength >= uiIndexPos );
                    uiCurrBinSize = uiIndexPos + 2 * uiBinSize > vChromosomes[ uiI ].uiLength
                                        ? vChromosomes[ uiI ].uiLength - uiIndexPos
                                        : uiBinSize;
                    break;
                case 3:
                case 4:
                case 5:
                    uiCurrBinSize = uiBinSize;
                    break;
                default:
                    throw std::logic_error( "unknown iSmallerBins value" );
            }
            bool bAddBin;
            switch( iSmallerBins )
            {
                case 0:
                case 2:
                case 5:
                    bAddBin = true;
                    break;
                case 1:
                    bAddBin = uiCurrBinSize == uiBinSize || uiIndexPos > 0;
                    break;
                case 3:
                case 4:
                    bAddBin = uiIndexPos + uiCurrBinSize <= uiChromosomeEndPos;
                    break;
                default:
                    throw std::logic_error( "unknown iSmallerBins value" );
            }
            if( bAddBin )
                vRet.push_back( AxisCoord{
                    .uiChromosome = uiI, //
                    .uiScreenPos = uiCurrScreenPos, //
                    .uiIndexPos = uiIndexPos, //
                    .uiScreenSize = uiCurrBinSize, //
                    .uiIndexSize = uiCurrBinSize, //
                    .uiRegionIdx = uiChr //
                } );
            bool bIncScreenPos;
            switch( iSmallerBins )
            {
                case 0:
                case 1:
                case 2:
                case 3:
                case 5:
                    bIncScreenPos = true;
                    break;
                case 4:
                    bIncScreenPos = uiIndexPos + uiCurrBinSize <= uiChromosomeEndPos;
                    break;
                default:
                    throw std::logic_error( "unknown iSmallerBins value" );
            }
            if( bIncScreenPos )
                uiCurrScreenPos += uiCurrBinSize;
            uiIndexPos += uiCurrBinSize;
        }

        if( vRet.size( ) > uiStartIdx )
        {
            assert( uiItrEndPos >= uiStartChromPos );
            vRet2.push_back( AxisRegion{
                {
                    .uiChromosome = uiI, //
                    .uiScreenPos = uiStartScreenPos, //
                    .uiIndexPos = uiStartChromPos, //
                    .uiScreenSize = uiCurrScreenPos - uiStartScreenPos, //
                    .uiIndexSize = uiItrEndPos - uiStartChromPos, //
                    .uiRegionIdx = uiChr //
                },
                .uiCoordStartIdx = uiStartIdx, //
                .uiNumCoords = vRet.size( ) - uiStartIdx //
            } );
        }

        uiChromosomeStartPos = uiChromosomeEndPos;

        if( uiScreenEndPos <= uiChromosomeStartPos )
            break;
        ++uiChr;
    }

    return std::make_pair( vRet, vRet2 );
}
size_t smaller_bin_to_num( std::string sVal )
{
    if( sVal == "smaller" )
        return 0;
    if( sVal == "smaller_if_fullsized_exists" )
        return 1;
    if( sVal == "larger" )
        return 2;
    if( sVal == "skip" )
        return 3;
    if( sVal == "fit_chrom_smaller" )
        return 4;
    if( sVal == "fit_chrom_larger" )
        return 5;
    throw std::logic_error( "unknown iSmallerBins value" );
}


size_t getEvenlyDividableByMaxTwoPowNIn( size_t uiFrom, size_t uiTo )
{
    size_t uiN = std::log2( uiTo - uiFrom ) + 1;
    while( uiN > 0 )
    {
        size_t uiX = 1 << uiN;
        if( uiFrom % uiX == 0 )
        {
            assert( uiFrom + uiX >= uiTo );
            return uiFrom;
        }
        size_t uiY = uiX * ( uiFrom / uiX ) + uiX;
        if( uiY < uiTo )
        {
            assert( uiY + uiX >= uiTo );
            assert( uiY >= uiFrom );
            return uiY;
        }
        --uiN;
    }
    return std::numeric_limits<size_t>::max( );
}

template <typename anno_t>
std::pair<std::vector<AxisCoord>, std::vector<AxisRegion>>
annoCoordsHelper( size_t uiBinSize, size_t uiScreenStartPos, size_t uiScreenEndPos, size_t /*iSmallerBins*/,
                  size_t iMultipleAnnosInBin, size_t iAnnoInMultipleBins, std::vector<ChromDesc> vChromosomes,
                  bool& bCancel, const json rJson, anno_t& rAnno )
{
    std::vector<AxisCoord> vRet;
    std::vector<AxisRegion> vRet2;
    vRet.reserve( 2 * ( uiScreenEndPos - uiScreenStartPos ) / uiBinSize );
    size_t uiCurrScreenPos = uiScreenStartPos;
    size_t uiChromosomeStartPos = 0;
    size_t uiChr = 0;
    for( size_t uiI = 0; uiI < vChromosomes.size( ); uiI++ )
    {
        if( bCancel )
            return std::make_pair( vRet, vRet2 );
        if( rJson.contains( vChromosomes[ uiI ].sName ) )
        {
            int64_t iDataSetId = rJson[ vChromosomes[ uiI ].sName ].get<int64_t>( );

            size_t uiChromSize;
            switch( iAnnoInMultipleBins )
            {
                case 0: // separate
                case 1: // stretch
                    uiChromSize = rAnno.totalIntervalSize( iDataSetId );
                    break;
                case 2: // squeeze
                    uiChromSize = rAnno.numIntervals( iDataSetId );
                    break;
                default:
                    throw std::logic_error( "unknown iAnnoInMultipleBins value" );
            }

            size_t uiChromosomeEndPos = uiChromosomeStartPos + uiChromSize;
            size_t uiItrEndPos = std::min( uiScreenEndPos, uiChromosomeEndPos );

            while( uiCurrScreenPos >= uiChromosomeStartPos && uiCurrScreenPos < uiItrEndPos )
            {
                if( bCancel )
                    return std::make_pair( vRet, vRet2 );
                auto xLower = rAnno.lowerBound( iDataSetId, uiCurrScreenPos - uiChromosomeStartPos,
                                                iAnnoInMultipleBins < 2, iAnnoInMultipleBins == 2 );
                auto xUpper = rAnno.upperBound( iDataSetId, uiCurrScreenPos - uiChromosomeStartPos + uiBinSize,
                                                iAnnoInMultipleBins < 2, iAnnoInMultipleBins == 2 );
                typename anno_t::interval_it_t xBegin, xPick;

                if( xLower >= xUpper )
                {
                    uiCurrScreenPos = uiChromosomeEndPos;
                    break;
                }

                while( xLower < xUpper )
                {
                    // uiDescId, uiAnnoStart, uiAnnoEnd
                    if( bCancel )
                        return std::make_pair( vRet, vRet2 );
                    // find bin coords
                    size_t uiIndexPos;
                    switch( iAnnoInMultipleBins )
                    {
                        case 0: // separate
                        case 1: // stretch
                            uiIndexPos = xLower->uiIntervalStart + uiCurrScreenPos - xLower->uiIntervalCoordsStart -
                                         uiChromosomeStartPos;
                            break;
                        case 2: // squeeze
                            uiIndexPos =
                                xLower->uiIntervalStart + uiCurrScreenPos - xLower->uiIntervalId - uiChromosomeStartPos;
                            break;
                        default:
                            throw std::logic_error( "unknown iAnnoInMultipleBins value" );
                    }
                    size_t uiCurrScreenSize;
                    size_t uiCurrIndexSize;
                    size_t uiAdd;
                    assert( xLower->uiIntervalEnd >= uiIndexPos );
                    assert( ( xUpper - 1 )->uiIntervalEnd >= uiIndexPos );
                    switch( iMultipleAnnosInBin )
                    {
                        case 0: // combine
                            assert( ( xUpper - 1 )->uiIntervalEnd >= uiIndexPos );
                            uiCurrIndexSize = ( xUpper - 1 )->uiIntervalEnd - uiIndexPos;
                            uiCurrScreenSize = ( xUpper - 1 )->uiIntervalCoordsEnd - xLower->uiIntervalCoordsStart;
                            break;
                        case 1: // first
                            assert( xLower->uiIntervalEnd >= uiIndexPos );
                            uiCurrIndexSize = xLower->uiIntervalEnd - uiIndexPos;
                            uiCurrScreenSize = ( xUpper - 1 )->uiIntervalCoordsEnd - xLower->uiIntervalCoordsStart;
                            break;
                        case 2: // max_fac_pow_two
                            xBegin = rAnno.begin( iDataSetId );
                            if( xLower->uiIntervalId < ( xUpper - 1 )->uiIntervalId )
                                uiAdd = getEvenlyDividableByMaxTwoPowNIn( xLower - xBegin, xUpper - xBegin );
                            else
                                uiAdd = xLower - xBegin;
                            if( uiAdd != std::numeric_limits<size_t>::max( ) )
                            {
                                xPick = xBegin + uiAdd;
                                switch( iAnnoInMultipleBins )
                                {
                                    case 0: // separate
                                        uiIndexPos =
                                            xPick->uiIntervalStart + uiCurrScreenPos - xPick->uiIntervalCoordsStart;
                                        uiCurrIndexSize = xPick->uiIntervalEnd - uiIndexPos;
                                        break;
                                    case 1: // stretch
                                    case 2: // squeeze
                                        uiIndexPos = xPick->uiIntervalStart;
                                        uiCurrIndexSize = xPick->uiIntervalEnd - xPick->uiIntervalStart;
                                        break;
                                    default:
                                        throw std::logic_error( "unknown iAnnoInMultipleBins value" );
                                }
                            }
                            else
                            {
                                uiIndexPos = 0;
                                uiCurrIndexSize = 0;
                            }

                            uiCurrScreenSize = ( xUpper - 1 )->uiIntervalCoordsEnd - xLower->uiIntervalCoordsStart;
                            break;
                        case 3: // force_separate
                            assert( xLower->uiIntervalEnd >= uiIndexPos );
                            uiCurrIndexSize = xLower->uiIntervalEnd - uiIndexPos;
                            uiCurrScreenSize = uiCurrIndexSize;
                            break;
                        default:
                            throw std::logic_error( "unknown iMultipleAnnosInBin value" );
                    }

                    switch( iAnnoInMultipleBins )
                    {
                        case 0: // separate
                        {
                            size_t uiAddGap = ( ( xUpper - 1 )->uiIntervalStart - xLower->uiIntervalEnd ) -
                                              ( ( xUpper - 1 )->uiIntervalCoordsStart - xLower->uiIntervalCoordsEnd );
                            uiCurrIndexSize = std::min( uiBinSize + uiAddGap, uiCurrIndexSize );
                            uiCurrScreenSize = std::min( uiBinSize, uiCurrScreenSize );
                        }
                        break;
                        case 1: // stretch
                            // do nothing
                            break;
                        case 2: // squeeze
                            if( iMultipleAnnosInBin == 3 )
                                uiCurrScreenSize = 1;
                            else
                                uiCurrScreenSize = xUpper - xLower;
                            break;
                        default:
                            throw std::logic_error( "unknown iAnnoInMultipleBins value" );
                    }

                    vRet.push_back( AxisCoord{
                        .uiChromosome = uiI, //
                        .uiScreenPos = uiCurrScreenPos, //
                        .uiIndexPos = uiIndexPos, //
                        .uiScreenSize = uiCurrScreenSize, //
                        .uiIndexSize = uiCurrIndexSize, //
                        .uiRegionIdx = uiChr //
                    } );
                    if( iAnnoInMultipleBins != 2 && vRet2.size( ) > 0 && vRet2.back( ).uiChromosome == uiI &&
                        vRet2.back( ).uiIndexPos + vRet2.back( ).uiIndexSize == uiIndexPos )
                    { // extend last AxisRegion
                        vRet2.back( ).uiScreenSize += uiCurrScreenSize;
                        vRet2.back( ).uiIndexSize += uiCurrIndexSize;
                        ++vRet2.back( ).uiNumCoords;
                    }
                    else // create new AxisRegion
                    {
                        vRet2.push_back( AxisRegion{
                            {
                                .uiChromosome = uiI, //
                                .uiScreenPos = uiCurrScreenPos, //
                                .uiIndexPos = uiIndexPos, //
                                .uiScreenSize = uiCurrScreenSize, //
                                .uiIndexSize = uiCurrIndexSize, //
                                .uiRegionIdx = uiChr //
                            },
                            .uiCoordStartIdx = vRet.size( ) - 1, //
                            .uiNumCoords = 1 //
                        } );
                    }

                    uiCurrScreenPos += uiCurrScreenSize;

                    // inc xLower
                    size_t uiPos;
                    switch( iMultipleAnnosInBin )
                    {
                        case 0: // combine
                        case 1: // first
                        case 2: // max_fac_pow_two
                            xLower = xUpper;
                            break;
                        case 3: // force_separate
                            uiPos = xLower->uiIntervalId;
                            while( xLower != xUpper && xLower->uiIntervalId == uiPos )
                                ++xLower;
                            break;
                        default:
                            throw std::logic_error( "unknown iMultipleAnnosInBin value" );
                    }
                }
            }
            ++uiChr;
            uiChromosomeStartPos += uiChromSize;
        }
    }

    return std::make_pair( vRet, vRet2 );
}

size_t multiple_anno( std::string sVal )
{
    if( sVal == "combine" )
        return 0;
    if( sVal == "first" )
        return 1;
    if( sVal == "max_fac_pow_two" )
        return 2;
    if( sVal == "force_separate" )
        return 3;
    throw std::logic_error( "unknown multiple_annos_in_bin value" );
}

std::vector<ChromDesc> activeChromList( json xChromLen, json xChromDisp )
{
    std::vector<ChromDesc> vRet;
    vRet.reserve( xChromDisp.size( ) );
    for( auto& xChrom : xChromDisp )
    {
        std::string sChrom = xChrom.get<std::string>( );
        vRet.emplace_back( ChromDesc{ .sName = sChrom, .uiLength = xChromLen[ sChrom ].get<size_t>( ) } );
    }
    return vRet;
}

bool PartialQuarry::setLCS( )
{
    uiLogestCommonSuffix = 0;

    if( getValue<json>( { "contigs", "list" } ).size( ) > 1 )
    {
        char cLast = '_';
        bool bContinue = true;
        while( bContinue )
        {
            ++uiLogestCommonSuffix;
            bool bFirst = true;
            for( auto& rChr : getValue<json>( { "contigs", "list" } ) )
            {
                CANCEL_RETURN;
                std::string sChr = rChr.get<std::string>( );
                if( sChr.size( ) < uiLogestCommonSuffix )
                {
                    bContinue = false;
                    break;
                }
                else if( bFirst )
                    cLast = sChr[ sChr.size( ) - uiLogestCommonSuffix ];
                else if( cLast != sChr[ sChr.size( ) - uiLogestCommonSuffix ] )
                {
                    bContinue = false;
                    break;
                }
                bFirst = false;
            }
        }
        --uiLogestCommonSuffix;
    }
    END_RETURN;
}

bool PartialQuarry::setCanvasSize( )
{
    for( size_t uiI = 0; uiI < 2; uiI++ )
    {

        size_t uiRunningStart = 0;
        std::string sCoords =
            getValue<std::string>( { "contigs", uiI == 0 ? "column_coordinates" : "row_coordinates" } );
        if( sCoords == "full_genome" )
        {
            for( ChromDesc& rDesc : this->vActiveChromosomes[ uiI ] )
            {
                CANCEL_RETURN;
                uiRunningStart += rDesc.uiLength;
            }
        }
        else
        {
            size_t iAnnoInMultipleBins =
                multiple_bins( getValue<std::string>( { "settings", "filters", "anno_in_multiple_bins" } ) );
            auto rJson = getValue<json>( { "annotation", "by_name", sCoords } );

            uiRunningStart = 0;
            for( auto xChr : this->vActiveChromosomes[ uiI ] )
            {
                CANCEL_RETURN;

                if( rJson.contains( xChr.sName ) )
                {
                    int64_t iDataSetId = rJson[ xChr.sName ].get<int64_t>( );

                    switch( iAnnoInMultipleBins )
                    {
                        case 0: // separate
                        case 1: // stretch
                            uiRunningStart += xIndices.vAnno.totalIntervalSize( iDataSetId );
                            break;
                        case 2: // squeeze
                            uiRunningStart += xIndices.vAnno.numIntervals( iDataSetId );
                    }
                }
            }
        }
        vCanvasSize[ uiI ] = uiRunningStart;
    }
    END_RETURN;
}

bool PartialQuarry::setTicks( )
{
    using namespace pybind11::literals;
    pybind11::gil_scoped_acquire acquire;
    size_t uiDividend = getValue<size_t>( { "dividend" } );
    const bool bSqueeze = getValue<std::string>( { "settings", "filters", "anno_in_multiple_bins" } ) == "squeeze";


    for( size_t uiI = 0; uiI < 2; uiI++ )
    {
        pybind11::list vNames;
        pybind11::list vStartPos;
        pybind11::list vFullList;

        size_t uiRunningStart = 0;
        std::string sCoords =
            getValue<std::string>( { "contigs", uiI == 0 ? "column_coordinates" : "row_coordinates" } );
        if( sCoords == "full_genome" )
        {
            for( ChromDesc& rDesc : this->vActiveChromosomes[ uiI ] )
            {
                CANCEL_RETURN;
                vNames.append( substringChr( rDesc.sName ) );
                vStartPos.append( uiRunningStart );
                vFullList.append( uiRunningStart );
                uiRunningStart += rDesc.uiLength;
            }
        }
        else
        {
            auto rJson = getValue<json>( { "annotation", "by_name", sCoords } );
            for( AxisRegion& xRegion : vAxisRegions[ uiI ] )
            {
                CANCEL_RETURN;
                std::string sChromName = vActiveChromosomes[ uiI ][ xRegion.uiChromosome ].sName;
                if( rJson.contains( sChromName ) )
                {
                    int64_t iDataSetId = rJson[ sChromName ].get<int64_t>( );
                    auto xFirst = xIndices.vAnno.lowerBound( iDataSetId, xRegion.uiIndexPos, !bSqueeze, bSqueeze );

                    std::vector<std::string> vSplit =
                        splitString<std::vector<std::string>>( xIndices.vAnno.desc( *xFirst ), '\n' );
                    std::string sId = "n/a";
                    const std::string csID = "ID=";
                    for( std::string sX : vSplit )
                        if( sX.substr( 0, csID.size( ) ) == csID )
                        {
                            size_t uiX = sX.rfind( ":" );
                            if( uiX != std::string::npos )
                                sId = sX.substr( csID.size( ), uiX - csID.size( ) );
                            else
                                sId = sX.substr( csID.size( ) );
                        }

                    vNames.append( substringChr( sChromName ) + " - " + sId );
                    vStartPos.append( xRegion.uiScreenPos );
                }
            }

            for( auto xChr : this->vActiveChromosomes[ uiI ] )
            {
                CANCEL_RETURN;
                if( rJson.contains( xChr.sName ) )
                {
                    int64_t iDataSetId = rJson[ xChr.sName ].get<int64_t>( );
                    vFullList.append( uiRunningStart );

                    if( bSqueeze )
                        uiRunningStart += xIndices.vAnno.numIntervals( iDataSetId );
                    else
                        uiRunningStart += xIndices.vAnno.totalIntervalSize( iDataSetId );
                }
            }
        }
        assert( uiRunningStart == vCanvasSize[ uiI ] );
        vFullList.append( uiRunningStart );

        xTicksCDS[ uiI ] = pybind11::dict( "contig_starts"_a = vStartPos,
                                           "genome_end"_a = vCanvasSize[ uiI ],
                                           "dividend"_a = uiDividend,
                                           "contig_names"_a = vNames );
        vTickLists[ uiI ] = vFullList;
        vTickLists2[ uiI ] = vStartPos;
    }
    END_RETURN;
}

const pybind11::dict PartialQuarry::getTicks( bool bXAxis )
{
    update( NodeNames::Ticks );
    return xTicksCDS[ bXAxis ? 0 : 1 ];
}

const pybind11::list PartialQuarry::getTickList( bool bXAxis )
{
    update( NodeNames::Ticks );
    return vTickLists[ bXAxis ? 0 : 1 ];
}

const pybind11::list PartialQuarry::getTickList2( bool bXAxis )
{
    update( NodeNames::Ticks );
    return vTickLists2[ bXAxis ? 0 : 1 ];
}

const std::array<size_t, 2> PartialQuarry::getCanvasSize( )
{
    update( NodeNames::CanvasSize );
    return vCanvasSize;
}

bool PartialQuarry::setActiveChrom( )
{
    for( bool bX : { true, false } )
        this->vActiveChromosomes[ bX ? 0 : 1 ] =
            activeChromList( getValue<json>( { "contigs", "lengths" } ),
                             getValue<json>( { "contigs", bX ? "displayed_on_x" : "displayed_on_y" } ) );
    END_RETURN;
}

bool PartialQuarry::setAxisCoords( )
{
    for( bool bX : { true, false } )
    {
        std::pair<std::vector<AxisCoord>, std::vector<AxisRegion>> xRet;
        if( getValue<std::string>( { "contigs", bX ? "column_coordinates" : "row_coordinates" } ) == "full_genome" )
            xRet = axisCoordsHelper(
                bX ? this->uiBinWidth : this->uiBinHeight,
                bX ? std::max( 0l, this->iStartX ) : std::max( 0l, this->iStartY ),
                bX ? std::max( 0l, this->iEndX ) : std::max( 0l, this->iEndY ),
                smaller_bin_to_num( getValue<std::string>( { "settings", "filters", "cut_off_bin" } ) ),
                this->vActiveChromosomes[ bX ? 0 : 1 ],
                this->bCancel );
        else
            xRet = annoCoordsHelper<SpsInterface<false>::anno_t>(
                bX ? this->uiBinWidth : this->uiBinHeight,
                bX ? std::max( 0l, this->iStartX ) : std::max( 0l, this->iStartY ),
                bX ? std::max( 0l, this->iEndX ) : std::max( 0l, this->iEndY ),
                smaller_bin_to_num( getValue<std::string>( { "settings", "filters", "cut_off_bin" } ) ),
                multiple_anno( getValue<std::string>( { "settings", "filters", "multiple_annos_in_bin" } ) ),
                multiple_bins( getValue<std::string>( { "settings", "filters", "anno_in_multiple_bins" } ) ),
                this->vActiveChromosomes[ bX ? 0 : 1 ], this->bCancel,
                getValue<json>(
                    { "annotation", "by_name",
                      getValue<std::string>( { "contigs", bX ? "column_coordinates" : "row_coordinates" } ) } ),
                xIndices.vAnno );
        this->vAxisCords[ bX ? 0 : 1 ] = xRet.first;
        this->vAxisRegions[ bX ? 0 : 1 ] = xRet.second;
    }
    CANCEL_RETURN;
    END_RETURN;
}

bool PartialQuarry::setGridSeqCoords( )
{
    if( getValue<std::string>( { "settings", "normalization", "normalize_by" } ) == "grid-seq" )
    {
        size_t uiGenomeSize = 0;
        for( ChromDesc& rDesc : this->vActiveChromosomes[ 0 ] )
        {
            CANCEL_RETURN;
            uiGenomeSize += rDesc.uiLength;
        }

        auto xRet = annoCoordsHelper<SpsInterface<false>::anno_t>(
            uiGenomeSize / getValue<size_t>( { "settings", "normalization", "grid_seq_slice_samples" } ), 0,
            uiGenomeSize, 0 /*<- unused */, 1 /* use first annotation in bin */,
            1 /* stretch bin over entire annotation*/, this->vActiveChromosomes[ 0 ], this->bCancel,
            getValue<json>( { "annotation", "by_name",
                              getValue<std::string>( { "settings", "normalization", "grid_seq_anno_type" } ) } ),
            xIndices.vAnno );
        vGridSeqCoords = xRet.first;
        vGridSeqRegions = xRet.second;
    }
    else
    {
        vGridSeqCoords.clear( );
        vGridSeqRegions.clear( );
    }
    CANCEL_RETURN;
    END_RETURN;
}


bool PartialQuarry::setAnnoFilters( )
{
    const json& rList = getValue<json>( { "annotation", "list" } );
    if( getValue<json>( { "annotation", "filter" } ).is_null( ) )
    {
        uiFromAnnoFilter = 0;
        uiToAnnoFilter = rList.size( ) * 3 + 2;
    }
    else
    {
        const std::array<std::array<size_t, 2>, 2> vvAF{
            /*x false*/ std::array<size_t, 2>{ /*y false*/ 0, /*y true*/ 0 },
            /*x true */ std::array<size_t, 2>{ /*y false*/ 1, /*y true*/ 1 } };
        const std::array<std::array<size_t, 2>, 2> vvAT{
            /*x false*/ std::array<size_t, 2>{ /*y false*/ 4, /*y true*/ 2 },
            /*x true */ std::array<size_t, 2>{ /*y false*/ 3, /*y true*/ 2 } };
        const size_t uiI = rList.find( getValue<std::string>( { "annotation", "filter" } ) ) - rList.begin( );
        const bool bFilterRow = getValue<bool>( { "annotation", "filter_row" } );
        const bool bFilterCol = getValue<bool>( { "annotation", "filter_col" } );
        uiFromAnnoFilter = uiI * 3 + vvAF[ bFilterRow ][ bFilterCol ];
        uiToAnnoFilter = uiI * 3 + vvAT[ bFilterRow ][ bFilterCol ];
    }

    END_RETURN;
}

bool PartialQuarry::setSymmetry( )
{
    std::string sSymmetry = getValue<std::string>( { "settings", "filters", "symmetry" } );
    if( sSymmetry == "all" )
        uiSymmetry = 0;
    else if( sSymmetry == "sym" )
        uiSymmetry = 1;
    else if( sSymmetry == "asym" )
        uiSymmetry = 2;
    else if( sSymmetry == "topToBot" )
        uiSymmetry = 3;
    else if( sSymmetry == "botToTop" )
        uiSymmetry = 4;
    else
        throw std::logic_error( "unknown symmetry setting" );
    END_RETURN;
}

bool PartialQuarry::setMappingQuality( )
{
    size_t uiMapQMin = getValue<size_t>( { "settings", "filters", "mapping_q", "val_max" } );
    size_t uiMapQMax = getValue<size_t>( { "settings", "filters", "mapping_q", "val_min" } );

    bool bIncomplAlignment = getValue<bool>( { "settings", "filters", "incomplete_alignments" } );

    if( !( uiMapQMin == 0 && bIncomplAlignment ) )
        ++uiMapQMin;
    ++uiMapQMax;

    uiMapQMin = 255 - uiMapQMin;
    uiMapQMax = 255 - uiMapQMax;
    END_RETURN;
}

template <typename out_t, typename in_t>
std::array<out_t, 2> PartialQuarry::binObjFromCoords( const in_t& xX, const in_t& xY )
{
    switch( uiSymmetry )
    {
        case 0:
            return { out_t{ BinCoordBase{ .uiChromosomeX = xX.uiChromosome,
                                          .uiChromosomeY = xY.uiChromosome,

                                          .uiScreenX = xX.uiScreenPos,
                                          .uiScreenY = xY.uiScreenPos,

                                          .uiIndexX = xX.uiIndexPos,
                                          .uiIndexY = xY.uiIndexPos,

                                          .uiScreenW = xX.uiScreenSize,
                                          .uiScreenH = xY.uiScreenSize,

                                          .uiIndexW = xX.uiIndexSize,
                                          .uiIndexH = xY.uiIndexSize } },
                     out_t{} };
            break;
        case 1:
        case 2:
            return { out_t{ BinCoordBase{ .uiChromosomeX = xX.uiChromosome,
                                          .uiChromosomeY = xY.uiChromosome,

                                          .uiScreenX = xX.uiScreenPos,
                                          .uiScreenY = xY.uiScreenPos,

                                          .uiIndexX = xX.uiIndexPos,
                                          .uiIndexY = xY.uiIndexPos,

                                          .uiScreenW = xX.uiScreenSize,
                                          .uiScreenH = xY.uiScreenSize,

                                          .uiIndexW = xX.uiIndexSize,
                                          .uiIndexH = xY.uiIndexSize } },
                     out_t{ BinCoordBase{ .uiChromosomeX = xY.uiChromosome,
                                          .uiChromosomeY = xX.uiChromosome,

                                          .uiScreenX = xX.uiScreenPos,
                                          .uiScreenY = xY.uiScreenPos,

                                          .uiIndexX = xY.uiIndexPos,
                                          .uiIndexY = xX.uiIndexPos,

                                          .uiScreenW = xY.uiScreenSize,
                                          .uiScreenH = xX.uiScreenSize,

                                          .uiIndexW = xY.uiIndexSize,
                                          .uiIndexH = xX.uiIndexSize } } };
            break;
        case 3:
            // @todo double query for bins overlapping diagonal (also fix binRegion function below)
            if( xY.uiRegionIdx > xX.uiRegionIdx ||
                ( xY.uiRegionIdx == xX.uiRegionIdx && xY.uiIndexPos + xY.uiIndexSize > xX.uiIndexPos ) )
            {
                size_t uiDiagIntersect1Y = std::max( xX.uiIndexPos, xY.uiIndexPos );
                size_t uiDiagIntersect1X = xX.uiIndexPos;

                size_t uiDiagIntersect2X = std::min( xX.uiIndexPos + xX.uiIndexSize, xY.uiIndexPos + xY.uiIndexSize );
                size_t uiDiagIntersect2Y = xY.uiIndexPos + xY.uiIndexSize;

                size_t uiIndexW = uiDiagIntersect2X - uiDiagIntersect1X;
                size_t uiIndexH = uiDiagIntersect2Y - uiDiagIntersect1Y;

                return { out_t{ BinCoordBase{ .uiChromosomeX = xY.uiChromosome,
                                              .uiChromosomeY = xX.uiChromosome,

                                              .uiScreenX = xX.uiScreenPos,
                                              .uiScreenY = xY.uiScreenPos,

                                              .uiIndexX = uiDiagIntersect1Y,
                                              .uiIndexY = uiDiagIntersect1X,

                                              .uiScreenW = xX.uiScreenSize,
                                              .uiScreenH = xY.uiScreenSize,

                                              .uiIndexW = uiIndexH,
                                              .uiIndexH = uiIndexW } },

                         out_t{} };
            }
            else
                return { out_t{ BinCoordBase{ .uiChromosomeX = xX.uiChromosome,
                                              .uiChromosomeY = xY.uiChromosome,

                                              .uiScreenX = xX.uiScreenPos,
                                              .uiScreenY = xY.uiScreenPos,

                                              .uiIndexX = xX.uiIndexPos,
                                              .uiIndexY = xY.uiIndexPos,

                                              .uiScreenW = xX.uiScreenSize,
                                              .uiScreenH = xY.uiScreenSize,

                                              .uiIndexW = xX.uiIndexSize,
                                              .uiIndexH = xY.uiIndexSize } },
                         out_t{} };
            break;
        case 4:

            // @todo double query for bins overlapping diagonal
            if( xX.uiRegionIdx > xY.uiRegionIdx ||
                ( xY.uiRegionIdx == xX.uiRegionIdx && xX.uiIndexPos /*+ xX.uiIndexSize*/ > xY.uiIndexPos ) )
                // size_t uiDiagIntersect1X = std::max( xY.uiIndexPos, xX.uiIndexPos );
                // size_t uiDiagIntersect1Y = xY.uiIndexPos;

                // size_t uiDiagIntersect2X =
                //     std::max( xX.uiIndexPos + xX.uiIndexSize, xY.uiIndexPos + xY.uiIndexSize );
                // size_t uiDiagIntersect2Y = xY.uiIndexPos + xY.uiIndexSize;

                // size_t uiIndexW = uiDiagIntersect2X - uiDiagIntersect1X;
                // size_t uiIndexH = uiDiagIntersect2Y - uiDiagIntersect1Y;

                return { out_t{ BinCoordBase{ .uiChromosomeX = xY.uiChromosome,
                                              .uiChromosomeY = xX.uiChromosome,

                                              .uiScreenX = xX.uiScreenPos,
                                              .uiScreenY = xY.uiScreenPos,

                                              .uiIndexX = xY.uiIndexPos,
                                              .uiIndexY = xX.uiIndexPos,

                                              .uiScreenW = xX.uiScreenSize,
                                              .uiScreenH = xY.uiScreenSize,

                                              .uiIndexW = xY.uiIndexSize,
                                              .uiIndexH = xX.uiIndexSize } },
                         out_t{} };
            else
                return { out_t{ BinCoordBase{ .uiChromosomeX = xX.uiChromosome,
                                              .uiChromosomeY = xY.uiChromosome,

                                              .uiScreenX = xX.uiScreenPos,
                                              .uiScreenY = xY.uiScreenPos,

                                              .uiIndexX = xX.uiIndexPos,
                                              .uiIndexY = xY.uiIndexPos,

                                              .uiScreenW = xX.uiScreenSize,
                                              .uiScreenH = xY.uiScreenSize,

                                              .uiIndexW = xX.uiIndexSize,
                                              .uiIndexH = xY.uiIndexSize } },
                         out_t{} };
            break;
        default:
            throw std::logic_error( "unknown symmetry setting" );
            break;
    }
}


std::array<BinCoordRegion, 2> PartialQuarry::binRegionObjFromCoords( const AxisRegion& xX, const AxisRegion& xY )
{
    std::array<BinCoordRegion, 2> xRet = binObjFromCoords<BinCoordRegion, AxisRegion>( xX, xY );

    switch( uiSymmetry )
    {
        case 0:
            xRet[ 0 ].uiNumCoordsX = xX.uiNumCoords;
            xRet[ 0 ].uiNumCoordsY = xY.uiNumCoords;

            xRet[ 0 ].uiCoordStartIdxX = xX.uiCoordStartIdx;
            xRet[ 0 ].uiCoordStartIdxY = xY.uiCoordStartIdx;
            break;

        case 1:
        case 2:
            xRet[ 0 ].uiNumCoordsX = xX.uiNumCoords;
            xRet[ 0 ].uiNumCoordsY = xY.uiNumCoords;

            xRet[ 0 ].uiCoordStartIdxX = xX.uiCoordStartIdx;
            xRet[ 0 ].uiCoordStartIdxY = xY.uiCoordStartIdx;


            xRet[ 1 ].uiNumCoordsX = xY.uiNumCoords;
            xRet[ 1 ].uiNumCoordsY = xX.uiNumCoords;

            xRet[ 1 ].uiCoordStartIdxX = xY.uiCoordStartIdx;
            xRet[ 1 ].uiCoordStartIdxY = xX.uiCoordStartIdx;
            break;

        case 3:
            // @todo this prob does not work
            if( xY.uiRegionIdx > xX.uiRegionIdx ||
                ( xY.uiRegionIdx == xX.uiRegionIdx && xY.uiIndexPos + xY.uiIndexSize > xX.uiIndexPos ) )
            {
                xRet[ 0 ].uiNumCoordsX = xX.uiNumCoords;
                xRet[ 0 ].uiNumCoordsY = xY.uiNumCoords;

                xRet[ 0 ].uiCoordStartIdxX = xX.uiCoordStartIdx;
                xRet[ 0 ].uiCoordStartIdxY = xY.uiCoordStartIdx;
            }
            else
            {
                xRet[ 0 ].uiNumCoordsX = xY.uiNumCoords;
                xRet[ 0 ].uiNumCoordsY = xX.uiNumCoords;

                xRet[ 0 ].uiCoordStartIdxX = xY.uiCoordStartIdx;
                xRet[ 0 ].uiCoordStartIdxY = xX.uiCoordStartIdx;
            }
            break;

        case 4:
            if( xX.uiRegionIdx > xY.uiRegionIdx ||
                ( xY.uiRegionIdx == xX.uiRegionIdx && xX.uiIndexPos /*+ xX.uiIndexSize*/ > xY.uiIndexPos ) )
            {
                xRet[ 0 ].uiNumCoordsX = xX.uiNumCoords;
                xRet[ 0 ].uiNumCoordsY = xY.uiNumCoords;

                xRet[ 0 ].uiCoordStartIdxX = xX.uiCoordStartIdx;
                xRet[ 0 ].uiCoordStartIdxY = xY.uiCoordStartIdx;
            }
            else
            {
                xRet[ 0 ].uiNumCoordsX = xY.uiNumCoords;
                xRet[ 0 ].uiNumCoordsY = xX.uiNumCoords;

                xRet[ 0 ].uiCoordStartIdxX = xY.uiCoordStartIdx;
                xRet[ 0 ].uiCoordStartIdxY = xX.uiCoordStartIdx;
            }
            break;

        default:
            throw std::logic_error( "unknown symmetry setting" );
            break;
    }

    return xRet;
}

bool PartialQuarry::setBinCoords( )
{
    size_t uiManhattenDist = 1000 * getValue<size_t>( { "settings", "filters", "min_diag_dist", "val" } ) /
                             getValue<size_t>( { "dividend" } );

    vBinCoords.clear( );
    vBinCoords.reserve( vAxisCords[ 0 ].size( ) * vAxisCords[ 1 ].size( ) );
    for( const AxisCoord& xX : vAxisCords[ 0 ] )
        for( const AxisCoord& xY : vAxisCords[ 1 ] )
        {
            CANCEL_RETURN;
            // @todo move filters to after the queries
            if( xX.uiChromosome != xY.uiChromosome ||
                (size_t)std::abs( (int64_t)xX.uiIndexPos - (int64_t)xY.uiIndexPos ) >= uiManhattenDist )
                vBinCoords.push_back( binObjFromCoords<BinCoord, AxisCoord>( xX, xY ) );
            else
                vBinCoords.push_back( { BinCoord{ }, BinCoord{} } );
        }

    vBinRegions.clear( );
    vBinRegions.reserve( vAxisRegions[ 0 ].size( ) * vAxisRegions[ 1 ].size( ) );
    for( const AxisRegion& xX : vAxisRegions[ 0 ] )
        for( const AxisRegion& xY : vAxisRegions[ 1 ] )
        {
            CANCEL_RETURN;
            vBinRegions.push_back( binRegionObjFromCoords( xX, xY ) );
        }
    END_RETURN;
}

bool PartialQuarry::setDecayCoords( )
{
    vDistDepDecCoords.clear( );
    if( getValue<bool>( { "settings", "normalization", "ddd" } ) ||
        getValue<bool>( { "settings", "normalization", "ddd_show" } ) )
    {
        vDistDepDecCoords.reserve( vBinCoords.size( ) );

        std::map<std::array<DecayCoord, 2>, size_t> vPtr;
        for( std::array<BinCoord, 2>& vCoords : vBinCoords )
        {
            CANCEL_RETURN;

            std::array<DecayCoord, 2> vKey;
            for( size_t uiI = 0; uiI < 2; uiI++ )
            {
                int64_t iA =
                    (int64_t)( vCoords[ uiI ].uiIndexX + vCoords[ uiI ].uiIndexW ) - (int64_t)vCoords[ uiI ].uiIndexY;

                int64_t iB =
                    (int64_t)vCoords[ uiI ].uiIndexX - (int64_t)( vCoords[ uiI ].uiIndexY + vCoords[ uiI ].uiIndexH );

                int64_t iS = std::min( iA, iB );
                int64_t iE = std::max( std::max( iA, iB ), iS + 1 );

                vKey[ uiI ] = DecayCoord{ vCoords[ uiI ].uiChromosomeX, vCoords[ uiI ].uiChromosomeY, iS, iE };
            }

            if( vPtr.count( vKey ) == 0 )
            {
                vPtr[ vKey ] = vDistDepDecCoords.size( );
                vDistDepDecCoords.push_back( vKey );
            }

            vCoords[ 0 ].uiDecayCoordIndex = vPtr[ vKey ];
            vCoords[ 1 ].uiDecayCoordIndex = vCoords[ 0 ].uiDecayCoordIndex;
        }
    }

    END_RETURN;
}

const std::vector<std::array<BinCoord, 2>>& PartialQuarry::getBinCoords( )
{
    update( NodeNames::BinCoords );
    return vBinCoords;
}

const pybind11::list PartialQuarry::getAnnotationList( bool bXAxis )
{
    update( NodeNames::ActiveChrom );

    pybind11::gil_scoped_acquire acquire;
    pybind11::list vRet;
    for( auto& xChr : this->vActiveChromosomes[ bXAxis ? 0 : 1 ] )
        vRet.append( xChr.sName );

    return vRet;
}


const std::vector<AxisCoord>& PartialQuarry::getAxisCoords( bool bXAxis )
{
    update( NodeNames::AxisCoords );
    return vAxisCords[ bXAxis ? 0 : 1 ];
}


void PartialQuarry::regCoords( )
{
    registerNode( NodeNames::LCS,
                  ComputeNode{ .sNodeName = "longest_common_substring",
                               .fFunc = &PartialQuarry::setLCS,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "contigs", "list" } },
                               .vSessionsIncomingInPrevious = {} } );

    registerNode( NodeNames::ActiveChrom,
                  ComputeNode{ .sNodeName = "active_chroms",
                               .fFunc = &PartialQuarry::setActiveChrom,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "contigs", "displayed_on_x" },
                                                     { "contigs", "displayed_on_y" },
                                                     { "contigs", "lengths" } },
                               .vSessionsIncomingInPrevious = {} } );

    registerNode(
        NodeNames::Ticks,
        ComputeNode{ .sNodeName = "ticks",
                     .fFunc = &PartialQuarry::setTicks,
                     .vIncomingFunctions = { NodeNames::LCS, NodeNames::AnnotationValues, NodeNames::CanvasSize },
                     .vIncomingSession = { },
                     .vSessionsIncomingInPrevious = { { "contigs", "column_coordinates" },
                                                      { "contigs", "row_coordinates" },
                                                      { "settings", "filters", "anno_in_multiple_bins" },
                                                      { "annotation", "by_name" },
                                                      { "dividend" } } } );

    registerNode( NodeNames::CanvasSize,
                  ComputeNode{ .sNodeName = "canvas_size",
                               .fFunc = &PartialQuarry::setCanvasSize,
                               .vIncomingFunctions = { NodeNames::ActiveChrom },
                               .vIncomingSession = { { "contigs", "column_coordinates" },
                                                     { "contigs", "row_coordinates" },
                                                     { "settings", "filters", "anno_in_multiple_bins" },
                                                     { "annotation", "by_name" } },
                               .vSessionsIncomingInPrevious = {} } );

    registerNode( NodeNames::AxisCoords,
                  ComputeNode{ .sNodeName = "axis_coords",
                               .fFunc = &PartialQuarry::setAxisCoords,
                               .vIncomingFunctions = { NodeNames::ActiveChrom, NodeNames::RenderArea },
                               .vIncomingSession = { { "settings", "filters", "cut_off_bin" },
                                                     { "settings", "filters", "multiple_annos_in_bin" },
                                                     { "settings", "filters", "anno_in_multiple_bins" },
                                                     { "contigs", "column_coordinates" },
                                                     { "contigs", "row_coordinates" },
                                                     { "annotation", "by_name" } },
                               .vSessionsIncomingInPrevious = {} } );

    registerNode( NodeNames::Symmetry, ComputeNode{ .sNodeName = "symmetry_setting",
                                                    .fFunc = &PartialQuarry::setSymmetry,
                                                    .vIncomingFunctions = { },
                                                    .vIncomingSession = { { "settings", "filters", "symmetry" } },
                                                    .vSessionsIncomingInPrevious = {} } );

    registerNode( NodeNames::MappingQuality,
                  ComputeNode{ .sNodeName = "mapping_quality_setting",
                               .fFunc = &PartialQuarry::setMappingQuality,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "settings", "filters", "mapping_q", "val_min" },
                                                     { "settings", "filters", "mapping_q", "val_max" },
                                                     { "settings", "filters", "incomplete_alignments" } },
                               .vSessionsIncomingInPrevious = {} } );

    registerNode( NodeNames::BinCoords,
                  ComputeNode{ .sNodeName = "bin_coords",
                               .fFunc = &PartialQuarry::setBinCoords,
                               .vIncomingFunctions = { NodeNames::AxisCoords, NodeNames::AnnoFilters,
                                                       NodeNames::IntersectionType, NodeNames::Symmetry },
                               .vIncomingSession = { { "settings", "filters", "min_diag_dist", "val" } },
                               .vSessionsIncomingInPrevious = { { "dividend" } } } );

    registerNode( NodeNames::AnnoFilters,
                  ComputeNode{ .sNodeName = "anno_filters",
                               .fFunc = &PartialQuarry::setAnnoFilters,
                               .vIncomingFunctions = { NodeNames::ActiveChrom },
                               .vIncomingSession = { { "annotation", "filterable" },
                                                     { "annotation", "filter" },
                                                     { "annotation", "filter_row" },
                                                     { "annotation", "filter_col" },
                                                     { "annotation", "by_name" },
                                                     { "contigs", "row_coordinates" },
                                                     { "contigs", "column_coordinates" },
                                                     { "annotation", "list" } },
                               .vSessionsIncomingInPrevious = {} } );

    registerNode( NodeNames::DecayCoords,
                  ComputeNode{ .sNodeName = "decay_coords",
                               .fFunc = &PartialQuarry::setDecayCoords,
                               .vIncomingFunctions = { NodeNames::BinCoords },
                               .vIncomingSession = { { "settings", "normalization", "ddd" },
                                                     { "settings", "normalization", "ddd_show" } },
                               .vSessionsIncomingInPrevious = {} } );

    registerNode( NodeNames::GridSeqCoords,
                  ComputeNode{ .sNodeName = "grid_seq_coords",
                               .fFunc = &PartialQuarry::setGridSeqCoords,
                               .vIncomingFunctions = { NodeNames::ActiveChrom },
                               .vIncomingSession = { { "settings", "normalization", "grid_seq_slice_samples" },
                                                     { "settings", "normalization", "grid_seq_anno_type" },
                                                     { "annotation", "by_name" },
                                                     { "settings", "normalization", "normalize_by" } },
                               .vSessionsIncomingInPrevious = {} } );
}

} // namespace cm