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
                case 6:
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
                case 6:
                    bAddBin = true;
                    break;
                case 1:
                    bAddBin = uiCurrBinSize == uiBinSize || uiIndexPos > 0;
                    break;
                case 3:
                case 4:
                    bAddBin = ( uiCurrScreenPos + uiCurrBinSize ) <= uiChromosomeEndPos;
                    break;
                default:
                    throw std::logic_error( "unknown iSmallerBins value" );
            }
            if( bAddBin )
                vRet.push_back( AxisCoord{
                    //{
                    /*.uiChromosome =*/uiI, //
                    /*.uiIndexPos =*/uiIndexPos, //
                    /*.uiIndexSize =*/uiCurrBinSize, //
                    //},
                    /*.uiScreenPos =*/uiCurrScreenPos, //
                    /*.uiScreenSize =*/uiCurrBinSize, //
                    /*.uiRegionIdx =*/uiChr //
                } );
            bool bIncScreenPos;
            switch( iSmallerBins )
            {
                case 0:
                case 1:
                case 2:
                case 3:
                case 5:
                case 6:
                    bIncScreenPos = true;
                    break;
                case 4:
                    bIncScreenPos = uiCurrScreenPos + uiCurrBinSize <= uiChromosomeEndPos;
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
                //{é
                //{
                /*.uiChromosome =*/uiI, //
                /*.uiIndexPos =*/uiStartChromPos, //
                /*.uiIndexSize =*/uiItrEndPos - uiStartChromPos, //
                //},
                /*.uiScreenPos =*/uiStartScreenPos, //
                /*.uiScreenSize =*/uiCurrScreenPos - uiStartScreenPos, //
                /*.uiRegionIdx =*/uiChr, //
                //},
                /*.uiCoordStartIdx =*/uiStartIdx, //
                /*.uiNumCoords =*/vRet.size( ) - uiStartIdx //
            } );
        }

        // when making contig ends larger don't skip the beginning of the next contig because of it
        if( iSmallerBins == 5 && uiCurrScreenPos >= uiChromosomeStartPos && uiCurrScreenPos > uiChromosomeEndPos )
            uiCurrScreenPos = uiChromosomeEndPos;

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
    if( sVal == "cover_multiple" )
        return 6;
    throw std::logic_error( "unknown iSmallerBins value" );
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
                        /*{*/
                        /*.uiChromosome =*/uiI, //
                        /*.uiIndexPos =*/uiIndexPos, //
                        /*.uiIndexSize =*/uiCurrIndexSize, //
                        /*},*/
                        /*.uiScreenPos =*/uiCurrScreenPos, //
                        /*.uiScreenSize =*/uiCurrScreenSize, //
                        /*.uiRegionIdx =*/uiChr //
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
                            /*{*/
                            /*{*/
                            /*.uiChromosome =*/uiI, //
                            /*.uiIndexPos =*/uiIndexPos, //
                            /*.uiIndexSize =*/uiCurrIndexSize, //
                            /*},*/
                            /*.uiScreenPos =*/uiCurrScreenPos, //
                            /*.uiScreenSize =*/uiCurrScreenSize, //
                            /*.uiRegionIdx =*/uiChr, //
                            /*},*/
                            /*.uiCoordStartIdx =*/vRet.size( ) - 1, //
                            /*.uiNumCoords =*/1 //
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
        vRet.emplace_back( ChromDesc{ /*.sName =*/sChrom, /*.uiUnadjustedLength =*/xChromLen[ sChrom ].get<size_t>( ),
                                      /*uiLength =*/0 } );
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

        const bool bAnnoCoords =
            getValue<bool>( { "settings", "filters", uiI == 0 ? "anno_coords_col" : "anno_coords_row" } );
        if( !bAnnoCoords )
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
            auto rJson = getValue<json>(
                { "annotation", "by_name", getValue<std::string>( { "contigs", "annotation_coordinates" } ) } );

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
                            uiRunningStart += pIndices->vAnno.totalIntervalSize( iDataSetId );
                            break;
                        case 2: // squeeze
                            uiRunningStart += pIndices->vAnno.numIntervals( iDataSetId );
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

        size_t uiRunningStart = 0;
        const bool bAnnoCoords =
            getValue<bool>( { "settings", "filters", uiI == 0 ? "anno_coords_col" : "anno_coords_row" } );
        auto rJson = getValue<json>(
            { "annotation", "by_name", getValue<std::string>( { "contigs", "annotation_coordinates" } ) } );
        for( ChromDesc& rDesc : this->vActiveChromosomes[ uiI ] )
        {
            CANCEL_RETURN;
            vNames.append( substringChr( rDesc.sName ) );
            vStartPos.append( uiRunningStart );
            if( !bAnnoCoords )
                uiRunningStart += rDesc.uiLength;
            else
            {
                int64_t iDataSetId = rJson[ rDesc.sName ].get<int64_t>( );
                if( rJson.contains( rDesc.sName ) )
                {
                    if( bSqueeze )
                        uiRunningStart += pIndices->vAnno.numIntervals( iDataSetId );
                    else
                        uiRunningStart += pIndices->vAnno.totalIntervalSize( iDataSetId );
                }
            }
        }

        assert( uiRunningStart == vCanvasSize[ uiI ] );
        vStartPos.append( uiRunningStart );

        xContigTicksCDS[ uiI ] = pybind11::dict( "contig_starts"_a = vStartPos, "genome_end"_a = vCanvasSize[ uiI ],
                                                 "contig_names"_a = vNames );
        vTickLists[ uiI ] = vStartPos;

        pybind11::list vScreenStart;
        pybind11::list vIndexStart;
        for( AxisRegion& xRegion : vAxisRegions[ uiI ] )
        {
            vScreenStart.append( xRegion.uiScreenPos );
            vIndexStart.append( xRegion.uiIndexPos );
        }

        xTicksCDS[ uiI ] = pybind11::dict( "screen_starts"_a = vScreenStart,
                                           "index_starts"_a = vIndexStart,
                                           "dividend"_a = uiDividend,
                                           "genome_end"_a = vCanvasSize[ uiI ] );
    }
    END_RETURN;
}

const pybind11::dict PartialQuarry::getTicks( bool bXAxis, const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::Ticks, fPyPrint );
    return xTicksCDS[ bXAxis ? 0 : 1 ];
}

const pybind11::dict PartialQuarry::getContigTicks( bool bXAxis,
                                                    const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::Ticks, fPyPrint );
    return xContigTicksCDS[ bXAxis ? 0 : 1 ];
}

const pybind11::list PartialQuarry::getTickList( bool bXAxis,
                                                 const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::Ticks, fPyPrint );
    return vTickLists[ bXAxis ? 0 : 1 ];
}


const std::array<size_t, 2> PartialQuarry::getCanvasSize( const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::CanvasSize, fPyPrint );
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

bool PartialQuarry::setActiveChromLength( )
{
    const size_t uiSmallerVal = smaller_bin_to_num( getValue<std::string>( { "settings", "filters", "cut_off_bin" } ) );
    for( bool bX : { true, false } )
    {
        const bool bGenomeCoords =
            !getValue<bool>( { "settings", "filters", bX ? "anno_coords_col" : "anno_coords_row" } );
        const size_t uiBinSize = bX ? uiBinWidth : uiBinHeight;
        for( auto& rChr : this->vActiveChromosomes[ bX ? 0 : 1 ] )
            if( uiSmallerVal == 4 && bGenomeCoords ) // fit_chrom_smaller
                rChr.uiLength = uiBinSize * ( rChr.uiUnadjustedLength / uiBinSize );
            else if( uiSmallerVal == 5 && bGenomeCoords ) // fit_chrom_larger
                rChr.uiLength = uiBinSize * ( ( ( rChr.uiUnadjustedLength - 1 ) / uiBinSize ) + 1 );
            else
                rChr.uiLength = rChr.uiUnadjustedLength;
    }
    END_RETURN;
}

bool PartialQuarry::setAxisCoords( )
{
    for( bool bX : { true, false } )
    {
        const bool bAnnoCoords =
            getValue<bool>( { "settings", "filters", bX ? "anno_coords_col" : "anno_coords_row" } );
        CANCEL_RETURN;
        std::pair<std::vector<AxisCoord>, std::vector<AxisRegion>> xRet;
        if( !bAnnoCoords )
            xRet = axisCoordsHelper(
                bX ? this->uiBinWidth : this->uiBinHeight,
                bX ? std::max( (int64_t)0, this->iStartX ) : std::max( (int64_t)0, this->iStartY ),
                bX ? std::max( (int64_t)0, this->iEndX ) : std::max( (int64_t)0, this->iEndY ),
                smaller_bin_to_num( getValue<std::string>( { "settings", "filters", "cut_off_bin" } ) ),
                this->vActiveChromosomes[ bX ? 0 : 1 ],
                this->bCancel );
        else
            xRet = annoCoordsHelper<SpsInterface<false>::anno_t>(
                bX ? this->uiBinWidth : this->uiBinHeight,
                bX ? std::max( (int64_t)0, this->iStartX ) : std::max( (int64_t)0, this->iStartY ),
                bX ? std::max( (int64_t)0, this->iEndX ) : std::max( (int64_t)0, this->iEndY ),
                smaller_bin_to_num( getValue<std::string>( { "settings", "filters", "cut_off_bin" } ) ),
                multiple_anno( getValue<std::string>( { "settings", "filters", "multiple_annos_in_bin" } ) ),
                multiple_bins( getValue<std::string>( { "settings", "filters", "anno_in_multiple_bins" } ) ),
                this->vActiveChromosomes[ bX ? 0 : 1 ], this->bCancel,
                getValue<json>(
                    { "annotation", "by_name", getValue<std::string>( { "contigs", "annotation_coordinates" } ) } ),
                pIndices->vAnno );
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
            pIndices->vAnno );
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

size_t PartialQuarry::getAnnoListId( std::string sAnno )
{
    const json& rList = getValue<json>( { "annotation", "list" } );
    size_t uiI = 0;
    for( const auto& rEle : rList )
    {
        if( rEle.get<std::string>( ) == sAnno )
            break;
        uiI += 1;
    }
    return uiI;
}

bool PartialQuarry::setAnnoFilters( )
{
    const std::string sCoordAnno = getValue<std::string>( { "contigs", "annotation_coordinates" } );
    const bool bAnnoCoordsRow = getValue<bool>( { "settings", "filters", "anno_coords_row" } );
    const bool bAnnoCoordsCol = getValue<bool>( { "settings", "filters", "anno_coords_col" } );

    const std::string sFilterAnno = getValue<std::string>( { "annotation", "filter" } );
    const bool bFilterAnnoRow = getValue<bool>( { "settings", "filters", "anno_filter_row" } );
    const bool bFilterAnnoCol = getValue<bool>( { "settings", "filters", "anno_filter_col" } );

    const bool bCoords = bAnnoCoordsRow || bAnnoCoordsCol;
    const bool bFilter = ( bFilterAnnoRow || bFilterAnnoCol ) && ( !bCoords || sCoordAnno == sFilterAnno );

    const bool bOneOf = bCoords || bFilter;
    const bool bRow = bOneOf && ( bAnnoCoordsRow || ( bFilter && bFilterAnnoRow ) );
    const bool bCol = bOneOf && ( bAnnoCoordsCol || ( bFilter && bFilterAnnoCol ) );

    if( bOneOf )
    {
        const size_t uiAnnotationFilterIdx = getAnnoListId( bCoords ? sCoordAnno : sFilterAnno );
        const std::array<std::array<size_t, 2>, 2> vvAF{
            /*x false*/ std::array<size_t, 2>{ /*y false*/ 0, /*y true*/ 0 },
            /*x true */ std::array<size_t, 2>{ /*y false*/ 1, /*y true*/ 1 } };
        const std::array<std::array<size_t, 2>, 2> vvAT{
            /*x false*/ std::array<size_t, 2>{ /*y false*/ 3, /*y true*/ 2 },
            /*x true */ std::array<size_t, 2>{ /*y false*/ 3, /*y true*/ 2 } };
        uiFromAnnoFilter = uiAnnotationFilterIdx * 3 + vvAF[ bRow ][ bCol ];
        uiToAnnoFilter = uiAnnotationFilterIdx * 3 + vvAT[ bRow ][ bCol ];
        std::cout << "uiFromAnnoFilter: " << uiFromAnnoFilter << " uiToAnnoFilter: " << uiToAnnoFilter
                  << " bCoords ? sCoordAnno : sFilterAnno:" << ( bCoords ? sCoordAnno : sFilterAnno )
                  << " uiAnnotationFilterIdx: " << uiAnnotationFilterIdx << std::endl;
    }
    else
    {
        uiFromAnnoFilter = 0;
        uiToAnnoFilter = getValue<json>( { "annotation", "list" } ).size( ) * 3 + 2;
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
    else if( sSymmetry == "mirror" )
        uiSymmetry = 3;
    else
        throw std::logic_error( "unknown symmetry setting" );
    END_RETURN;
}

bool PartialQuarry::setMappingQuality( )
{
    uiMapQMin = getValue<size_t>( { "settings", "filters", "mapping_q", "val_max" } );
    uiMapQMax = getValue<size_t>( { "settings", "filters", "mapping_q", "val_min" } );

    bool bIncomplAlignment = getValue<bool>( { "settings", "filters", "incomplete_alignments" } );

    if( !( uiMapQMax == 0 && bIncomplAlignment ) )
        ++uiMapQMax;
    ++uiMapQMin;

    uiMapQMin = 255 - uiMapQMin;
    uiMapQMax = 255 - uiMapQMax;
    END_RETURN;
}

bool PartialQuarry::setDirectionality( )
{
    const std::string sDir = getValue<std::string>( { "settings", "filters", "directionality" } );
    if( sDir == "all" )
    {
        uiFromSameStrandFilter = 0;
        uiToSameStrandFilter = 2;

        uiFromYStrandFilter = 0;
        uiToYStrandFilter = 2;

        ui1DFromStrandFilter = 0;
        ui1DToStrandFilter = 2;
    }
    else if( sDir == "same" )
    {
        uiFromSameStrandFilter = 0;
        uiToSameStrandFilter = 1;

        uiFromYStrandFilter = 0;
        uiToYStrandFilter = 2;

        ui1DFromStrandFilter = 0;
        ui1DToStrandFilter = 2;
    }
    else if( sDir == "oppo" )
    {
        uiFromSameStrandFilter = 1;
        uiToSameStrandFilter = 2;

        uiFromYStrandFilter = 0;
        uiToYStrandFilter = 2;

        ui1DFromStrandFilter = 0;
        ui1DToStrandFilter = 2;
    }
    else if( sDir == "forw" )
    {
        uiFromSameStrandFilter = 0;
        uiToSameStrandFilter = 1;

        uiFromYStrandFilter = 0;
        uiToYStrandFilter = 1;

        ui1DFromStrandFilter = 0;
        ui1DToStrandFilter = 1;
    }
    else if( sDir == "rev" )
    {
        uiFromSameStrandFilter = 0;
        uiToSameStrandFilter = 1;

        uiFromYStrandFilter = 1;
        uiToYStrandFilter = 2;

        ui1DFromStrandFilter = 1;
        ui1DToStrandFilter = 2;
    }
    else
        throw std::logic_error( "unknown directionality setting" );
    END_RETURN;
}


template <typename out_t, typename in_t> const out_t makeSymBin( const in_t& xX, const in_t& xY )
{
    return out_t{ BinCoordBase{ /*.uiChromosomeX =*/xX.uiChromosome,
                                /*.uiChromosomeY =*/xY.uiChromosome,

                                /*.uiScreenX =*/xX.uiScreenPos,
                                /*.uiScreenY =*/xY.uiScreenPos,

                                /*.uiIndexX =*/xX.uiIndexPos,
                                /*.uiIndexY =*/xY.uiIndexPos,

                                /*.uiScreenW =*/xX.uiScreenSize,
                                /*.uiScreenH =*/xY.uiScreenSize,

                                /*.uiIndexW =*/xX.uiIndexSize,
                                /*.uiIndexH =*/xY.uiIndexSize } };
}

template <typename out_t, typename in_t> const out_t makeAsymBin( const in_t& xX, const in_t& xY )
{
    return out_t{ BinCoordBase{ /*.uiChromosomeX =*/xY.uiChromosome,
                                /*.uiChromosomeY =*/xX.uiChromosome,

                                /*.uiScreenX =*/xX.uiScreenPos,
                                /*.uiScreenY =*/xY.uiScreenPos,

                                /*.uiIndexX =*/xY.uiIndexPos,
                                /*.uiIndexY =*/xX.uiIndexPos,

                                /*.uiScreenW =*/xX.uiScreenSize,
                                /*.uiScreenH =*/xY.uiScreenSize,

                                /*.uiIndexW =*/xY.uiIndexSize,
                                /*.uiIndexH =*/xX.uiIndexSize } };
}

template <typename out_t, typename in_t>
std::array<out_t, 2> PartialQuarry::binObjFromCoords( const in_t& xX, const in_t& xY )
{
    switch( uiSymmetry )
    {
        case 0: // all
            return { makeSymBin<out_t, in_t>( xX, xY ), out_t{} };
            break;
        case 1: // sym
        case 2: // asym
        case 3: // mirror
            return { makeSymBin<out_t, in_t>( xX, xY ), makeAsymBin<out_t, in_t>( xX, xY ) };
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
        case 0: // all
            xRet[ 0 ].uiNumCoordsX = xX.uiNumCoords;
            xRet[ 0 ].uiNumCoordsY = xY.uiNumCoords;

            xRet[ 0 ].uiCoordStartIdxX = xX.uiCoordStartIdx;
            xRet[ 0 ].uiCoordStartIdxY = xY.uiCoordStartIdx;
            break;

        case 1: // sym
        case 2: // asym
        case 3: // mirror
            xRet[ 0 ].uiNumCoordsX = xX.uiNumCoords;
            xRet[ 0 ].uiNumCoordsY = xY.uiNumCoords;

            xRet[ 0 ].uiCoordStartIdxX = xX.uiCoordStartIdx;
            xRet[ 0 ].uiCoordStartIdxY = xY.uiCoordStartIdx;


            xRet[ 1 ].uiNumCoordsX = xY.uiNumCoords;
            xRet[ 1 ].uiNumCoordsY = xX.uiNumCoords;

            xRet[ 1 ].uiCoordStartIdxX = xY.uiCoordStartIdx;
            xRet[ 1 ].uiCoordStartIdxY = xX.uiCoordStartIdx;
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
                if( vCoords[ uiI ].uiChromosomeX != std::numeric_limits<size_t>::max( ) )
                {
                    assert( vCoords[ uiI ].uiChromosomeX < vActiveChromosomes[ 0 ].size( ) );
                    assert( vCoords[ uiI ].uiChromosomeY < vActiveChromosomes[ 1 ].size( ) );

                    int64_t iA = (int64_t)( vCoords[ uiI ].uiIndexX + vCoords[ uiI ].uiIndexW ) -
                                 (int64_t)vCoords[ uiI ].uiIndexY;

                    int64_t iB = (int64_t)vCoords[ uiI ].uiIndexX -
                                 (int64_t)( vCoords[ uiI ].uiIndexY + vCoords[ uiI ].uiIndexH );

                    int64_t iS = std::min( iA, iB );
                    int64_t iE = std::max( std::max( iA, iB ), iS + 1 );

                    vKey[ uiI ] = DecayCoord{ vCoords[ uiI ].uiChromosomeX, vCoords[ uiI ].uiChromosomeY, iS, iE };
                }
                else
                    vKey[ uiI ] = { };
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

const std::vector<std::array<BinCoord, 2>>&
PartialQuarry::getBinCoords( const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::BinCoords, fPyPrint );
    return vBinCoords;
}

const pybind11::list PartialQuarry::getAnnotationList( bool bXAxis,
                                                       const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::ActiveChrom, fPyPrint );

    pybind11::gil_scoped_acquire acquire;
    pybind11::list vRet;
    for( auto& xChr : this->vActiveChromosomes[ bXAxis ? 0 : 1 ] )
        vRet.append( xChr.sName );

    return vRet;
}


const std::vector<AxisCoord>& PartialQuarry::getAxisCoords( bool bXAxis,
                                                            const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::AxisCoords, fPyPrint );
    return vAxisCords[ bXAxis ? 0 : 1 ];
}

size_t PartialQuarry::getAxisSize( bool bXAxis, const std::function<void( const std::string& )>& fPyPrint )
{
    update( NodeNames::AxisCoords, fPyPrint );
    return vAxisCords[ bXAxis ? 0 : 1 ].size( );
}


void PartialQuarry::regCoords( )
{
    registerNode( NodeNames::LCS,
                  ComputeNode{ /*.sNodeName =*/"longest_common_substring",
                               /*.fFunc =*/&PartialQuarry::setLCS,
                               /*.vIncomingFunctions =*/{ },
                               /*.vIncomingSession =*/{ { "contigs", "list" } },
                               /*.vSessionsIncomingInPrevious =*/{} } );

    registerNode(
        NodeNames::ActiveChrom,
        ComputeNode{ /*.sNodeName =*/"active_chroms",
                     /*.fFunc =*/&PartialQuarry::setActiveChrom,
                     /*.vIncomingFunctions =*/{ },
                     /*.vIncomingSession =*/
                     { { "contigs", "displayed_on_x" }, { "contigs", "displayed_on_y" }, { "contigs", "lengths" } },
                     /*.vSessionsIncomingInPrevious =*/{} } );

    registerNode( NodeNames::ActiveChromLength,
                  ComputeNode{ /*.sNodeName =*/"active_chroms_length",
                               /*.fFunc =*/&PartialQuarry::setActiveChromLength,
                               /*.vIncomingFunctions =*/{ NodeNames::ActiveChrom, NodeNames::BinSize },
                               /*.vIncomingSession =*/
                               { { "settings", "filters", "cut_off_bin" },
                                 { "contigs", "annotation_coordinates" },
                                 { "settings", "filters", "anno_coords_row" },
                                 { "settings", "filters", "anno_coords_col" } },
                               /*.vSessionsIncomingInPrevious =*/{} } );

    registerNode(
        NodeNames::Ticks,
        ComputeNode{ /*.sNodeName =*/"ticks",
                     /*.fFunc =*/&PartialQuarry::setTicks,
                     /*.vIncomingFunctions =*/{ NodeNames::LCS, NodeNames::AnnotationValues, NodeNames::CanvasSize },
                     /*.vIncomingSession =*/{ },
                     /*.vSessionsIncomingInPrevious =*/
                     { { "contigs", "annotation_coordinates" },
                       { "settings", "filters", "anno_coords_row" },
                       { "settings", "filters", "anno_coords_col" },
                       { "settings", "filters", "anno_in_multiple_bins" },
                       { "annotation", "by_name" },
                       { "dividend" } } } );

    registerNode( NodeNames::CanvasSize,
                  ComputeNode{ /*.sNodeName =*/"canvas_size",
                               /*.fFunc =*/&PartialQuarry::setCanvasSize,
                               /*.vIncomingFunctions =*/{ NodeNames::ActiveChromLength },
                               /*.vIncomingSession =*/
                               { { "contigs", "annotation_coordinates" },
                                 { "settings", "filters", "anno_coords_row" },
                                 { "settings", "filters", "anno_coords_col" },
                                 { "settings", "filters", "anno_in_multiple_bins" },
                                 { "annotation", "by_name" } },
                               /*.vSessionsIncomingInPrevious =*/{} } );

    registerNode( NodeNames::AxisCoords,
                  ComputeNode{ /*.sNodeName =*/"axis_coords",
                               /*.fFunc =*/&PartialQuarry::setAxisCoords,
                               /*.vIncomingFunctions =*/{ NodeNames::ActiveChromLength, NodeNames::RenderArea },
                               /*.vIncomingSession =*/
                               { { "settings", "filters", "cut_off_bin" },
                                 { "settings", "filters", "multiple_annos_in_bin" },
                                 { "settings", "filters", "anno_in_multiple_bins" },
                                 { "contigs", "annotation_coordinates" },
                                 { "settings", "filters", "anno_coords_row" },
                                 { "settings", "filters", "anno_coords_col" },
                                 { "annotation", "by_name" } },
                               /*.vSessionsIncomingInPrevious =*/{} } );

    registerNode( NodeNames::Symmetry, ComputeNode{ /*.sNodeName =*/"symmetry_setting",
                                                    /*.fFunc =*/&PartialQuarry::setSymmetry,
                                                    /*.vIncomingFunctions =*/{ },
                                                    /*.vIncomingSession =*/{ { "settings", "filters", "symmetry" } },
                                                    /*.vSessionsIncomingInPrevious =*/{} } );

    registerNode( NodeNames::MappingQuality,
                  ComputeNode{ /*.sNodeName =*/"mapping_quality_setting",
                               /*.fFunc =*/&PartialQuarry::setMappingQuality,
                               /*.vIncomingFunctions =*/{ },
                               /*.vIncomingSession =*/
                               { { "settings", "filters", "mapping_q", "val_min" },
                                 { "settings", "filters", "mapping_q", "val_max" },
                                 { "settings", "filters", "incomplete_alignments" } },
                               /*.vSessionsIncomingInPrevious =*/{} } );

    registerNode( NodeNames::Directionality,
                  ComputeNode{ /*.sNodeName =*/"directionality_setting",
                               /*.fFunc =*/&PartialQuarry::setDirectionality,
                               /*.vIncomingFunctions =*/{ },
                               /*.vIncomingSession =*/{ { "settings", "filters", "directionality" } },
                               /*.vSessionsIncomingInPrevious =*/{} } );

    registerNode( NodeNames::BinCoords,
                  ComputeNode{ /*.sNodeName =*/"bin_coords",
                               /*.fFunc =*/&PartialQuarry::setBinCoords,
                               /*.vIncomingFunctions =*/
                               { NodeNames::AxisCoords, NodeNames::AnnoFilters, NodeNames::IntersectionType,
                                 NodeNames::Symmetry },
                               /*.vIncomingSession =*/{ { "settings", "filters", "min_diag_dist", "val" } },
                               /*.vSessionsIncomingInPrevious =*/{ { "dividend" } } } );

    registerNode( NodeNames::AnnoFilters,
                  ComputeNode{ /*.sNodeName =*/"anno_filters",
                               /*.fFunc =*/&PartialQuarry::setAnnoFilters,
                               /*.vIncomingFunctions =*/{ NodeNames::ActiveChromLength },
                               /*.vIncomingSession =*/
                               { { "annotation", "filter" },
                                 { "settings", "filters", "anno_filter_row" },
                                 { "settings", "filters", "anno_filter_col" },
                                 { "contigs", "annotation_coordinates" },
                                 { "settings", "filters", "anno_coords_row" },
                                 { "settings", "filters", "anno_coords_col" },
                                 { "annotation", "list" } },
                               /*.vSessionsIncomingInPrevious =*/{} } );

    registerNode( NodeNames::DecayCoords,
                  ComputeNode{ /*.sNodeName =*/"decay_coords",
                               /*.fFunc =*/&PartialQuarry::setDecayCoords,
                               /*.vIncomingFunctions =*/{ NodeNames::BinCoords },
                               /*.vIncomingSession =*/
                               { { "settings", "normalization", "ddd" }, { "settings", "normalization", "ddd_show" } },
                               /*.vSessionsIncomingInPrevious =*/{} } );

    registerNode( NodeNames::GridSeqCoords,
                  ComputeNode{ /*.sNodeName =*/"grid_seq_coords",
                               /*.fFunc =*/&PartialQuarry::setGridSeqCoords,
                               /*.vIncomingFunctions =*/{ NodeNames::ActiveChromLength },
                               /*.vIncomingSession =*/
                               { { "settings", "normalization", "grid_seq_slice_samples" },
                                 { "settings", "normalization", "grid_seq_anno_type" },
                                 { "annotation", "by_name" },
                                 { "settings", "normalization", "normalize_by" } },
                               /*.vSessionsIncomingInPrevious =*/{} } );
}

} // namespace cm