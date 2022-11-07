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
    for( auto xChr : vChromosomes )
    {
        size_t uiChromosomeEndPos = uiChromosomeStartPos + xChr.uiLength;
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
                    uiCurrBinSize = std::min( uiBinSize, xChr.uiLength - uiIndexPos );
                    break;
                case 2:
                    uiCurrBinSize = uiIndexPos + 2 * uiBinSize > xChr.uiLength ? xChr.uiLength - uiIndexPos : uiBinSize;
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
                    .sChromosome = xChr.sName, //
                    .uiScreenPos = uiCurrScreenPos, //
                    .uiIndexPos = uiIndexPos, //
                    .uiSize = uiCurrBinSize //
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

        if(uiStartScreenPos >= uiChromosomeStartPos)
            vRet2.push_back( AxisRegion{
                {
                    .sChromosome = xChr.sName, //
                    .uiScreenPos = uiStartScreenPos, //
                    .uiIndexPos = uiStartChromPos, //
                    .uiSize = uiItrEndPos - uiStartChromPos, //
                },
                .uiCoordStartIdx = uiStartIdx, //
                .uiNumCoords = vRet.size( ) - uiStartIdx //
            } );

        uiChromosomeStartPos = uiChromosomeEndPos;

        if( uiScreenEndPos <= uiChromosomeStartPos )
            break;
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

std::vector<ChromDesc> activeChromList( json& xChromLen, json& xChromDisp, std::string sExistenceKey )
{
    std::vector<ChromDesc> vRet;
    vRet.reserve( xChromDisp[ sExistenceKey ].size( ) );
    for( auto& xChrom : xChromDisp[ sExistenceKey ] )
    {
        std::string sChrom = xChrom.get<std::string>( );
        vRet.emplace_back( ChromDesc{ .sName = sChrom, .uiLength = xChromLen[ sChrom ].get<size_t>( ) } );
    }
    return vRet;
}

bool PartialQuarry::setLCS( )
{
    uiLogestCommonSuffix = 0;

    if( this->xSession[ "contigs" ][ "list" ].size( ) > 1 )
    {
        char cLast = '_';
        bool bContinue = true;
        while( bContinue )
        {
            ++uiLogestCommonSuffix;
            bool bFirst = true;
            for( auto& rChr : this->xSession[ "contigs" ][ "list" ] )
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
            }
        }
        --uiLogestCommonSuffix;
    }
    END_RETURN;
}

bool PartialQuarry::setTicks( )
{
    using namespace pybind11::literals;
    pybind11::gil_scoped_acquire acquire;


    for( size_t uiI = 0; uiI < 2; uiI++ )
    {
        pybind11::list vNames;
        pybind11::list vStartPos;
        pybind11::list vFullList;

        size_t uiRunningStart = 0;
        for( ChromDesc& rDesc : this->vActiveChromosomes[ uiI ] )
        {
            CANCEL_RETURN;
            vNames.append( rDesc.sName.substr( 0, rDesc.sName.size( ) - uiLogestCommonSuffix ) );
            vStartPos.append( uiRunningStart );
            vFullList.append( uiRunningStart );
            uiRunningStart += rDesc.uiLength;
        }
        vFullList.append( uiRunningStart );
        vCanvasSize[ uiI ] = uiRunningStart;

        xTicksCDS[ uiI ] = pybind11::dict( "contig_starts"_a = vStartPos,
                                           "genome_end"_a = uiRunningStart,
                                           "dividend"_a = this->xSession[ "dividend" ].get<size_t>( ),
                                           "contig_names"_a = vNames );
        vTickLists[ uiI ] = vFullList;
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

const std::array<size_t, 2> PartialQuarry::getCanvasSize( )
{
    update( NodeNames::Ticks );
    return vCanvasSize;
}

bool PartialQuarry::setActiveChrom( )
{
    for( bool bX : { true, false } )
        this->vActiveChromosomes[ bX ? 0 : 1 ] = activeChromList( this->xSession[ "contigs" ][ "lengths" ],
                                                                  this->xSession[ "contigs" ],
                                                                  bX ? "displayed_on_x" : "displayed_on_y" );
    END_RETURN;
}

bool PartialQuarry::setAxisCoords( )
{
    for( bool bX : { true, false } )
    {
        auto xRet = axisCoordsHelper(
            bX ? this->uiBinWidth : this->uiBinHeight,
            bX ? std::max( 0l, this->iStartX ) : std::max( 0l, this->iStartY ),
            bX ? std::max( 0l, this->iEndX ) : std::max( 0l, this->iEndY ),
            smaller_bin_to_num( this->xSession[ "settings" ][ "filters" ][ "cut_off_bin" ].get<std::string>( ) ),
            this->vActiveChromosomes[ bX ? 0 : 1 ],
            this->bCancel );
        this->vAxisCords[ bX ? 0 : 1 ] = xRet.first;
        this->vAxisRegions[ bX ? 0 : 1 ] = xRet.second;
    }
    CANCEL_RETURN;
    END_RETURN;
}

bool PartialQuarry::setSymmetry( )
{
    std::string sSymmetry = this->xSession[ "settings" ][ "filters" ][ "symmetry" ].get<std::string>( );
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
bool PartialQuarry::setBinCoords( )
{
    size_t uiManhattenDist = 1000 *
                             this->xSession[ "settings" ][ "filters" ][ "min_diag_dist" ][ "val" ].get<size_t>( ) /
                             this->xSession[ "dividend" ].get<size_t>( );

    vBinCoords.clear( );
    vBinCoords.reserve( vAxisCords[ 0 ].size( ) * vAxisCords[ 1 ].size( ) );
    for( AxisCoord& xX : vAxisCords[ 0 ] )
        for( AxisCoord& xY : vAxisCords[ 1 ] )
        {
            CANCEL_RETURN;
            if( xX.sChromosome != xY.sChromosome ||
                (size_t)std::abs( (int64_t)xX.uiIndexPos - (int64_t)xY.uiIndexPos ) >= uiManhattenDist )
                switch( uiSymmetry )
                {
                    case 0:
                        vBinCoords.push_back( { BinCoord{ .sChromosomeX = xX.sChromosome,
                                                          .sChromosomeY = xY.sChromosome,

                                                          .uiScreenX = xX.uiScreenPos,
                                                          .uiScreenY = xY.uiScreenPos,

                                                          .uiIndexX = xX.uiIndexPos,
                                                          .uiIndexY = xY.uiIndexPos,

                                                          .uiW = xX.uiSize,
                                                          .uiH = xY.uiSize },
                                                BinCoord{} } );
                        break;
                    case 1:
                    case 2:
                        vBinCoords.push_back( { BinCoord{ .sChromosomeX = xX.sChromosome,
                                                          .sChromosomeY = xY.sChromosome,

                                                          .uiScreenX = xX.uiScreenPos,
                                                          .uiScreenY = xY.uiScreenPos,

                                                          .uiIndexX = xX.uiIndexPos,
                                                          .uiIndexY = xY.uiIndexPos,

                                                          .uiW = xX.uiSize,
                                                          .uiH = xY.uiSize },
                                                BinCoord{ .sChromosomeX = xY.sChromosome,
                                                          .sChromosomeY = xX.sChromosome,

                                                          .uiScreenX = xX.uiScreenPos,
                                                          .uiScreenY = xY.uiScreenPos,

                                                          .uiIndexX = xY.uiIndexPos,
                                                          .uiIndexY = xX.uiIndexPos,

                                                          .uiW = xY.uiSize,
                                                          .uiH = xX.uiSize } } );
                        break;
                    case 3:
                        if( xY.uiIndexPos + xY.uiSize > xX.uiIndexPos )
                        {
                            size_t uiDiagIntersect1Y = std::max( xX.uiIndexPos, xY.uiIndexPos );
                            size_t uiDiagIntersect1X = xX.uiIndexPos;

                            size_t uiDiagIntersect2X = std::min( xX.uiIndexPos + xX.uiSize, xY.uiIndexPos + xY.uiSize );
                            size_t uiDiagIntersect2Y = xY.uiIndexPos + xY.uiSize;

                            size_t uiW = uiDiagIntersect2X - uiDiagIntersect1X;
                            size_t uiH = uiDiagIntersect2Y - uiDiagIntersect1Y;

                            vBinCoords.push_back( { BinCoord{ .sChromosomeX = xX.sChromosome,
                                                              .sChromosomeY = xY.sChromosome,

                                                              .uiScreenX = xX.uiScreenPos,
                                                              .uiScreenY = xY.uiScreenPos,

                                                              .uiIndexX = xX.uiIndexPos,
                                                              .uiIndexY = xY.uiIndexPos,

                                                              .uiW = xX.uiSize,
                                                              .uiH = xY.uiSize },

                                                    BinCoord{ .sChromosomeX = xY.sChromosome,
                                                              .sChromosomeY = xX.sChromosome,

                                                              .uiScreenX = xX.uiScreenPos,
                                                              .uiScreenY = xY.uiScreenPos,

                                                              .uiIndexX = uiDiagIntersect1Y,
                                                              .uiIndexY = uiDiagIntersect1X,

                                                              .uiW = uiH,
                                                              .uiH = uiW } } );
                        }
                        else
                            vBinCoords.push_back( { BinCoord{ .sChromosomeX = xX.sChromosome,
                                                              .sChromosomeY = xY.sChromosome,

                                                              .uiScreenX = xX.uiScreenPos,
                                                              .uiScreenY = xY.uiScreenPos,

                                                              .uiIndexX = xX.uiIndexPos,
                                                              .uiIndexY = xY.uiIndexPos,

                                                              .uiW = xX.uiSize,
                                                              .uiH = xY.uiSize },
                                                    BinCoord{} } );
                        break;
                    case 4:
                        if( xX.uiIndexPos + xX.uiSize > xY.uiIndexPos )
                        {
                            size_t uiDiagIntersect1X = std::max( xY.uiIndexPos, xX.uiIndexPos );
                            size_t uiDiagIntersect1Y = xY.uiIndexPos;

                            size_t uiDiagIntersect2X = std::max( xX.uiIndexPos + xX.uiSize, xY.uiIndexPos + xY.uiSize );
                            size_t uiDiagIntersect2Y = xY.uiIndexPos + xY.uiSize;

                            size_t uiW = uiDiagIntersect2X - uiDiagIntersect1X;
                            size_t uiH = uiDiagIntersect2Y - uiDiagIntersect1Y;

                            vBinCoords.push_back( { BinCoord{ .sChromosomeX = xX.sChromosome,
                                                              .sChromosomeY = xY.sChromosome,

                                                              .uiScreenX = xX.uiScreenPos,
                                                              .uiScreenY = xY.uiScreenPos,

                                                              .uiIndexX = xX.uiIndexPos,
                                                              .uiIndexY = xY.uiIndexPos,

                                                              .uiW = xX.uiSize,
                                                              .uiH = xY.uiSize },
                                                    BinCoord{ .sChromosomeX = xY.sChromosome,
                                                              .sChromosomeY = xX.sChromosome,

                                                              .uiScreenX = xX.uiScreenPos,
                                                              .uiScreenY = xY.uiScreenPos,

                                                              .uiIndexX = uiDiagIntersect1Y,
                                                              .uiIndexY = uiDiagIntersect1X,

                                                              .uiW = uiH,
                                                              .uiH = uiW } } );
                        }
                        else
                            vBinCoords.push_back( { BinCoord{ .sChromosomeX = xX.sChromosome,
                                                              .sChromosomeY = xY.sChromosome,

                                                              .uiScreenX = xX.uiScreenPos,
                                                              .uiScreenY = xY.uiScreenPos,

                                                              .uiIndexX = xX.uiIndexPos,
                                                              .uiIndexY = xY.uiIndexPos,

                                                              .uiW = xX.uiSize,
                                                              .uiH = xY.uiSize },
                                                    BinCoord{} } );
                        break;
                    default:
                        throw std::logic_error( "unknown symmetry setting" );
                        break;
                }
            else
                vBinCoords.push_back( { BinCoord{ }, BinCoord{} } );
        }
    END_RETURN;
}


const std::vector<std::array<BinCoord, 2>>& PartialQuarry::getBinCoords( )
{
    update( NodeNames::BinCoords );
    return vBinCoords;
}


const std::vector<AxisCoord>& PartialQuarry::getAxisCoords( bool bXAxis )
{
    update( NodeNames::AxisCoords );
    return vAxisCords[ bXAxis ? 0 : 1 ];
}


void PartialQuarry::regCoords( )
{
    registerNode( NodeNames::LCS,
                  ComputeNode{ .sNodeName = "active_chroms",
                               .fFunc = &PartialQuarry::setLCS,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "contigs", "list" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::ActiveChrom,
                  ComputeNode{ .sNodeName = "active_chroms",
                               .fFunc = &PartialQuarry::setActiveChrom,
                               .vIncomingFunctions = { },
                               .vIncomingSession = { { "contigs", "displayed_on_x" }, { "contigs", "displayed_on_y" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::Ticks,
                  ComputeNode{ .sNodeName = "ticks",
                               .fFunc = &PartialQuarry::setTicks,
                               .vIncomingFunctions = { NodeNames::ActiveChrom, NodeNames::LCS },
                               .vIncomingSession = { },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::AxisCoords,
                  ComputeNode{ .sNodeName = "axis_coords",
                               .fFunc = &PartialQuarry::setAxisCoords,
                               .vIncomingFunctions = { NodeNames::ActiveChrom, NodeNames::RenderArea },
                               .vIncomingSession = { { "settings", "filters", "cut_off_bin" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::Symmetry, ComputeNode{ .sNodeName = "symmetry_setting",
                                                    .fFunc = &PartialQuarry::setSymmetry,
                                                    .vIncomingFunctions = { },
                                                    .vIncomingSession = { { "settings", "filters", "symmetry" } },
                                                    .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::BinCoords,
                  ComputeNode{ .sNodeName = "bin_coords",
                               .fFunc = &PartialQuarry::setBinCoords,
                               .vIncomingFunctions = { NodeNames::AxisCoords, NodeNames::Symmetry },
                               .vIncomingSession = { { "settings", "filters", "min_diag_dist", "val" } },
                               .uiLastUpdated = uiCurrTime } );
}

} // namespace cm