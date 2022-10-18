#include "cm/computation.h"
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
std::vector<AxisCoord> axisCoordsHelper( size_t uiBinSize, size_t uiScreenStartPos, size_t uiScreenEndPos,
                                         bool bEmptyBinForChromBorder, size_t iSmallerBins,
                                         std::vector<ChromDesc> vChromosomes )
{
    std::vector<AxisCoord> vRet;
    vRet.reserve( 2 * ( uiScreenEndPos - uiScreenStartPos ) / uiBinSize );
    size_t uiCurrScreenPos = uiScreenStartPos;
    size_t uiChromosomeStartPos = 0;
    for( auto xChr : vChromosomes )
    {
        if( uiChromosomeStartPos > 0 && bEmptyBinForChromBorder )
            vRet.push_back( AxisCoord{
                .sChromosome = "", //
                .uiScreenPos = uiCurrScreenPos, //
                .uiIndexPos = 0, //
                .uiSize = 0 //
            } );
        size_t uiChromosomeEndPos = uiChromosomeStartPos + xChr.uiLength;
        size_t uiIndexPos = 0;
        size_t uiItrEndPos = std::min( uiScreenEndPos, uiChromosomeEndPos );
        while( uiCurrScreenPos < uiItrEndPos )
        {
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

        uiChromosomeStartPos = uiChromosomeEndPos;

        if( uiScreenEndPos <= uiChromosomeStartPos )
            break;
    }

    return vRet;
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

std::vector<ChromDesc> activeChromList( json& xChromList, json& xChromLen, json& xChromDisp, std::string sExistenceKey )
{
    std::vector<ChromDesc> vRet;
    vRet.reserve( xChromList.size( ) );
    for( auto& xChrom : xChromList )
    {
        std::string sChrom = xChrom.get<std::string>( );
        if( xChromDisp[ sChrom ][ sExistenceKey ].get<bool>( ) )
            vRet.emplace_back( ChromDesc{ .sName = sChrom, .uiLength = xChromLen[ sChrom ].get<size_t>( ) } );
    }
    return vRet;
}

void ContactMapping::setActiveChrom( )
{
    for( bool bX : { true, false } )
        this->vActiveChromosomes[ bX ? 0 : 1 ] = activeChromList( this->xSession[ "contigs" ][ "list" ],
                                                                  this->xSession[ "contigs" ][ "lengths" ],
                                                                  this->xSession[ "contigs" ][ "displayed" ],
                                                                  bX ? "displayed_on_x" : "displayed_on_y" );
}

void ContactMapping::setAxisCoords( )
{
    for( bool bX : { true, false } )
        this->vAxisCords[ bX ? 0 : 1 ] = axisCoordsHelper(
            bX ? this->uiBinWidth : this->uiBinHeight,
            bX ? std::max( 0l, this->iStartX ) : std::max( 0l, this->iStartY ),
            bX ? std::max( 0l, this->iEndX ) : std::max( 0l, this->iEndY ),
            false,
            smaller_bin_to_num( this->xRenderSettings[ "filters" ][ "cut_off_bin" ].get<std::string>( ) ),
            this->vActiveChromosomes[ bX ? 0 : 1 ] );
}

void ContactMapping::setSymmetry( )
{
    std::string sSymmetry = this->xRenderSettings[ "filters" ][ "symmetry" ].get<std::string>( );
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
}
void ContactMapping::setBinCoords( )
{
    size_t uiManhattenDist = 1000 * this->xRenderSettings[ "filters" ][ "min_diag_dist" ][ "val" ].get<size_t>( ) /
                             this->xSession[ "dividend" ].get<size_t>( );

    vBinCoords.reserve( vAxisCords[ 0 ].size( ) * vAxisCords[ 1 ].size( ) );
    for( AxisCoord& xX : vAxisCords[ 0 ] )
        for( AxisCoord& xY : vAxisCords[ 1 ] )
        {
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
}

void ContactMapping::regCoords( )
{
    registerNode( NodeNames::ActiveChrom, ComputeNode{ .sNodeName = "active_chroms",
                                                       .fFunc = &ContactMapping::setActiveChrom,
                                                       .vIncomingFunctions = { },
                                                       .vIncomingSession = { { "contigs", "displayed" } },
                                                       .vIncomingRender = { },
                                                       .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::AxisCoords,
                  ComputeNode{ .sNodeName = "axis_coords",
                               .fFunc = &ContactMapping::setAxisCoords,
                               .vIncomingFunctions = { NodeNames::ActiveChrom, NodeNames::RenderArea },
                               .vIncomingSession = { },
                               .vIncomingRender = { { "filters", "cut_off_bin" } },
                               .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::Symmetry, ComputeNode{ .sNodeName = "symmetry",
                                                    .fFunc = &ContactMapping::setSymmetry,
                                                    .vIncomingFunctions = { },
                                                    .vIncomingSession = { },
                                                    .vIncomingRender = { { "filters", "symmetry" } },
                                                    .uiLastUpdated = uiCurrTime } );

    registerNode( NodeNames::BinCoords,
                  ComputeNode{ .sNodeName = "bin_coords",
                               .fFunc = &ContactMapping::setBinCoords,
                               .vIncomingFunctions = { NodeNames::AxisCoords, NodeNames::Symmetry },
                               .vIncomingSession = { },
                               .vIncomingRender = { { "filters", "min_diag_dist", "val" } },
                               .uiLastUpdated = uiCurrTime } );
}

} // namespace cm