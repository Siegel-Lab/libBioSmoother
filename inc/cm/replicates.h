#include "cm/computation.h"
#include <cmath>

#pragma once

namespace cm
{

void Computation::setActiveReplicates( )
{
    for( auto& xRepl : this->xSession[ "replicates" ][ "list" ] )
    {
        std::string sRepl = xRepl.get<std::string>( );
        if( this->xSession[ "replicates" ][ "by_name" ][ sRepl ][ "in_group_a" ].get<bool>( ) ||
            this->xSession[ "replicates" ][ "by_name" ][ sRepl ][ "in_group_b" ].get<bool>( ) )
            vActiveReplicates.push_back( sRepl );
    }
}

void Computation::setBinValues( )
{
    vvBinValues.reserve( vActiveReplicates.size( ) );

    size_t uiSymmetry = this->getSymmetry( );

    std::string sRenderSettings = this->xRenderSettings[ "filters" ][ "ambiguous_mapping" ].get<std::string>( );
    sps::IntersectionType xIntersect;
    if( sRenderSetting == "enclosed" )
        xIntersect = sps::IntersectionType::enclosed;
    else if( sRenderSetting == "encloses" )
        xIntersect = sps::IntersectionType::encloses;
    else if( sRenderSetting == "overlaps" )
        xIntersect = sps::IntersectionType::overlaps;
    else if( sRenderSetting == "first" )
        xIntersect = sps::IntersectionType::first;
    else if( sRenderSetting == "last" )
        xIntersect = sps::IntersectionType::last;
    else if( sRenderSetting == "points_only" )
        xIntersect = sps::IntersectionType::points_only;
    else
        throw std::runtime_error( "unknown ambiguous_mapping value" );

    for( std::string& sRep : vActiveReplicates )
    {
        vvBinValues.emplace_back( );
        vvBinValues.back( ).reserve( vBinCoords.size( ) );
        auto& xRep = this->xSession[ "replicates" ][ "by_name" ][ sRep ];

        bool bHasMapQ = xRep[ "has_map_q" ];
        bool bHasMultiMap = xRep[ "has_multimapping" ];

        size_t uiMapQMin = 255 - this->xRenderSettings[ "filters" ][ "mapping_q" ][ "val_min" ].get<size_t>( );
        size_t uiMapQMax = 255 - this->xRenderSettings[ "filters" ][ "mapping_q" ][ "val_max" ].get<size_t>( );

        for( std::array<BinCoord, 2>& vCoords : vBinCoords )
        {
            std::array<size_t, 2> vVals;
            for( size_t uiI = 0; uiI < 2; uiI++ )
                if( vCoords[ uiI ].sChromosomeX != "" )
                {
                    int64_t iDataSetId =
                        xRep[ "ids" ][ vCoords[ uiI ].sChromosomeX ][ vCoords[ uiI ].sChromosomeY ].get<int64_t>( );

                    if( iDataSetId != -1 )
                    {
                        if( bHasMapQ && bHasMultiMap )
                            vVals[ uiI ] = xIndices.getIndex<3, 2>( ).count(
                                iDatasetId,
                                { vCoords[ uiI ].uiIndexY, vCoords[ uiI ].uiIndexX, uiMapQMax },
                                { vCoords[ uiI ].uiIndexY + vCoords[ uiI ].uiH,
                                  vCoords[ uiI ].uiIndexX + vCoords[ uiI ].uiW, uiMapQMin },
                                xIntersect,
                                0 );
                        else if( !bHasMapQ && bHasMultiMap )
                            vVals[ uiI ] =
                                xIndices.getIndex<2, 2>( ).count( iDatasetId,
                                                                  { vCoords[ uiI ].uiIndexY, vCoords[ uiI ].uiIndexX },
                                                                  { vCoords[ uiI ].uiIndexY + vCoords[ uiI ].uiH,
                                                                    vCoords[ uiI ].uiIndexX + vCoords[ uiI ].uiW },
                                                                  xIntersect,
                                                                  0 );
                        else if( bHasMapQ && !bHasMultiMap )
                            vVals[ uiI ] = xIndices.getIndex<3, 0>( ).count(
                                iDatasetId,
                                { vCoords[ uiI ].uiIndexY, vCoords[ uiI ].uiIndexX, uiMapQMax },
                                { vCoords[ uiI ].uiIndexY + vCoords[ uiI ].uiH,
                                  vCoords[ uiI ].uiIndexX + vCoords[ uiI ].uiW, uiMapQMin },
                                xIntersect,
                                0 );
                        else // if(!bHasMapQ && !bHasMultiMap)
                            vVals[ uiI ] =
                                xIndices.getIndex<2, 0>( ).count( iDatasetId,
                                                                  { vCoords[ uiI ].uiIndexY, vCoords[ uiI ].uiIndexX },
                                                                  { vCoords[ uiI ].uiIndexY + vCoords[ uiI ].uiH,
                                                                    vCoords[ uiI ].uiIndexX + vCoords[ uiI ].uiW },
                                                                  xIntersect,
                                                                  0 );
                    }
                    else
                        vVals[ uiI ] = 0;
                }

            switch( uiSymmetry )
            {
                case 0:
                    vvBinValues.back( ).push_back( vVals[ 0 ] );
                    break;
                case 1:
                    vvBinValues.back( ).push_back( std::min( vVals[ 0 ], vVals[ 1 ] ) );
                case 2:
                    vvBinValues.back( ).push_back( (size_t)std::abs( vVals[ 0 ] - (int64_t)vVals[ 1 ] ) );
                    break;
                case 3:
                case 4:
                    vvBinValues.back( ).push_back( vVals[ 0 ] + vVals[ 1 ] );
                    break;
                default:
                    throw std::runtime_error( "unknown symmetry setting" );
                    break;
            }
        }
    }
}

/*

    def flatten_bins(self, bins):
        if len(bins) > 0:
            self.render_step_log("flatten_bins", 0, len(bins[0]))
            group_a = range(len(self.group_a))
            group_b = range(len(self.group_a), len(self.group_a) + len(self.group_b))
            ret = [[], []]
            for idx, _ in enumerate(bins[0]):
                self.render_step_log("flatten_bins", idx, len(bins[0]))
                a = []
                b = []
                for idx_2 in group_a:
                    a.append(bins[idx_2][idx])
                for idx_2 in group_b:
                    b.append(bins[idx_2][idx])
                if self.settings['replicates']['in_group'] == "min":
                    aa = min(a) if len(a) > 0 else 0
                    bb = min(b) if len(b) > 0 else 0
                elif self.settings['replicates']['in_group'] == "sum":
                    aa = sum(a)
                    bb = sum(b)
                elif self.settings['replicates']['in_group'] == "dif":
                    aa = sum(abs(x-y) for x in a for y in a)
                    bb = sum(abs(x-y) for x in b for y in b)
                else:
                    raise RuntimeError("Unknown in group value")
                ret[0].append(aa)
                ret[1].append(bb)
            return ret
        return [[], []]

*/


void Computation::setFlatValues( )
{
    if( vvBinValues.size( ) > 0 )
    {
        vvFlatValues.reserve( vvBinValues[ 0 ].size( ) );
        std::array<std::vector<size_t>, 2> vInGroup;
        vInGroup[0].reserve( vActiveReplicates.size( ) );
        vInGroup[1].reserve( vActiveReplicates.size( ) );
        for(size_t uiI = 0; uiI < vActiveReplicates.size(); uiI++)
        {
            if(this->xSession[ "replicates" ][ "by_name" ][ vActiveReplicates[uiI] ][ "in_group_a" ].get<bool>( ))
                vInGroup[0].push_back( uiI );
            if(this->xSession[ "replicates" ][ "by_name" ][ vActiveReplicates[uiI] ][ "in_group_b" ].get<bool>( ))
                vInGroup[1].push_back( uiI );
        }

        std::string sInGroupSetting = this->xRenderSettings[ "replicates" ][ "in_group" ].get<std::string>( );
        size_t iInGroupSetting;
        if(sInGroupSetting == "min")
            iInGroupSetting = 0;
        else if(sInGroupSetting == "sum")
            iInGroupSetting = 1;
        else if(sInGroupSetting == "dif")
            iInGroupSetting = 2;
        else if(sInGroupSetting == "max")
            iInGroupSetting = 3;
        else
            throw std::runtime_error("invalid value for in_group");

        for(size_t uiI = 0; uiI < vvFlatValues.size(); uiI++)
            for(size_t uiJ = 0; uiJ < 2; uiJ++)
            {
                std::vector<size_t> vCollected;
                vCollected.reserve(vInGroup[uiJ].size());
                for( size_t uiX : vInGroup[uiJ] )
                    vCollected.push_back(vvBinValues[uiX][uiI]);
                
                size_t uiVal = 0;
                if(iInGroupSetting == 0 && vCollected.size() > 0)
                    uiVal = std::numeric_limits<size_t>::max();

                for(size_t uiC : vCollected)
                    switch(iInGroupSetting)
                    {
                        case 0:
                            uiVal = std::min(uiVal, uiC);
                            break;
                        case 1:
                            uiVal += uiC;
                            break;
                        case 2:
                            for(size_t uiC2 : vCollected)
                                uiVal += (size_t)std::abs(uiC - (int64_t)uiC2);
                            break;
                        case 3:
                            uiVal = std::max(uiVal, uiC);
                            break;
                        default:
                            throw std::runtime_error("invalid value for in_group");
                    }

            }
    }
    else
        vvFlatValues.clear( );
}

} // namespace cm