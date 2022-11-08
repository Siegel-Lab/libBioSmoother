#include "sps/desc.h"
#include <algorithm>
#include <cassert>
#include <set>
#include <tuple>

#pragma once

namespace cm
{

template <template <typename> typename vec_gen_t> class AnnotationDescIndex
{
    struct Interval
    {
        size_t uiIntervalId;
        size_t uiIntervalStart;
        size_t uiIntervalEnd;
        size_t uiAnnoStart;
        size_t uiAnnoEnd;
        size_t uiDescId;
        bool bForwStrnd;
        size_t uiIntervalCoordsStart;
        size_t uiIntervalCoordsEnd;
        size_t uiAnnoCoordsStart;
        size_t uiAnnoCoordsEnd;
    };

    struct Dataset
    {
        size_t uiStart;
        size_t uiEnd;
    };

    typename sps::DescImpl<vec_gen_t> vDesc;
    vec_gen_t<Interval> xIntervalVecGen = vec_gen_t<Interval>( );
    typename vec_gen_t<Interval>::file_t xIntervalFile;
    typename vec_gen_t<Interval>::vec_t vIntervals;

    vec_gen_t<Dataset> xDatasetVecGen = vec_gen_t<Dataset>( );
    typename vec_gen_t<Dataset>::file_t xDatasetFile;
    typename vec_gen_t<Dataset>::vec_t vDatasets;

  public:
    AnnotationDescIndex( std::string sPrefix, bool bWrite )
        : vDesc( sPrefix, bWrite ),
          xIntervalFile( xIntervalVecGen.file( sPrefix + ".intervals", bWrite ) ),
          vIntervals( xIntervalVecGen.vec( xIntervalFile ) ),
          xDatasetFile( xDatasetVecGen.file( sPrefix + ".datasets", bWrite ) ),
          vDatasets( xDatasetVecGen.vec( xDatasetFile ) )
    {}

    AnnotationDescIndex( )
    {}

    size_t addIntervals( std::vector<std::tuple<size_t, size_t, std::string, bool>> vIntervalsIn,
                         size_t /*uiVerbosity = 0*/ )
    {
        size_t uiStartSize = vIntervals.size( );
        std::vector<std::pair<size_t, size_t>> vLineSweepPos;
        std::vector<size_t> vDescID;
        vDescID.resize( vIntervalsIn.size( ) );
        for( size_t uiI = 0; uiI < vIntervalsIn.size( ); uiI++ )
        {
            vLineSweepPos.push_back( std::make_pair( std::get<0>( vIntervalsIn[ uiI ] ), uiI ) );
            vLineSweepPos.push_back( std::make_pair( std::get<1>( vIntervalsIn[ uiI ] ) + 1, uiI ) );
            vDescID[ uiI ] = vDesc.add( std::get<2>( vIntervalsIn[ uiI ] ) );
        }
        std::sort( vLineSweepPos.begin( ), vLineSweepPos.end( ) );

        size_t uiIntervalCoordPos = 0;
        size_t uiIntervalId = 0;
        std::set<size_t> xActiveIntervals;
        for( size_t uiJ = 0; uiJ < vLineSweepPos.size( ); )
        {
            size_t uiCurrPos = vLineSweepPos[ uiJ ].first;
            while( uiJ < vLineSweepPos.size( ) && vLineSweepPos[ uiJ ].first == uiCurrPos )
            {
                if( std::get<0>( vIntervalsIn[ vLineSweepPos[ uiJ ].second ] ) == uiCurrPos )
                    // interval start
                    xActiveIntervals.insert( vLineSweepPos[ uiJ ].second );
                else /*if(std::get<1>(vIntervalsIn[vLineSweepPos[uiJ].second]) + 1 == uiCurrPos)*/
                {
                    assert( xActiveIntervals.count( vLineSweepPos[ uiJ ].second ) > 0 );
                    // interval end
                    xActiveIntervals.erase( vLineSweepPos[ uiJ ].second );
                }
                uiJ++;
            }

            assert( uiJ < vLineSweepPos.size( ) || xActiveIntervals.size( ) == 0 );

            if( xActiveIntervals.size( ) > 0 )
            {
                size_t uiNextPos = vLineSweepPos[ uiJ + 1 ].first;

                for( size_t uiActive : xActiveIntervals )
                    vIntervals.push_back( Interval{
                        .uiIntervalId = uiIntervalId,

                        .uiIntervalStart = uiCurrPos,
                        .uiIntervalEnd = uiNextPos,

                        .uiAnnoStart = std::get<0>( vIntervalsIn[ uiActive ] ),
                        .uiAnnoEnd = std::get<1>( vIntervalsIn[ uiActive ] ),

                        .uiDescId = vDescID[ uiActive ],
                        .bForwStrnd = std::get<3>( vIntervalsIn[ uiActive ] ),

                        .uiIntervalCoordsStart = uiIntervalCoordPos,
                        .uiIntervalCoordsEnd = uiIntervalCoordPos + ( uiNextPos - uiCurrPos ),

                        .uiAnnoCoordsStart = uiIntervalCoordPos + std::get<0>( vIntervalsIn[ uiActive ] ) - uiCurrPos,
                        .uiAnnoCoordsEnd = uiIntervalCoordPos + std::get<1>( vIntervalsIn[ uiActive ] ) - uiCurrPos,
                    } );
                uiIntervalCoordPos += ( uiNextPos - uiCurrPos );
                ++uiIntervalId;
            }
        }
        vDatasets.push_back( Dataset{ .uiStart = uiStartSize, .uiEnd = vIntervals.size( ) } );
        return vDatasets.size( ) - 1;
    }

    void iterate( size_t uiDatasetId, size_t uiFrom, size_t uiTo,
                  std::function<bool( std::tuple<size_t, size_t, std::string, bool> )> fYield,
                  bool bIntervalCoords = false )
    {
        std::set<size_t> xActiveIntervals;
        auto xEnd = vIntervals.begin( ) + vDatasets[ uiDatasetId ].uiEnd;

        auto xStart = std::lower_bound(
            vIntervals.begin( ) + vDatasets[ uiDatasetId ].uiStart,
            xEnd,
            uiFrom,
            bIntervalCoords ? []( const Interval& rI, size_t uiFrom ) { return rI.uiIntervalCoordsEnd < uiFrom; }
                            : []( const Interval& rI, size_t uiFrom ) { return rI.uiIntervalEnd < uiFrom; } );

        while( xStart < xEnd && xStart->uiIntervalStart < uiTo )
        {
            if( xActiveIntervals.count( xStart->uiDescId ) == 0 )
            {
                xActiveIntervals.insert( xStart->uiDescId );
                if( !fYield( std::make_tuple( bIntervalCoords ? xStart->uiAnnoCoordsStart : xStart->uiAnnoStart,
                                              bIntervalCoords ? xStart->uiAnnoCoordsEnd : xStart->uiAnnoEnd,
                                              vDesc.get( xStart->uiDescId ),
                                              xStart->bForwStrnd ) ) )
                    return;
            }
            ++xStart;
        }
    }

    std::vector<std::tuple<size_t, size_t, std::string, bool>> query( size_t uiDatasetId, size_t uiFrom, size_t uiTo,
                                                                      bool bIntervalCoords = false )
    {
        std::vector<std::tuple<size_t, size_t, std::string, bool>> vRet;

        iterate(
            uiDatasetId, uiFrom, uiTo,
            [ & ]( std::tuple<size_t, size_t, std::string, bool> xTup ) {
                vRet.push_back( xTup );
                return true;
            },
            bIntervalCoords );

        std::sort( vRet.begin( ), vRet.end( ) );

        return vRet;
    }

    std::array<typename vec_gen_t<Interval>::vec_t::iterator, 2> getRange( size_t uiDatasetId, size_t uiFrom,
                                                                           size_t uiTo, bool bIntervalCoords = false )
    {
        auto xBegin = vIntervals.begin( ) + vDatasets[ uiDatasetId ].uiStart;
        auto xEnd = vIntervals.begin( ) + vDatasets[ uiDatasetId ].uiEnd;

        auto xStart = std::lower_bound(
            xBegin, xEnd, uiFrom,
            bIntervalCoords ? []( const Interval& rI, size_t uiFrom ) { return rI.uiIntervalCoordsEnd < uiFrom; }
                            : []( const Interval& rI, size_t uiFrom ) { return rI.uiIntervalEnd < uiFrom; } );
        auto xStop = std::upper_bound(
            xBegin, xEnd, uiTo,
            bIntervalCoords ? []( size_t uiTo, const Interval& rI ) { return uiTo < rI.uiIntervalCoordsStart; }
                            : []( size_t uiTo, const Interval& rI ) { return uiTo < rI.uiIntervalStart; } );

        return std::array<typename vec_gen_t<Interval>::vec_t::iterator, 2>{ xStart, xStop };
    }

    size_t count( size_t uiDatasetId, size_t uiFrom, size_t uiTo, bool bIntervalCoords = false )
    {
        auto xRange = getRange( uiDatasetId, uiFrom, uiTo, bIntervalCoords );

        return xRange[ 1 ] - xRange[ 0 ];
    }

    size_t totalIntervalSize( size_t uiDatasetId )
    {
        auto xBegin = vIntervals.begin( ) + vDatasets[ uiDatasetId ].uiStart;
        auto xEnd = vIntervals.begin( ) + vDatasets[ uiDatasetId ].uiEnd - 1;


        return xEnd->uiIntervalCoordsEnd - xBegin->uiIntervalCoordsStart;
    }

    size_t numIntervals( size_t uiDatasetId )
    {
        auto xEnd = vIntervals.begin( ) + vDatasets[ uiDatasetId ].uiEnd - 1;

        return xEnd->uiIntervalId + 1;
    }
};

} // namespace cm
