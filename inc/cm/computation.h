#include "cm/sps_interface.h"
#include <nlohmann/json.hpp>

#pragma once

using json = nlohmann::json;

namespace cm
{

struct ChromDesc
{
    std::string sName;
    size_t uiLength;
};

struct AxisCoord
{
    std::string sChromosome;
    size_t uiScreenPos;
    size_t uiIndexPos;
    size_t uiSize;
};

struct BinCoord
{
    std::string sChromosomeX, sChromosomeY;
    size_t uiScreenX, uiScreenY;
    size_t uiIndexX, uiIndexY;
    size_t uiW, uiH;
};

class ContactMapping
{
  private:
    SpsInterface<false> xIndices;

    json xRenderSettings, xSession;

    size_t uiBinWidth, uiBinHeight;
    int64_t iStartX, iStartY, iEndX, iEndY;

    std::array<std::vector<ChromDesc>, 2> vActiveChromosomes;
    std::array<std::vector<AxisCoord>, 2> vAxisCords;

    std::vector<std::string> vActiveReplicates;

    sps::IntersectionType xIntersect;

    std::vector<std::array<BinCoord, 2>> vBinCoords;

    std::array<std::vector<std::pair<std::string, bool>>, 2> vActiveCoverage;
    size_t uiSymmetry;
    size_t iInGroupSetting;

    std::vector<std::vector<size_t>> vvBinValues;
    std::vector<std::array<size_t, 2>> vvFlatValues;
    std::vector<std::array<double, 2>> vvNormalized;

    std::array<std::vector<std::vector<size_t>>, 2> vvCoverageValues;
    std::array<std::vector<size_t>, 2> vvFlatCoverageValues;

    // bin_size.h
    size_t nextEvenNumber( double fX );

    // bin_size.h
    void setBinSize( );
    // bin_size.h
    void setRenderArea( );
    // coords.h
    void setActiveChrom( );
    // coords.h
    void setAxisCoords( );

    // coords.h
    void setSymmetry( );
    // coords.h
    void setBinCoords( );

    // replicates.h
    void setIntersectionType( );

    // replicates.h
    void setActiveReplicates( );

    // coverage.h
    void setActiveCoverage( );

    // coverage.h
    void setCoverageValues( );

    // coverage.h
    void setFlatCoverageValues( );

    // replicates.h
    size_t symmetry( size_t uiA, size_t uiB );

    // replicates.h
    void setBinValues( );

    // replicates.h
    size_t getFlatValue( std::vector<size_t> vCollected );

    // replicates.h
    void setFlatValues( );

    // replicates.h
    void setInGroup( );

    // normalization.h
    void doNotNormalize( );
    // normalization.h
    void normalizeSize( size_t uiSize );
    // normalization.h
    void normalizeMaxBin( );
    // normalization.h
    void normalize( );

  public:
    ContactMapping( std::string sPrefix ) : xIndices( sPrefix )
    {}

    void computeAll( /*json xRenderSettings, json xSession*/ )
    {
        //this->xRenderSettings = xRenderSettings;
        //this->xSession = xSession;
        setBinSize( );
        setRenderArea( );

        setActiveChrom( );
        setAxisCoords( );

        setBinCoords( );

        setActiveReplicates( );
        setIntersectionType( );

        setSymmetry( );
        setInGroup( );

        setActiveCoverage( );
        setCoverageValues( );
        setFlatCoverageValues( );

        setBinValues( );
        setFlatValues( );
    }
};

} // namespace cm
