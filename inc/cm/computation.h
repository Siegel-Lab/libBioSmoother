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
    std::vector<std::array<BinCoord, 2>> vBinCoords;

    std::vector<std::vector<size_t>> vvBinValues;
    std::vector << std::array < size_t, 2 >> vvFlatValues;

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
    size_t getSymmetry( );
    // coords.h
    void setBinCoords( );

    // replicates.h
    void setActiveReplicates( );

    // replicates.h
    void setBinValues( );

    // replicates.h
    void setFlatValues( );

  public:
    ContactMapping( std::string sPrefix ) : xIndices( sPrefix )
    {}

    void computeAll( json xRenderSettings, json xSession )
    {
        this->xRenderSettings = xRenderSettings;
        this->xSession = xSession;
        setBinSize( );
        setRenderArea( );

        setActiveChrom( );
        setAxisCoords( );

        setBinCoords( );

        setActiveReplicates( );

        setBinValues( );
        setFlatValues( );
    }
};

} // namespace cm
