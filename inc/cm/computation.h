#include <nlohmann/json.hpp>

#pragma once

using json = nlohmann::json;

namespace cm
{

struct AxisCoord{
    std::string sChromosome;
    size_t uiRenderPos;
    size_t uiSize;
    size_t uiChromPos;
};

class Computation
{
  public:
    json xRenderSettings, xSession;

    size_t uiBinWidth, uiBinHeight;
    int64_t iStartX, iStartY, iEndX, iEndY;

    std::array<std::vector<>, 2> vAxisCords;

    Computation( json xRenderSettings, json xSession ) : xRenderSettings( xRenderSettings ), xSession( xSession )
    {}

  private:
    // bin_size.h
    size_t nextEvenNumber( double fX );

  public:
    // bin_size.h
    void setBinSize( );
    // bin_size.h
    void setRenderArea( );
    // coords.h
    void setAxisCoords( bool bX );

    void computeAll( )
    {
        setBinSize( );
        setRenderArea( );
        setAxisCoords(true);
        setAxisCoords(false);
    }
};

} // namespace cm
