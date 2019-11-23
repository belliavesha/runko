#pragma once

#include "../../pic/tile.h"
#include "../../definitions.h"

namespace pic {

/// Reflecting piston wall
template<size_t D>
class Piston
{

  public:

  Piston() {};

  double walloc = 0.0;     // x location of the wall
  double gammawall = 1.0;  // gamma of the moving wall
  double betawall  = 0.0;  // beta of the moving wall

  /// \brief interpolate electromagnetic fields to particle locations
  void solve(pic::Tile<D>& );
    
  /// \brief apply conducting boundaries behind piston head
  void field_bc(pic::Tile<D>& );

  /// Small current deposition routine for individual particles 
  void zigzag(pic::Tile<D>&, 
      real_long x0, real_long y0, real_long z0, 
      real_long x1, real_long y1, real_long z1, 
      real_long q);
};


} // end of namespace pic
