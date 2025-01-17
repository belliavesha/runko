#include "conductor.h"

#include <cmath> 
#include <cassert>

using std::min;
using std::max;


// General dipole formula in cartesian coordinates
template<>
double fields::Conductor<3>::dipole(
        double x,
        double y,
        double z,
        int dim) {

  //non-rotating magnetic moment vector components
  double p1 = sin(chi), p2 = 0.0, p3 = cos(chi);

  // final rotating magnetic moment with phase included
  double mux = p1*cos(phase) - p2*sin(phase);
  double muy = p1*sin(phase) + p2*cos(phase);
  double muz = p3;

  double rad = std::sqrt(x*x + y*y + z*z);

  // mu . r
  double mudotr = mux*x + muy*y + muz*z;

   if     (dim == 0) return 3.0*x*mudotr/( pow(rad,5) + EPS) - mux/(pow(rad,3) + EPS); //x
   else if(dim == 1) return 3.0*y*mudotr/( pow(rad,5) + EPS) - muy/(pow(rad,3) + EPS); //y
   else if(dim == 2) return 3.0*z*mudotr/( pow(rad,5) + EPS) - muz/(pow(rad,3) + EPS); //z
   return 0.0;
}


real_short shape(real_short r, real_short r0, real_short delta) 
{ 
  return 0.5 * (1.0 - tanh((r - r0) / delta)); 
}


class StaggeredSphericalCoordinates
{
  double cx, cy, cz;
  double r;

  public:

  StaggeredSphericalCoordinates( double cenx, double ceny, double cenz, double radius) 
      : cx(cenx), cy(ceny), cz(cenz),r(radius)
  {}

  double x(double i, double stg) { return (i - cx + stg)/r; }
  double y(double j, double stg) { return (j - cy + stg)/r; }
  double z(double k, double stg) { return (k - cz + stg)/r; }
};




template<>
void fields::Conductor<3>::insert_em(
    fields::Tile<3>& tile)
{

  // Tile limits
  auto mins = tile.mins;
  //auto maxs = tile.maxs;

  auto& yee = tile.get_yee();

  real_short bxd, byd, bzd, exd, eyd, ezd;
  real_short iglob, jglob, kglob;
  real_short xr,yr,zr;
  real_short vx,vy;

  // angular velocity
  real_short Omega = 2.0*PI/period;

  if(period < EPS) Omega = 0.0; // reality check

  // helper class for staggered grid positions
  StaggeredSphericalCoordinates coord(cenx,ceny,cenz,1.0);


  // set transverse directions to zero to make this conductor
  for(int k=-3; k<static_cast<int>(tile.mesh_lengths[2])+3; k++) 
  for(int j=-3; j<static_cast<int>(tile.mesh_lengths[1])+3; j++) 
  for(int i=-3; i<static_cast<int>(tile.mesh_lengths[0])+3; i++) {

    //-----------
    // global grid coordinates
    iglob = static_cast<real_short>(i) + mins[0];
    jglob = static_cast<real_short>(j) + mins[1];
    kglob = static_cast<real_short>(k) + mins[2];


    //--------------------------------------------------
    // magnetic field
    
    // x coord staggering for B: 1,0,0
    xr = coord.x(iglob, 0.5);
    yr = coord.y(jglob, 0.0);
    zr = coord.z(kglob, 0.0);
    bxd = B0*dipole(xr,yr,zr,0);

    // y coord staggering for B: 0,1,0
    xr = coord.x(iglob, 0.0);
    yr = coord.y(jglob, 0.5);
    zr = coord.z(kglob, 0.0);
    byd = B0*dipole(xr,yr,zr,1);

    // z coord staggering for B: 0,0,1
    xr = coord.x(iglob, 0.0);
    yr = coord.y(jglob, 0.0);
    zr = coord.z(kglob, 0.5);
    bzd = B0*dipole(xr,yr,zr,2);

    yee.bx(i,j,k) = bxd;
    yee.by(i,j,k) = byd;
    yee.bz(i,j,k) = bzd;


    //--------------------------------------------------
    // electric field

    // x coord staggering for E: 0,1,1
    
    xr = coord.x(iglob, 0.0);
    yr = coord.y(jglob, 0.5);
    zr = coord.z(kglob, 0.5);
    bzd = B0*dipole(xr,yr,zr,2);
    
    vx = -Omega*yr;
    vy = +Omega*xr;
    exd = -vy*bzd;

    // y coord staggering for E: 1,0,1
    xr = coord.x(iglob, 0.5);
    yr = coord.y(jglob, 0.0);
    zr = coord.z(kglob, 0.5);
    bzd = B0*dipole(xr,yr,zr,2);

    vx = -Omega*yr;
    vy = +Omega*xr;
    eyd = vx*bzd;

    // z coord staggering for E: 1,1,0
    xr = coord.x(iglob, 0.5);
    yr = coord.y(jglob, 0.5);
    zr = coord.z(kglob, 0.0);
    bxd = B0*dipole(xr,yr,zr,0);
    byd = B0*dipole(xr,yr,zr,1);

    vx = -Omega*yr;
    vy = +Omega*xr;
    ezd = -vx*byd + vy*bxd;

    yee.ex(i,j,k) = exd;
    yee.ey(i,j,k) = eyd;
    yee.ez(i,j,k) = ezd;

    //yee.jx(i,j,k) = 0.0;
    //yee.jy(i,j,k) = 0.0;
    //yee.jz(i,j,k) = 0.0;
  }

}

template<>
void fields::Conductor<3>::update_b(
    fields::Tile<3>& tile)
{

  real_short bxd, byd, bzd;
  real_short bxi, byi, bzi;
  real_short bxnew, bynew, bznew;

  real_short iglob, jglob, kglob;

  real_short xr0,yr0,zr0;
  real_short xr,yr,zr;

  real_short s,r;

  // helper class for staggered grid positions
  StaggeredSphericalCoordinates coord(cenx,ceny,cenz,1.0);


  // smoothing scales for different components
  //real_short delta = 1.0;  //from conductor.h
  //
  //real_short delta_erad   = 1.0*delta; 
  //real_short delta_eperp  = 0.5*delta;
  //real_short delta_brad   = 1.0*delta;
  //real_short delta_bperp  = 1.0*delta;

  // Tile limits
  auto mins = tile.mins;
  //auto maxs = tile.maxs;

  auto& yee = tile.get_yee();

  // set transverse directions to zero to make this conductor
  for(int k=-1; k<static_cast<int>(tile.mesh_lengths[2])+1; k++) 
  for(int j=-1; j<static_cast<int>(tile.mesh_lengths[1])+1; j++) 
  for(int i=-1; i<static_cast<int>(tile.mesh_lengths[0])+1; i++) {

    //-----------
    // global grid coordinates
    iglob = static_cast<real_short>(i) + mins[0];
    jglob = static_cast<real_short>(j) + mins[1];
    kglob = static_cast<real_short>(k) + mins[2];

    // spherical coordinates 
    xr0 = coord.x(iglob, 0.5);
    yr0 = coord.y(jglob, 0.5);
    zr0 = coord.z(kglob, 0.5);

    // modify fields if we are inside 2R_* radius
    if(std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0) < 2.0*radius) {

      //--------------------------------------------------
      // Bx: x coord staggering for B: 1,0,0
      xr = coord.x(iglob, 0.5);
      yr = coord.y(jglob, 0.0);
      zr = coord.z(kglob, 0.0);
      r = std::sqrt( xr*xr + yr*yr + zr*zr );

      // interior diple field
      bxd = B0*dipole(xr,yr,zr,0);
      //byd = B0*dipole(xr,yr,zr,1);
      //bzd = B0*dipole(xr,yr,zr,2);

      bxi = yee.bx(i,j,k);
      //byi = yee.by(i,j,k); // interpolate to this location
      //bzi = yee.bz(i,j,k); // interpolate to this location

      s = shape(r, radius, delta);
      bxnew = s*bxd  + (1-s)*bxi;

      // Bx radial  component 
      //s = shape(r, radius-offs_brad, delta_brad);
      //bxrad =     s*(bxd - (bxd*xr + byd*yr + bzd*zr)*xr/r/r)
      //        (1-s)*(bxi - (bxi*xr + byi*yr + bzi*zr)*xr/r/r);
      //
      //s = shape(r, radius-offs_bperp, delta_bperp);
      // Bx perp  component 
      //bxperp =      s*(bxd*xr + byd*yr + bzd*zr)*xr/r/r 
      //        + (1-s)*(bxi*xr + byi*yr + bzi*zr)*xr/r/r;


      //--------------------------------------------------
      // By: y coord staggering for B: 0,1,0
      xr = coord.x(iglob, 0.0);
      yr = coord.y(jglob, 0.5);
      zr = coord.z(kglob, 0.0);
      r = std::sqrt( xr*xr + yr*yr + zr*zr );

      // interior diple field
      //bxd = B0*dipole(xr,yr,zr,0);
      byd = B0*dipole(xr,yr,zr,1);
      //bzd = B0*dipole(xr,yr,zr,2);
      byi = yee.by(i,j,k);
      s = shape(r, radius, delta);
      bynew = s*byd  + (1-s)*byi;


      //--------------------------------------------------
      // Bz: z coord staggering for B: 0,0,1
      xr = coord.x(iglob, 0.0);
      yr = coord.y(jglob, 0.0);
      zr = coord.z(kglob, 0.5);
      r = std::sqrt( xr*xr + yr*yr + zr*zr );

      // interior diple field
      //bxd = B0*dipole(xr,yr,zr,0);
      //byd = B0*dipole(xr,yr,zr,1);
      bzd = B0*dipole(xr,yr,zr,2);
      bzi = yee.bz(i,j,k);
      s = shape(r, radius, delta);
      bznew = s*bzd  + (1-s)*bzi;


      //--------------------------------------------------
      yee.bx(i,j,k) = bxnew;
      yee.by(i,j,k) = bynew;
      yee.bz(i,j,k) = bznew;
    }
  }
}


template<>
void fields::Conductor<3>::update_e(
    fields::Tile<3>& tile)
{

  real_short exd, eyd, ezd;
  real_short bxd, byd, bzd;
  real_short exi, eyi, ezi;
  real_short exnew, eynew, eznew;

  real_short iglob, jglob, kglob;

  real_short xr0,yr0,zr0;
  real_short xr,yr,zr;

  real_short vx,vy;
  real_short s,r;


  // angular velocity
  real_short Omega = 2.0*PI/period;
  if(period < EPS) Omega = 0.0; // reality check

  // helper class for staggered grid positions
  StaggeredSphericalCoordinates coord(cenx,ceny,cenz,1.0);


  // smoothing scales for Brad/Bperp components
  //real_short delta = 2.0;  // from conductor.h
  //
  //real_short delta_erad   = 1.0*delta; 
  //real_short delta_eperp  = 0.5*delta;
  //real_short delta_brad   = 1.0*delta;
  //real_short delta_bperp  = 1.0*delta;

  // Tile limits
  auto mins = tile.mins;
  //auto maxs = tile.maxs;

  auto& yee = tile.get_yee();

  // set transverse directions to zero to make this conductor
  for(int k=-1; k<static_cast<int>(tile.mesh_lengths[2])+1; k++) 
  for(int j=-1; j<static_cast<int>(tile.mesh_lengths[1])+1; j++) 
  for(int i=-1; i<static_cast<int>(tile.mesh_lengths[0])+1; i++) {

    //-----------
    // global grid coordinates
    iglob = static_cast<real_short>(i) + mins[0];
    jglob = static_cast<real_short>(j) + mins[1];
    kglob = static_cast<real_short>(k) + mins[2];

    // spherical coordinates
    xr0 = coord.x(iglob, 0.5);
    yr0 = coord.y(jglob, 0.5);
    zr0 = coord.z(kglob, 0.5);

    // modify fields if we are inside 2R_* radius
    if(std::sqrt(xr0*xr0 + yr0*yr0 + zr0*zr0) < 2.0*radius) {

      //std::cout << 
      //    "ig" << iglob 
      //    "jg" << jglob 
      //    "kg" << kglob 
      //    "xr" << xr0 
      //    "yr" << yr0 
      //    "zr" << zr0 
      //    << std::endl;

      //-------------------------------------------------- 
      // ex
      xr = coord.x(iglob, 0.0);
      yr = coord.y(jglob, 0.5);
      zr = coord.z(kglob, 0.5);
      r = std::sqrt( xr*xr + yr*yr + zr*zr );

      bzd = B0*dipole(xr,yr,zr,2);
      
      //vx = -Omega*yr;
      vy = +Omega*xr;
      exd = -vy*bzd;
      exi = yee.ex(i,j,k);

      s = shape(r, radius, delta);
      exnew = s*exd  + (1-s)*exi;


      //-------------------------------------------------- 
      // ey

      // y coord staggering for E: 1,0,1
      xr = coord.x(iglob, 0.5);
      yr = coord.y(jglob, 0.0);
      zr = coord.z(kglob, 0.5);
      r = std::sqrt( xr*xr + yr*yr + zr*zr );

      bzd = B0*dipole(xr,yr,zr,2);

      vx = -Omega*yr;
      //vy = +Omega*xr;
      eyd = vx*bzd;
      eyi = yee.ey(i,j,k);

      s = shape(r, radius, delta);
      eynew = s*eyd  + (1-s)*eyi;


      //-------------------------------------------------- 
      // ez

      // z coord staggering for E: 1,1,0
      xr = coord.x(iglob, 0.5);
      yr = coord.y(jglob, 0.5);
      zr = coord.z(kglob, 0.0);
      r = std::sqrt( xr*xr + yr*yr + zr*zr );

      bxd = B0*dipole(xr,yr,zr,0);
      byd = B0*dipole(xr,yr,zr,1);

      vx = -Omega*yr;
      vy = +Omega*xr;
      ezd = -vx*byd + vy*bxd;
      ezi = yee.ez(i,j,k);

      s = shape(r, radius, delta);
      eznew = s*ezd  + (1-s)*ezi;


      //--------------------------------------------------
      yee.ex(i,j,k) = exnew;
      yee.ey(i,j,k) = eynew;
      yee.ez(i,j,k) = eznew;
    }
  }
}

//--------------------------------------------------
// explicit template instantiation

//template class fields::Conductor<2>; // 2D
template class fields::Conductor<3>; // 3D

