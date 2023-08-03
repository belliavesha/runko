#include "rgca.h"

#include <cmath> 
#include "../../tools/signum.h"

using toolbox::sign;



inline real_long _lerp(
      real_long c000,
      real_long c100,
      real_long c010,
      real_long c110,
      real_long c001,
      real_long c101,
      real_long c011,
      real_long c111,
      real_long dx, real_long dy, real_long dz
      ) 
{
      real_long c00 = c000 * (1.0-dx) + c100 * dx;
      real_long c10 = c010 * (1.0-dx) + c110 * dx;
      real_long c0  = c00  * (1.0-dy) + c10  * dy;
      real_long c01 = c001 * (1.0-dx) + c101 * dx;
      real_long c11 = c011 * (1.0-dx) + c111 * dx;
      real_long c1  = c01  * (1.0-dy) + c11  * dy;
      real_long c   = c0   * (1.0-dz) + c1   * dz;

      
      real_long dc_dx = (c100 - c000) * (1.0 - dy) * (1.0 - dz) + (c110 - c010) * dy * (1.0 - dz) + (c101 - c001) * (1.0 - dy) * dz + (c111 - c011) * dy * dz;
      real_long dc_dy = (c010 - c000) * (1.0 - dx) * (1.0 - dz) + (c110 - c100) * dx * (1.0 - dz) + (c011 - c001) * (1.0 - dx) * dz + (c111 - c101) * dx * dz;
      real_long dc_dz = (c001 - c000) * (1.0 - dx) * (1.0 - dy) + (c101 - c100) * dx * (1.0 - dy) + (c011 - c010) * (1.0 - dx) * dy + (c111 - c110) * dx * dy;

      return {c, dc_dx, dc_dy, dc_dz};
}


//-------------------------------------------------- 
inline auto ExB_drift( 
            real_long  ex,  real_long  ey,  real_long  ez,
            real_long  bx,  real_long  by,  real_long  bz
                     ) -> std::tuple<real_long, real_long, real_long, real_long, real_long>
{
    const real_long b = sqrt( bx*bx + by*by + bz*bz );

    real_long vex = (ey*bz - ez*by)/(b*b + EPS);
    real_long vey = (ez*bx - ex*bz)/(b*b + EPS);
    real_long vez = (ex*by - ey*bx)/(b*b + EPS);

    real_long ve2 = vex*vex + vey*vey + vez*vez; //|u|^2
    real_long kappa = 1.0/(sqrt(1. - ve2) + EPS); // gamma factor

    return {vex, vey, vez, kappa, 0.};
}
    

//ExB in units of c with E.B != correction
inline auto ExB_drift_rel( 
            real_long  ex,  real_long  ey,  real_long  ez,
            real_long  bx,  real_long  by,  real_long  bz
                     ) -> std::tuple<real_long, real_long, real_long, real_long, real_long>
{
    const real_long b = sqrt( bx*bx + by*by + bz*bz );
    const real_long e = sqrt( ex*ex + ey*ey + ez*ez );

    real_long vex = (ey*bz - ez*by)/(b*b + e*e + EPS);
    real_long vey = (ez*bx - ex*bz)/(b*b + e*e + EPS);
    real_long vez = (ex*by - ey*bx)/(b*b + e*e + EPS);
    real_long we2  = vex*vex + vey*vey + vez*vez; //|u|^2

    //// we -> ve
    real_long ginv = (1. - sqrt(1. - 4.*we2))/(2.*we2 + EPS); // eq4
    vex *= ginv;
    vey *= ginv;
    vez *= ginv;

    real_long ve2 = vex*vex + vey*vey + vez*vez; //|u|^2
    // if (ve2>1.0) {
    //	    std::cout << "ve2 = "<<ve2<<" ; x :"<<vex<<" y : "<<vey<<" z : "<<vez<<"\n\n";
    //	    std::cout << std::flush;
    // }

		    
    real_long kappa = 1.0/(sqrt(1. - ve2) + EPS); // gamma factor

    return {vex, vey, vez, kappa, we2};
}


// b: normal unit vector 
inline auto mag_unit_vec( 
            real_long  bx,  real_long  by,  real_long  bz
                     ) -> std::tuple< real_long, real_long, real_long>
{
    const real_long b = sqrt( bx*bx + by*by + bz*bz );

    real_long bnx = bx/(b + EPS); 
    real_long bny = by/(b + EPS); 
    real_long bnz = bz/(b + EPS); 

    return {bnx, bny, bnz};
}


// b*: exact, E.B corrected "relativistic" unit B field vector
inline auto mag_unit_vec_rel( 
            real_long  ex,  real_long  ey,  real_long  ez,
            real_long  bx,  real_long  by,  real_long  bz,
            real_long we2
                     ) -> std::tuple< real_long, real_long, real_long>
{
    const real_long b = sqrt( bx*bx + by*by + bz*bz );
    const real_long e = sqrt( ex*ex + ey*ey + ez*ez );

    real_long edotb = ex*bx + ey*by + ez*bz;
    real_long eperpx = ex - edotb*bx/(b*b + EPS);
    real_long eperpy = ey - edotb*by/(b*b + EPS);
    real_long eperpz = ez - edotb*bz/(b*b + EPS);
    real_long eperp2 = eperpx*eperpx + eperpy*eperpy + eperpz*eperpz;

    real_long bp = sqrt(0.5*(b*b - e*e + (e*e + b*b)*sqrt(1.0 - 4.*we2))); // eq6
    real_long ep = edotb/(bp + EPS); //eq 5

    real_long psi = ep*(b*b - bp*bp)/(bp*eperp2 + EPS);
    real_long eta = 1.0/(b*sqrt( psi*psi*eperp2/(b*b + EPS) + 1.) + EPS);
    real_long zeta = psi*eta;

    // rotation giving b* 
    real_long bnx = zeta*eperpx + eta*bx;
    real_long bny = zeta*eperpy + eta*by;
    real_long bnz = zeta*eperpz + eta*bz;

    // normalize to unit vector
    real_long bn = sqrt( bnx*bnx + bny*bny + bnz*bnz );

    // bnx *= 1.0/(bn + EPS);
    // bny *= 1.0/(bn + EPS);
    // bnz *= 1.0/(bn + EPS);

    bnx = bx/(b + EPS); 
    bny = by/(b + EPS); 
    bnz = bz/(b + EPS); 

    return {bnx, bny, bnz, bn};
}



template<size_t D, size_t V>
void pic::rGCAPusher<D,V>::push_container(
    pic::ParticleContainer<D>& container, 
    pic::Tile<D>& tile
    )
{
  int nparts = container.size();

  // initialize pointers to particle arrays
  real_prtcl* loc[3];
  for( int i=0; i<3; i++) loc[i] = &( container.loc(i,0) );

  real_prtcl* vel[3];
  for( int i=0; i<3; i++) vel[i] = &( container.vel(i,0) );


  real_long ex0 = 0.0, ey0 = 0.0, ez0 = 0.0;
  real_long bx0 = 0.0, by0 = 0.0, bz0 = 0.0;

  const real_long u_upper_bound = 7e2 ; //1e14;

  // make sure E and B tmp arrays are of correct size
  if(container.Epart.size() != (size_t)3*nparts)
    container.Epart.resize(3*nparts);
  if(container.Bpart.size() != (size_t)3*nparts)
    container.Bpart.resize(3*nparts);

  // fields at prtcl loc
  real_prtcl *exP, *eyP, *ezP, *bxP, *byP, *bzP;
  exP = &( container.Epart[0*nparts] );
  eyP = &( container.Epart[1*nparts] );
  ezP = &( container.Epart[2*nparts] );

  bxP = &( container.Bpart[0*nparts] );
  byP = &( container.Bpart[1*nparts] );
  bzP = &( container.Bpart[2*nparts] );

  const int Nx = tile.mesh_lengths[0];
  const int Ny = tile.mesh_lengths[1];
  const int Nz = tile.mesh_lengths[2];

  // fields at grid
  auto& yee = tile.get_yee();

  auto& exM = yee.ex;
  auto& eyM = yee.ey;
  auto& ezM = yee.ez;

  auto& bxM = yee.bx;
  auto& byM = yee.by;
  auto& bzM = yee.bz;

  // loop over particles
  int n1 = 0;
  int n2 = nparts;

  real_long c = tile.cfl;
  real_long cinv = 1.0/c; 

  // half charge-to-mass ratio (sign only because fields are in units of q)
  real_long qm = sign(container.q)/container.m;
  real_long me = container.m;

  real_long loc0n, loc1n, loc2n;
  real_long vel0n, vel1n, vel2n;

  //real_long edotb, eperpx, eperpy, eperpz, eperp, bp, ep, psi, eta, zeta;
  real_long ugx, ugy, ugz, ug2, ug2n;
  real_long G0, G1, G01;
  real_long mu;

  // work variables
  real_long bx1, by1, bz1;
  real_long ex1, ey1, ez1;

  real_long exx1, exy1, exz1; 
  real_long eyx1, eyy1, eyz1; 
  real_long ezx1, ezy1, ezz1; 
  real_long bxx1, bxy1, bxz1; 
  real_long byx1, byy1, byz1; 
  real_long bzx1, bzy1, bzz1; 

  real_long vn1x, vn1y, vn1z;
  real_long un1x, un1y, un1z, un12;

  real_long v_cur_x, v_cur_y, v_cur_z;
  real_long v_cur_x1, v_cur_y1, v_cur_z1;
  real_long v_pol_x, v_pol_y, v_pol_z;
  real_long v_mir_x, v_mir_y, v_mir_z;

  real_long bgradb_x, bgradb_y, bgradb_z;
  real_long wegradb_x, wegradb_y, wegradb_z;

  real_long exc000, exc100, exc010, exc110, exc001, exc101, exc011, exc111;
  real_long eyc000, eyc100, eyc010, eyc110, eyc001, eyc101, eyc011, eyc111;
  real_long ezc000, ezc100, ezc010, ezc110, ezc001, ezc101, ezc011, ezc111;
  real_long bxc000, bxc100, bxc010, bxc110, bxc001, bxc101, bxc011, bxc111;
  real_long byc000, byc100, byc010, byc110, byc001, byc101, byc011, byc111;
  real_long bzc000, bzc100, bzc010, bzc110, bzc001, bzc101, bzc011, bzc111;

  real_long bnx000, bnx100, bnx010, bnx110, bnx001, bnx101, bnx011, bnx111;
  real_long bny000, bny100, bny010, bny110, bny001, bny101, bny011, bny111;
  real_long bnz000, bnz100, bnz010, bnz110, bnz001, bnz101, bnz011, bnz111;

  real_long bpr000, bpr100, bpr010, bpr110, bpr001, bpr101, bpr011, bpr111;
  
  real_long wex000, wex100, wex010, wex110, wex001, wex101, wex011, wex111;
  real_long wey000, wey100, wey010, wey110, wey001, wey101, wey011, wey111;
  real_long wez000, wez100, wez010, wez110, wez001, wez101, wez011, wez111;


  real_long wex1, wexx1, wexy1, wexz1; 
  real_long wey1, weyx1, weyy1, weyz1; 
  real_long wez1, wezx1, wezy1, wezz1; 
  real_long bnx1, bnxx1, bnxy1, bnxz1; 
  real_long bny1, bnyx1, bnyy1, bnyz1; 
  real_long bnz1, bnzx1, bnzy1, bnzz1; 

  real_long bpr1, bprx1, bpry1, bprz1; 

  // mesh sizes for 1D indexing
  const size_t iy = D >= 2 ? yee.ex.indx(0,1,0) - yee.ex.indx(0,0,0) : 0;
  const size_t iz = D >= 3 ? yee.ex.indx(0,0,1) - yee.ex.indx(0,0,0) : 0;
  auto mins = tile.mins;

  real_long dx=0.0, dy=0.0, dz=0.0;
  int i=0, j=0, k=0;

  // loop over prtcls
  for(int n=n1; n<n2; n++) {
      bool crash_flag = false;

    loc0n = static_cast<real_long>( loc[0][n] );
    loc1n = static_cast<real_long>( loc[1][n] );
    loc2n = static_cast<real_long>( loc[2][n] );

    vel0n = static_cast<real_long>( vel[0][n] );
    vel1n = static_cast<real_long>( vel[1][n] );
    vel2n = static_cast<real_long>( vel[2][n] );

    // read particle-specific fields
    ex0 = static_cast<real_long>( (exP[n] + this->get_ex_ext(0,0,0))*cinv );
    ey0 = static_cast<real_long>( (eyP[n] + this->get_ey_ext(0,0,0))*cinv );
    ez0 = static_cast<real_long>( (ezP[n] + this->get_ez_ext(0,0,0))*cinv );
    bx0 = static_cast<real_long>( (bxP[n] + this->get_bx_ext(0,0,0))*cinv );
    by0 = static_cast<real_long>( (byP[n] + this->get_by_ext(0,0,0))*cinv );
    bz0 = static_cast<real_long>( (bzP[n] + this->get_bz_ext(0,0,0))*cinv );

    //-------------------------------------------------- 
    // iterate: step0
      
    const real_long R0x = loc0n;
    const real_long R0y = loc1n;
    const real_long R0z = loc2n;

    const real_long b0 = sqrt( bx0*bx0 + by0*by0 + bz0*bz0 );
    const real_long e0 = sqrt( ex0*ex0 + ey0*ey0 + ez0*ez0 );

    {

        if(D >= 1) i  = static_cast<int>(floor(R0x));
        if(D >= 2) j  = static_cast<int>(floor(R0y));
        if(D >= 3) k  = static_cast<int>(floor(R0z));

        if(D >= 1) dx = R0x - i;
        if(D >= 2) dy = R0y - j;
        if(D >= 3) dz = R0z - k;

        // normalize to tile units
        if(D >= 1) i -= mins[0];
        if(D >= 2) j -= mins[1];
        if(D >= 3) k -= mins[2];

        //extrapolate if on tile boundary
        // min left side of tile
        if(D >= 1) { if(i <= -2 ) { dx -= i+2; i = -2; } }
        if(D >= 2) { if(j <= -2 ) { dy -= j+2; j = -2; } }
        if(D >= 3) { if(k <= -2 ) { dz -= k+2; k = -2; } }
        
        // max right side of tile
        if(D >= 1) { if(i >= Nx+1 ) { dx += i-Nx-1; i = Nx+1; } }
        if(D >= 2) { if(j >= Ny+1 ) { dy += j-Ny-1; j = Ny+1; } }
        if(D >= 3) { if(k >= Nz+1 ) { dz += k-Nz-1; k = Nz+1; } }


        const size_t ind = yee.ex.indx(i,j,k);

        //ex =
        exc000 = cinv* 0.5*(exM(ind       ) +exM(ind-1      ));
        exc100 = cinv* 0.5*(exM(ind       ) +exM(ind+1      ));
        exc010 = cinv* 0.5*(exM(ind+iy    ) +exM(ind-1+iy   ));
        exc110 = cinv* 0.5*(exM(ind+iy    ) +exM(ind+1+iy   ));
        exc001 = cinv* 0.5*(exM(ind+iz    ) +exM(ind-1+iz   ));
        exc101 = cinv* 0.5*(exM(ind+iz    ) +exM(ind+1+iz   ));
        exc011 = cinv* 0.5*(exM(ind+iy+iz ) +exM(ind-1+iy+iz));
        exc111 = cinv* 0.5*(exM(ind+iy+iz ) +exM(ind+1+iy+iz));

        //ey
        eyc000 = cinv* 0.5*(eyM(ind      ) +eyM(ind-iy     ));
        eyc100 = cinv* 0.5*(eyM(ind+1    ) +eyM(ind+1-iy   ));
        eyc010 = cinv* 0.5*(eyM(ind      ) +eyM(ind+iy     ));
        eyc110 = cinv* 0.5*(eyM(ind+1    ) +eyM(ind+1+iy   ));
        eyc001 = cinv* 0.5*(eyM(ind+iz   ) +eyM(ind-iy+iz  ));
        eyc101 = cinv* 0.5*(eyM(ind+1+iz ) +eyM(ind+1-iy+iz));
        eyc011 = cinv* 0.5*(eyM(ind+iz   ) +eyM(ind+iy+iz  ));
        eyc111 = cinv* 0.5*(eyM(ind+1+iz ) +eyM(ind+1+iy+iz));

        //ez
        ezc000 = cinv* 0.5*(ezM(ind      ) + ezM(ind-iz     ));
        ezc100 = cinv* 0.5*(ezM(ind+1    ) + ezM(ind+1-iz   ));
        ezc010 = cinv* 0.5*(ezM(ind+iy   ) + ezM(ind+iy-iz  ));
        ezc110 = cinv* 0.5*(ezM(ind+1+iy ) + ezM(ind+1+iy-iz));
        ezc001 = cinv* 0.5*(ezM(ind      ) + ezM(ind+iz     ));
        ezc101 = cinv* 0.5*(ezM(ind+1    ) + ezM(ind+1+iz   ));
        ezc011 = cinv* 0.5*(ezM(ind+iy   ) + ezM(ind+iy+iz  ));
        ezc111 = cinv* 0.5*(ezM(ind+1+iy ) + ezM(ind+1+iy+iz));

        //-------------------------------------------------- 
        // bx
        bxc000 = cinv* 0.25*( bxM(ind)+   bxM(ind-iy)+   bxM(ind-iz)+      bxM(ind-iy-iz));
        bxc100 = cinv* 0.25*( bxM(ind+1)+ bxM(ind+1-iy)+ bxM(ind+1-iz)+    bxM(ind+1-iy-iz));
        bxc001 = cinv* 0.25*( bxM(ind)+   bxM(ind+iz)+   bxM(ind-iy)+      bxM(ind-iy+iz));
        bxc101 = cinv* 0.25*( bxM(ind+1)+ bxM(ind+1+iz)+ bxM(ind+1-iy)+    bxM(ind+1-iy+iz));
        bxc010 = cinv* 0.25*( bxM(ind)+   bxM(ind+iy)+   bxM(ind-iz)+      bxM(ind+iy-iz));
        bxc110 = cinv* 0.25*( bxM(ind+1)+ bxM(ind+1-iz)+ bxM(ind+1+iy-iz)+ bxM(ind+1+iy));
        bxc011 = cinv* 0.25*( bxM(ind)+   bxM(ind+iy)+   bxM(ind+iy+iz)+   bxM(ind+iz));
        bxc111 = cinv* 0.25*( bxM(ind+1)+ bxM(ind+1+iy)+ bxM(ind+1+iy+iz)+ bxM(ind+1+iz));

        // by
        byc000 = cinv* 0.25*( byM(ind-1-iz)+    byM(ind-1)+       byM(ind-iz)+      byM(ind));
        byc100 = cinv* 0.25*( byM(ind-iz)+      byM(ind)+         byM(ind+1-iz)+    byM(ind+1));
        byc001 = cinv* 0.25*( byM(ind-1)+       byM(ind-1+iz)+    byM(ind)+         byM(ind+iz));
        byc101 = cinv* 0.25*( byM(ind)+         byM(ind+iz)+      byM(ind+1)+       byM(ind+1+iz));
        byc010 = cinv* 0.25*( byM(ind-1+iy-iz)+ byM(ind-1+iy)+    byM(ind+iy-iz)+   byM(ind+iy));
        byc110 = cinv* 0.25*( byM(ind+iy-iz)+   byM(ind+iy)+      byM(ind+1+iy-iz)+ byM(ind+1+iy));
        byc011 = cinv* 0.25*( byM(ind-1+iy)+    byM(ind-1+iy+iz)+ byM(ind+iy)+      byM(ind+iy+iz));
        byc111 = cinv* 0.25*( byM(ind+iy)+      byM(ind+iy+iz)+   byM(ind+1+iy)+    byM(ind+1+iy+iz));

        // bz
        bzc000 = cinv* 0.25*( bzM(ind-1-iy)+    bzM(ind-1)+       bzM(ind-iy)+      bzM(ind));
        bzc100 = cinv* 0.25*( bzM(ind-iy)+      bzM(ind)+         bzM(ind+1-iy)+    bzM(ind+1));
        bzc001 = cinv* 0.25*( bzM(ind-1-iy+iz)+ bzM(ind-1+iz)+    bzM(ind-iy+iz)+   bzM(ind+iz));
        bzc101 = cinv* 0.25*( bzM(ind-iy+iz)+   bzM(ind+iz)+      bzM(ind+1-iy+iz)+ bzM(ind+1+iz));
        bzc010 = cinv* 0.25*( bzM(ind-1)+       bzM(ind-1+iy)+    bzM(ind)+         bzM(ind+iy));
        bzc110 = cinv* 0.25*( bzM(ind)+         bzM(ind+iy)+      bzM(ind+1)+       bzM(ind+1+iy));
        bzc011 = cinv* 0.25*( bzM(ind-1+iz)+    bzM(ind-1+iy+iz)+ bzM(ind+iz)+      bzM(ind+iy+iz));
        bzc111 = cinv* 0.25*( bzM(ind+iz)+      bzM(ind+iy+iz)+   bzM(ind+1+iz)+    bzM(ind+1+iy+iz));
        

        wex000, wey000, wez000, kappa0, we2 = ExB_drift_rel( exc000, ey000, ez000, bx000, by000, bz000 );
        bnx000, bny000, bnz000, bpr000 = mag_unit_vec_rel( exc000, ey000, ez000, bx000, by000, bz000, we2);
        wex001, wey001, wez001, kappa0, we2 = ExB_drift_rel( exc001, ey001, ez001, bx001, by001, bz001 );
        bnx001, bny001, bnz001, bpr001 = mag_unit_vec_rel( exc001, ey001, ez001, bx001, by001, bz001, we2);
        wex010, wey010, wez010, kappa0, we2 = ExB_drift_rel( exc010, ey010, ez010, bx010, by010, bz010 );
        bnx010, bny010, bnz010, bpr010 = mag_unit_vec_rel( exc010, ey010, ez010, bx010, by010, bz010, we2);
        wex011, wey011, wez011, kappa0, we2 = ExB_drift_rel( exc011, ey011, ez011, bx011, by011, bz011 );
        bnx011, bny011, bnz011, bpr011 = mag_unit_vec_rel( exc011, ey011, ez011, bx011, by011, bz011, we2);
        wex100, wey100, wez100, kappa0, we2 = ExB_drift_rel( exc100, ey100, ez100, bx100, by100, bz100 );
        bnx100, bny100, bnz100, bpr100 = mag_unit_vec_rel( exc100, ey100, ez100, bx100, by100, bz100, we2);
        wex101, wey101, wez101, kappa0, we2 = ExB_drift_rel( exc101, ey101, ez101, bx101, by101, bz101 );
        bnx101, bny101, bnz101, bpr101 = mag_unit_vec_rel( exc101, ey101, ez101, bx101, by101, bz101, we2);
        wex110, wey110, wez110, kappa0, we2 = ExB_drift_rel( exc110, ey110, ez110, bx110, by110, bz110 );
        bnx110, bny110, bnz110, bpr110 = mag_unit_vec_rel( exc110, ey110, ez110, bx110, by110, bz110, we2);
        wex111, wey111, wez111, kappa0, we2 = ExB_drift_rel( exc111, ey111, ez111, bx111, by111, bz111 );
        bnx111, bny111, bnz111, bpr111 = mag_unit_vec_rel( exc111, ey111, ez111, bx111, by111, bz111, we2);

        wex1, wexx1, wexy1, wexz1 = _lerp(wex000, wex100, wex010, wex110, wex001, wex101, wex011, wex111, dx, dy, dz); // components of the bulk velocity and their gradients 
        wey1, weyx1, weyy1, weyz1 = _lerp(wey000, wey100, wey010, wey110, wey001, wey101, wey011, wey111, dx, dy, dz); 
        wez1, wezx1, wezy1, wezz1 = _lerp(wez000, wez100, wez010, wez110, wez001, wez101, wez011, wez111, dx, dy, dz);

        bnx1, bnxx1, bnxy1, bnxz1 = _lerp(bnx000, bnx100, bnx010, bnx110, bnx001, bnx101, bnx011, bnx111, dx, dy, dz); // components of the unit b vecot and their gradients 
        bny1, bnyx1, bnyy1, bnyz1 = _lerp(bny000, bny100, bny010, bny110, bny001, bny101, bny011, bny111, dx, dy, dz);
        bnz1, bnzx1, bnzy1, bnzz1 = _lerp(bnz000, bnz100, bnz010, bnz110, bnz001, bnz101, bnz011, bnz111, dx, dy, dz);

        bpr1, bprx1, bpry1, bprz1 = _lerp(bpr000, bpr100, bpr010, bpr110, bpr001, bpr101, bpr011, bpr111, dx, dy, dz); // B' and its gradient 

        ex1, exx1, exy1, exz1 = _lerp(exc000, exc100, exc010, exc110, exc001, exc101, exc011, exc111, dx, dy, dz); // electric field components and their gradients
        ey1, eyx1, eyy1, eyz1 = _lerp(eyc000, eyc100, eyc010, eyc110, eyc001, eyc101, eyc011, eyc111, dx, dy, dz);
        ez1, ezx1, ezy1, ezz1 = _lerp(ezc000, ezc100, ezc010, ezc110, ezc001, ezc101, ezc011, ezc111, dx, dy, dz);
        ex1 += this->get_ex_ext(0,0,0);
        ey1 += this->get_ey_ext(0,0,0);
        ez1 += this->get_ez_ext(0,0,0);
        
        bx1, bxx1, bxy1, bxz1 = _lerp(bxc000, bxc100, bxc010, bxc110, bxc001, bxc101, bxc011, bxc111, dx, dy, dz); // magnetic field components and their gradients
        by1, byx1, byy1, byz1 = _lerp(byc000, byc100, byc010, byc110, byc001, byc101, byc011, byc111, dx, dy, dz);
        bz1, bzx1, bzy1, bzz1 = _lerp(bzc000, bzc100, bzc010, bzc110, bzc001, bzc101, bzc011, bzc111, dx, dy, dz);
        bx1 += this->get_bx_ext(0,0,0);
        by1 += this->get_by_ext(0,0,0);
        bz1 += this->get_bz_ext(0,0,0);


    }






    //-------------------------------------------------- 
    //ExB in units of c

    // non-rel / rel ExB drift velocity
    // auto [vex0, vey0, vez0, kappa0, we2] = ExB_drift( ex0, ey0, ez0, bx0, by0, bz0 );
    auto [vex0, vey0, vez0, kappa0, we2] = ExB_drift_rel( ex0, ey0, ez0, bx0, by0, bz0 );

    //-------------------------------------------------- 
    // magnetic field unit vector b
      
    // auto [bnx0, bny0, bnz0] = mag_unit_vec(bx0, by0, bz0);
    auto [bnx0, bny0, bnz0, bpr0] = mag_unit_vec_rel( ex0, ey0, ez0, bx0, by0, bz0, we2);

    //--------------------------------------------------
    // epar = e.b
    real_long epar = ex0*bnx0 + ey0*bny0 + ez0*bnz0;

    // project velocity to the b field; upar at t_n-1/2 via u_par = (u . b)b
    real_long upar01  = vel0n*bnx0 + vel1n*bny0 + vel2n*bnz0;
    if (false) {  // upar01> 1.0 ){
	std::cout << "upar > 1 = " << upar01 
	    <<" vel0n: "<< vel0n << " vel1n: :"<< vel1n <<"vel2n: "<< vel2n<<"  ... "
    	    <<" bnx0: "<< bnx0 <<" bny0: "<< bny0 <<" bnz0: "<<bnz0
	    <<"\n\n";
	std::cout << std::flush;
	assert(false);	
    }



    //--------------------------------------------------
    // gyro motion four-velocity
      
    // NOTE: assuming full gamma factor here from previous step
    G0 = sqrt(1.0 + vel0n*vel0n + vel1n*vel1n + vel2n*vel2n );
    
    bgradb_x = upar01*upar01/G0*(bnx1*bnxx1 + bny1*bnxy1 + bny1*bnxy1); // u_par^2/Gamma * (b .\grad) b
    bgradb_y = upar01*upar01/G0*(bnx1*bnyx1 + bny1*bnyy1 + bny1*bnyy1); 
    bgradb_z = upar01*upar01/G0*(bnx1*bnzx1 + bny1*bnzy1 + bny1*bnzy1); 

    wegradb_x = upar01*(wex1*bnxx1 + wey1*bnxy1 + wey1*bnxy1); // u_par * (v_E .\grad) b
    wegradb_y = upar01*(wex1*bnyx1 + wey1*bnyy1 + wey1*bnyy1); 
    wegradb_z = upar01*(wex1*bnzx1 + wey1*bnzy1 + wey1*bnzy1); 

    v_cur_x = (kappa0*kappa0/qm/b0)*(bny*(bgradb_z + wegradb_z) - bnz*(bgradb_y + wegradb_y)); // bacchini '20, (8)
    v_cur_y = (kappa0*kappa0/qm/b0)*(bnz*(bgradb_x + wegradb_x) - bnx*(bgradb_z + wegradb_z)); // curvature drift
    v_cur_z = (kappa0*kappa0/qm/b0)*(bnx*(bgradb_y + wegradb_y) - bny*(bgradb_x + wegradb_x));


    ugx = vel0n - upar01*bnx0 - vex0*G0 - v_cur_x*G0;
    ugy = vel1n - upar01*bny0 - vey0*G0 - v_cur_y*G0;
    ugz = vel2n - upar01*bnz0 - vez0*G0 - v_cur_z*G0;
    ug2 = ugx*ugx + ugy*ugy + ugz*ugz;

    // magnetic moment = m u_g^2/2 B_0 \gamma
    // FIXME: this is taken as kappa0 in paper eqs; correct?
    //mu = me*ug2/(2.0*b0*G);
    mu = me*ug2/(2.0*b0*kappa0);


    v_mir_x = (mu/qm/me/G0)*(bny*(bprz1) - bnz*(bpry1)); // similar to Ripperda '18
    v_mir_y = (mu/qm/me/G0)*(bnz*(bprx1) - bnx*(bprz1)); // mirror drift
    v_mir_z = (mu/qm/me/G0)*(bnx*(bpry1) - bny*(bprx1));


    //--------------------------------------------------
    // upar at t_n+1/2
    // increase velocity with parallel E field: u += q/m dt Epar
    // NOTE: standard c * vel + qm*epar changed to vel + qm*epar/c
    // NOTE: cinv is multiplied to b0 in the beginning
    // FIXME: or multiply cinv here?
    // upar01 += qm*epar;
    const real_long k0 = sqrt(1.0 + upar01*upar01 + ug2 );     // gamma
    upar01 += qm*epar; //SWAP?
  

    //--------------------------------------------------
    // step1
    real_long R1x = R0x;
    real_long R1y = R0y;
    real_long R1z = R0z;


    for(size_t iter=0; iter<5; iter++){

      // re-use previous interpolation for first step
      if(iter == 0) { 
        // ex1 = ex0; ey1 = ey0; ez1 = ez0;
        // bx1 = bx0; by1 = by0; bz1 = bz0;
        // exx1 = 0.0; exy1 = 0.0; exz1 = 0.0; 
        // eyx1 = 0.0; eyy1 = 0.0; eyz1 = 0.0; 
        // ezx1 = 0.0; ezy1 = 0.0; ezz1 = 0.0; 
        // bxx1 = 0.0; bxy1 = 0.0; bxz1 = 0.0; 
        // byx1 = 0.0; byy1 = 0.0; byz1 = 0.0; 
        // bzx1 = 0.0; bzy1 = 0.0; bzz1 = 0.0; 
      } else {

        if(D >= 1) i  = static_cast<int>(floor(R1x));
        if(D >= 2) j  = static_cast<int>(floor(R1y));
        if(D >= 3) k  = static_cast<int>(floor(R1z));

        if(D >= 1) dx = R1x - i;
        if(D >= 2) dy = R1y - j;
        if(D >= 3) dz = R1z - k;

        // normalize to tile units
        if(D >= 1) i -= mins[0];
        if(D >= 2) j -= mins[1];
        if(D >= 3) k -= mins[2];

        //extrapolate if on tile boundary
        // min left side of tile
        if(D >= 1) { if(i <= -2 ) { dx -= i+2; i = -2; } }
        if(D >= 2) { if(j <= -2 ) { dy -= j+2; j = -2; } }
        if(D >= 3) { if(k <= -2 ) { dz -= k+2; k = -2; } }
        
        // max right side of tile
        if(D >= 1) { if(i >= Nx+1 ) { dx += i-Nx-1; i = Nx+1; } }
        if(D >= 2) { if(j >= Ny+1 ) { dy += j-Ny-1; j = Ny+1; } }
        if(D >= 3) { if(k >= Nz+1 ) { dz += k-Nz-1; k = Nz+1; } }


        const size_t ind = yee.ex.indx(i,j,k);

        //ex =
        exc000 = cinv* 0.5*(exM(ind       ) +exM(ind-1      ));
        exc100 = cinv* 0.5*(exM(ind       ) +exM(ind+1      ));
        exc010 = cinv* 0.5*(exM(ind+iy    ) +exM(ind-1+iy   ));
        exc110 = cinv* 0.5*(exM(ind+iy    ) +exM(ind+1+iy   ));
        exc001 = cinv* 0.5*(exM(ind+iz    ) +exM(ind-1+iz   ));
        exc101 = cinv* 0.5*(exM(ind+iz    ) +exM(ind+1+iz   ));
        exc011 = cinv* 0.5*(exM(ind+iy+iz ) +exM(ind-1+iy+iz));
        exc111 = cinv* 0.5*(exM(ind+iy+iz ) +exM(ind+1+iy+iz));

        //ey
        eyc000 = cinv* 0.5*(eyM(ind      ) +eyM(ind-iy     ));
        eyc100 = cinv* 0.5*(eyM(ind+1    ) +eyM(ind+1-iy   ));
        eyc010 = cinv* 0.5*(eyM(ind      ) +eyM(ind+iy     ));
        eyc110 = cinv* 0.5*(eyM(ind+1    ) +eyM(ind+1+iy   ));
        eyc001 = cinv* 0.5*(eyM(ind+iz   ) +eyM(ind-iy+iz  ));
        eyc101 = cinv* 0.5*(eyM(ind+1+iz ) +eyM(ind+1-iy+iz));
        eyc011 = cinv* 0.5*(eyM(ind+iz   ) +eyM(ind+iy+iz  ));
        eyc111 = cinv* 0.5*(eyM(ind+1+iz ) +eyM(ind+1+iy+iz));

        //ez
        ezc000 = cinv* 0.5*(ezM(ind      ) + ezM(ind-iz     ));
        ezc100 = cinv* 0.5*(ezM(ind+1    ) + ezM(ind+1-iz   ));
        ezc010 = cinv* 0.5*(ezM(ind+iy   ) + ezM(ind+iy-iz  ));
        ezc110 = cinv* 0.5*(ezM(ind+1+iy ) + ezM(ind+1+iy-iz));
        ezc001 = cinv* 0.5*(ezM(ind      ) + ezM(ind+iz     ));
        ezc101 = cinv* 0.5*(ezM(ind+1    ) + ezM(ind+1+iz   ));
        ezc011 = cinv* 0.5*(ezM(ind+iy   ) + ezM(ind+iy+iz  ));
        ezc111 = cinv* 0.5*(ezM(ind+1+iy ) + ezM(ind+1+iy+iz));

        //-------------------------------------------------- 
        // bx
        bxc000 = cinv* 0.25*( bxM(ind)+   bxM(ind-iy)+   bxM(ind-iz)+      bxM(ind-iy-iz));
        bxc100 = cinv* 0.25*( bxM(ind+1)+ bxM(ind+1-iy)+ bxM(ind+1-iz)+    bxM(ind+1-iy-iz));
        bxc001 = cinv* 0.25*( bxM(ind)+   bxM(ind+iz)+   bxM(ind-iy)+      bxM(ind-iy+iz));
        bxc101 = cinv* 0.25*( bxM(ind+1)+ bxM(ind+1+iz)+ bxM(ind+1-iy)+    bxM(ind+1-iy+iz));
        bxc010 = cinv* 0.25*( bxM(ind)+   bxM(ind+iy)+   bxM(ind-iz)+      bxM(ind+iy-iz));
        bxc110 = cinv* 0.25*( bxM(ind+1)+ bxM(ind+1-iz)+ bxM(ind+1+iy-iz)+ bxM(ind+1+iy));
        bxc011 = cinv* 0.25*( bxM(ind)+   bxM(ind+iy)+   bxM(ind+iy+iz)+   bxM(ind+iz));
        bxc111 = cinv* 0.25*( bxM(ind+1)+ bxM(ind+1+iy)+ bxM(ind+1+iy+iz)+ bxM(ind+1+iz));

        // by
        byc000 = cinv* 0.25*( byM(ind-1-iz)+    byM(ind-1)+       byM(ind-iz)+      byM(ind));
        byc100 = cinv* 0.25*( byM(ind-iz)+      byM(ind)+         byM(ind+1-iz)+    byM(ind+1));
        byc001 = cinv* 0.25*( byM(ind-1)+       byM(ind-1+iz)+    byM(ind)+         byM(ind+iz));
        byc101 = cinv* 0.25*( byM(ind)+         byM(ind+iz)+      byM(ind+1)+       byM(ind+1+iz));
        byc010 = cinv* 0.25*( byM(ind-1+iy-iz)+ byM(ind-1+iy)+    byM(ind+iy-iz)+   byM(ind+iy));
        byc110 = cinv* 0.25*( byM(ind+iy-iz)+   byM(ind+iy)+      byM(ind+1+iy-iz)+ byM(ind+1+iy));
        byc011 = cinv* 0.25*( byM(ind-1+iy)+    byM(ind-1+iy+iz)+ byM(ind+iy)+      byM(ind+iy+iz));
        byc111 = cinv* 0.25*( byM(ind+iy)+      byM(ind+iy+iz)+   byM(ind+1+iy)+    byM(ind+1+iy+iz));

        // bz
        bzc000 = cinv* 0.25*( bzM(ind-1-iy)+    bzM(ind-1)+       bzM(ind-iy)+      bzM(ind));
        bzc100 = cinv* 0.25*( bzM(ind-iy)+      bzM(ind)+         bzM(ind+1-iy)+    bzM(ind+1));
        bzc001 = cinv* 0.25*( bzM(ind-1-iy+iz)+ bzM(ind-1+iz)+    bzM(ind-iy+iz)+   bzM(ind+iz));
        bzc101 = cinv* 0.25*( bzM(ind-iy+iz)+   bzM(ind+iz)+      bzM(ind+1-iy+iz)+ bzM(ind+1+iz));
        bzc010 = cinv* 0.25*( bzM(ind-1)+       bzM(ind-1+iy)+    bzM(ind)+         bzM(ind+iy));
        bzc110 = cinv* 0.25*( bzM(ind)+         bzM(ind+iy)+      bzM(ind+1)+       bzM(ind+1+iy));
        bzc011 = cinv* 0.25*( bzM(ind-1+iz)+    bzM(ind-1+iy+iz)+ bzM(ind+iz)+      bzM(ind+iy+iz));
        bzc111 = cinv* 0.25*( bzM(ind+iz)+      bzM(ind+iy+iz)+   bzM(ind+1+iz)+    bzM(ind+1+iy+iz));
        

        wex000, wey000, wez000, kappa0, we2 = ExB_drift_rel( exc000, ey000, ez000, bx000, by000, bz000 );
        bnx000, bny000, bnz000, bpr000 = mag_unit_vec_rel( exc000, ey000, ez000, bx000, by000, bz000, we2);
        wex001, wey001, wez001, kappa0, we2 = ExB_drift_rel( exc001, ey001, ez001, bx001, by001, bz001 );
        bnx001, bny001, bnz001, bpr001 = mag_unit_vec_rel( exc001, ey001, ez001, bx001, by001, bz001, we2);
        wex010, wey010, wez010, kappa0, we2 = ExB_drift_rel( exc010, ey010, ez010, bx010, by010, bz010 );
        bnx010, bny010, bnz010, bpr010 = mag_unit_vec_rel( exc010, ey010, ez010, bx010, by010, bz010, we2);
        wex011, wey011, wez011, kappa0, we2 = ExB_drift_rel( exc011, ey011, ez011, bx011, by011, bz011 );
        bnx011, bny011, bnz011, bpr011 = mag_unit_vec_rel( exc011, ey011, ez011, bx011, by011, bz011, we2);
        wex100, wey100, wez100, kappa0, we2 = ExB_drift_rel( exc100, ey100, ez100, bx100, by100, bz100 );
        bnx100, bny100, bnz100, bpr100 = mag_unit_vec_rel( exc100, ey100, ez100, bx100, by100, bz100, we2);
        wex101, wey101, wez101, kappa0, we2 = ExB_drift_rel( exc101, ey101, ez101, bx101, by101, bz101 );
        bnx101, bny101, bnz101, bpr101 = mag_unit_vec_rel( exc101, ey101, ez101, bx101, by101, bz101, we2);
        wex110, wey110, wez110, kappa0, we2 = ExB_drift_rel( exc110, ey110, ez110, bx110, by110, bz110 );
        bnx110, bny110, bnz110, bpr110 = mag_unit_vec_rel( exc110, ey110, ez110, bx110, by110, bz110, we2);
        wex111, wey111, wez111, kappa0, we2 = ExB_drift_rel( exc111, ey111, ez111, bx111, by111, bz111 );
        bnx111, bny111, bnz111, bpr111 = mag_unit_vec_rel( exc111, ey111, ez111, bx111, by111, bz111, we2);

        wex1, wexx1, wexy1, wexz1 = _lerp(wex000, wex100, wex010, wex110, wex001, wex101, wex011, wex111, dx, dy, dz);
        wey1, weyx1, weyy1, weyz1 = _lerp(wey000, wey100, wey010, wey110, wey001, wey101, wey011, wey111, dx, dy, dz);
        wez1, wezx1, wezy1, wezz1 = _lerp(wez000, wez100, wez010, wez110, wez001, wez101, wez011, wez111, dx, dy, dz);

        bnx1, bnxx1, bnxy1, bnxz1 = _lerp(bnx000, bnx100, bnx010, bnx110, bnx001, bnx101, bnx011, bnx111, dx, dy, dz);
        bny1, bnyx1, bnyy1, bnyz1 = _lerp(bny000, bny100, bny010, bny110, bny001, bny101, bny011, bny111, dx, dy, dz);
        bnz1, bnzx1, bnzy1, bnzz1 = _lerp(bnz000, bnz100, bnz010, bnz110, bnz001, bnz101, bnz011, bnz111, dx, dy, dz);

        bpr1, bprx1, bpry1, bprz1 = _lerp(bpr000, bpr100, bpr010, bpr110, bpr001, bpr101, bpr011, bpr111, dx, dy, dz);

        ex1, exx1, exy1, exz1 = _lerp(exc000, exc100, exc010, exc110, exc001, exc101, exc011, exc111, dx, dy, dz);
        ey1, eyx1, eyy1, eyz1 = _lerp(eyc000, eyc100, eyc010, eyc110, eyc001, eyc101, eyc011, eyc111, dx, dy, dz);
        ez1, ezx1, ezy1, ezz1 = _lerp(ezc000, ezc100, ezc010, ezc110, ezc001, ezc101, ezc011, ezc111, dx, dy, dz);
        ex1 += this->get_ex_ext(0,0,0);
        ey1 += this->get_ey_ext(0,0,0);
        ez1 += this->get_ez_ext(0,0,0);
        
        bx1, bxx1, bxy1, bxz1 = _lerp(bxc000, bxc100, bxc010, bxc110, bxc001, bxc101, bxc011, bxc111, dx, dy, dz);
        by1, byx1, byy1, byz1 = _lerp(byc000, byc100, byc010, byc110, byc001, byc101, byc011, byc111, dx, dy, dz);
        bz1, bzx1, bzy1, bzz1 = _lerp(bzc000, bzc100, bzc010, bzc110, bzc001, bzc101, bzc011, bzc111, dx, dy, dz);
        bx1 += this->get_bx_ext(0,0,0);
        by1 += this->get_by_ext(0,0,0);
        bz1 += this->get_bz_ext(0,0,0);

        bgradb_x = upar01*upar01/G0*(bnx1*bnxx1 + bny1*bnxy1 + bny1*bnxy1); // u_par^2/Gamma * (b .\grad) b
        bgradb_y = upar01*upar01/G0*(bnx1*bnyx1 + bny1*bnyy1 + bny1*bnyy1); 
        bgradb_z = upar01*upar01/G0*(bnx1*bnzx1 + bny1*bnzy1 + bny1*bnzy1); 

        wegradb_x = upar01*(wex1*bnxx1 + wey1*bnxy1 + wey1*bnxy1); // u_par * (v_E .\grad) b
        wegradb_y = upar01*(wex1*bnyx1 + wey1*bnyy1 + wey1*bnyy1); 
        wegradb_z = upar01*(wex1*bnzx1 + wey1*bnzy1 + wey1*bnzy1); 

        v_cur_x1 = (kappa0*kappa0/qm/b0)*(bny*(bgradb_z + wegradb_z) - bnz*(bgradb_y + wegradb_y)); // bacchini '20, (8)
        v_cur_y1 = (kappa0*kappa0/qm/b0)*(bnz*(bgradb_x + wegradb_x) - bnx*(bgradb_z + wegradb_z)); // curvature drift
        v_cur_z1 = (kappa0*kappa0/qm/b0)*(bnx*(bgradb_y + wegradb_y) - bny*(bgradb_x + wegradb_x));

        // ex1 *= cinv;
        // ey1 *= cinv;
        // ez1 *= cinv;
        // bx1 *= cinv;
        // by1 *= cinv;
        // bz1 *= cinv;
      }

      //-------------------------------------------------- 

      //-------------------------------------------------- 
      //ExB in units of c at new location

      // non-rel / rel ExB drift velocity at the new location
      // auto [vex1, vey1, vez1, kappa1, we2] = ExB_drift( ex1, ey1, ez1, bx1, by1, bz1 );
      auto [vex1, vey1, vez1, kappa1, we2] = ExB_drift_rel( ex1, ey1, ez1, bx1, by1, bz1 );

      //-------------------------------------------------- 
      // magnetic field unit vector b at new location

      // auto [bnx1, bny1, bnz1] = mag_unit_vec(bx1, by1, bz1);
      auto [bnx1, bny1, bnz1, bpr1] = mag_unit_vec_rel( ex1, ey1, ez1, bx1, by1, bz1, we2);

      //-------------------------------------------------- 
      // location update
      // NOTE: assuming mu -> 0 zero during the time step

      // Construct new gyrovelocity from conservation of magnetic moment
      // magnetic moment = m u_g^2/2 B_0 \gamma
        
      // mu = 0.0; // NOTE: synchrotron losses are assumed to bring mag. mom. to zero
      real_long b1 = sqrt( bx1*bx1 + by1*by1 + bz1*bz1 );
      ug2n = mu*2.0*b1*kappa1/container.m;

      const real_long k1 = sqrt(1.0 + upar01*upar01 + ug2n);     // gamma
      EPkappa1 =  1.0/(sqrt(1. - 0.5*(vex1*vex1 + vey1*vey1 + vez1*vez1 + vex0*vex0 + vey0*vey0 + vez0*vez0)) + EPS); 

      G0 = k0*kappa0; // inv Gamma at t = n FIXME: or previous G0?
      G1 = k1*kappa1; // inv Gamma at t = n+1


      // GCA coordinate-velocity at v_n+1/2 
      // vn1x = 0.5*upar01*( bnx0/G0 + bnx1/G1 ) + 0.5*(vex0 + vex1);
      // vn1y = 0.5*upar01*( bny0/G0 + bny1/G1 ) + 0.5*(vey0 + vey1);
      // vn1z = 0.5*upar01*( bnz0/G0 + bnz1/G1 ) + 0.5*(vez0 + vez1);
      // vn1x = 0.5*upar01*( bnx0/G1 + bnx1/G1 ) + 0.5*(vex0 + vex1);
      // vn1y = 0.5*upar01*( bny0/G1 + bny1/G1 ) + 0.5*(vey0 + vey1);
      // vn1z = 0.5*upar01*( bnz0/G1 + bnz1/G1 ) + 0.5*(vez0 + vez1);
      vn1x = 0.5*upar01*( bnx0/G1 + bnx1/G1 ) + 0.5*(vex0 + vex1) + 0.5*(v_cur_x + v_cur_x1);
      vn1y = 0.5*upar01*( bny0/G1 + bny1/G1 ) + 0.5*(vey0 + vey1) + 0.5*(v_cur_y + v_cur_y1);
      vn1z = 0.5*upar01*( bnz0/G1 + bnz1/G1 ) + 0.5*(vez0 + vez1) + 0.5*(v_cur_z + v_cur_z1);
      G01 = 1.0/sqrt(1.0 - vn1x*vn1x - vn1y*vn1y - vn1z*vn1z);
      //-------------------------------------------------- 
      // newest GCA four velocity at u_n+1 (for velocity update)
      // gyromotion is taken to be same direction as in the previous step (hence ug_i/sqrt(ug2))
      // un1x = upar01*bnx1 + vex1*G1 + sqrt(ug2n)*ugx/sqrt(ug2);
      // un1y = upar01*bny1 + vey1*G1 + sqrt(ug2n)*ugy/sqrt(ug2);
      // un1z = upar01*bnz1 + vez1*G1 + sqrt(ug2n)*ugz/sqrt(ug2);
      // un1x = 0.5*vn1x*(G0+G1) + sqrt(ug2n)*ugx/sqrt(ug2);
      // un1y = 0.5*vn1y*(G0+G1) + sqrt(ug2n)*ugy/sqrt(ug2);
      // un1z = 0.5*vn1z*(G0+G1) + sqrt(ug2n)*ugz/sqrt(ug2);
      // un1x = vn1x*(G01) + sqrt(ug2n)*ugx/sqrt(ug2);
      // un1y = vn1y*(G01) + sqrt(ug2n)*ugy/sqrt(ug2);
      // un1z = vn1z*(G01) + sqrt(ug2n)*ugz/sqrt(ug2);
      un1x = vn1x*(G1)  + sqrt(ug2n)*ugx/sqrt(ug2);
      un1y = vn1y*(G1)  + sqrt(ug2n)*ugy/sqrt(ug2);
      un1z = vn1z*(G1)  + sqrt(ug2n)*ugz/sqrt(ug2);
      un12 = un1x*un1x + un1y*un1y + un1z*un1z; 


      // real_long gradb2x = bx1*bxx + by*byx + bz*bzx;
      // real_long gradb2y = bx1*bxy + by*byy + bz*bzy;
      // real_long gradb2z = bx1*bxz + by*byz + bz*bzz;

      // real_long gradbxx = (bxx - bx1*gradb2x/(b1*b1))/b1;
      // real_long gradbxy = (bxy - bx1*gradb2y/(b1*b1))/b1;
      // real_long gradbxz = (bxz - bx1*gradb2z/(b1*b1))/b1;
      // real_long gradbyx = (byx - by1*gradb2x/(b1*b1))/b1;
      // real_long gradbyy = (byy - by1*gradb2y/(b1*b1))/b1;
      // real_long gradbyz = (byz - by1*gradb2z/(b1*b1))/b1;
      // real_long gradbzx = (bzx - bz1*gradb2x/(b1*b1))/b1;
      // real_long gradbzy = (bzy - bz1*gradb2y/(b1*b1))/b1;
      // real_long gradbzz = (bzz - bz1*gradb2z/(b1*b1))/b1;
      //-------------------------------------------------- 
      // location error
      real_long Hx = R1x - (R0x + c*vn1x);
      real_long Hy = R1y - (R0y + c*vn1y);
      real_long Hz = R1z - (R0z + c*vn1z);
      real_long H = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);

      //-------------------------------------------------- 
      // guiding center location update
      if(D>=1) R1x = R0x + vn1x*c;
      if(D>=2) R1y = R0y + vn1y*c;
      if(D>=3) R1z = R0z + vn1z*c;


      //-------------------------------------------------- 
      if(false) {
      //if(n == 0) {
      //if(mu > 1.0) {
        //crash_flag = true;
          
        real_long b1 = sqrt( bx1*bx1 + by1*by1 + bz1*bz1 );
        real_long e1 = sqrt( ex1*ex1 + ey1*ey1 + ez1*ez1 );

        // E.B violation; i.e .how bad is the pure ExB drift assumption
        real_long EB0 = (ex0*bx0 + ey0*by0 + ez0*bz0)/b0/b0;
        real_long EB1 = (ex1*bx1 + ey1*by1 + ez1*bz1)/b1/b1;
        real_long bn0 = bnx0*bnx0 + bny0*bny0 + bnz0*bnz0;
        real_long bn1 = bnx1*bnx1 + bny1*bny1 + bnz1*bnz1;

        real_long egtb0 = e0/b0;
        real_long egtb1 = e1/b1;

        std::cout 
        << "iter:" << iter << " H:" << H << " E.B:" << EB0 << " " << EB1
        << " e/b0:" << egtb0 << " e/b1:" << egtb1
        << "\n"
        << " R0x:" << R0x << " R0y:" << R0y << " R0z:" << R0z 
        << "\n"
        << " R1x:" << R1x << " R1y:" << R1y << " R1z:" << R1z  
        << "\n"
        << " Hx:"  << Hx  << " Hy:"  << Hy  << " Hz:"  << Hz
        << "\n"
        //<< " ex0:" << ex0 << " ey0:" << ey0 << " ez0:" << ez0
        //<< " bx0:" << bx0 << " by0:" << by0 << " bz0:" << bz0
        << " bnx0:" << bnx0 << " bny:" << bny0 << " bnz:" << bnz0 << " b: " << bn0
        //<< "\n"
        //<< " *bnx0:" << bnx0_2 << " bny:" << bny0_2 << " bnz:" << bnz0_2
        << "\n"
        << " bnx1:" << bnx1 << " bny:" << bny1 << " bnz:" << bnz1 << " b: " << bn1
        << "\n"
        //<< " eparx:" << eparx << " epary:" << epary << " eparz:" << eparz
        //<< " uparx:" << uparx01 << " upary:" << upary01 << " uparz:" << uparz01 << 
        " upar:" << upar01
        << " k0:" << k0 
        << " k1:" << k1 
        << "\n"
        << " mu:" << mu 
        << " ug:" << sqrt(ug2) << " ugn:" << sqrt(ug2n)
        << "\n"
        << " vex0:" << vex0 << " vey:" << vey0 << " vez:" << vez0 << " kappa0:" << kappa0
        << "\n"
        << " vex1:" << vex1 << " vey:" << vey1 << " vez:" << vez1 << " kappa1:" << kappa1
        << "\n\n";
      std::cout << std::flush;

      }

      // exit if converged
      if(H < 1e-5) break;

    }//end of iteration

    if(crash_flag) assert(false);

    un1x = std::min(u_upper_bound, un1x); 
    un1y = std::min(u_upper_bound, un1y); 
    un1z = std::min(u_upper_bound, un1z); 
    un1x = std::max(-u_upper_bound, un1x); 
    un1y = std::max(-u_upper_bound, un1y); 
    un1z = std::max(-u_upper_bound, un1z); 
    vel[0][n] = static_cast<real_prtcl>( un1x );
    vel[1][n] = static_cast<real_prtcl>( un1y );
    vel[2][n] = static_cast<real_prtcl>( un1z );

    // position update from iteration, new location is following gyro center position
    if(D>=1) loc[0][n] = R1x;
    if(D>=2) loc[1][n] = R1y;
    if(D>=3) loc[2][n] = R1z;  

    // store also the field values at the new point 
    exP[n] = static_cast<real_prtcl>( ex1 );
    eyP[n] = static_cast<real_prtcl>( ey1 );
    ezP[n] = static_cast<real_prtcl>( ez1 );
    bxP[n] = static_cast<real_prtcl>( bx1 );
    byP[n] = static_cast<real_prtcl>( by1 );
    bzP[n] = static_cast<real_prtcl>( bz1 );

    bool debug_flag = /// (un12 > 1.0e80) ||
    std::isinf(vel[0][n]) ||
    std::isinf(vel[1][n]) ||
    std::isinf(vel[2][n]) ||
    std::isnan(vel[0][n]) ||
    std::isnan(vel[1][n]) ||
    std::isnan(vel[2][n]) ||
    std::isinf(loc[0][n]) ||
    std::isinf(loc[1][n]) ||
    std::isinf(loc[2][n]) ||
    std::isnan(loc[0][n]) ||
    std::isnan(loc[1][n]) ||
    std::isnan(loc[2][n]);   

    //if(1./kinv01 > 30.0) debug_flag = true;
    if(debug_flag){
      std::cout 
        << " loc0n:" << loc[0][n] << " loc1n:" << loc[1][n] << " loc2n:" << loc[2][n]
        << " un1x:" << un1x << " un1y:" << un1y << " un1z:" << un1z
      	<< " un12: " << un12
        << " vel0n:" << vel[0][n] << " vel1n:" << vel[1][n] << " vel2n:" << vel[2][n]
        << " G1:" << G1 << " ugx:" << ugx << " ug2:" << ug2
        << " ex0:" << ex0 << " ey0:" << ey0 << " ez0:" << ez0
        << " bx0:" << bx0 << " by0:" << by0 << " bz0:" << bz0
        << " bnx0:" << bnx0 << " bny:" << bny0 << " bnz:" << bnz0
        //<< " bnx1:" << bnx1 << " bny:" << bny1 << " bnz:" << bnz1
        //<< " eparx:" << eparx << " epary:" << epary << " eparz:" << eparz
        //<< " uparx:" << uparx01 << " upary:" << upary01 << " uparz:" << uparz01
        << " upar:" << upar01 << " qm:"<<qm <<" epar:"<<epar
        << " k0:" << k0 
        //<< " k1:" << k1 
        << " vex0:" << vex0 << " vey:" << vey0 << " vez:" << vez0 << " kappa0:" << kappa0
        //<< " vex1:" << vex1 << " vey:" << vey1 << " vez:" << vez1 << " kappa1:" << kappa1
        << "\n";
      std::cout << std::flush;
      assert(false);
    }
  }

}

/*
template<size_t D, size_t V>
void pic::rGCAPusher<D,V>::calc_gradients(
    pic::Tile<D>& tile
    )
{

    const int Nx = tile.mesh_lengths[0];
    const int Ny = tile.mesh_lengths[1];
    const int Nz = tile.mesh_lengths[2];

    // fields at grid
    auto& yee = tile.get_yee();

    auto& exM = yee.ex;
    auto& eyM = yee.ey;
    auto& ezM = yee.ez;

    auto& bxM = yee.bx;
    auto& byM = yee.by;
    auto& bzM = yee.bz;


    const size_t tot = Nx*Ny*Nz //?
    // mesh sizes for 1D indexing
    const size_t iy = D >= 2 ? yee.ex.indx(0,1,0) - yee.ex.indx(0,0,0) : 0;
    const size_t iz = D >= 3 ? yee.ex.indx(0,0,1) - yee.ex.indx(0,0,0) : 0;

    real_long* dbxx = new real_long[];
    real_long dbxxN[Nx*Ny*Nz];
    real_long dbxyN[Nx*Ny*Nz];
    real_long dbxzN[Nx*Ny*Nz];
    real_long dbyyN[Nx*Ny*Nz];
    real_long dbyzN[Nx*Ny*Nz];
    real_long dbzzN[Nx*Ny*Nz];
    real_long dexxN[Nx*Ny*Nz];
    real_long dexyN[Nx*Ny*Nz];
    real_long dexzN[Nx*Ny*Nz];
    real_long deyyN[Nx*Ny*Nz];
    real_long deyzN[Nx*Ny*Nz];
    real_long dezzN[Nx*Ny*Nz];


    return;
}
*/



//--------------------------------------------------
// explicit template instantiation

template class pic::rGCAPusher<1,3>; // 1D3V
template class pic::rGCAPusher<2,3>; // 2D3V
template class pic::rGCAPusher<3,3>; // 3D3V



