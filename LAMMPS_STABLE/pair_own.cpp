/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "pair_own.h"                                                                  // yukawa->Own
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairOwn::PairOwn(LAMMPS *lmp) : Pair(lmp)                             //pairyukawa->pairOwn
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairOwn::~PairOwn()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

//    memory->destroy(rad);							// org
    memory->destroy(cut);
    memory->destroy(sigma);
    memory->destroy(alpha);
    memory->destroy(offset);                                                    //use Uc
  }
}

/* ---------------------------------------------------------------------- */

void PairOwn::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,rinv,r2inv,screening,forceOwn,factor;                               //forceyukwa->force
  double cons1,var;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;                                       // 1-2,1-3,1-4 prefactor for lj ?
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor = special_lj[sbmask(j)];                                           // input units -> LJ units
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
//      delx = delx - int(delx);
//      dely = dely - int(dely);
//      delz = delz - int(delz);
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
//        r2inv = 1.0/rsq;
        r = sqrt(rsq);
        rinv = 1.0/r;
	var = 1.0 - r/sigma[itype][jtype];
	cons1 = alpha[itype][jtype] - 1.0;
//        screening = exp(-kappa*r);
//        force = a[itype][jtype] * screening * (kappa + rinv);
	  forceOwn = pow(var,cons1)/sigma[itype][jtype];

        fpair = factor * forceOwn * rinv;                                           

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = pow(var,alpha[itype][jtype]) / alpha[itype][jtype] - offset[itype][jtype];                                         // energy except coulobmic
          evdwl *= factor;                                                      //factor
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);                   //updata virial term
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairOwn::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

//  memory->create(rad,n+1,"pair:rad");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(alpha,n+1,n+1,"pair:alpha");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairOwn::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");                //2->1

//  kappa = force->numeric(FLERR,arg[0]);
  cut_global = force->numeric(FLERR,arg[0]);                                    //1->0  set cut

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
---------------------------------------------------------------------- --- */

void PairOwn::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);                                   //if d(i) is not same, define two type and sigma,
  force->bounds(arg[1],atom->ntypes,jlo,jhi);                                   // look lj/cut.cpp

  double sigma_one = force->numeric(FLERR,arg[2]);
  double alpha_one = force->numeric(FLERR,arg[3]);

  double cut_one = cut_global;
  if (narg == 5) cut_one = force->numeric(FLERR,arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      sigma[i][j] = sigma_one;
      alpha[i][j] = alpha_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;                                                        // ?
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairOwn::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    alpha[i][j] = mix_energy(alpha[i][i],alpha[j][j],1.0,1.0);                              // where?
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);                              // where?
  } 

  if (offset_flag) {                                                            // when use restart
    offset[i][j] = 0.0;
  } else offset[i][j] = 0.0;

  sigma[j][i] = sigma[i][j];
  alpha[j][i] = alpha[i][j];
  cut[j][i] = cut[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */    // not use now

void PairOwn::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&alpha[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOwn::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&alpha[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&alpha[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOwn::write_restart_settings(FILE *fp)
{
//  fwrite(&kappa,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOwn::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
//    fread(&kappa,sizeof(double),1,fp);
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
//  MPI_Bcast(&kappa,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairOwn::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,sigma[i][i],alpha[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairOwn::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",i,j,sigma[i][j],alpha[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairOwn::single(int i, int j, int itype, int jtype, double rsq,
                          double factor_coul, double factor_lj,
                          double &fforce)
{
  double r2inv,r,rinv,forceOwn,phi;
  double var2, cons;

//  r2inv = 1.0/rsq;
  r = sqrt(rsq);
  rinv = 1.0/r;
  var2 = 1.0 - r/sigma[itype][jtype];
  cons = alpha[itype][jtype] - 1.0;
  forceOwn = pow(var2,cons);
  fforce = factor_lj*forceOwn * rinv;                                             //factor_lj ?

  phi = pow(var2,alpha[itype][jtype]) / alpha[itype][jtype];
  return factor_lj*phi;                                                        // factor_lj ?
}
