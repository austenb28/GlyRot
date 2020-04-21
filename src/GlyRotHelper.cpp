#include <ctype.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "gromacs/commandline.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/smalloc.h"

#include "gromacs/gmxpreprocess/readir.h"

using namespace gmx;

static int sudorun(int argc, char *argv[]);
static atom_id my_ind(t_blocka *grps, int i, int j);
static void rotate_glycan(t_blocka *grps, int *isize, int gnum, t_state *state, double dtheta, int x, int y, rvec *xp);
static void rotate_atom(int atom, t_state *state, double dtheta, int x, int y, rvec *xp);
static double my_dist(rvec **xp, int x, int y);

class GlyRotModule : public CommandLineModuleInterface
{
	public:
		GlyRotModule(){}
		virtual const char *name() const{return "GlyRot";}
		virtual const char *shortDescription() const{return "Glycan Rotation Module";}
		virtual void init(CommandLineModuleSettings *settings){}
		virtual int run(int argc, char *argv[]);
		virtual void writeHelp(const CommandLineHelpContext &context) const;
};

void GlyRotModule::writeHelp(const CommandLineHelpContext &context) const
{
	writeCommandLineHelpCMain(context, name(), &sudorun );
}

int GlyRotModule::run(int argc, char *argv[])
{
	const char *desc[] = {
        "[THISMODULE] reads a [REF].tpr[ref] of a dry glycoprotein and",
	"a properly defined index file.  [THISMODULE] rotates N-glycans",
	" about three bonds such that a minimum energy rotamer is obtained."};
	int NDESC = (int)(sizeof(desc)/sizeof(desc[0]));
	
	t_filenm fnm[] = {
                          { efTPR, "-s", "topol", ffREAD },
	                  { efNDX, "-n", "index", ffREAD },
	                  { efXTC, "-to", "trx_out", ffWRITE }
	                 };
	int NFILE = (int)(sizeof(fnm)/sizeof(fnm[0]));
	static int maxwarn = 0;

	static real dtheta = 90; //degrees
	static real theta1 = 0;
	static real theta2 = 0;
	static real theta3 = 0;
	t_pargs pa[] = {
                        { "-dtheta", FALSE, etREAL, {&dtheta}, 
	                  "Angular rotation resolution"}
	               };
	int NPARG = (int)(sizeof(pa)/sizeof(pa[0]));

	output_env_t oenv = NULL;

	char buf[256];
	t_topology top;
	int ePBC;
	rvec *xp;
	matrix box = {{0}};

	t_inputrec *ir;
	t_state *state = NULL;
        gmx_mtop_t *mtop = NULL;

	t_blocka *grps;
	char    **gnames;
	int gnum;
	int *isize;
	
	int *ND2_list, *HD21_list, *CG_list, *CA_list, *CB_list, 
	    *HB2_list, *HB3_list, *OD1_list;
	int nND2, nHD21, nCG, nCA, nCB, nHB2, nHB3, nOD1;
	int j,k,l,m,n,o,p,q;
	atom_id C1  = -1,
			CG  = -1, 
			ND2 = -1, 
			HD21= -1, 
			CA  = -1, 
			CB  = -1,   
			HB2 = -1, 
			HB3 = -1, 
			OD1 = -1;

	t_trxstatus *tstatus = NULL;
	t_trxframe fr;

	if (!parse_common_args(&argc, argv, 0, NFILE, fnm, NPARG, pa, NDESC, desc,
                           0, NULL, &oenv))
	{
		return 0;
	}

	init_warning(TRUE, maxwarn);
	
	/* code below is not executed for help call */

	snew(state, 1);
	snew(ir, 1);
	snew(mtop, 1);

	/* read structure/topology (two times, inefficient)*/
	read_tps_conf(opt2fn("-s", NFILE, fnm), buf, &top, &ePBC, &xp,
                             NULL, box, TRUE);

	read_tpx_state(opt2fn("-s", NFILE, fnm), ir, state, NULL, mtop);

	/*read index from ndx*/
	grps = init_index(ftp2fn(efNDX, NFILE, fnm), &gnames);
	
	snew(isize, grps->nr);

	fprintf(stderr, "\n");
	for (int i = 0; (i < grps->nr); i++)
	{
		isize[i] = grps->index[i+1]-grps->index[i];
		fprintf(stderr, "Group %5d (%15s) has %5d elements\n", i, gnames[i],
		isize[i]);
	}

	gnum = 1;

	nND2 = 0;
	nHD21 = 0;
	nCG = 0;
	nOD1 = 0;
	nCB = 0;
	nCA = 0;
	nHB2 = 0;
	nHB3 = 0;
	for(int i = 0; i < isize[0]; i++)
	{
		if(!strcmp(*(top.atoms.atomname[i]),"ND2"))
		{
			nND2++;
		}
		if(!strcmp(*(top.atoms.atomname[i]),"HD21"))
		{
			nHD21++;
		}
		if(!strcmp(*(top.atoms.atomname[i]),"CG"))
		{
			nCG++;
		}
		if(!strcmp(*(top.atoms.atomname[i]),"OD1"))
		{
			nOD1++;
		}
		if(!strcmp(*(top.atoms.atomname[i]),"CB"))
		{
			nCB++;
		}
		if(!strcmp(*(top.atoms.atomname[i]),"CA"))
		{
			nCA++;
		}
		if(!strcmp(*(top.atoms.atomname[i]),"HB2"))
		{
			nHB2++;
		}
		if(!strcmp(*(top.atoms.atomname[i]),"HB3"))
		{
			nHB3++;
		}
	}
	snew(ND2_list,nND2);
	snew(HD21_list,nHD21);
	snew(CG_list,nCG);
	snew(OD1_list,nOD1);
	snew(CB_list,nCB);
	snew(CA_list,nCA);
	snew(HB2_list,nHB2);
	snew(HB3_list,nHB3);
	/*initialize a list of ND2s and HD21s
	and CGs (a little inefficient)*/
	j=0;k=0;l=0;m=0;n=0;o=0;p=0;q=0;
	for(int i = 0; i < isize[0]; i++)
	{
		if(!strcmp(*(top.atoms.atomname[i]),"ND2"))
		{
			ND2_list[j] = i;
			j++;
		}
		if(!strcmp(*(top.atoms.atomname[i]),"HD21"))
		{
			HD21_list[k] = i;
			k++;
		}
		if(!strcmp(*(top.atoms.atomname[i]),"CG"))
		{
			CG_list[l] = i;
			l++;
		}
		if(!strcmp(*(top.atoms.atomname[i]),"OD1"))
		{
			OD1_list[m] = i;
			m++;
		}
		if(!strcmp(*(top.atoms.atomname[i]),"CB"))
		{
			CB_list[n] = i;
			n++;
		}
		if(!strcmp(*(top.atoms.atomname[i]),"CA"))
		{
			CA_list[o] = i;
			o++;
		}
		if(!strcmp(*(top.atoms.atomname[i]),"HB2"))
		{
			HB2_list[p] = i;
			p++;
		}
		if(!strcmp(*(top.atoms.atomname[i]),"HB3"))
		{
			HB3_list[q] = i;
			q++;
		}
	}

	/*find C1 of first sugar of glycan*/

	for(int i = 0; i < isize[gnum]; i++)
	{
		C1 = my_ind(grps,gnum,i);
		if(!strcmp(*(top.atoms.atomname[C1]),"C1"))
		{
			break;
		}
	}

	/*find nearest ND2 to C1*/
	for(int i = 0; i < nND2; i++)
	{
		ND2 = ND2_list[i];
		if(my_dist(&xp,ND2,C1) < 0.16)
		{
			break;
		}
	}

	/*find nearest HD21 to ND2 (a bit inefficient)*/
	for(int i = 0; i < nHD21; i++)
	{
		HD21 = HD21_list[i];
		if(my_dist(&xp,HD21,ND2) < 0.16)
		{
			break;
		}
	}
	/*find nearest CG to ND2 (a bit inefficient)*/
	for(int i = 0; i < nCG; i++)
	{
		CG = CG_list[i];
		if(my_dist(&xp,CG,ND2) < 0.16)
		{
			break;
		}
	}

	/*find nearest CB to CG (a bit inefficient)*/
	for(int i = 0; i < nCB; i++)
	{
		CB = CB_list[i];
		if(my_dist(&xp,CB,CG) < 0.16)
		{
			break;
		}
	}

	/*find nearest OD1 to CG (a bit inefficient)*/
	for(int i = 0; i < nOD1; i++)
	{
		OD1 = OD1_list[i];
		if(my_dist(&xp,OD1,CG) < 0.16)
		{
			break;
		}
	}
	/*find nearest CA to CB (a bit inefficient)*/
	for(int i = 0; i < nCA; i++)
	{
		CA = CA_list[i];
		if(my_dist(&xp,CA,CB) < 0.16)
		{
			break;
		}
	}
	/*find nearest HB2 to CB (a bit inefficient)*/
	for(int i = 0; i < nHB2; i++)
	{
		HB2 = HB2_list[i];
		if(my_dist(&xp,HB2,CB) < 0.16)
		{
			break;
		}
	}
	/*find nearest HB3 to CB (a bit inefficient)*/
	for(int i = 0; i < nHB3; i++)
	{
		HB3 = HB3_list[i];
		if(my_dist(&xp,HB3,CB) < 0.16)
		{
			break;
		}
	}

	/*write first output frame*/
	clear_trxframe(&fr, TRUE);	
        fr.bAtoms = TRUE;
        fr.atoms  = &(mtop->moltype[0].atoms);
        fr.bX     = TRUE;
        fr.bBox   = TRUE;
	fr.natoms = state->natoms;
	fr.time = 0;
	fr.x = state->x;
	fr.step = 0;
	for(int i = 0; i <= 2; i++)
	{
		for(int j = 0; j <= 2; j++)
		{
			fr.box[i][j] = box[i][j];
		}
	}
	tstatus = open_trx(opt2fn("-to", NFILE, fnm), "w");
	write_trxframe(tstatus, &fr, NULL);
	fr.step++;

	dtheta = dtheta*3.14159265359/180;
	k = 1;
	while(true)
	{
		while(true)
		{
			while(true)
			{
				/*rotate glycan along ND2 - C1 bond*/
				theta3 += dtheta;
				if(theta3 >= 3.14159265359*2 - dtheta/2)
				{
					break;
				}
				rotate_glycan(grps, isize, gnum, state, theta3, ND2, C1, xp);
				fr.time = k++;
				write_trxframe(tstatus, &fr, NULL);
				fr.step++;
			}

			/*reset and rotate glycan and other atoms along CB - CG bond, in xp!*/
			theta2 += dtheta;
			if(theta2 >= 3.14159265359*2 - dtheta/2)
			{
				break;
			}
			theta3 = 0;
			rotate_glycan(grps, isize, gnum, state, theta3, ND2, C1, xp);
			rotate_glycan(grps, isize, gnum, state, dtheta, CB, CG, xp);
			rotate_atom(ND2, state, dtheta, CB, CG, xp);
			rotate_atom(HD21, state, dtheta, CB, CG, xp);
			rotate_atom(OD1, state, dtheta, CB, CG, xp);

			fr.time = k++;
			write_trxframe(tstatus, &fr, NULL);
			fr.step++;

			for(int j = 0; j < isize[gnum]; j++)
			{
				for(int i = 0; i < 3; i++)
				{
					xp[my_ind(grps,gnum,j)][i] = state->x[my_ind(grps,gnum,j)][i];
				}
			}
			for(int i = 0; i < 3; i++)
			{
				xp[ND2][i] = state->x[ND2][i];
				xp[HD21][i] = state->x[HD21][i];
				xp[OD1][i] = state->x[OD1][i];
			}
		}
		/*reset and rotate glycan and other atoms along CA - CB bond, in xp!*/
		theta1 += dtheta;
		if(theta1 >= 3.14159265359*2 - dtheta/2)
		{
			break;
		}
		theta2 = 0;
		theta3 = 0;
		rotate_glycan(grps, isize, gnum, state, theta3, ND2, C1, xp);
		rotate_glycan(grps, isize, gnum, state, dtheta, CB, CG, xp);
		rotate_atom(ND2, state,  theta2, CB, CG, xp);
		rotate_atom(HD21, state, theta2, CB, CG, xp);
		rotate_atom(OD1, state,  theta2, CB, CG, xp);
		rotate_glycan(grps, isize, gnum, state, dtheta, CA, CB, xp);
		rotate_atom(ND2, state,  dtheta, CA, CB, xp);
		rotate_atom(HD21, state, dtheta, CA, CB, xp);
		rotate_atom(OD1, state,  dtheta, CA, CB, xp);
		rotate_atom(CG , state,  dtheta, CA, CB, xp);
		rotate_atom(HB2, state,  dtheta, CA, CB, xp);
		rotate_atom(HB3, state,  dtheta, CA, CB, xp);
		fr.time = k++;
		write_trxframe(tstatus, &fr, NULL);
		fr.step++;
		for(int j = 0; j < isize[gnum]; j++)
		{
			for(int i = 0; i < 3; i++)
			{
				xp[my_ind(grps,gnum,j)][i] = state->x[my_ind(grps,gnum,j)][i];
			}
		}
		for(int i = 0; i < 3; i++)
		{
			xp[ND2][i] = state->x[ND2][i];
			xp[HD21][i] = state->x[HD21][i];
			xp[OD1][i] = state->x[OD1][i];
			xp[CG][i] = state->x[CG][i];
			xp[HB2][i] = state->x[HB2][i];
			xp[HB3][i] = state->x[HB3][i];
		}
	}

	close_trx(tstatus);

	return 0;
}

static double sumprod(rvec a, rvec b)
{
	double val = 0;
	for(int i = 0; i < 3; i++)
	{
		val = val + a[i]*b[i];
	}
	return val;
}

static double my_dist(rvec **xp, int x, int y)
{
	double dist = 0;

	for(int j = 0; j < 3; j++)
	{
		dist = dist + ((*xp)[x][j]-(*xp)[y][j])*((*xp)[x][j]-(*xp)[y][j]);
	}
	dist = sqrt(dist);
	return dist;
}

static void rotate_atom(int atom, t_state *state, double dtheta, int x, int y, rvec *xp)
{
	rvec k,z,l,m,p,yy;
	double t,a1,a2,b1,spl,spk,spyk,c1,cdtheta;

	for(int i = 0; i < 3; i++)
	{
		yy[i] = xp[y][i];
		k[i] = yy[i] - xp[x][i];
	}
	spk = sumprod(k,k);
	spyk = sumprod(yy,k);
	c1 = sin(dtheta)/sqrt(spk);
	cdtheta = cos(dtheta);
	for(int i = 0; i < 3; i++)
	{
		p[i] = xp[atom][i];
	}
	t = (sumprod(p,k) - spyk)/spk;

	for(int i = 0; i < 3; i++)
	{
		z[i] = yy[i] + t*k[i];
		l[i] = p[i] - z[i];
	}
	spl = sumprod(l,l);
	if(spl != 0)
	{
		a1 = spl*c1;
		a2 = a1*k[1];
		a1 = a1*k[0];
		b1 = spl*cdtheta;
		m[2] = (a1*l[1]+b1*l[2]-a2*l[0])/spl;
		m[1] = (l[1]*m[2]-a1)/l[2];
		m[0] = (a2+l[0]*m[2])/l[2];
		for(int i = 0; i < 3; i++)
		{
			state->x[atom][i] = m[i]+z[i];
		}
	}
}


static void rotate_glycan(t_blocka *grps, int *isize, int gnum, t_state *state, 
                          double dtheta, int x, int y, rvec *xp)
{
	rvec k,z,l,m,p,yy;
	double t,a1,a2,b1,spl,spk,spyk,c1,cdtheta;

	for(int i = 0; i < 3; i++)
	{
		yy[i] = xp[y][i];
		k[i] = yy[i] - xp[x][i];
	}
	spk = sumprod(k,k);
	spyk = sumprod(yy,k);
	c1 = sin(dtheta)/sqrt(spk);
	cdtheta = cos(dtheta);
	for(int j = 0; j < isize[gnum]; j++)
	{
		for(int i = 0; i < 3; i++)
		{
			p[i] = xp[my_ind(grps,gnum,j)][i];
		}
		t = (sumprod(p,k) - spyk)/spk;

		for(int i = 0; i < 3; i++)
		{
			z[i] = yy[i] + t*k[i];
			l[i] = p[i] - z[i];
		}
		spl = sumprod(l,l);
		if(spl != 0)
		{
			a1 = spl*c1;
			a2 = a1*k[1];
			a1 = a1*k[0];
			b1 = spl*cdtheta;
			m[2] = (a1*l[1]+b1*l[2]-a2*l[0])/spl;
			m[1] = (l[1]*m[2]-a1)/l[2];
			m[0] = (a2+l[0]*m[2])/l[2];
			for(int i = 0; i < 3; i++)
			{
				state->x[my_ind(grps,gnum,j)][i] = m[i]+z[i];
			}
		}
	}
}

static atom_id my_ind(t_blocka *grps, int i, int j)
{
	return grps->a[grps->index[i]+j];
}

static int sudorun(int argc, char *argv[])
{
	GlyRotModule my_mod;
	return my_mod.run(argc, argv);
}

int main(int argc,char *argv[])
{
	GlyRotModule my_mod;

	return runCommandLineModule(argc, argv, &my_mod);
}

