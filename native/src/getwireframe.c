#include "common.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "vmath.h"
#include "raytrace.h"
#include "bu.h"

struct wc_solid_geometry
{
  struct bu_list l;
  char *solid_name;
  int dashflag;
  int color[3];
  struct bu_list vhead;
};

struct wc_display
{
  struct wc_solid_geometry *solids;
  char title[1024];
  int use_fixed_color;
  int allow_dash;
  int number_solids;
  int fixed_color[3];
};

int
rgb_int(int r, int green, int blue)
{
  return 0x00ff00;
}

union tree *
leaf_func(struct db_tree_state *tsp, const struct db_full_path *path, struct rt_db_internal *ip, genptr_t data)
{
  union tree *curtree;
  struct wc_display *display = (struct wc_display *) data;
  struct wc_solid_geometry *solid;

  RT_CK_DB_INTERNAL(ip);
  RT_CK_TESS_TOL(tsp->ts_ttol);
  BN_CK_TOL(tsp->ts_tol);
  RT_CK_RESOURCE(tsp->ts_resp);

  for (BU_LIST_FOR(solid, wc_solid_geometry, (struct bu_list * )display->solids))
  {
    if (strcmp(solid->solid_name, path->fp_names[path->fp_len - 1]->d_namep) == 0)
    {
      RT_GET_TREE(curtree, tsp->ts_resp);
      curtree->tr_op = OP_NOP;

      return curtree;
    }
  }

  display->number_solids += 1;
  BU_GET(solid, struct wc_solid_geometry);

  BU_LIST_INIT(&solid->vhead);

  if (!ip->idb_meth->ft_plot || ip->idb_meth->ft_plot(&solid->vhead, ip, tsp->ts_ttol, tsp->ts_tol, NULL ) < 0)
  {
    bu_log("%s: plot failure\n", DB_FULL_PATH_CUR_DIR(path)->d_namep);
    return TREE_NULL; /* ERROR */
  }

  solid->solid_name = path->fp_names[path->fp_len - 1]->d_namep;
  if (display->allow_dash)
  {
    solid->dashflag = (tsp->ts_sofar & (TS_SOFAR_MINUS | TS_SOFAR_INTER));
  } else
  {
    solid->dashflag = 0;
  }
  if (display->use_fixed_color)
  {
    solid->color[0] = display->fixed_color[0];
    solid->color[1] = display->fixed_color[1];
    solid->color[2] = display->fixed_color[2];
  } else
  {
    solid->color[0] = (int) (tsp->ts_mater.ma_color[0] * 255.);
    solid->color[1] = (int) (tsp->ts_mater.ma_color[1] * 255.);
    solid->color[2] = (int) (tsp->ts_mater.ma_color[2] * 255.);
  }

  BU_LIST_PUSH(&(display->solids->l), &(solid->l));

  /* Indicate success by returning something other than TREE_NULL */
  RT_GET_TREE(curtree, tsp->ts_resp);
  curtree->tr_op = OP_NOP;

  return curtree;
}

int main(int argc, char **argv)
{
  static struct rt_i *rtip;
  char *av[2];
  struct model *the_model;
  struct rt_tess_tol ttol;
  struct bn_tol tol;
  struct db_tree_state tree_state;
  struct wc_display display;
  struct wc_solid_geometry *solid;
  struct bn_vlist *vp;
  int i;
  int start_line;
  int position;
  int current_count;
  float min[3], max[3];

  min[0] = min[1] = min[2] = INFINITY;
  max[0] = max[1] = max[2] = -INFINITY;

  if (argc < 3)
  {
    bu_exit(1, "Usage: %s model.g objects...\n", argv[0]);
  }

  display.use_fixed_color = 0;
  display.fixed_color[0] = 255;
  display.fixed_color[1] = 0;
  display.fixed_color[2] = 0;
  display.allow_dash = 1;
  display.number_solids = 0;

  BU_GET(display.solids, struct wc_solid_geometry);
  BU_LIST_INIT(&(display.solids->l));

  strcpy(display.title, "");
  rtip = rt_dirbuild(argv[1], display.title, sizeof(display.title));
  if (rtip == RTI_NULL )
  {
    bu_exit(2, "Building the database directory for [%s] FAILED\n", argv[1]);
  }

  while (argc > 2)
  {
    printf("load (%s)\n", argv[2]);
    if (rt_gettree(rtip, argv[2]) < 0)
      printf("Loading the geometry for [%s] FAILED\n", argv[2]);
    av[0] = argv[2];
    argc--;
    argv++;
  }

  tree_state = rt_initial_tree_state;
  tree_state.ts_tol = &tol;
  tree_state.ts_ttol = &ttol;
  tree_state.ts_m = &the_model;

  ttol.magic = RT_TESS_TOL_MAGIC;
  ttol.abs = 0.0;
  ttol.rel = 0.01;
  ttol.norm = 0.0;
  tol.magic = BN_TOL_MAGIC;
  tol.dist = 0.0005;
  tol.dist_sq = tol.dist * tol.dist;
  tol.perp = 1e-6;
  tol.para = 1 - tol.perp;

  /* make empty NMG model */
  the_model = nmg_mm();

  db_walk_tree(rtip->rti_dbip, 1, (const char **) av, 1, &tree_state, NULL, NULL, leaf_func, &display);

  if (display.title[0])
  {
    printf("Title: %s\n", display.title);
  }

  printf("number of solids (%d)\n", display.number_solids);

  for (BU_LIST_FOR(solid, wc_solid_geometry, (struct bu_list * )display.solids))
  {
    start_line = 0;
    for (BU_LIST_FOR(vp, bn_vlist, &solid->vhead))
    {
      for (i = 0; i < vp->nused; i++)
      {
	if (vp->pt[i][0] < min[0]) min[0] = vp->pt[i][0];
	if (vp->pt[i][1] < min[1]) min[1] = vp->pt[i][1];
	if (vp->pt[i][2] < min[2]) min[2] = vp->pt[i][2];
	if (vp->pt[i][0] > max[0]) max[0] = vp->pt[i][0];
	if (vp->pt[i][1] > max[1]) max[1] = vp->pt[i][1];
	if (vp->pt[i][2] > max[2]) max[2] = vp->pt[i][2];
	printf("%f, %f, %f,", vp->pt[i][0], vp->pt[i][1], vp->pt[i][2]);
	start_line += 1;
	if (start_line > 3)
	{
	  printf("\n");
	  start_line = 0;
	}
      }
    }
    printf("\n");
  }
  printf("\n");
  position = 0;
  for (BU_LIST_FOR(solid, wc_solid_geometry, (struct bu_list * )display.solids))
  {
    current_count = 0;
    for (BU_LIST_FOR(vp, bn_vlist, &solid->vhead))
    {
      for (i = 0; i < vp->nused; i++)
      {
	if (vp->cmd[i] == 0)
	{
	  if (current_count > 0)
	  {
	    printf(", %d],\n", current_count);
	    current_count = 0;
	  }
	  printf("  [%d", position);
	}
	position += 1;
	current_count += 1;
      }
    }
    printf(", %d],\n", current_count);
  }
  for (BU_LIST_FOR(solid, wc_solid_geometry, (struct bu_list * )display.solids))
  {
    start_line = 0;
    for (BU_LIST_FOR(vp, bn_vlist, &solid->vhead))
    {
      for (i = 0; i < vp->nused; i++)
      {
	printf("%f, %f, %f, 1.0, ", (solid->color[0]/255.0), (solid->color[1]/255.0), (solid->color[2]/255.0));
	start_line += 1;
	if (start_line > 3)
	{
	  printf("\n");
	  start_line = 0;
	}
      }
    }
    printf("\n");
  }
  printf("\n");

  printf("min (%f, %f, %f)\n", min[0], min[1], min[2]);
  printf("max (%f, %f, %f)\n", max[0], max[1], max[2]);
  printf("diff (%f, %f, %f)\n", max[0]-min[0], max[1]-min[1], max[2]-min[2]);
  printf("center (%f, %f, %f)\n", min[0]+((max[0]-min[0])/2), min[1]+((max[1]-min[1])/2), min[2]+((max[2]-min[2])/2));
  return 0;
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
