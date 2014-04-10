/*
    This file is part of darktable,
    copyright (c) 2009--2011 johannes hanika.
    copyright (c) 2014 Stéphane Gimenez

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "bauhaus/bauhaus.h"
#include "common/darktable.h"
#include "common/interpolation.h"
#include "develop/imageop.h"
#include "dtgtk/slider.h"
#include "gui/gtk.h"
#include <gtk/gtk.h>
#include <stdlib.h>

#define CACORRECT_DEBUG 0

// this is the version of the modules parameters,
// and includes version information about compile-time dt
DT_MODULE_INTROSPECTION(1, dt_iop_cacorrect_params_t)

typedef struct dt_iop_cacorrect_params_t
{
  int keep;
}
dt_iop_cacorrect_params_t;

typedef struct dt_iop_cacorrect_gui_data_t
{
}
dt_iop_cacorrect_gui_data_t;

typedef struct dt_iop_cacorrect_global_data_t
{
}
dt_iop_cacorrect_global_data_t;

typedef struct dt_iop_cacorrect_data_t
{
}
dt_iop_cacorrect_data_t;

// this returns a translatable name
const char *name()
{
  // make sure you put all your translatable strings into _() !
  return _("chromatic aberrations");
}

int
groups ()
{
  return IOP_GROUP_CORRECT;
}

int flags ()
{
  return IOP_FLAGS_ONE_INSTANCE;
}

static int get_legacy()
{
  return dt_conf_get_int("plugins/darkroom/cacorrect/legacy");
}

static int get_displayed()
{
  return dt_conf_get_int("plugins/darkroom/cacorrect/display");
}

static int get_visuals()
{
  return dt_conf_get_int("plugins/darkroom/cacorrect/visuals");
}

int
output_bpp(dt_iop_module_t *module, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  return sizeof(float);
}

/** modify regions of interest (optional, per pixel ops don't need this) */
// void modify_roi_out(struct dt_iop_module_t *self, struct dt_dev_pixelpipe_iop_t *piece, dt_iop_roi_t *roi_out, const dt_iop_roi_t *roi_in);

void modify_roi_in(struct dt_iop_module_t *self, struct dt_dev_pixelpipe_iop_t *piece, const dt_iop_roi_t *roi_out, dt_iop_roi_t *roi_in)
{
/*
  roi_in->width  = roi_out->width * 1.5;
  roi_in->height = roi_out->height * 1.5;
  roi_in->x = roi_out->x - roi_out->width * 0.25;
  roi_in->y = roi_out->y - roi_out->width * 0.25;
  if (roi_in->x < 0) { roi_in->width -= roi_in->x; roi_in->x = 0; }
  if (roi_in->y < 0) { roi_in->height -= roi_in->y; roi_in->y = 0; }
  roi_in->width = MIN(piece->pipe->image.width, roi_in->width);
  roi_in->height = MIN(piece->pipe->image.height, roi_in->height);
*/
  roi_in->x = 0;
  roi_in->y = 0;
  roi_in->width = piece->pipe->image.width;
  roi_in->height = piece->pipe->image.height;
}

static int
FC(const int row, const int col, const unsigned int filters)
{
  return filters >> (((row << 1 & 14) + (col & 1)) << 1) & 3;
}

/*==============================================================================
 * begin raw therapee code, hg checkout of june 05, 2013 branch master.
 *============================================================================*/

////////////////////////////////////////////////////////////////
//
//		Chromatic Aberration Auto-correction
//
//		copyright (c) 2008-2010  Emil Martinec <ejmartin@uchicago.edu>
//
//
// code dated: November 26, 2010
//
//	CA_correct_RT.cc is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
static int
LinEqSolve(int nDim, double* pfMatr, double* pfVect, double* pfSolution)
{
//==============================================================================
// return 1 if system not solving, 0 if system solved
// nDim - system dimension
// pfMatr - matrix with coefficients
// pfVect - vector with free members
// pfSolution - vector with system solution
// pfMatr becames trianglular after function call
// pfVect changes after function call
//
// Developer: Henry Guennadi Levkin
//
//==============================================================================

  double fMaxElem;
  double fAcc;

  int i, j, k, m;

  for(k=0; k<(nDim-1); k++)  // base row of matrix
  {
    // search of line with max element
    fMaxElem = fabs( pfMatr[k*nDim + k] );
    m = k;
    for (i=k+1; i<nDim; i++)
    {
      if(fMaxElem < fabs(pfMatr[i*nDim + k]) )
      {
        fMaxElem = pfMatr[i*nDim + k];
        m = i;
      }
    }

    // permutation of base line (index k) and max element line(index m)
    if(m != k)
    {
      for(i=k; i<nDim; i++)
      {
        fAcc               = pfMatr[k*nDim + i];
        pfMatr[k*nDim + i] = pfMatr[m*nDim + i];
        pfMatr[m*nDim + i] = fAcc;
      }
      fAcc = pfVect[k];
      pfVect[k] = pfVect[m];
      pfVect[m] = fAcc;
    }

    if( pfMatr[k*nDim + k] == 0.)
    {
      //linear system has no solution
      return 1; // needs improvement !!!
    }

    // triangulation of matrix with coefficients
    for(j=(k+1); j<nDim; j++)  // current row of matrix
    {
      fAcc = - pfMatr[j*nDim + k] / pfMatr[k*nDim + k];
      for(i=k; i<nDim; i++)
      {
        pfMatr[j*nDim + i] = pfMatr[j*nDim + i] + fAcc*pfMatr[k*nDim + i];
      }
      pfVect[j] = pfVect[j] + fAcc*pfVect[k]; // free member recalculation
    }
  }

  for(k=(nDim-1); k>=0; k--)
  {
    pfSolution[k] = pfVect[k];
    for(i=(k+1); i<nDim; i++)
    {
      pfSolution[k] -= (pfMatr[k*nDim + i]*pfSolution[i]);
    }
    pfSolution[k] = pfSolution[k] / pfMatr[k*nDim + k];
  }

  return 0;
}
//end of linear equation solver
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "caold.c"

/*==============================================================================
 * end raw therapee code
 *============================================================================*/

//todo: can we cache the analyse phase somehow?

static void
CA_correct(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const float *const in, float *out, const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out)
{
  const int width  = roi_in->width;
  const int height = roi_in->height;

#if CACORRECT_DEBUG
  printf("cacorrect: i.w=%d i.h=%d i.x=%d i.y=%d\n",
         roi_in->width, roi_in->height, roi_in->x, roi_in->y);
  printf("cacorrect: o.w=%d o.h=%d o.x=%d o.y=%d\n",
         roi_out->width, roi_out->height, roi_out->x, roi_out->y);
#endif

  const uint32_t filters = dt_image_flipped_filter(&piece->pipe->image);

  const int sl = 8; // max allowed CA shift, must be even
  const int srad = 22; // size of samples, must be even
  const int subs = 22; // spacing between samples, must be even
  const int deg = 4; // order of 2-variables polynomial fit

  const int degnum = ((deg + 1) * (deg + 2)) / 2;

  const double gamma_analyse = 2.2; // experimental gammas
  const double gamma_shift = 2.2;

  const int TS = (width > 2024 && height > 2024) ? 256 : 128;
  const int border = 2 + srad + sl;
  const int ncoeff = 2;

  int visuals = get_visuals();
  if (visuals == 1 && piece->pipe->type != DT_DEV_PIXELPIPE_EXPORT)
    visuals = 0;

  // polynomial fit coefficients
  double fitmat[2][ncoeff][degnum*degnum];
  double fitvec[2][ncoeff][degnum];
  double fitsol[2][ncoeff][degnum];

  // these small sizes aren't safe to run through the code below.
  if (height < 3 * border || width < 3 * border) return;

  const struct dt_interpolation *ip_trans =
    dt_interpolation_new(DT_INTERPOLATION_USERPREF);

  // initialize fit arrays
  for (int color = 0; color < 2; color++)
    for (int coeff = 0; coeff < ncoeff; coeff++)
    {
      for (int i = 0; i < degnum*degnum; i++)
        fitmat[color][coeff][i] = 0;
      for (int i = 0; i < degnum; i++)
        fitvec[color][coeff][i] = 0;
    }

#if CACORRECT_DEBUG
  printf("pre-analysis\n");
#endif

  // pre-analysis (todo: allow a second pass to ditch out the outliers)

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int top = 0; top < height; top += TS - 2*border)
    for (int left = 0; left < width; left += TS - 2*border)
    {
      int rm = top + TS > height ? height - top : TS;
      int cm = left + TS > width ? width - left : TS;

      float rgb[TS][TS];
      float la[TS][TS][2];

      for (int r = 0; r < rm; r++)
        for (int c = 0; c < cm; c++)
        {
          int inr = top + r;
          int inc = left + c;
          rgb[r][c] = in[inr*width+inc];
          rgb[r][c] = pow(rgb[r][c], 1/gamma_analyse);
        }

      // compute derivatives
      for (int r = 2; r < rm - 2; r++)
        for (int c = 2; c < cm - 2; c++)
          {
            la[r][c][0] = (rgb[r+2][c] - rgb[r-2][c])/4;
            la[r][c][1] = (rgb[r][c+2] - rgb[r][c-2])/4;
          }

      // thread-local components to be summed up
      double tfitmat[2][ncoeff][degnum*degnum];
      double tfitvec[2][ncoeff][degnum];

      for (int color = 0; color < 2; color++)
        for (int coeff = 0; coeff < ncoeff; coeff++)
        {
          for (int i = 0; i < degnum*degnum; i++)
            tfitmat[color][coeff][i] = 0;
          for (int i = 0; i < degnum; i++)
            tfitvec[color][coeff][i] = 0;
        }

      // compute centered offsets for sampling areas
      int rsuboff = (height/4 - border/2) % subs < subs/2 ? 0 : subs/2;
      int csuboff = (width/4 - border/2) % subs < subs/2 ? 0 : subs/2;
      int roff = (height/4 + rsuboff - (top + border)/2) % subs;
      int coff = (width/4 + csuboff - (left + border)/2) % subs;
      if (roff < 0) roff += subs;
      if (coff < 0) coff += subs;

      for (int r = border + 2*roff; r < rm - border; r += 2*subs)
        for (int c = border + 2*coff; c < cm - border; c += 2*subs)
        {
          int inr = top + r;
          int inc = left + c;
          if (inr < border || inr > height - 1 - border ||
              inc < border || inc > width - 1 - border)
            continue;

          double y = 2.0 * (2.0 * (double)inr/height - 1);
          double x = 2.0 * (2.0 * (double)inc/width - 1);

          for (int color = 0; color < 2; color++)
          {

            // compute variation map
            double t[2*sl+1][2*sl+1];
            for (int i = 0; i < 2*sl+1; i++)
              for (int j = 0; j < 2*sl+1; j++)
                t[i][j] = 0.;

            for (int gi = -srad; gi < srad; gi++)
              for (int gj = -srad; gj < srad; gj++)
              {
                if (FC(inr+gi, inc+gj, filters) != 1) continue;
                double dgv = la[r+gi][c+gj][0];
                double dgh = la[r+gi][c+gj][1];
                for (int di = -sl; di <= sl; di++)
                  for (int dj = -sl; dj <= sl; dj++)
                    if (FC(inr+gi+di, inc+gj+dj, filters) == 2*color)
                      if (di*di + dj*dj < (int)((sl+0.5)*(sl+0.5))) // mask
                      {
                        double dcv = la[r+gi+di][c+gj+dj][0];
                        double dch = la[r+gi+di][c+gj+dj][1];
                        double kv = dcv - dgv;
                        double kh = dch - dgh;
                        t[di+sl][dj+sl] += kv * kv + kh * kh;
                      }
              }

            // compute averages and locate indexed minimum
            int I = 0, J = 0;
            double min = INFINITY;
            for (int di = -sl+1; di <= sl-1; di++)
              for (int dj = -sl+1+(~di&1); dj <= sl-1; dj += 2)
                if (di*di + dj*dj < (int)((sl-0.5)*(sl-0.5))) // mask
                {
                  double sum = t[di+sl+1][dj+sl] + t[di+sl-1][dj+sl]
                    + t[di+sl][dj+sl+1] + t[di+sl][dj+sl-1];
                  t[di+sl][dj+sl] = sum/4;
                  if (di*di + dj*dj < (int)((sl-2.5)*(sl-2.5))) // mask
                    if (sum < min) { min = sum; I = di; J = dj; }
                }

            // adjustments
            double bc[2] = { 0 };
            double bw[2] = { 0 };
            {
              int i = I + sl;
              int j = J + sl;
              double dv = (t[i+2][j] - t[i-2][j])/4;
              double dh = (t[i][j+2] - t[i][j-2])/4;
              double d2v = (t[i+2][j] + t[i-2][j] - 2*t[i][j])/8;
              double d2h = (t[i][j+2] + t[i][j-2] - 2*t[i][j])/8;
              double d2m =
                (t[i+1][j+1] - t[i-1][j+1] - t[i+1][j-1] + t[i-1][j-1])/4;
              double div = 4*d2v*d2h - d2m*d2m;
              if (div > 0)
              {
                double sqrtD = sqrt((d2v-d2h)*(d2v-d2h) + d2m*d2m);
                /* double sin2th = d2m / sqrtD; */
                double cos2th = (d2v - d2h) / sqrtD;
                double sqsinth = (1 - cos2th)/2;
                double sqcosth = (1 + cos2th)/2;
                double lb = (d2v+d2h + sqrtD)/2;
                double ls = (d2v+d2h - sqrtD)/2;
                double pv = -(2*dv*d2h - dh*d2m) / div;
                double ph = -(2*dh*d2v - dv*d2m) / div;
                pv = isnan(pv) ? 0 : CLAMP(pv, -2.5, +2.5);
                ph = isnan(ph) ? 0 : CLAMP(ph, -2.5, +2.5);
                bc[0] = I + pv;
                bc[1] = J + ph;
                // weights (deviations are in 1/l^2)
                bw[0] = 1/(sqrt(sqcosth/(lb*lb) + sqsinth/(ls*ls)));
                bw[1] = 1/(sqrt(sqsinth/(lb*lb) + sqcosth/(ls*ls)));
              }
            }

            // show weights
            if (visuals)
            {
              int outr = inr + roi_in->y - roi_out->y;
              int outc = inc + roi_in->x - roi_out->x;
              int off = FC(inr, inc, filters)&1;
              int outi = outr*roi_out->width+outc;
              if (outr >= 4 && outr < roi_out->height - 4
                  && outc >= 4 && outc < roi_out->width - 4)
              {
                double s = 25*sqrt(bw[0]*bw[0]+bw[1]*bw[1]);
                for (int di = -6; di <= 6; di++)
                  for (int dj = -6; dj <= 6; dj++)
                    if (FC(inr+di, inc+dj+off, filters) == 1)
                      out[outi+di*roi_out->width+dj+off] +=
                        MIN(1,s*exp(1.5*(-di*di-dj*dj)));
              }
            }

#if CACORRECT_DEBUG
            if (y > -0.1 && y < 0.1 && x > -0.1 && x < 0.1)
            {
              printf("x=%+.4lf y=%+.4lf g=%.2lf ", x, y, rgb[r][c]);
              printf("%+.5f (%.5f) %+.5f (%.5f)\n", bc[0], bw[0], bc[1], bw[1]);
            }
#endif

            double piy = 1.; // pow_i(y)
            int iobs = 0;
            for (int i = 0; i <= deg; i++, piy *= y)
            {
              double pjx = 1.; // pow_j(x)
              for (int j = 0; j <= deg - i; j++, pjx *= x)
              {
                int ideg = 0;
                double pmy = 1.; // pow_m(y)
                double o = piy*pjx;
                for (int m = 0; m <= deg; m++, pmy *= y)
                {
                  double pnx = 1.; // pow_n(x)
                  for (int n = 0; n <= deg - m; n++, pnx *= x)
                  {
                    for (int coeff = 0; coeff < ncoeff; coeff++)
                    {
                      double v = o*bw[coeff]*pmy*pnx;
                      tfitmat[color][coeff][iobs*degnum+ideg] += v;
                    }
                    ideg++;
                  }
                }
                for (int coeff = 0; coeff < ncoeff; coeff++)
                {
                  double v = o*bw[coeff]*bc[coeff];
                  tfitvec[color][coeff][iobs] += v;
                }
                iobs++;
              }
            }
          }
        }

#ifdef _OPENMP
#pragma omp critical
#endif
      {
        for (int color = 0; color < 2; color++)
          for (int coeff = 0; coeff < ncoeff; coeff++)
          {
            for (int i = 0; i < degnum*degnum; i++)
              fitmat[color][coeff][i] += tfitmat[color][coeff][i];
            for (int i = 0; i < degnum; i++)
              fitvec[color][coeff][i] += tfitvec[color][coeff][i];
          }
      }
    }

  // end of pre-analysis

#if CACORRECT_DEBUG
  printf("solve\n");
#endif

  // fit parameters to blockshifts
  for (int color = 0; color < 2; color++)
    for (int coeff = 0; coeff < ncoeff; coeff++)
    {
      int res =
        LinEqSolve(degnum, fitmat[color][coeff], fitvec[color][coeff],
                   fitsol[color][coeff]);
      if (res)
      {
        printf("Automatic color correction pass failed.\n");
        printf("Can't solve linear equations, color=%d coeff=%d.\n",
               color, coeff);
        return;
      }
    }

  // shifting red and blue layers

#if CACORRECT_DEBUG
  printf("shift\n");
#endif

  const int margin = sl + 8; // 8 for interpolation

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int top = -margin; top < height - margin; top += TS - 2*margin)
    for (int left = -margin; left < width - margin; left += TS - 2*margin)
    {
      int rm = top + TS > height + margin ? height + margin - top : TS;
      int cm = left + TS > width + margin ? width + margin - left : TS;

      float rgb[TS][TS];

      for (int r = 0; r < rm; r++)
        for (int c = 0; c < cm; c++)
          {
            int inr = top + r;
            int inc = left + c;
            if (FC(inr, inc, filters) == 1) continue;
            if (inr < 0 || inr >= height || inc < 0 || inc >= width)
            {
              rgb[r][c] = 0./0.; // nan values will be extrapolated
            }
            else
            {
              rgb[r][c] = in[inr*width+inc];
              rgb[r][c] = pow(rgb[r][c], 1/gamma_shift);
            }
          }

      for (int r = margin; r < rm - margin; r++)
        for (int c = margin; c < cm - margin; c++)
        {
          int inr = top + r;
          int inc = left + c;
          if (FC(inr, inc, filters) == 1) continue;
          int color = FC(inr, inc, filters)/2;

          int outr = inr + roi_in->y - roi_out->y;
          int outc = inc + roi_in->x - roi_out->x;
          if (outr < 0 || outr >= roi_out->height
              || outc < 0 || outc >= roi_out->width)
            continue;

          double bc[2] = { 0 };
          double y = 2.0 * (2.0 * (double)inr/height - 1);
          double x = 2.0 * (2.0 * (double)inc/width - 1);

          // compute polynomial
          int ideg = 0;
          double pmy = 1.; // pow_m(y)
          for (int m = 0; m <= deg; m++, pmy *= y)
          {
            double pnx = 1.; // pow_n(x)
            for (int n = 0; n <= deg - m; n++, pnx *= x)
            {
              for (int coeff = 0; coeff < ncoeff; coeff++)
              {
                double o = pmy*pnx;
                bc[coeff] += o*fitsol[color][coeff][ideg];
              }
              ideg++;
            }
          }

          double bv = CLAMPS(bc[0], -sl, sl);
          double bh = CLAMPS(bc[1], -sl, sl);

          // shift by interpolation
          int xmin = MAX(-inc, -margin);
          int ymin = MAX(-inr, -margin);
          int xmax = MIN(width-1-inc, +margin);
          int ymax = MIN(height-1-inr, +margin);
          double v =
            dt_interpolation_compute_extrapolated_sample(
              ip_trans, &rgb[r][c], bh/2, bv/2,
              xmin/2, ymin/2, xmax/2, ymax/2, 2, 2*TS);

          v = MAX(0., v);

          // show shift lines
          if (visuals)
          {
            double norm = sqrt(bv*bv + bh*bh);
            for (double l = 0; l <= sl; l += 0.5)
            {
              double d = norm - l;
              if (fabs(d) < 0.01) {
                v *= d >= 0 ? 1.75 : 1.5;
                v += d > 0 ? 0.2 : 0.075;
              }
            }
          }

          v = pow(v, gamma_shift);

          out[outr*roi_out->width+outc] = v;

        }
    }

#if CACORRECT_DEBUG
  printf("end\n");
#endif

}

/** process, all real work is done here. */
void process (struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, void *i, void *o, const dt_iop_roi_t *roi_in, const dt_iop_roi_t *roi_out)
{
  float *in = (float *)i;
  float *out = (float *)o;
  for (int i = 0; i < roi_out->height; i++)
    for (int j = 0; j < roi_out->width; j++)
      out[i*roi_out->width+j] =
        in[(i+roi_out->y-roi_in->y)*roi_in->width+(j+roi_out->x-roi_in->x)];
  int type = get_legacy();
  if (type == 0)
    CA_correct(self, piece, in, out, roi_in, roi_out);
  else if (type == 1)
    CA_correct_old(self, piece, in, out, roi_in, roi_out);
}

void reload_defaults(dt_iop_module_t *module)
{
  // init defaults:
  dt_iop_cacorrect_params_t tmp = { 50 };
  memcpy(module->params, &tmp, sizeof(dt_iop_cacorrect_params_t));
  memcpy(module->default_params, &tmp, sizeof(dt_iop_cacorrect_params_t));

  // can't be switched on for non-raw images:
  if(dt_image_is_raw(&module->dev->image_storage))
    module->hide_enable_button = 0;
  else
    module->hide_enable_button = 1;
  module->default_enabled = 0;
}

/** init, cleanup, commit to pipeline */
void init(dt_iop_module_t *module)
{
  module->params = malloc(sizeof(dt_iop_cacorrect_params_t));
  module->default_params = malloc(sizeof(dt_iop_cacorrect_params_t));
  module->data = NULL;
  module->gui_data = NULL;
  // our module is disabled by default
  module->default_enabled = 0;

  // we come just before demosaicing.
  module->priority = 70; // module order created by iop_dependencies.py, do not edit!
  module->params_size = sizeof(dt_iop_cacorrect_params_t);
}

void cleanup(dt_iop_module_t *module)
{
  free(module->params);
  module->params = NULL;
  module->gui_data = NULL;
  module->data = NULL;
}

/** commit is the synch point between core and gui, so it copies params to pipe data. */
void commit_params(struct dt_iop_module_t *self, dt_iop_params_t *params, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  // preview pipe doesn't have mosaiced data either:
  if (!(pipe->image.flags & DT_IMAGE_RAW)
     || dt_dev_pixelpipe_uses_downsampled_input(pipe))
    piece->enabled = 0;
  if (piece->pipe->type == DT_DEV_PIXELPIPE_THUMBNAIL)
    piece->enabled = 0;
  if (get_displayed() == 0 && piece->pipe->type != DT_DEV_PIXELPIPE_EXPORT)
    piece->enabled = 0;

}

void init_pipe(struct dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  piece->data = malloc(sizeof(dt_iop_cacorrect_data_t));
  self->commit_params(self, self->default_params, pipe, piece);
}

void cleanup_pipe(struct dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  free(piece->data);
}

void gui_update(dt_iop_module_t *self)
{
}

void gui_init(dt_iop_module_t *self)
{
  self->gui_data = malloc(sizeof(dt_iop_cacorrect_gui_data_t));
  self->widget = gtk_vbox_new(TRUE, DT_BAUHAUS_SPACE);
}

void gui_cleanup(dt_iop_module_t *self)
{
  free(self->gui_data);
  self->gui_data = NULL;
}

// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.sh
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-space on;
