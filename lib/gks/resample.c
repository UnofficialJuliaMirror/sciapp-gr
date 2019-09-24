#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "gkscore.h"
#include "gks.h"

#ifndef INFINITY
#define INFINITY (1.0 / 0.0)
#endif
#ifndef M_PI
#define M_PI (3.141592653589793)
#endif

static int resample_method = GKS_K_RESAMPLE_DEFAULT;

/*!
 * Set the resample method for up and downscaling. The default method is nearest neighbour.
 *
 * \param[in] flag Resample method, valid options are GKS_K_RESAMPLE_DEFAULT, GKS_K_RESAMPLE_LINEAR,
 * GKS_K_RESAMPLE_NEAREST and GKS_K_RESAMPLE_LANCZOS
 */
void gks_set_resample_method(int flag)
{
  if (flag == GKS_K_RESAMPLE_NEAREST || flag == GKS_K_RESAMPLE_LANCZOS || flag == GKS_K_RESAMPLE_LINEAR)
    {
      resample_method = flag;
    }
  else
    {
      fprintf(stderr, "Only GKS_K_RESAMPLE_DEFAULT, GKS_K_RESAMPLE_LINEAR, GKS_K_RESAMPLE_NEAREST and "
                      "GKS_K_RESAMPLE_LANCZOS are valid values for the flag parameter\n");
    }
}

/*!
 * Inquire the resample flag status.
 *
 * \returns Resample flag
 */
void gks_inq_resample_method(int *flag)
{
  flag[0] = resample_method;
}

static double *linspace(double a, double b, int n, double u[])
{
  double c;
  int i;

  /* make sure points and array numbers are valid */
  if (n < 2 || u == 0) return (void *)0;

  /* step size */
  c = (b - a) / (n - 1);

  /* fill vector */
  for (i = 0; i < n - 1; ++i)
    {
      u[i] = a + i * c;
    }

  /* fix last entry to b */
  u[n - 1] = b;
  /* done */
  return u;
}


static double lanczos(double x, int a)
{
  if (x == 0.0)
    {
      return 1.0;
    }
  else if (-a < x && x < a)
    {
      return (a * sin(M_PI * x) * sin(M_PI * x / a)) / (x * x * M_PI * M_PI);
    }
  return 0.0;
}


static double integrate_box(double left, double right, int width)
{
  if (left > width / 2.0)
    {
      return 0.0;
    }
  if (right < -width / 2.0)
    {
      return 0.0;
    }
  if (left < -width / 2.0)
    {
      left = -width / 2.0;
    }
  if (right > width / 2.0)
    {
      right = width / 2.0;
    }
  return (right - left) / width;
}


static size_t get_index(double x, size_t height, size_t width, double y)
{
  double i_x = x;
  double i_y = y;
  if (x < 0)
    {
      i_x = 0;
    }
  else if (x >= width)
    {
      i_x = width - 1;
    }
  if (y < 0)
    {
      i_y = 0;
    }
  else if (y >= height)
    {
      i_y = height - 1;
    }
  return (size_t)(i_y * width + i_x);
}

/*!
 * Use linear resampling to rescale the data image.
 *
 * \param[in] source_img Data array of the original image
 * \param[out] result_img Data array with resampled values
 * \param[in] width Size of a row in source_img
 * \param[in] height Number of rows in source_img
 * \param[in] new_w Size of a row in result_img
 * \param[in] new_h Number of rows in result_img
 * \param[in] stride Stride of source_img
 * \param[in] w Width of linear filter
 *
 */
static void resample_rgba(const unsigned char *source_img, unsigned char *result_img, size_t width, size_t height,
                          size_t new_w, size_t new_h, size_t stride, int w)
{
  if (w <= 0)
    {
      gks_fatal_error("w greater than 0 required!\n");
    }
  int i, j, h, l;
  double *one_dir_img = (double *)gks_malloc(4 * (int)height * (int)new_w * sizeof(double));
  /* x coordinates before resampling */
  double *x = (double *)gks_malloc(sizeof(double) * new_w);
  double *y = (double *)gks_malloc(sizeof(double) * new_h);
  /* array to store filter values once to increase efficiency */
  double *horizontal_values = (double *)gks_malloc(sizeof(double) * new_w * 2 * w);
  double *vertical_values = (double *)gks_malloc(sizeof(double) * new_h * 2 * w);
  size_t index;

  linspace(0, width - 1, (int)new_w, x);
  linspace(0, height - 1, (int)new_h, y);

  /* resample width */
  for (i = 0; i < height; i++)
    {
      for (j = 0; j < new_w; j++)
        {
          for (h = (int)floor(-w / 2.0 - 0.5 + x[j]) + 1; h < w / 2.0 + 0.5 + x[j]; h++)
            {
              if (i == 0)
                {
                  /* compute linear filter integral values once and store the resulting data */
                  horizontal_values[j * 2 * w + h - (int)floor(-w / 2.0 - 0.5 + x[j]) + 1] =
                      integrate_box(-0.5 + h - x[j], 0.5 + h - x[j], w);
                }
              index = get_index(h, height, stride, i) * 4;
              for (l = 0; l < 4; l++)
                {
                  /* calculate the influence of the image points */
                  one_dir_img[(i * new_w + j) * 4 + l] +=
                      source_img[index + l] * horizontal_values[j * 2 * w + h - (int)floor(-w / 2.0 - 0.5 + x[j]) + 1];
                }
            }
        }
    }

  /* resample height */
  for (i = 0; i < new_w; i++)
    {
      for (j = 0; j < new_h; j++)
        {
          double rgba[4] = {0};
          for (h = (int)floor(-w / 2.0 - 0.5 + y[j]) + 1; h < w / 2.0 + 0.5 + y[j]; h++)
            {
              if (i == 0)
                {
                  /* compute linear filter integral values ones and store the resulting data */
                  vertical_values[j * 2 * w + h - (int)floor(-w / 2.0 - 0.5 + y[j]) + 1] =
                      integrate_box(-0.5 + h - y[j], 0.5 + h - y[j], w);
                }
              index = get_index(i, height, new_w, h) * 4;
              for (l = 0; l < 4; l++)
                {
                  /* calculate the influence of image points */
                  rgba[l] +=
                      one_dir_img[index + l] * vertical_values[j * 2 * w + h - (int)floor(-w / 2.0 - 0.5 + y[j]) + 1];
                }
            }
          result_img[(j * new_w + i) * 4 + 0] = (unsigned char)round(rgba[0]);
          result_img[(j * new_w + i) * 4 + 1] = (unsigned char)round(rgba[1]);
          result_img[(j * new_w + i) * 4 + 2] = (unsigned char)round(rgba[2]);
          result_img[(j * new_w + i) * 4 + 3] = (unsigned char)round(rgba[3]);
        }
    }
  gks_free(horizontal_values);
  gks_free(vertical_values);
  gks_free(x);
  gks_free(y);
  gks_free(one_dir_img);
}

/*!
 * Use lanczos resampling to rescale the data image.
 *
 * \param[in] source_img Data array of the original image
 * \param[out] result_img Data array with resampled values
 * \param[in] width Size of a row in source_img
 * \param[in] height Number of rows in source_img
 * \param[in] new_w Size of a row in result_img
 * \param[in] new_h Number of rows in result_img
 * \param[in] stride Stride of source_img
 * \param[in] a Half width of lanczos filter
 * \param[in] min_val Lower border for colour values
 * \param[in] max_val Upperborder for colour values
 * \param[in] swapx True if x coordinates are swapped in source_img
 * \param[in] swapy True if y coordinates are swapped in source_img
 *
 * Best results will be generated when a = 2 or a = 3.
 */
static void resample_rgba_lanczos(const unsigned char *source_img, unsigned char *result_img, size_t width,
                                  size_t height, size_t new_w, size_t new_h, size_t stride, int a, int min_val,
                                  int max_val, int swapx, int swapy)
{
  if (a <= 0)
    {
      gks_fatal_error("a greater than 0 required!\n");
    }
  int i, j, h, l, hy, jx;
  double *one_dir_img = (double *)gks_malloc(4 * (int)height * (int)new_w * sizeof(double));
  double *horizontal_values = (double *)gks_malloc(sizeof(double) * new_w * a * 2);
  double *vertical_values = (double *)gks_malloc(sizeof(double) * new_h * a * 2);

  /* precompute horizontal lanczos filter values*/
  for (j = 0; j < new_w; j++)
    {
      jx = j;
      if (swapx)
        {
          jx = (int)new_w - 1 - j;
        }
      double destination_position = jx / (double)(new_w - 1) * (double)width - 0.5;
      double sum = 0;
      for (i = 0; i < 2 * a; i++)
        {
          int source_position = (int)floor(destination_position - (a - 1)) + i;
          if (source_position < 0)
            {
              continue;
            }
          if (source_position > width - 1)
            {
              break;
            }
          double lanczos_factor = lanczos(source_position - destination_position, a);
          sum += lanczos_factor;
          horizontal_values[jx * 2 * a + i] = lanczos_factor;
        }
      for (i = 0; i < 2 * a; i++)
        {
          horizontal_values[jx * 2 * a + i] /= sum;
        }
    }

  /* precompute vertical lanczos filter values */
  for (h = 0; h < new_h; h++)
    {
      hy = h;
      if (swapy)
        {
          hy = (int)new_h - 1 - h;
        }
      double destination_position = hy / (double)(new_h - 1) * (double)height - 0.5;
      double sum = 0;
      for (i = 0; i < 2 * a; i++)
        {
          int source_position = (int)floor(destination_position - (a - 1)) + i;
          if (source_position < 0)
            {
              continue;
            }
          if (source_position > height - 1)
            {
              break;
            }
          double lanczos_factor = lanczos(source_position - destination_position, a);
          sum += lanczos_factor;
          vertical_values[hy * 2 * a + i] = lanczos_factor;
        }
      for (i = 0; i < 2 * a; i++)
        {
          vertical_values[hy * 2 * a + i] /= sum;
        }
    }

  /* resample width */
  for (h = 0; h < height; h++)
    {
      hy = h;
      if (swapy)
        {
          hy = (int)height - 1 - h;
        }
      for (j = 0; j < new_w; j++)
        {
          jx = j;
          if (swapx)
            {
              jx = (int)new_w - 1 - j;
            }
          /* linspace between -0.5 and width-0.5 */
          double destination_position = jx / (double)(new_w - 1) * (double)width - 0.5;
          int source_position_offset = (int)floor(destination_position - (a - 1));
          double result[4] = {0};
          /* filter values are not zero */
          for (i = 0; i < 2 * a; i++)
            {
              int source_position = i + source_position_offset;
              if (source_position < 0)
                {
                  continue;
                }
              if (source_position > width - 1)
                {
                  break;
                }
              double lanczos_factor = horizontal_values[jx * 2 * a + i];
              for (l = 0; l < 4; l++)
                {
                  result[l] += lanczos_factor * source_img[(hy * width + source_position) * 4 + l];
                }
            }
          for (l = 0; l < 4; l++)
            {
              one_dir_img[(hy * new_w + jx) * 4 + l] = result[l];
            }
        }
    }

  /* resample height*/
  for (j = 0; j < new_w; j++)
    {
      jx = j;
      if (swapx)
        {
          jx = (int)new_w - 1 - j;
        }
      for (h = 0; h < new_h; h++)
        {
          hy = h;
          if (swapy)
            {
              hy = (int)new_h - 1 - h;
            }
          double destination_position = hy / (double)(new_h - 1) * (double)height - 0.5;
          int source_position_offset = (int)floor(destination_position - (a - 1));
          double result[4] = {0};
          /* where filter is not null */
          for (i = 0; i < 2 * a; i++)
            {
              int source_position = i + source_position_offset;
              if (source_position < 0)
                {
                  continue;
                }
              if (source_position > height - 1)
                {
                  break;
                }
              double lanczos_factor = vertical_values[hy * 2 * a + i];
              for (l = 0; l < 4; l++)
                {
                  result[l] += lanczos_factor * one_dir_img[(source_position * new_w + jx) * 4 + l];
                }
            }
          for (l = 0; l < 4; l++)
            {
              /* clipping */
              if (result[l] > max_val)
                {
                  result[l] = max_val;
                }
              else if (result[l] < min_val)
                {
                  result[l] = min_val;
                }
              result_img[(h * new_w + j) * 4 + l] = (unsigned char)round(result[l]);
            }
        }
    }
  gks_free(horizontal_values);
  gks_free(vertical_values);
  gks_free(one_dir_img);
}

/*!
 * Method to resample the image faster. For this the image is going to be divided into nine regions.
 *
 * \param[in] source_img Data array of the original image
 * \param[out] result_img Data array with resampled values
 * \param[in] width Size of a row in source_img
 * \param[in] height Number of rows in source_img
 * \param[in] new_w Size of a row in result_img
 * \param[in] new_h Number of rows in result_img
 * \param[in] stride Stride of source_img
 * \param[in] w Width of linear filter
 * \param[in] swapx True if x coordinates are swapped in source_img
 * \param[in] swapy True if y coordinates are swapped in source_img
 */
static void resample_rgba_linear(const unsigned char *source_img, unsigned char *result_img, size_t width,
                                 size_t height, size_t new_w, size_t new_h, size_t stride, int w, int swapx, int swapy)
{
  int i, j, ix, jy;
  size_t index;
  /* lower vertical border */
  size_t ver_low = new_h / (2 * height);
  size_t ver_up = new_h * (2 * height - 1) / (2 * height);
  size_t ver_cal = new_h * (2 * height - 2) / (2 * height);
  /* left horizontal border*/
  size_t hori_left = new_w / (2 * width);
  size_t hori_right = new_w * (2 * width - 1) / (2 * width);
  size_t hori_cal = new_w * (2 * width - 2) / (2 * width);
  unsigned char *result = (unsigned char *)gks_malloc((int)(4 * hori_cal * ver_cal));
  unsigned char *left = (unsigned char *)gks_malloc((int)(8 * (width - 1) * (height - 1) * hori_left * ver_low));
  unsigned char *top = (unsigned char *)gks_malloc((int)(8 * (width - 1) * (height - 1) * hori_left * ver_low));
  unsigned char *right = (unsigned char *)gks_malloc((int)(8 * (width - 1) * (height - 1) * hori_left * ver_low));
  unsigned char *bot = (unsigned char *)gks_malloc((int)(8 * (width - 1) * (height - 1) * hori_left * ver_low));

  resample_rgba(source_img, result, width, height, hori_cal, ver_cal, stride, 1);
  resample_rgba(source_img, left, 1, height, hori_left, ver_cal, stride, 1);
  resample_rgba(&source_img[4 * (width - 1)], right, 1, height, hori_left, ver_cal + 1, stride, 1);
  resample_rgba(source_img, top, width, 1, hori_cal, ver_low + 1, stride, 1);
  resample_rgba(&source_img[4 * height * (width - 1)], bot, width, 1, hori_cal, ver_low, stride, 1);
  for (i = 0; i < new_w; i++)
    {
      ix = i;
      if (swapx)
        {
          ix = (int)new_w - 1 - i;
        }
      for (j = 0; j < new_h; j++)
        {
          jy = j;
          if (swapy)
            {
              jy = (int)new_h - 1 - j;
            }
          index = (j * new_w + i) * 4;
          if (ix < hori_left)
            {
              if (jy < ver_low)
                {
                  /* top left corner color like the pixel color in this corner */
                  result_img[index + 0] = source_img[0];
                  result_img[index + 1] = source_img[1];
                  result_img[index + 2] = source_img[2];
                  result_img[index + 3] = source_img[3];
                }
              if (jy >= ver_low)
                {
                  /* bottom left corner color like the pixel color in this corner */
                  result_img[index + 0] = source_img[(height - 1) * 4 * stride + 0];
                  result_img[index + 1] = source_img[(height - 1) * 4 * stride + 1];
                  result_img[index + 2] = source_img[(height - 1) * 4 * stride + 2];
                  result_img[index + 3] = source_img[(height - 1) * 4 * stride + 3];
                }
              if (jy >= ver_low && jy < ver_up)
                {
                  /* left side between top and bottom corner is going to be resampled over the pixels in that line */
                  result_img[index + 0] = left[((jy - ver_low) * hori_left + ix) * 4 + 0];
                  result_img[index + 1] = left[((jy - ver_low) * hori_left + ix) * 4 + 1];
                  result_img[index + 2] = left[((jy - ver_low) * hori_left + ix) * 4 + 2];
                  result_img[index + 3] = left[((jy - ver_low) * hori_left + ix) * 4 + 3];
                }
            }
          if (ix >= hori_right)
            {
              if (jy < ver_low)
                {
                  /* top right corner color like the pixel color in this corner */
                  result_img[index + 0] = source_img[4 * (height - 1) + 0];
                  result_img[index + 1] = source_img[4 * (height - 1) + 1];
                  result_img[index + 2] = source_img[4 * (height - 1) + 2];
                  result_img[index + 3] = source_img[4 * (height - 1) + 3];
                }
              if (jy >= ver_up)
                {
                  /* bottom right corner color like the pixel color in this corner */
                  result_img[index + 0] = source_img[height * 4 * stride - 4];
                  result_img[index + 1] = source_img[height * 4 * stride - 3];
                  result_img[index + 2] = source_img[height * 4 * stride - 2];
                  result_img[index + 3] = source_img[height * 4 * stride - 1];
                }
              if (jy >= ver_low && jy < ver_up)
                {
                  /* right side between top and bottom corner is going to be resampled over the pixels in that line */
                  result_img[index + 0] = right[((jy - ver_low) * hori_left + ix - hori_right) * 4 + 0];
                  result_img[index + 1] = right[((jy - ver_low) * hori_left + ix - hori_right) * 4 + 1];
                  result_img[index + 2] = right[((jy - ver_low) * hori_left + ix - hori_right) * 4 + 2];
                  result_img[index + 3] = right[((jy - ver_low) * hori_left + ix - hori_right) * 4 + 3];
                }
            }
          if (ix >= hori_left && ix < hori_right)
            {
              if (jy >= ver_low && jy < ver_up)
                {
                  /* rest of the image is going to be resampled linear above the complete source_img */
                  result_img[index + 0] = result[((jy - ver_low) * hori_cal + ix - hori_left) * 4 + 0];
                  result_img[index + 1] = result[((jy - ver_low) * hori_cal + ix - hori_left) * 4 + 1];
                  result_img[index + 2] = result[((jy - ver_low) * hori_cal + ix - hori_left) * 4 + 2];
                  result_img[index + 3] = result[((jy - ver_low) * hori_cal + ix - hori_left) * 4 + 3];
                }
              if (jy < ver_low)
                {
                  /* top side between left and right corner is going to be resampled over the pixels in that line */
                  result_img[index + 0] = top[(jy * hori_cal + ix - hori_left) * 4 + 0];
                  result_img[index + 1] = top[(jy * hori_cal + ix - hori_left) * 4 + 1];
                  result_img[index + 2] = top[(jy * hori_cal + ix - hori_left) * 4 + 2];
                  result_img[index + 3] = top[(jy * hori_cal + ix - hori_left) * 4 + 3];
                }
              if (jy >= ver_up)
                {
                  /* bottom side between left and right corner is going to be resampled over the pixels in that line */
                  result_img[index + 0] = bot[((jy - ver_up) * hori_cal + ix - hori_left) * 4 + 0];
                  result_img[index + 1] = bot[((jy - ver_up) * hori_cal + ix - hori_left) * 4 + 1];
                  result_img[index + 2] = bot[((jy - ver_up) * hori_cal + ix - hori_left) * 4 + 2];
                  result_img[index + 3] = bot[((jy - ver_up) * hori_cal + ix - hori_left) * 4 + 3];
                }
            }
        }
    }
}

static void resample_rgba_nearest(const unsigned char *source_img, unsigned char *result_img, size_t width,
                                  size_t height, size_t new_w, size_t new_h, size_t stride, int swapx, int swapy)
{
  int red, green, blue, alpha;
  int i, j, ix, iy, rgb;
  for (j = 0; j < new_h; j++)
    {
      iy = (int)height * j / (int)new_h;
      if (swapy)
        {
          iy = (int)height - 1 - iy;
        }
      for (i = 0; i < new_w; i++)
        {
          ix = (int)width * i / (int)new_w;
          if (swapx)
            {
              ix = (int)width - 1 - ix;
            }
          rgb = source_img[iy * 4 * new_w + ix];
          red = (rgb & 0xff);
          green = (rgb & 0xff00) >> 8;
          blue = (rgb & 0xff0000) >> 16;
          alpha = (rgb & 0xff000000) >> 24;
          result_img[j * 4 * new_w + i * 4 + 0] = (unsigned char)red;
          result_img[j * 4 * new_w + i * 4 + 1] = (unsigned char)green;
          result_img[j * 4 * new_w + i * 4 + 2] = (unsigned char)blue;
          result_img[j * 4 * new_w + i * 4 + 3] = (unsigned char)alpha;
        }
    }
}

/*!
 * Method switches between different resample filters depending on the resample flag status.
 *
 * \param[in] source_img Data array of the original image
 * \param[out] result_img Data array with resampled values
 * \param[in] width Size of a row in source_img
 * \param[in] height Number of rows in source_img
 * \param[in] new_w Size of a row in result_img
 * \param[in] new_h Number of rows in result_img
 * \param[in] stride Stride of source_img
 * \param[in] swapx True if x coordinates are swapped in source_img
 * \param[in] swapy True if y coordinates are swapped in source_img
 *
 * Calculate new image data with linear or lanczos interpolation to receive smoother results.
 * This way the source_img can be up or downscaled. It is an alternative for the nearest neighbour resampling.
 */
void gks_resample(const unsigned char *source_img, unsigned char *result_img, size_t width, size_t height, size_t new_w,
                  size_t new_h, size_t stride, int swapx, int swapy)
{
  int method[1] = {0};
  gks_inq_resample_method(method);
  if (method[0] == GKS_K_RESAMPLE_LINEAR)
    {
      resample_rgba_linear(source_img, result_img, width, height, new_w, new_h, stride, 1, swapx, swapy);
    }
  else if (method[0] == GKS_K_RESAMPLE_LANCZOS)
    {
      resample_rgba_lanczos(source_img, result_img, width, height, new_w, new_h, stride, 3, 0, 255, swapx, swapy);
    }
  else if (method[0] == GKS_K_RESAMPLE_NEAREST)
    {
      resample_rgba_nearest(source_img, result_img, width, height, new_w, new_h, stride, swapx, swapy);
    }
  else
    {
      gks_fatal_error(
          "Only GKS_K_RESAMPLE_LINEAR, GKS_K_RESAMPLE_NEAREST and GKS_K_RESAMPLE_LANCZOS are valid for resample\n");
    }
}
