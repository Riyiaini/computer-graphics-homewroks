#include "rasterizer.h"
#include <chrono>
#include "texture.h"
#include <immintrin.h>

using namespace std;

namespace CGL
{

  void Coords_computer::get_denom(float x0, float y0, float x1, float y1, float x2, float y2)
  {
    this->x0 = x0;
    this->y0 = y0;

    v0 = Vector2D(x1 - x0, y1 - y0);
    v1 = Vector2D(x2 - x0, y2 - y0);

    d00 = v0.x * v0.x + v0.y * v0.y; // v0·v0
    d01 = v0.x * v1.x + v0.y * v1.y; // v0·v1
    d11 = v1.x * v1.x + v1.y * v1.y; // v1·v1

    denom = d00 * d11 - d01 * d01;
  }

  Coords Coords_computer::compute_coords(float x, float y)
  {
    Vector2D v2(x - x0, y - y0);

    float d20 = v2.x * v0.x + v2.y * v0.y; // v2·v0
    float d21 = v2.x * v1.x + v2.y * v1.y; // v2·v1

    if (std::abs(denom) < 1e-10f)
    {
      return {0, 0, 0};
    }

    struct Coords coords;
    coords.beta = (d11 * d20 - d01 * d21) / denom;
    coords.gamma = (d00 * d21 - d01 * d20) / denom;
    coords.alpha = 1.0f - coords.beta - coords.gamma;

    return coords;
  }

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
                               size_t width, size_t height,
                               unsigned int sample_rate)
  {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;
    switch (sample_rate)
    {
    case 4:
      this->rate = 2;
      break;
    case 9:
      this->rate = 3;
      break;
    case 16:
      this->rate = 4;
      break;
    default:
      this->rate = 1;
    }

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c)
  {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)

    if (rate == 1)
    {
      sample_buffer[y * width + x] = c;
    }
    else
    {
      // point in super sampling space
      /* int pixel_x = x / rate;
      int pixel_y = y / rate;

      // offset in the super sampling box
      int sample_x = x % rate;
      int sample_y = y % rate;

      int pixel_index = pixel_y * width + pixel_x;
      int sample_index = sample_y * rate + sample_x;
      sample_buffer[pixel_index * sample_rate + sample_index] = c; */

      int ry = y % rate;
      int rx = x % rate;
      int ny = y - ry;
      int nx = x - rx;
      sample_buffer[(ny * width + nx + ry) * rate + rx] = c;
    }
  }

  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c, bool buffer[])
  {
    int pixel_index = y * width + x;
    int pixel_offset = pixel_index * sample_rate;

    if (buffer == NULL)
    {
      for (int i = 0; i < sample_rate; i++)
      {
        sample_buffer[pixel_offset + i] = c;
      }
    }
    else
    {
      for (int i = 0; i < sample_rate; i++)
      {
        if (buffer[i])
        {
          sample_buffer[pixel_offset + i] = c;
        }
      }
    }
  }

  void RasterizerImp::fill_pixel(size_t x, size_t y, Color colors[], bool buffer[])
  {

    int pixel_index = y * width + x;
    int pixel_offset = pixel_index * sample_rate;

    if (buffer == NULL)
    {
      for (int i = 0; i < sample_rate; i++)
      {
        sample_buffer[pixel_offset + i] = colors[i];
      }
    }
    else
    {
      for (int i = 0; i < sample_rate; i++)
      {
        if (buffer[i])
        {
          sample_buffer[pixel_offset + i] = colors[i];
        }
      }
    }
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color)
  {
    // fill in the nearest pixel
    /* int sx = (int)floor(x * rate);
    int sy = (int)floor(y * rate); */
    int sx = (int)floor(x * rate);
    int sy = (int)floor(y * rate);

    // check bounds
    if (sx < 0 || sx >= width * rate)
      return;
    if (sy < 0 || sy >= height * rate)
      return;

    fill_pixel(sx, sy, color);
  }

  void RasterizerImp::rasterize_pixel(float x, float y, Color color, bool buffer[])
  {
    // fill in the nearest pixel
    /* int sx = (int)floor(x * rate);
    int sy = (int)floor(y * rate); */
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width)
      return;
    if (sy < 0 || sy >= height)
      return;

    fill_pixel(sx, sy, color, buffer);
  }

  void RasterizerImp::rasterize_pixel(float x, float y, Color colors[], bool buffer[])
  {
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width)
      return;
    if (sy < 0 || sy >= height)
      return;

    fill_pixel(sx, sy, colors, buffer);
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
                                     float x1, float y1,
                                     Color color)
  {
    if (x0 > x1)
    {
      swap(x0, x1);
      swap(y0, y1);
    }

    float pt[] = {x0, y0};
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = {1, m};
    int steep = abs(m) > 1;
    if (steep)
    {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0))
    {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0];
      pt[1] += dpt[1];
    }
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
                                         float x1, float y1,
                                         float x2, float y2,
                                         Color color)
  {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling

    // auto start = std::chrono::high_resolution_clock::now();

    int miny = clamp(std::floor(std::min(std::min(y0, y1), y2)), 0.0f, static_cast<float>(height - 1));
    int maxy = clamp(std::ceil(std::max(std::max(y0, y1), y2)), 0.0f, static_cast<float>(height - 1));

    if (sample_rate == 1)
    {

      for (int y = miny; y <= maxy; y++)
      {
        float minx = std::numeric_limits<float>::max();
        float maxx = std::numeric_limits<float>::lowest();
        if (find_intersection(y + 0.5, x0, y0, x1, y1, x2, y2, minx, maxx))
        {
          int start_x = clamp(static_cast<int>(std::round(minx)), 0, (int)width - 1);
          int end_x = clamp(static_cast<int>(std::round(maxx)), 0, (int)width - 1);
          for (int x = start_x; x < end_x; x++)
          {
            rasterize_point(x, y, color);
          }
        }
      }
      return;
    }

    // TODO: Task 2: Update to implement super-sampled rasterization

    float delta = 1.0f / rate;
    float half_delta = delta / 2.0;

    for (int pixel_y = miny; pixel_y <= maxy; pixel_y++)
    {
      float start_x[rate], end_x[rate];
      int lminx = width;
      int lmaxx = 0;
      int rminx = width;
      int rmaxx = 0;

      int index = 0;
      bool find = false;
      for (float y = pixel_y; y < pixel_y + 1; y += delta)
      {
        float minx = width;
        float maxx = 0;
        if (find_intersection(y + half_delta, x0, y0, x1, y1, x2, y2, minx, maxx))
        {
          find = true;
          start_x[index] = minx;
          lminx = std::min(lminx, static_cast<int>(std::floor(minx)));
          lmaxx = std::max(lmaxx, static_cast<int>(std::ceil(minx)));
          end_x[index++] = maxx;
          rminx = std::min(rminx, static_cast<int>(std::floor(maxx)));
          rmaxx = std::max(rmaxx, static_cast<int>(std::ceil(maxx)));
        }
        else
        {
          start_x[index] = std::numeric_limits<float>::max();
          end_x[index++] = std::numeric_limits<float>::lowest();
        }
      }
      // process the left edge
      if (find)
      {
        lminx = clamp(lminx, 0, (int)width - 1);
        lmaxx = clamp(lmaxx, 0, (int)width - 1);
        rminx = clamp(rminx, 0, (int)width - 1);
        rmaxx = clamp(rmaxx, 0, (int)width - 1);
        for (int pixel_x = lminx; pixel_x < lmaxx; pixel_x++)
        {
          bool buffer[sample_rate];
          memset(buffer, 0x00, sizeof(buffer));
          for (int i = 0; i < rate; i++)
          {
            float startx = start_x[i], endx = end_x[i];

            int start = std::max(0, (int)(std::floor((startx - pixel_x - half_delta) * rate)) + 1);
            int end = std::min((int)rate, (int)(std::floor((endx - pixel_x - half_delta) * rate)) + 1);
            for (int j = start; j < end; j++)
            {
              buffer[i * rate + j] = true;
            }
          }
          rasterize_pixel(pixel_x, pixel_y, color, buffer);
        }
        for (int pixel_x = lmaxx; pixel_x < rminx; pixel_x++)
        {
          bool buffer[sample_rate];
          std::fill(buffer, buffer + sample_rate, true);
          rasterize_pixel(pixel_x, pixel_y, color, buffer);
        }
        // process the right edge
        for (int pixel_x = rminx; pixel_x < rmaxx; pixel_x++)
        {
          bool buffer[sample_rate];
          memset(buffer, 0x00, sizeof(buffer));
          for (int i = 0; i < rate; i++)
          {
            float startx = start_x[i], endx = end_x[i];

            int start = std::max(0, (int)(std::floor((startx - pixel_x - half_delta) * rate)) + 1);
            int end = std::min((int)rate, (int)(std::floor((endx - pixel_x - half_delta) * rate)) + 1);
            for (int j = start; j < end; j++)
            {
              buffer[i * rate + j] = true;
            }
            
          }
          rasterize_pixel(pixel_x, pixel_y, color, buffer);
        }
      }
    }
    /* auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "duration: " << duration << " ms" << std::endl; */
  }

  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
                                                            float x1, float y1, Color c1,
                                                            float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle

    int miny = clamp(std::floor(std::min(std::min(y0, y1), y2)), 0.0f, static_cast<float>(height - 1));
    int maxy = clamp(std::ceil(std::max(std::max(y0, y1), y2)), 0.0f, static_cast<float>(height - 1));

    struct Coords_computer cct;
    cct.get_denom(x0, y0, x1, y1, x2, y2);

    if (sample_rate == 1)
    {

      for (int y = miny; y <= maxy; y++)
      {
        float minx = std::numeric_limits<float>::max();
        float maxx = std::numeric_limits<float>::lowest();
        if (find_intersection(y + 0.5, x0, y0, x1, y1, x2, y2, minx, maxx))
        {
          int start_x = clamp(static_cast<int>(std::round(minx)), 0, (int)width - 1);
          int end_x = clamp(static_cast<int>(std::round(maxx)), 0, (int)width - 1);
          for (int x = start_x; x < end_x; x++)
          {
            Coords coords = cct.compute_coords(x + 0.5, y + 0.5);
            Color color = c0 * coords.alpha + c1 * coords.beta + c2 * coords.gamma;
            rasterize_point(x, y, color);
          }
        }
      }
      return;
    }

    float delta = 1.0f / rate;
    float half_delta = delta / 2.0;

    for (int pixel_y = miny; pixel_y <= maxy; pixel_y++)
    {
      float start_x[rate], end_x[rate];
      int lminx = width;
      int lmaxx = 0;
      int rminx = width;
      int rmaxx = 0;

      int index = 0;
      bool found = false;
      for (float y = pixel_y; y < pixel_y + 1; y += delta)
      {
        float minx = width;
        float maxx = 0;
        if (find_intersection(y + half_delta, x0, y0, x1, y1, x2, y2, minx, maxx))
        {
          found = true;
          start_x[index] = minx;
          lminx = std::min(lminx, static_cast<int>(std::floor(minx)));
          lmaxx = std::max(lmaxx, static_cast<int>(std::ceil(minx)));
          end_x[index++] = maxx;
          rminx = std::min(rminx, static_cast<int>(std::floor(maxx)));
          rmaxx = std::max(rmaxx, static_cast<int>(std::ceil(maxx)));
        }
        else
        {
          start_x[index] = std::numeric_limits<float>::max();
          end_x[index++] = std::numeric_limits<float>::lowest();
        }
      }
      // process the left edge
      if (found)
      {
        lminx = clamp(lminx, 0, (int)width - 1);
        lmaxx = clamp(lmaxx, 0, (int)width - 1);
        rminx = clamp(rmaxx, 0, (int)width - 1);
        rmaxx = clamp(rminx, 0, (int)width - 1);
        for (int pixel_x = lminx; pixel_x < lmaxx; pixel_x++)
        {
          bool buffer[sample_rate];
          Color colors[sample_rate];
          memset(buffer, 0x00, sizeof(buffer));
          for (int i = 0; i < rate; i++)
          {
            float startx = start_x[i], endx = end_x[i];

            int start = std::max(0, (int)(std::floor((startx - pixel_x - half_delta) * rate)) + 1);
            int end = std::min((int)rate, (int)(std::floor((endx - pixel_x - half_delta) * rate)) + 1);
            for (int j = start; j < end; j++)
            {
              buffer[i * rate + j] = true;
              Coords coords = cct.compute_coords(pixel_x + j * delta + half_delta, pixel_y + i * delta + half_delta);
              colors[i * rate + j] = c0 * coords.alpha + c1 * coords.beta + c2 * coords.gamma;
            }
            
          }
          rasterize_pixel(pixel_x, pixel_y, colors, buffer);
        }
        for (int pixel_x = lmaxx; pixel_x < rminx; pixel_x++)
        {
          bool buffer[sample_rate];
          std::fill(buffer, buffer + sample_rate, true);
          Color colors[sample_rate];
          for (int i = 0; i < rate; i++)
          {
            for (int j = 0; j < rate; j++)
            {
              buffer[i * rate + j] = true;
              Coords coords = cct.compute_coords(pixel_x + j * delta + half_delta, pixel_y + i * delta + half_delta);
              colors[i * rate + j] = c0 * coords.alpha + c1 * coords.beta + c2 * coords.gamma;
            }
          }
          rasterize_pixel(pixel_x, pixel_y, colors, buffer);
        }
        // process the right edge
        for (int pixel_x = rminx; pixel_x < rmaxx; pixel_x++)
        {
          bool buffer[sample_rate];
          Color colors[sample_rate];
          memset(buffer, 0x00, sizeof(buffer));
          for (int i = 0; i < rate; i++)
          {
            float startx = start_x[i], endx = end_x[i];

            int start = std::max(0, (int)(std::floor((startx - pixel_x - half_delta) * rate)) + 1);
            int end = std::min((int)rate, (int)(std::floor((endx - pixel_x - half_delta) * rate)) + 1);
            for (int j = start; j < end; j++)
            {
              buffer[i * rate + j] = true;
              Coords coords = cct.compute_coords(pixel_x + j * delta + half_delta, pixel_y + i * delta + half_delta);
              colors[i * rate + j] = c0 * coords.alpha + c1 * coords.beta + c2 * coords.gamma;
            }

          }
          rasterize_pixel(pixel_x, pixel_y, colors, buffer);
        }
      }
    }
  }

  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
                                                  float x1, float y1, float u1, float v1,
                                                  float x2, float y2, float u2, float v2,
                                                  Texture &tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.

    int miny = clamp(std::floor(std::min(std::min(y0, y1), y2)), 0.0f, static_cast<float>(height - 1));
    int maxy = clamp(std::ceil(std::max(std::max(y0, y1), y2)), 0.0f, static_cast<float>(height - 1));
    int minx = clamp(std::floor(std::min(std::min(x0, x1), x2)), 0.0f, static_cast<float>(width - 1));
    int maxx = clamp(std::ceil(std::max(std::max(x0, x1), x2)), 0.0f, static_cast<float>(width - 1));
    /* int miny = std::floor(std::min(std::min(y0, y1), y2));
    int maxy = std::ceil(std::max(std::max(y0, y1), y2));
    int minx = std::floor(std::min(std::min(x0, x1), x2));
    int maxx = std::ceil(std::max(std::max(x0, x1), x2));
 */
    if (miny == maxy || minx == maxx)
    {
      return;
    }

    struct Coords_computer cct;
    cct.get_denom(x0, y0, x1, y1, x2, y2);

    if (sample_rate == 1)
    {

      Vector2D uv_buffer[(maxx - minx + 1) * (maxy - miny + 1)];
      bool filled_buffer[(maxx - minx + 1) * (maxy - miny + 1)];
      int buffer_width = maxx - minx + 1;
      int buffer_size = buffer_width * (maxy - miny + 1);
      memset(uv_buffer, 0x00, sizeof(uv_buffer));
      memset(filled_buffer, 0x00, sizeof(filled_buffer));

      for (int y = miny; y < maxy; y++)
      {
        float leftx = std::numeric_limits<float>::max();
        float rightx = std::numeric_limits<float>::lowest();
        if (find_intersection(y + 0.5, x0, y0, x1, y1, x2, y2, leftx, rightx))
        {
          int start_x = clamp(static_cast<int>(std::round(leftx)), 0, (int)width - 1);
          int end_x = clamp(static_cast<int>(std::round(rightx)), 0, (int)width - 1);
          for (int x = start_x; x < end_x; x++)
          {
            int base_index = (y - miny) * buffer_width + (x - minx);
            int indices[3] = {
                base_index,
                base_index + 1,
                base_index + buffer_width};
            constexpr float offsets[3][2] = {
                {0.5f, 0.5f}, // middle
                {1.5f, 0.5f}, // right
                {0.5f, 1.5f}  // down
            };

            for (int i = 0; i < 3; i++)
            {
              if (filled_buffer[indices[i]])
              {
                continue;
              }
              float xx = x + offsets[i][0];
              float yy = y + offsets[i][1];
              Coords coords = cct.compute_coords(xx, yy);
              Vector2D uv;
              uv.x = coords.alpha * u0 + coords.beta * u1 + coords.gamma * u2;
              uv.y = coords.alpha * v0 + coords.beta * v1 + coords.gamma * v2;
              uv_buffer[indices[i]] = uv;
              filled_buffer[indices[i]] = true;
            }

            SampleParams params;
            params.psm = psm;
            params.lsm = lsm;
            params.p_uv = uv_buffer[indices[0]];
            params.p_dx_uv = uv_buffer[indices[1]];
            params.p_dy_uv = uv_buffer[indices[2]];
            params.rate = 1;

            Color color = tex.sample(params);

            rasterize_point(x, y, color);
          }
        }
      }
      return;
    }

    // supser sampling
    float delta = 1.0f / rate;
    float half_delta = delta / 2.0;

    for (int pixel_y = miny; pixel_y < maxy; pixel_y++)
    {
      float start_x[rate], end_x[rate];
      int lminx = width;
      int lmaxx = 0;
      int rminx = width;
      int rmaxx = 0;

      int index = 0;
      bool found = false;
      for (float y = pixel_y; y < pixel_y + 1; y += delta)
      {
        float minx = width;
        float maxx = 0;
        if (find_intersection(y + half_delta, x0, y0, x1, y1, x2, y2, minx, maxx))
        {
          found = true;
          start_x[index] = minx;
          lminx = std::min(lminx, static_cast<int>(std::floor(minx)));
          lmaxx = std::max(lmaxx, static_cast<int>(std::ceil(minx)));
          end_x[index++] = maxx;
          rminx = std::min(rminx, static_cast<int>(std::floor(maxx)));
          rmaxx = std::max(rmaxx, static_cast<int>(std::ceil(maxx)));
        }
        else
        {
          start_x[index] = width;
          end_x[index++] = 0;
        }
      }

      // process the left edge
      if (found)
      {
        lminx = clamp(lminx, 0, (int)width - 1);
        lmaxx = clamp(lmaxx, 0, (int)width - 1);
        rminx = clamp(rmaxx, 0, (int)width - 1);
        rmaxx = clamp(rminx, 0, (int)width - 1);

        bool buffer[sample_rate];                    // which samples are covered
        Color colors[sample_rate];                   // color for each sample
        bool filled_buffer[(rate + 1) * (rate + 1)]; // whether the uv is computed
        Vector2D uv_buffer[(rate + 1) * (rate + 1)]; // uv for each sample and its right and down neighbor

        for (int pixel_x = lminx; pixel_x < lmaxx; pixel_x++)
        {

          memset(buffer, 0x00, sizeof(buffer));
          memset(filled_buffer, 0x00, sizeof(filled_buffer));
          memset(uv_buffer, 0x00, sizeof(uv_buffer));

          for (int i = 0; i < rate; i++)
          {
            float startx = start_x[i], endx = end_x[i];
            int start = std::max(0, (int)(std::floor((startx - pixel_x - half_delta) * rate)) + 1);
            int end = std::min((int)rate, (int)(std::floor((endx - pixel_x - half_delta) * rate)) + 1);
            for (int j = start; j < end; j++)
            {
              buffer[i * rate + j] = true;
              int base_index = i * (rate + 1) + j;
              int indices[3] = {
                  base_index,                  // current
                  base_index + 1,              // right
                  base_index + ((int)rate + 1) // down
              };
              for (int k = 0; k < 3; k++)
              {
                if (filled_buffer[indices[k]])
                {
                  continue;
                }
                float xx = pixel_x + j * delta + half_delta + (k == 1 ? delta : 0);
                float yy = pixel_y + i * delta + half_delta + (k == 2 ? delta : 0);
                Coords coords = cct.compute_coords(xx, yy);
                Vector2D uv;
                uv.x = coords.alpha * u0 + coords.beta * u1 + coords.gamma * u2;
                uv.y = coords.alpha * v0 + coords.beta * v1 + coords.gamma * v2;
                uv_buffer[indices[k]] = uv;
                filled_buffer[indices[k]] = true;
              }
              SampleParams params;
              params.psm = psm;
              params.lsm = lsm;
              params.p_uv = uv_buffer[indices[0]];
              params.p_dx_uv = uv_buffer[indices[1]];
              params.p_dy_uv = uv_buffer[indices[2]];
              params.rate = rate;
              colors[i * rate + j] = tex.sample(params);
            }
          }
          rasterize_pixel(pixel_x, pixel_y, colors, buffer);
        }

        for (int pixel_x = lmaxx; pixel_x < rminx; pixel_x++)
        {
          std::fill(buffer, buffer + sample_rate, true); // fill all samples
          memset(filled_buffer, 0x00, sizeof(filled_buffer));
          memset(uv_buffer, 0x00, sizeof(uv_buffer));

          for (int i = 0; i < rate; ++i)
          {
            for (int j = 0; j < rate; ++j)
            {
              int base_index = i * (rate + 1) + j;
              int indices[3] = {
                  base_index,                  // current
                  base_index + 1,              // right
                  base_index + ((int)rate + 1) // down
              };
              for (int k = 0; k < 3; ++k)
              {
                if (filled_buffer[indices[k]])
                {
                  continue;
                }
                float xx = pixel_x + j * delta + half_delta + (k == 1 ? delta : 0);
                float yy = pixel_y + i * delta + half_delta + (k == 2 ? delta : 0);
                Coords coords = cct.compute_coords(xx, yy);
                Vector2D uv;
                uv.x = coords.alpha * u0 + coords.beta * u1 + coords.gamma * u2;
                uv.y = coords.alpha * v0 + coords.beta * v1 + coords.gamma * v2;
                uv_buffer[indices[k]] = uv;
                filled_buffer[indices[k]] = true;
              }

              SampleParams params;
              params.psm = psm;
              params.lsm = lsm;
              params.p_uv = uv_buffer[indices[0]];
              params.p_dx_uv = uv_buffer[indices[1]];
              params.p_dy_uv = uv_buffer[indices[2]];
              params.rate = rate;

              Color color = tex.sample(params);

              colors[i * rate + j] = color;
            }
          }

          rasterize_pixel(pixel_x, pixel_y, colors, buffer);
        }

        // process the right edge
        for (int pixel_x = rminx; pixel_x < rmaxx; pixel_x++)
        {
          memset(buffer, 0x00, sizeof(buffer));
          memset(filled_buffer, 0x00, sizeof(filled_buffer));
          memset(uv_buffer, 0x00, sizeof(uv_buffer));

          for (int i = 0; i < rate; i++)
          {
            float startx = start_x[i], endx = end_x[i];
  
            int start = std::max(0, (int)(std::floor((startx - pixel_x - half_delta) * rate)) + 1);
            int end = std::min((int)rate, (int)(std::floor((endx - pixel_x - half_delta) * rate)) + 1);
            for (int j = start; j < end; j++)
            {
              buffer[i * rate + j] = true;
              int base_index = i * (rate + 1) + j;
              int indices[3] = {
                  base_index,                  // current
                  base_index + 1,              // right
                  base_index + ((int)rate + 1) // down
              };
              for (int k = 0; k < 3; k++)
              {
                if (filled_buffer[indices[k]])
                {
                  continue;
                }
                float xx = pixel_x + j * delta + half_delta + (k == 1 ? delta : 0);
                float yy = pixel_y + i * delta + half_delta + (k == 2 ? delta : 0);
                Coords coords = cct.compute_coords(xx, yy);
                Vector2D uv;
                uv.x = coords.alpha * u0 + coords.beta * u1 + coords.gamma * u2;
                uv.y = coords.alpha * v0 + coords.beta * v1 + coords.gamma * v2;
                uv_buffer[indices[k]] = uv;
                filled_buffer[indices[k]] = true;
              }
              SampleParams params;
              params.psm = psm;
              params.lsm = lsm;
              params.p_uv = uv_buffer[indices[0]];
              params.p_dx_uv = uv_buffer[indices[1]];
              params.p_dy_uv = uv_buffer[indices[2]];
              params.rate = rate;
              colors[i * rate + j] = tex.sample(params);
            }
          }
          rasterize_pixel(pixel_x, pixel_y, colors, buffer);
        }
      }
    }

    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle
  }

  void RasterizerImp::set_sample_rate(unsigned int rate)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;
    switch (rate)
    {
    case 4:
      this->rate = 2;
      break;
    case 9:
      this->rate = 3;
      break;
    case 16:
      this->rate = 4;
      break;
    default:
      this->rate = 1;
    }

    this->sample_buffer.resize(sample_rate * width * height, Color::White);
  }

  void RasterizerImp::set_framebuffer_target(unsigned char *rgb_framebuffer,
                                             size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;

    this->sample_buffer.resize(sample_rate * width * height, Color::White);
  }

  void RasterizerImp::clear_buffers()
  {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }

  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer()
  {
    // TODO: Task 2: You will likely want to update this function for supersampling support
    auto start = std::chrono::high_resolution_clock::now();

    int BLOCKSIZE = 4; // blocking to impove cache hit rate

    if (rate == 1)
    {
      for (int start_y = 0; start_y < height; start_y += BLOCKSIZE)
      {
        for (int start_x = 0; start_x < width; start_x += BLOCKSIZE)
        {
          int end_y = std::min(start_y + BLOCKSIZE, (int)height);
          int end_x = std::min(start_x + BLOCKSIZE, (int)width);
          for (int off_y = start_y; off_y < end_y; off_y++)
          {
            for (int off_x = start_x; off_x < end_x; off_x++)
            {
              int pixel_index = off_y * width + off_x;
              Color col = sample_buffer[pixel_index];

              for (int k = 0; k < 3; ++k)
              {
                this->rgb_framebuffer_target[3 * pixel_index + k] = (&col.r)[k] * 255;
              }
            }
          }
        }
      }
    }
    else if (rate > 1)
    {
      for (int start_y = 0; start_y < height; start_y += BLOCKSIZE)
      {
        for (int start_x = 0; start_x < width; start_x += BLOCKSIZE)
        {
          int end_y = std::min(start_y + BLOCKSIZE, (int)height);
          int end_x = std::min(start_x + BLOCKSIZE, (int)width);
          for (int off_y = start_y; off_y < end_y; off_y++)
          {
            for (int off_x = start_x; off_x < end_x; off_x++)
            {
              int y = off_y;
              int x = off_x;

              int pixel_index = y * width + x;
              int sample_start = pixel_index * sample_rate;
              float r = 0.0f, g = 0.0f, b = 0.0f;
              const Color *col = &sample_buffer[sample_start];
              for (int off = 0; off < sample_rate; ++off)
              {
                r += col[off].r;
                g += col[off].g;
                b += col[off].b;
              }
              this->rgb_framebuffer_target[3 * pixel_index] = r * 255.0f / static_cast<float>(sample_rate);
              this->rgb_framebuffer_target[3 * pixel_index + 1] = g * 255.0f / static_cast<float>(sample_rate);
              this->rgb_framebuffer_target[3 * pixel_index + 2] = b * 255.0f / static_cast<float>(sample_rate);
            }
          }
        }
      }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "resolve duration: " << duration << " ms" << std::endl;
  }

  bool RasterizerImp::find_intersection(float y_m, float x0, float y0, float x1, float y1, float x2, float y2, float &minx, float &maxx)
  {
    bool found_intersection = false;
    if (y_m >= std::min(y0, y1) && y_m <= std::max(y0, y1) && y1 != y0)
    {
      float x = (x0 * (y1 - y_m) + x1 * (y_m - y0)) / (y1 - y0);
      minx = std::min(minx, x);
      maxx = std::max(maxx, x);
      found_intersection = true;
    }
    if (y_m >= std::min(y1, y2) && y_m <= std::max(y1, y2) && y2 != y1)
    {
      float x = (x1 * (y2 - y_m) + x2 * (y_m - y1)) / (y2 - y1);
      minx = std::min(minx, x);
      maxx = std::max(maxx, x);
      found_intersection = true;
    }
    if (y_m >= std::min(y2, y0) && y_m <= std::max(y2, y0) && y0 != y2)
    {
      float x = (x2 * (y0 - y_m) + x0 * (y_m - y2)) / (y0 - y2);
      minx = std::min(minx, x);
      maxx = std::max(maxx, x);
      found_intersection = true;
    }
    return found_intersection;
  }

  Rasterizer::~Rasterizer() {}

} // CGL
