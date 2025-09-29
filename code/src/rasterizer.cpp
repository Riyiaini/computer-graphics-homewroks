#include "rasterizer.h"
#include <chrono>

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)


    sample_buffer[y * width + x] = c;
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling

    // auto start = std::chrono::high_resolution_clock::now();

    /* solution 1 */
   /*  Vector2D v1{x2 - x1, y2 - y1};
    Vector2D v2{x0 - x2, y0 - y2};

    if(cross(v1, v2) > 0) {
      std::swap(x1, x2);
      std::swap(y1, y2);
    }
    int minx = std::floor(std::min(std::min(x0, x1), x2));
    int miny = std::floor(std::min(std::min(y0, y1), y2));
    int maxx = std::floor(std::max(std::max(x0, x1), x2));
    int maxy = std::floor(std::max(std::max(y0, y1), y2));

    bool find = false;
    for (int y = miny; y <= maxy; y++) {
      for (int x = minx; x <= maxx; x++) {
      
        if (inside_triangle(x0, y0, x1, y1, x2, y2, x+0.5f, y+0.5f)) {
          find = true;
          rasterize_point(x, y, color);
        } else {
        } 
      }
    } */

    /* solution 2*/
    int miny = std::floor(std::min(std::min(y0, y1), y2));
    int maxy = std::ceil(std::max(std::max(y0, y1), y2)); 

    for (int y = miny; y <= maxy; y++) {
      float minx = std::numeric_limits<float>::max();
      float maxx = std::numeric_limits<float>::lowest();
      bool found_intersection = false;
      float y_m = y + 0.5;
      if (y_m >= std::min(y0, y1) && y_m <= std::max(y0, y1) && y1 != y0) {
        float x = (x0 * (y1 - y_m) + x1 * (y_m - y0)) / (y1 - y0);
        minx = std::min(minx, x);
        maxx = std::max(maxx, x);
        found_intersection = true;
      }
      if (y_m >= std::min(y1, y2) && y_m <= std::max(y1, y2) && y2 != y1) {
        float x = (x1 * (y2 - y_m) + x2 * (y_m - y1)) / (y2 - y1);
        minx = std::min(minx, x);
        maxx = std::max(maxx, x);
        found_intersection = true;
      }
      if (y_m >= std::min(y2, y0) && y_m <= std::max(y2, y0) && y0 != y2) {
        float x = (x2 * (y0 - y_m) + x0 * (y_m - y2)) / (y0 - y2);
        minx = std::min(minx, x);
        maxx = std::max(maxx, x);
        found_intersection = true;
      }
      if (found_intersection) {
        int start_x = std::round(minx);
        int end_x = std::round(maxx);
        for (int x = start_x; x < end_x; x++) {
            rasterize_point(x, y, color);
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



  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle




  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;


    this->sample_buffer.resize(width * height, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize(width * height, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support


    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color col = sample_buffer[y * width + x];

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }

  }

  bool RasterizerImp::inside_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    float px, float py)
  {
    float cp1 = (x1 - x0) * (py - y0) - (y1 - y0) * (px - x0);
    float cp2 = (x2 - x1) * (py - y1) - (y2 - y1) * (px - x1);
    float cp3 = (x0 - x2) * (py - y2) - (y0 - y2) * (px - x2);

    if ((x1 > x0 && cp1 > 0) || cp1 >= 0) {
      return false;
    }
    if ((x2 > x1 && cp2 > 0) || cp2 >= 0) {
      return false;
    }
    if ((x0 > x2 && cp3 > 0) || cp3 >= 0) {
      return false;
    }

    return true;
  }


  Rasterizer::~Rasterizer() { }


}// CGL
