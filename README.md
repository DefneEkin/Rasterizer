# Rasterizer

A C++ software rasterizer implementing key stages of the **Forward Rendering Pipeline**.  
This project renders 3D scenes defined in an XML-like input file by applying modeling transformations, camera projection, triangle rasterization, interpolation, clipping, and depth buffering.

The rasterizer supports both **wireframe** and **solid** rendering modes, backface culling, and multiple camera configurations.

---

## Features

- **Modeling Transformations:** translation, scaling, rotation (applied per mesh)
- **Viewing Transformation:** camera setup using position, gaze, up vector
- **Projection:** supports both **perspective** and **orthographic** projection
- **Triangle Rasterization:**
  - Midpoint algorithm for edge drawing (wireframe)
  - Barycentric interpolation for solid shading
- **Color Interpolation:** per-vertex RGB colors smoothly interpolated over surfaces
- **Clipping (wireframe mode):** implement Liang–Barsky
- **Backface Culling:** optional, configurable in input
- **Depth Buffering:** ensures correct visibility when triangles overlap
- **Multiple Cameras:** renders one output PPM file per camera


---

## Building the Executable

The project includes a `Makefile` for easy compilation.

```bash
make
```

This command will create an executable named **`rasterizer`**.

---

## Running the Rasterizer

The executable takes one argument: the path to the input scene file.


```bash
./rasterizer <scene_file.xml>
```


Each camera defined in the file produces its own `.ppm` output image.

---

## Input File Overview

Input files are structured similarly to XML and include:

- **BackgroundColor** — default color when no geometry is drawn  
- **Culling** — enables/disables backface culling  
- **Cameras** — position, gaze, up, projection type, planes, resolution, output name  
- **Vertices** — ID, position, and RGB color  
- **Transformations:**  
  - Translation (tx, ty, tz)  
  - Scaling (sx, sy, sz)  
  - Rotation (angle and axis)  
- **Meshes:**  
  - wireframe / solid mode  
  - list and order of transformations  
  - list of triangles defined by vertex IDs  

---

## Rendering Pipeline Summary

1. **Modeling Transformations**  
   Apply translation, scaling, and rotations to mesh vertices.

2. **Viewing Transformation**  
   Convert coordinates into the camera’s coordinate system.

3. **Projection**  
   Project 3D coordinates to the image plane (orthographic or perspective).

4. **Clipping (wireframe only)**  
   Clip lines using a parametric or region-based algorithm.

5. **Rasterization**  
   - Wireframe: draw triangle edges  
   - Solid: fill triangle surfaces using barycentric coordinates  

6. **Depth Test**  
   Use a depth buffer to keep only the nearest fragment.

7. **Color Output**  
   Write final image as a `.ppm` file.
