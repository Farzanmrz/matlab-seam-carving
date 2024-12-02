# MATLAB Seam Carving for Image Resizing and Energy Optimization

This repository contains a MATLAB implementation of seam carving for image resizing and energy optimization using neural network-inspired bilinear transformations.

## Objective

This project explores the application of seam carving techniques to resize images by iteratively removing pixels with the least energy. The energy optimization step enhances the quality of the resized image, while neural network-inspired bilinear interpolation is used for smooth resizing.

## Methods Used

### Seam Carving

The project employs seam carving, an advanced content-aware image resizing technique. Energy maps are computed to identify the least important pixels, which are removed iteratively.

### Bilinear and Nearest Neighbor Interpolation

- Nearest Neighbor Interpolation: A basic resampling method applied during image resizing.
- Bilinear Interpolation: Used for smoother resizing by averaging pixel values.

### Energy Optimization

Energy optimization involves computing gradient magnitudes to identify important features in the image, ensuring quality is preserved during resizing.

### Video Output

The project also visualizes the seam carving process by generating an `.mp4` video that highlights the gradual removal of seams.

## Files

- `main.m`: The main MATLAB script for processing images, resizing, and video generation.
- `img.jpeg` and `img2.jpeg`: Example input images for testing.
- `LICENSE`: MIT license.
- `README.md`: Project description and instructions.

## Contact

* Farzan Mirza: [farzan.mirza@drexel.edu](mailto:farzan.mirza@drexel.edu) | [LinkedIn](https://www.linkedin.com/in/farzan-mirza13/)
