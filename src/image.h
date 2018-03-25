#ifndef __IMAGE_H__
#define __IMAGE_H__

void saveBMP(const char* filename, const unsigned char* image, int width, int height);
void initilize_img(unsigned char* image, int img_width, int img_height);
void update_img(unsigned char* image, double *p, int img_width, int img_height);

#endif