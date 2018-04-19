#include "define.h"
#include <assert.h>


void saveBMP(const char* filename, const unsigned char* result, int w, int h){
    //printf("saving image...\n");
	FILE *f;
	unsigned char *img = NULL;
	int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int

	unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
	unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
	unsigned char bmppad[3] = {0,0,0};

	bmpfileheader[ 2] = (unsigned char)(filesize    );
	bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
	bmpfileheader[ 4] = (unsigned char)(filesize>>16);
	bmpfileheader[ 5] = (unsigned char)(filesize>>24);

	bmpinfoheader[ 4] = (unsigned char)(       w    );
	bmpinfoheader[ 5] = (unsigned char)(       w>> 8);
	bmpinfoheader[ 6] = (unsigned char)(       w>>16);
	bmpinfoheader[ 7] = (unsigned char)(       w>>24);
	bmpinfoheader[ 8] = (unsigned char)(       h    );
	bmpinfoheader[ 9] = (unsigned char)(       h>> 8);
	bmpinfoheader[10] = (unsigned char)(       h>>16);
	bmpinfoheader[11] = (unsigned char)(       h>>24);

	f = fopen(filename,"wb");
	fwrite(bmpfileheader,1,14,f);
	fwrite(bmpinfoheader,1,40,f);


	img = (unsigned char *)malloc(3*w);
	assert(img);

	int i,j;
	for(j=0; j<h; j++)
	{
	    for(i=0; i<w; i++)
		{
            img[i*3+0] = result[(j*w+i)*3+0];
            img[i*3+1] = result[(j*w+i)*3+1];
            img[i*3+2] = result[(j*w+i)*3+2];
		}
		fwrite(img,3,w,f);
	    fwrite(bmppad,1,(4-(w*3)%4)%4,f);
	}
	fclose(f);
    //printf("finished saving!\n");
}

void initilize_img(unsigned char* image, int img_width, int img_height) {
	int i, j;

	for (i = 0; i < img_height; ++i) {
		for (j = 0; j < img_width; ++j) {
			image[(i * img_width + j) * 3 + 0] = (unsigned char)0;
			image[(i * img_width + j) * 3 + 1] = (unsigned char)0;
			image[(i * img_width + j) * 3 + 2] = (unsigned char)0;
		}
	}
}

void paint_square(unsigned char* image, int i, int j, int img_width, int img_height,
				int filterRad, unsigned char r, unsigned char g, unsigned char b) {
	int ii, jj;

	int min_j = (j - filterRad < 0) ? 0 : j - filterRad;
    int max_j = (j + filterRad >= img_width) ? img_width - 1 : j + filterRad;
    int min_i = (i - filterRad < 0) ? 0 : i - filterRad;
    int max_i = (i + filterRad >= img_height) ? img_height - 1 : i + filterRad;
    for (ii = min_i; ii <= max_i; ++ii) {
    	for (jj = min_j; jj <= max_j; ++jj) {
				image[(ii * img_width + jj) * 3 + 0] = b;
				image[(ii * img_width + jj) * 3 + 1] = g;
				image[(ii * img_width + jj) * 3 + 2] = r;
    	}
    }
}

int update_img(unsigned char* image, double *p, int img_width, int img_height) {
	// POS_MAX_X POS_MIN_X POS_MAX_Y POS_MIN_Y
	int x = (int)(img_width * p[POS_X_COL] / (double)X_RNG);
	int y = (int)(img_height * p[POS_Y_COL] / (double)Y_RNG);
	if (x < 0 || y < 0 || x >= img_width || y >= img_height) {
		// Out of range
		return 0;
	}

	unsigned char r, g, b;
	r = p[TYPE_COL] == LIGHT ? LIGHT_R :
		p[TYPE_COL] == MEDIUM ? MEDIUM_R :
		p[TYPE_COL] == HEAVY ? HEAVY_R : 0;
	g = p[TYPE_COL] == LIGHT ? LIGHT_G :
		p[TYPE_COL] == MEDIUM ? MEDIUM_G :
		p[TYPE_COL] == HEAVY ? HEAVY_G : 0;
	b = p[TYPE_COL] == LIGHT ? LIGHT_B :
		p[TYPE_COL] == MEDIUM ? MEDIUM_B :
		p[TYPE_COL] == HEAVY ? HEAVY_B : 0;

	paint_square(image, y, x, img_width, img_height, PARTICLE_R, r, g, b);
	//image[(y * img_width + x) * 3 + 0] = r;
	//image[(y * img_width + x) * 3 + 1] = g;
	//image[(y * img_width + x) * 3 + 2] = b;
	return 1;
}
