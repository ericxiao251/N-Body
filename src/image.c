#include "define.h"
#include <assert.h>

void saveBMP(const char* filename, const unsigned char* result, int w, int h){
    printf("saving image...\n");
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
    printf("finished saving!\n");
}


void initilize_img(unsigned char* image, int img_width, int img_height) {
	int i, j;

	for (i = 0; i < img_height; ++i) {
		for (j = 0; j < img_width; ++j) {
			image[(j * img_width + i) * 3 + 0] = (unsigned char)0;
			image[(j * img_width + i) * 3 + 1] = (unsigned char)0;
			image[(j * img_width + i) * 3 + 2] = (unsigned char)0;
		}
	}
}

void update_img(unsigned char* image, double *p, int img_width, int img_height) {
	// POS_MAX_X POS_MIN_X POS_MAX_Y POS_MIN_Y
	int x = (int)(X_RNG * p[POS_X_COL] / (double)img_width);
	int y = (int)(Y_RNG * p[POS_Y_COL] / (double)img_height);
	if (x < 0 || y < 0 || x >= img_width || y >= img_height) {
		// Out of range
		return;
	}

	unsigned char r, g, b;
	if (massLightMin <= p[WEIGHT_COL] <= massLightMax) {
		r = LIGHT_R;
		g = LIGHT_G;
		b = LIGHT_B;
	} else if (massMediumMin <= p[WEIGHT_COL] <= massMediumMax) {
		r = MEDIUM_R;
		g = MEDIUM_G;
		b = MEDIUM_B;
	} else {
		r = HEAVY_R;
		g = HEAVY_G;
		b = HEAVY_B;
	}

	image[(y * img_width + x) * 3 + 0] = r;
	image[(y * img_width + x) * 3 + 1] = g;
	image[(y * img_width + x) * 3 + 2] = b;
}