#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

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
