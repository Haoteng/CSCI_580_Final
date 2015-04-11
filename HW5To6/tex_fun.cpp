/* Texture functions for cs580 GzLib	*/
#include	<stdlib.h>
#include	<stdio.h>
#include	<math.h>
#include	"gz.h"
#include "PooyanToolset.h"

GzColor *image;
int xs, ys;
int reset = 1;

/* Image texture function */
int tex_fun(float u, float v, GzColor color) {

	unsigned char pixel[3];
	unsigned char dummy;
	char foo[8];
	int i;
	FILE *fd;

	if (reset) { /* open and load texture file */
		fd = fopen("texture", "rb");
		if (fd == NULL) {
			fprintf(stderr, "texture file not found\n");
			return GZ_FAILURE;
		}
		fscanf(fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
		image = (GzColor*) malloc(sizeof(GzColor) * (xs + 1) * (ys + 1));
		if (image == NULL) {
			fprintf(stderr, "malloc for texture image failed\n");
			return GZ_FAILURE;
		}

		for (i = 0; i < xs * ys; i++) { /* create array of GzColor values */
			fread(pixel, sizeof(pixel), 1, fd);
			image[i][RED] = (float) ((int) pixel[RED]) * (1.0 / 255.0);
			image[i][GREEN] = (float) ((int) pixel[GREEN]) * (1.0 / 255.0);
			image[i][BLUE] = (float) ((int) pixel[BLUE]) * (1.0 / 255.0);
		}

		reset = 0; /* init is done */
		fclose(fd);
	}

	/* bounds-test u,v to make sure nothing will overflow image array bounds */
	float zero = 0, one = 1;
	PooyanToolset::fitValue(u, zero, one);
	PooyanToolset::fitValue(v, zero, one);
	/* determine texture cell corner values and perform bilinear interpolation */
	int x = u * (xs - 1);
	int y = v * (ys - 1);
	int x2 = x + ((x + 1) < xs ? 1 : 0);
	int y2 = y + ((y + 1) < ys ? 1 : 0);
	GzColor A;
	PooyanToolset::copyFromVoid(A, image[x + y * xs]);
	GzColor B;
	PooyanToolset::copyFromVoid(B, image[x2 + y * xs]);
	GzColor C;
	PooyanToolset::copyFromVoid(C, image[x2 + y2 * xs]);
	GzColor D;
	PooyanToolset::copyFromVoid(D, image[x + y2 * xs]);

	/* set color to interpolated GzColor value and return */
	for (int i = 0; i < 3; i++)
		color[i] = PooyanToolset::bilinear(u * (xs - 1) - x, v * (ys - 1) - y,
				A[i], B[i], C[i], D[i]);
	return GZ_SUCCESS;
}


/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color) {
	//Julia set
	float x_r=u;
	float x_i=v;
	float c_i=-.2;
	float c_r=.8;
	for (int n = 0; n < 10; n++) {
		float r= x_r*x_r -x_i*x_i+c_r;
		float i= 2*x_r*x_i + c_i;
		x_r=r;
		x_i=i;
		if (i*i+r*r>2)
			break ;
	}
	float l = sqrt(x_r * x_r + x_i * x_i)/2;
	float zero=0,one=1;
	PooyanToolset::fitValue(l,zero,one);
	for (int i=0;i<3;i++)
		color[i]=l*2/(4-i);
	return GZ_SUCCESS;
}

