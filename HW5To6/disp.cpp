/* 
*  disp.cpp -- definition file for Display
*  USC csci 580 
*/

#include "gz.h"
#include "disp.h"
#include <stdio.h>
#include "PooyanToolset.h"


#ifndef MIN_INTENSITY
#define MIN_INTENSITY 0
#endif

#ifndef MAX_INTENSITY
#define MAX_INTENSITY 4095
#endif

#ifndef MAX_INT
#define MAX_INT 2147483647
#endif


//This function checks if w and h are in the max range defined in disp.h
int isWidthHeightinMAXRange(int width,int height){
    return PooyanToolset::isInRange(width,0,MAXXRES)
        && PooyanToolset::isInRange(height,0,MAXYRES) ;
}

//This function checks if a pixel (x,y) is in a display (frame)
int isInDisplay(GzDisplay *display, int x , int y) {
    return PooyanToolset::isInRange(x,0,display->xres-1)
        && PooyanToolset::isInRange(y,0,display->yres-1);
}

// This funciton fits the value for RGBA in the frame
void fitIntensity(GzIntensity& x){
    PooyanToolset::fitValue(x,(GzIntensity)MIN_INTENSITY, (GzIntensity)MAX_INTENSITY);
}

int GzNewFrameBuffer(char** framebuffer, int width, int height)
{
/* create a framebuffer:
 -- allocate memory for framebuffer : (sizeof)GzPixel x width x height
 -- pass back pointer
*/
    if (! isWidthHeightinMAXRange(width,height) ) // check if the heigh or width are out of boundries
        return GZ_FAILURE;
    (*framebuffer)= new char[sizeof(GzPixel) * width * height];
    return GZ_SUCCESS;
}

int GzNewDisplay(GzDisplay	**display, GzDisplayClass dispClass, int xRes, int yRes)
{

/* create a display:
  -- allocate memory for indicated class and resolution
  -- pass back pointer to GzDisplay object in display
*/
    if (! isWidthHeightinMAXRange(xRes,yRes) )
        return GZ_FAILURE;
    (*display)= new GzDisplay () ;
    (*display)->xres=xRes;
    (*display)->yres=yRes;
    (*display)->dispClass=dispClass;
    (*display)->fbuf=new GzPixel[xRes*yRes];
	return GZ_SUCCESS;
}

int GzFreeDisplay(GzDisplay	*display)
{
    /* clean up, free memory */
    if (display){
        if (display->fbuf)
            delete [] display->fbuf;
        delete display;
    }
	return GZ_SUCCESS;
}

int GzGetDisplayParams(GzDisplay *display, int *xRes, int *yRes, GzDisplayClass	*dispClass)
{
    /* pass back values for an open display */
    (*dispClass) = display->dispClass;
    (*xRes) = display->xres;
    (*yRes) = display->yres;
	return GZ_SUCCESS;
}

int GzInitDisplay(GzDisplay	*display)
{
    /* set everything to some default values - start a new frame */
    for (int i=0; i< display->xres;i++)
        for (int j=0; j< display->yres;j++)
            GzPutDisplay(display, i, j, MAX_INTENSITY, MAX_INTENSITY, MAX_INTENSITY, 1, MAX_INT);// in fact we are refilling the display with a default pixel
    return GZ_SUCCESS;
}

int GzResetDisplay(GzDisplay	*display)
{
    /* set everything to some default values - start a new frame */
    for (int i=0; i< display->xres;i++)
        for (int j=0; j< display->yres;j++)
            GzPutDisplay(display, i, j, MIN_INTENSITY, MIN_INTENSITY, MIN_INTENSITY, 0, MAX_INT);// in fact we are refilling the display with a default pixel
    return GZ_SUCCESS;
}



int GzPutDisplay(GzDisplay *display, int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
    /* write pixel values into the display */
    if (!isInDisplay(display,i,j))
        return GZ_FAILURE;
    fitIntensity(r);
    fitIntensity(g);
    fitIntensity(b);
    fitIntensity(a);

	GzPixel * pixel = &(display->fbuf[ARRAY(i, j)]);
	pixel->blue = b;
	pixel->green = g;
	pixel->red = r;
	pixel->alpha = a;
	pixel->z = z;


    return GZ_SUCCESS;
}


int GzGetDisplay(GzDisplay *display, int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
	/* pass back pixel value in the display */
	/* check display class to see what vars are valid */
    if (!isInDisplay(display,i,j))
        return GZ_FAILURE;

    GzPixel * pixel = &(display->fbuf[ARRAY(i,j)]);
    (*r)=pixel->red;
    (*g)=pixel->green;
    (*b)=pixel->blue;
    (*r)=pixel->red;
    (*a)=pixel->alpha;
    (*z)=pixel->z;

    return GZ_SUCCESS;
}


int GzFlushDisplay2File(FILE* outfile, GzDisplay *display)
{
    /* write pixels to ppm file based on display class -- "P6 %d %d 255\r" */

    fprintf(outfile,"P3 %d %d %d\r" , display->xres, display->yres,MAX_INTENSITY);
    for (int i=0;i<display->xres*display->yres;i++){
    	if ((i+1)%display->yres ==0 )
    		fprintf (outfile,"\r");
    	GzPixel pixel = display->fbuf[i];
    	fprintf (outfile,"%d %d %d " ,pixel.red,pixel.green,pixel.blue);
    }
	return GZ_SUCCESS;
}


int GzFlushDisplay2FrameBuffer(char* framebuffer, GzDisplay *display)
{

	/* write pixels to framebuffer:
		- Put the pixels into the frame buffer
		- Caution: store the pixel to the frame buffer as the order of blue, green, and red
		- Not red, green, and blue !!!
	*/
	for ( int i=0 ; i<display->xres*display->yres; i++){
		GzPixel pixel = display->fbuf[i];
		framebuffer[i*3]=pixel.blue;
		framebuffer[i*3+1]=pixel.green;
		framebuffer[i*3+2]=pixel.red;
	}
	return GZ_SUCCESS;
}

int GzAddDisplay(GzDisplay **base,GzDisplay* second, float weight){

	int ps = (*base)->xres*(*base)->yres;
	for (int c=0;c<ps;c++){
		(*base)->fbuf[c].red+=second->fbuf[c].red*weight;
		(*base)->fbuf[c].green+=second->fbuf[c].green*weight;
		(*base)->fbuf[c].blue+=second->fbuf[c].blue*weight;
		(*base)->fbuf[c].alpha+=second->fbuf[c].alpha*weight;
		(*base)->fbuf[c].z+=second->fbuf[c].z*weight;
	}
	return GZ_SUCCESS;
}
